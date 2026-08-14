//! WGS84 geodetic -> ECEF -> local ENU tangent plane.
//!
//! Tracking is done in metres in a local East/North/Up frame rather than in degrees.
//! Degrees are not a metric space: a 0.01 deg error is 1.1 km of northing but only
//! 0.87 km of easting at Athens, and any filter run directly on lat/lon inherits that
//! distortion. Converting once at the sensor boundary keeps the estimator linear.

pub const WGS84_A: f64 = 6_378_137.0;
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
pub const WGS84_E2: f64 = WGS84_F * (2.0 - WGS84_F);
/// Semi-minor axis.
pub const WGS84_B: f64 = WGS84_A * (1.0 - WGS84_F);
/// Second eccentricity squared, used by the inverse transform.
pub const WGS84_EP2: f64 = (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_B * WGS84_B);

pub const FT_TO_M: f64 = 0.3048;
pub const KT_TO_MS: f64 = 0.514_444_4;
pub const MS_TO_KT: f64 = 1.0 / KT_TO_MS;
pub const FPM_TO_MS: f64 = FT_TO_M / 60.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Geodetic {
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub alt_m: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Enu {
    pub e: f64,
    pub n: f64,
    pub u: f64,
}

impl std::ops::Sub for Enu {
    type Output = Enu;
    fn sub(self, other: Enu) -> Enu {
        Enu {
            e: self.e - other.e,
            n: self.n - other.n,
            u: self.u - other.u,
        }
    }
}

impl std::ops::Add for Enu {
    type Output = Enu;
    fn add(self, other: Enu) -> Enu {
        Enu {
            e: self.e + other.e,
            n: self.n + other.n,
            u: self.u + other.u,
        }
    }
}

/// Scaling a displacement by seconds is how a velocity becomes an extrapolation.
impl std::ops::Mul<f64> for Enu {
    type Output = Enu;
    fn mul(self, k: f64) -> Enu {
        Enu {
            e: self.e * k,
            n: self.n * k,
            u: self.u * k,
        }
    }
}

impl Enu {
    /// Horizontal magnitude, ignoring altitude.
    pub fn horiz(self) -> f64 {
        self.e.hypot(self.n)
    }

    pub fn norm(self) -> f64 {
        (self.e * self.e + self.n * self.n + self.u * self.u).sqrt()
    }

    /// Compass bearing from the frame origin, degrees true.
    pub fn bearing_deg(self) -> f64 {
        let b = self.e.atan2(self.n).to_degrees();
        if b < 0.0 {
            b + 360.0
        } else {
            b
        }
    }
}

/// Earth-centred Earth-fixed position, metres.
#[derive(Clone, Copy, Debug)]
struct Ecef {
    x: f64,
    y: f64,
    z: f64,
}

fn to_ecef(p: Geodetic) -> Ecef {
    let (lat, lon) = (p.lat_deg.to_radians(), p.lon_deg.to_radians());
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = lon.sin_cos();
    // Radius of curvature in the prime vertical.
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    Ecef {
        x: (n + p.alt_m) * cos_lat * cos_lon,
        y: (n + p.alt_m) * cos_lat * sin_lon,
        z: (n * (1.0 - WGS84_E2) + p.alt_m) * sin_lat,
    }
}

/// ECEF back to geodetic, by Bowring's closed form.
///
/// Needed because the ENU "Up" axis is *not* altitude. It is height above the tangent
/// plane at the frame origin, and the Earth curves away from that plane as the square of
/// the range: about 17 km of droop at 460 km. An aircraft at the edge of a 250 NM
/// picture, cruising at 36,000 ft, sits roughly 19,000 ft *below* the Athens tangent
/// plane. Reporting `u` as altitude is therefore wrong by tens of thousands of feet at
/// long range, and wrong in a way that looks plausible at short range.
///
/// This does not change the tracking frame. ENU remains the estimation and screening
/// frame, where the curvature error is common-mode between two nearby aircraft and
/// cancels in their separation. The inverse exists only so the picture can report an
/// altitude that means what it says.
fn from_ecef(q: Ecef) -> Geodetic {
    let p = q.x.hypot(q.y);
    let lon = q.y.atan2(q.x);

    // Bowring's parametric latitude, then one closed-form refinement. Accurate to well
    // under a millimetre for any altitude an aircraft will ever occupy.
    let theta = (q.z * WGS84_A).atan2(p * WGS84_B);
    let (sin_t, cos_t) = theta.sin_cos();
    let lat =
        (q.z + WGS84_EP2 * WGS84_B * sin_t.powi(3)).atan2(p - WGS84_E2 * WGS84_A * cos_t.powi(3));

    let (sin_lat, cos_lat) = lat.sin_cos();
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();

    // `p / cos(lat)` degenerates at the poles, where `p` goes to zero. Nothing in this
    // application flies there, but a transform that silently returns infinity for a
    // legal input is a trap for whoever reuses it next.
    let alt = if cos_lat.abs() > 1e-10 {
        p / cos_lat - n
    } else {
        q.z.abs() / sin_lat.abs() - n * (1.0 - WGS84_E2)
    };

    Geodetic {
        lat_deg: lat.to_degrees(),
        lon_deg: lon.to_degrees(),
        alt_m: alt,
    }
}

/// A local tangent plane anchored at the sensor site.
#[derive(Clone, Copy, Debug)]
pub struct Frame {
    origin: Ecef,
    sin_lat: f64,
    cos_lat: f64,
    sin_lon: f64,
    cos_lon: f64,
    pub site: Geodetic,
}

impl Frame {
    pub fn new(site: Geodetic) -> Self {
        let (sin_lat, cos_lat) = site.lat_deg.to_radians().sin_cos();
        let (sin_lon, cos_lon) = site.lon_deg.to_radians().sin_cos();
        Self {
            origin: to_ecef(site),
            sin_lat,
            cos_lat,
            sin_lon,
            cos_lon,
            site,
        }
    }

    pub fn to_enu(&self, p: Geodetic) -> Enu {
        let q = to_ecef(p);
        let (dx, dy, dz) = (
            q.x - self.origin.x,
            q.y - self.origin.y,
            q.z - self.origin.z,
        );
        Enu {
            e: -self.sin_lon * dx + self.cos_lon * dy,
            n: -self.sin_lat * self.cos_lon * dx - self.sin_lat * self.sin_lon * dy
                + self.cos_lat * dz,
            u: self.cos_lat * self.cos_lon * dx
                + self.cos_lat * self.sin_lon * dy
                + self.sin_lat * dz,
        }
    }

    /// The exact inverse of `to_enu`: a local ENU position back to latitude, longitude
    /// and true altitude above the ellipsoid.
    ///
    /// Used by the picture, not by the tracker. See `from_ecef` for why a separate
    /// altitude is needed at all.
    pub fn to_geodetic(&self, p: Enu) -> Geodetic {
        // Transpose of the forward rotation. For an orthonormal matrix the transpose is
        // the inverse, so this really is exact rather than an approximation of it.
        let dx = -self.sin_lon * p.e - self.sin_lat * self.cos_lon * p.n
            + self.cos_lat * self.cos_lon * p.u;
        let dy = self.cos_lon * p.e - self.sin_lat * self.sin_lon * p.n
            + self.cos_lat * self.sin_lon * p.u;
        let dz = self.cos_lat * p.n + self.sin_lat * p.u;

        from_ecef(Ecef {
            x: self.origin.x + dx,
            y: self.origin.y + dy,
            z: self.origin.z + dz,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn athens() -> Frame {
        Frame::new(Geodetic {
            lat_deg: 37.9838,
            lon_deg: 23.7275,
            alt_m: 0.0,
        })
    }

    #[test]
    fn origin_maps_to_zero() {
        let f = athens();
        let p = f.to_enu(f.site);
        assert!(p.norm() < 1e-6, "origin should be the ENU zero, got {p:?}");
    }

    #[test]
    fn one_degree_north_is_about_111km() {
        let f = athens();
        let p = f.to_enu(Geodetic {
            lat_deg: 38.9838,
            lon_deg: 23.7275,
            alt_m: 0.0,
        });
        assert!(p.n > 110_000.0 && p.n < 112_000.0, "northing was {}", p.n);
        assert!(p.e.abs() < 1.0, "easting should be ~0, was {}", p.e);
    }

    #[test]
    fn easting_is_shorter_than_northing_at_this_latitude() {
        // The whole reason we leave degrees at the sensor boundary.
        let f = athens();
        let north = f
            .to_enu(Geodetic {
                lat_deg: 38.9838,
                lon_deg: 23.7275,
                alt_m: 0.0,
            })
            .n;
        let east = f
            .to_enu(Geodetic {
                lat_deg: 37.9838,
                lon_deg: 24.7275,
                alt_m: 0.0,
            })
            .e;
        assert!(east < north * 0.85, "east {east} vs north {north}");
    }

    #[test]
    fn altitude_becomes_up() {
        let f = athens();
        let p = f.to_enu(Geodetic {
            lat_deg: 37.9838,
            lon_deg: 23.7275,
            alt_m: 10_000.0,
        });
        assert!((p.u - 10_000.0).abs() < 1.0, "up was {}", p.u);
    }

    #[test]
    fn enu_round_trips_back_to_geodetic() {
        let f = athens();
        for (lat, lon, alt) in [
            (37.9838, 23.7275, 0.0),
            (38.5, 24.1, 11_000.0),
            (41.0, 20.0, 12_496.8),
            (34.2, 27.9, 300.0),
        ] {
            let g = Geodetic {
                lat_deg: lat,
                lon_deg: lon,
                alt_m: alt,
            };
            let back = f.to_geodetic(f.to_enu(g));
            assert!(
                (back.lat_deg - lat).abs() < 1e-9,
                "lat {lat} -> {}",
                back.lat_deg
            );
            assert!(
                (back.lon_deg - lon).abs() < 1e-9,
                "lon {lon} -> {}",
                back.lon_deg
            );
            assert!(
                (back.alt_m - alt).abs() < 1e-3,
                "alt {alt} -> {}",
                back.alt_m
            );
        }
    }

    #[test]
    fn up_is_not_altitude_at_long_range() {
        // The defect this inverse exists for, pinned as a test. An aircraft at the edge
        // of the configured 250 NM picture, cruising at 36,089 ft, sits about 19,000 ft
        // *below* the Athens tangent plane. Reporting `u` as altitude was wrong by
        // roughly 55,000 ft — and wrong in the direction that reads as "below sea
        // level", which is at least obviously broken rather than quietly plausible.
        let f = athens();
        let truth_m = 11_000.0;
        let far = Geodetic {
            lat_deg: 37.9838 + 4.16,
            lon_deg: 23.7275,
            alt_m: truth_m,
        };
        let p = f.to_enu(far);

        assert!(p.horiz() > 460_000.0, "range was {}", p.horiz() / 1000.0);
        assert!(
            p.u < -5_000.0,
            "the tangent-plane height should be far below zero here, was {}",
            p.u
        );
        assert!(
            (f.to_geodetic(p).alt_m - truth_m).abs() < 0.01,
            "the inverse must recover true altitude, got {}",
            f.to_geodetic(p).alt_m
        );
    }

    #[test]
    fn curvature_error_cancels_between_a_nearby_pair() {
        // Why the screening math was left alone. The tangent-plane droop is common-mode
        // for two aircraft near each other, so their ENU vertical *separation* is
        // right even where their individual ENU heights are badly wrong. At 460 km a
        // true 305 m separation reads within a few metres — negligible against the
        // 305 m minimum the screen actually tests.
        let f = athens();
        // Same 463 km edge-of-picture geometry as the test above, where the droop is
        // at its worst for this configuration.
        let (lat, lon) = (37.9838 + 4.16, 23.7275);
        let lower = f.to_enu(Geodetic {
            lat_deg: lat,
            lon_deg: lon,
            alt_m: 11_000.0,
        });
        let upper = f.to_enu(Geodetic {
            lat_deg: lat,
            lon_deg: lon + 0.114,
            alt_m: 11_305.0,
        });

        assert!(lower.u < -5_000.0, "both should be far below the plane");
        let enu_separation = upper.u - lower.u;
        assert!(
            (enu_separation - 305.0).abs() < 10.0,
            "pair separation drifted to {enu_separation} m; the screening frame is no \
             longer safe to use unconverted"
        );
    }

    #[test]
    fn bearing_is_compass_convention() {
        let east = Enu {
            e: 1000.0,
            n: 0.0,
            u: 0.0,
        };
        let north = Enu {
            e: 0.0,
            n: 1000.0,
            u: 0.0,
        };
        assert!((east.bearing_deg() - 90.0).abs() < 1e-9);
        assert!(north.bearing_deg().abs() < 1e-9);
    }
}
