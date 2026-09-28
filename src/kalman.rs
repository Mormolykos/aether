//! Constant-velocity Kalman filter, one scalar axis.
//!
//! A 3-D constant-velocity target with diagonal process and measurement noise has a
//! block-diagonal covariance: the East, North and Up blocks never exchange information.
//! Running three independent 2-state filters is therefore numerically identical to one
//! 6-state filter, with no 6x6 matrix inversion, no heap allocation and a fixed
//! instruction count per update. That matters when the same code has to hold a hard
//! deadline on a constrained target.
//!
//! State: [position, velocity]. Measurement: position only.
//!
//! # Covariance in square-root form
//!
//! The covariance is carried as its lower-triangular Cholesky factor `L`, with
//! `P = L Lᵀ`, and never as `P` itself. `L Lᵀ` is symmetric and positive semidefinite for
//! every real `L`, so no rounding can leave the filter holding an invalid covariance. This
//! is square-root filtering in the sense of Potter (1963) and Bierman (1977), reduced to
//! the two-state case.
//!
//! The textbook update `P⁺ = (I − K H) P` subtracts nearly equal numbers whenever a
//! measurement is much more precise than the prediction it corrects. With no process
//! noise, `R = 1e-20` and a 9 ms step, the velocity variance should come out near
//! 1.2e-16. It came out as 250000 − 249999.99…, which rounds to −2.9e-11: one unit in the
//! last place of 250000, with the wrong sign. The position variance came out as exactly 0,
//! because `1 − K₀` rounded to zero once `R` fell below the spacing of `P₀₀`. That second
//! failure needs no exotic tuning: `R = 1e-8` and a 60 s gap produce it too. The factored
//! forms below contain no such subtraction.
//!
//! * **Predict.** The pre-array `[F L | L_Q]` is re-triangularised by Givens rotations.
//!   Rotations are orthogonal, so the result satisfies `L' L'ᵀ = F P Fᵀ + Q` exactly in
//!   exact arithmetic, and to backward-stable rounding otherwise.
//! * **Update.** For a position-only measurement the new factor is the old one with its
//!   first column scaled by `√(R/S)` and its second column unchanged. The derivation is
//!   in [`Kf1D::update`].

#[derive(Clone, Copy, Debug)]
pub struct Kf1D {
    pub x: f64,
    pub v: f64,
    /// Lower-triangular Cholesky factor of the covariance, `[[l00, 0], [l10, l11]]`, so
    /// that `P = [[l00², l00·l10], [l00·l10, l10² + l11²]]`. `l00` is kept non-negative.
    l00: f64,
    l10: f64,
    l11: f64,
    /// Continuous white-noise-acceleration PSD, m^2/s^3.
    q: f64,
    /// Measurement variance, m^2.
    r: f64,
}

/// Rotate row `k` into row `i` so that the first entry of row `k` becomes zero. A Givens
/// rotation: orthogonal, so it leaves `Σ rowᵀ row` unchanged.
fn rotate(rows: &mut [(f64, f64); 4], i: usize, k: usize) {
    let ((xi, yi), (xk, yk)) = (rows[i], rows[k]);
    if xk == 0.0 {
        return;
    }
    let r = xi.hypot(xk);
    let (c, s) = (xi / r, xk / r);
    rows[i] = (r, c * yi + s * yk);
    rows[k] = (0.0, c * yk - s * yi);
}

/// Result of a measurement update, kept for gating and track quality reporting.
#[derive(Clone, Copy, Debug)]
pub struct Innovation {
    /// Measurement minus prediction.
    pub y: f64,
    /// Innovation variance.
    pub s: f64,
}

impl Innovation {
    /// Normalised innovation squared. Should average ~1.0 on a consistent filter.
    pub fn nis(&self) -> f64 {
        self.y * self.y / self.s
    }

    pub fn sigma(&self) -> f64 {
        self.y.abs() / self.s.sqrt()
    }

    /// A finite residual over a finite, positive variance: the only kind a gate can
    /// judge. Anything else makes `sigma` NaN or infinite, and a NaN compares false
    /// against every threshold, so it would pass a gate asking "is it too large?".
    pub fn is_usable(&self) -> bool {
        self.y.is_finite() && self.s.is_finite() && self.s > 0.0
    }
}

impl Kf1D {
    /// Initialise on a first measurement: position is known to `r`, velocity is not
    /// known at all, so it gets a deliberately loose prior rather than a guess.
    pub fn new(x0: f64, v0: f64, q: f64, r: f64) -> Self {
        Self {
            x: x0,
            v: v0,
            l00: r.sqrt(),
            l10: 0.0,
            l11: 500.0,
            q,
            r,
        }
    }

    pub fn predict(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }
        self.x += self.v * dt;

        // P' = F P Fᵀ + Q with F = [[1, dt], [0, 1]], formed as M Mᵀ for the 2×4 pre-array
        // M = [F L | L_Q]. Below are the rows of Mᵀ, which are the columns of M:
        //   F L = [[l00 + dt·l10, dt·l11], [l10, l11]]
        //   L_Q = Cholesky factor of Q = q·[[dt³/3, dt²/2], [dt²/2, dt]]
        //       = √(q·dt) · [[dt/√3, 0], [√3/2, 1/2]], in closed form.
        let (a, b, c) = (self.l00, self.l10, self.l11);
        let g = (self.q * dt).sqrt();
        let root3 = 3f64.sqrt();
        let mut rows = [
            (a + dt * b, b),
            (dt * c, c),
            (g * dt / root3, g * root3 / 2.0),
            (0.0, g / 2.0),
        ];
        // Mᵀ = Q R by rotations: zero the first entry of rows 1–3 against row 0. What is
        // left in rows 1–3 is a single column, whose norm is the last entry of R. Then
        // M Mᵀ = Rᵀ R, so L' = Rᵀ.
        for k in 1..4 {
            rotate(&mut rows, 0, k);
        }
        let (l00, l10) = rows[0];
        // A factor's columns may change sign without changing L Lᵀ; keep l00 ≥ 0 so the
        // position sigma can be read off it directly.
        let sign = if l00 < 0.0 { -1.0 } else { 1.0 };
        self.l00 = sign * l00;
        self.l10 = sign * l10;
        self.l11 = rows[1].1.hypot(rows[2].1).hypot(rows[3].1);
    }

    /// Innovation for a candidate measurement, without applying it. Used for gating.
    pub fn innovation(&self, z: f64) -> Innovation {
        Innovation {
            y: z - self.x,
            s: self.l00 * self.l00 + self.r,
        }
    }

    /// Apply a position measurement.
    ///
    /// With `L = [[a, 0], [b, c]]`, `S = a² + R` and gain `K = [a², a·b] / S`,
    ///
    /// ```text
    /// P⁺ = P − K H P = [[a²·R/S,  a·b·R/S], [a·b·R/S,  c² + b²·R/S]]
    ///    = L⁺ L⁺ᵀ   with   L⁺ = [[a·f, 0], [b·f, c]],   f = √(R/S),
    /// ```
    ///
    /// because `S − a² = R` exactly. So the update is two multiplications by `f`, and
    /// `f` is formed as `√R / √S` so that it does not underflow when `R/S` would.
    pub fn update(&mut self, z: f64) -> Innovation {
        let innov = self.innovation(z);
        let (a, b) = (self.l00, self.l10);
        self.x += a * a / innov.s * innov.y;
        self.v += a * b / innov.s * innov.y;

        let f = self.r.sqrt() / innov.s.sqrt();
        self.l00 = a * f;
        self.l10 = b * f;
        innov
    }

    /// Every state and covariance-factor term is finite. A finite factor is a valid
    /// covariance: symmetric and positive semidefinite by construction.
    ///
    /// The filter itself does not validate what it is given: it is the arithmetic, and
    /// `TrackStore` is the boundary that decides what reaches it and what it keeps.
    pub fn is_finite(&self) -> bool {
        [self.x, self.v, self.l00, self.l10, self.l11]
            .iter()
            .all(|t| t.is_finite())
    }

    /// The covariance `[[P00, P01], [P01, P11]]`, formed from the factor. Symmetric by
    /// construction, and never stored.
    pub fn covariance(&self) -> [[f64; 2]; 2] {
        let (a, b, c) = (self.l00, self.l10, self.l11);
        let p01 = a * b;
        [[a * a, p01], [p01, b * b + c * c]]
    }

    /// One-sigma position uncertainty, metres. Read from the factor: no variance is
    /// ever negative, so nothing needs clamping before the square root.
    pub fn pos_sigma(&self) -> f64 {
        self.l00.abs()
    }

    /// One-sigma velocity uncertainty, metres per second.
    pub fn vel_sigma(&self) -> f64 {
        self.l10.hypot(self.l11)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converges_on_a_constant_velocity_target() {
        // Truth: starts at 0, moves at 200 m/s. Measurements every second, 30 m noise,
        // deterministic sawtooth so the test cannot flake.
        let mut kf = Kf1D::new(0.0, 0.0, 4.0, 900.0);
        let truth_v = 200.0;
        for step in 1..=60 {
            let t = step as f64;
            let noise = if step % 2 == 0 { 30.0 } else { -30.0 };
            kf.predict(1.0);
            kf.update(truth_v * t + noise);
        }
        assert!(
            (kf.v - truth_v).abs() < 10.0,
            "velocity estimate was {}",
            kf.v
        );
        assert!(
            kf.pos_sigma() < 30.0,
            "position sigma did not shrink: {}",
            kf.pos_sigma()
        );
        assert!(
            kf.vel_sigma() < 20.0,
            "velocity sigma did not shrink: {}",
            kf.vel_sigma()
        );
    }

    #[test]
    fn uncertainty_grows_when_coasting() {
        let mut kf = Kf1D::new(0.0, 100.0, 4.0, 900.0);
        kf.update(0.0);
        let before = kf.pos_sigma();
        kf.predict(30.0);
        assert!(
            kf.pos_sigma() > before,
            "coasting must widen the error ellipse"
        );
    }

    /// Covariance validity, checked from what the filter reports. Both one-sigma values
    /// must be positive: a zero sigma on a filter that has only ever seen noisy
    /// measurements is a variance that went negative and was clamped, or one that was
    /// rounded away.
    fn assert_valid(kf: &Kf1D, what: &str) {
        let (ps, vs) = (kf.pos_sigma(), kf.vel_sigma());
        assert!(
            ps.is_finite() && ps > 0.0 && vs.is_finite() && vs > 0.0,
            "{what}: pos sigma {ps}, vel sigma {vs}, {kf:?}"
        );
    }

    #[test]
    fn a_precise_measurement_cannot_make_a_variance_negative() {
        // Audit regression, the exact case: no process noise, R = 1e-20, two identical
        // stationary plots 9 ms apart. The subtractive update computed the velocity
        // variance as 250000 - 249999.99..., rounded to -2.9e-11, and the position
        // variance as exactly 0.
        let mut kf = Kf1D::new(0.0, 0.0, 0.0, 1e-20);
        kf.predict(0.009);
        kf.update(0.0);
        assert_valid(&kf, "after the 9 ms update");
        for step in 1..=5 {
            kf.predict(1.0);
            let innov = kf.innovation(0.0);
            assert!(innov.is_usable(), "step {step}: {innov:?}");
            kf.update(0.0);
            assert_valid(&kf, &format!("stationary update {step}"));
        }
    }

    #[test]
    fn tiny_measurement_variances_keep_a_valid_covariance() {
        // Neighbouring scales, each step sizes from microseconds to minutes, repeated
        // stationary updates, and process noise both off and on.
        for r in [1e-8, 1e-12, 1e-16, 1e-20, 1e-40, 1e-100, 1e-300, 5e-324] {
            for q in [0.0, 1e-20, 4.0] {
                for dt in [1e-6, 1e-3, 0.009, 1.0, 60.0] {
                    let mut kf = Kf1D::new(0.0, 0.0, q, r);
                    for step in 0..200 {
                        kf.predict(dt);
                        kf.update(0.0);
                        assert_valid(&kf, &format!("r {r} q {q} dt {dt} step {step}"));
                    }
                }
            }
        }
    }

    #[test]
    fn covariance_stays_symmetric_and_positive_definite_over_long_runs() {
        // 60 generated runs of 2,000 steps: step sizes from a microsecond to two minutes,
        // noisy measurements, and tunings across the admitted domain, from the smallest
        // normal double to 1e9. After every step, on the factor the filter actually holds:
        //   * every term finite, l00 > 0 and l11 > 0, so P = L Lᵀ has determinant
        //     (l00·l11)² > 0: positive definite, and neither sigma can read zero.
        // And on the covariance formed from it, wherever that is representable:
        //   * symmetric: the off-diagonal terms are one stored number, so exactly equal;
        //   * P01² ≤ P00·P11, to within 5ε relative. In exact arithmetic
        //     P01² = a²b² ≤ a²(b² + c²) = P00·P11. With u = ε/2, the computed left side
        //     is at most (1 + u)³ times its exact value: three roundings, in a·b and in
        //     squaring it. The computed right side is at least (1 − u)⁵ times its exact
        //     value before the allowance: one rounding in a·a, two in b·b + c·c (both
        //     terms positive), one in the product and one in applying the allowance. So
        //     the check can fail spuriously only if (1 + τ) < (1 + u)³ / (1 − u)⁵, which
        //     is 1 + 4ε plus second-order terms. τ = 5ε clears that with a margin of ε;
        //     it is the bound, not a tuning knob. The
        //     bound assumes normal numbers, so it is checked where P00·P11 is one. At the
        //     smallest variances the covariance entries themselves fall below the smallest
        //     double while the factor, which holds their square roots, does not.
        let mut seed = 0x5ca1_ab1e_u64;
        let mut unit = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let qs = [0.0, 1e-20, 4.0, 1e9];
        let rs = [f64::MIN_POSITIVE, 1e-20, 1e-8, 900.0, 1e9];
        let (mut steps, mut matrix_checked) = (0, 0);
        for run in 0..60 {
            let (q, r) = (qs[run % qs.len()], rs[run % rs.len()]);
            let mut kf = Kf1D::new(0.0, 0.0, q, r);
            let mut truth = 0.0;
            for step in 0..2_000 {
                let dt = 10f64.powf(unit() * 8.0 - 6.0).min(120.0);
                truth += 200.0 * dt;
                kf.predict(dt);
                kf.update(truth + (unit() - 0.5) * 2.0 * r.sqrt().min(1e4));
                let [[p00, p01], [p10, p11]] = kf.covariance();
                let at = format!(
                    "run {run} (q {q:e}, r {r:e}) step {step}: l00 {:e} l10 {:e} l11 {:e} \
                     P00 {p00:e} P01 {p01:e} P11 {p11:e}",
                    kf.l00, kf.l10, kf.l11
                );
                assert!(kf.is_finite() && kf.l00 > 0.0 && kf.l11 > 0.0, "{at}");
                assert_eq!(p01.to_bits(), p10.to_bits(), "{at}");
                steps += 1;
                if p00 * p11 >= f64::MIN_POSITIVE {
                    matrix_checked += 1;
                    assert!(p00 > 0.0 && p11 > 0.0, "{at}");
                    assert!(p01 * p01 <= p00 * p11 * (1.0 + 5.0 * f64::EPSILON), "{at}");
                }
            }
        }
        // The representable case must be the common one, or the matrix check is vacuous.
        assert!(
            matrix_checked * 10 > steps * 7,
            "{matrix_checked} of {steps}"
        );
    }

    #[test]
    fn the_factored_filter_matches_the_textbook_filter_where_that_is_accurate() {
        // Validity alone does not pin the arithmetic: a factor that dropped part of Q
        // would still be valid. At ordinary tunings the textbook forms lose nothing, so
        // the factored filter must reproduce them, state and covariance, step by step.
        // The two round differently and agree only to rounding. 1e-9 relative leaves
        // seven orders of magnitude above double precision and is still far below what
        // a change of formula does to the result.
        for (q, r) in [(4.0, 900.0), (0.5, 25.0), (50.0, 1e4), (4.0, 1.0)] {
            let mut kf = Kf1D::new(10.0, 200.0, q, r);
            let (mut x, mut v) = (10.0, 200.0);
            let (mut p00, mut p01, mut p11) = (r, 0.0, 250_000.0);
            for step in 1..=50 {
                let dt = [0.5, 1.0, 2.0, 7.0][step % 4];
                x += v * dt;
                (p00, p01, p11) = (
                    p00 + 2.0 * dt * p01 + dt * dt * p11 + q * dt.powi(3) / 3.0,
                    p01 + dt * p11 + q * dt * dt / 2.0,
                    p11 + q * dt,
                );
                kf.predict(dt);

                let z = 260.0 * step as f64 + if step % 2 == 0 { 17.0 } else { -23.0 };
                let s = p00 + r;
                let (k0, k1) = (p00 / s, p01 / s);
                let y = z - x;
                (x, v) = (x + k0 * y, v + k1 * y);
                (p00, p01, p11) = ((1.0 - k0) * p00, (1.0 - k0) * p01, p11 - k1 * p01);
                kf.update(z);

                let [[g00, g01], [_, g11]] = kf.covariance();
                let scale = (p00 * p11).sqrt();
                let close = |got: f64, want: f64, of: f64| (got - want).abs() <= 1e-9 * of;
                assert!(
                    close(kf.x, x, x.abs())
                        && close(kf.v, v, v.abs())
                        && close(g00, p00, p00)
                        && close(g01, p01, scale)
                        && close(g11, p11, p11),
                    "q {q} r {r} step {step}: factored x {} v {} P [{g00:e} {g01:e} {g11:e}], \
                     textbook x {x} v {v} P [{p00:e} {p01:e} {p11:e}]",
                    kf.x,
                    kf.v
                );
            }
        }
    }

    #[test]
    fn gating_flags_an_implausible_jump() {
        let mut kf = Kf1D::new(0.0, 0.0, 4.0, 900.0);
        for _ in 0..20 {
            kf.predict(1.0);
            kf.update(0.0);
        }
        // A 50 km jump on a stationary track is a decoded-garbage report, not a manoeuvre.
        assert!(kf.innovation(50_000.0).sigma() > 5.0);
        assert!(kf.innovation(40.0).sigma() < 5.0);
    }
}
