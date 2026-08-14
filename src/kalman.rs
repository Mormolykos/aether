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

#[derive(Clone, Copy, Debug)]
pub struct Kf1D {
    pub x: f64,
    pub v: f64,
    /// Covariance [[p00, p01], [p10, p11]], symmetric.
    p00: f64,
    p01: f64,
    p10: f64,
    p11: f64,
    /// Continuous white-noise-acceleration PSD, m^2/s^3.
    q: f64,
    /// Measurement variance, m^2.
    r: f64,
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
}

impl Kf1D {
    /// Initialise on a first measurement: position is known to `r`, velocity is not
    /// known at all, so it gets a deliberately loose prior rather than a guess.
    pub fn new(x0: f64, v0: f64, q: f64, r: f64) -> Self {
        Self {
            x: x0,
            v: v0,
            p00: r,
            p01: 0.0,
            p10: 0.0,
            p11: 500.0 * 500.0,
            q,
            r,
        }
    }

    pub fn predict(&mut self, dt: f64) {
        if dt <= 0.0 {
            return;
        }
        self.x += self.v * dt;

        // P = F P F^T + Q, with F = [[1, dt], [0, 1]].
        let p00 = self.p00 + dt * (self.p01 + self.p10) + dt * dt * self.p11;
        let p01 = self.p01 + dt * self.p11;
        let p10 = self.p10 + dt * self.p11;
        let p11 = self.p11;

        // Continuous white-noise acceleration discretised over dt.
        let (dt2, dt3) = (dt * dt, dt * dt * dt);
        self.p00 = p00 + self.q * dt3 / 3.0;
        self.p01 = p01 + self.q * dt2 / 2.0;
        self.p10 = p10 + self.q * dt2 / 2.0;
        self.p11 = p11 + self.q * dt;
    }

    /// Innovation for a candidate measurement, without applying it. Used for gating.
    pub fn innovation(&self, z: f64) -> Innovation {
        Innovation {
            y: z - self.x,
            s: self.p00 + self.r,
        }
    }

    pub fn update(&mut self, z: f64) -> Innovation {
        let innov = self.innovation(z);
        let k0 = self.p00 / innov.s;
        let k1 = self.p10 / innov.s;

        self.x += k0 * innov.y;
        self.v += k1 * innov.y;

        // P = (I - K H) P, H = [1, 0]. Old values on the right-hand side.
        let (p00, p01, p10, p11) = (self.p00, self.p01, self.p10, self.p11);
        self.p00 = (1.0 - k0) * p00;
        self.p01 = (1.0 - k0) * p01;
        self.p10 = p10 - k1 * p00;
        self.p11 = p11 - k1 * p01;

        innov
    }

    /// One-sigma position uncertainty, metres.
    pub fn pos_sigma(&self) -> f64 {
        self.p00.max(0.0).sqrt()
    }

    /// One-sigma velocity uncertainty, metres per second.
    pub fn vel_sigma(&self) -> f64 {
        self.p11.max(0.0).sqrt()
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
