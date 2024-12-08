use std::f64;

use xilem_web::svg::kurbo;

use kurbo::{common::GAUSS_LEGENDRE_COEFFS_32, Affine, Point, Vec2};
use nalgebra::{ArrayStorage, Matrix5, Vector5};

use crate::utils::*;

use super::*;

impl<D: DualNum<f64> + Copy> HyperbezParams<D> {
    pub fn kappa_extrema(&self) -> ArrayVec<D, 2> {
        let quad_a = self.a * self.c * 4.;
        let quad_b = self.a * self.d + self.b * self.c * 6.;
        let quad_c = self.b * self.d * 3. - self.a * 2.;
        let det = (quad_b.powi(2) - quad_a * quad_c * 4.).sqrt();
        if !det.re().is_finite() {
            return ArrayVec::new();
        }
        let p = (-quad_b + det) / quad_a * 0.5;
        let m = (-quad_b - det) / quad_a * 0.5;
        [(m.re() < 0.).then_some(m), (p.re() > 1.).then_some(p)]
            .into_iter()
            .flatten()
            .collect()
    }
}

const PI_2: f64 = f64::consts::FRAC_PI_2;

impl HyperbezParams<f64> {
    pub fn from_control(p1: Point, p2: Point) -> Self {
        // let pts = this.cubic.pts;
        // let chord = pts[3].minus(pts[0]).hypot();
        let chord = 1.;
        // let dx0 = pts[1].x - pts[0].x;
        // let dy0 = pts[1].y - pts[0].y;
        // let dx1 = pts[2].x - pts[3].x;
        // let dy1 = pts[2].y - pts[3].y;
        let dp0 = p1.to_vec2();
        let dp1 = p2 - Point::new(1., 0.);
        // let th0 = Math.atan2(-dy0, dx0);
        // let th1 = Math.atan2(-dy1, -dx1);
        let th0 = -(dp0.atan2());
        let th1 = (-dp1).atan2();
        let tens0 = dp0.hypot() / chord * 1.5 * (th0.cos() + 1.);
        let tens1 = dp1.hypot() / chord * 1.5 * (th1.cos() + 1.);
        let d0 = dp0.hypot() / chord;
        let d1 = dp1.hypot() / chord;
        let mut k0 = k_for_tension(tens0);
        let mut k1 = k_for_tension(tens1);
        let cbr = (k0 / k1).powf(1. / 3.);
        fn soft(x: f64) -> f64 {
            (0.5 * (1. + x * x)).sqrt()
        }
        k1 /= soft(cbr);
        k0 /= soft(1. / cbr);

        let dc = (p2 - p1).hypot() / chord;
        let kmid = dc.powf(1.5);
        let ratio = (d0 / d1).powf(1.5);
        let blend = 0.5 + 0.5 * (3. - 10. * dc).tanh();
        k0 += blend * (kmid / ratio - k0);
        k1 += blend * (kmid * ratio - k1);
        let [c, d] = quadratic_for_endk(k0, k1);
        //console.log('dc', dc, 'c', cd.c, 'd', cd.d);
        // let endk = endk_for_quadratic(c, d);
        // console.log(dc, k0, k1, blend);
        let [a, b] = solve_thetas(th0, th1, c, d, 1.);
        HyperbezParams::new(a, b, c, d, 1.)
    }

    fn integrate_any<R: std::ops::Add<Output = R> + std::ops::Mul<f64, Output = R>>(
        &self,
        f: impl Fn(f64) -> R,
        init: R,
        t: f64,
    ) -> R {
        // TODO: improve accuracy by subdividing in near-cusp cases
        let mut xy = init;
        let u0 = 0.5 * t;
        for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
            let u = u0 + u0 * xi;
            xy = xy + f(u) * *wi;
        }
        xy * u0
    }

    fn denom(&self) -> f64 {
        2. / (4. * self.c - self.d.powi(2))
    }

    pub fn dtheta_da(&self, t: f64) -> f64 {
        let beta1 = -2. * self.denom();
        -beta1 + (beta1 - 2. * self.d * t) / self.q(t).sqrt()
    }

    pub fn dtheta_db(&self, t: f64) -> f64 {
        let beta0 = self.d * self.denom();
        -beta0 + (beta0 - 4. * self.c * t) / self.q(t).sqrt()
    }

    pub fn dtheta_dc(&self, t: f64) -> f64 {
        let denom = self.denom();
        let term1 = 2. * self.num0 * denom;
        let term2 = -self.int_helper(t) * t.powi(2) / (2. * self.q(t));
        let term3 = 2. * denom * (-self.int_helper(t) + self.b * t / self.q(t).sqrt());
        term1 + term2 + term3
    }

    pub fn dtheta_dd(&self, t: f64) -> f64 {
        let denom = self.denom();
        let term1 = -self.d * self.num0 * denom;
        let term2 = -self.b * denom;
        let term3 = -self.int_helper(t) * t / (2. * self.q(t));
        let term4 =
            denom * (self.d * self.int_helper(t) + (self.b - self.a * t) / self.q(t).sqrt());
        term1 + term2 + term3 + term4
    }

    pub fn dx_da(&self, t: f64) -> f64 {
        -self.theta(t).sin() * self.dtheta_da(t)
    }

    pub fn dx_db(&self, t: f64) -> f64 {
        -self.theta(t).sin() * self.dtheta_db(t)
    }

    pub fn dx_dc(&self, t: f64) -> f64 {
        -self.theta(t).sin() * self.dtheta_dc(t)
    }

    pub fn dx_dd(&self, t: f64) -> f64 {
        -self.theta(t).sin() * self.dtheta_dd(t)
    }

    pub fn dy_da(&self, t: f64) -> f64 {
        self.theta(t).cos() * self.dtheta_da(t)
    }

    pub fn dy_db(&self, t: f64) -> f64 {
        self.theta(t).cos() * self.dtheta_db(t)
    }

    pub fn dy_dc(&self, t: f64) -> f64 {
        self.theta(t).cos() * self.dtheta_dc(t)
    }

    pub fn dy_dd(&self, t: f64) -> f64 {
        self.theta(t).cos() * self.dtheta_dd(t)
    }

    // pub fn dp_da(&self, t: f64) -> Vec2 {
    //     self.integrate_any(|t| Vec2::new(self.dx_da(t), self.dy_da(t)), Vec2::ZERO, t)
    // }

    // pub fn dp_db(&self, t: f64) -> Vec2 {
    //     self.integrate_any(|t| Vec2::new(self.dx_db(t), self.dy_db(t)), Vec2::ZERO, t)
    // }

    // pub fn dp_dc(&self, t: f64) -> Vec2 {
    //     self.integrate_any(|t| Vec2::new(self.dx_dc(t), self.dy_da(t)), Vec2::ZERO, t)
    // }

    // pub fn dp_dd(&self, t: f64) -> Vec2 {
    //     self.integrate_any(|t| Vec2::new(self.dx_dd(t), self.dy_da(t)), Vec2::ZERO, t)
    // }

    pub fn dp_da(&self, t: f64) -> Vec2 {
        self.integrate_any(
            |t| self.dtheta_da(t) * Vec2::from_angle(self.theta(t) + PI_2),
            Vec2::ZERO,
            t,
        )
    }

    pub fn dp_db(&self, t: f64) -> Vec2 {
        self.integrate_any(
            |t| self.dtheta_db(t) * Vec2::from_angle(self.theta(t) + PI_2),
            Vec2::ZERO,
            t,
        )
    }

    pub fn dp_dc(&self, t: f64) -> Vec2 {
        self.integrate_any(
            |t| self.dtheta_dc(t) * Vec2::from_angle(self.theta(t) + PI_2),
            Vec2::ZERO,
            t,
        )
    }

    pub fn dp_dd(&self, t: f64) -> Vec2 {
        self.integrate_any(
            |t| self.dtheta_dd(t) * Vec2::from_angle(self.theta(t) + PI_2),
            Vec2::ZERO,
            t,
        )
    }

    fn system_for_solving(
        p0_5_i: kurbo::Point,
        phi0_5_i: f64,
        theta1_i: f64,
        p1_angle_i: f64,
    ) -> impl Fn([f64; 5]) -> ([f64; 5], [[f64; 5]; 5]) {
        move |[a, b, c, d, t]: [f64; 5]| -> ([f64; 5], [[f64; 5]; 5]) {
            let hb = HyperbezParams::new(a, b, c, d, 1.);

            let p1 = as_vec2(hb.integrate(1.));
            let p1_hypot = p1.hypot();
            let p0_5 = hb.integrate(t) / p1_hypot;
            let phi0_5 = hb.theta(t);
            let theta1 = hb.theta(1.);
            let p1_angle = p1.atan2();

            let dp1_da = hb.dp_da(1.);
            let dp1_db = hb.dp_db(1.);
            let dp1_dc = hb.dp_dc(1.);
            let dp1_dd = hb.dp_dd(1.);

            let p1_norm = p1.normalize();
            let dp1_hypot_da = p1_norm.dot(dp1_da);
            let dp1_hypot_db = p1_norm.dot(dp1_db);
            let dp1_hypot_dc = p1_norm.dot(dp1_dc);
            let dp1_hypot_dd = p1_norm.dot(dp1_dd);

            let dp0_5_da = hb.dp_da(t);
            let dp0_5_db = hb.dp_db(t);
            let dp0_5_dc = hb.dp_dc(t);
            let dp0_5_dd = hb.dp_dd(t);
            let dp0_5_dt = Vec2::from_angle(hb.theta(t)) / p1_hypot;

            let p1_hypot2 = p1.hypot2();
            let dp0_5_da_x =
                Vec2::new(p1_hypot, p0_5.x).cross(Vec2::new(dp1_hypot_da, dp0_5_da.x)) / p1_hypot2;
            let dp0_5_db_x =
                Vec2::new(p1_hypot, p0_5.x).cross(Vec2::new(dp1_hypot_db, dp0_5_db.x)) / p1_hypot2;
            let dp0_5_dc_x =
                Vec2::new(p1_hypot, p0_5.x).cross(Vec2::new(dp1_hypot_dc, dp0_5_dc.x)) / p1_hypot2;
            let dp0_5_dd_x =
                Vec2::new(p1_hypot, p0_5.x).cross(Vec2::new(dp1_hypot_dd, dp0_5_dd.x)) / p1_hypot2;

            let dp0_5_da_y =
                Vec2::new(p1_hypot, p0_5.y).cross(Vec2::new(dp1_hypot_da, dp0_5_da.y)) / p1_hypot2;
            let dp0_5_db_y =
                Vec2::new(p1_hypot, p0_5.y).cross(Vec2::new(dp1_hypot_db, dp0_5_db.y)) / p1_hypot2;
            let dp0_5_dc_y =
                Vec2::new(p1_hypot, p0_5.y).cross(Vec2::new(dp1_hypot_dc, dp0_5_dc.y)) / p1_hypot2;
            let dp0_5_dd_y =
                Vec2::new(p1_hypot, p0_5.y).cross(Vec2::new(dp1_hypot_dd, dp0_5_dd.y)) / p1_hypot2;

            let dphi0_5_da = hb.dtheta_da(t);
            let dphi0_5_db = hb.dtheta_db(t);
            let dphi0_5_dc = hb.dtheta_dc(t);
            let dphi0_5_dd = hb.dtheta_dd(t);
            let dphi0_5_dt = hb.kappa(t);

            let dtheta1_da = hb.dtheta_da(1.);
            let dtheta1_db = hb.dtheta_db(1.);
            let dtheta1_dc = hb.dtheta_dc(1.);
            let dtheta1_dd = hb.dtheta_dd(1.);
            let dtheta1_dt = 0.;

            let dp1_angle_da = p1.cross(dp1_da) / p1_hypot2;
            let dp1_angle_db = p1.cross(dp1_db) / p1_hypot2;
            let dp1_angle_dc = p1.cross(dp1_dc) / p1_hypot2;
            let dp1_angle_dd = p1.cross(dp1_dd) / p1_hypot2;
            let dp1_angle_dt = 0.;

            let err = [
                p0_5.x - p0_5_i.x,
                p0_5.y - p0_5_i.y,
                norm_radians(phi0_5 - phi0_5_i),
                norm_radians(theta1 - theta1_i),
                norm_radians(p1_angle - p1_angle_i),
            ];

            let jac = [
                [dp0_5_da_x, dp0_5_db_x, dp0_5_dc_x, dp0_5_dd_x, dp0_5_dt.x],
                [dp0_5_da_y, dp0_5_db_y, dp0_5_dc_y, dp0_5_dd_y, dp0_5_dt.y],
                [dphi0_5_da, dphi0_5_db, dphi0_5_dc, dphi0_5_dd, dphi0_5_dt],
                [dtheta1_da, dtheta1_db, dtheta1_dc, dtheta1_dd, dtheta1_dt],
                [
                    dp1_angle_da,
                    dp1_angle_db,
                    dp1_angle_dc,
                    dp1_angle_dd,
                    dp1_angle_dt,
                ],
            ];

            (err, jac)
        }
    }
}

#[must_use]
fn solve_iterate_once(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta1: f64,
    p1_angle: f64,
    guess: [f64; 5],
) -> ([f64; 5], Option<[f64; 5]>) {
    let (err, jac) = HyperbezParams::system_for_solving(p0_5, phi0_5, theta1, p1_angle)(guess);

    let guess_v = Vector5::from_data(ArrayStorage([guess]));
    let err_v = Vector5::from_data(ArrayStorage([err]));
    let mut jac_v = Matrix5::from_data(ArrayStorage(jac)).transpose();
    let new_guess_v = (err_v.norm_squared().is_finite() && jac_v.try_inverse_mut())
        .then(|| guess_v - jac_v * err_v);

    (err, new_guess_v.map(|v| v.data.0[0]))
}

pub fn solve_for_params_exact(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta0: f64,
    theta1: f64,
    mut guess: [f64; 5],
    threshold: f64,
    n_iter: usize,
) -> SolveResult<5> {
    if radian_in_line(theta0) && radian_in_line(theta1) {
        return Ok(Solution {
            params: Vector5::new(0., 0., 1., -1., 0.5),
            err: Vector5::zeros(),
            iter: 0,
        });
    }

    let p1_angle = -theta0;
    let p0_5 = Affine::rotate(p1_angle) * p0_5;
    let phi0_5 = phi0_5 + p1_angle;
    let theta1 = theta1 + p1_angle;

    let mut err = [f64::INFINITY; 5];
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once(p0_5, phi0_5, theta1, p1_angle, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity {
                guess: Vector5::from_data(ArrayStorage([guess])),
                err: Vector5::from_data(ArrayStorage([new_err])),
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() <= threshold) {
            return Ok(Solution {
                params: Vector5::from_data(ArrayStorage([new_guess])),
                err: Vector5::from_data(ArrayStorage([new_err])),
                iter: i,
            });
        }
        guess = new_guess;
        err = new_err;
    }
    Err(SolveError::OutOfIteration {
        guess: Vector5::from_data(ArrayStorage([guess])),
        err: Vector5::from_data(ArrayStorage([err])),
    })
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn conrtol_test_sym1() {
        dbg!(HyperbezParams::from_control(
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1() {
        solve_helper(
            kurbo::Point::new(0.1, 0.4),
            kurbo::Point::new(0.3, 0.4),
            solve_with_guess(solve_for_params_exact, [-1., -1., -1., 1.]),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        solve_helper(
            kurbo::Point::new(0.1, 0.3),
            kurbo::Point::new(0.3, 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test3() {
        solve_helper(
            kurbo::Point::new(0., 0.3),
            kurbo::Point::new(1., 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test4() {
        solve_helper(
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [0., -1., -1., 1.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test5() {
        solve_helper(
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }
}
