use std::f64;

use approx::AbsDiffEq as _;
use nalgebra::{Vector2, Vector5};
use num_dual::jacobian;
use xilem_web::svg::kurbo::{self, Affine};

use kurbo::{common::GAUSS_LEGENDRE_COEFFS_32, ParamCurve, ParamCurveDeriv};

use crate::num_dual_ext::*;

#[derive(Clone, Copy, Debug)]
pub struct HyperbezParams<D> {
    a: D,
    b: D,
    c: D,
    d: D,
    e: D,

    num0: D,
    num1: D,
}

impl<D: DualNumExt<f64>> HyperbezParams<D> {
    /// Create a new hyperbezier with the given parameters.
    pub fn new(a: D, b: D, c: D, d: D, e: D) -> Self {
        let denom = D::from(2.) / (c * 4. - d * d);
        let beta0 = d * denom;
        let num0 = a * (-d * beta0 * 0.5 - 1.) / c + b * beta0;
        let num1 = (b * c * 2. - d * a) * denom;

        HyperbezParams {
            a,
            b,
            c,
            d,
            e,
            num0,
            num1,
        }
    }

    fn int_helper(&self, t: D) -> D {
        // assumes self.e = 1
        let q = self.c * t * t + self.d * t + 1.;
        (self.num0 + self.num1 * t) / q.sqrt()
    }

    /// Determine the angle for the given parameter.
    ///
    /// This can be interpreted as a Whewell representation of the
    /// curve. The `t` parameter ranges from 0 to 1, and the returned
    /// value is 0 for `t = 0`.
    pub fn theta(&self, t: D) -> D {
        self.int_helper(t) - self.num0
    }

    /// Returns [q, κ]
    fn eval_q(&self, t: D) -> [D; 2] {
        let q = self.c * t * t + self.d * t + self.e;
        let k = (self.a * t + self.b) / (q * q.sqrt());
        [q, k]
    }

    pub fn kappa(&self, t: D) -> D {
        self.eval_q(t)[1]
    }

    pub fn q(&self, t: D) -> D {
        self.eval_q(t)[0]
    }

    fn report_endpoints(&self) -> [D; 2] {
        let mut sum = D::from(0.);
        for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
            // for (let i = 0; i < co.length; i += 2) {
            // let xi = co[i + 1];
            // let wi = co[i];
            let t = 0.5 + 0.5 * xi;
            let q = self.q(D::from(t));
            sum += q.powf(-1.5) * *wi;
        }
        let integral = sum * 0.5;
        let q0 = self.q(D::from(0.));
        let q1 = self.q(D::from(1.));
        [q0.powf(-1.5) / integral, q1.powf(-1.5) / integral]
    }

    /// Evaluate the position of the raw curve.
    ///
    /// This is simply the integral of the Whewell representation,
    /// so that the total arc length is unit, and the initial tangent
    /// is horizontal.
    pub fn integrate(&self, t: D) -> Vector2<D> {
        // TODO: improve accuracy by subdividing in near-cusp cases
        let mut xy = Vector2::new(D::from(0.), D::from(0.));
        let u0 = t * 0.5;
        for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
            let u = u0 + u0 * *xi;
            let (y, x) = self.theta(u).sin_cos();
            xy += Vector2::new(x, y) * D::from(*wi);
        }
        xy * u0
    }

    pub fn make(a: D, b: D, c0: D, c1: D) -> Self {
        let c = -(c0 + c1) * (c0 + c1);
        let d = c0 * (c0 + c1) * 2.;
        let e = -c0 * c0 + 1.;
        HyperbezParams::new(a, b, c, d, e)
    }

    // pub fn from_control(p1: Point, p2: Point) -> Self {
    //     // let pts = this.cubic.pts;
    //     // let chord = pts[3].minus(pts[0]).hypot();
    //     let chord = 1.;
    //     // let dx0 = pts[1].x - pts[0].x;
    //     // let dy0 = pts[1].y - pts[0].y;
    //     // let dx1 = pts[2].x - pts[3].x;
    //     // let dy1 = pts[2].y - pts[3].y;
    //     let dp0 = p1.to_vec2();
    //     let dp1 = p2 - Point::new(1., 0.);
    //     // let th0 = Math.atan2(-dy0, dx0);
    //     // let th1 = Math.atan2(-dy1, -dx1);
    //     let th0 = -(dp0.atan2());
    //     let th1 = (-dp1).atan2();
    //     let tens0 = dp0.hypot() / chord * 1.5 * (th0.cos() + 1.);
    //     let tens1 = dp1.hypot() / chord * 1.5 * (th1.cos() + 1.);
    //     let d0 = dp0.hypot() / chord;
    //     let d1 = dp1.hypot() / chord;
    //     let mut k0 = k_for_tension(tens0);
    //     let mut k1 = k_for_tension(tens1);
    //     let cbr = (k0 / k1).powf(1. / 3.);
    //     fn soft(x: f64) -> f64 {
    //         (0.5 * (1. + x * x)).sqrt()
    //     }
    //     k1 /= soft(cbr);
    //     k0 /= soft(1. / cbr);

    //     let dc = (p2 - p1).hypot() / chord;
    //     let kmid = dc.powf(1.5);
    //     let ratio = (d0 / d1).powf(1.5);
    //     let blend = 0.5 + 0.5 * (3. - 10. * dc).tanh();
    //     k0 += blend * (kmid / ratio - k0);
    //     k1 += blend * (kmid * ratio - k1);
    //     let [c, d] = quadratic_for_endk(k0, k1);
    //     //console.log('dc', dc, 'c', cd.c, 'd', cd.d);
    //     let endk = endk_for_quadratic(c, d);
    //     // console.log(dc, k0, k1, blend);
    //     let [a, b] = solve_thetas(th0, th1, c, d, 1.);
    //     HyperbezParams::new(a, b, c, d, 1.)
    // }

    /// Returns [θ0, θ1]
    pub fn calc_thetas(&self) -> [D; 2] {
        let p = self.integrate(D::from(1.));
        let th0 = p.y.atan2(p.x);
        let th1 = self.theta(D::from(1.)) - th0;
        [th0, th1]
    }

    pub fn a(&self) -> D {
        self.a
    }
    pub fn b(&self) -> D {
        self.b
    }
    pub fn c(&self) -> D {
        self.c
    }
    pub fn d(&self) -> D {
        self.d
    }
    pub fn e(&self) -> D {
        self.e
    }
    pub fn num1(&self) -> D {
        self.num1
    }
    pub fn num0(&self) -> D {
        self.num0
    }
}

fn system_for_solving<D: DualNumExt<f64>>(
    p0_5_i: kurbo::Point,
    phi0_5_i: f64,
    theta1_i: f64,
    p1_angle_i: f64,
) -> impl Fn(Vector5<D>) -> Vector5<D> {
    move |guess: Vector5<D>| -> Vector5<D> {
        let a = guess.x;
        let b = guess.y;
        let c = guess.z;
        let d = guess.w;
        let t = guess.a;

        let hb = HyperbezParams::new(a, b, c, d, D::from(1.));

        let p1 = hb.integrate(D::from(1.));
        let p1_hypot = (p1.y.powi(2) + p1.x.powi(2)).sqrt();
        let p0_5 = hb.integrate(t) / p1_hypot;
        let phi0_5 = hb.theta(t);
        let theta1 = hb.theta(D::from(1.));
        let p1_angle = p1.y.atan2(p1.x);

        // dbg!(
        Vector5::new(
            p0_5.x - p0_5_i.x,
            p0_5.y - p0_5_i.y,
            norm_radians(phi0_5 - phi0_5_i),
            norm_radians(theta1 - theta1_i),
            norm_radians(p1_angle - p1_angle_i),
        )
        // )
    }
}

#[must_use]
fn solve_iterate_once(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta1: f64,
    p1_angle: f64,
    guess: Vector5<f64>,
) -> (Vector5<f64>, Option<Vector5<f64>>) {
    let (f, mut jac) = jacobian(system_for_solving(p0_5, phi0_5, theta1, p1_angle), guess);
    // jac.transpose_mut();
    // dbg!(&f, &jac);
    let new_guess =
        (f.norm_squared().is_finite() && jac.try_inverse_mut()).then(|| guess + jac * f);
    (f, new_guess)
}

#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub params: Vector5<f64>,
    pub err: Vector5<f64>,
    pub iter: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum SolveError {
    Singularity {
        guess: Vector5<f64>,
        err: Vector5<f64>,
        iter: usize,
    },
    OutOfIteration {
        guess: Vector5<f64>,
        err: Vector5<f64>,
    },
}

pub fn solve_for_params_exact(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta0: f64,
    theta1: f64,
    guess: [f64; 5],
    threshold: f64,
    n_iter: usize,
) -> Result<Solution, SolveError> {
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

    let mut guess = Vector5::from_data(nalgebra::ArrayStorage([guess]));
    let mut err = Vector5::from_data(nalgebra::ArrayStorage([[f64::INFINITY; 5]]));
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once(p0_5, phi0_5, theta1, p1_angle, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity {
                guess,
                err: new_err,
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() < threshold) {
            return Ok(Solution {
                params: new_guess,
                err: new_err,
                iter: i,
            });
        }
        guess = new_guess;
        err = new_err;
    }
    Err(SolveError::OutOfIteration { guess, err })
}

const EPSILON: f64 = 0.01;
pub fn radian_in_line(theta: f64) -> bool {
    (theta % f64::consts::PI).abs_diff_eq(&0., EPSILON)
}

pub fn solve_for_cubic(
    cubicbez: kurbo::CubicBez,
    guess: [f64; 4],
    threshold: f64,
    n_iter: usize,
) -> Result<Solution, SolveError> {
    let p0_5 = cubicbez.eval(0.5);
    let res = solve_for_params_exact(
        p0_5,
        cubicbez.deriv().eval(0.5).to_vec2().atan2(),
        cubicbez.p1.to_vec2().atan2(),
        (cubicbez.p3 - cubicbez.p2).atan2(),
        [guess[0], guess[1], guess[2], guess[3], p0_5.x],
        threshold,
        n_iter,
    );
    dbg!(p0_5);
    dbg!(res)
}

mod tests {
    use super::*;

    #[test]
    #[allow(unused_must_use)]
    fn test1() {
        let cubicbez = kurbo::CubicBez::new(
            kurbo::Point::ZERO,
            kurbo::Point::new(0.1, 0.4),
            kurbo::Point::new(0.3, 0.4),
            kurbo::Point::new(1., 0.),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
        solve_for_cubic(cubicbez, [-1., -1., 1., 0.], 1e-2, 11);
    }

    #[test]
    #[allow(unused_must_use)]
    fn test2() {
        let cubicbez = kurbo::CubicBez::new(
            kurbo::Point::ZERO,
            kurbo::Point::new(0.1, 0.3),
            kurbo::Point::new(0.3, 0.3),
            kurbo::Point::new(1., 0.),
        );
        solve_for_cubic(cubicbez, [-1., -1., 1., 0.], 1e-2, 11);
    }

    #[test]
    #[allow(unused_must_use)]
    fn test3() {
        let cubicbez = kurbo::CubicBez::new(
            kurbo::Point::ZERO,
            kurbo::Point::new(0., 0.3),
            kurbo::Point::new(1., 0.3),
            kurbo::Point::new(1., 0.),
        );
        solve_for_cubic(cubicbez, [-1., -1., 1., 0.], 1e-2, 11);
    }

    #[test]
    #[allow(unused_must_use)]
    fn test4() {
        let cubicbez = kurbo::CubicBez::new(
            kurbo::Point::ZERO,
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
            kurbo::Point::new(1., 0.),
        );
        solve_for_cubic(cubicbez, [-1., -1., 1., 0.], 1e-3, 11);
    }
}
