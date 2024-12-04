//! Representation and computation of hyperbeziers.

use std::f64;

use xilem_web::svg::kurbo;

use arrayvec::ArrayVec;
use kurbo::{
    common::GAUSS_LEGENDRE_COEFFS_32, CurveFitSample, ParamCurve, ParamCurveFit, Point, Vec2,
};
use nalgebra::Vector2;
use num_dual::DualNum;

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

impl<D: DualNum<f64> + Copy> HyperbezParams<D> {
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

    // fn report_endpoints(&self) -> [D; 2] {
    //     let mut sum = D::from(0.);
    //     for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
    //         // for (let i = 0; i < co.length; i += 2) {
    //         // let xi = co[i + 1];
    //         // let wi = co[i];
    //         let t = 0.5 + 0.5 * xi;
    //         let q = self.q(D::from(t));
    //         sum += q.powf(-1.5) * *wi;
    //     }
    //     let integral = sum * 0.5;
    //     let q0 = self.q(D::from(0.));
    //     let q1 = self.q(D::from(1.));
    //     [q0.powf(-1.5) / integral, q1.powf(-1.5) / integral]
    // }

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

#[derive(Clone, Copy, Debug)]
pub struct Hyperbezier {
    params: HyperbezParams<f64>,
    p0: Point,
    p1: Point,
    scale_rot: Vec2,
}

impl Hyperbezier {
    /// Create a new hyperbezier curve with given parameters and end points.
    pub fn from_points_params(params: HyperbezParams<f64>, p0: Point, p1: Point) -> Self {
        let uv = params.integrate(1.0);
        let uv = Vec2::new(uv.x, uv.y);
        let uv_scaled = uv / uv.length_squared();
        let d = p1 - p0;
        let scale_rot = Vec2::new(uv_scaled.dot(d), uv_scaled.cross(d));
        Hyperbezier {
            params,
            p0,
            p1,
            scale_rot,
        }
    }

    pub fn params(&self) -> &HyperbezParams<f64> {
        &self.params
    }

    pub fn scale_rot(&self) -> Vec2 {
        self.scale_rot
    }

    pub fn theta(&self, t: f64) -> f64 {
        self.params.theta(t) + self.scale_rot.angle()
    }

    pub fn kappa(&self, t: f64) -> f64 {
        self.params.kappa(t) / self.scale_rot.length()
    }
}

impl ParamCurve for Hyperbezier {
    fn eval(&self, t: f64) -> Point {
        if t == 1.0 {
            self.p1
        } else {
            let s = self.scale_rot;
            let uv = self.params.integrate(t);
            self.p0 + Vec2::new(s.x * uv.x - s.y * uv.y, s.x * uv.y + s.y * uv.x)
        }
    }

    fn start(&self) -> Point {
        self.p0
    }

    fn end(&self) -> Point {
        self.p1
    }

    fn subsegment(&self, range: std::ops::Range<f64>) -> Self {
        let (t0, t1) = (range.start, range.end);
        let dt = t1 - t0;
        let a = self.params.a * dt;
        let b = self.params.b + self.params.a * t0;
        let c = self.params.c * dt * dt;
        let d = (self.params.d + 2. * self.params.c * t0) * dt;
        let e = self.params.c * t0 * t0 + self.params.d * t0 + 1.;
        let s = 1. / e;
        let ps = dt * s * s.sqrt();
        let params = HyperbezParams::new(a * ps, b * ps, c * s, d * s, 1.);
        let p0 = self.eval(t0);
        let p1 = self.eval(t1);
        Hyperbezier::from_points_params(params, p0, p1)
    }
}

impl ParamCurveFit for Hyperbezier {
    fn sample_pt_tangent(&self, t: f64, _sign: f64) -> CurveFitSample {
        let (p, tangent) = self.sample_pt_deriv(t);
        CurveFitSample { p, tangent }
    }

    fn sample_pt_deriv(&self, t: f64) -> (Point, Vec2) {
        let p = self.eval(t);
        let uv = Vec2::from_angle(self.params.theta(t));
        let s = self.scale_rot;
        let d = Vec2::new(s.x * uv.x - s.y * uv.y, s.x * uv.y + s.y * uv.x);
        (p, d)
    }

    fn break_cusp(&self, _: std::ops::Range<f64>) -> Option<f64> {
        None
    }
}

// impl ParamCurveNearest for Hyperbezier {
//     fn nearest(&self, p: Point, accuracy: f64) -> kurbo::Nearest {
//         let p_local = Affine::translate(self.p0.to_vec2())
//             .then_rotate(-self.scale_rot.angle())
//             .then_scale(1. / self.scale_rot.length())
//             * p;

//         // 1. if theta1 < 2pi, check if p_local is in the sweep region between normal0 and normal1
//         //   1a. if true, then subdivide and repeat from 1 for each half
//         //   1b. if false, then either s = 0 or s = 1 is nearest, check and tell
//         // 2. otherwise sibdivide for theta1 = 2pi and repeat from 1

//         if self.params.theta(1.) >= std::f64::consts::TAU {}

//         todo!()
//     }
// }

/// For curve κ = 1/(1-x^2)^2, in range -x..x, what is the expected height?
pub fn forward_scale(x: f64) -> f64 {
    if x == 0.0 {
        return 1.0;
    }
    if x.abs() == 1.0 {
        return 0.0;
    }
    //const a = (1 - x * x) * Math.atanh(x);
    //return 2 * a / (x + a);
    let u = 1. - x * x;
    2. * (u.sqrt() - u) / (x * x)
}

/// Simple inverse of forward_scale
pub fn inv_scale(y: f64) -> f64 {
    ((4. - 4. * y) / (y * y - 4. * y + 4.)).sqrt()
    // let a = 0;
    // let b = 1;
    // for (let i = 0; i < 20; i++) {
    //     let m = (a + b) * 0.5;
    //     if (forward_scale(m) > y) {
    //         a = m;
    //     } else {
    //         b = m;
    //     }
    // }
    // const c = Math.sqrt((4 - 4 * y)/(y * y - 4 * y + 4));
    // console.log(c, (a + b) * 0.5);
    // return (a + b) * 0.5;
}

/// Compute slope / value ratio for endpoint
pub fn inv_scale_slope(y: f64) -> f64 {
    16. * (1. - 1. * y) / (y * y)
}

pub fn k_for_tension(t: f64) -> f64 {
    let b = 0.25;
    1. / (t * t * (b + (1. - b) * t))
}

/// Returns [k0, k1]
pub fn endk_for_quadratic(c: f64, d: f64) -> [f64; 2] {
    let dis = 4. * c - d * d;
    //console.log('evaluating', c, d, 'dis', dis);
    let integral = (4. * c + 2. * d) / (dis * (c + d + 1.).sqrt()) - (2. * d / dis);
    let k0 = 1. / integral;
    let k1 = k0 * (c + d + 1.).powf(-1.5);
    // let gamma = (k1 / k0).cbrt();
    // let beta = 1. / (gamma * gamma) - 1.;
    // let inv_k0 = ((-2. * gamma - 2.) * d + 4. * beta * gamma) / (-d * d - 4. * d + 4. * beta);
    //console.log('kk', inv_k0, 1/k0);
    //console.log('qu', 1/k0 * d * d + (-2 * gamma - 2 + 4/k0) * d + 4 * beta * (gamma - 1/k0));
    //console.log('dd', (-2 * gamma - 2) * d + 4 * beta * gamma);
    [k0, k1]
}

// pub fn copysign(x, y) {
//     const a = Math.abs(x);
//     return y < 0 ? -a : a;
// }

pub fn solve_quadratic(c0: f64, c1: f64, c2: f64) -> ArrayVec<f64, 2> {
    let sc0 = c0 / c2;
    let sc1 = c1 / c2;
    if !(sc0.is_finite() && sc1.is_finite()) {
        let root = -c0 / c1;
        return if root.is_finite() {
            [root].into_iter().collect()
        } else if c0 == 0. && c1 == 0. {
            [0.].into_iter().collect()
        } else {
            ArrayVec::new()
        };
    }
    let arg = sc1 * sc1 - 4. * sc0;
    let root1 = if arg.is_finite() {
        if arg < 0. {
            return ArrayVec::new();
        } else if arg == 0. {
            return [-0.5 * sc1].into_iter().collect();
        }
        -0.5 * (sc1 + arg.sqrt().copysign(sc1))
    } else {
        -sc1
    };
    let root2 = sc0 / root1;
    if root2.is_finite() {
        return if root2 > root1 {
            [root1, root2]
        } else {
            [root2, root1]
        }
        .into_iter()
        .collect();
    }
    [root1].into_iter().collect()
}

/// Solve quadratic parameters for given k0's normalized to integral
/// Returns [c, d]
pub fn quadratic_for_endk(k0: f64, k1: f64) -> [f64; 2] {
    let gamma = (k1 / k0).cbrt();
    let beta = 1. / (gamma * gamma) - 1.;
    let d = (gamma * k0 - 1.) * (2. * gamma + 2.) / gamma;
    /*
    let c0 = 4 * beta * (gamma - ik0);
    let c1 = -2 * gamma - 2 + 4 * ik0;
    let c2 = ik0;
    let roots = solve_quadratic(c0, c1, c2);
    let mut best_err = 1e9;
    let mut best_d = 0;
    for (var cand_d of roots) {
        let mut endk = endk_for_quadratic(beta - cand_d, cand_d);
        let err = Math.pow(endk.k0 - k0, 2) + Math.pow(endk.k1 - k1, 2);
        if (err < best_err) {
            best_d = cand_d;
            best_err = err;
        }
    }
    let d = best_d;
    let c = beta - d;
    let wrong_root = 4 * beta * gamma / (2 * gamma + 2);
    let other_root = c0 / (c2 * wrong_root);
    console.log('other root', other_root, fast);
    */
    [beta - d, d]
}

/// Returns [a, b]
pub fn solve_thetas(th0: f64, th1: f64, c: f64, d: f64, e: f64) -> [f64; 2] {
    let mut a = 0.;
    let mut b = 0.;
    for _ in 0..10 {
        let ths = HyperbezParams::new(a, b, c, d, e).calc_thetas();
        let dth0 = ths[0] - th0;
        let dth1 = ths[1] - th1;
        if dth0.abs() < 1e-9 && dth1.abs() < 1e-1 {
            break;
        }
        const EPS: f64 = 1e-6;
        const IEPS: f64 = 0.5 / EPS;
        let tap = HyperbezParams::new(a + EPS, b, c, d, e).calc_thetas();
        let tam = HyperbezParams::new(a - EPS, b, c, d, e).calc_thetas();
        let tbp = HyperbezParams::new(a, b + EPS, c, d, e).calc_thetas();
        let tbm = HyperbezParams::new(a, b - EPS, c, d, e).calc_thetas();
        let dth0da = IEPS * (tap[0] - tam[0]);
        let dth1da = IEPS * (tap[1] - tam[1]);
        let dth0db = IEPS * (tbp[0] - tbm[0]);
        let dth1db = IEPS * (tbp[1] - tbm[1]);
        let det = dth0da * dth1db - dth1da * dth0db;
        let da = (dth0 * dth1db - dth1 * dth0db) / det;
        let db = (dth0da * dth1 - dth1da * dth0) / det;
        a -= da;
        b -= db;
        //console.log(dth0, dth1, a, b);
    }
    [a, b]
}

pub mod solver {
    use nalgebra::{SVector, Vector5};

    use super::*;

    #[derive(Debug, Clone, Copy)]
    pub struct Solution<const N: usize> {
        pub params: Vector5<f64>,
        pub err: SVector<f64, N>,
        pub iter: usize,
    }

    #[derive(Debug, Clone, Copy)]
    pub enum SolveError<const N: usize> {
        Singularity {
            guess: Vector5<f64>,
            err: SVector<f64, N>,
            iter: usize,
        },
        OutOfIteration {
            guess: Vector5<f64>,
            err: SVector<f64, N>,
        },
    }

    pub type SolveResult<const N: usize> = Result<Solution<N>, SolveError<N>>;

    pub mod analytic;
    pub mod dual;
}
