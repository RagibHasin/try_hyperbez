//! Representation and computation of hyperbeziers.

use core::f64;

use xilem_web::svg::kurbo;

use arrayvec::ArrayVec;
use kurbo::{
    common::GAUSS_LEGENDRE_COEFFS_32, Affine, CurveFitSample, ParamCurve, ParamCurveFit, Point,
    Vec2,
};

#[derive(Clone, Copy, Debug)]
pub struct HyperbezParams {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,

    num0: f64,
    num1: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Hyperbezier {
    params: HyperbezParams,
    p0: Point,
    p1: Point,
    scale_rot: Vec2,
}

impl HyperbezParams {
    /// Create a new hyperbezier with the given parameters.
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let denom = 2. / (4. * c - d * d);
        let beta0 = d * denom;
        let num0 = a * (-1. - 0.5 * d * beta0) / c + b * beta0;
        let num1 = (2. * b * c - d * a) * denom;

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

    fn int_helper(&self, t: f64) -> f64 {
        // assumes self.e = 1
        let q = self.c * t * t + self.d * t + 1.;
        (self.num0 + self.num1 * t) / q.sqrt()
    }

    /// Determine the angle for the given parameter.
    ///
    /// This can be interpreted as a Whewell representation of the
    /// curve. The `t` parameter ranges from 0 to 1, and the returned
    /// value is 0 for `t = 0`.
    pub fn theta(&self, t: f64) -> f64 {
        self.int_helper(t) - self.num0
    }

    /// Returns [q, κ]
    fn eval_q(&self, t: f64) -> [f64; 2] {
        let q = self.c * t * t + self.d * t + self.e;
        let k = (self.a * t + self.b) / (q * q.sqrt());
        [q, k]
    }

    pub fn kappa(&self, t: f64) -> f64 {
        self.eval_q(t)[1]
    }

    pub fn q(&self, t: f64) -> f64 {
        self.eval_q(t)[0]
    }

    // fn report_endpoints(&self) -> [f64; 2] {
    //     let mut sum = 0.;
    //     for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
    //         // for (let i = 0; i < co.length; i += 2) {
    //         // let xi = co[i + 1];
    //         // let wi = co[i];
    //         let t = 0.5 + 0.5 * xi;
    //         let q = self.eval_q(t)[0];
    //         sum += wi * q.powf(-1.5);
    //     }
    //     let integral = 0.5 * sum;
    //     let q0 = self.q(0.);
    //     let q1 = self.q(1.);
    //     [q0.powf(-1.5) / integral, q1.powf(-1.5) / integral]
    // }

    /// Evaluate the position of the raw curve.
    ///
    /// This is simply the integral of the Whewell representation,
    /// so that the total arc length is unit, and the initial tangent
    /// is horizontal.
    pub fn integrate(&self, t: f64) -> Vec2 {
        // TODO: improve accuracy by subdividing in near-cusp cases
        let mut xy = Vec2::ZERO;
        let u0 = 0.5 * t;
        for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
            let u = u0 + u0 * xi;
            xy += *wi * Vec2::from_angle(self.theta(u));
        }
        u0 * xy
    }

    pub fn make(a: f64, b: f64, c0: f64, c1: f64) -> Self {
        let c = -(c0 + c1) * (c0 + c1);
        let d = 2. * c0 * (c0 + c1);
        let e = 1. - c0 * c0;
        HyperbezParams::new(a, b, c, d, e)
    }

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

    /// Returns [θ0, θ1]
    pub fn calc_thetas(&self) -> [f64; 2] {
        let p = self.integrate(1.);
        let th0 = p.atan2();
        let th1 = self.theta(1.) - th0;
        [th0, th1]
    }

    pub fn a(&self) -> f64 {
        self.a
    }
    pub fn b(&self) -> f64 {
        self.b
    }
    pub fn c(&self) -> f64 {
        self.c
    }
    pub fn d(&self) -> f64 {
        self.d
    }
    pub fn e(&self) -> f64 {
        self.e
    }
    pub fn num1(&self) -> f64 {
        self.num1
    }
    pub fn num0(&self) -> f64 {
        self.num0
    }
}

impl Hyperbezier {
    /// Create a new hyperbezier curve with given parameters and end points.
    pub fn from_points_params(params: HyperbezParams, p0: Point, p1: Point) -> Self {
        let uv = params.integrate(1.0);
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

    pub fn params(&self) -> &HyperbezParams {
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

impl HyperbezParams {
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

            let p1 = hb.integrate(1.);
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

const PI_2: f64 = f64::consts::FRAC_PI_2;

fn norm_radians(mut theta: f64) -> f64 {
    theta = theta.rem_euclid(f64::consts::TAU);
    if theta > f64::consts::PI {
        theta -= f64::consts::TAU;
    } else if theta < -f64::consts::PI {
        theta += f64::consts::TAU;
    }
    theta
}

use nalgebra as na;

use crate::utils::radian_in_line;

#[must_use]
fn solve_iterate_once(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta1: f64,
    p1_angle: f64,
    guess: [f64; 5],
) -> ([f64; 5], Option<[f64; 5]>) {
    let (err, jac) = HyperbezParams::system_for_solving(p0_5, phi0_5, theta1, p1_angle)(guess);

    let guess_v = na::Vector5::from_data(na::ArrayStorage([guess]));
    let err_v = na::Vector5::from_data(na::ArrayStorage([err]));
    let mut jac_v = na::Matrix5::from_data(na::ArrayStorage(jac)).transpose();
    let new_guess_v = (err_v.norm_squared().is_finite() && jac_v.try_inverse_mut())
        .then(|| guess_v - jac_v * err_v);

    (err, new_guess_v.map(|v| v.data.0[0]))
}

#[derive(Debug, Clone, Copy)]
pub struct Solution {
    pub params: [f64; 5],
    pub err: [f64; 5],
    pub iter: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum SolveError {
    Singularity {
        guess: [f64; 5],
        err: [f64; 5],
        iter: usize,
    },
    OutOfIteration {
        guess: [f64; 5],
        err: [f64; 5],
    },
}

pub fn solve_for_params_exact(
    p0_5: kurbo::Point,
    phi0_5: f64,
    theta0: f64,
    theta1: f64,
    mut guess: [f64; 5],
    threshold: f64,
    n_iter: usize,
) -> Result<Solution, SolveError> {
    if radian_in_line(theta0) && radian_in_line(theta1) {
        return Ok(Solution {
            params: [0., 0., 1., -1., 0.5],
            err: [0.; 5],
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
                guess,
                err: new_err,
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() <= threshold) {
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

#[allow(unused_must_use)]
mod tests {
    use super::*;
    use crate::utils::*;

    #[test]
    fn test1() {
        solve_helper(
            kurbo::Point::new(0.1, 0.4),
            kurbo::Point::new(0.3, 0.4),
            solve_with_guess(solve_for_params_exact, [-1., -1., -1., 1.]),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    fn test2() {
        solve_helper(
            kurbo::Point::new(0.1, 0.3),
            kurbo::Point::new(0.3, 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    fn test3() {
        solve_helper(
            kurbo::Point::new(0., 0.3),
            kurbo::Point::new(1., 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    fn test4() {
        solve_helper(
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [0., -1., -1., 1.]),
        );
    }

    #[test]
    fn test5() {
        solve_helper(
            kurbo::Point::new(0.3, 0.15),
            kurbo::Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    fn listy() {
        use std::io::Write;
        let mut file = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .open(r#"D:\Projects\Long\Playground\try_hyperbez\listy_big.csv"#)
                .unwrap(),
        );
        for b in (0..=200).map(|i| -20. + i as f64 * 0.2) {
            for a in (0..=200).map(|i| -20. + i as f64 * 0.2) {
                let hb = HyperbezParams::new(a, b, -1., 1., 1.);
                let theta1 = hb.theta(1.);
                let p1_arg = hb.integrate(1.).atan2();
                writeln!(&mut file, "{a:.1},{b:.1},{p1_arg:.3},{theta1:.3}");
            }
        }
    }
}

pub mod solver;
