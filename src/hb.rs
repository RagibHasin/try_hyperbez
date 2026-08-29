//! Representation and computation of hyperbeziers.

use std::f64;

use xilem_web::svg::kurbo;

use arrayvec::ArrayVec;
use kurbo::{
    CurveFitSample, Nearest, ParamCurve, ParamCurveFit, ParamCurveNearest, Point, Vec2,
    common::GAUSS_LEGENDRE_COEFFS_32,
};
use nalgebra::Vector2;
use num_dual::DualNum;

use crate::utils::*;

#[derive(Clone, Copy, Debug)]
pub struct HyperbezParams<D> {
    a: D,
    b: D,
    c: D,
    d: D,

    num0: D,
    num1: D,
}

impl<D: DualNum<Primitive = f64> + Copy> HyperbezParams<D> {
    /// Create a new hyperbezier with the given parameters.
    pub fn new(a: D, b: D, c: D, d: D) -> Self {
        let denom = D::from(2.) / (c * 4. - d * d);
        let beta0 = d * denom;
        let num0 = a * (-d * beta0 * 0.5 - 1.) / c + b * beta0;
        let num1 = (b * c * 2. - d * a) * denom;

        HyperbezParams {
            a,
            b,
            c,
            d,
            num0,
            num1,
        }
    }

    pub fn new_with_endk_ln(a: D, b: D, l0: D, l1: D) -> Self {
        let [c, d] = Self::quadratic_for_endk(l0, l1);
        Self::new(a, b, c, d)
    }

    fn int_helper(&self, t: D) -> D {
        if self.c.is_zero() && self.d.is_zero() {
            self.a * t.powi(2) * 0.5 + self.b * t
        } else {
            (self.num0 + self.num1 * t) / self.q(t).sqrt()
        }
    }

    /// Determine the angle for the given parameter.
    ///
    /// This can be interpreted as a Whewell representation of the
    /// curve. The `t` parameter ranges from 0 to 1, and the returned
    /// value is 0 for `t = 0`.
    pub fn theta(&self, t: D) -> D {
        if self.c.is_zero() && self.d.is_zero() {
            self.int_helper(t)
        } else {
            self.int_helper(t) - self.num0
        }
    }

    pub fn kappa(&self, t: D) -> D {
        let q = self.q(t);
        (self.a * t + self.b) / (q * q.sqrt())
    }

    pub fn deriv_kappa(&self, t: D) -> D {
        let q = self.q(t);
        let q_sqrt = q.sqrt();
        self.a / (q * q_sqrt)
            - (self.a * t + self.b) * (self.c * t * 2. + self.d) * 1.5 / (q.powi(2) * q_sqrt)
    }

    pub fn q(&self, t: D) -> D {
        self.c * t * t + self.d * t + D::from(1.)
    }

    pub fn kappa_extrema(&self) -> ArrayVec<D, 2> {
        fn filter_and_collect<D: DualNum<Primitive = f64> + Copy, const N: usize>(
            res: [D; N],
        ) -> ArrayVec<D, 2> {
            let mut res: ArrayVec<D, 2> = res
                .into_iter()
                .filter(|t| (0.0..1.).contains(&t.re()))
                .collect();
            res.sort_by(|a, b| a.re().total_cmp(&b.re()));
            res
        }

        let a = self.a;
        let b = self.b;
        let c = self.c;
        let d = self.d;

        match (a.re() == 0., c.re() == 0.) {
            (true, true) => ArrayVec::new(),
            (true, false) => filter_and_collect([-d / c * 0.5]),
            (false, true) => filter_and_collect([(a * 2. - b * d * 3.) / (a * d)]),
            (false, false) => {
                let quad_a = -a * c * 4.;
                let quad_b = -a * d - b * c * 6.;
                let quad_c = a * 2. - b * d * 3.;
                let det = (quad_b.powi(2) - quad_a * quad_c * 4.).sqrt();
                if !det.re().is_finite() {
                    return ArrayVec::new();
                }
                filter_and_collect([
                    (-quad_b - det) / quad_a * 0.5,
                    (-quad_b + det) / quad_a * 0.5,
                ])
            }
        }
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

    fn integrate_any<R: std::ops::Add<Output = R> + std::ops::Mul<D, Output = R>>(
        f: impl Fn(D) -> R,
        init: R,
        t: D,
    ) -> R {
        let mut xy = init;
        let u0 = t * 0.5;
        for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
            let u = u0 + u0 * *xi;
            xy = xy + f(u) * D::from(*wi);
        }
        xy * u0
    }

    /// Evaluate the position of the raw curve.
    ///
    /// This is simply the integral of the Whewell representation,
    /// so that the total arc length is unit, and the initial tangent
    /// is horizontal.
    pub fn integrate(&self, t: D) -> Vector2<D> {
        let zero = Vector2::new(D::from(0.), D::from(0.));

        // Self::integrate_any(
        //     |u| {
        //         let (y, x) = self.theta(u).sin_cos();
        //         Vector2::new(x, y)
        //     },
        //     zero,
        //     t,
        // )

        let integrate_self = |end_t| {
            Self::integrate_any(
                |u| {
                    let (y, x) = self.theta(u).sin_cos();
                    Vector2::new(x, y)
                },
                zero,
                end_t,
            )
        };

        let integrate_subseg = |w0, w1, end_t: Option<_>| {
            let theta_w0 = self.theta(w0);
            let subseg = self.subsegment(w0..w1);
            let dw = w1 - w0;
            Self::integrate_any(
                |u| {
                    let (y, x) = (theta_w0 + subseg.theta(u)).sin_cos();
                    Vector2::new(x, y)
                },
                zero,
                end_t.unwrap_or((t - w0) / dw),
            ) * dw
        };

        match *self.kappa_extrema().as_slice() {
            [] => integrate_self(t),
            [w0, ..] if t.re() <= w0.re() => integrate_self(t),
            [w0] => integrate_subseg(w0, D::from(1.), None) + integrate_self(w0),
            [w0, w1] if (w0.re()..w1.re()).contains(&t.re()) => {
                integrate_subseg(w0, w1, None) + integrate_self(w0)
            }
            [w0, w1] => {
                integrate_subseg(w1, D::from(1.), None)
                    + integrate_subseg(w0, w1, Some(D::from(1.)))
                    + integrate_self(w0)
            }
            _ => unreachable!(),
        }
    }

    pub fn make(a: D, b: D, c0: D, c1: D) -> Self {
        let c = -(c0 + c1) * (c0 + c1);
        let d = c0 * (c0 + c1) * 2.;
        let e = -c0 * c0 + 1.;
        let e32 = e * e.sqrt();
        HyperbezParams::new(a / e32, b / e32, c / e, d / e)
    }

    /// Returns [θ0, θ1]
    pub fn calc_thetas(&self) -> [D; 2] {
        let p = self.integrate(D::from(1.));
        let th0 = p.y.atan2(p.x);
        let th1 = self.theta(D::from(1.)) - th0;
        [th0, th1]
    }

    pub fn subsegment(&self, range: std::ops::Range<D>) -> Self {
        let (t0, t1) = (range.start, range.end);
        let dt = t1 - t0;
        let a = self.a * dt;
        let b = self.b + self.a * t0;
        let c = self.c * dt * dt;
        let d = (self.d + self.c * t0 * 2.) * dt;
        let e = self.c * t0 * t0 + self.d * t0 + 1.;
        let s = D::from(1.) / e;
        let ps = dt * s * s.sqrt();
        HyperbezParams::new(a * ps, b * ps, c * s, d * s)
    }

    pub fn k_for_tension(t: D) -> D {
        let b = 0.25;
        (t * t * (t * (1. - b) + b)).recip()
    }

    /// Returns [k0, k1]
    pub fn endk_for_quadratic(c: D, d: D) -> [D; 2] {
        // let dis = c * 4. - d.powi(2);
        // let integral = (c * 4. + d * 2.) / (dis * (c + d + 1.).sqrt()) - (d * 2. / dis);
        // let k0 = integral.recip();
        // let k1 = k0 * (c + d + 1.).powf(-1.5);

        let beta = c + d;
        let gamma = (beta + 1.).sqrt().recip();
        let k0 = (d * gamma / (gamma * 2. + 2.) + 1.) / gamma;
        let k1 = gamma.powi(3) * k0;

        [k0, k1]
    }

    /// Solve quadratic parameters for given k0's normalized to integral
    /// Returns [c, d, e]
    pub fn quadratic_for_endk_ln(l0: D, l1: D) -> [D; 3] {
        let e_ln = -(l0 * 2. + l1) / 3. - f64::consts::LN_2;
        let e = e_ln.exp();
        let gamma_ln = (l1 - l0) / 3.;
        let beta = (gamma_ln * -2. + e_ln).exp() - e;
        let d = (((l0 * 2. + l1) / 3.).exp() - 1.)
            * ((-gamma_ln + e_ln + f64::consts::LN_2).exp() + e * 2.);

        [beta - d, d, e]
    }

    /// Solve quadratic parameters for given k0's normalized to integral
    /// Returns [c, d]
    pub fn quadratic_for_endk(k0: D, k1: D) -> [D; 2] {
        let gamma = (k1 / k0).cbrt();
        let beta = gamma.powi(-2) - 1.;
        let d = (gamma * k0 - 1.) * (gamma * 2. + 2.) / gamma;

        [beta - d, d]
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
    pub fn num0(&self) -> D {
        self.num0
    }
    pub fn num1(&self) -> D {
        self.num1
    }
    pub fn endk(&self) -> [D; 2] {
        Self::endk_for_quadratic(self.c, self.d)
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
        let uv = as_vec2(params.integrate(1.0));
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
        let params = self.params.subsegment(range.clone());
        let p0 = self.eval(range.start);
        let p1 = self.eval(range.end);
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

impl ParamCurveNearest for Hyperbezier {
    fn nearest(&self, p: Point, accuracy: f64) -> Nearest {
        let p = kurbo::Affine::translate(self.p0.to_vec2())
            .then_rotate(-self.scale_rot.angle())
            .then_scale(1. / self.scale_rot.length())
            * p;

        let g = |s| {
            (
                s,
                (p.to_vec2() - as_vec2(self.params.integrate(s)))
                    .dot(Vec2::from_angle(self.params.theta(s))),
            )
        };
        let extrema = [0.]
            .into_iter()
            .chain(self.params.kappa_extrema())
            .chain([1.])
            .collect::<ArrayVec<_, 4>>();
        let samples = extrema
            .array_windows()
            .copied()
            .flat_map(|[a, b]| [a, a + (b - a) / 3., b - (b - a) / 3.])
            .chain([1.])
            .map(g)
            .collect::<ArrayVec<_, 10>>();

        let max_iter = (-accuracy.log2()).ceil() as u8;
        let hit = |s| Nearest {
            distance_sq: self.scale_rot.hypot2()
                * (p - as_vec2(self.params.integrate(s))).to_vec2().hypot2(),
            t: s,
        };

        let filtered = samples.array_windows().copied().flat_map(
            |[(mut a_s, mut a_g), (mut b_s, mut b_g)]| {
                if a_g.signum() == b_g.signum() {
                    return None;
                }

                for _ in 0..max_iter {
                    if b_s - a_s <= accuracy {
                        return Some((b_s + a_s) / 2.);
                    }

                    let (c_s, c_g) = g(a_s + a_g * (b_s - a_s) / (a_g - b_g));
                    if c_g.signum() == a_g.signum() {
                        (a_s, a_g) = (c_s, c_g);
                    } else {
                        (b_s, b_g) = (c_s, c_g);
                    }
                }
                None
            },
        );

        [0.].into_iter()
            .chain(filtered)
            .chain([1.])
            .map(hit)
            .min_by(|a, b| a.distance_sq.total_cmp(&b.distance_sq))
            .unwrap()
    }
}

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

/// Returns [a, b]
pub fn solve_thetas(th0: f64, th1: f64, c: f64, d: f64) -> [f64; 2] {
    let mut a = 0.;
    let mut b = 0.;
    for _ in 0..10 {
        let ths = HyperbezParams::new(a, b, c, d).calc_thetas();
        let dth0 = ths[0] - th0;
        let dth1 = ths[1] - th1;
        if dth0.abs() < 1e-9 && dth1.abs() < 1e-1 {
            break;
        }
        const EPS: f64 = 1e-6;
        const IEPS: f64 = 0.5 / EPS;
        let tap = HyperbezParams::new(a + EPS, b, c, d).calc_thetas();
        let tam = HyperbezParams::new(a - EPS, b, c, d).calc_thetas();
        let tbp = HyperbezParams::new(a, b + EPS, c, d).calc_thetas();
        let tbm = HyperbezParams::new(a, b - EPS, c, d).calc_thetas();
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

pub mod render;

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

    pub mod coprop_dual;
    pub mod ladder_alt;
    pub mod ptan_analytic;
    pub mod ptan_dual;
    pub mod ptan_dual_endk;
    pub mod q0q1;
    pub mod theta_kappa;
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1() {
        let hb = HyperbezParams::new(0., -1., -1., 1.);
        let seg0 = hb.subsegment(0.0..0.5);
        let seg1 = hb.subsegment(0.5..1.);
        tracing::trace!(?seg0, ?seg1);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        let hb = HyperbezParams::new(
            0.,
            -0.004213609035097804,
            15.592889809259216,
            -7.892748222945392,
        );
        let extrema = hb.kappa_extrema();
        let seg0 = hb.subsegment(0.0..0.5);
        let seg1 = hb.subsegment(0.5..1.);
        tracing::trace!(?extrema, ?seg0, ?seg1);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test3() {
        let hb = HyperbezParams::new(8.2, -6., 3.4, -3.4);
        let extrema = hb.kappa_extrema();
        let seg0 = hb.subsegment(0.0..0.5);
        let seg1 = hb.subsegment(0.5..1.);
        tracing::trace!(?extrema, ?seg0, ?seg1);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test4() {
        let hb = HyperbezParams::new(0., -3f64.sqrt() / 2., -0.9, 0.9);
        let extrema = hb.kappa_extrema();
        tracing::trace!(?extrema);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test_hit0() {
        let hb = Hyperbezier::from_points_params(
            HyperbezParams::new(0., -1.5, 0.18, -0.18),
            Point::ZERO,
            Point::new(500., 0.),
        );
        let test_p = [
            Point::new(-1., 0.),
            Point::new(501., 0.),
            Point::new(250., 0.),
            Point::new(150., 0.),
            Point::new(350., 0.),
            Point::new(250., 200.),
            Point::new(150., 200.),
            Point::new(350., 200.),
        ];
        let accuracy = 0.001953125;
        for p in test_p {
            let result = hb.nearest(p, accuracy);
            tracing::trace!(%p, ?result);
        }
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test_hit1() {
        let hb = Hyperbezier::from_points_params(
            HyperbezParams::new(32.790, -6.350, 25.950, -9.476),
            Point::ZERO,
            Point::new(500., 0.),
        );
        let test_p = [
            Point::new(325., -100.),
            Point::new(350., -80.),
            Point::new(390., -70.),
            Point::new(420., -50.),
            Point::new(445., -35.),
            Point::new(550., -100.),
        ];
        let accuracy = 0.005;
        for p in test_p {
            let result = hb.nearest(p, accuracy);
            tracing::trace!(%p, ?result);
        }
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn kootikoo() {
        for l0 in -10..=10 {
            for l1 in -10..=10 {
                let l0 = l0 as f64;
                let l1 = l1 as f64;
                let [ln_c, ln_d, ln_e] = HyperbezParams::<f64>::quadratic_for_endk_ln(l0, l1);
                let ln_c_n = ln_c / ln_e;
                let ln_d_n = ln_d / ln_e;
                let k0 = l0.exp();
                let k1 = l1.exp();
                let [j_c, j_d] = HyperbezParams::<f64>::quadratic_for_endk(k0, k1);
                tracing::trace!(l0, l1, k0, k1, ln_c, ln_d, ln_e, ln_c_n, ln_d_n, j_c, j_d);
            }
        }
    }
}
