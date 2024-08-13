use core::f64;
use std::ops::Range;

use approx::AbsDiffEq;
use nalgebra::{allocator::Allocator, Complex, DefaultAllocator, Dim, RealField, Vector2, Vector3};
use num_dual::{hessian, jacobian, Dual2Vec, DualNum, DualVec};

use num_traits::NumOps;
use xilem_web::svg::kurbo::{common::GAUSS_LEGENDRE_COEFFS_32, ParamCurve, Vec2};

pub use spline::hyperbezier::*;

pub fn d_limit(a: f64, b: f64, c: f64) -> Range<f64> {
    if c == 0. {
        return -1.0..10.;
    }
    let end = 2. * c.abs().sqrt();
    let start = -end / (b * c - a).abs().log10().max(4. / 3.);
    start..end
}

pub fn d_limit_rounded(a: f64, b: f64, c: f64) -> Range<f64> {
    let d_limit = d_limit(a, b, c);
    ((d_limit.start * 10.).ceil() / 10.)..(((d_limit.end - f64::EPSILON) * 10.).floor() / 10.)
}

pub fn base(hb: &Hyperbezier) -> Vec2 {
    hb.end() - hb.start()
}

pub fn tension(hb: &Hyperbezier, t: f64) -> f64 {
    let base = base(hb).length();
    let kappa = hb.params().kappa(t);
    base * kappa_to_tau(kappa)
}

pub fn kappa_to_tau(kappa: f64) -> f64 {
    // ((kappa.abs() + 1.).log2() + 1.) / 5.
    let a = 16.;
    let b = 2.;
    let c = -1.9;
    (kappa.abs() + a).log2() / b + c
}

pub fn tau_to_kappa(tau: f64) -> f64 {
    // (tau * 5. - 1.).exp2() - 1.
    let a = 16.;
    let b = 2.;
    let c = -1.9;
    ((tau - c) * b).exp2() - a
}

pub fn integrate_theta(hb: &HyperbezParams, t: f64) -> f64 {
    let c = hb.c();
    let d = hb.d();

    let th_a = hb.th_a();
    let th_b = hb.th_b();

    let q = c * t * t + d * t + 1.;

    let c_sqrt = Complex::from(c).sqrt();
    let q_sqrt = Complex::from(q).sqrt();
    (th_a * q_sqrt / c
        - th_a / c
        - t * th_b
        - (d * th_a - 2. * c * th_b) / (2. * c * c_sqrt)
            * (((d + 2. * c * t) / (2. * c_sqrt * q_sqrt)).atanh() - (d / (2. * c_sqrt)).atanh()))
    .re
}

pub fn hypo_of_qoppa(hb: &Hyperbezier, t: f64) -> f64 {
    let params = hb.params();
    let qoppa = qoppa(params, t);
    let r = 2. / (qoppa * (2. * hb.theta(t)).sin()).abs().sqrt();
    base(hb).length() * r.atan() * f64::consts::FRAC_2_PI
}

pub fn qoppa(hb: &HyperbezParams, t: f64) -> f64 {
    let a = hb.a();
    let b = hb.b();
    let c = hb.c();
    let d = hb.d();

    let q = c * t * t + d * t + 1.;

    a / q.powf(1.5) + 1.5 * (a * t + b) * (2. * c * t + d) / q.powf(2.5)
}

pub fn dust_of_ka_qo(hb: &Hyperbezier, t: f64) -> f64 {
    let params = hb.params();
    let ka_qo = frac_kappa_qoppa(params, t).abs() + 0.125;
    base(hb).length() * ka_qo.atan() * f64::consts::FRAC_2_PI
}

pub fn frac_kappa_qoppa(hb: &HyperbezParams, t: f64) -> f64 {
    let a = hb.a();
    let b = hb.b();
    let c = hb.c();
    let d = hb.d();

    let q = c * t * t + d * t + 1.;

    2. * q * (a * t + b) / (3. * b * (2. * c * t + d) + a * (8. * q - 3. * d * t - 6.))
}

pub fn abs_max_theta(hb: &HyperbezParams) -> f64 {
    let c = hb.c();
    let d = hb.d();

    let th_a = hb.th_a();
    let th_b = hb.th_b();
    let t = (-2. * th_a + d * th_b) / (d * th_a - 2. * c * th_b);

    hb.theta(if (0. ..=1.).contains(&t) { t } else { 1. })
}

pub trait DualNumExt<T: Copy>: DualNum<T> + Copy {
    fn re_mut(&mut self) -> &mut T;
}

impl<T, F, D> DualNumExt<T> for DualVec<T, F, D>
where
    T: DualNum<F> + RealField + Copy,
    D: Dim,
    DualVec<T, F, D>: DualNum<T> + Copy,
    DefaultAllocator: Allocator<D>,
{
    fn re_mut(&mut self) -> &mut T {
        &mut self.re
    }
}

impl<T, F, D> DualNumExt<T> for Dual2Vec<T, F, D>
where
    T: DualNum<F> + RealField + Copy,
    D: Dim,
    Dual2Vec<T, F, D>: DualNum<T> + Copy,
    DefaultAllocator: Allocator<D, D> + Allocator<nalgebra::Const<1>, D>,
    num_dual::Derivative<T, F, nalgebra::Const<1>, D>: Copy,
{
    fn re_mut(&mut self) -> &mut T {
        &mut self.re
    }
}

pub trait Guess {
    type A: Copy;
    type CD: DualNumExt<f64> + NumOps<Self::A>;

    fn from_parts(a: Self::CD, c: Self::CD, d: Self::CD) -> Self;
    fn into_parts(self) -> (Self::A, Self::CD, Self::CD);
}

impl<D: DualNumExt<f64>> Guess for Vector3<D> {
    type A = D;
    type CD = D;

    fn from_parts(a: Self::CD, c: Self::CD, d: Self::CD) -> Self {
        Vector3::new(a, c, d)
    }

    fn into_parts(self) -> (Self::A, Self::CD, Self::CD) {
        (self.x, self.y, self.z)
    }
}

impl<D: DualNumExt<f64>> Guess for (f64, Vector2<D>) {
    type A = f64;
    type CD = D;

    fn from_parts(a: Self::CD, c: Self::CD, d: Self::CD) -> Self {
        (a.re(), Vector2::new(c, d))
    }

    fn into_parts(self) -> (Self::A, Self::CD, Self::CD) {
        (self.0, self.1.x, self.1.y)
    }
}

fn norm_radians<D: DualNumExt<f64>>(mut theta: D) -> D {
    let re = theta.re_mut();
    *re = re.rem_euclid(f64::consts::TAU);
    if *re > f64::consts::PI {
        *re -= f64::consts::TAU;
    } else if *re < -f64::consts::PI {
        *re += f64::consts::TAU;
    }
    theta
}

// pub fn system_for_solving<D: DualNumExt<f64> + NumOps<<G as Guess>::A>, G: Guess<CD = D>>(
pub fn system_for_solving<G: Guess<CD: DualNumExt<f64>>>(
    theta0_i: f64,
    theta1_i: f64,
    kappa0_i: f64,
    kappa1_i: f64,
) -> impl Fn(G) -> G {
    move |guess: G| -> G {
        let (a, c, d) = guess.into_parts();
        let b = kappa0_i;

        // from `HyperbezParams::new`
        let denom = G::CD::from(2.) / (c * 4. - d * d);
        let dd = d * denom;
        let th_a = (c * b * 2. - d * a) * denom;
        let th_b = dd * b - (d * dd * 0.5 + 1.) * a / c;

        // from `HyperbezParams::theta`
        let q1 = c + d + 1.;
        let q1_sqrt = q1.sqrt();
        let theta1 = (th_a + th_b) / q1_sqrt - th_b;
        let kappa1 = th_a / q1_sqrt - ((c * 2. + d) * (th_a + th_b)) / q1 / q1_sqrt * 0.5;

        let theta = |t: f64| {
            let q = c * t * t + d * t + 1.;
            (th_a * t + th_b) / q.sqrt() - th_b
        };

        // from `HyperbezParams::integrate`
        let arg_uv = {
            // TODO: improve accuracy by subdividing in near-cusp cases
            let mut x = G::CD::from(0.);
            let mut y = G::CD::from(0.);
            let u0 = 0.5;
            for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
                let u = u0 + u0 * xi;
                let (dy, dx) = theta(u).sin_cos();
                x += dx * *wi;
                y += dy * *wi;
            }
            y.atan2(x)
        };

        // dbg!(theta0_i, theta1_i, kappa0_i, kappa1_i,);
        // dbg!(theta1, kappa0, kappa1, arg_uv);

        G::from_parts(
            norm_radians(theta1 - theta1_i + theta0_i),
            kappa1 - kappa1_i,
            norm_radians(arg_uv + theta0_i),
        )
    }
}

pub const EPSILON: f64 = 9765625e-10;

pub fn radian_in_line(theta: f64) -> bool {
    (theta % f64::consts::PI).abs_diff_eq(&0., EPSILON)
}

#[must_use]
fn solve_iterate_once(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    guess: Vector3<f64>,
) -> (f64, Option<Vector3<f64>>) {
    if kappa1.abs_diff_ne(&0., EPSILON) {
        let (f, mut jac) = jacobian(system_for_solving(theta0, theta1, kappa0, kappa1), guess);
        let err = f.norm_squared();
        let new_guess = (err.is_finite() && jac.try_inverse_mut()).then(|| guess - jac * f);
        (err, new_guess)
    } else {
        let guess = Vector2::new(guess.y, guess.z);
        let (f, mut jac) = jacobian(
            |guess| {
                let (_, new_guess) =
                    system_for_solving(theta0, theta1, kappa0, 0.)((-kappa0, guess));
                new_guess
            },
            guess,
        );
        let err = f.norm_squared();
        let new_guess = (err.is_finite() && jac.try_inverse_mut())
            .then(|| guess - jac * f)
            .map(|new_guess| Vector3::new(-kappa1, new_guess.x, new_guess.y));
        (err, new_guess)
    }
}

#[must_use]
fn optim_iterate_once(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    guess: Vector3<f64>,
) -> (f64, Option<Vector3<f64>>) {
    if kappa1.abs_diff_ne(&0., EPSILON) {
        let (f, g, mut h) = hessian(
            |guess| {
                system_for_solving(theta0, theta1, kappa0, kappa1)(guess)
                    .into_iter()
                    .map(|f| f.powi(2))
                    .reduce(|a, b| a + b)
                    .unwrap()
            },
            guess,
        );
        let new_guess = (f.is_finite() && h.try_inverse_mut()).then(|| guess - h * g);
        (f, new_guess)
    } else {
        let guess = Vector2::new(guess.y, guess.z);
        let (f, g, mut h) = hessian(
            |guess| {
                let (_, new_guess) =
                    system_for_solving(theta0, theta1, kappa0, 0.)((-kappa0, guess));
                new_guess
                    .into_iter()
                    .map(|f| f.powi(2))
                    .reduce(|a, b| a + b)
                    .unwrap()
            },
            guess,
        );
        let new_guess = (f.is_finite() && h.try_inverse_mut())
            .then(|| guess - h * g)
            .map(|new_guess| Vector3::new(-kappa1, new_guess.x, new_guess.y));
        (f, new_guess)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum SolveError {
    Singularity(f64, usize),
    OutOfIteration([f64; 3], f64),
}

pub fn solve_for_params_exact(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    [a, c, d]: [f64; 3],
    threshold: f64,
    n_iter: usize,
) -> Result<([f64; 4], f64, usize), SolveError> {
    if radian_in_line(theta0) && radian_in_line(theta1) {
        return Ok(([0.; 4], 0., 0));
    }

    let mut guess = Vector3::new(a, c, d);
    let mut err = f64::INFINITY;
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once(theta0, theta1, kappa0, kappa1, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity(new_err, i));
        };
        if new_err <= threshold {
            return Ok(([new_guess.x, kappa0, new_guess.y, new_guess.z], new_err, i));
        }
        guess = new_guess;
        err = new_err;
    }
    Err(SolveError::OutOfIteration([guess.x, guess.y, guess.z], err))
}

pub fn infer_guess(theta0: f64, theta1: f64, kappa0: f64, kappa1: f64, guessed_d: f64) -> [f64; 3] {
    let b = kappa0;
    let d = guessed_d; // guess

    let d_theta = theta1 - theta0;
    let d_theta_2 = d_theta.powi(2);
    let d_theta_3 = d_theta.powi(3);

    if kappa1.abs_diff_ne(&0., EPSILON) {
        let kappa1_2 = kappa1.powi(2);

        // 4 b k1 T^2 + 4 k1 T^3 + 2 d k1 T^3 + T^4
        let determinant =
            4. * b * kappa1 * d_theta_2 + 2. * kappa1 * (2. + d) * d_theta_3 + d_theta.powi(4);

        // 2 b k1 - 2 (1 + d) k1^2 + (2 + d) k1 T + T^2
        let adj =
            2. * b * kappa1 - 2. * (1. + d) * kappa1_2 + (2. + d) * kappa1 * d_theta + d_theta_2;

        let c = (adj + determinant.sqrt()) / (2. * kappa1_2);

        // -b + (1 + c + d)^(3/2) k1
        let a = -b + (1. + c + d).powf(1.5) * kappa1;

        [a, c, d]
    } else {
        let a = -b;

        // d^2/4 + (b (2 + d))/T + b^2/T^2
        let c = d.powi(2) / 4. + b * (2. + d) / d_theta + b.powi(2) / d_theta_2;

        [a, c, d]
    }
}

pub fn solve_for_params_exact_inferred(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    threshold: f64,
    n_iter: usize,
) -> Result<([f64; 4], f64, usize), SolveError> {
    solve_for_params_exact(
        theta0,
        theta1,
        kappa0,
        kappa1,
        [2.8284271238 * kappa1 - kappa0, 1., 0.],
        threshold,
        n_iter,
    )
}

pub fn solve_for_params(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    guess: [f64; 3],
    threshold: f64,
    n_iter: usize,
) -> [Result<HyperbezParams, SolveError>; 2] {
    let params0 =
        solve_for_params_exact(theta0, theta1, -kappa0, -kappa1, guess, threshold, n_iter);
    let params1 = solve_for_params_exact(theta0, theta1, kappa0, -kappa1, guess, threshold, n_iter);

    let hb_from_params =
        |([a, b, c, d], ..): ([f64; 4], f64, usize)| HyperbezParams::new(a, b, c, d);

    let params0 = params0.map(hb_from_params);
    let params1 = params1.map(hb_from_params);

    let f_base = |p: &HyperbezParams| p.integrate(1.).length_squared();

    if params0.as_ref().map(f_base).ok() >= params1.as_ref().map(f_base).ok() {
        [params0, params1]
    } else {
        [params1, params0]
    }
}

pub fn solve_for_params_inferred(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    threshold: f64,
    n_iter: usize,
) -> [Result<HyperbezParams, SolveError>; 4] {
    let hb_from_params =
        |([a, b, c, d], ..): ([f64; 4], f64, usize)| HyperbezParams::new(a, b, c, d);

    let mut params = [0, 1, 2, 3].map(|i| {
        solve_for_params_exact_inferred(
            theta0,
            theta1,
            if i & 0b01 == 0 { -kappa0 } else { kappa0 },
            if i & 0b10 == 0 { -kappa1 } else { kappa1 },
            threshold,
            n_iter,
        )
        .map(hb_from_params)
    });

    let f_base = |p: &HyperbezParams| p.integrate(1.).length_squared();

    params.sort_unstable_by(|a, b| {
        b.as_ref()
            .map(f_base)
            .ok()
            .partial_cmp(&a.as_ref().map(f_base).ok())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    params
}

pub fn optim_for_params_exact(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    [a, c, d]: [f64; 3],
    threshold: f64,
    n_iter: usize,
) -> Result<([f64; 4], f64, usize), SolveError> {
    if radian_in_line(theta0) && radian_in_line(theta1) {
        return Ok(([0.; 4], 0., 0));
    }

    let mut guess = Vector3::new(a, c, d);
    let mut err = f64::INFINITY;
    for i in 0..n_iter {
        let (new_err, new_guess) = optim_iterate_once(theta0, theta1, kappa0, kappa1, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity(new_err, i));
        };
        if new_err <= threshold {
            return Ok(([new_guess.x, kappa0, new_guess.y, new_guess.z], new_err, i));
        }
        guess = new_guess;
        err = new_err;
    }
    Err(SolveError::OutOfIteration([guess.x, guess.y, guess.z], err))
}

pub fn optim_for_params_exact_inferred(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    threshold: f64,
    n_iter: usize,
) -> Result<([f64; 4], f64, usize), SolveError> {
    optim_for_params_exact(
        theta0,
        theta1,
        kappa0,
        kappa1,
        [2. * kappa1 - kappa0, 1., 0.],
        threshold,
        n_iter,
    )
}

pub fn optim_for_params(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    guess: [f64; 3],
    threshold: f64,
    n_iter: usize,
) -> [Result<HyperbezParams, SolveError>; 2] {
    let params0 =
        optim_for_params_exact(theta0, theta1, -kappa0, -kappa1, guess, threshold, n_iter);
    let params1 = optim_for_params_exact(theta0, theta1, kappa0, -kappa1, guess, threshold, n_iter);

    let hb_from_params =
        |([a, b, c, d], ..): ([f64; 4], f64, usize)| HyperbezParams::new(a, b, c, d);

    let params0 = params0.map(hb_from_params);
    let params1 = params1.map(hb_from_params);

    let f_base = |p: &HyperbezParams| p.integrate(1.).length_squared();

    if params0.as_ref().map(f_base).ok() >= params1.as_ref().map(f_base).ok() {
        [params0, params1]
    } else {
        [params1, params0]
    }
}

pub fn optim_for_params_inferred(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    threshold: f64,
    n_iter: usize,
) -> [Result<HyperbezParams, SolveError>; 4] {
    let hb_from_params =
        |([a, b, c, d], ..): ([f64; 4], f64, usize)| HyperbezParams::new(a, b, c, d);

    let mut params = [0, 1, 2, 3].map(|i| {
        optim_for_params_exact_inferred(
            theta0,
            theta1,
            if i & 0b01 == 0 { -kappa0 } else { kappa0 },
            if i & 0b10 == 0 { -kappa1 } else { kappa1 },
            threshold,
            n_iter,
        )
        .map(hb_from_params)
    });

    let f_base = |p: &HyperbezParams| p.integrate(1.).length_squared();

    params.sort_unstable_by(|a, b| {
        b.as_ref()
            .map(f_base)
            .ok()
            .partial_cmp(&a.as_ref().map(f_base).ok())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    params
}

pub struct DebugHbParams<'a>(pub &'a HyperbezParams);

impl std::fmt::Debug for DebugHbParams<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let arg_uv = self.0.integrate(1.).angle();
        f.debug_struct("HyperbezParams")
            .field("theta0", &-arg_uv.to_degrees())
            .field("theta1", &(-arg_uv + self.0.theta(1.)).to_degrees())
            .field("kappa0", &self.0.kappa(0.))
            .field("kappa1", &self.0.kappa(1.))
            .finish()
    }
}

#[test]
fn list_kappas_for_taus() {
    for tau in ((10..20).step_by(2))
        .chain((20..40).step_by(5))
        .chain((40..=400).step_by(10))
        .map(|i| i as f64 / 100.)
    {
        let kappa = tau_to_kappa(tau);
        println!("{tau:>5.2} => {kappa:>5.3}");
    }
}

#[test]
#[expect(unused_must_use)]
fn test_solve() {
    dbg!(solve_for_params(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        [-1., 1., 0.],
        1e-8,
        50,
    ));

    dbg!(solve_for_params(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        [0., 1., 0.],
        1e-4,
        5,
    ));

    dbg!(solve_for_params(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        [-0.1, 1., 0.],
        1e-4,
        5,
    ));

    dbg!(solve_for_params(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        [1., 2., 1.],
        1e-4,
        5,
    ));

    dbg!(solve_for_params_exact(
        1.598187631519324,
        // -78.40181236848068,
        -72.1186270613,
        40.0,
        -113.1370849898476,
        [-19., 1., -1.5],
        1e-4,
        5,
    ));

    dbg!(solve_for_params_exact(
        1.598187631519324,
        // -78.40181236848068,
        -72.1186270613,
        40.0,
        -113.1370849898476,
        [-1., 1., 0.],
        1e-4,
        5,
    ));
}

#[test]
#[expect(unused_must_use)]
fn test_solve_inferred() {
    dbg!(solve_for_params_inferred(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        1e-8,
        50,
    ));

    dbg!(solve_for_params_inferred(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        1e-4,
        5,
    ));

    dbg!(solve_for_params_inferred(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        1e-4,
        5,
    ));

    dbg!(solve_for_params_inferred(
        0.5334841996861125,
        -0.4665158003138875,
        1.0,
        0.7071067811865475,
        1e-4,
        5,
    ));

    dbg!(solve_for_params_exact_inferred(
        1.598187631519324,
        // -78.40181236848068,
        -72.1186270613,
        40.0,
        -113.1370849898476,
        1e-4,
        5,
    ));

    dbg!(solve_for_params_exact_inferred(
        1.598187631519324,
        // -78.40181236848068,
        -72.1186270613,
        40.0,
        -113.1370849898476,
        1e-4,
        5,
    ));
}

#[test]
#[expect(unused_must_use)]
fn test_solve_list() {
    let guess = infer_guess(
        1.598187631519324,
        // -78.40181236848068,
        -72.1186270613,
        40.0,
        -113.1370849898476,
        -1.,
    );

    dbg!(guess);

    let list = [-250.];
    for a in list {
        let init = [a, 1., 0.];
        let soln = solve_for_params_exact(
            1.598187631519324,
            // -78.40181236848068,
            -72.1186270613,
            40.0,
            -113.1370849898476,
            guess,
            1e-4,
            10,
        );
        let opted = optim_for_params_exact(
            1.598187631519324,
            // -78.40181236848068,
            -72.1186270613,
            40.0,
            -113.1370849898476,
            guess,
            1e-4,
            10,
        );
        dbg!(soln, opted);
    }
}
