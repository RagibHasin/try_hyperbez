use std::f64;

use xilem_web::svg::kurbo;

use kurbo::{Affine, Point};
use nalgebra::{Vector2, Vector3, Vector5};
use num_dual::{jacobian, DualNum, DualVec64};

use crate::utils::*;

use super::*;

fn system_for_solving<D: DualNum<f64> + Copy>(
    p0_5_i: Point,
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
    p0_5: Point,
    phi0_5: f64,
    theta1: f64,
    p1_angle: f64,
    guess: Vector5<f64>,
) -> (Vector5<f64>, Option<Vector5<f64>>) {
    let (f, mut jac) = jacobian(system_for_solving(p0_5, phi0_5, theta1, p1_angle), guess);
    // jac.transpose_mut();
    // dbg!(&f, &jac);
    let new_guess =
        (f.norm_squared().is_finite() && jac.try_inverse_mut()).then(|| guess - jac * f);
    (f, new_guess)
}

pub fn solve_for_params_exact(
    p0_5: Point,
    phi0_5: f64,
    theta0: f64,
    theta1: f64,
    guess: [f64; 5],
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

    let mut guess = Vector5::from_data(nalgebra::ArrayStorage([guess]));
    let mut err = Vector5::from_data(nalgebra::ArrayStorage([[f64::INFINITY; 5]]));
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once(p0_5, phi0_5, theta1, p1_angle, guess);
        tracing::trace!(i, ?new_err, ?new_guess);
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

#[must_use]
fn solve_iterate_once_for_ab(
    theta1: f64,
    p1_angle: f64,
    guess: Vector5<f64>,
) -> (Vector2<f64>, Option<Vector2<f64>>) {
    let (f, mut jac) = jacobian(
        |guess_i: Vector2<DualVec64<nalgebra::U2>>| {
            let result = system_for_solving(Point::ZERO, 0., theta1, p1_angle)(Vector5::new(
                guess_i.x,
                guess_i.y,
                guess.z.into(),
                guess.w.into(),
                guess.a.into(),
            ));
            Vector2::new(result.w, result.a)
        },
        guess.xy(),
    );
    let new_guess =
        (f.norm_squared().is_finite() && jac.try_inverse_mut()).then(|| guess.xy() - jac * f);
    (f, new_guess)
}

pub fn solve_for_ab_exact(
    theta0: f64,
    theta1: f64,
    guess: [f64; 5],
    threshold: f64,
    n_iter: usize,
) -> SolveResult<2> {
    if radian_in_line(theta0) && radian_in_line(theta1) {
        return Ok(Solution {
            params: Vector5::new(0., 0., guess[2], guess[3], guess[4]),
            err: Vector2::zeros(),
            iter: 0,
        });
    }

    let p1_angle = -theta0;
    let theta1 = theta1 + p1_angle;

    let mut guess = Vector5::from_data(nalgebra::ArrayStorage([guess]));
    let mut err = Vector2::from_data(nalgebra::ArrayStorage([[f64::INFINITY; 2]]));
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once_for_ab(theta1, p1_angle, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity {
                guess,
                err: new_err,
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() < threshold) {
            return Ok(Solution {
                params: Vector5::new(new_guess.x, new_guess.y, guess[2], guess[3], guess[4]),
                err: new_err,
                iter: i,
            });
        }
        guess.x = new_guess.x;
        guess.y = new_guess.y;
        err = new_err;
    }
    Err(SolveError::OutOfIteration { guess, err })
}

#[must_use]
fn solve_iterate_once_for_cdt(
    p0_5: Point,
    phi0_5: f64,
    guess: Vector5<f64>,
) -> (Vector3<f64>, Option<Vector3<f64>>) {
    let guess_o = Vector3::new(guess.z, guess.w, guess.a);
    let (f, mut jac) = jacobian(
        |guess_i: Vector3<DualVec64<nalgebra::U3>>| {
            let [[p0_5_x_r, p0_5_y_r, phi0_5_r, _, p1_angle_r]] =
                system_for_solving(Point::ZERO, 0., 0., 0.)(Vector5::new(
                    guess.x.into(),
                    guess.y.into(),
                    guess_i.x,
                    guess_i.y,
                    guess_i.z,
                ))
                .data
                .0;
            let theta0_r = -p1_angle_r;
            let p0_5_x_o = p0_5_x_r * theta0_r.cos() - p0_5_y_r * theta0_r.sin() - p0_5.x;
            let p0_5_y_o = p0_5_x_r * theta0_r.sin() + p0_5_y_r * theta0_r.cos() - p0_5.y;
            let phi0_5_o = norm_radians(phi0_5_r + theta0_r - phi0_5);
            // tracing::trace!(
            //     ?p0_5_x_r,
            //     ?p0_5_y_r,
            //     ?phi0_5_r,
            //     ?p1_angle_r,
            //     ?p0_5_x_o,
            //     ?p0_5_y_o,
            //     ?phi0_5_o,
            // );

            Vector3::new(p0_5_x_o, p0_5_y_o, phi0_5_o)
        },
        guess_o,
    );
    let new_guess =
        (f.norm_squared().is_finite() && jac.try_inverse_mut()).then(|| guess_o - jac * f);
    (f, new_guess)
}

pub fn solve_for_cdt_exact(
    p0_5: Point,
    phi0_5: f64,
    guess: [f64; 5],
    threshold: f64,
    n_iter: usize,
) -> SolveResult<3> {
    let mut guess = Vector5::from_data(nalgebra::ArrayStorage([guess]));
    let mut err = Vector3::from_data(nalgebra::ArrayStorage([[f64::INFINITY; 3]]));
    for i in 0..n_iter {
        let (new_err, new_guess) = solve_iterate_once_for_cdt(p0_5, phi0_5, guess);
        let Some(new_guess) = new_guess else {
            return Err(SolveError::Singularity {
                guess,
                err: new_err,
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() < threshold) {
            return Ok(Solution {
                params: Vector5::new(guess[0], guess[1], new_guess.x, new_guess.y, new_guess.z),
                err: new_err,
                iter: i,
            });
        }
        guess.z = new_guess.x;
        guess.w = new_guess.y;
        guess.a = new_guess.z;
        err = new_err;
    }
    Err(SolveError::OutOfIteration { guess, err })
}

#[tracing::instrument(ret(level = tracing::Level::TRACE))]
fn arm_limit_from_thetas(a: f64, b: f64) -> f64 {
    fn f_normally(a: f64, b: f64) -> f64 {
        a.sin() / (a + b).sin()
    }

    fn f_regular(a: f64, b: f64) -> f64 {
        let normally = f_normally(a, b);
        let coeff_ab_max = ((a + b) * f64::consts::FRAC_1_PI).powi(36);
        let ang_3deg = 3f64.to_radians();
        // mtor = mitigator
        let mted_a = a * 14. / 15. + ang_3deg;
        let mted_b = b * 14. / 15. + ang_3deg;
        let mtor_ab_max = f_normally(mted_a, mted_b);
        let mted_ab_max = (1. - coeff_ab_max) * normally + coeff_ab_max * mtor_ab_max;

        mted_ab_max + (-20. * mted_ab_max - 1.).exp2()
    }

    let a = a.abs();
    let b = b.abs();

    let ab = a + b;
    use f64::consts::PI;
    if ab == 0. {
        0.5
    } else if ab <= PI {
        f_regular(a, b)
    } else {
        f_regular(a * PI / ab, b * PI / ab).powi(2) / f_regular(PI - a, PI - b)
    }
}

fn kappa_from_k(k: f64) -> f64 {
    k.exp() - 1. / (k + 1.)
}

/// (k, κ)
#[tracing::instrument(ret(level = tracing::Level::TRACE))]
fn kappa_from_arm_limit(arm: Vec2, limit: f64) -> (f64, f64) {
    let theta = arm.angle().abs();
    let k_circle = 0.5 * theta.sin() / limit - 1.;
    let a = -2. * theta / kappa_from_k(k_circle);
    let k = arm.length() / limit - 1.;
    let kappa = a * kappa_from_k(k);
    tracing::trace!(theta, k_circle, a, k, kappa);
    (k, kappa.abs())
}

pub fn make_hyperbez(cb: kurbo::CubicBez) -> HyperbezParams<f64> {
    let c0 = cb.p1.to_vec2();
    let c1 = cb.p3 - cb.p2;
    let theta0 = c0.atan2();
    let theta1 = c1.atan2();

    // // TODO: special-case zero components of control arms
    // let a0 = (c0.x * c0.y).abs();
    // let a1 = (c1.x * c1.y).abs();

    let (a0, a1) = if c0.x == 0. || c1.x == 0. {
        (c0.y.abs() * c0.length(), c1.y.abs() * c1.length())
    } else if c0.y == 0. || c1.y == 0. {
        (c0.x.abs() * c0.length(), c1.x.abs() * c1.length())
    } else {
        ((c0.x * c0.y).abs(), (c1.x * c1.y).abs())
    };
    tracing::trace!(a0, a1);

    let inner_theta1 = -theta1;
    let lim0 = arm_limit_from_thetas(inner_theta1, theta0);
    let lim1 = arm_limit_from_thetas(theta0, inner_theta1);
    let (k0, curv0) = kappa_from_arm_limit(c0, lim0);
    let (k1, curv1) = kappa_from_arm_limit(c1, lim1);

    // both handle on the same side of base
    let same_sided = theta0.is_sign_positive() == theta1.is_sign_negative();
    let loopy = same_sided && k0 > 0. && k1 > 0.;
    let kappa0 = curv0.copysign(theta0 * if loopy { 1. } else { -1. });
    let kappa1 = curv1.copysign(kappa0 * if same_sided { 1. } else { -1. });
    let traversed_theta = theta1 - theta0 + f64::from(loopy) * f64::consts::TAU.copysign(kappa0);

    tracing::trace!(
        theta0,
        theta1,
        lim0,
        lim1,
        k0,
        k1,
        kappa0,
        kappa1,
        traversed_theta,
        same_sided,
        loopy
    );

    // we know, s_critical = -d / (2 c) ; extrema of curvature denominator
    // approximating: s_critical = abs area under control arm 0 / abs area under control arm 0 and 1
    // and assume: d = m c => s_critical = -m / 2
    // therefore: m = -2 s_critical
    let m = -2. * a0 / (a0 + a1);

    if (m + 1.).abs() < 1e-6 {
        let a = 0.;
        let b = kappa0;
        let c = 4. - 2. * (kappa0 + kappa1) / traversed_theta;
        let d = -c;
        return HyperbezParams::new(a, b, c, d, 1.);
    }

    let guess_c = {
        let a = (traversed_theta * m - 2. * (m + 1.) * kappa1).powi(2);
        let b = -4.
            * (traversed_theta.powi(2)
                + kappa0 * traversed_theta * m
                + kappa1 * traversed_theta * (2. + 3. * m)
                - 2. * (1. + m) * kappa1 * (kappa0 + kappa1));
        let c = 4. * (kappa0 + kappa1) * (kappa0 + kappa1 - 2. * traversed_theta);
        let det = (b.powi(2) - 4. * a * c).sqrt();
        [-b + det, -b - det].map(|x| x * 0.5 / a)
    };
    let guess_a = guess_c.map(|c| -kappa0 + kappa1 * (c * (m + 1.) + 1.).powf(1.5));
    let guess_d = guess_c.map(|c| m * c);

    tracing::trace!(m, ?guess_a, ?guess_c, ?guess_d);

    let [a, b, c, d] = [guess_a[0], kappa0, guess_c[0], guess_d[0]];
    HyperbezParams::new(a, b, c, d, 1.)
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.1, 0.5),
            Point::new(0.9, 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0., 0.5),
            Point::new(1., 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test3() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.0, 0.5),
            Point::new(1., 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test4() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.004, 0.5),
            Point::new(1.004, 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test5() {
        // fails
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.02, 0.5),
            Point::new(1.02, 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test6() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.2, 0.05),
            Point::new(1.2, 0.05),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test7() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.2, 0.1),
            Point::new(1.2, 0.1),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test8() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.51, 0.5),
            Point::new(0.49, 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test9() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.02, 0.5),
            Point::new(1.02, 0.5),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test10() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(-0.25, 0.02),
            Point::new(1.25, 0.02),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test11() {
        let arm = arm_limit_from_thetas(f64::consts::FRAC_PI_2, f64::consts::FRAC_PI_2);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test12() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.1, -0.1),
            Point::new(0.9, -0.1),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test13() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.1, 0.1),
            Point::new(0.9, 0.1),
            Point::new(1., 0.),
        ));
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test14() {
        let hb = make_hyperbez(kurbo::CubicBez::new(
            Point::ZERO,
            Point::new(0.2, 0.5),
            Point::new(0.8, -0.5),
            Point::new(1., 0.),
        ));
    }
}
