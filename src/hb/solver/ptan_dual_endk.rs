use std::f64;

use xilem_web::svg::kurbo;

use kurbo::{Affine, ParamCurve, ParamCurveDeriv, Point};
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

pub fn solve_inferring_full(cb: kurbo::CubicBez, threshold: f64, n_iter: usize) -> [f64; 5]{
    let p0_5 = cb.eval(0.5);
    let phi0_5 = cb.deriv().eval(0.5).to_vec2().atan2();
    let c0 = cb.p1.to_vec2();
    let c1 = cb.p3 - cb.p2;
    let theta0 = c0.atan2();
    let theta1 = c1.atan2();

    let p1_angle_i = -theta0;
    let p0_5_i = Affine::rotate(p1_angle_i) * p0_5;
    let phi0_5_i = phi0_5 + p1_angle_i;
    let theta1_i = theta1 + p1_angle_i;

    let [guess_c, guess_d] = {
        let th0 = -theta0;
        let th1 = theta1;
        let d0 = c0.hypot();
        let d1 = c1.hypot();
        let tens0 = d0 * 1.5 * (th0.cos() + 1.);
        let tens1 = d1 * 1.5 * (th1.cos() + 1.);
        let mut k0 = HyperbezParams::<f64>::k_for_tension(tens0);
        let mut k1 = HyperbezParams::<f64>::k_for_tension(tens1);
        let cbr = (k0 / k1).powf(1. / 3.);
        // tracing::trace!(th0, th1, d0, d1, tens0, tens1, k0, k1, cbr);
        fn soft(x: f64) -> f64 {
            (0.5 * (1. + x * x)).sqrt()
        }
        k1 /= soft(cbr);
        k0 /= soft(1. / cbr);
        // tracing::trace!(k0, k1);

        let dc = (cb.p2 - cb.p1).hypot();
        let kmid = dc.powf(1.5);
        let ratio = (d0 / d1).powf(1.5);
        let blend = 0.5 + 0.5 * (3. - 10. * dc).tanh();
        k0 += blend * (kmid / ratio - k0);
        k1 += blend * (kmid * ratio - k1);
        // tracing::trace!(dc, kmid, ratio, blend, k0, k1);
        HyperbezParams::<f64>::quadratic_for_endk(k0, k1)
    };

    let guess_b = make_guess_b(guess_c, guess_d, theta1, theta0);
    let guess_t = (p0_5.x * theta1.abs() + 0.5 * theta0.abs()) / (theta1.abs() + theta0.abs());
    let mut guess = [0., guess_b, guess_c, guess_d, guess_t];
    tracing::trace!(?guess);

    let mut err = f64::INFINITY;
    // for i in 0..n_iter {
    for i in 0..1 {
        let ab = solve_for_ab_exact(theta0, theta1, guess, threshold, n_iter / 2);
        tracing::trace!(?ab);

        if let Ok(Solution { params, .. }) = ab {
            guess = params.data.0[0];
        };
        tracing::trace!(?guess);

        let cdt = solve_for_cdt_exact(p0_5, phi0_5, guess, 1e-2, 5);
        tracing::trace!(?cdt);

        if let Ok(Solution { params, .. }) = cdt {
            guess = params.data.0[0];
        };
        // guess[1] = make_guess_b(guess[2], guess[3], theta1, theta0);
        tracing::trace!(?guess);

        let hb = HyperbezParams::new(guess[0], guess[1], guess[2], guess[3], 1.);
        let p1_r = hb.integrate(1.);
        let p1_angle_r = p1_r.y.atan2(p1_r.x);
        let theta1_r = hb.theta(1.);
        let p0_5_r = hb.integrate(guess[4]);
        let phi0_5_r = hb.theta(guess[4]);

        let new_err = ([
            p0_5_r.x - p0_5_i.x,
            p0_5_r.y - p0_5_i.y,
            phi0_5_r - phi0_5_i,
            theta1_r - theta1_i,
            p1_angle_r - p1_angle_i,
        ]
        .into_iter()
        .map(|e| e.powi(2))
        .sum::<f64>()
            / 5.)
            .sqrt();

        tracing::trace!(i, ?err, ?new_err);

        if new_err < threshold || new_err > err + threshold.powi(2) {
            break;
        }

        err = new_err;
    }

    let abcdt = solve_for_params_exact(p0_5, phi0_5, theta0, theta1, guess, 1e-2, 7);
    tracing::trace!(?abcdt);

    if let Ok(Solution { params, .. }) = abcdt {
        guess = params.data.0[0];
    };
    tracing::trace!(?guess);

    guess
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1() {
        solve_helper(
            Point::new(0.1, 0.4),
            Point::new(0.3, 0.4),
            // solve_with_guess(solve_for_params_exact, [-1., -1., -1., 1.]),
            solve_with_guess(solve_for_params_exact, [3.4, -3., 1., -1.]),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1_1() {
        solve_helper(
            Point::new(0.1, 0.4),
            Point::new(0.3, 0.4),
            // solve_with_guess(solve_for_params_exact, [-1., -1., -1., 1.]),
            // solve_with_guess(solve_for_params_exact, [5.3, -5., -1., 1.]),
            solve_with_guess(solve_for_params_exact, [2.8, -2.2, 2.2, -2.2]),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        solve_helper(
            Point::new(0.1, 0.3),
            Point::new(0.3, 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test3() {
        solve_helper(
            Point::new(0., 0.3),
            Point::new(1., 0.3),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test4() {
        solve_helper(
            Point::new(0.3, 0.15),
            Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [0., -1., -1., 1.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test5() {
        solve_helper(
            Point::new(0.3, 0.15),
            Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [-1., -1., 1., 0.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test6() {
        solve_helper(
            Point::new(0.3, 0.15),
            Point::new(0.7, 0.15),
            solve_with_guess(solve_for_params_exact, [0., -0.7, 1., -1.]),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test7() {
        solve_helper(
            Point::new(0.3, 0.15),
            Point::new(0.7, 0.15),
            solve_inferring(solve_for_params_exact),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test8() {
        solve_helper(
            Point::new(0.3, 0.25),
            Point::new(0.7, 0.25),
            solve_inferring(solve_for_params_exact),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test9() {
        solve_helper(
            Point::new(0.1, 0.4),
            Point::new(0.3, 0.4),
            solve_inferring(solve_for_params_exact),
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test10() {
        solve_helper(
            Point::new(0.65, 0.5),
            Point::new(0.35, 0.5),
            solve_inferring(solve_for_params_exact),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test11() {
        solve_helper(
            Point::new(-0.5, 0.5),
            Point::new(1.5, 0.5),
            solve_inferring(solve_for_params_exact),
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test12() {
        solve_helper(
            Point::new(0.1, 0.4),
            Point::new(0.3, 0.4),
            solve_inferring_full,
        );
        // solve_for_cubic(cubicbez, [0., 0., 1., -1., 0.5]);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test13() {
        solve_helper(
            Point::new(0.65, 0.5),
            Point::new(0.35, 0.5),
            solve_inferring_full,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test14() {
        solve_helper(
            Point::new(-0.5, 0.5),
            Point::new(1.5, 0.5),
            solve_inferring_full,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test15() {
        solve_helper(
            Point::new(0.1, 0.4),
            Point::new(0.9, 0.4),
            solve_inferring_full,
        );
    }
}
