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

pub fn make_hyperbez(
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    loopy: bool,
) -> HyperbezParams<f64> {
    if theta0.abs() < f64::EPSILON && theta1.abs() < f64::EPSILON {
        return HyperbezParams::new(0., 0., -1., 1., 1.);
    }

    let a0 = theta0 / kappa0.powi(2);
    let a1 = -theta1 / kappa1.powi(2);
    tracing::trace!(a0, a1);

    // both handle on the same side of base
    let same_sided = theta0.is_sign_positive() == theta1.is_sign_negative();
    let traversed_theta = theta1 - theta0 + f64::from(loopy) * f64::consts::TAU.copysign(kappa0);

    tracing::trace!(traversed_theta, same_sided, loopy);

    // we know, s_critical = -d / (2 c) ; extrema of curvature denominator
    // approximating: s_critical = abs area under control arm 0 / abs area under control arm 0 and 1
    // and assume: d = m c => s_critical = -m / 2
    // therefore: m = -2 s_critical
    let m = if (a0 + a1).abs() <= 1e-6 {
        -1.
    } else {
        -2. * a0 / (a0 + a1)
    };

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
        let _ = make_hyperbez(0., -0., -1., -1., false);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        let _ = make_hyperbez(0.0174532925, -0.0174532925, -1., -1., false);
        let _ = make_hyperbez(-0.0174532925, 0.0174532925, 1., 1., false);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test3() {
        let _ = make_hyperbez(-0.0174532925, 0.0174532925, 1., 1., false);
    }
}
