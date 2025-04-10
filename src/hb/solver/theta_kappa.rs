use std::f64;

use num_traits::Signed;
use xilem_web::svg::kurbo;

use kurbo::{Affine, Point};
use nalgebra::{Vector2, Vector3, Vector5};
use num_dual::{jacobian, DualNum, DualSVec64, DualVec, DualVec64};

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

#[derive(Debug, Clone, Copy)]
struct Solution<const N: usize> {
    params: SVector<f64, N>,
    err: SVector<f64, N>,
    iter: usize,
}

#[derive(Debug, Clone, Copy)]
enum SolveError<const N: usize> {
    Singularity {
        guess: SVector<f64, N>,
        err: SVector<f64, N>,
        iter: usize,
    },
    OutOfIteration {
        guess: SVector<f64, N>,
        err: SVector<f64, N>,
    },
}

type SolveResult<const N: usize> = Result<Solution<N>, SolveError<N>>;

fn solve<const ORDER: usize>(
    f: impl Fn([DualSVec64<ORDER>; ORDER]) -> [DualSVec64<ORDER>; ORDER],
    u: impl Fn([f64; ORDER]) -> [f64; ORDER],
    guess: [f64; ORDER],
    threshold: f64,
    n_iter: usize,
) -> SolveResult<ORDER> {
    let mut guess = SVector::from_data(nalgebra::ArrayStorage([guess]));
    let mut err = guess;
    err.fill(f64::INFINITY);
    for i in 0..n_iter {
        let (new_err, mut jac) = jacobian(
            |guess| SVector::from_data(nalgebra::ArrayStorage([f(guess.data.0[0])])),
            guess,
        );

        tracing::trace!(?guess, ?new_err, ?jac);

        if new_err.norm_squared().is_infinite() || !jac.try_inverse_mut() {
            return Err(SolveError::Singularity {
                guess,
                err: new_err,
                iter: i,
            });
        }

        if new_err.iter().all(|e| e.abs() < threshold) {
            return Ok(Solution {
                params: guess,
                err: new_err,
                iter: i,
            });
        }

        guess = SVector::from_data(nalgebra::ArrayStorage([u(
            (guess - jac * new_err).data.0[0]
        )]));
        err = new_err;
        tracing::trace!(?guess, ?jac);
    }
    Err(SolveError::OutOfIteration { guess, err })
}

#[tracing::instrument(level = "trace", skip_all, ret(level = "trace"))]
pub fn make_hyperbez(
    theta0: f64,
    inner_theta1: f64,
    kappa0: f64,
    kappa1: f64,
    loopy: bool,
) -> HyperbezParams<f64> {
    tracing::trace!(theta0, inner_theta1, kappa0, kappa1, loopy);

    // straight-line
    if theta0.abs() < f64::EPSILON && inner_theta1.abs() < f64::EPSILON {
        return HyperbezParams::new(0., 0., -1., 1., 1.);
    }

    let theta1 = -inner_theta1;
    // both handle on the same side of base
    let same_sided = theta0.is_sign_positive() == theta1.is_sign_negative();
    let traversed_theta = theta1 - theta0 + f64::from(loopy) * f64::consts::TAU.copysign(kappa0);

    tracing::trace!(traversed_theta, same_sided);

    const EPS: f64 = 1e-6;

    // symmetric
    if (theta0 - inner_theta1).abs() < EPS && (kappa0 - kappa1).abs() < EPS {
        let a = 0.;
        let b = kappa0;
        let c = 4. - 2. * (kappa0 + kappa1) / traversed_theta;
        let d = -c;
        tracing::trace!(c, "Symmetric");
        return HyperbezParams::new(a, b, c, d, 1.);
    }

    // antisymmetric
    if (theta0 + inner_theta1).abs() < EPS && (kappa0 + kappa1).abs() < EPS {
        let a = 2. * kappa1;
        let b = kappa0;
        let theta0_5 = 3.4 * theta0 * theta0.sin() / kappa0 - theta0;
        let c = 2. * (1. + (-kappa0 - (theta0_5 * (theta0_5 + 2. * kappa0)).sqrt()) / theta0_5);
        // let c = 2. + (kappa0 + 2. * (theta0 * (theta0 - kappa0)).sqrt()) / theta0;
        tracing::trace!(theta0_5, c, "Antisymmetric");
        let soln = solve(
            |[c]| {
                let d = -c;

                let dual_1 = DualVec::from(1.);
                let hb = HyperbezParams::new(DualVec::from(a), DualVec::from(b), c, d, dual_1);

                let p1 = hb.integrate(dual_1);
                let p1_angle = p1.y.atan2(p1.x);

                [norm_radians(p1_angle + DualVec::from(theta0))]
            },
            |[c]| [c.min(3.99999)],
            [c],
            1e-6,
            10,
        );
        tracing::trace!(?soln);
        let c = match soln {
            Ok(Solution { params: guess, .. })
            | Err(
                SolveError::OutOfIteration { guess, .. } | SolveError::Singularity { guess, .. },
            ) => guess.x,
        };
        let d = -c;
        let hb = HyperbezParams::new(a, b, c, d, 1.);
        let theta0_5 = (theta0_5 + theta0).to_degrees();
        let theta0_5_act = (hb.theta(0.5) + theta0).to_degrees();
        tracing::trace!(theta0_5, theta0_5_act);
        return hb;
    }

    fn criticality_bias(theta: f64, kappa: f64) -> f64 {
        let sin_2theta = (2. * theta).sin();
        (sin_2theta + f64::from(sin_2theta.is_sign_negative())) * (-kappa.abs()).exp2()
    }
    fn criticality_bias_alt(theta: f64, kappa: f64) -> f64 {
        let cos_theta = -theta.cos();
        (cos_theta + f64::from(cos_theta.is_sign_negative())) * (-kappa.abs()).exp2()
    }

    let (a0, a1) = if (theta0 - f64::consts::FRAC_PI_2).abs() < EPS
        || (inner_theta1 - f64::consts::FRAC_PI_2).abs() < EPS
    {
        (
            criticality_bias_alt(theta0, kappa0),
            criticality_bias_alt(inner_theta1, kappa1),
        )
    } else {
        (
            criticality_bias(theta0, kappa0),
            criticality_bias(inner_theta1, kappa1),
        )
    };
    // let a0 = theta0 / kappa0.powi(2);
    // let a1 = -theta1 / kappa1.powi(2);

    // we know, s_critical = -d / (2 c) ; extrema of curvature denominator
    // approximating: s_critical = abs area under control arm 0 / abs area under control arm 0 and 1
    // and assume: d = m c => s_critical = -m / 2
    // therefore: m = -2 s_critical
    let m = if (a0 + a1).abs() <= EPS {
        -1.
    } else {
        -2. * a0 / (a0 + a1)
    };

    tracing::trace!(a0, a1, m);

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

    tracing::trace!(?guess_a, ?guess_c, ?guess_d);

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

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test5() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            -f64::consts::FRAC_PI_4,
            -f64::consts::SQRT_2,
            f64::consts::SQRT_2,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test6() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_3,
            -f64::consts::FRAC_PI_3,
            -f64::consts::SQRT_2,
            f64::consts::SQRT_2,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test7() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            -f64::consts::FRAC_PI_4,
            -2. * f64::consts::SQRT_2,
            2. * f64::consts::SQRT_2,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test8() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_3,
            -f64::consts::FRAC_PI_3,
            -1.9,
            1.9,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test9() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            -f64::consts::FRAC_PI_4,
            -1.,
            1.,
            false,
        );
    }
}
