use std::{convert::identity, f64};

use num_dual::{jacobian, DualNum, DualSVec64, DualVec};
use num_traits::Signed;

use crate::utils::*;

use super::*;

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
    threshold: [f64; ORDER],
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

        if !new_err.norm_squared().is_finite() || !jac.try_inverse_mut() {
            return Err(SolveError::Singularity {
                guess,
                err: new_err,
                iter: i,
            });
        }

        if new_err.iter().zip(threshold).all(|(e, eps)| e.abs() < eps) {
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

    let b = kappa0;
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
        // TODO: refine this guess
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

                [norm_radians(p1_angle + theta0)]
            },
            identity,
            [c],
            [1e-2],
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
        (theta / kappa).abs()
    }
    let a0 = criticality_bias(theta0, kappa0);
    let a1 = criticality_bias(inner_theta1, kappa1);

    let s_critical = if a0.is_infinite() {
        1.
    } else if a1.is_infinite() {
        0.
    } else {
        a0 / (a0 + a1)
    };
    tracing::trace!(a0, a1, s_critical);

    // let guess_a_alt = 0.;
    // let (guess_c_alt, guess_d_alt) = {
    //     let kappa012 = kappa0 * kappa1.powi(2);
    //     let l0 = traversed_theta.powi(3) + kappa012;
    //     let l1 = (l0 * kappa012 * traversed_theta.powi(6)).abs();
    //     let l2 =
    //         l1 - kappa012 * traversed_theta.powi(9) + kappa012.powi(2) * traversed_theta.powi(6);
    //     let l2_t = l2.cbrt();
    //     let cbrt2 = 2f64.cbrt();
    //     let d = -2.
    //         + 2. * kappa0 / traversed_theta
    //         + (cbrt2.powi(2) * l2_t) / (traversed_theta.powi(3) * kappa1)
    //         - (2. * cbrt2 * traversed_theta.powi(2) * kappa0 * kappa1) / l2_t;
    //     let c = (cbrt2.powi(2) * l1 * traversed_theta * kappa012
    //         - cbrt2 * l2_t * traversed_theta.powi(9) * kappa012
    //         + 4. * cbrt2 * l2_t * traversed_theta.powi(8) * kappa0 * kappa1.powi(3)
    //         + 2. * cbrt2.powi(2) * traversed_theta.powi(12) * kappa0 * kappa1.powi(3)
    //         - 4. * l2_t.powi(2) * traversed_theta.powi(2) * kappa012.powi(2)
    //         + cbrt2.powi(2) * traversed_theta.powi(10) * kappa012.powi(2)
    //         + cbrt2.powi(2) * traversed_theta.powi(7) * kappa012.powi(3)
    //         + cbrt2 * l1 * (l2_t - 2. * cbrt2 * kappa0 * kappa1.powi(3))
    //         - 4. * l2_t
    //             * traversed_theta.powi(5)
    //             * kappa012
    //             * (l2_t - cbrt2 * kappa0 * kappa1.powi(3))
    //         - 2. * traversed_theta.powi(3)
    //             * (cbrt2.powi(2) * l1 * kappa1 - l2_t.powi(2) * kappa0 * kappa1.powi(4))
    //         + traversed_theta.powi(6)
    //             * (2. * l2_t.powi(2) * kappa1.powi(2)
    //                 - cbrt2 * l2_t * kappa012.powi(2)
    //                 - 2. * cbrt2.powi(2) * kappa0.powi(3) * kappa1.powi(7)))
    //         / (2. * l0 * l2_t.powi(2) * traversed_theta.powi(3) * kappa1.powi(2));
    //     (c, d)
    // };

    let forge = |a: f64| {
        let k012 = (a + kappa0) * kappa1.powi(2);
        let l0 = traversed_theta.powi(3) + k012;
        let l1 = (l0 * k012 * traversed_theta.powi(6)).abs();
        let l2 = l1 + traversed_theta.powi(6) * k012 * (k012 - traversed_theta.powi(3));
        let l2_t = l2.cbrt();
        let l3 = traversed_theta * (a + kappa0) * kappa1;
        let l4 = l3 - 2. * k012;
        let cbrt2 = 2f64.cbrt();
        let d = -2.
            + 2. * kappa0 / traversed_theta
            + (cbrt2.powi(2) * l2_t) / (traversed_theta.powi(3) * kappa1)
            - (2. * cbrt2 * traversed_theta * l3) / l2_t;
        let c = (l1
            * (cbrt2 * l2_t + cbrt2.powi(2) * (l4 - 2. * traversed_theta.powi(3)) * kappa1)
            + traversed_theta.powi(2)
                * (k012 + traversed_theta.powi(3))
                * kappa1.powi(2)
                * (2. * l2_t.powi(2) * (traversed_theta - 2. * kappa0)
                    - cbrt2
                        * l2_t
                        * traversed_theta.powi(3)
                        * (a + kappa0)
                        * (traversed_theta - 4. * kappa1)
                    + cbrt2.powi(2)
                        * l3
                        * traversed_theta.powi(3)
                        * (l4 + 2. * traversed_theta.powi(3))))
            / (2. * l0 * l2_t.powi(2) * traversed_theta.powi(3) * kappa1.powi(2));
        tracing::trace!(a, c, d);
        (c, d)
    };

    let guess_a_alt = if kappa0.is_sign_positive() == kappa1.is_sign_positive() {
        0.
    } else {
        -2. * kappa0
    };
    let guess_a_alt = 2. * kappa0 * (theta0 + theta1).sin();
    let guess_a_alt = kappa1 - kappa0;
    let (guess_c_alt, guess_d_alt) = forge(guess_a_alt);
    // let a = 2. * (theta0 - inner_theta1).sin() * kappa0;
    // let guess_a_alt = -3. * b * (2. * guess_c_alt * s_critical + guess_d_alt)
    //     / (4. * guess_c_alt * s_critical.powi(2) + guess_d_alt * s_critical - 2.);
    // let (guess_c_alt, guess_d_alt) = forge(guess_a_alt);

    fn calc_c<D: DualNum<f64> + Copy>(a: D, b: D, d: D, traversed_theta: f64) -> D {
        let denom = ((d + 2.) * traversed_theta - b * 2.)
            * (-(d - 2.) * traversed_theta * b * 4.
                + b.powi(2) * 4.
                + (a * 16. + (d + 2.).powi(2) * traversed_theta) * traversed_theta)
                .sqrt();
        let avgr = -d * traversed_theta * b * 4.
            + b.powi(2) * 4.
            + (a * 8. + (d.powi(2) - d * 4. - 4.) * traversed_theta) * traversed_theta;
        (avgr - denom) / (8. * traversed_theta.powi(2))
    }

    let soln = solve(
        |[a, d]| {
            let dual_1 = DualVec::from(1.);

            let b = b.into();
            let c = calc_c(a, b, d, traversed_theta);
            let hb = HyperbezParams::new(a, b, c, d, dual_1);

            let p1 = hb.integrate(dual_1);
            let p1_angle = p1.y.atan2(p1.x);

            let kappa1_g = hb.kappa(dual_1);

            [norm_radians(p1_angle + theta0), kappa1_g - kappa1]
        },
        identity,
        [guess_a_alt, guess_d_alt],
        [1e-2, 1e-3],
        10,
    );
    tracing::trace!(?soln);
    let [a, d] = match soln {
        Ok(Solution { params: guess, .. })
        | Err(SolveError::OutOfIteration { guess, .. } | SolveError::Singularity { guess, .. }) => {
            guess.data.0[0]
        }
    };
    let c = calc_c(a, b, d, traversed_theta);
    // let [a, b, c, d] = [guess_a[0], kappa0, guess_c[0], guess_d[0]];
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

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test10() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            f64::consts::PI,
            -f64::consts::SQRT_2,
            -f64::consts::SQRT_2,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test11() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            46f64.to_radians(),
            -f64::consts::SQRT_2,
            -f64::consts::SQRT_2,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test12() {
        let _ = make_hyperbez(
            f64::consts::FRAC_PI_4,
            f64::consts::FRAC_PI_3,
            // 46f64.to_radians(),
            -f64::consts::FRAC_PI_2,
            -f64::consts::FRAC_PI_3 * 2.,
            // -1.7,
            false,
        );
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test13() {
        let _ = make_hyperbez(60f64.to_radians(), -15f64.to_radians(), -1., 1. / 3., false);
    }
}
