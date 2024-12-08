use std::f64;

use xilem_web::svg::kurbo::{CubicBez, ParamCurve, ParamCurveDeriv, Point, Vec2};

use crate::hb::{k_for_tension, quadratic_for_endk, solver};

pub fn norm_radians<D: num_dual::DualNum<f64> + Copy>(theta: D) -> D {
    let mut re = theta.re().rem_euclid(f64::consts::TAU);
    if re > f64::consts::PI {
        re -= f64::consts::TAU;
    }
    theta + (re - theta.re())
}

pub fn radian_in_line(theta: f64) -> bool {
    (theta % f64::consts::PI).abs() <= 0.01
}

pub fn as_vec2(v: nalgebra::Vector2<f64>) -> Vec2 {
    Vec2::new(v.x, v.y)
}

// MARK: Solver

pub fn solve_helper<R>(p1: Point, p2: Point, solve: impl Fn(CubicBez, f64, usize) -> R) -> R {
    let cubicbez = CubicBez::new(Point::ZERO, p1, p2, Point::new(1., 0.));
    solve(cubicbez, 1e-2, 11)
}

pub fn solve_helper_ext<R>(
    p1: Point,
    p2: Point,
    threshold: f64,
    n_iter: usize,
    solve: impl Fn(CubicBez, f64, usize) -> R,
) -> R {
    let cubicbez = CubicBez::new(Point::ZERO, p1, p2, Point::new(1., 0.));
    solve(cubicbez, threshold, n_iter)
}

fn solve_for_cubic<R: std::fmt::Debug>(
    solve_exact: impl Fn(Point, f64, f64, f64, [f64; 5], f64, usize) -> R,
    cb: CubicBez,
    guess: [f64; 4],
    threshold: f64,
    n_iter: usize,
) -> R {
    let p0_5 = cb.eval(0.5);
    let phi0_5 = cb.deriv().eval(0.5).to_vec2().atan2();
    let theta0 = cb.p1.to_vec2().atan2();
    let theta1 = (cb.p3 - cb.p2).atan2();
    let guess_t = (p0_5.x * theta1.abs() + 0.5 * theta0.abs()) / (theta1.abs() + theta0.abs());
    tracing::trace!(guess_t);
    let res = solve_exact(
        p0_5,
        phi0_5,
        theta0,
        theta1,
        [guess[0], guess[1], guess[2], guess[3], guess_t],
        threshold,
        n_iter,
    );
    dbg!(p0_5);
    dbg!(res)
}

pub fn solve_with_guess<R: std::fmt::Debug>(
    solve_exact: impl Copy + Fn(Point, f64, f64, f64, [f64; 5], f64, usize) -> R,
    guess: [f64; 4],
) -> impl Fn(CubicBez, f64, usize) -> R {
    move |cb, threshold, n_iter| solve_for_cubic(solve_exact, cb, guess, threshold, n_iter)
}

pub fn solve_inferring<R: std::fmt::Debug>(
    solve_exact: impl Copy + Fn(Point, f64, f64, f64, [f64; 5], f64, usize) -> R,
) -> impl Fn(CubicBez, f64, usize) -> R {
    move |cb, threshold, n_iter| {
        let p0_5 = cb.eval(0.5);
        let phi0_5 = cb.deriv().eval(0.5).to_vec2().atan2();
        let c0 = cb.p1.to_vec2();
        let c1 = cb.p3 - cb.p2;
        let theta0 = c0.atan2();
        let theta1 = c1.atan2();

        let [guess_c, guess_d] = {
            let th0 = -theta0;
            let th1 = theta1;
            let d0 = c0.hypot();
            let d1 = c1.hypot();
            let tens0 = d0 * 1.5 * (th0.cos() + 1.);
            let tens1 = d1 * 1.5 * (th1.cos() + 1.);
            let mut k0 = k_for_tension(tens0);
            let mut k1 = k_for_tension(tens1);
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
            quadratic_for_endk(k0, k1)
        };

        let guess_b = make_guess_b(guess_c, guess_d, theta1, theta0);
        let guess_t = (p0_5.x * theta1.abs() + 0.5 * theta0.abs()) / (theta1.abs() + theta0.abs());
        let guess = [0., guess_b, guess_c, guess_d, guess_t];
        tracing::trace!(?guess);

        let cdt = solver::dual::solve_for_cdt_exact(p0_5, phi0_5, guess, 1e-2, 5);
        tracing::trace!(?cdt);

        let [.., guess_c, guess_d, guess_t] = if let Ok(solver::Solution { params, .. }) = cdt {
            params.data.0[0]
        } else {
            guess
        };
        let guess_b = make_guess_b(guess_c, guess_d, theta1, theta0);
        tracing::trace!(guess_b, guess_c, guess_d, guess_t);

        let s = solver::dual::solve_for_ab_exact(
            theta0,
            theta1,
            [0., guess_b, guess_c, guess_d, guess_t],
            1e-2,
            5,
        );
        tracing::trace!(?s);

        // let a_b = 1.5625 + theta1 / (theta0 + theta1 * 0.34);
        // let b = (theta1 - theta0) * 5. / (2. * a_b + 4.);
        // let a = a_b * b;
        // let guess = [0., theta1 * 2.5, -1., 1.];

        let guess = if let Ok(solver::Solution { params, .. }) = s {
            [params.x, params.y, params.z, params.w]
        } else {
            // [a, b, -1., 1.]
            [0., guess_b, guess_c, guess_d]
        };
        dbg!(guess);

        solve_for_cubic(solve_exact, cb, guess, threshold, n_iter)
    }
}

pub fn make_guess_b(guess_c: f64, guess_d: f64, theta1: f64, theta0: f64) -> f64 {
    // (-4 c T + d^2 T)/(2 d - (4 c)/Sqrt[1 + c + d] - (2 d)/Sqrt[1 + c + d])
    let guess_q1_sqrt = (guess_c + guess_d + 1.).sqrt();
    (theta1 - theta0) * (-4. * guess_c + guess_d.powi(2))
        / (2. * guess_d - (4. * guess_c + 2. * guess_d) / guess_q1_sqrt)
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test_norm_radians() {
        let eps = 1e-3;
        for th in [
            eps,
            -eps,
            f64::consts::PI + eps,
            f64::consts::PI - eps,
            2. * f64::consts::PI + eps,
            2. * f64::consts::PI - eps,
            3. * f64::consts::PI + eps,
            3. * f64::consts::PI - eps,
            7. * f64::consts::PI + eps,
            7. * f64::consts::PI - eps,
            11. * f64::consts::PI + eps,
            11. * f64::consts::PI - eps,
            23. * f64::consts::PI + eps,
            23. * f64::consts::PI - eps,
            -f64::consts::PI + eps,
            -f64::consts::PI - eps,
            -2. * f64::consts::PI + eps,
            -2. * f64::consts::PI - eps,
            -3. * f64::consts::PI + eps,
            -3. * f64::consts::PI - eps,
            -7. * f64::consts::PI + eps,
            -7. * f64::consts::PI - eps,
            -11. * f64::consts::PI + eps,
            -11. * f64::consts::PI - eps,
            -23. * f64::consts::PI + eps,
            -23. * f64::consts::PI - eps,
        ] {
            let norm = norm_radians(th);
            tracing::trace!(th, norm);
        }
    }
}
