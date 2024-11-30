use std::f64;

use xilem_web::svg::kurbo::{CubicBez, ParamCurve, ParamCurveDeriv, Point};

use crate::hb_extra::solver;

pub fn norm_radians<D: num_dual::DualNum<f64> + Copy>(theta: D) -> D {
    let mut re = theta.re().rem_euclid(f64::consts::TAU);
    if re > f64::consts::PI {
        re -= f64::consts::TAU;
    }
    theta - (re - theta.re())
}

pub fn radian_in_line(theta: f64) -> bool {
    (theta % f64::consts::PI).abs() <= 0.01
}

pub fn solve_helper<R>(p1: Point, p2: Point, solve: impl Fn(CubicBez, f64, usize) -> R) -> R {
    let cubicbez = CubicBez::new(Point::ZERO, p1, p2, Point::new(1., 0.));
    solve(cubicbez, 1e-4, 11)
}

fn solve_for_cubic<R: std::fmt::Debug>(
    solve_exact: impl Fn(Point, f64, f64, f64, [f64; 5], f64, usize) -> R,
    cb: CubicBez,
    guess: [f64; 4],
    threshold: f64,
    n_iter: usize,
) -> R {
    let p0_5 = cb.eval(0.5);
    let res = solve_exact(
        p0_5,
        cb.deriv().eval(0.5).to_vec2().atan2(),
        cb.p1.to_vec2().atan2(),
        (cb.p3 - cb.p2).atan2(),
        [guess[0], guess[1], guess[2], guess[3], p0_5.x],
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
        let theta0 = cb.p1.to_vec2().atan2();
        let theta1 = (cb.p3 - cb.p2).atan2();

        let guess_t = (p0_5.x * theta1.abs() + 0.5 * theta0.abs()) / (theta1.abs() + theta0.abs());
        let cdt = solver::solve_for_cdt_exact(p0_5, phi0_5, [0., -1., -1., 1., guess_t], 1e-2, 5);
        tracing::trace!(?cdt);

        let [guess_c, guess_d, guess_t] = if let Ok(solver::Solution { params, .. }) = cdt {
            params.data.0[0]
        } else {
            [-1., 1., guess_t]
        };
        dbg!(guess_c, guess_d, guess_t);

        let guess_q1_sqrt = (guess_c + guess_d + 1.).sqrt();
        let guess_b = (theta1 - theta0) * (-4. * guess_c + guess_d.powi(2))
            / (2. * guess_d - (4. * guess_c + 2. * guess_d) / guess_q1_sqrt);
        // let guess_b = (-4 c T + d^2 T)/(2 d - (4 c)/Sqrt[1 + c + d] - (2 d)/Sqrt[1 + c + d]);

        let s = solver::solve_for_ab_exact(
            theta0,
            theta1,
            [0., guess_b, guess_c, guess_d, guess_t],
            1e-2,
            5,
        );
        tracing::trace!(?s);

        let a_b = 1.5625 + theta1 / (theta0 + theta1 * 0.34);
        let b = (theta1 - theta0) * 5. / (2. * a_b + 4.);
        let a = a_b * b;
        // let guess = [0., theta1 * 2.5, -1., 1.];

        let guess = if let Ok(solver::Solution { params, .. }) = s {
            [params.x, params.y, -1., 1.]
        } else {
            [a, b, -1., 1.]
        };
        dbg!(guess);

        solve_for_cubic(solve_exact, cb, guess, threshold, n_iter)
    }
}
