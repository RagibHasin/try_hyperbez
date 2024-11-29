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
        let theta0 = cb.p1.to_vec2().atan2();
        let theta1 = (cb.p3 - cb.p2).atan2();

        let s = solver::solve_for_ab_exact(theta0, theta1);
        dbg!(s);

        // let xho = 1. / (theta0 / theta1 + 0.34) + 1. / 0.64;
        let a_b = 1.5625 + theta1 / (theta0 + theta1 * 0.34);
        let b = (theta1 - theta0) * 5. / (2. * a_b + 4.);
        let a = a_b * b;
        // let guess = [0., theta1 * 2.5, -1., 1.];

        let guess = if let Ok(solver::GuessAB { params, .. }) = s {
            [params.x, params.y, -1., 1.]
        } else {
            [a, b, -1., 1.]
        };
        dbg!(guess);

        solve_for_cubic(solve_exact, cb, guess, threshold, n_iter)
    }
}
