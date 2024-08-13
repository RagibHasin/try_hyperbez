use core::f64;

use hyperbez_toy::{hb, lut::*};

#[allow(unused_must_use)]
fn main() {
    for i in 0..N_THETA {
        for j in 0..N_THETA {
            for k in 5..N_KAPPA - 5 {
                for l in 5..N_KAPPA - 5 {
                    // let [theta0, theta1, kappa0, kappa1] = result_from_index(index_from_roll(roll));
                    let [theta0, theta1, kappa0, kappa1] = result_from_index([i, j, k, l]);

                    if hb::radian_in_line(theta0) && hb::radian_in_line(theta1) {
                        println!();
                        continue;
                    }

                    for guessed_d in [-1., -0.5, 0., 0.5, 1.] {
                        let guess = hb::infer_guess(theta0, theta1, kappa0, kappa1, guessed_d);

                        let op = |([a, b, c, d], ..): ([f64; 4], f64, usize)| {
                            let params = hb::HyperbezParams::new(a, b, c, d);

                            let arg_uv = params.integrate(1.).angle();

                            let theta0_s = -arg_uv;
                            let theta1_s = -arg_uv + params.theta(1.);
                            let kappa1_s = params.kappa(1.);

                            let d_theta0 = (theta0_s - theta0) % f64::consts::PI;
                            let d_theta1 = (theta1_s - theta1) % f64::consts::PI;
                            let d_kappa1 = kappa1_s - kappa1;

                            (
                                params,
                                [theta0_s, theta1_s, kappa1_s],
                                [d_theta0, d_theta1, d_kappa1],
                            )
                        };

                        let solved = hb::solve_for_params_exact(
                            theta0, theta1, kappa0, kappa1, guess, 1e-4, 20,
                        )
                        .map(op);

                        let optimd = hb::optim_for_params_exact(
                            theta0, theta1, kappa0, kappa1, guess, 1e-4, 20,
                        )
                        .map(op);

                        dbg!(theta0, theta1, kappa0, kappa1, guess, solved, optimd);
                        println!();
                    }
                }
            }
        }
    }
}
