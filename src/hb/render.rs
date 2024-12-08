use super::*;

use kurbo::{Affine, BezPath, CubicBez, ParamCurveDeriv};
use nalgebra::{ArrayStorage, Matrix4, Vector4};

impl Hyperbezier {
    pub fn render(&self, accuracy: f64) -> BezPath {
        Affine::rotate(self.scale_rot.atan2()).then_scale(self.scale_rot.hypot())
            * self.params.render(accuracy)
    }
}

impl HyperbezParams<f64> {
    pub fn render(&self, accuracy: f64) -> BezPath {
        let deflection = self.theta(1.).abs().ceil();
        let step_size = deflection.recip();
        let mut path = BezPath::new();
        path.move_to((0., 0.));

        for i in 0..deflection as usize {
            let start = i as f64 * step_size;
            let end = start + step_size;
            let subseg = self.subsegment(start..end);
            let p0 = as_vec2(self.integrate(start));
            let p1 = as_vec2(self.integrate(end));
            let uv = p1 - p0;
            let rendered_subseg = subseg.render_to_cubic(accuracy).expect("not to hit error"); // TODO: fix it
            let CubicBez { p1, p2, p3, .. } = Affine::rotate(uv.atan2())
                .then_scale(uv.hypot())
                .then_translate(p0)
                * rendered_subseg;
            path.curve_to(p1, p2, p3);
        }

        path
    }

    pub fn render_to_cubic(&self, accuracy: f64) -> RenderResult<CubicBez> {
        let p1_3 = as_vec2(self.integrate(1. / 3.));
        let p2_3 = as_vec2(self.integrate(2. / 3.));
        let p1 = as_vec2(self.integrate(1.));
        let theta1 = self.theta(1.);

        let rc0 = p1_3.dot(p1) / p1.hypot();
        let rc1 = p2_3.dot(p1) / p1.hypot();

        let res = render_for_params_exact(
            p1_3,
            p2_3,
            p1,
            theta1,
            [rc0, rc1, 1. / 3., 2. / 3.],
            accuracy,
            accuracy.log2().abs().ceil() as usize + 2,
        );

        tracing::trace!(?res);

        res.map(|Render { params, .. }| {
            let [rc0, rc1, ..] = params.data.0[0];
            let c1 = rc1 * Vec2::from_angle(theta1);
            Affine::rotate(-p1.atan2()).then_scale(p1.hypot().recip())
                * CubicBez::new(
                    Point::ZERO,
                    Point::new(rc0, 0.),
                    (p1 - c1).to_point(),
                    p1.to_point(),
                )
        })
    }
}

fn system_for_rendering(
    p1_3_i: Vec2,
    p2_3_i: Vec2,
    p1: Vec2,
    theta1: f64,
) -> impl Fn([f64; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    move |[rc0, rc1, t1_3, t2_3]: [f64; 4]| -> ([f64; 4], [[f64; 4]; 4]) {
        let c1 = rc1 * Vec2::from_angle(theta1);
        let cb = CubicBez::new(
            Point::ZERO,
            Point::new(rc0, 0.),
            (p1 - c1).to_point(),
            p1.to_point(),
        );

        let p1_3 = cb.eval(t1_3);
        let p2_3 = cb.eval(t2_3);

        let dp1_3_dt1_3 = cb.deriv().eval(t1_3);

        let dp1_3_x_drc0 = 3. * (1. - t1_3).powi(2) * t1_3;
        let dp1_3_x_drc1 = -3. * (1. - t1_3) * t1_3.powi(2) * theta1.cos();
        let dp1_3_x_dt1_3 = dp1_3_dt1_3.x;
        let dp1_3_x_dt2_3 = 0.;

        let dp1_3_y_drc0 = 0.;
        let dp1_3_y_drc1 = -3. * (1. - t1_3) * t1_3.powi(2) * theta1.sin();
        let dp1_3_y_dt1_3 = dp1_3_dt1_3.y;
        let dp1_3_y_dt2_3 = 0.;

        let dp2_3_dt2_3 = cb.deriv().eval(t2_3);

        let dp2_3_x_drc0 = 3. * (1. - t2_3).powi(2) * t2_3;
        let dp2_3_x_drc1 = -3. * (1. - t2_3) * t2_3.powi(2) * theta1.cos();
        let dp2_3_x_dt1_3 = 0.;
        let dp2_3_x_dt2_3 = dp2_3_dt2_3.x;

        let dp2_3_y_drc0 = 0.;
        let dp2_3_y_drc1 = -3. * (1. - t2_3) * t2_3.powi(2) * theta1.sin();
        let dp2_3_y_dt1_3 = 0.;
        let dp2_3_y_dt2_3 = dp2_3_dt2_3.y;

        let err = [
            p1_3.x - p1_3_i.x,
            p1_3.y - p1_3_i.y,
            p2_3.x - p2_3_i.x,
            p2_3.y - p2_3_i.y,
        ];

        let jac = [
            [dp1_3_x_drc0, dp1_3_x_drc1, dp1_3_x_dt1_3, dp1_3_x_dt2_3],
            [dp1_3_y_drc0, dp1_3_y_drc1, dp1_3_y_dt1_3, dp1_3_y_dt2_3],
            [dp2_3_x_drc0, dp2_3_x_drc1, dp2_3_x_dt1_3, dp2_3_x_dt2_3],
            [dp2_3_y_drc0, dp2_3_y_drc1, dp2_3_y_dt1_3, dp2_3_y_dt2_3],
        ];

        (err, jac)
    }
}

#[must_use]
fn render_iterate_once(
    p1_3: Vec2,
    p2_3: Vec2,
    p1: Vec2,
    theta1: f64,
    guess: [f64; 4],
) -> ([f64; 4], Option<[f64; 4]>) {
    let (err, jac) = system_for_rendering(p1_3, p2_3, p1, theta1)(guess);

    let guess_v = Vector4::from_data(ArrayStorage([guess]));
    let err_v = Vector4::from_data(ArrayStorage([err]));
    let mut jac_v = Matrix4::from_data(ArrayStorage(jac)).transpose();
    let new_guess_v = (err_v.norm_squared().is_finite() && jac_v.try_inverse_mut())
        .then(|| guess_v - jac_v * err_v);

    (err, new_guess_v.map(|v| v.data.0[0]))
}

#[derive(Debug, Clone, Copy)]
pub struct Render {
    pub params: Vector4<f64>,
    pub err: Vector4<f64>,
    pub iter: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum RenderError {
    Singularity {
        guess: Vector4<f64>,
        err: Vector4<f64>,
        iter: usize,
    },
    OutOfIteration {
        guess: Vector4<f64>,
        err: Vector4<f64>,
    },
}

pub type RenderResult<T = Render> = Result<T, RenderError>;

pub fn render_for_params_exact(
    p1_3: Vec2,
    p2_3: Vec2,
    p1: Vec2,
    theta1: f64,
    mut guess: [f64; 4],
    threshold: f64,
    n_iter: usize,
) -> RenderResult {
    if p1_3.y == 0. && p2_3.y == 0. && p1.y == 0. {
        return Ok(Render {
            params: Vector4::new(1. / 3., 2. / 3., 1. / 3., 2. / 3.),
            err: Vector4::zeros(),
            iter: 0,
        });
    }

    // let p1_angle = -theta0;
    // let p0_5 = Affine::rotate(p1_angle) * p0_5;
    // let phi0_5 = phi0_5 + p1_angle;
    // let theta1 = theta1 + p1_angle;

    let mut err = [f64::INFINITY; 4];
    for i in 0..n_iter {
        let (new_err, new_guess) = render_iterate_once(p1_3, p2_3, p1, theta1, guess);
        let Some(new_guess) = new_guess else {
            return Err(RenderError::Singularity {
                guess: Vector4::from_data(ArrayStorage([guess])),
                err: Vector4::from_data(ArrayStorage([new_err])),
                iter: i,
            });
        };
        if new_err.iter().all(|e| e.abs() <= threshold) {
            return Ok(Render {
                params: Vector4::from_data(ArrayStorage([new_guess])),
                err: Vector4::from_data(ArrayStorage([new_err])),
                iter: i,
            });
        }
        guess = new_guess;
        err = new_err;
    }
    Err(RenderError::OutOfIteration {
        guess: Vector4::from_data(ArrayStorage([guess])),
        err: Vector4::from_data(ArrayStorage([err])),
    })
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test1() {
        let hb = HyperbezParams::new(0., -1., -1., 1., 1.);
        let res = hb.render_to_cubic(1e-3);
        tracing::trace!(?res);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2() {
        let hb = HyperbezParams::new(0., -2., -1., 1., 1.);
        let res = hb.render_to_cubic(1e-3);
        tracing::trace!(?res);
    }

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test2_2() {
        let hb = HyperbezParams::new(0., -2., -1., 1., 1.);
        let path = hb.render(1e-3).to_svg();
        tracing::trace!(%path);
    }
}
