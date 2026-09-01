use std::f64;

use xilem_web::svg::kurbo;

use crate::utils::*;

use super::*;

/// Given the denominator quadratic (c,d), solve (a,b) so the curve interpolates the
/// endpoint tangents (G1). theta is linear in (a,b), so the net turn
/// theta(1) = a*F(1) + b*G(1) = th1 - th0 eliminates b exactly -- and
/// pins the winding class, excluding stray full loops. What remains is a
/// scalar Newton on a: the residual is the start-tangent error, read off
/// the argument of the position integral I1.
///
/// Init 1 uses the ghost cubic's midpoint tangent as an estimate of
/// theta(1/2); with the net-turn condition that is a closed-form 2x2
/// (again by linearity). Init 2 is the small-angle/Euler formula. A
/// converged solution whose interior theta excursion leaves the principal
/// band (a hidden +-2pi round trip) is rejected and the next init tried.
fn fit_ab(c: f64, d: f64, l0: Vec2, l1: Vec2) -> Result<[f64; 2], String> {
    let f = HyperbezParams::new(1., 0., c, d);
    let g = HyperbezParams::new(0., 1., c, d);
    let f_1 = f.theta(1.);
    let g_1 = g.theta(1.);

    let [th0, th1] = [l0, l1].map(Vec2::angle);
    let dth = th1 - th0;

    if g_1.abs() < 1e-30 {
        return Err("g(1) vanished".into());
    }

    let b = |a| (dth - a * f_1) / g_1;

    let residual = |a: f64| {
        num_dual::first_derivative(
            |a| {
                let b = (-a * f_1 + dth) / g_1;
                let p1 = HyperbezParams::new(a, b, c.into(), d.into()).integrate(1f64.into());
                if p1.x == 0. || p1.y == 0. {
                    f64::NAN.into()
                } else {
                    norm_radians(-p1.y.atan2(p1.x) - th0)
                }
            },
            a,
        )
    };

    let solve_step = |mut a| {
        let mut residue = 1.;
        let mut derivative;
        // if residue.is_nan() {
        //     return Err("integral vanished".to_string());
        // }

        'newt: for i in 0..30 {
            (residue, derivative) = residual(a);
            if residue.abs() < 1e-11 {
                tracing::trace!(i, "converged");
                break;
            }
            if derivative == 0. {
                residue = 1.;
                tracing::trace!(i, "diverged");
                break;
            }
            let newton_step = -residue / derivative;
            for j in 0..12 {
                let damping = 2f64.powi(-j);
                let new_a = a + damping * newton_step;
                let (new_residue, _) = residual(new_a);
                if !new_residue.is_nan() && new_residue.abs() < residue.abs() {
                    a = new_a;
                    residue = new_residue;
                    continue 'newt;
                }
            }
            tracing::trace!(i, "didn't make progress");
            break;
        }

        if residue.abs() >= 1e-8 {
            return Err(format!("residue: {residue}"));
        }

        let b = b(a);
        // reject interior-excursion loopers (net turn right, round trip inside)
        if a.abs() > 1e-14 {
            let tx = -b / a;
            if tx > 0. && tx < 1. {
                let ex = a * f.theta(tx) + b * g.theta(tx);
                if !(dth.min(0.) - f64::consts::PI - 0.3 < ex
                    && ex < dth.max(0.) + f64::consts::PI + 0.3)
                {
                    return Err("loop".into());
                }
            }
        }
        Ok([a, b])
    };

    // init ladder: midpoint-tangent informed (closed-form 2x2), then Euler
    let dr = Vec2::new(2., 0.) - l0.normalize() - l1.normalize();
    let adj = if dr == Vec2::ZERO {
        0.
    } else {
        norm_radians(-dr.atan2() - th0 - 0.5 * dth)
    };
    let th_h = 0.5 * dth + adj;
    let f_h = f.theta(0.5);
    let g_h = g.theta(0.5);
    let det = f_h * g_1 - g_h * f_1;

    if det.abs() > 1e-300 {
        solve_step((th_h * g_1 - g_h * dth) / det)
    } else {
        Err("determinant vanished".into())
    }
}

/// The control-polygon -> hyperbezier map.
///
///   u, v : arms normalized by l*. Written as l*1.5*(1+cos th) = l/l*
///          with a guard so u -> 0 continuously as a tangent approaches
///          +-180 deg. That fade is what makes rotating a control point
///          through the back direction continuous: the +-2pi winding
///          ambiguity collapses into a sub-pixel curl at the endpoint.
///
///   ln q = K*ln u + CK*max(ln u, 0)^3, then sc(), a smooth clamp
///          (p-norm min with radius L): 5.75..8 above (corner headroom),
///          15 below (deep enough that the residual seam at +-180 deg is
///          ~3e-5 of the chord). Note sc(0) = 0, so u = v = 1 maps to
///          q0 = q1 = 1: an exact Euler spiral or circle at every angle.
///
///   qm   : this representation fixes the gauge freedom by normalizing
///          integral_0^1 q^(-3/2) ds = 1 (rather than q(0) = 1), with
///          (q0, qm, q1) the Bernstein weights of q. The antiderivative
///          of q^(-3/2) is elementary and the condition telescopes to a
///          closed form for qm.
///
///   Then convert to the q(0) = 1 gauge used above and solve the angles.
pub fn solve(
    cb: kurbo::CubicBez,
    param_u0: f64,
    param_eps: f64,
    param_l: f64,
    param_k: f64,
    param_c: f64,
) -> Result<HyperbezParams<f64>, String> {
    let l0 = cb.p1.to_vec2();
    let l1 = cb.p3 - cb.p2;

    let mut err = String::new();
    for i in 0..14 {
        let loosening = 0.98.powi(i);

        let [q0, q1] = [l0, l1].map(|l| {
            let u = param_u0 + l.length() * param_c * (1. + param_eps.exp2() + l.angle().cos());
            let lg_u = u.log2();
            let x = param_k * lg_u * loosening;
            let lg_q = x / (1. + (x / param_l).powi(4)).powf(0.25);
            lg_q.exp2()
        });

        let qm = 1. / q0.sqrt() + 1. / q1.sqrt() - (q0 * q1).sqrt();

        let c = (q0 - 2. * qm + q1) / q0;
        let d = 2. * (qm - q0) / q0;

        match fit_ab(c, d, l0, l1) {
            Ok([a, b]) => {
                tracing::trace!(i, "fitted");
                return Ok(HyperbezParams::new(a, b, c, d));
            }
            Err(e) => err = e,
        }
    }

    Err(err)
}
