use std::f64;

use xilem_web::svg::kurbo;

use crate::utils::*;

use super::*;

/// Cut list for piecewise integration: curvature extrema, plus a geometric
/// ladder of cuts around a deep minimum of q (near-corner), where the
/// q^(-3/2) factor makes theta vary on a scale ~ w0 = sqrt(q_min/c).
fn cut_points(a: f64, b: f64, c: f64, d: f64, t: f64) -> Vec<f64> {
    let mut cuts = HyperbezParams::new(a, b, c, d)
        .kappa_extrema()
        .into_iter()
        .filter(|w| *w < t)
        .collect::<Vec<_>>();

    if c > 1e-12 {
        let tv = -0.5 * d / c;
        if tv > 0. && tv < t {
            let qv = (c * tv + d) * tv + 1.;
            if qv < 0.05 {
                let w0 = (qv.max(1e-300) / c).sqrt();
                cuts.push(tv);
                let mut k = 3.;
                while k * w0 < 0.5 {
                    if tv - k * w0 > 0. && tv - k * w0 < t {
                        cuts.push(tv - k * w0);
                    }
                    if tv + k * w0 > 0. && tv + k * w0 < t {
                        cuts.push(tv + k * w0);
                    }
                    k *= 3.;
                }
            }
        }
    }

    cuts.extend([0., t]);
    cuts.sort_by(f64::total_cmp);
    cuts
}

/// Position integral z(t) = integral_0^t e^(i*theta) ds. GL-32 per cut
/// interval; intervals not starting at 0 are re-expressed via subseg so
/// each gets a fresh, well-conditioned chart.
fn integrate(a: f64, b: f64, c: f64, d: f64, t: f64, theta: impl Fn(f64) -> f64) -> Vec2 {
    //   theta = theta || makeTheta(a, b, c, d);
    let pts = cut_points(a, b, c, d, t);
    let mut z = Vec2::ZERO;

    for [w0, w1] in pts.array_windows().copied() {
        if w1 - w0 < 1e-15 {
            continue;
        }
        if w0 == 0. {
            let mut s = Vec2::ZERO;
            for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
                let u = 0.5 * w1 * (*xi + 1.);
                let th = theta(u);
                s += *wi * Vec2::from_angle(th);
            }
            z += 0.5 * w1 * s;
        } else {
            let sp = HyperbezParams::new(a, b, c, d).subsegment(w0..w1);
            let off = theta(w0);
            let mut s = Vec2::ZERO;
            for (wi, xi) in GAUSS_LEGENDRE_COEFFS_32 {
                let u = 0.5 * w1 * (*xi + 1.);
                let th = off + sp.theta(u);
                s += *wi * Vec2::from_angle(th);
            }
            z += 0.5 * (w1 - w0) * s;
        }
    }
    z
}

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
fn fit_ab(c: f64, d: f64, l0: Vec2, l1: Vec2) -> (Option<[f64; 2]>, String) {
    let f = HyperbezParams::new(1., 0., c, d);
    let g = HyperbezParams::new(0., 1., c, d);
    let f_1 = f.theta(1.);
    let g_1 = g.theta(1.);

    let [th0, th1] = [l0, l1].map(Vec2::angle);
    let dth = th1 - th0;

    if g_1.abs() < 1e-30 {
        return (None, "g(1) vanished".into());
    }

    let b = |a| (dth - a * f_1) / g_1;

    let residual = |a| {
        let b = b(a);
        let theta = |t| a * f.theta(t) + b * g.theta(t);
        let integral = integrate(a, b, c, d, 1., theta);
        (integral.hypot() >= 1e-12).then(|| norm_radians(-integral.atan2() - th0))
    };

    let solve_step = |mut a| {
        let mut residue = residual(a).ok_or("integral vanished")?;

        for _ in 0..30 {
            if residue.abs() < 1e-11 {
                break;
            }
            let h = 1e-7 * a.abs().max(1.);
            let Some(residue_next) = residual(a + h) else {
                residue = 1.;
                break;
            };
            let dr = (residue_next - residue) / h;
            if dr == 0. {
                residue = 1.;
                break;
            }
            let step = -residue / dr;
            let mut lam = 1.;
            let mut ok = false;
            for _ in 0..12 {
                if let Some(rn) = residual(a + lam * step)
                    && rn.abs() < residue.abs()
                {
                    a += lam * step;
                    residue = rn;
                    ok = true;
                    break;
                }
                lam *= 0.5;
            }
            if !ok {
                break;
            }
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
    let raw = if dr.hypot() > 1e-9 {
        -dr.atan2()
    } else {
        th0 + 0.5 * dth
    };
    let estg = 0.5 * dth + norm_radians(raw - th0 - 0.5 * dth);
    let f_h = f.theta(0.5);
    let g_h = g.theta(0.5);
    let det = f_h * g_1 - g_h * f_1;

    if det.abs() > 1e-10 * (f_h * g_1).abs().max((g_h * f_1).abs()).max(1e-300) {
        match solve_step((estg * g_1 - g_h * dth) / det) {
            Ok(ab) => (Some(ab), String::new()),
            Err(e) => match solve_step(6. * (th0 + th1)) {
                Ok(ab) => (Some(ab), format!("midpoint-tangent failed: {e}")),
                Err(e2) => (
                    None,
                    format!("midpoint-tangent failed: {e}, euler approx failed: {e2}"),
                ),
            },
        }
    } else {
        (None, "determinant vanished".into())
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
) -> (Option<HyperbezParams<f64>>, String) {
    let l0 = cb.p1.to_vec2();
    let l1 = cb.p3 - cb.p2;

    let [q0, q1] = [l0, l1].map(|l| {
        let u = param_u0 + l.length() * 1.5 * (1. + param_eps.exp2() + l.angle().cos());
        let ln_u = u.ln();
        let x = param_k * ln_u + param_c * ln_u.max(0.).powi(3);
        let ln_q = x / (1. + (x / param_l).powi(4)).powf(0.25);
        ln_q.exp()
    });

    let qm = 1. / q0.sqrt() + 1. / q1.sqrt() - (q0 * q1).sqrt();

    let c = (q0 - 2. * qm + q1) / q0;
    let d = 2. * (qm - q0) / q0;

    let (ab, comment) = fit_ab(c, d, l0, l1);

    (ab.map(|[a, b]| HyperbezParams::new(a, b, c, d)), comment)
}
