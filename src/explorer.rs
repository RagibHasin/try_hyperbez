use std::ops::Range;

use xilem_web::{
    elements::{
        html::{self, button, div},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, ParamCurve, Point, Vec2},
    DomView,
};

use hyperbez_toy::*;

use crate::components::*;

pub(crate) struct AppState {
    a: f64,
    b: f64,
    c: f64,
    d: f64,

    optimize: bool,
    accuracy_order: f64,

    plots: plots::State,
    sheet: sheet::State,
}

const BASE_WIDTH: f64 = 500.;

impl Default for AppState {
    fn default() -> Self {
        Self {
            a: 0.,
            b: -1.,
            c: -1.,
            d: 1.,
            optimize: false,
            accuracy_order: 1.,
            plots: Default::default(),
            sheet: Default::default(),
        }
    }
}

pub fn d_limit(a: f64, b: f64, c: f64) -> Range<f64> {
    if c == 0. {
        return -1.0..10.;
    }
    let end = 2. * c.abs().sqrt();
    let start = -end / (b * c - a).abs().log10().max(4. / 3.);
    start..end
}

pub fn d_limit_rounded(a: f64, b: f64, c: f64) -> Range<f64> {
    let d_limit = d_limit(a, b, c);
    ((d_limit.start * 10.).ceil() / 10.)..(((d_limit.end - f64::EPSILON) * 10.).floor() / 10.)
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let d_limit = d_limit_rounded(state.a, state.b, state.c);
    // state.d = state.d.clamp(d_limit.start, d_limit.end);

    let hyperbez = hb::Hyperbezier::from_points_params(
        hb::HyperbezParams::new(state.a, state.b, state.c, state.d, 1.),
        Point::ZERO,
        Point::new(BASE_WIDTH, 0.),
    );

    // {
    //     let arg_uv = hyperbez.params().integrate(1.).angle().to_degrees();
    //     let int_theta = hb::integrate_theta(hyperbez.params(), 1.).to_degrees();
    //     tracing::trace!(?arg_uv, ?int_theta);
    // }

    let accuracy = 10f64.powf(-state.accuracy_order);
    let path = Affine::FLIP_Y
        * if state.optimize {
            kurbo::fit_to_bezpath_opt(&hyperbez, accuracy)
        } else {
            kurbo::fit_to_bezpath(&hyperbez, accuracy)
        };
    let arclen = hyperbez.scale_rot().length() / BASE_WIDTH;
    let (theta, kappa): (Vec<_>, Vec<_>) = (0..=1000)
        .map(|i| i as f64 * 1e-3)
        .map(|t| {
            (
                hyperbez.theta(t).to_degrees(),
                hyperbez.kappa(t) * hyperbez.scale_rot().length(),
            )
        })
        .unzip();
    let theta0 = *theta.first().unwrap();
    let theta1 = *theta.last().unwrap();
    let kappa0 = *kappa.first().unwrap();
    let kappa1 = *kappa.last().unwrap();

    const NODE_RADIUS: f64 = 5.;
    let points = path
        .elements()
        .iter()
        .filter_map(|e| e.end_point())
        .map(|p| Circle::new(p, NODE_RADIUS))
        .collect::<Vec<_>>();

    let mut hovered_point = None;
    let mut hovered_theta = None;
    let mut hovered_kappa = None;
    let mut hover_mark = None;
    if let Some(s) = state.plots.hovered_x() {
        hovered_point = Some(hyperbez.eval(s));
        let i = (s * 1e3) as usize;
        (hovered_theta, hovered_kappa) = (Some(theta[i]), Some(kappa[i]));
        let (theta, kappa) = (theta[i].to_radians(), kappa[i]);

        let p = hyperbez.eval(s);
        let r_curv = hyperbez.scale_rot().length() / kappa;
        let tangent_half = 0.25 * BASE_WIDTH * Vec2::from_angle(theta);
        let tangent = Affine::FLIP_Y * Line::new(p - tangent_half, p + tangent_half);
        let r_curv = Affine::FLIP_Y
            * Line::new(
                p,
                p + r_curv * Vec2::from_angle(theta + std::f64::consts::FRAC_PI_2),
            );

        hover_mark = Some(
            g((
                tangent.id("tangent"),
                r_curv.id("curvature"),
                Circle::new(Affine::FLIP_Y * p, 3.),
            ))
            .class("hover"),
        );
    };

    let frag_arclen = labeled_valued("S / b: ", (), format!("{:.3}", arclen));
    let frag_theta0 = labeled_valued("θ₀: ", (), format!("{:.1}°", theta0));
    let frag_theta1 = labeled_valued("θ₁: ", (), format!("{:.1}°", theta1));
    let frag_kappa0 = labeled_valued("κ₀: ", (), format!("{:.3}", kappa0));
    let frag_kappa1 = labeled_valued("κ₁: ", (), format!("{:.3}", kappa1));
    let frag_n_points = labeled_valued("n: ", (), points.len());

    let empty = "~".to_string();
    let frag_hovered_s = labeled_valued(
        "s: ",
        (),
        state
            .plots
            .hovered_x()
            .map_or(empty.clone(), |s| format!("{:.3}", s)),
    );
    let frag_hovered_p_x = labeled_valued(
        ("P", html::sub("x"), "(s): "),
        (),
        hovered_point.map_or(empty.clone(), |v| format!("{:.2}", v.x)),
    );
    let frag_hovered_p_y = labeled_valued(
        ("P", html::sub("x"), "(s): "),
        (),
        hovered_point.map_or(empty.clone(), |v| format!("{:.2}", v.y)),
    );
    let frag_hovered_theta = labeled_valued(
        "θ(s): ",
        (),
        hovered_theta.map_or(empty.clone(), |v| format!("{:.1}°", v)),
    );
    let frag_hovered_kappa = labeled_valued(
        "κ(s): ",
        (),
        hovered_kappa.map_or(empty.clone(), |v| format!("{:.2}", v)),
    );

    let frag_results = (
        div((
            frag_arclen,
            spacer(),
            frag_theta0,
            frag_theta1,
            spacer(),
            frag_kappa0,
            frag_kappa1,
            spacer(),
            frag_n_points,
        ))
        .class("results"),
        div((
            frag_hovered_s,
            frag_hovered_p_x,
            frag_hovered_p_y,
            frag_hovered_theta,
            frag_hovered_kappa,
        ))
        .class("results"),
    );

    let frag_plots = state
        .plots
        .view(&theta, &kappa)
        .map_state(|state: &mut AppState| &mut state.plots);

    let frag_svg = state
        .sheet
        .view((path.id("hyperbez"), hover_mark, g(points).id("nodes")))
        .adapt(move |state: &mut AppState, thunk| thunk.call(&mut state.sheet).map(|_| {}));

    let frag_a = labeled_valued(
        "a: ",
        slider(state.a, -20., 20., 0.2).map_state(move |state: &mut AppState| &mut state.a),
        textbox(state.a).map_state(move |state: &mut AppState| &mut state.a),
    );
    let frag_b = labeled_valued(
        "b: ",
        slider(state.b, -20., 20., 0.2).map_state(move |state: &mut AppState| &mut state.b),
        textbox(state.b).map_state(move |state: &mut AppState| &mut state.b),
    );
    let frag_c = labeled_valued(
        "c: ",
        slider(state.c, 0.1, 10., 0.1).adapt(move |state: &mut AppState, thunk| {
            let mut temp_c = state.c;
            let r = thunk.call(&mut temp_c);
            state.c = if temp_c == 0. {
                if state.c > 0. {
                    -0.1
                } else {
                    0.1
                }
            } else {
                temp_c
            };
            r
        }),
        textbox(state.c).map_state(move |state: &mut AppState| &mut state.c),
    );
    let frag_d = labeled_valued(
        "d: ",
        slider(state.d, d_limit.start, d_limit.end, 0.1)
            .map_state(move |state: &mut AppState| &mut state.d),
        textbox(state.d).map_state(move |state: &mut AppState| &mut state.d),
    );

    let frag_optimize = button({
        let this = &state;
        if this.optimize {
            "Do not optimize"
        } else {
            "Optimize"
        }
    })
    .on_click(|state: &mut AppState, _| state.optimize = !state.optimize);
    let frag_accuracy = labeled_valued(
        "log(α): ",
        slider(state.accuracy_order, 1., 4., 1.)
            .map_state(move |state: &mut AppState| &mut state.accuracy_order),
        state.accuracy_order,
    );

    let frag_options = div((
        frag_a,
        frag_b,
        frag_c,
        frag_d,
        spacer(),
        frag_optimize,
        frag_accuracy,
    ))
    .id("options");

    div((
        div((div((frag_options, frag_results)).id("ui"), frag_plots)).id("pane-left"),
        div(frag_svg).id("render-sheet"),
    ))
    .id("app-root")
}
