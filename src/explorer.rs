use std::{ops::Range, rc::Rc};

use wasm_bindgen::JsCast;
use xilem_web::{
    elements::{
        html::{self, div, option, select},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, ParamCurve, Point, Vec2},
    AnyDomView, DomView,
};

use hyperbez_toy::*;

use crate::components::*;

#[derive(Default)]
pub(crate) struct AppState {
    data: AppData,
    memo: Memoized<AppData, MemoizedState>,

    plots: plots::State,
    sheet: sheet::State,
}

#[derive(Debug, PartialEq, Clone, Copy)]
struct AppData {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,

    render_method: RenderMethod,
    accuracy_order: f64,
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum RenderMethod {
    UnoptimizedCurveFit,
    OptimizedCurveFit,
    SubdivisionSolve,
}
struct MemoizedState {
    hyperbez: hb::Hyperbezier,
    theta: Rc<[f64]>,
    kappa: Rc<[f64]>,

    frag_path: Rc<AnyDomView<sheet::State, sheet::DragAction>>,
    frag_points: Rc<AnyDomView<sheet::State, sheet::DragAction>>,
    frag_results: Rc<AnyDomView<AppState>>,
    frag_options: Rc<AnyDomView<AppData>>,
}

const BASE_WIDTH: f64 = 500.;

impl Default for AppData {
    fn default() -> Self {
        Self {
            a: 0.,
            b: -1.,
            c: -1.,
            d: 1.,
            e: 1.,
            render_method: RenderMethod::UnoptimizedCurveFit,
            accuracy_order: 1.,
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

fn memoized_app_logic(data: &AppData) -> MemoizedState {
    let d_limit = d_limit_rounded(data.a, data.b, data.c);
    // state.d = state.d.clamp(d_limit.start, d_limit.end);

    let hyperbez = hb::Hyperbezier::from_points_params(
        hb::HyperbezParams::new(data.a, data.b, data.c, data.d, data.e),
        Point::ZERO,
        Point::new(BASE_WIDTH, 0.),
    );

    // {
    //     let arg_uv = hyperbez.params().integrate(1.).angle().to_degrees();
    //     let int_theta = hb::integrate_theta(hyperbez.params(), 1.).to_degrees();
    //     tracing::trace!(?arg_uv, ?int_theta);
    // }

    let accuracy = 10f64.powf(-data.accuracy_order);
    let path = Affine::FLIP_Y
        * match data.render_method {
            RenderMethod::UnoptimizedCurveFit => kurbo::fit_to_bezpath(&hyperbez, accuracy),
            RenderMethod::OptimizedCurveFit => kurbo::fit_to_bezpath_opt(&hyperbez, accuracy),
            RenderMethod::SubdivisionSolve => hyperbez.render(accuracy),
        };
    let arclen = hyperbez.scale_rot().length() / BASE_WIDTH;
    let plot_points = (0..=1000).map(|i| i as f64 * 1e-3);
    let theta = plot_points
        .clone()
        .map(|t| hyperbez.theta(t).to_degrees())
        .collect::<Rc<_>>();
    let kappa = plot_points
        .map(|t| hyperbez.kappa(t) * hyperbez.scale_rot().length())
        .collect::<Rc<_>>();
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
    let n_points = points.len();

    let frag_path = path.id("hyperbez");
    let frag_points = g(points).id("nodes");

    let frag_arclen = labeled_valued("S / b: ", (), format!("{:.3}", arclen));
    let frag_theta0 = labeled_valued("θ₀: ", (), format!("{:.1}°", theta0));
    let frag_theta1 = labeled_valued("θ₁: ", (), format!("{:.1}°", theta1));
    let frag_kappa0 = labeled_valued("κ₀: ", (), format!("{:.3}", kappa0));
    let frag_kappa1 = labeled_valued("κ₁: ", (), format!("{:.3}", kappa1));

    let [k0, k1] = hyperbez.params().endk();
    let frag_k0 = labeled_valued("k₀: ", (), format!("{k0:.3}"));
    let frag_k1 = labeled_valued("k₁: ", (), format!("{k1:.3}"));

    let frag_n_points = labeled_valued("n: ", (), n_points);

    let frag_results = div((
        frag_arclen,
        spacer(),
        frag_theta0,
        frag_theta1,
        spacer(),
        frag_kappa0,
        frag_kappa1,
        spacer(),
        spacer(),
        frag_k0,
        frag_k1,
        spacer(),
        frag_n_points,
    ))
    .class("results");

    let frag_a = labeled_valued(
        "a: ",
        slider(data.a, -20., 20., 0.2).map_state(move |data: &mut AppData| &mut data.a),
        textbox(data.a).map_state(move |data: &mut AppData| &mut data.a),
    );
    let frag_b = labeled_valued(
        "b: ",
        slider(data.b, -20., 20., 0.2).map_state(move |data: &mut AppData| &mut data.b),
        textbox(data.b).map_state(move |data: &mut AppData| &mut data.b),
    );
    let frag_c = labeled_valued(
        "c: ",
        slider(data.c, 0.1, 10., 0.1).adapt(move |data: &mut AppData, thunk| {
            let mut temp_c = data.c;
            let r = thunk.call(&mut temp_c);
            data.c = if temp_c == 0. {
                if data.c > 0. {
                    -0.1
                } else {
                    0.1
                }
            } else {
                temp_c
            };
            r
        }),
        textbox(data.c).map_state(move |data: &mut AppData| &mut data.c),
    );
    let frag_d = labeled_valued(
        "d: ",
        slider(data.d, d_limit.start, d_limit.end, 0.1)
            .map_state(move |data: &mut AppData| &mut data.d),
        textbox(data.d).map_state(move |data: &mut AppData| &mut data.d),
    );
    let frag_e = labeled_valued(
        "e: ",
        slider(data.e, 0.1, 20., 0.1).map_state(move |data: &mut AppData| &mut data.e),
        textbox(data.e).map_state(move |data: &mut AppData| &mut data.e),
    );

    let frag_render_method = select((
        option("Unoptimized Curve Fitting")
            .value("UnoptimizedCurveFit")
            .selected(matches!(
                data.render_method,
                RenderMethod::UnoptimizedCurveFit
            )),
        option("Optimized Curve Fitting")
            .value("OptimizedCurveFit")
            .selected(matches!(
                data.render_method,
                RenderMethod::OptimizedCurveFit
            )),
        option("Subdivision Solve")
            .value("SubdivisionSolve")
            .selected(matches!(data.render_method, RenderMethod::SubdivisionSolve)),
    ))
    .on_change(move |data: &mut AppData, e| {
        match e
            .target()
            .unwrap()
            .unchecked_into::<web_sys::HtmlSelectElement>()
            .value()
            .as_ref()
        {
            "UnoptimizedCurveFit" => data.render_method = RenderMethod::UnoptimizedCurveFit,
            "OptimizedCurveFit" => data.render_method = RenderMethod::OptimizedCurveFit,
            "SubdivisionSolve" => data.render_method = RenderMethod::SubdivisionSolve,
            _ => {}
        }
    });
    let frag_accuracy = labeled_valued(
        "log(α): ",
        slider(data.accuracy_order, 1., 4., 1.)
            .map_state(move |data: &mut AppData| &mut data.accuracy_order),
        data.accuracy_order,
    );

    let frag_options = div((
        frag_a,
        frag_b,
        frag_c,
        frag_d,
        frag_e,
        spacer(),
        frag_render_method,
        frag_accuracy,
    ))
    .id("options");

    MemoizedState {
        hyperbez,
        theta,
        kappa,
        frag_results: Rc::new(frag_results),
        frag_path: Rc::new(frag_path),
        frag_points: Rc::new(frag_points),
        frag_options: Rc::new(frag_options),
    }
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let MemoizedState {
        hyperbez,
        theta,
        kappa,
        frag_path,
        frag_points,
        frag_results,
        frag_options,
    } = state.memo.update(state.data, memoized_app_logic);

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
        frag_results.clone(),
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
        .view(theta, kappa)
        .map_state(|state: &mut AppState| &mut state.plots);

    let frag_svg = state
        .sheet
        .view((frag_path.clone(), hover_mark, frag_points.clone()))
        .adapt(move |state: &mut AppState, thunk| thunk.call(&mut state.sheet).map(|_| {}));

    div((
        div((
            div((
                frag_options
                    .clone()
                    .map_state(|state: &mut AppState| &mut state.data),
                frag_results,
            ))
            .id("ui"),
            frag_plots,
        ))
        .id("pane-left"),
        div(frag_svg).id("render-sheet"),
    ))
    .id("app-root")
}
