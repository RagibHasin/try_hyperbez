use std::{ops::Range, rc::Rc};

use wasm_bindgen::JsCast;
use xilem_web::{
    AnyDomView, DomView,
    core::View,
    elements::{
        html::{self, div, option, select},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, Point},
};

use hyperbez_toy::{utils::parse_param, *};

use crate::common::*;
use crate::components::*;

#[derive(Default)]
pub(crate) struct AppState {
    pub data: AppData,
    memo: Memoized<AppData, MemoizedState>,

    plots: plots::State,
    sheet: sheet::State,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub(crate) struct AppData {
    a: f64,
    b: f64,
    c: f64,
    d_m: f64,

    old_c: f64,

    is_d: bool,
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

impl From<AppData> for AppState {
    fn from(data: AppData) -> Self {
        AppState {
            data,
            ..Default::default()
        }
    }
}

impl std::str::FromStr for AppData {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let params = &mut s.split(",");
        let a = parse_param(params, "a")?;
        let b = parse_param(params, "b")?;
        let c = parse_param(params, "c")?;
        let d = parse_param(params, "d")?;
        let render_method = match params.next().ok_or("not enough params: render_option")? {
            "unopt" => RenderMethod::UnoptimizedCurveFit,
            "opt" => RenderMethod::OptimizedCurveFit,
            "subdiv" => RenderMethod::SubdivisionSolve,
            _ => unreachable!(),
        };
        let accuracy_order = parse_param(params, "accuracy_order")?;
        Ok(AppData {
            a,
            b,
            c,
            d_m: d,
            old_c: c,
            is_d: true,
            render_method,
            accuracy_order,
        })
    }
}

impl std::fmt::Display for AppData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let AppData {
            a,
            b,
            c,
            d_m,
            is_d,
            render_method,
            accuracy_order,
            ..
        } = *self;
        let render = match render_method {
            RenderMethod::UnoptimizedCurveFit => "unopt",
            RenderMethod::OptimizedCurveFit => "opt",
            RenderMethod::SubdivisionSolve => "subdiv",
        };
        let d = d_m * if is_d { 1. } else { c };
        write!(f, "{a},{b},{c},{d},{render},{accuracy_order}",)
    }
}

impl Default for AppData {
    fn default() -> Self {
        Self {
            a: 0.,
            b: -1.,
            c: -1.,
            d_m: 1.,
            old_c: -1.,
            is_d: true,
            render_method: RenderMethod::UnoptimizedCurveFit,
            accuracy_order: 1.,
        }
    }
}

impl AppData {
    fn toggle_is_d(&mut self) {
        self.d_m *= if self.is_d { self.c.recip() } else { self.c };
        self.is_d = !self.is_d;
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

    let d = data.d_m * if data.is_d { 1. } else { data.c };
    let hyperbez = hb::Hyperbezier::from_points_params(
        hb::HyperbezParams::new(data.a, data.b, data.c, data.d_m),
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

    let points = path
        .elements()
        .iter()
        .filter_map(|e| e.end_point())
        .map(|p| Circle::new(p, NODE_RADIUS))
        .collect::<Vec<_>>();
    let n_points = points.len();

    let frag_path = path.id("hyperbez");
    let frag_points = g(points).id("nodes");

    let frag_m_d = if data.is_d {
        labeled_valued("m: ", (), format!("{:.3}", data.d_m / data.c))
    } else {
        labeled_valued("d: ", (), format!("{:.3}", d))
    };
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
        frag_m_d,
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

    let frag_a = labeled_slider!("a: ", AppData::data.a, -20., 20., 0.2);
    let frag_b = labeled_slider!("b: ", AppData::data.b, -20., 20., 0.2);
    let frag_c = labeled_slider!("c: ", AppData::data.c, -10., 10., 0.1);
    let frag_d_m = labeled_slider!(
        html::span(if data.is_d { "d: " } else { "m: " })
            .on_dblclick(|data: &mut AppData, _| data.toggle_is_d()),
        AppData::data.d_m,
        if data.is_d { d_limit.start } else { -3. },
        if data.is_d { d_limit.end } else { -0.1 },
        0.1,
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
    .on_change(|data: &mut AppData, e| {
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
    let frag_accuracy = labeled_slider!("log(α): ", AppData::data.accuracy_order, 1., 4., 1.);

    let frag_options = div((
        frag_a,
        frag_b,
        frag_c,
        frag_d_m,
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

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> + use<> {
    if state.data.c == 0. {
        state.data.c = state.data.old_c;
    } else {
        state.data.old_c = state.data.c;
    };
    let memo = state.memo.update(state.data, memoized_app_logic);

    app_view!(state, memo, (memo.frag_results.clone()), (), |_, _| ())
}
