use std::rc::Rc;

use xilem_web::{
    AnyDomView, DomView,
    core::View,
    elements::{
        html::{self, div},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, Point, Shape},
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
    p1: Point,
    p2: Point,

    param_u0: f64,
    param_eps: f64,
    param_l: f64,
    param_k: f64,
}

type SheetElement = AnyDomView<sheet::State, sheet::DragAction>;

struct MemoizedState {
    hyperbez: hb::Hyperbezier,
    theta: Rc<[f64]>,
    kappa: Rc<[f64]>,

    frag_controls: Rc<SheetElement>,
    frag_cubicbez: Rc<SheetElement>,
    frag_path: Rc<SheetElement>,
    frag_points: Rc<SheetElement>,
    frag_results_1: Rc<AnyDomView<AppState>>,
    frag_results_2: Rc<AnyDomView<AppState>>,
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
        let p1_x = parse_param(params, "p1_x")?;
        let p1_y = parse_param(params, "p1_y")?;
        let p2_x = parse_param(params, "p2_x")?;
        let p2_y = parse_param(params, "p2_y")?;
        let param_u0 = parse_param(params, "param_u0")?;
        let param_eps = parse_param(params, "param_eps")?;
        let param_l = parse_param(params, "param_l")?;
        let param_k = parse_param(params, "param_k")?;
        Ok(AppData {
            p1: Point::new(p1_x, p1_y),
            p2: Point::new(p2_x, p2_y),
            param_u0,
            param_eps,
            param_l,
            param_k,
        })
    }
}

impl std::fmt::Display for AppData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let AppData {
            p1,
            p2,
            param_u0,
            param_eps,
            param_l,
            param_k,
        } = *self;
        write!(
            f,
            "{},{},{},{},{param_u0},{param_eps},{param_l},{param_k}",
            p1.x, p1.y, p2.x, p2.y
        )
    }
}

impl Default for AppData {
    fn default() -> Self {
        Self {
            p1: Point::new(140., 140.),
            p2: Point::new(360., 140.),
            param_u0: 0.1,
            param_eps: -4.,
            param_l: 10.,
            param_k: 1.,
        }
    }
}

fn memoized_app_logic(data: &AppData) -> MemoizedState {
    let scale_down = kurbo::TranslateScale::scale(1. / BASE_WIDTH);
    let cubicbez = kurbo::CubicBez::new(P0, data.p1, data.p2, P3);

    let params = hb::solver::ladder_alt::solve(
        scale_down * cubicbez,
        data.param_u0,
        data.param_eps,
        data.param_l,
        data.param_k,
    )
    .unwrap_or_else(|comment| {
        tracing::trace!(comment);
        hb::HyperbezParams::from_control(scale_down * data.p1, scale_down * data.p2)
    });
    let hyperbez = hb::Hyperbezier::from_points_params(params, P0, P3);

    let cubicbez = Affine::FLIP_Y * cubicbez.into_path(0.);
    let accuracy = 0.1;
    let optimize = false;
    let path = Affine::FLIP_Y
        * if optimize {
            kurbo::fit_to_bezpath_opt(&hyperbez, accuracy)
        } else {
            kurbo::fit_to_bezpath(&hyperbez, accuracy)
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

    let control0 = Affine::FLIP_Y * data.p1;
    let control0 = (
        Line::new((0., 0.), control0),
        Circle::new(control0, NODE_RADIUS).on_mousedown(|state: &mut sheet::State, e| {
            state.set_drag_handle(Some(sheet::Handle::C0));
            e.stop_propagation();
        }),
    );

    let control1 = Affine::FLIP_Y * data.p2;
    let control1 = (
        Line::new((BASE_WIDTH, 0.), control1),
        Circle::new(control1, NODE_RADIUS).on_mousedown(|state: &mut sheet::State, e| {
            state.set_drag_handle(Some(sheet::Handle::C1));
            e.stop_propagation();
        }),
    );

    let frag_controls = g((control0, control1)).class("control");
    let frag_cubicbez = cubicbez.id("cubicbez");
    let frag_path = path.id("hyperbez");
    let frag_points = g(points).id("nodes");

    let frag_a = labeled_valued("a: ", (), format!("{:.3}", params.a()));
    let frag_b = labeled_valued("b: ", (), format!("{:.3}", params.b()));
    let frag_c = labeled_valued("c: ", (), format!("{:.3}", params.c()));
    let frag_d = labeled_valued("d: ", (), format!("{:.3}", params.d()));

    let frag_q1 = labeled_valued("q₁: ", (), format!("{:.3}", params.q(1.)));

    let frag_arclen = labeled_valued("S / b: ", (), format!("{:.3}", arclen));
    let frag_theta0 = labeled_valued("θ₀: ", (), format!("{:.1}°", theta0));
    let frag_theta1 = labeled_valued("θ₁: ", (), format!("{:.1}°", theta1));
    let frag_kappa0 = labeled_valued("κ₀: ", (), format!("{:.3}", kappa0));
    let frag_kappa1 = labeled_valued("κ₁: ", (), format!("{:.3}", kappa1));
    let frag_n_points = labeled_valued("n: ", (), n_points);

    let frag_results_1 = div((frag_a, frag_b, frag_c, frag_d, spacer(), frag_q1)).class("results");
    let frag_results_2 = div((
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
    .class("results");

    let frag_p1_x = labeled_textbox!(("P1", html::sub("x"), ": "), AppData::data.p1.x);
    let frag_p1_y = labeled_textbox!(("P1", html::sub("y"), ": "), AppData::data.p1.y);
    let frag_p2_x = labeled_textbox!(("P2", html::sub("x"), ": "), AppData::data.p2.x);
    let frag_p2_y = labeled_textbox!(("P2", html::sub("y"), ": "), AppData::data.p2.y);

    let frag_param_u0 = labeled_slider!("u₀: ", AppData::data.param_u0, 0., 1., 0.02);
    let frag_param_eps = labeled_slider!("ϵ: ", AppData::data.param_eps, -40., -1., 1.);
    let frag_param_l = labeled_slider!("λ: ", AppData::data.param_l, 5., 20., 1.);
    let frag_param_k = labeled_slider!("K: ", AppData::data.param_k, 1., 3., 0.02);

    let frag_options = div((
        frag_p1_x,
        frag_p1_y,
        frag_p2_x,
        frag_p2_y,
        frag_param_u0,
        frag_param_eps,
        frag_param_l,
        frag_param_k,
    ))
    .id("options");

    MemoizedState {
        hyperbez,
        theta,
        kappa,
        frag_results_1: Rc::new(frag_results_1),
        frag_results_2: Rc::new(frag_results_2),
        frag_controls: Rc::new(frag_controls),
        frag_cubicbez: Rc::new(frag_cubicbez),
        frag_path: Rc::new(frag_path),
        frag_points: Rc::new(frag_points),
        frag_options: Rc::new(frag_options),
    }
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> + use<> {
    let memo = state.memo.update(state.data, memoized_app_logic);

    app_view!(
        state,
        memo,
        (memo.frag_results_1.clone(), memo.frag_results_2.clone()),
        (memo.frag_controls.clone(), memo.frag_cubicbez.clone()),
        |state: &mut AppState, sheet::DragAction { handle, point }| {
            *match handle {
                sheet::Handle::C0 => &mut state.data.p1,
                sheet::Handle::C1 => &mut state.data.p2,
            } = point;
        }
    )
}
