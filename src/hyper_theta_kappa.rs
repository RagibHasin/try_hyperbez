use std::{f64, fmt::Write, rc::Rc};

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

use hyperbez_toy::{utils::parse_param, *};

use crate::components::*;

#[derive(Default)]
pub(crate) struct AppState {
    pub data: AppData,
    memo: Memoized<AppData, MemoizedState>,

    plots: plots::State,
    sheet: sheet::State<Handle>,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub(crate) struct AppData {
    theta0: f64,
    theta1: f64,
    kappa0: f64,
    kappa1: f64,
    loopy: bool,

    running: bool,
    arrangement: Arrangement,
}

#[derive(Debug, PartialEq, Clone, Copy)]
enum Arrangement {
    Symmetric,
    Antisymmetric,
    Free,
}

type SheetElement = AnyDomView<sheet::State<Handle>, sheet::DragAction<Handle>>;

struct MemoizedState {
    hyperbez: hb::Hyperbezier,
    theta: Rc<[f64]>,
    kappa: Rc<[f64]>,

    frag_controls: Rc<SheetElement>,
    frag_path: Rc<SheetElement>,
    frag_points: Rc<SheetElement>,
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

const BASE_WIDTH: f64 = 500.;

impl std::str::FromStr for AppData {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let params = &mut s.split(",");
        let theta0 = parse_param::<f64>(params, "theta0")?.to_radians();
        let kappa0 = parse_param::<f64>(params, "kappa0")?;
        let (arrangement, theta1, kappa1) = match params
            .next()
            .ok_or("not enough params: arrangement, theta1, kappa1")?
        {
            "sym" => (Arrangement::Symmetric, theta0, kappa0),
            "anti" => (Arrangement::Antisymmetric, -theta0, -kappa0),
            _ => (
                Arrangement::Free,
                parse_param::<f64>(params, "theta1")?.to_radians(),
                parse_param(params, "kappa1")?,
            ),
        };
        let loopy = match params.next().ok_or("not enough params: loopy")? {
            "T" => true,
            "F" => false,
            _ => unreachable!(),
        };
        Ok(AppData {
            theta0,
            theta1,
            kappa0,
            kappa1,
            loopy,
            running: true,
            arrangement,
        })
    }
}

impl std::fmt::Display for AppData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let AppData {
            theta0,
            theta1,
            kappa0,
            kappa1,
            loopy,
            arrangement,
            ..
        } = *self;
        let theta0 = theta0.to_degrees();
        let theta1 = theta1.to_degrees();
        write!(f, "{theta0},{kappa0},")?;
        match arrangement {
            Arrangement::Symmetric => f.write_str("sym,"),
            Arrangement::Antisymmetric => f.write_str("anti,"),
            Arrangement::Free => write!(f, "{theta1},{kappa1},"),
        }?;
        f.write_char(if loopy { 'T' } else { 'F' })
    }
}

const P0: Point = Point::ZERO;
const P3: Point = Point::new(BASE_WIDTH, 0.);

impl Default for AppData {
    fn default() -> Self {
        Self {
            theta0: f64::consts::FRAC_PI_4,
            theta1: f64::consts::FRAC_PI_4,
            kappa0: -f64::consts::FRAC_PI_2,
            kappa1: -f64::consts::FRAC_PI_2,
            loopy: false,
            running: true,
            arrangement: Arrangement::Symmetric,
        }
    }
}

impl AppData {
    fn maintain_arrangement(&mut self, handle: Handle) {
        match self.arrangement {
            Arrangement::Symmetric => match handle {
                Handle::Theta0 => {
                    self.theta1 = self.theta0;
                    self.kappa1 = self.kappa0;
                }
                Handle::Theta1 => {
                    self.theta0 = self.theta1;
                    self.kappa0 = self.kappa1;
                }
            },
            Arrangement::Antisymmetric => match handle {
                Handle::Theta0 => {
                    self.theta1 = -self.theta0;
                    self.kappa1 = -self.kappa0;
                }
                Handle::Theta1 => {
                    self.theta0 = -self.theta1;
                    self.kappa0 = -self.kappa1;
                }
            },
            Arrangement::Free => {}
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Handle {
    Theta0,
    Theta1,
}

fn memoized_app_logic(data: &AppData, memo: Option<&mut MemoizedState>) -> Option<MemoizedState> {
    if !data.running {
        let memo = memo?;
        memo.frag_controls = Rc::new(frag_controls(data));
        memo.frag_options = Rc::new(frag_options(data));
        return None;
    }

    let params = hb::solver::theta_kappa::make_hyperbez(
        data.theta0,
        data.theta1,
        data.kappa0,
        data.kappa1,
        data.loopy,
    );

    let hyperbez = hb::Hyperbezier::from_points_params(params, P0, P3);

    let path = Affine::FLIP_Y * kurbo::fit_to_bezpath(&hyperbez, 0.1);
    let arclen = hyperbez.scale_rot().length() / BASE_WIDTH;
    let plot_points = (0..=1000).map(|i| i as f64 * 1e-3);
    let theta = plot_points
        .clone()
        .map(|t| hyperbez.theta(t).to_degrees())
        .collect::<Rc<_>>();
    let kappa = plot_points
        .map(|t| hyperbez.kappa(t) * hyperbez.scale_rot().length())
        .collect::<Rc<_>>();

    const NODE_RADIUS: f64 = 5.;
    let points = path
        .elements()
        .iter()
        .filter_map(|e| e.end_point())
        .map(|p| Circle::new(p, NODE_RADIUS))
        .collect::<Vec<_>>();
    let n_points = points.len();

    fn frag_controls(
        data: &AppData,
    ) -> impl DomView<sheet::State<Handle>, sheet::DragAction<Handle>> {
        const CONTROL_LENGTH: f64 = 100.;
        let control0 = Affine::FLIP_Y * (CONTROL_LENGTH * Vec2::from_angle(data.theta0)).to_point();
        let control0 = (
            Line::new((0., 0.), control0),
            Circle::new(control0, NODE_RADIUS).on_mousedown(
                |state: &mut sheet::State<Handle>, e| {
                    state.set_drag(Some(Handle::Theta0));
                    e.stop_propagation();
                },
            ),
        );

        let control1 = Affine::FLIP_Y * (P3 - CONTROL_LENGTH * Vec2::from_angle(-data.theta1));
        let control1 = (
            Line::new(P3, control1),
            Circle::new(control1, NODE_RADIUS).on_mousedown(
                |state: &mut sheet::State<Handle>, e| {
                    state.set_drag(Some(Handle::Theta1));
                    e.stop_propagation();
                },
            ),
        );

        g((control0, control1)).class("control")
    }

    let frag_controls = frag_controls(data);
    let frag_path = path.id("hyperbez");
    let frag_points = g(points).id("nodes");

    let frag_a = labeled_valued("a: ", (), format!("{:.3}", params.a()));
    let frag_b = labeled_valued("b: ", (), format!("{:.3}", params.b()));
    let frag_c = labeled_valued("c: ", (), format!("{:.3}", params.c()));
    let frag_d = labeled_valued("d: ", (), format!("{:.3}", params.d()));
    let frag_e = labeled_valued("e: ", (), format!("{:.3}", params.e()));

    let frag_arclen = labeled_valued("S / b: ", (), format!("{:.3}", arclen));

    let [k0, k1] = hyperbez.params().endk();
    let frag_k0 = labeled_valued("k₀: ", (), format!("{k0:.3}"));
    let frag_k1 = labeled_valued("k₁: ", (), format!("{k1:.3}"));

    let frag_n_points = labeled_valued("n: ", (), n_points);

    let frag_results = div((
        frag_arclen,
        spacer(),
        frag_a,
        frag_b,
        frag_c,
        frag_d,
        frag_e,
        spacer(),
        frag_k0,
        frag_k1,
        spacer(),
        frag_n_points,
    ))
    .class("results");

    fn frag_options(data: &AppData) -> impl DomView<AppData> {
        let frag_theta0 = labeled_valued(
            "θ₀",
            div(()),
            textbox(data.theta0.to_degrees())
                .map_state(|data: &mut AppData| &mut data.theta0)
                .map_action(move |data: &mut AppData, _| {
                    data.theta0 = data.theta0.to_radians();
                    data.maintain_arrangement(Handle::Theta0);
                }),
        );
        let frag_theta1 = labeled_valued(
            "θ₁",
            div(()),
            textbox(data.theta1.to_degrees())
                .map_state(|data: &mut AppData| &mut data.theta1)
                .map_action(move |data: &mut AppData, _| {
                    data.theta1 = data.theta1.to_radians();
                    data.maintain_arrangement(Handle::Theta1);
                }),
        );
        let frag_kappa0 = labeled_valued(
            "κ₀",
            div(()),
            textbox(data.kappa0)
                .map_state(|data: &mut AppData| &mut data.kappa0)
                .map_action(move |data: &mut AppData, _| data.maintain_arrangement(Handle::Theta0)),
        );
        let frag_kappa1 = labeled_valued(
            "κ₁",
            div(()),
            textbox(data.kappa1)
                .map_state(|data: &mut AppData| &mut data.kappa1)
                .map_action(move |data: &mut AppData, _| data.maintain_arrangement(Handle::Theta1)),
        );

        let frag_loopy = html::button(if data.loopy {
            "Remove loop"
        } else {
            "Create loop"
        })
        .on_click(|data: &mut AppData, _| data.loopy = !data.loopy);

        let frag_arrangement = select((
            option("Symmetric")
                .value("Symmetric")
                .selected(matches!(data.arrangement, Arrangement::Symmetric)),
            option("Antisymmetric")
                .value("Antisymmetric")
                .selected(matches!(data.arrangement, Arrangement::Antisymmetric)),
            option("Free")
                .value("Free")
                .selected(matches!(data.arrangement, Arrangement::Free)),
        ))
        .on_change(move |data: &mut AppData, e| {
            match e
                .target()
                .unwrap()
                .unchecked_into::<web_sys::HtmlSelectElement>()
                .value()
                .as_ref()
            {
                "Symmetric" => data.arrangement = Arrangement::Symmetric,
                "Antisymmetric" => data.arrangement = Arrangement::Antisymmetric,
                "Free" => data.arrangement = Arrangement::Free,
                _ => {}
            }
            data.maintain_arrangement(Handle::Theta0);
        });

        let frag_running = html::button(if data.running { "Pause" } else { "Resume" })
            .on_click(|data: &mut AppData, _| data.running = !data.running);

        div((
            frag_theta0,
            frag_theta1,
            frag_kappa0,
            frag_kappa1,
            spacer(),
            frag_loopy,
            spacer(),
            frag_arrangement,
            frag_running,
        ))
        .id("options")
    }
    let frag_options = frag_options(data);

    Some(MemoizedState {
        hyperbez,
        theta,
        kappa,
        frag_results: Rc::new(frag_results),
        frag_controls: Rc::new(frag_controls),
        frag_path: Rc::new(frag_path),
        frag_points: Rc::new(frag_points),
        frag_options: Rc::new(frag_options),
    })
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let MemoizedState {
        hyperbez,
        theta,
        kappa,
        frag_controls,
        frag_path,
        frag_points,
        frag_results,
        frag_options,
    } = state.memo.try_update(state.data, memoized_app_logic);

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
        .view((
            frag_controls.clone(),
            frag_path.clone(),
            hover_mark,
            frag_points.clone(),
        ))
        .map_state(|state: &mut AppState| &mut state.sheet)
        .map_action(
            move |state: &mut AppState, sheet::DragAction { data, event }| {
                let p = Affine::FLIP_Y
                    * Affine::scale(state.sheet.zoom())
                        .then_translate(state.sheet.origin().to_vec2())
                    * Point::new(event.offset_x() as f64, event.offset_y() as f64);

                match data {
                    Handle::Theta0 => state.data.theta0 = p.to_vec2().angle(),
                    Handle::Theta1 => {
                        state.data.theta1 =
                            (Point::new(BASE_WIDTH, 0.) - Affine::FLIP_Y * p).angle()
                    }
                };
                state.data.maintain_arrangement(data);
            },
        );

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
