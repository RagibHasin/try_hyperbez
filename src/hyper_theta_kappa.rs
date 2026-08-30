use std::{f64, rc::Rc};

use wasm_bindgen::JsCast;
use xilem_web::{
    AnyDomView, DomView,
    core::View,
    elements::{
        html::{self, button, div, option, select},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, Point, Vec2},
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

type SheetElement = AnyDomView<sheet::State, sheet::DragAction>;

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
        f.write_str(if loopy { "T" } else { "F" })
    }
}

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
    fn maintain_arrangement(&mut self, handle: sheet::Handle) {
        match self.arrangement {
            Arrangement::Symmetric => match handle {
                sheet::Handle::C0 => {
                    self.theta1 = self.theta0;
                    self.kappa1 = self.kappa0;
                }
                sheet::Handle::C1 => {
                    self.theta0 = self.theta1;
                    self.kappa0 = self.kappa1;
                }
            },
            Arrangement::Antisymmetric => match handle {
                sheet::Handle::C0 => {
                    self.theta1 = -self.theta0;
                    self.kappa1 = -self.kappa0;
                }
                sheet::Handle::C1 => {
                    self.theta0 = -self.theta1;
                    self.kappa0 = -self.kappa1;
                }
            },
            Arrangement::Free => {}
        }
    }
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

    fn frag_controls(data: &AppData) -> impl DomView<sheet::State, sheet::DragAction> + use<> {
        const CONTROL_LENGTH: f64 = 100.;
        let control0 = Affine::FLIP_Y * (CONTROL_LENGTH * Vec2::from_angle(data.theta0)).to_point();
        let control0 = (
            Line::new((0., 0.), control0),
            Circle::new(control0, NODE_RADIUS).on_mousedown(|state: &mut sheet::State, e| {
                state.set_drag_handle(Some(sheet::Handle::C0));
                e.stop_propagation();
            }),
        );

        let control1 = Affine::FLIP_Y * (P3 - CONTROL_LENGTH * Vec2::from_angle(-data.theta1));
        let control1 = (
            Line::new(P3, control1),
            Circle::new(control1, NODE_RADIUS).on_mousedown(|state: &mut sheet::State, e| {
                state.set_drag_handle(Some(sheet::Handle::C1));
                e.stop_propagation();
            }),
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
        spacer(),
        frag_k0,
        frag_k1,
        spacer(),
        frag_n_points,
    ))
    .class("results");

    fn frag_options(data: &AppData) -> impl DomView<AppData> + use<> {
        let frag_theta0 = labeled_valued(
            "θ₀",
            div(()),
            textbox(data.theta0.to_degrees())
                .map_state(|data: &mut AppData| &mut data.theta0)
                .map_message_result(|data: &mut AppData, r| {
                    data.theta0 = data.theta0.to_radians();
                    data.maintain_arrangement(sheet::Handle::C0);
                    r
                }),
        );
        let frag_theta1 = labeled_valued(
            "θ₁",
            div(()),
            textbox(data.theta1.to_degrees())
                .map_state(|data: &mut AppData| &mut data.theta1)
                .map_message_result(|data: &mut AppData, r| {
                    data.theta1 = data.theta1.to_radians();
                    data.maintain_arrangement(sheet::Handle::C1);
                    r
                }),
        );
        let frag_kappa0 = labeled_valued(
            "κ₀",
            div(()),
            textbox(data.kappa0)
                .map_state(|data: &mut AppData| &mut data.kappa0)
                .map_message_result(|data: &mut AppData, r| {
                    data.maintain_arrangement(sheet::Handle::C0);
                    r
                }),
        );
        let frag_kappa1 = labeled_valued(
            "κ₁",
            div(()),
            textbox(data.kappa1)
                .map_state(|data: &mut AppData| &mut data.kappa1)
                .map_message_result(|data: &mut AppData, r| {
                    data.maintain_arrangement(sheet::Handle::C1);
                    r
                }),
        );

        let frag_loopy = button(if data.loopy {
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
            data.maintain_arrangement(sheet::Handle::C0);
        });

        let frag_running = button(if data.running { "Pause" } else { "Resume" })
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

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> + use<> {
    let memo = state.memo.try_update(state.data, memoized_app_logic);

    app_view!(
        state,
        memo,
        (memo.frag_results.clone()),
        (memo.frag_controls.clone()),
        |state: &mut AppState, sheet::DragAction { handle, point }| {
            match handle {
                sheet::Handle::C0 => state.data.theta0 = point.to_vec2().angle(),
                sheet::Handle::C1 => {
                    state.data.theta1 =
                        (Point::new(BASE_WIDTH, 0.) - Affine::FLIP_Y * point).angle()
                }
            };
            state.data.maintain_arrangement(handle);
        }
    )
}
