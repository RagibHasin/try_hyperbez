use xilem_web::{
    elements::{
        html::{self, button, div},
        svg::g,
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, ParamCurve, Point, Shape, Vec2},
    DomView,
};

use hyperbez_toy::*;

use crate::components::*;

#[derive(Debug, PartialEq)]
pub(crate) struct AppState {
    p1: Point,
    p2: Point,

    symmetric: bool,

    plots: plots::State,
    sheet: sheet::State<Handle>,
}

const BASE_WIDTH: f64 = 500.;

impl AppState {
    fn maintain_symmetry(&mut self, handle: Handle) {
        if self.symmetric {
            let (dest, src) = match handle {
                Handle::P1 => (&mut self.p2, self.p1),
                Handle::P2 => (&mut self.p1, self.p2),
            };
            *dest = Point::new(BASE_WIDTH - src.x, src.y);
        }
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            p1: Point::new(50., 200.),
            p2: Point::new(450., 200.),
            symmetric: true,
            plots: Default::default(),
            sheet: Default::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Handle {
    P1,
    P2,
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let p0 = Point::ZERO;
    let p3 = Point::new(BASE_WIDTH, 0.);
    let scale_down = kurbo::TranslateScale::scale(1. / BASE_WIDTH);
    let cubicbez = kurbo::CubicBez::new(p0, state.p1, state.p2, p3);

    // let (params, raw_params, n_iter, err, err_type) = match utils::solve_helper(
    //     scale_down * state.p1,
    //     scale_down * p2,
    //     utils::solve_inferring(hb_extra::solver::solve_for_params_exact),
    // ) {
    //     Ok(s) => (
    //         hb_extra::HyperbezParams::new(s.params.x, s.params.y, s.params.z, s.params.w, 1.),
    //         s.params.data.0[0],
    //         s.iter,
    //         [s.err.x, s.err.y, s.err.z, s.err.w],
    //         "None",
    //     ),
    //     Err(hb_extra::solver::SolveError::Singularity { guess, err, iter }) => (
    //         hb_extra::HyperbezParams::new(0., -1., -1., 1., 1.),
    //         guess.data.0[0],
    //         iter,
    //         [err.x, err.y, err.z, err.w],
    //         "Singularity",
    //     ),
    //     Err(hb_extra::solver::SolveError::OutOfIteration { guess, err }) => (
    //         hb_extra::HyperbezParams::new(0., -1., -1., 1., 1.),
    //         guess.data.0[0],
    //         11,
    //         [err.x, err.y, err.z, err.w],
    //         "OOI",
    //     ),
    // };

    // tracing::trace!(?raw_params, n_iter, ?err, err_type);

    // let params = hb_extra::HyperbezParams::from_control(scale_down * state.p1, scale_down * p2);

    let raw_params = utils::solve_helper_ext(
        scale_down * state.p1,
        scale_down * state.p2,
        1e-2,
        11,
        hb::solver::dual::solve_inferring_full,
    );
    let params = hb::HyperbezParams::new(
        raw_params[0],
        raw_params[1],
        raw_params[2],
        raw_params[3],
        1.,
    );
    let hyperbez = hb::Hyperbezier::from_points_params(params, p0, p3);

    // {
    //     let arg_uv = hyperbez.params().integrate(1.).angle().to_degrees();
    //     let int_theta = hb_extra::integrate_theta(hyperbez.params(), 1.).to_degrees();
    //     tracing::trace!(?arg_uv, ?int_theta);
    // }

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

    let control0 = Affine::FLIP_Y * state.p1;
    let control0 = (
        Line::new((0., 0.), control0),
        Circle::new(control0, NODE_RADIUS).on_mousedown(|state: &mut sheet::State<Handle>, e| {
            state.set_drag(Some(Handle::P1));
            e.stop_propagation();
        }),
    );

    let control1 = Affine::FLIP_Y * state.p2;
    let control1 = (
        Line::new((BASE_WIDTH, 0.), control1),
        Circle::new(control1, NODE_RADIUS).on_mousedown(|state: &mut sheet::State<Handle>, e| {
            state.set_drag(Some(Handle::P2));
            e.stop_propagation();
        }),
    );

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

    let frag_a = labeled_valued("a: ", (), format!("{:.3}", params.a()));
    let frag_b = labeled_valued("b: ", (), format!("{:.3}", params.b()));
    let frag_c = labeled_valued("c: ", (), format!("{:.3}", params.c()));
    let frag_d = labeled_valued("d: ", (), format!("{:.3}", params.d()));
    let frag_e = labeled_valued("e: ", (), format!("{:.3}", params.e()));

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
        div((frag_a, frag_b, frag_c, frag_d, frag_e)).class("results"),
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
        .view((
            g((control0, control1)).class("control"),
            cubicbez.id("cubicbez"),
            path.id("hyperbez"),
            hover_mark,
            g(points).id("nodes"),
        ))
        .adapt(move |state: &mut AppState, thunk| {
            thunk
                .call(&mut state.sheet)
                .map(|sheet::DragAction { data, event }| {
                    let p = Affine::FLIP_Y
                        * Affine::scale(state.sheet.zoom())
                            .then_translate(state.sheet.origin().to_vec2())
                        * Point::new(event.offset_x() as f64, event.offset_y() as f64);

                    *match data {
                        Handle::P1 => &mut state.p1,
                        Handle::P2 => &mut state.p2,
                    } = p;
                    state.maintain_symmetry(data);
                })
        });

    let frag_p1_x = labeled_valued(
        ("P1", html::sub("x"), ": "),
        div(()),
        textbox(state.p1.x).adapt(move |state: &mut AppState, thunk| {
            let res = thunk.call(&mut state.p1.x);
            state.maintain_symmetry(Handle::P1);
            res
        }),
    );
    let frag_p1_y = labeled_valued(
        ("P1", html::sub("y"), ": "),
        div(()),
        textbox(state.p1.y).adapt(move |state: &mut AppState, thunk| {
            let res = thunk.call(&mut state.p1.y);
            state.maintain_symmetry(Handle::P1);
            res
        }),
    );
    let frag_p2_x = labeled_valued(
        ("P2", html::sub("x"), ": "),
        div(()),
        textbox(state.p2.x).adapt(move |state: &mut AppState, thunk| {
            let res = thunk.call(&mut state.p2.x);
            state.maintain_symmetry(Handle::P2);
            res
        }),
    );
    let frag_p2_y = labeled_valued(
        ("P2", html::sub("y"), ": "),
        div(()),
        textbox(state.p2.y).adapt(move |state: &mut AppState, thunk| {
            let res = thunk.call(&mut state.p2.y);
            state.maintain_symmetry(Handle::P2);
            res
        }),
    );

    let frag_symmetric = button(if state.symmetric {
        "Make asymmetric"
    } else {
        "Make symmetric"
    })
    .on_click(|state: &mut AppState, _| state.symmetric = !state.symmetric);

    let frag_options = div((
        frag_p1_x,
        frag_p1_y,
        frag_p2_x,
        frag_p2_y,
        spacer(),
        frag_symmetric,
    ))
    .id("options");

    div((
        div((div((frag_options, frag_results)).id("ui"), frag_plots)).id("pane-left"),
        div(frag_svg).id("render-sheet"),
    ))
    .id("app-root")
}
