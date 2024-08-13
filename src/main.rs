// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use web_sys::wasm_bindgen::JsCast;
use xilem_web::{
    elements::{
        html::{button, div},
        svg::{g, svg},
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, ParamCurve, Point, Rect, Size, Vec2},
    App, DomView,
};

use hyperbez_toy::*;

mod components;
use components::*;

struct AppState {
    a: f64,
    b: f64,
    c: f64,
    d: f64,

    optimize: bool,
    accuracy_order: f64,

    hovered_s: Option<f64>,
    drag: ControlDrag,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            a: -1.,
            b: -1.,
            c: 1.,
            d: 0.,
            optimize: false,
            accuracy_order: 1.,
            hovered_s: None,
            drag: ControlDrag::None,
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
enum ControlDrag {
    #[default]
    None,
    P0,
    P1,
}

fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let d_limit = hb::d_limit_rounded(state.a, state.b, state.c);
    state.d = state.d.clamp(d_limit.start, d_limit.end);

    let base_width = 500.;
    let hyperbez = hb::Hyperbezier::from_points_params(
        hb::HyperbezParams::new(state.a, state.b, state.c, state.d),
        Point::ZERO,
        Point::new(base_width, 0.),
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
    let arclen = hyperbez.scale_rot().length() / base_width;
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

    let qoppa0 = hb::qoppa(hyperbez.params(), 0.);
    let qoppa1 = hb::qoppa(hyperbez.params(), 1.);
    let ka_qo0 = hb::frac_kappa_qoppa(hyperbez.params(), 0.);
    let ka_qo1 = hb::frac_kappa_qoppa(hyperbez.params(), 1.);

    const NODE_RADIUS: f64 = 5.;
    let points = path
        .elements()
        .iter()
        .filter_map(|e| e.end_point())
        .map(|p| Circle::new(p, NODE_RADIUS))
        .collect::<Vec<_>>();

    let tension0 = hb::tension(&hyperbez, 0.);
    let control0 = Affine::FLIP_Y * (tension0 * Vec2::from_angle(theta0.to_radians())).to_point();
    let control0 = (
        Line::new((0., 0.), control0),
        Circle::new(control0, NODE_RADIUS)
            .on_mousedown(|state: &mut AppState, _| state.drag = ControlDrag::P0),
    );

    let tension1 = hb::tension(&hyperbez, 1.);
    let control1 = Affine::FLIP_Y
        * (Point::new(base_width, 0.) - tension1 * Vec2::from_angle(theta1.to_radians()));
    let control1 = (
        Line::new((base_width, 0.), control1),
        Circle::new(control1, NODE_RADIUS)
            .on_mousedown(|state: &mut AppState, _| state.drag = ControlDrag::P1),
    );

    let trigon0 = {
        let apex = Affine::FLIP_Y
            * (hb::hypo_of_qoppa(&hyperbez, 0.) * Vec2::from_angle(theta0.to_radians())).to_point();
        let mut trigon = kurbo::BezPath::new();
        trigon.move_to((0., 0.));
        trigon.line_to((apex.x, 0.));
        trigon.line_to(apex);
        trigon.class("area-qoppa")
    };

    let trigon1 = {
        let base1 = Point::new(base_width, 0.);
        let apex = Affine::FLIP_Y
            * (base1 - hb::hypo_of_qoppa(&hyperbez, 1.) * Vec2::from_angle(theta1.to_radians()));
        let mut trigon = kurbo::BezPath::new();
        trigon.move_to(base1);
        trigon.line_to((apex.x, 0.));
        trigon.line_to(apex);
        trigon.class("area-qoppa")
    };

    let horton0 = {
        let apex = Affine::FLIP_Y
            * (hb::dust_of_ka_qo(&hyperbez, 0.) * Vec2::from_angle(theta0.to_radians())).to_point();
        let mut horton = kurbo::BezPath::new();
        horton.move_to((0., 0.));
        horton.line_to((apex.x, 0.));
        horton.line_to(apex);
        horton.class("area-kappa-qoppa")
    };

    let horton1 = {
        let base1 = Point::new(base_width, 0.);
        let apex = Affine::FLIP_Y
            * (base1 - hb::dust_of_ka_qo(&hyperbez, 1.) * Vec2::from_angle(theta1.to_radians()));
        let mut horton = kurbo::BezPath::new();
        horton.move_to(base1);
        horton.line_to((apex.x, 0.));
        horton.line_to(apex);
        horton.class("area-kappa-qoppa")
    };

    let mut hovered_theta = None;
    let mut hovered_kappa = None;
    let mut hover_mark = None;
    if let Some(s) = state.hovered_s {
        let i = (s * 1e3) as usize;
        (hovered_theta, hovered_kappa) = (Some(theta[i]), Some(kappa[i]));
        let (theta, kappa) = (theta[i].to_radians(), kappa[i]);

        let p = hyperbez.eval(s);
        let r_curv = hyperbez.scale_rot().length() / kappa;
        let tangent_half = 0.25 * base_width * Vec2::from_angle(theta);
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

    let frag_qoppa0 = labeled_valued("ϟ₀: ", (), format!("{:.3}", qoppa0));
    let frag_qoppa1 = labeled_valued("ϟ₁: ", (), format!("{:.3}", qoppa1));
    let frag_ka_qo0 = labeled_valued("κ₀/ϟ₀: ", (), format!("{:.3}", ka_qo0));
    let frag_ka_qo1 = labeled_valued("κ₁/ϟ₁: ", (), format!("{:.3}", ka_qo1));

    let empty = "~".to_string();
    let frag_hovered_s = labeled_valued(
        "s: ",
        (),
        state
            .hovered_s
            .map_or(empty.clone(), |s| format!("{:.3}", s)),
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
            frag_qoppa0,
            frag_qoppa1,
            spacer(),
            frag_ka_qo0,
            frag_ka_qo1,
            spacer(),
            frag_hovered_s,
            frag_hovered_theta,
            frag_hovered_kappa,
        ))
        .class("results"),
    );

    let plot_size = Size::new(760., 270.);
    let frag_plots = div((
        plot(&theta, plot_size, state.hovered_s, "θ (°)")
            .map_state(|state: &mut AppState| &mut state.hovered_s),
        plot(&kappa, plot_size, state.hovered_s, "κ")
            .map_state(|state: &mut AppState| &mut state.hovered_s),
    ))
    .id("plots");

    let sheet_rect = Rect::from_origin_size((-350., -450.), (1200., 900.));
    let frag_svg = svg((
        g((trigon0, trigon1, horton0, horton1, control0, control1)).class("control"),
        path.id("hyperbez"),
        hover_mark,
        g(points).id("nodes"),
    ))
    .attr(
        "viewBox",
        format!(
            "{} {} {} {}",
            sheet_rect.x0,
            sheet_rect.y0,
            sheet_rect.width(),
            sheet_rect.height(),
        ),
    )
    .on_mouseup(|state: &mut AppState, _| state.drag = ControlDrag::None)
    .on_mousemove(move |state: &mut AppState, e| {
        let sheet = e
            .current_target()
            .unwrap()
            .unchecked_into::<web_sys::Element>();
        let scale = sheet_rect.width() / sheet.client_width() as f64;
        let p = Affine::FLIP_Y
            * Affine::scale(scale)
                .then_translate(sheet_rect.origin().to_vec2())
                .then_scale(1. / base_width)
            * Point::new(e.offset_x() as f64, e.offset_y() as f64);
        if let ControlDrag::P0 = state.drag {
            let control0 = p.to_vec2();
            let theta0 = control0.angle();
            let kappa0 = hb::tau_to_kappa(control0.hypot());
            let theta1 = theta1.to_radians();
            tracing::trace!(?control0, theta0, kappa0, theta1, kappa1);
            let soln = hb::solve_for_params_inferred(theta0, theta1, kappa0, kappa1, 1e-4, 10);
            tracing::trace!(?soln);
            if let [Ok(params), ..] = soln {
                tracing::trace!(?params);

                state.a = params.a();
                state.b = params.b();
                state.c = params.c();
                state.d = params.d();
            }
        } else if let ControlDrag::P1 = state.drag {
            let control1 = Vec2::new(1., 0.) - p.to_vec2();
            let theta1 = control1.angle();
            let kappa1 = hb::tau_to_kappa(control1.hypot());
            let theta0 = theta0.to_radians();
            tracing::trace!(?control1, theta0, kappa0, theta1, kappa1);
            let soln = hb::solve_for_params_inferred(theta0, theta1, kappa0, kappa1, 1e-4, 10);
            tracing::trace!(?soln);
            if let [Ok(params), ..] = soln {
                tracing::trace!(?params);

                state.a = params.a();
                state.b = params.b();
                state.c = params.c();
                state.d = params.d();
            }
        }
    });

    let frag_a = labeled_valued(
        "a: ",
        slider(state.a, -20., 20., 0.2).map_state(move |state: &mut AppState| &mut state.a),
        format!("{:.1}", state.a),
    );
    let frag_b = labeled_valued(
        "b: ",
        slider(state.b, -20., 20., 0.2).map_state(move |state: &mut AppState| &mut state.b),
        format!("{:.1}", state.b),
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
        format!("{:.1}", state.c),
    );
    let frag_d = labeled_valued(
        "d: ",
        slider(state.d, d_limit.start, d_limit.end, 0.1)
            .map_state(move |state: &mut AppState| &mut state.d),
        format!("{:.1}", state.d),
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

    div((
        div((
            div((
                div((
                    frag_a,
                    frag_b,
                    frag_c,
                    frag_d,
                    spacer(),
                    frag_optimize,
                    frag_accuracy,
                ))
                .id("options"),
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

pub fn main() {
    use tracing_subscriber::{fmt::format::Pretty, prelude::*};
    use tracing_web::{performance_layer, MakeWebConsoleWriter};

    console_error_panic_hook::set_once();

    let fmt_layer = tracing_subscriber::fmt::layer()
        .pretty()
        .with_ansi(false)
        .without_time()
        .with_writer(MakeWebConsoleWriter::new());
    let perf_layer = performance_layer().with_details_from_fields(Pretty::default());

    tracing_subscriber::registry()
        .with(fmt_layer)
        .with(perf_layer)
        .init();

    let state = AppState::default();
    App::new(xilem_web::document_body(), state, app_logic).run();
}
