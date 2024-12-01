use web_sys::wasm_bindgen::JsCast;
use xilem_web::{
    elements::{
        html::div,
        svg::{g, svg},
    },
    interfaces::*,
    svg::kurbo::{self, Affine, Circle, Line, ParamCurve, Point, Shape, Size, Vec2},
    DomView,
};

use hyperbez_toy::*;

use crate::*;

use components::*;

#[derive(Debug, PartialEq)]
pub(crate) struct AppState {
    p1: Point,
    p2: Point,

    sheet_origin: Point,
    sheet_zoom: f64,
    hovered_s: Option<f64>,
    drag: DragElement,
}

impl AppState {
    fn sheet_scale(&self, e: &web_sys::MouseEvent) -> f64 {
        let sheet = e
            .current_target()
            .unwrap()
            .unchecked_into::<web_sys::Element>();
        self.sheet_zoom * 1200. / sheet.client_width() as f64
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            p1: Point::new(50., 200.),
            p2: Point::new(150., 200.),
            sheet_origin: Point::new(-350., -450.),
            sheet_zoom: 1.,
            hovered_s: None,
            drag: DragElement::None,
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum DragElement {
    None,
    Sheet,
    P0,
    P1,
}

pub(crate) fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let base_width = 500.;

    let p0 = Point::ZERO;
    let p3 = Point::new(base_width, 0.);
    let scale_down = kurbo::TranslateScale::scale(1. / base_width);
    let cubicbez = kurbo::CubicBez::new(p0, state.p1, state.p2, p3);

    let params =
        hb_extra::HyperbezParams::from_control(scale_down * state.p1, scale_down * state.p2);
    let hyperbez = hb_extra::Hyperbezier::from_points_params(params, p0, p3);

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
        Circle::new(control0, NODE_RADIUS).on_mousedown(|state: &mut AppState, e| {
            state.drag = DragElement::P0;
            e.stop_propagation();
        }),
    );

    let control1 = Affine::FLIP_Y * state.p2;
    let control1 = (
        Line::new((base_width, 0.), control1),
        Circle::new(control1, NODE_RADIUS).on_mousedown(|state: &mut AppState, e| {
            state.drag = DragElement::P1;
            e.stop_propagation();
        }),
    );

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
        div((frag_hovered_s, frag_hovered_theta, frag_hovered_kappa)).class("results"),
    );

    let plot_size = Size::new(760., 270.);
    let frag_plots = div((
        plot(&theta, plot_size, state.hovered_s, "θ (°)")
            .map_state(|state: &mut AppState| &mut state.hovered_s),
        plot(&kappa, plot_size, state.hovered_s, "κ")
            .map_state(|state: &mut AppState| &mut state.hovered_s),
    ))
    .id("plots");

    let sheet_size = state.sheet_zoom * Size::new(1200., 900.);
    let frag_svg = svg((
        g((control0, control1)).class("control"),
        cubicbez.id("cubicbez"),
        path.id("hyperbez"),
        hover_mark,
        g(points).id("nodes"),
    ))
    .attr(
        "viewBox",
        format!(
            "{} {} {} {}",
            state.sheet_origin.x, state.sheet_origin.y, sheet_size.width, sheet_size.height,
        ),
    )
    .on_mousedown(|state: &mut AppState, _| state.drag = DragElement::Sheet)
    .on_mouseup(|state: &mut AppState, _| state.drag = DragElement::None)
    .on_mousemove(move |state: &mut AppState, e| {
        if let DragElement::None = state.drag {
            return;
        };

        let scale = state.sheet_scale(&e);
        let p = Affine::FLIP_Y
            * Affine::scale(scale).then_translate(state.sheet_origin.to_vec2())
            * Point::new(e.offset_x() as f64, e.offset_y() as f64);
        tracing::trace!(?p);

        match state.drag {
            DragElement::P0 => state.p1 = p,
            DragElement::P1 => state.p2 = p,
            DragElement::Sheet => {
                let delta = scale * Vec2::new(e.movement_x() as f64, e.movement_y() as f64);
                tracing::trace!(?delta);
                state.sheet_origin -= delta
            }
            DragElement::None => unreachable!(),
        }
    })
    .on_wheel(|state: &mut AppState, e| {
        e.prevent_default();

        let factor = if e.delta_y() > 0. {
            1.25
        } else if e.delta_y() < 0. {
            0.8
        } else {
            1.
        };

        let origin_delta = (factor - 1.)
            * state.sheet_scale(&e)
            * Vec2::new(e.offset_x() as f64, e.offset_y() as f64);
        tracing::trace!(factor, ?origin_delta);
        state.sheet_origin -= origin_delta;
        state.sheet_zoom *= factor;
    })
    .passive(false);

    div((
        div((div(frag_results).id("ui"), frag_plots)).id("pane-left"),
        div(frag_svg).id("render-sheet"),
    ))
    .id("app-root")
}
