use xilem_web::svg::kurbo::Point;

pub(crate) const BASE_WIDTH: f64 = 500.;

pub(crate) const P0: Point = Point::ZERO;
pub(crate) const P3: Point = Point::new(BASE_WIDTH, 0.);

macro_rules! app_view {
    ($state:ident, $memo:ident, ($($frag_results:expr),+), ($($frag_sheet:expr),*), $map_action:expr) => {{
        let mut hovered_point = None;
        let mut hovered_theta = None;
        let mut hovered_kappa = None;
        let mut hover_mark = None;
        if let Some(s) = $state.plots.hovered_x {
            let p = kurbo::ParamCurve::eval(&$memo.hyperbez, s);

            hovered_point = Some(p);
            let i = (s * 1e3) as usize;
            (hovered_theta, hovered_kappa) = (Some($memo.theta[i]), Some($memo.kappa[i]));
            let (theta, kappa) = ($memo.theta[i].to_radians(), $memo.kappa[i]);

            let r_curv = $memo.hyperbez.scale_rot().length() / kappa;
            let tangent_half = 0.25 * BASE_WIDTH * kurbo::Vec2::from_angle(theta);
            let tangent = Affine::FLIP_Y * Line::new(p - tangent_half, p + tangent_half);
            let r_curv = Affine::FLIP_Y
                * Line::new(
                    p,
                    p + r_curv * kurbo::Vec2::from_angle(theta + std::f64::consts::FRAC_PI_2),
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
            $state
                .plots
                .hovered_x
                .map_or(empty.clone(), |s| format!("{:.3}", s)),
        );
        let frag_hovered_p_x = labeled_valued(
            ("P", html::sub("x"), "(s): "),
            (),
            hovered_point.map_or(empty.clone(), |v| format!("{:.2}", v.x)),
        );
        let frag_hovered_p_y = labeled_valued(
            ("P", html::sub("y"), "(s): "),
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
            $($frag_results,)*
            div((
                frag_hovered_s,
                frag_hovered_p_x,
                frag_hovered_p_y,
                frag_hovered_theta,
                frag_hovered_kappa,
            ))
            .class("results"),
        );

        let frag_plots = $state
            .plots
            .view(&*$memo.theta, &*$memo.kappa)
            .map_state(|state: &mut AppState| &mut state.plots);

        let frag_svg = $state
            .sheet
            .view((
                $($frag_sheet,)*
                $memo.frag_path.clone(),
                hover_mark,
                $memo.frag_points.clone(),
            ))
            .map_state(|state: &mut AppState| &mut state.sheet)
            .map_action($map_action)
            .map_message_result(|state: &mut AppState, r| {
                if let Some(p) = state.sheet.hovered_pt() {
                    let hit = xilem_web::svg::kurbo::ParamCurveNearest::nearest(&state.memo.value().hyperbez, p, 0.005);
                    tracing::trace!(%p, ?hit, zoom = state.sheet.zoom());
                    state.plots.hovered_x =
                        (hit.distance_sq * state.sheet.zoom().powi(-2) < 25.).then_some(hit.t);
                }
                r
            });

        div((
            div((
                div((
                    $memo.frag_options
                        .clone()
                        .map_state(|state: &mut AppState| &mut state.data),
                    frag_results,
                ))
                .id("ui"),
                frag_plots,
            ))
            .id("pane-left"),
            frag_svg,
        ))
        .id("app-root")
    }};
}

pub(crate) use app_view;
