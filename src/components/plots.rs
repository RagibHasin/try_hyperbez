use wasm_bindgen::JsCast;
use xilem_web::{
    elements::{
        html::div,
        svg::{g, svg, text},
    },
    interfaces::Element,
    svg::kurbo::{Affine, BezPath, Circle, Line, PathEl, Point, Size},
    DomView,
};

#[derive(Debug, PartialEq)]
pub struct State {
    size: Size,
    hovered_x: Option<f64>,
}

impl Default for State {
    fn default() -> Self {
        Self {
            size: Size::new(760., 540.),
            hovered_x: None,
        }
    }
}

impl State {
    pub fn view(&mut self, theta: &[f64], kappa: &[f64]) -> impl DomView<Self> {
        let mut plot_size = 1.5 * self.size;
        plot_size.height /= 2.;
        div((
            plot(theta, plot_size, self.hovered_x, "θ (°)")
                .map_state(|state: &mut Self| &mut state.hovered_x),
            plot(kappa, plot_size, self.hovered_x, "κ")
                .map_state(|state: &mut Self| &mut state.hovered_x),
        ))
        .on_resize(|state: &mut Self, e| {
            state.size.width = e.content_rect().width();
            state.size.height = e.content_rect().height();
        })
        .id("plots")
    }

    pub fn hovered_x(&self) -> Option<f64> {
        self.hovered_x
    }
}

pub fn plot(
    values: &[f64],
    size: Size,
    hover: Option<f64>,
    caption: &str,
) -> impl DomView<Option<f64>> {
    let Size { width, height } = size;
    let axis_width = width - 70.;
    let axis_height = height - 40.;

    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let trace = BezPath::from_iter(
        std::iter::once(PathEl::MoveTo(Point::new(0., values[0]))).chain(
            values
                .iter()
                .inspect(|v| {
                    min = min.min(**v);
                    max = max.max(**v);
                })
                .enumerate()
                .skip(1)
                .map(|(i, v)| PathEl::LineTo(Point::new(i as f64 * 1e-3, *v))),
        ),
    );

    let y_tick_interval = make_nice((max - min) * 5. / axis_height);
    let y_tick_resolution = n_dec_digits((y_tick_interval * 10.).fract(), 3);

    let min_i = (min / y_tick_interval).floor() as i32;
    let max_i = (max / y_tick_interval).ceil() as i32;
    let min = min_i as f64 * y_tick_interval;
    let max = max_i as f64 * y_tick_interval;

    let x_axis_y = 0f64.clamp(min, max);

    let axis_height = axis_height + if x_axis_y == 0. { 20. } else { 0. };
    let scale_y = axis_height / (max - min);
    let tx = Affine::FLIP_Y * Affine::scale_non_uniform(axis_width, scale_y);

    let x_axis = Vec::from_iter((0..=100).map(|i| {
        let x = i as f64 * 1e-2 * axis_width;
        let y = -x_axis_y * scale_y;
        let dy = if i % 10 == 0 { 5. } else { 3. };
        g((
            Line::new((x, y - dy), (x, y + dy)),
            (i % 10 == 0 && i != 0).then(|| {
                text(format!("{:.1}", i as f64 * 1e-2))
                    .attr("x", x)
                    .attr("y", y + 10.)
            }),
        ))
        .class("tick")
    }));
    let x_axis = g(((tx * Line::new((0., x_axis_y), (1., x_axis_y))), x_axis)).class(["axis", "x"]);

    let y_axis = Vec::from_iter((min_i..=max_i).map(|i| {
        let u = i as f64 * y_tick_interval;
        let y = u * -scale_y;
        let dx = if i % 10 == 0 { 5. } else { 3. };
        g((
            Line::new((-dx, y), (dx, y)),
            (i % 10 == 0).then(|| {
                text(format!("{:.*}", y_tick_resolution, u))
                    .attr("x", -10.)
                    .attr("y", y)
            }),
        ))
        .class("tick")
    }));
    let y_axis = g((
        (tx * Line::new((0., min), (0., max))),
        y_axis,
        text(caption.to_string())
            .attr("x", 10.)
            .attr("y", max_i as f64 * y_tick_interval * -scale_y)
            .class("caption"),
    ))
    .class(["axis", "y"]);

    let trace = (tx * trace).class("trace");

    let marker = hover.map(|s| {
        let v = values[(s * 1e3) as usize];
        g((
            (tx * Line::new((s, min), (s, max))), // x marker
            (tx * Line::new((0., v), (1., v))),   // y marker
            Circle::new(tx * Point::new(s, v), 3.),
        ))
        .class("hover")
    });

    svg((x_axis, y_axis, trace, marker))
        .on_mousemove(move |state: &mut Option<f64>, e| {
            let scale = width
                / e.current_target()
                    .unwrap()
                    .unchecked_into::<web_sys::Element>()
                    .client_width() as f64;
            let s = (e.offset_x() as f64 * scale - 50.) / axis_width;
            *state = (0.0..=1.).contains(&s).then_some(s);
        })
        .on_mouseleave(|state: &mut Option<f64>, _| *state = None)
        .attr(
            "viewBox",
            format!("-50 {} {width} {height}", -max * scale_y - 10.),
        )
        .class("plot")
}

fn n_dec_digits(mut fract: f64, limit: usize) -> usize {
    for i in 0..limit {
        if fract == 0. {
            return i;
        } else {
            fract = (fract * 10.).fract()
        }
    }
    limit
}

fn make_nice(v: f64) -> f64 {
    [10f64.powf(v.log10().floor()), 5f64.powf(v.log(5.).floor())]
        .into_iter()
        .chain((v > 0.25).then(|| 2f64.powf(v.log2().floor())))
        .map(|l| (l, (v - l).abs()))
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .unwrap()
        .0
}
