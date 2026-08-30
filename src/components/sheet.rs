// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use xilem_web::{
    Action, DomFragment, DomView,
    elements::{
        html::{button, div},
        svg::svg,
    },
    interfaces::Element,
    svg::kurbo::{Affine, Point, Size, Vec2},
};

#[derive(Debug, PartialEq)]
pub struct State {
    size: Size,
    origin: Point,
    zoom: f64,
    drag: DragElement,
    hovered_pt: Option<Point>,
}

impl Default for State {
    fn default() -> Self {
        Self {
            size: Size::new(1200., 900.),
            origin: Point::new(-350., -450.),
            zoom: 1.5,
            drag: DragElement::None,
            hovered_pt: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DragElement {
    None,
    Sheet,
    Handle(Handle),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Handle {
    C0,
    C1,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DragAction {
    pub handle: Handle,
    pub point: Point,
}
impl Action for DragAction {}

impl State {
    pub fn view<Children: DomFragment<Self, DragAction>>(
        &mut self,
        children: Children,
    ) -> impl DomView<Self, DragAction> + use<Children> {
        let sheet_size = self.zoom * self.size;
        let sheet = svg(children)
            .attr(
                "viewBox",
                format!(
                    "{} {} {} {}",
                    self.origin.x, self.origin.y, sheet_size.width, sheet_size.height,
                ),
            )
            .on_mousedown(|state: &mut Self, _| state.drag = DragElement::Sheet)
            .on_mouseup(|state: &mut Self, _| state.drag = DragElement::None)
            .on_mousemove(move |state: &mut Self, e| {
                let p = Affine::FLIP_Y
                    * Affine::scale(state.zoom).then_translate(state.origin.to_vec2())
                    * Point::new(e.offset_x() as f64, e.offset_y() as f64);
                state.hovered_pt = Some(p);
                match state.drag {
                    DragElement::Handle(data) => Some(DragAction {
                        handle: data,
                        point: p,
                    }),
                    DragElement::Sheet => {
                        let delta =
                            state.zoom * Vec2::new(e.movement_x() as f64, e.movement_y() as f64);
                        state.origin -= delta;
                        None
                    }
                    DragElement::None => None,
                }
            })
            .on_mouseleave(|state: &mut Self, _| state.hovered_pt = None)
            .on_wheel(|state: &mut Self, e| {
                e.prevent_default();

                let factor = if e.delta_y() > 0. {
                    1.25
                } else if e.delta_y() < 0. {
                    0.8
                } else {
                    1.
                };

                let origin_delta = (factor - 1.)
                    * state.zoom
                    * Vec2::new(e.offset_x() as f64, e.offset_y() as f64);
                state.origin -= origin_delta;
                state.zoom *= factor;
            })
            .passive(false)
            .on_resize(|state: &mut Self, e| {
                let new_size = Size::new(e.content_rect().width(), e.content_rect().height());
                state.origin -= (new_size * state.zoom - state.size).to_vec2() / 2.;
                state.size = new_size;
            });

        div((
            sheet,
            div(button("Reset View").on_click(|state: &mut Self, _| {
                state.origin = Point::new(250., 0.) - state.size.to_vec2() * 0.75;
                state.zoom = 1.5;
            }))
            .id("toolbar"),
        ))
        .id("render-sheet")
    }

    pub fn size(&self) -> Size {
        self.size
    }

    pub fn origin(&self) -> Point {
        self.origin
    }

    pub fn zoom(&self) -> f64 {
        self.zoom
    }

    pub fn hovered_pt(&self) -> Option<Point> {
        self.hovered_pt
    }

    pub fn set_drag_handle(&mut self, handle: Option<Handle>) {
        self.drag = if let Some(handle) = handle {
            DragElement::Handle(handle)
        } else {
            DragElement::None
        }
    }
}
