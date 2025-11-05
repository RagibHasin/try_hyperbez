// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use web_sys::MouseEvent;
use xilem_web::{
    core::Edit,
    elements::svg::svg,
    interfaces::Element,
    svg::kurbo::{Point, Size, Vec2},
    Action, DomFragment, DomView,
};

#[derive(Debug, PartialEq)]
pub struct State<DragData = NoData> {
    size: Size,
    origin: Point,
    zoom: f64,
    drag: DragElement<DragData>,
}

impl<O> Default for State<O> {
    fn default() -> Self {
        Self {
            size: Size::new(1200., 900.),
            origin: Point::new(-350., -450.),
            zoom: 1.5,
            drag: DragElement::None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DragElement<O = NoData> {
    None,
    Sheet,
    Other(O),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DragAction<O = NoData> {
    pub data: O,
    pub event: MouseEvent,
}
impl<O> Action for DragAction<O> {}

impl<DragData: Copy + 'static> State<DragData> {
    pub fn view<Children: DomFragment<Edit<Self>, DragAction<DragData>>>(
        &mut self,
        children: Children,
    ) -> impl DomView<Edit<Self>, DragAction<DragData>> {
        let sheet_size = self.zoom * self.size;
        svg(children)
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
                if let DragElement::None = state.drag {
                    return None;
                };

                match state.drag {
                    DragElement::Other(data) => Some(DragAction { data, event: e }),
                    DragElement::Sheet => {
                        let delta =
                            state.zoom * Vec2::new(e.movement_x() as f64, e.movement_y() as f64);
                        state.origin -= delta;
                        None
                    }
                    DragElement::None => unreachable!(),
                }
            })
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
                let new_width = e.content_rect().width();
                let new_height = e.content_rect().height();
                state.origin.x -= (new_width * state.zoom - state.size.width) / 2.;
                state.origin.y -= (new_height * state.zoom - state.size.height) / 2.;
                state.size.width = new_width;
                state.size.height = new_height;
            })
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

    pub fn set_drag(&mut self, e: Option<DragData>) {
        self.drag = if let Some(e) = e {
            DragElement::Other(e)
        } else {
            DragElement::None
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum NoData {}
