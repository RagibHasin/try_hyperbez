// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use wasm_bindgen::JsCast;
use xilem_web::{
    elements::html::{div, input, span},
    interfaces::*,
    DomFragment, DomView,
};

#[derive(Debug)]
pub struct Memoized<K, V> {
    key: K,
    value: Option<V>,
}

impl<K: PartialEq, V> Memoized<K, V> {
    pub fn update(&mut self, data: K, logic: impl FnOnce(&K) -> V) -> &V {
        self.try_update(data, |data, _| Some(logic(data)))
    }

    pub fn try_update(
        &mut self,
        data: K,
        logic: impl FnOnce(&K, Option<&mut V>) -> Option<V>,
    ) -> &V {
        if (self.key != data) || self.value.is_none() {
            if let Some(value) = logic(&data, self.value.as_mut()) {
                self.value = Some(value);
            }
            self.key = data;
        }
        self.value()
    }

    pub fn value(&self) -> &V {
        self.value.as_ref().unwrap()
    }
}

impl<K: Default, V> Default for Memoized<K, V> {
    fn default() -> Self {
        Self {
            key: Default::default(),
            value: None,
        }
    }
}

pub fn slider(value: f64, min: f64, max: f64, step: f64) -> impl DomView<f64> {
    input(())
        .attr("value", value)
        .attr("type", "range")
        .attr("min", min)
        .attr("max", max)
        .attr("step", step)
        .on_input(|state: &mut f64, e| {
            if let Ok(val_f64) = e
                .target()
                .unwrap()
                .unchecked_into::<web_sys::HtmlInputElement>()
                .value()
                .parse::<f64>()
            {
                *state = val_f64;
            }
        })
}

pub fn textbox(value: f64) -> impl DomView<f64> {
    input(())
        .attr("value", format!("{:.1}", value))
        .attr("type", "text")
        .on_change(|state: &mut f64, e| {
            if let Ok(val_f64) = e
                .target()
                .unwrap()
                .unchecked_into::<web_sys::HtmlInputElement>()
                .value()
                .parse::<f64>()
            {
                *state = val_f64;
            }
        })
}

pub fn labeled_valued<T: 'static>(
    label: impl DomFragment<T>,
    edit: impl DomFragment<T>,
    value: impl DomFragment<T>,
) -> impl DomFragment<T> {
    (span(label).class("label"), edit, span(value).class("label"))
}

pub fn spacer<T: 'static>() -> impl DomView<T> {
    div(()).class("spacer")
}

pub mod plots;
pub mod sheet;
