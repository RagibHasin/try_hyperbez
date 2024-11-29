// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use wasm_bindgen::JsCast;
use xilem_web::{
    core::one_of::OneOf3,
    elements::html::{div, option, select},
    interfaces::{Element, HtmlOptionElement},
    App, DomFragment, DomView,
};

pub mod components;

mod jsport;
mod old;
mod ptan;

enum AppState {
    Old(old::AppState),
    JsPort(jsport::AppState),
    PointTangent(ptan::AppState),
}

fn app_logic(state: &mut AppState) -> impl DomFragment<AppState> {
    let app = match state {
        AppState::Old(state) => {
            OneOf3::A(old::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::Old(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
        AppState::JsPort(state) => {
            OneOf3::B(jsport::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::JsPort(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
        AppState::PointTangent(state) => {
            OneOf3::C(ptan::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::PointTangent(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
    };

    let toolbar = div(select((
        option("Explorer")
            .value("old")
            .selected(matches!(state, AppState::Old(_))),
        option("JavaScript port")
            .value("jsport")
            .selected(matches!(state, AppState::JsPort(_))),
        option("Point · Tangent")
            .value("ptan")
            .selected(matches!(state, AppState::PointTangent(_))),
    ))
    .on_change(move |state: &mut AppState, e| {
        match e
            .target()
            .unwrap()
            .unchecked_into::<web_sys::HtmlSelectElement>()
            .value()
            .as_ref()
        {
            "old" => *state = AppState::Old(old::AppState::default()),
            "jsport" => *state = AppState::JsPort(jsport::AppState::default()),
            "ptan" => *state = AppState::PointTangent(ptan::AppState::default()),
            _ => {}
        }
    }))
    .id("toolbar");

    (app, toolbar)
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

    App::new(
        xilem_web::document_body(),
        AppState::PointTangent(Default::default()),
        app_logic,
    )
    .run();
}
