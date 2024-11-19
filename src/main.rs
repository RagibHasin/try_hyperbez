// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use xilem_web::{core::one_of::Either, App, DomView};

pub mod components;

mod jsport;
mod old;

enum AppState {
    Old(old::AppState),
    JsPort(jsport::AppState),
}

fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    match state {
        AppState::Old(state) => {
            Either::A(old::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::Old(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
        AppState::JsPort(state) => {
            Either::B(jsport::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::JsPort(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
    }
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

    let state = AppState::JsPort(Default::default());
    App::new(xilem_web::document_body(), state, app_logic).run();
}
