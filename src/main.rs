// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use wasm_bindgen::JsCast;
use xilem_web::{
    core::one_of::OneOf5,
    elements::html::{div, option, select},
    interfaces::{Element, HtmlOptionElement},
    App, DomFragment, DomView,
};

pub mod components;

mod coprop;
mod euler_approx;
mod explorer;
mod explorer_theta_kappa;
mod ptan;

enum AppState {
    Explorer(explorer::AppState),
    ExplorerThetaKappa(explorer_theta_kappa::AppState),
    EulerApprox(euler_approx::AppState),
    PointTangent(ptan::AppState),
    Coproportional(coprop::AppState),
}

fn app_logic(state: &mut AppState) -> impl DomFragment<AppState> {
    let app = match state {
        AppState::Explorer(state) => OneOf5::A(explorer::app_logic(state).map_state(
            |state: &mut AppState| {
                if let AppState::Explorer(state) = state {
                    state
                } else {
                    unreachable!()
                }
            },
        )),
        AppState::ExplorerThetaKappa(state) => OneOf5::B(
            explorer_theta_kappa::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::ExplorerThetaKappa(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }),
        ),
        AppState::EulerApprox(state) => OneOf5::C(euler_approx::app_logic(state).map_state(
            |state: &mut AppState| {
                if let AppState::EulerApprox(state) = state {
                    state
                } else {
                    unreachable!()
                }
            },
        )),
        AppState::PointTangent(state) => {
            OneOf5::D(ptan::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::PointTangent(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
        AppState::Coproportional(state) => {
            OneOf5::E(coprop::app_logic(state).map_state(|state: &mut AppState| {
                if let AppState::Coproportional(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
    };

    let toolbar = div(select((
        option("Explorer")
            .value("explorer")
            .selected(matches!(state, AppState::Explorer(_))),
        option("Explorer θ⋅κ")
            .value("explorer_theta_kappa")
            .selected(matches!(state, AppState::ExplorerThetaKappa(_))),
        option("Euler Approximation")
            .value("euler_approx")
            .selected(matches!(state, AppState::EulerApprox(_))),
        option("Point · Tangent")
            .value("ptan")
            .selected(matches!(state, AppState::PointTangent(_))),
        option("Coproportional")
            .value("coprop")
            .selected(matches!(state, AppState::Coproportional(_))),
    ))
    .on_change(move |state: &mut AppState, e| {
        match e
            .target()
            .unwrap()
            .unchecked_into::<web_sys::HtmlSelectElement>()
            .value()
            .as_ref()
        {
            "explorer" => *state = AppState::Explorer(explorer::AppState::default()),
            "explorer_theta_kappa" => {
                *state = AppState::ExplorerThetaKappa(explorer_theta_kappa::AppState::default())
            }
            "euler_approx" => *state = AppState::EulerApprox(euler_approx::AppState::default()),
            "ptan" => *state = AppState::PointTangent(ptan::AppState::default()),
            "coprop" => *state = AppState::Coproportional(coprop::AppState::default()),
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
