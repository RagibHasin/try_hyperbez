// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use wasm_bindgen::prelude::*;
use web_sys::window;
use xilem_web::{
    concurrent::task,
    core::{fork, one_of::OneOf5},
    elements::html::{div, option, select},
    interfaces::{Element, HtmlOptionElement},
    App, DomFragment, DomView,
};

pub mod components;

mod coprop;
mod euler_approx;
mod hyper_theta_kappa;
mod hyperparams;
mod ptan;

enum Explorer {
    HyperParams(hyperparams::AppState),
    ThetaKappa(hyper_theta_kappa::AppState),
    EulerApprox(euler_approx::AppState),
    PointTangent(ptan::AppState),
    Coproportional(coprop::AppState),
}

impl std::str::FromStr for Explorer {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (explorer, params) = s.split_once(';').ok_or("invalid format")?;
        Ok(match explorer {
            "#hyperparams" => Explorer::HyperParams(params.parse::<hyperparams::AppData>()?.into()),
            "#hyper_theta_kappa" => {
                Explorer::ThetaKappa(params.parse::<hyper_theta_kappa::AppData>()?.into())
            }
            "#euler_approx" => {
                Explorer::EulerApprox(params.parse::<euler_approx::AppData>()?.into())
            }
            "#ptan" => Explorer::PointTangent(params.parse::<ptan::AppData>()?.into()),
            "#coprop" => Explorer::Coproportional(params.parse::<coprop::AppData>()?.into()),
            _ => unreachable!(),
        })
    }
}

impl std::fmt::Display for Explorer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Explorer::HyperParams(e) => write!(f, "#hyperparams;{}", e.data),
            Explorer::ThetaKappa(e) => write!(f, "#hyper_theta_kappa;{}", e.data),
            Explorer::EulerApprox(e) => write!(f, "#euler_approx;{}", e.data),
            Explorer::PointTangent(e) => write!(f, "#ptan;{}", e.data),
            Explorer::Coproportional(e) => write!(f, "#coprop;{}", e.data),
        }
    }
}

impl Explorer {
    fn try_update(&mut self, s: &str) -> Result<(), Box<dyn std::error::Error>> {
        let (explorer, params) = s.split_once(';').ok_or("invalid format")?;
        match self {
            Explorer::HyperParams(e) if explorer == "#hyperparams" => e.data = params.parse()?,
            Explorer::ThetaKappa(e) if explorer == "#hyper_theta_kappa" => {
                e.data = params.parse()?
            }
            Explorer::EulerApprox(e) if explorer == "#euler_approx" => e.data = params.parse()?,
            Explorer::PointTangent(e) if explorer == "#ptan" => e.data = params.parse()?,
            Explorer::Coproportional(e) if explorer == "#coprop" => e.data = params.parse()?,
            _ => *self = s.parse()?,
        }
        Ok(())
    }
}

struct AppState {
    explorer: Explorer,
    data_fragment: String,
}

impl From<Explorer> for AppState {
    fn from(explorer: Explorer) -> Self {
        Self {
            data_fragment: explorer.to_string(),
            explorer,
        }
    }
}

fn explorer_app(state: &mut Explorer) -> impl DomView<Explorer> {
    match state {
        Explorer::HyperParams(state) => OneOf5::A(hyperparams::app_logic(state).map_state(
            |state: &mut Explorer| {
                if let Explorer::HyperParams(state) = state {
                    state
                } else {
                    unreachable!()
                }
            },
        )),
        Explorer::ThetaKappa(state) => OneOf5::B(hyper_theta_kappa::app_logic(state).map_state(
            |state: &mut Explorer| {
                if let Explorer::ThetaKappa(state) = state {
                    state
                } else {
                    unreachable!()
                }
            },
        )),
        Explorer::EulerApprox(state) => OneOf5::C(euler_approx::app_logic(state).map_state(
            |state: &mut Explorer| {
                if let Explorer::EulerApprox(state) = state {
                    state
                } else {
                    unreachable!()
                }
            },
        )),
        Explorer::PointTangent(state) => {
            OneOf5::D(ptan::app_logic(state).map_state(|state: &mut Explorer| {
                if let Explorer::PointTangent(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
        Explorer::Coproportional(state) => {
            OneOf5::E(coprop::app_logic(state).map_state(|state: &mut Explorer| {
                if let Explorer::Coproportional(state) = state {
                    state
                } else {
                    unreachable!()
                }
            }))
        }
    }
}

fn app_logic(state: &mut AppState) -> impl DomFragment<AppState> {
    let app = explorer_app(&mut state.explorer)
        .map_state(|state: &mut AppState| &mut state.explorer)
        .map_action(|state: &mut AppState, _| {
            state.data_fragment = state.explorer.to_string();
            window()
                .unwrap_throw()
                .history()
                .unwrap_throw()
                .replace_state_with_url(&JsValue::NULL, "", Some(&state.data_fragment))
                .ok();
        });

    let toolbar = div(select((
        option("HyperParams")
            .value("hyperparams")
            .selected(matches!(state.explorer, Explorer::HyperParams(_))),
        option("θ⋅κ → HyperParams")
            .value("hyper_theta_kappa")
            .selected(matches!(state.explorer, Explorer::ThetaKappa(_))),
        option("Euler Approximation")
            .value("euler_approx")
            .selected(matches!(state.explorer, Explorer::EulerApprox(_))),
        option("Point · Tangent")
            .value("ptan")
            .selected(matches!(state.explorer, Explorer::PointTangent(_))),
        option("Coproportional")
            .value("coprop")
            .selected(matches!(state.explorer, Explorer::Coproportional(_))),
    ))
    .on_change(move |state: &mut AppState, e| {
        *state = AppState::from(
            match e
                .target()
                .unwrap_throw()
                .unchecked_into::<web_sys::HtmlSelectElement>()
                .value()
                .as_ref()
            {
                "hyperparams" => Explorer::HyperParams(hyperparams::AppState::default()),
                "hyper_theta_kappa" => Explorer::ThetaKappa(hyper_theta_kappa::AppState::default()),
                "euler_approx" => Explorer::EulerApprox(euler_approx::AppState::default()),
                "ptan" => Explorer::PointTangent(ptan::AppState::default()),
                "coprop" => Explorer::Coproportional(coprop::AppState::default()),
                _ => unreachable!(),
            },
        )
    }))
    .id("toolbar");

    (
        fork(
            app,
            task(
                |proxy, _| async {
                    let callback: Closure<dyn Fn(web_sys::Event)> =
                        Closure::new(move |_| proxy.send_message(()));
                    window()
                        .unwrap_throw()
                        .set_onhashchange(Some(callback.as_ref().unchecked_ref()));
                    std::mem::forget(callback);
                },
                |state: &mut AppState, _: ()| {
                    if let Ok(url_fragment) = window().unwrap_throw().location().hash() {
                        if state.data_fragment != url_fragment
                            && state.explorer.try_update(&url_fragment).is_ok()
                        {
                            state.data_fragment = url_fragment;
                        }
                    }
                },
            ),
        ),
        toolbar,
    )
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

    let app_state = AppState::from(
        window()
            .unwrap_throw()
            .location()
            .hash()
            .ok()
            .as_deref()
            .and_then(|f| f.parse().ok())
            .unwrap_or_else(|| Explorer::PointTangent(Default::default())),
    );

    App::new(
        xilem_web::document_body(),
        app_state,
        // AppState::from(Explorer::PointTangent(Default::default())),
        app_logic,
    )
    .run();
}

#[allow(unused_must_use)]
#[cfg(test)]
mod tests {
    use super::*;
    use test_log::test;

    #[test]
    #[test_log(default_log_filter = "trace")]
    fn test_explorer_from_str() {
        if let Err(e) = "#ptan;50,240,sym".parse::<Explorer>() {
            panic!("{e}");
        }
    }
}
