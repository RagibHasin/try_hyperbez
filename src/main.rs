// Copyright 2024 Muhammad Ragib Hasin
// SPDX-License-Identifier: Apache-2.0

use wasm_bindgen::prelude::*;
use web_sys::window;
use xilem_web::{
    App, DomFragment, DomView,
    concurrent::task,
    core::{View, fork, one_of::OneOf7},
    elements::html::{button, div, option, select},
    interfaces::{Element, HtmlOptionElement},
};

pub mod common;
pub mod components;

mod coprop;
mod euler_approx;
mod hyper_theta_kappa;
mod hyperparams;
mod ladder_alt;
mod ptan;
mod q0q1;

enum Explorer {
    HyperParams(hyperparams::AppState),
    ThetaKappa(hyper_theta_kappa::AppState),
    EulerApprox(euler_approx::AppState),
    PointTangent(ptan::AppState),
    Coproportional(coprop::AppState),
    Q0Q1(q0q1::AppState),
    LadderAlt(ladder_alt::AppState),
}

impl std::str::FromStr for Explorer {
    type Err = Box<dyn std::error::Error>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (tag, params) = s.split_once(';').ok_or("invalid format")?;
        Ok(match tag {
            "#hyperparams" => Explorer::HyperParams(params.parse::<hyperparams::AppData>()?.into()),
            "#hyper_theta_kappa" => {
                Explorer::ThetaKappa(params.parse::<hyper_theta_kappa::AppData>()?.into())
            }
            "#euler_approx" => {
                Explorer::EulerApprox(params.parse::<euler_approx::AppData>()?.into())
            }
            "#ptan" => Explorer::PointTangent(params.parse::<ptan::AppData>()?.into()),
            "#coprop" => Explorer::Coproportional(params.parse::<coprop::AppData>()?.into()),
            "#q0q1" => Explorer::Q0Q1(params.parse::<q0q1::AppData>()?.into()),
            "#ladder_alt" => Explorer::LadderAlt(params.parse::<ladder_alt::AppData>()?.into()),
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
            Explorer::Q0Q1(e) => write!(f, "#q0q1;{}", e.data),
            Explorer::LadderAlt(e) => write!(f, "#ladder_alt;{}", e.data),
        }
    }
}

macro_rules! state_mapper {
    ($var:ident) => {
        |state: &mut Explorer| {
            if let Explorer::$var(state) = state {
                state
            } else {
                unreachable!()
            }
        }
    };
}

impl Explorer {
    fn try_update(&mut self, s: &str) -> Result<(), Box<dyn std::error::Error>> {
        let (tag, params) = s.split_once(';').ok_or("invalid format")?;
        match self {
            Explorer::HyperParams(e) if tag == "#hyperparams" => e.data = params.parse()?,
            Explorer::ThetaKappa(e) if tag == "#hyper_theta_kappa" => e.data = params.parse()?,
            Explorer::EulerApprox(e) if tag == "#euler_approx" => e.data = params.parse()?,
            Explorer::PointTangent(e) if tag == "#ptan" => e.data = params.parse()?,
            Explorer::Coproportional(e) if tag == "#coprop" => e.data = params.parse()?,
            Explorer::Q0Q1(e) if tag == "#q0q1" => e.data = params.parse()?,
            Explorer::LadderAlt(e) if tag == "#ladder_alt" => e.data = params.parse()?,
            _ => *self = s.parse()?,
        }
        Ok(())
    }

    fn refresh(&mut self) {
        match self {
            Explorer::HyperParams(e) => e.data = hyperparams::AppData::default(),
            Explorer::ThetaKappa(e) => e.data = hyper_theta_kappa::AppData::default(),
            Explorer::EulerApprox(e) => e.data = euler_approx::AppData::default(),
            Explorer::PointTangent(e) => e.data = ptan::AppData::default(),
            Explorer::Coproportional(e) => e.data = coprop::AppData::default(),
            Explorer::Q0Q1(e) => e.data = q0q1::AppData::default(),
            Explorer::LadderAlt(e) => e.data = ladder_alt::AppData::default(),
        }
    }

    fn view(&mut self) -> impl DomView<Self> + use<> {
        match self {
            Explorer::HyperParams(state) => {
                OneOf7::A(hyperparams::app_logic(state).map_state(state_mapper!(HyperParams)))
            }
            Explorer::ThetaKappa(state) => {
                OneOf7::B(hyper_theta_kappa::app_logic(state).map_state(state_mapper!(ThetaKappa)))
            }
            Explorer::EulerApprox(state) => {
                OneOf7::C(euler_approx::app_logic(state).map_state(state_mapper!(EulerApprox)))
            }
            Explorer::PointTangent(state) => {
                OneOf7::D(ptan::app_logic(state).map_state(state_mapper!(PointTangent)))
            }
            Explorer::Coproportional(state) => {
                OneOf7::E(coprop::app_logic(state).map_state(state_mapper!(Coproportional)))
            }
            Explorer::Q0Q1(state) => {
                OneOf7::F(q0q1::app_logic(state).map_state(state_mapper!(Q0Q1)))
            }
            Explorer::LadderAlt(state) => {
                OneOf7::G(ladder_alt::app_logic(state).map_state(state_mapper!(LadderAlt)))
            }
        }
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

fn app_logic(state: &mut AppState) -> impl DomFragment<AppState> + use<> {
    let app = state
        .explorer
        .view()
        .map_state(|state: &mut AppState| &mut state.explorer)
        .map_message_result(|state: &mut AppState, r| {
            state.data_fragment = state.explorer.to_string();
            window()
                .unwrap_throw()
                .history()
                .unwrap_throw()
                .replace_state_with_url(&JsValue::NULL, "", Some(&state.data_fragment))
                .ok();
            r
        });

    let toolbar = div((
        select((
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
            option("q₀ q₁")
                .value("q0q1")
                .selected(matches!(state.explorer, Explorer::Q0Q1(_))),
            option("Ladder (Alt)")
                .value("ladder_alt")
                .selected(matches!(state.explorer, Explorer::LadderAlt(_))),
        ))
        .on_change(|state: &mut AppState, e| {
            *state = AppState::from(
                match e
                    .target()
                    .unwrap_throw()
                    .unchecked_into::<web_sys::HtmlSelectElement>()
                    .value()
                    .as_str()
                {
                    "hyperparams" => Explorer::HyperParams(hyperparams::AppState::default()),
                    "hyper_theta_kappa" => {
                        Explorer::ThetaKappa(hyper_theta_kappa::AppState::default())
                    }
                    "euler_approx" => Explorer::EulerApprox(euler_approx::AppState::default()),
                    "ptan" => Explorer::PointTangent(ptan::AppState::default()),
                    "coprop" => Explorer::Coproportional(coprop::AppState::default()),
                    "q0q1" => Explorer::Q0Q1(q0q1::AppState::default()),
                    "ladder_alt" => Explorer::LadderAlt(ladder_alt::AppState::default()),
                    _ => unreachable!(),
                },
            )
        }),
        button("🔄︎").on_click(|state: &mut AppState, _| {
            state.explorer.refresh();
            state.data_fragment = state.explorer.to_string();
        }),
    ))
    .id("toolbar");

    (
        fork(
            app,
            task::<_, _, _, AppState, _, _>(
                |proxy, _| async {
                    let callback: Closure<dyn Fn(web_sys::Event)> =
                        Closure::new(move |_| proxy.send_message(()));
                    window()
                        .unwrap_throw()
                        .set_onhashchange(Some(callback.as_ref().unchecked_ref()));
                    std::mem::forget(callback);
                },
                |state: &mut AppState, ()| {
                    if let Ok(url_fragment) = window().unwrap_throw().location().hash()
                        && state.data_fragment != url_fragment
                        && state.explorer.try_update(&url_fragment).is_ok()
                    {
                        state.data_fragment = url_fragment;
                    }
                },
            ),
        ),
        toolbar,
    )
}

pub fn main() {
    use tracing_subscriber::{fmt::format::Pretty, prelude::*};
    use tracing_web::{MakeWebConsoleWriter, performance_layer};

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

    App::new(xilem_web::document_body(), app_state, app_logic).run();
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
