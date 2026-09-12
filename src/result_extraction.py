"""Extract numerical figure inputs from a solved CISTAR flowsheet."""

from __future__ import annotations

import numbers
from typing import Any

import pandas as pd
from idaes.core.util.tables import arcs_to_stream_dict, create_stream_table_dataframe
from pyomo.environ import value

from src.utility_minimization_1d import gen_curves, heat_ex_data, return_HX_results


def dataframe_payload(frame: pd.DataFrame) -> dict[str, Any]:
    """Return a JSON-safe, orientation-preserving DataFrame representation."""

    def normalize(value_: Any) -> Any:
        if value_ is None:
            return None
        if isinstance(value_, bool):
            return value_
        if isinstance(value_, numbers.Integral):
            return int(value_)
        if isinstance(value_, numbers.Real):
            return None if pd.isna(value_) else float(value_)
        if isinstance(value_, str):
            return value_
        return str(value_)

    return {
        "columns": [str(column) for column in frame.columns],
        "index": [str(index) for index in frame.index],
        "data": [
            [normalize(value_) for value_ in row]
            for row in frame.itertuples(index=False, name=None)
        ],
    }


def collect_figure_data(model: Any) -> dict[str, Any]:
    """Collect the numerical series used by the published analysis figures."""
    heating = [model.fs.H101, model.fs.H103, model.fs.R101]
    cooling = [
        model.fs.H102,
        model.fs.H104,
        model.fs.H105,
        model.fs.H106,
        model.fs.R102,
    ]
    curve_data = heat_ex_data(model.fs, heating, cooling)
    hot_temperature, hot_heat = gen_curves(
        curve_data.Cooling_Tin,
        curve_data.Cooling_Tout,
        curve_data.Cooling_Q,
    )
    cold_temperature, cold_heat = gen_curves(
        curve_data.Heating_Tin,
        curve_data.Heating_Tout,
        -curve_data.Heating_Q,
    )
    cold_heat = cold_heat + sum(curve_data.Heating_Q) + value(curve_data.Qw)

    component_flow = {
        component: value(
            model.fs.F102.liq_outlet.flow_mol_phase_comp[0, "Liq", component]
        )
        for component in model.fs.liquid
    }
    total_liquid_flow = sum(component_flow.values())

    return {
        "upstream_emissions_kg_co2e_per_gj": value(model.fs.upstream_emissions),
        "liquid_product": {
            "total_flow_mol_per_s": total_liquid_flow,
            "component_flow_mol_per_s": component_flow,
            "component_mole_percent": {
                component: 100.0 * flow / total_liquid_flow
                for component, flow in component_flow.items()
            },
            "component_lhv_contribution_mj_per_s": {
                component: value(
                    model.fs.LHV_per_component_per_stream[
                        "liq_outlet", component
                    ]
                )
                * total_liquid_flow
                for component in model.fs.liquid
            },
        },
        "composite_curves": {
            "hot": {
                "temperature_k": hot_temperature.tolist(),
                "cumulative_heat_gj_per_hour": (-hot_heat).tolist(),
            },
            "cold": {
                "temperature_k": cold_temperature.tolist(),
                "cumulative_heat_gj_per_hour": cold_heat.tolist(),
            },
        },
        "heat_exchanger_table": dataframe_payload(
            return_HX_results(model.fs, heating + cooling)
        ),
        "stream_table": dataframe_payload(
            create_stream_table_dataframe(
                arcs_to_stream_dict(model, descend_into=True)
            )
        ),
    }
