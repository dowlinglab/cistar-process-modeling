# adding grandparent directory to path
import sys
sys.path.append("./../")

from pyomo.environ import (Constraint,
                           ConstraintList,
                           Var,
                           ConcreteModel,
                           Expression,
                           Param,
                           Set,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value,
                           minimize)

from pyomo.network import Arc

# Import plotting functions
import matplotlib.pyplot as plt

# Import numpy library 
import numpy as np

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock, MaterialFlowBasis

# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale

from pyomo.environ import units as pyunits

from idaes.models.properties.modular_properties.state_definitions import FcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal

from idaes.models.properties.modular_properties.reactions.equilibrium_constant import \
    van_t_hoff
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import \
    power_law_equil

# First, create a thermophsyical property definition that will be used
# with the reactions
import idaes.models.properties.modular_properties.pure.Perrys as Perrys

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom, large_residuals_set
from idaes.core.util.initialization import propagate_state

# Import the Generic Parameter Block
from idaes.models.properties.modular_properties.base.generic_property import (
        GenericParameterBlock)

from idaes.models.properties.modular_properties.base.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm

# Import compressor unit model from the model library
from idaes.models.unit_models.pressure_changer import PressureChanger,ThermodynamicAssumption

# from packed_bed_reactor import PBR

from idaes.models.unit_models import (
                               Heater,
                               StoichiometricReactor,
                               Translator,
                               Flash,
                               Mixer,
                               Separator,
                               PFR)

# Import compressor unit model from the model library
from idaes.models.unit_models.pressure_changer import PressureChanger,ThermodynamicAssumption
from idaes.models.unit_models.separator import SplittingType

# Import python path
import os

# Import idaes model serializer to store initialized model
from idaes.core.util import model_serializer as ms

import warnings

# Import idaes logger to set output levels
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale

from idaes.core.util.initialization import propagate_state

# Stream component properties
from src.state_properties.properties_vap import configuration as configuration_vapor
from src.state_properties.properties_vap_H2_permeate import configuration as configuration_H2_permeate
from src.state_properties.properties_VLE_FpcTP import configuration as configuration_VLE
from src.state_properties.properties_vap_post_flash_2 import configuration as configuration_vap_post_flash_2

# dehydrogenation simple kinetics reaction block
from src.dehydro_reactions import rxn_configuration

# oligomerization kinetics reaction block
from src.reaction_network_generator import return_reaction_network

# Emissions and costing functions
from src.emissions_calculations import calc_lhv_values, calculate_stream_energies, calculate_emissions, create_ghg_objective
from src.costing_function import add_costing,calculate_costs_for_objective
from src.utility_minimization_1d import (
    min_utility,
    PinchDataClass,
    heat_ex_data,
    gen_curves,
    print_HX_results,
    generate_curves,
    heat_data,
    pinch_calc,
    return_data
)

import pandas as pd

from src.utility_minimization_1d import heat_ex_data, return_HX_results
from src.plotting_functions import generate_and_save_composite_curves

def create_flowsheet(model_code):
    """
    Function to create and return ConcreteModel flowsheet
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    create_properties(m,model_code)
    return m


def create_properties(m,model_code):
    """
    Function to create stream properties
    Arguents:
        m: int
        Integer to identify ROK model. Possible values: 2,3,4,5
    """
    rate_reactions = return_reaction_network(model_code = model_code)
    # Add properties parameter blocks to the flowsheet with specifications

    # general vapor properties
    m.fs.Vap_props = GenericParameterBlock(**configuration_vapor)
    # vapor properties for components in H2-membrane permeate (H2, CH4, C2H6, C2H4)
    m.fs.H2_permeate_props = GenericParameterBlock(**configuration_H2_permeate)
    # general VLE properties
    m.fs.VLE_props = GenericParameterBlock(**configuration_VLE)
    # Vapor properties without H2 and CH4 for second flash vapor
    m.fs.Vap_props_post_flash_2 = GenericParameterBlock(**configuration_vap_post_flash_2)

    m.fs.liquid = Set(initialize=['ethane', 'propane', 'nbutane', 'ibutane', 'pentane', 'hexane','heptane',
                    'octane','ethylene', 'propene', 'butene', 'pentene', 'hexene', 'heptene', 'octene', 'nonene'])

    m.fs.vapor = Set(initialize=['hydrogen','methane','ethane', 'propane', 'nbutane', 'ibutane', 'pentane',
                    'hexane','heptane','octane', 'ethylene', 'propene', 'butene', 'pentene', 'hexene', 
                    'heptene', 'octene', 'nonene'])

    m.fs.vapor_post_flash_2 = Set(initialize=['ethane', 'propane', 'nbutane', 'ibutane', 'pentane', 'hexane', 
                                'heptane', 'octane', 'ethylene', 'propene', 'butene', 'pentene', 'hexene', 
                                'heptene', 'octene', 'nonene'])

    m.fs.vapor_H2_permeate = Set(initialize=['hydrogen','methane','ethane','ethylene'])

    m.fs.rxn_params = GenericReactionParameterBlock(
                        property_package= m.fs.Vap_props,
                        **rxn_configuration)

    m.fs.reactions = GenericReactionParameterBlock(
                    property_package= m.fs.Vap_props, 
                    base_units= {"time": pyunits.s,
                                "length": pyunits.m,
                                "mass": pyunits.kg,
                                "amount": pyunits.mol,
                                "temperature": pyunits.K},
                    rate_reactions= rate_reactions,
                    reaction_basis= MaterialFlowBasis.molar)


def define_models(m, catalyst_mass = 1167.003367):
    """
    Function to define unit models and properties
    Arguments:
        catalyst_mass: float
            Mass of catalyst loading in oligomerization reactor, kg
    """
    m.fs.M101 = Mixer(
                property_package=m.fs.Vap_props,
                inlet_list=["feed", "recycle"])

    m.fs.H101 = Heater(
                property_package=m.fs.Vap_props,
                has_pressure_change=True,
                has_phase_equilibrium=True)

    m.fs.R101 = StoichiometricReactor(
                property_package=m.fs.Vap_props,
                reaction_package=m.fs.rxn_params,
                has_heat_of_reaction=True,
                has_heat_transfer=True,
                has_pressure_change=False)

    m.fs.H102 = Heater(
                property_package=m.fs.Vap_props,
                has_pressure_change=False,
                has_phase_equilibrium=True)

    m.fs.S101 = Separator(
                property_package=m.fs.Vap_props,
                ideal_separation=False,
                outlet_list=["H2_purge", "outlet"],
                split_basis=SplittingType.componentFlow)
    
    m.fs.T101 = Translator(
                inlet_property_package=m.fs.Vap_props,
                outlet_property_package=m.fs.H2_permeate_props,
                has_phase_equilibrium=False,
                outlet_state_defined=True)

    m.fs.E101 = PressureChanger(
                dynamic=False,
                property_package=m.fs.H2_permeate_props,
                compressor=False,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic)

    m.fs.H103 = Heater(
                property_package=m.fs.Vap_props,
                has_pressure_change=True,
                has_phase_equilibrium=True)

    m.fs.R102 = PFR(
                property_package=m.fs.Vap_props,
                reaction_package=m.fs.reactions,
                has_equilibrium_reactions=False,
                has_heat_transfer=True,
                has_heat_of_reaction=False,
                has_pressure_change=False,
                )

    m.fs.H104 = Heater(
                property_package=m.fs.Vap_props,
                has_pressure_change=False,
                has_phase_equilibrium=True)

    m.fs.T102 = Translator(
                inlet_property_package=m.fs.Vap_props,
                outlet_property_package=m.fs.VLE_props,
                has_phase_equilibrium=False,
                outlet_state_defined=True)
    
    m.fs.H105 = Heater(
                property_package=m.fs.VLE_props,
                has_pressure_change=False,
                has_phase_equilibrium=True)

    m.fs.F101 = Flash(
                property_package=m.fs.VLE_props,
                has_heat_transfer=True,
                has_pressure_change=True)
    
    m.fs.H106 = Heater(
                property_package=m.fs.VLE_props,
                has_pressure_change=True,
                has_phase_equilibrium=True)
    
    m.fs.F102 = Flash(
                property_package=m.fs.VLE_props,
                has_heat_transfer=True,
                has_pressure_change=True)

    m.fs.T103 = Translator(
                inlet_property_package=m.fs.VLE_props,
                outlet_property_package=m.fs.Vap_props,
                has_phase_equilibrium=False,
                outlet_state_defined=True)

    m.fs.T104 = Translator(
                inlet_property_package=m.fs.VLE_props,
                outlet_property_package=m.fs.Vap_props_post_flash_2,
                has_phase_equilibrium=False,
                outlet_state_defined=True)

    m.fs.C101 = PressureChanger(
                dynamic=False,
                property_package=m.fs.Vap_props_post_flash_2,
                compressor=True,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic)
    
    m.fs.T105 = Translator(
                inlet_property_package=m.fs.Vap_props_post_flash_2,
                outlet_property_package=m.fs.Vap_props,
                has_phase_equilibrium=False,
                outlet_state_defined=True)
    
    m.fs.M102 = Mixer(
                property_package=m.fs.Vap_props,
                inlet_list=["vapor_recycle_F101", "vapor_recycle_F102"])

    m.fs.S102 = Separator(
                property_package=m.fs.Vap_props,
                ideal_separation=False,
                outlet_list=["purge", "recycle"])

    m.fs.C102 = PressureChanger(
                dynamic=False,
                property_package=m.fs.Vap_props,
                compressor=True,
                thermodynamic_assumption=ThermodynamicAssumption.isentropic)


def define_arcs(m):
    """
    Function to define connections using arcs
    """
    m.fs.s01 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s02 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s03 = Arc(source=m.fs.R101.outlet, destination=m.fs.H102.inlet)
    m.fs.s04 = Arc(source=m.fs.H102.outlet, destination=m.fs.S101.inlet)
    m.fs.s05 = Arc(source=m.fs.S101.outlet, destination=m.fs.H103.inlet)
    m.fs.s06 = Arc(source=m.fs.S101.H2_purge, destination=m.fs.T101.inlet)
    m.fs.s07 = Arc(source=m.fs.T101.outlet, destination=m.fs.E101.inlet)
    m.fs.s08 = Arc(source=m.fs.H103.outlet, destination=m.fs.R102.inlet)
    m.fs.s09 = Arc(source=m.fs.R102.outlet, destination=m.fs.H104.inlet)
    m.fs.s10 = Arc(source=m.fs.H104.outlet, destination=m.fs.T102.inlet)

    m.fs.s11 = Arc(source=m.fs.T102.outlet, destination=m.fs.H105.inlet)
    m.fs.s12 = Arc(source=m.fs.H105.outlet, destination=m.fs.F101.inlet)
    m.fs.s13 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.H106.inlet)
    m.fs.s14 = Arc(source=m.fs.H106.outlet, destination=m.fs.F102.inlet)
    
    m.fs.s15 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.T103.inlet)
    m.fs.s16 = Arc(source=m.fs.F102.vap_outlet, destination=m.fs.T104.inlet)
    m.fs.s17 = Arc(source=m.fs.T103.outlet, destination=m.fs.M102.vapor_recycle_F101)
    m.fs.s18 = Arc(source=m.fs.T104.outlet, destination=m.fs.C101.inlet)
    m.fs.s19 = Arc(source=m.fs.C101.outlet, destination=m.fs.T105.inlet)
    m.fs.s20 = Arc(source=m.fs.T105.outlet, destination=m.fs.M102.vapor_recycle_F102)
    m.fs.s21 = Arc(source=m.fs.M102.outlet, destination=m.fs.S102.inlet)
    m.fs.s22 = Arc(source=m.fs.S102.recycle, destination=m.fs.C102.inlet)
    m.fs.s23 = Arc(source=m.fs.C102.outlet, destination=m.fs.M101.recycle)


## Unit model constraints and variable fixing
def set_unit_model_variables(m,
                            model_code,
                            feed_flow_rate = 481.3888889, 
                            feed_temp = 308.0,
                            feed_pressure = 700000.0,
                            inlet_composition_dict={},
                            dehydro_conv_dict= {'ethane':0.3566,
                                                'propane':0.6632,
                                                'nbutane':0.5188}):
    """
    Function to set constraints and variable values for unit models
    Arguments:
        feed_flow_rate: float
            Fresh feed flow rate, mol/s
        feed_temp: float
            Fresh feed temperature, K
        feed_pressure: float
            Fresh feed pressure, Pa
        inlet_composition_dict={}
            Fresh feed composition in mole fraction
        dehydro_conv_dict: dict
            Conversions for ethane, propane, and nbutane
    """
    # Mixer M101 to mix recycle and feed
    m.fs.M101.feed.flow_mol.fix(feed_flow_rate)
    for k,v in inlet_composition_dict.items():
        m.fs.M101.feed.mole_frac_comp[0, k].fix(v)
    m.fs.M101.feed.temperature.fix(feed_temp)
    m.fs.M101.feed.pressure.fix(feed_pressure)

    for k,v in inlet_composition_dict.items():
        m.fs.M101.recycle.mole_frac_comp[0, k].unfix()


    ### Heat exchanger H101
    m.fs.H101.outlet.temperature.fix(1073.0)
    m.fs.H101.outlet.pressure.fix(600000.0)
    # H101 is a heater, heat_duty > 0
    m.fs.H101.heat_duty.setlb(0.0)


    ### Dehydrogenation reactor R101
    m.fs.R101.C2H6conversion = Var(initialize=dehydro_conv_dict['ethane'], bounds=(0, 1))
    m.fs.R101.C2H6conv_constraint = Constraint(
        expr=m.fs.R101.C2H6conversion*m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "ethane"] ==
        (m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "ethane"] -
        m.fs.R101.outlet.flow_mol[0]*
        m.fs.R101.outlet.mole_frac_comp[0, "ethane"]))

    m.fs.R101.C3H8conversion = Var(initialize=dehydro_conv_dict['propane'], bounds=(0, 1))
    m.fs.R101.C3H8conv_constraint = Constraint(
        expr=m.fs.R101.C3H8conversion*m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "propane"] ==
        (m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "propane"] -
        m.fs.R101.outlet.flow_mol[0]*
        m.fs.R101.outlet.mole_frac_comp[0, "propane"]))

    m.fs.R101.nC4H10conversion = Var(initialize=dehydro_conv_dict['nbutane'], bounds=(0, 1))
    m.fs.R101.nC4H10conv_constraint = Constraint(
        expr=m.fs.R101.nC4H10conversion*m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "nbutane"] ==
        (m.fs.R101.inlet.flow_mol[0]*
        m.fs.R101.inlet.mole_frac_comp[0, "nbutane"] -
        m.fs.R101.outlet.flow_mol[0]*
        m.fs.R101.outlet.mole_frac_comp[0, "nbutane"]))

    m.fs.R101.C2H6conversion.fix(dehydro_conv_dict['ethane'])
    m.fs.R101.C3H8conversion.fix(dehydro_conv_dict['propane'])
    m.fs.R101.nC4H10conversion.fix(dehydro_conv_dict['nbutane'])
    m.fs.R101.outlet.temperature.fix(1073.0)


    ### Heat exchanger H102
    m.fs.H102.outlet.temperature.fix(473.0)
    # H102 is a cooler, heat_duty < 0
    m.fs.H102.heat_duty.setub(0.0)


    ### Hydrogen Separation S101
    m.fs.S101.H2_recovery_frac = Param(initialize=0.57565,mutable=True)
    m.fs.S101.coproduct_leak_wrt_H2 = Param(initialize=7.5,mutable=True)

    m.fs.S101.H2_purge_split_fraction_constraint = ConstraintList()
    for i in m.fs.vapor:
        if i == 'hydrogen':
            m.fs.S101.H2_purge_split_fraction_constraint.add(expr = m.fs.S101.split_fraction[0.0, 'H2_purge', i] == m.fs.S101.H2_recovery_frac)
        elif i in ['methane','ethane','ethylene']:
            m.fs.S101.H2_purge_split_fraction_constraint.add(expr = m.fs.S101.split_fraction[0.0, 'H2_purge', i] * m.fs.S101.coproduct_leak_wrt_H2 == m.fs.S101.H2_recovery_frac)
        else:
            m.fs.S101.split_fraction[0.0, 'H2_purge', i].fix(1e-8)

    
    ### General Vapor to H2-permeate translator T101
    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T101.eq_temperature = Constraint(
        expr=m.fs.T101.outlet.temperature[0] ==
        m.fs.T101.inlet.temperature[0])

    # Add constraint: Outlet pressure = Inlet pressure
    m.fs.T101.eq_pressure = Constraint(
        expr=m.fs.T101.outlet.pressure[0] ==
        m.fs.T101.inlet.pressure[0])
    
    # Add constraint: Outlet flow rate = Inlet flow rate
    m.fs.T101.eq_flow = Constraint(
        expr=m.fs.T101.outlet.flow_mol[0] == 
        m.fs.T101.inlet.flow_mol[0])

    def return_initial_vap_constraint_T101(blk,i):
        return blk.outlet.flow_mol[0] * blk.outlet.mole_frac_comp[0, i] == blk.inlet.flow_mol[0] * blk.inlet.mole_frac_comp[0, i]
    m.fs.T101.eq_vap_frac_constraint_T101 = Constraint(m.fs.vapor_H2_permeate, rule=return_initial_vap_constraint_T101)
            

    ### Expander E101 to de-pressurize H2 outlet
    m.fs.E101.outlet.pressure.fix(100000.0)
    m.fs.E101.efficiency_isentropic.fix(0.85)
    

    ### Heat exchanger H103 to pressurize oligomerization reactor feed
    m.fs.H103.outlet.temperature.fix(573.0)
    m.fs.H103.outlet.pressure.fix(500000.0)
    # H103 needs heat, heat_duty > 0
    m.fs.H103.heat_duty.setlb(0.0)


    ### Oligomerization Reactor R102
    m.fs.R102.length.fix(1.0)
    m.fs.R102.volume.fix(m.fs.reactions.catalyst_mass/m.fs.reactions.bulk_density)

    for i in m.fs.R102.control_volume.length_domain:
        if i != 0:
            m.fs.R102.control_volume.properties[0.0,i].temperature.fix(m.fs.H103.outlet.temperature[0])


    ### Heat Exchanger H104
    m.fs.H104.outlet.temperature.fix(308.0)
    # H104 is a cooler, heat_duty < 0
    m.fs.H104.heat_duty.setub(0.0)
    

    ### Vapor to Vapor-Liquid translator T102
    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T102.eq_temperature = Constraint(
        expr=m.fs.T102.outlet.temperature[0] ==
        m.fs.T102.inlet.temperature[0])

    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T102.eq_pressure = Constraint(
        expr=m.fs.T102.outlet.pressure[0] ==
        m.fs.T102.inlet.pressure[0])

    def return_vap_outlet_constraint_T102(blk,i):
        return blk.outlet.flow_mol_phase_comp[0, 'Vap', i] == blk.inlet.flow_mol[0] * blk.inlet.mole_frac_comp[0, i]
    m.fs.T102.eq_vap_frac_constraint_T102 = Constraint(m.fs.vapor, rule=return_vap_outlet_constraint_T102)

    def return_liq_outlet_constraint_T102(blk,i):
        return blk.outlet.flow_mol_phase_comp[0,'Liq', i] == 1e-8
    m.fs.T102.eq_liq_flow_constraint_T102 = Constraint(m.fs.liquid, rule=return_liq_outlet_constraint_T102)


    ### Phase changer H105
    def return_liq_constraint_H105(blk):
        return sum(blk.outlet.flow_mol_phase_comp[0, 'Liq', i] for i in m.fs.liquid) >= 10.0
    m.fs.H105.eq_liq_frac_const = Constraint(rule=return_liq_constraint_H105)
    def return_out_liq_constraint_H105_lb(blk,i):
        return blk.outlet.flow_mol_phase_comp[0, 'Liq', i] >= 1e-4
    m.fs.H105.out_liq_constraint_H105_lb = Constraint(m.fs.liquid,expr = return_out_liq_constraint_H105_lb)
    # H105 temperature flexibility (- 2C)
    m.fs.H105.temperature_constraint_ub = Constraint(expr = m.fs.H105.outlet.temperature[0] <= m.fs.H105.inlet.temperature[0])
    m.fs.H105.temperature_constraint_lb = Constraint(expr = m.fs.H105.outlet.temperature[0] >= m.fs.H105.inlet.temperature[0]-2.0)
    # H105 is a cooler, heat_duty < 0
    m.fs.H105.heat_duty.setub(0.0)


    ### Flash column F101
    m.fs.F101.deltaP.setub(0.0)
    m.fs.F101.temperature_constraint_ub = Constraint(expr = m.fs.F101.control_volume.properties_out[0.0].temperature <= m.fs.F101.control_volume.properties_in[0.0].temperature)
    m.fs.F101.heat_duty.fix(0.0)


    ### Cooler H106 for expansion and cooling
    # constraints to ensure 2-phase flow exiting the expander valve
    # outlet vapor flow sum constraint
    def return_out_vap_sum_constraint_H106(blk):
        return sum(blk.outlet.flow_mol_phase_comp[0, 'Vap', i] for i in m.fs.vapor) >= 1.0
    m.fs.H106.eq_vap_frac_const = Constraint(rule=return_out_vap_sum_constraint_H106)
    # outlet vapor flow lower bounds
    def return_out_vap_constraint_H106_lb(blk,i):
        return blk.outlet.flow_mol_phase_comp[0, 'Vap', i] >= blk.inlet.flow_mol_phase_comp[0, 'Vap', i]
    m.fs.H106.out_vap_constraint_H106_lb = Constraint(m.fs.vapor,expr = return_out_vap_constraint_H106_lb)
    m.fs.H106.outlet.pressure.setlb(100000.0)
    m.fs.H106.outlet.pressure.setub(250000.0)
    m.fs.H106.outlet_temperature_ub = Constraint(expr=m.fs.H106.outlet.temperature[0] <= m.fs.H106.inlet.temperature[0])
    m.fs.H106.outlet.temperature.fix(298.0)
    m.fs.H106.outlet.temperature.setlb(295.0)
    # H106 is a cooler, heat_duty < 0
    m.fs.H106.heat_duty.setub(0.0)


    ### Flash column F102
    # outlet vapor flow sum constraint
    # outlet vapor flow lower and upper bounds
    def return_out_vap_constraint_F102_lb(blk,i):
        return blk.vap_outlet.flow_mol_phase_comp[0, 'Vap', i] >= blk.inlet.flow_mol_phase_comp[0, 'Vap', i]
    m.fs.F102.out_vap_constraint_F102_lb = Constraint(m.fs.vapor,expr = return_out_vap_constraint_F102_lb)
    m.fs.F102.deltaP.setub(0.0)
    m.fs.F102.pressure_drop_constraint_ub = Constraint(expr = m.fs.F102.control_volume.properties_out[0.0].pressure <= 150000.0)
    m.fs.F102.pressure_drop_constraint_lb = Constraint(expr = m.fs.F102.control_volume.properties_out[0.0].pressure >= 100000.0)
    m.fs.F102.outlet_temperature_lb = Constraint(expr = m.fs.F102.control_volume.properties_out[0.0].temperature >= 290.0)
    m.fs.F102.heat_duty.fix(0.0)


    ### Translator T103 for vapor outlet from F102
    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T103.eq_temperature_T103 = Constraint(
        expr=m.fs.T103.outlet.temperature[0] ==
        m.fs.T103.inlet.temperature[0])

    # Add constraint: Outlet pressure = Inlet pressure
    m.fs.T103.eq_pressure_T103 = Constraint(
        expr=m.fs.T103.outlet.pressure[0] ==
        m.fs.T103.inlet.pressure[0])

    m.fs.T103.flow_mol_constraint_T103 = Constraint(
        expr=m.fs.T103.outlet.flow_mol[0] == sum(m.fs.T103.inlet.flow_mol_phase_comp[0, 'Vap', i] for i in m.fs.vapor))
    def return_initial_vap_constraint_T103(blk,i):
        return blk.outlet.flow_mol[0] * blk.outlet.mole_frac_comp[0, i] == blk.inlet.flow_mol_phase_comp[0, 'Vap', i]
    m.fs.T103.eq_vap_frac_constraint_T103 = Constraint(m.fs.vapor, rule=return_initial_vap_constraint_T103)


    ### Translator T104 for vapor outlet from F102
    # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T104.eq_temperature_T104 = Constraint(
        expr=m.fs.T104.outlet.temperature[0] ==
        m.fs.T104.inlet.temperature[0])

    # Add constraint: Outlet pressure = Inlet pressure
    m.fs.T104.eq_pressure_T104 = Constraint(
        expr=m.fs.T104.outlet.pressure[0] ==
        m.fs.T104.inlet.pressure[0])

    m.fs.T104.flow_mol_constraint_T104 = Constraint(
        expr=m.fs.T104.outlet.flow_mol[0] == sum(m.fs.T104.inlet.flow_mol_phase_comp[0, 'Vap', i] for i in m.fs.vapor_post_flash_2))
    def return_initial_vap_constraint_T104(blk,i):
        return blk.outlet.flow_mol[0] * blk.outlet.mole_frac_comp[0, i] == blk.inlet.flow_mol_phase_comp[0, 'Vap', i]
    m.fs.T104.eq_vap_frac_constraint_T104 = Constraint(m.fs.vapor_post_flash_2, rule=return_initial_vap_constraint_T104)


    ### Compressor C101 to pressurize F102 vapor stream
    m.fs.C101.outlet.pressure.fix(400000.0)
    m.fs.C101.efficiency_isentropic.fix(0.85)


    ### Translator T105 for general vapor outlet
        # Add constraint: Outlet temperature = Inlet temperature
    m.fs.T105.eq_temperature = Constraint(
        expr=m.fs.T105.outlet.temperature[0] ==
        m.fs.T105.inlet.temperature[0])

    # Add constraint: Outlet pressure = Inlet pressure
    m.fs.T105.eq_pressure = Constraint(
        expr=m.fs.T105.outlet.pressure[0] ==
        m.fs.T105.inlet.pressure[0])
    
    # Add constraint: Outlet flow rate = Inlet flow rate
    m.fs.T105.eq_flow = Constraint(
        expr=m.fs.T105.outlet.flow_mol[0] == 
        m.fs.T105.inlet.flow_mol[0])

    def return_initial_vap_constraint_T105(blk,i):
        if i in ['hydrogen','methane']:
            return blk.outlet.mole_frac_comp[0, i] == 1e-8
        else:
            return blk.outlet.flow_mol[0] * blk.outlet.mole_frac_comp[0, i] == blk.inlet.flow_mol[0] * blk.inlet.mole_frac_comp[0, i]
    m.fs.T105.eq_vap_frac_constraint_T105 = Constraint(m.fs.vapor, rule=return_initial_vap_constraint_T105)


    ### Splitter S102 to purge from recycle stream
    m.fs.S102.split_fraction[0, "purge"].fix(0.01)


    ### Compressor C102 to pressurize recycle stream
    m.fs.C102.outlet.pressure.fix(700000)
    m.fs.C102.efficiency_isentropic.fix(0.85)

    TransformationFactory("network.expand_arcs").apply_to(m)


def replace_heater_heat_duty_constraint_with_bounds(m):
    """
    Helper function to check and remove inequality constraints on heater/cooler
    heat_duty and replace with bounds using setlb() and setub()
    """
    if hasattr(m.fs.H101,'heat_duty_constraint'):
        m.fs.H101.del_component(m.fs.H101.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H101)))
    if hasattr(m.fs.H102,'heat_duty_constraint'):
        m.fs.H102.del_component(m.fs.H102.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H102)))
    if hasattr(m.fs.H103,'heat_duty_constraint'):
        m.fs.H103.del_component(m.fs.H103.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H103)))
    if hasattr(m.fs.H104,'heat_duty_constraint'):
        m.fs.H104.del_component(m.fs.H104.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H104)))
    if hasattr(m.fs.H105,'heat_duty_constraint'):
        m.fs.H105.del_component(m.fs.H105.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H105)))
    if hasattr(m.fs.H106,'heat_duty_constraint'):
        m.fs.H106.del_component(m.fs.H106.heat_duty_constraint)
        print('{} heat duty inequality constraint replaced by variable bound.\n'.format(str(m.fs.H106)))
    # H101 is a heater, heat_duty > 0
    m.fs.H101.heat_duty.setlb(0.0)
    # H102 is a cooler, heat_duty < 0
    m.fs.H102.heat_duty.setub(0.0)
    # H103 needs heat, heat_duty > 0
    m.fs.H103.heat_duty.setlb(0.0)
    # H104 is a cooler, heat_duty < 0
    m.fs.H104.heat_duty.setub(0.0)
    # H105 is a cooler, heat_duty < 0
    m.fs.H105.heat_duty.setub(0.0)
    # H106 is a cooler, heat_duty < 0
    m.fs.H106.heat_duty.setub(0.0)


def update_model_for_optimization(m,obj_reset_flag=False):
    """
    Helper function to:
        - Add/update constraints and update DOF bounds for optimization
        - Add objective function for optimization problem
    Arguments:
        obj_reset_flag: boolean, default=False
            True: Delete and re-define objective function to ensure previous solve does not affect current solve
            False: Use existing objective function value
    """
    ### Add/update constraints and update DOF bounds
    print('Model updates:\n\n')
    # Upper bound on recycle ratio
    if hasattr(m.fs, 'recycle_ratio_ub'):
        pass
    else:
        print('Constraint added for recycle ratio <= 8.0.\n')
        m.fs.recycle_ratio_ub = Constraint(expr = m.fs.C102.outlet.flow_mol[0] <= 8.0 * m.fs.M101.feed.flow_mol[0])

    # H103 outlet temperature bound
    print('{} temperature lower and upper bounds added.\n'.format(str(m.fs.H103)))
    m.fs.H103.outlet.temperature.setlb(500.0)
    m.fs.H103.outlet.temperature.setub(700.0)

    # R102 isothermal operation constraint
    def R102_isothermal_constraint(blk,i):
        if i != 0:
            return blk.R102.control_volume.properties[0.0,i].temperature == blk.H103.outlet.temperature[0]
        else:
            return Constraint.Skip
    if hasattr(m.fs,'isothermal_R102'):
        if obj_reset_flag:
            m.fs.del_component(m.fs.isothermal_R102)
            for i in m.fs.R102.control_volume.length_domain:
                m.fs.R102.control_volume.properties[0.0,i].temperature.unfix()
            m.fs.isothermal_R102 = Constraint(m.fs.R102.control_volume.length_domain,rule=R102_isothermal_constraint)
            print('{} isothermal operation constraint re-defined.\n'.format(str(m.fs.R102)))
        else:
            pass
    else:
        print('{} isothermal operation constraint added.\n'.format(str(m.fs.R102)))
        for i in m.fs.R102.control_volume.length_domain:
            m.fs.R102.control_volume.properties[0.0,i].temperature.unfix()
        m.fs.isothermal_R102 = Constraint(m.fs.R102.control_volume.length_domain,rule=R102_isothermal_constraint)

    # H106 outlet temperature bound
    print('{} temperature lower bound updated.\n'.format(str(m.fs.H106)))
    m.fs.H106.outlet.temperature.setlb(290.0)

    # Update F102 outlet temperature lower bound
    if hasattr(m.fs.F102,'outlet_temperature_lb'):
        m.fs.F102.del_component(m.fs.F102.outlet_temperature_lb)
        print('{} temperature lower bound updated.\n'.format(str(m.fs.F102)))
    m.fs.F102.outlet_temperature_lb = Constraint(expr = m.fs.F102.control_volume.properties_out[0.0].temperature >= 288.0)

    ### Add objective function for optimization

    if hasattr(m.fs,'obj'):
        if obj_reset_flag:
            m.fs.del_component(m.fs.obj)
            m.fs.obj = Objective(expr= m.fs.min_sell_price, sense=minimize)
            print('Objective function deleted and re-defined.\n')
        else:
            pass
    else:
        print('Objective function added.\n')
        m.fs.obj = Objective(expr= m.fs.min_sell_price, sense=minimize)


def unfix_DOFs_pre_optimization(m):
    """
    Helper function to unfix optimization degrees of freedom before optimization
    """
    # heat utility
    m.fs.Qs.unfix()

    # H103/R102 temperature
    m.fs.H103.outlet.temperature.unfix()
    
    # H104 pressure
    m.fs.H104.outlet.temperature.unfix()
    
    # H105 temperature
    m.fs.H105.outlet.temperature.unfix()
    
    # F101 pressure drop
    m.fs.F101.deltaP.unfix()

    # H106 temperature and pressure
    m.fs.H106.outlet.temperature.unfix()
    m.fs.H106.outlet.pressure.unfix()
        
    # F102 pressure drop
    m.fs.F102.deltaP.unfix()

def fix_DOFs_post_optimization(m):
    """
    Helper function to fix optimization degrees of freedom post-optimization
    """
    # heat utility
    m.fs.Qs.fix()

    # H103/R102 temperature
    m.fs.H103.outlet.temperature.fix()
    
    # H104 pressure
    m.fs.H104.outlet.temperature.fix()
    
    # H105 temperature
    m.fs.H105.outlet.temperature.fix()
    
    # F101 pressure drop
    m.fs.F101.deltaP.fix()

    # H106 temperature and pressure
    m.fs.H106.outlet.temperature.fix()
    m.fs.H106.outlet.pressure.fix()
        
    # F102 pressure drop
    m.fs.F102.deltaP.fix()


def vapor_only_to_vapor_liquid_reformulate(blk, error_check=True, zero_liquid_outlet=True):
    """
    Function to check if translator block with vapor to vapor-liquid flow 
    has zero flow in vapor or liquid phase. If either phase flow is zero, 
    alter model equations to remove degenerate zero-flow phase composition equations 
    and fix values to 1e-8
    Arguments:
        blk: idaes Translator block model object
             to investigate and alter
        error_check: boolean
            True (default): ensure if all liquid or vapor flowrates are zero
            False: skip check, no need to alter model equations
        zero_liquid_outlet: boolean
            True (default): when translator outlet is supposed to have only liquid flow
            False: If translator outlet has zero vapor flow (provisional for future)
    Return:
        None
    """
    # standard case: zero liquid flow in outlet
    if zero_liquid_outlet:
        if error_check:
            liquid_flow = True # flag to indicate if block outlet has non-zero liquid-phase flow

            # assert liquid flow rate is zero
            assert sum(value(blk.properties_out[0.0].flow_mol_phase_comp[p,j]) for p,j 
                        in blk.properties_out[0.0].phase_component_set if p == 'Liq') < 1e-6,\
                        'Liquid phase flow is not zero. Translator block {} liquid phase composition equations will not be modified.\n'.format(str(blk))
            # check if degenerate constraints have been deactivated previously
            # if active constraints found, set liquid_flow to False
            # if no active constraints found, no further changes are neededm set liquid_flow to True
            if any(blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].active for p,j 
                in blk.properties_out[0.0].phase_component_set if p == 'Liq'):
                print('Liquid phase flow is zero. Translator block {} liquid phase composition equations are being modified...\n'.format(str(blk)))
                liquid_flow = False
            else:
                print('Translator block {} liquid phase composition equations already modified, no degenerate constraints remaining.\n'.format(str(blk)))
            
            # if liquid_flow flag is False, deactivate degenerate liquid phase composition constraints 
            # and fix liquid phase composition values to 1e-8
            if not liquid_flow:
                for p,j in blk.properties_out[0.0].phase_component_set:
                    if p == 'Liq':
                        # check if degenerate constraint is active
                        if blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].active:
                            # deactivate constraint
                            blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].deactivate()
                            # unfix variable
                            blk.properties_out[0.0].mole_frac_phase_comp[p,j].unfix()
                            # fix variable to 1e-8
                            blk.properties_out[0.0].mole_frac_phase_comp[p,j].fix(1e-8)
                print('Degenerate constraints removed and liquid phase composition values set to 1e-8\n')

    # If need is to modify translator block with vapor liquid phase outlet with zero vapor flow
    else:
        if error_check:
            vapor_flow = True # flag to indicate if block outlet has non-zero vapor-phase flow
            # assert vapor flow rate is zero
            assert sum(value(blk.properties_out[0.0].flow_mol_phase_comp[p,j]) for p,j 
                    in blk.properties_out[0.0].phase_component_set if p == 'Vap') < 1e-6, \
                    'Vapor phase flow is not zero. Translator block {} vapor phase composition equations not modified.\n'.format(str(blk))
            # check if degenerate constraints have been deactivated previously
            # if active constraints found, set vapor_flow to False
            # if no active constraints found, no further changes are neededm set vapor_flow to True
            if any(blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].active for p,j 
                in blk.properties_out[0.0].phase_component_set if p == 'Vap'):
                print('Vapor phase flow is zero. Translator block {} vapor phase composition equations are being modified...\n'.format(str(blk)))
                vapor_flow = False
            else:
                print('Translator block {} vapor phase composition equations already modified, no degenerate constraints remaining.\n'.format(str(blk)))
        
            # if vapor_flow flag is False, deactivate degenerate vapor phase composition constraints 
            # and fix vapor phase composition values to 1e-8
            if not vapor_flow:
                for p,j in blk.properties_out[0.0].phase_component_set:
                    if p == 'Vap':
                        # check if degenerate constraint is active
                        if blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].active:
                            # deactivate constraint
                            blk.properties_out[0.0].mole_frac_phase_comp_eq[p,j].deactivate()
                            # unfix variable
                            blk.properties_out[0.0].mole_frac_phase_comp[p,j].unfix()
                            # fix variable to 1e-8
                            blk.properties_out[0.0].mole_frac_phase_comp[p,j].fix(1e-8)
                print('Degenerate constraints removed and vapor phase composition values set to 1e-8\n')


def H106_inlet_vapor_reformulate(blk, error_check=True):
    """
    Function to check if unit block with vapor-liquid flow 
    has zero flow in vapor phase. If vapor phase flow is zero, 
    alter model equations to remove degenerate zero-flow phase composition equations 
    and fix values to 1e-8
    Arguments:
        blk: idaes unit block model object
             block to investigate and alter
        error_check: boolean
            True (default): ensure if all vapor flowrates are zero
            False: skip check, no need to alter model equations
    Return:
        None
    """
    # standard case: zero liquid flow in outlet
    if error_check:
        vapor_flow = True # flag to indicate if block outlet has non-zero vapor-phase flow
        # assert vapor flow rate is zero
        assert sum(value(blk.control_volume.properties_in[0.0].flow_mol_phase_comp[p,j]) for p,j 
                in blk.control_volume.properties_out[0.0].phase_component_set if p == 'Vap') < 1e-6, \
                'Vapor phase flow is not zero. Block {} vapor phase composition equations not modified.\n'.format(str(blk))
        # check if degenerate constraints have been deactivated previously
        # if active constraints found, set vapor_flow to False
        # if no active constraints found, no further changes are neededm set vapor_flow to True
        if any(blk.control_volume.properties_in[0.0].mole_frac_phase_comp_eq[p,j].active for p,j 
            in blk.control_volume.properties_in[0.0].phase_component_set if p == 'Vap'):
            print('Vapor phase flow is zero. Block {} vapor phase composition equations are being modified...\n'.format(str(blk)))
            vapor_flow = False
        else:
            print('Block {} vapor phase composition equations already modified, no degenerate constraints remaining.\n'.format(str(blk)))
    
        # if vapor_flow flag is False, deactivate degenerate vapor phase composition constraints 
        # and fix vapor phase composition values to 1e-8
        if not vapor_flow:
            for p,j in blk.control_volume.properties_in[0.0].phase_component_set:
                if p == 'Vap':
                    # check if degenerate constraint is active
                    if blk.control_volume.properties_in[0.0].mole_frac_phase_comp_eq[p,j].active:
                        # deactivate constraint
                        blk.control_volume.properties_in[0.0].mole_frac_phase_comp_eq[p,j].deactivate()
                        # unfix variable
                        blk.control_volume.properties_in[0.0].mole_frac_phase_comp[p,j].unfix()
                        # fix variable to 1e-8
                        blk.control_volume.properties_in[0.0].mole_frac_phase_comp[p,j].fix(1e-8)
            print('Degenerate constraints removed and vapor phase composition values set to 1e-8\n')

def update_model_after_initialization(m):
    """
    Function to add additional constraints to model and unfix corresponding variables 
    after initialization
    """
    # set compressor C101 outlet pressure to F101 outlet pressure to reduce pressure loss in recycle
    m.fs.C101.outlet.pressure.unfix()
    m.fs.C101_outlet_pressure_constraint = Constraint(expr = m.fs.C101.outlet.pressure[0] == m.fs.F101.vap_outlet.pressure[0])



def create_model_and_return_optimization_results(init_file_name, converged_file_name, costing_file_name,
                                                optimal_file_name,model_code = 5, region = 'Bakken', c_tax_rate=0.0,
                                                M_catalyst=1167.003367):
    """
    Function to create an instance of the process model and return the model with optimization results
    Arguments:
        init_file_name: str
            file name for flowsheet initialization data
        converged_file_name: str
            file name for converged flowsheet data
        costing_file_name: str
            file name for converged flowsheet with costing data
        optimal_file_name: str
            file name for optimal flowsheet data
        model_code: int in 2,3,4,5
            model code for ROK model
            Default: 5
        region: str
            feed composition region
            Default: 'Bakken'
        c_tax_rate: float
            C-tax rate in $/kg CO2e
            Default: 0.0
        M_catalyst: float
            Catalyst loading for oligomerization in kg
    Returns:
        m: Pyomo model
            Flowsheet model with optimal design data
        
    """
    print("\t\t ***** Creating flowsheet model with ROK model M{}, C-tax rate = {} for region {}.... ***** \t\t".format(model_code,c_tax_rate, region))
    # define flowsheet
    m = create_flowsheet(model_code)
    # define unit models
    define_models(m, catalyst_mass = M_catalyst)
    #define connections
    define_arcs(m)

    inlet_df = pd.read_csv('data/NGL_compositions.csv')

    inlet_composition_dict = {}

    for col in inlet_df.columns:
        if col == region:
            for i,r in inlet_df.iterrows():
                if r[col] == 0.0:
                    inlet_composition_dict[r['Species']] = 1e-6
                else:
                    inlet_composition_dict[r['Species']] = round(r[col],4)

    inlet_flow_rate = 481.3888889

    inlet_composition_dict

    dehydro_conv_dict = {'ethane':0.3566,
                        'propane':0.6632,
                        'nbutane':0.5188}
    
    ## Define constraints and set-points for equipment

    set_unit_model_variables(m, model_code=model_code, feed_flow_rate = inlet_flow_rate, 
                            feed_temp = 308.0, feed_pressure = 700000.0,
                            inlet_composition_dict = inlet_composition_dict,
                            dehydro_conv_dict = dehydro_conv_dict)

    ## Scale model components

    if model_code == 2 or model_code == 3:
        set_scaling_factors(m,flow_mol_scaling_factor = 1e-2, inlet_composition_dict = inlet_composition_dict)
    elif model_code == 4 or model_code == 5:
        set_scaling_factors(m,flow_mol_scaling_factor = 1e-3, inlet_composition_dict = inlet_composition_dict)
    else:
        pass

    ## Read-in initialization file
    ms.from_json(m, fname=init_file_name)
    print("Unit initialization data read-in successful.")

    ## Add post-initialization constraints
    update_model_after_initialization(m)
    vapor_only_to_vapor_liquid_reformulate(m.fs.T102)
    vapor_only_to_vapor_liquid_reformulate(m.fs.T102)

    ## Read-in converged flowsheet file
    ms.from_json(m, fname=converged_file_name)
    print("Converged flowsheet data read-in successful.")

    
    ## Equipment costing
    add_costing(m)

    ### Utility minimization
    min_utility(
        m.fs, [m.fs.H101, m.fs.H103, m.fs.R101], [m.fs.H102, m.fs.H104, m.fs.H105, m.fs.H106, m.fs.R102], 10.0
    )
    m.fs.Qs.fix()

    # Emissions calculations
    calc_lhv_values(m,region,'data/LHV.xlsx','data/NGL_compositions.csv','data/NGL_fraction.csv')
    calculate_stream_energies(m)
    calculate_emissions(m,region,'data/emissions_factor_by_region.csv')
    create_ghg_objective(m)
    calculate_costs_for_objective(m,c_tax_flag=True, c_tax_val = c_tax_rate)
    

    # Read-in initialization file
    ms.from_json(m, fname=costing_file_name)
    print('Pre-optimization QW:',m.fs.Qw())
    print('Pre-optimization Qs:',m.fs.Qs())
    print("Costing data read-in successful.")

    CDs_initialized = heat_ex_data(m.fs, [m.fs.H101, m.fs.H103, m.fs.R101], [m.fs.H102, m.fs.H104, m.fs.H105, m.fs.H106, m.fs.R102])
    cold_df_initial, hot_df_initial = generate_and_save_composite_curves(CDs_initialized,model_code, c_tax_rate, region,optimal_solution=False)
    HI_df_initialized = return_HX_results(m.fs,[m.fs.H101, m.fs.H103, m.fs.R101, m.fs.H102, m.fs.H104, m.fs.H105, m.fs.H106, m.fs.R102])


    # Post-costing updates
    replace_heater_heat_duty_constraint_with_bounds(m)
    update_model_for_optimization(m)
    ms.from_json(m, fname=optimal_file_name
                 )
    print('Post-optimization QW:',m.fs.Qw())
    print('Post-optimization Qs:',m.fs.Qs())
    print("Optimization results read-in successful.")

    CDs_optimal = heat_ex_data(m.fs, [m.fs.H101, m.fs.H103, m.fs.R101], [m.fs.H102, m.fs.H104, m.fs.H105, m.fs.H106, m.fs.R102])
    cold_df_optimal, hot_df_optimal = generate_and_save_composite_curves(CDs_optimal,model_code, c_tax_rate, region,optimal_solution=True)
    HI_df_optimal = return_HX_results(m.fs,[m.fs.H101, m.fs.H103, m.fs.R101, m.fs.H102, m.fs.H104, m.fs.H105, m.fs.H106, m.fs.R102])

    print("\t\t ***** Flowsheet model with ROK model M{} and C-tax rate = {} for region {} creation complete. ***** \t\t\n\n".format(model_code, c_tax_rate, region))
    return m, HI_df_initialized, HI_df_optimal, cold_df_optimal, hot_df_optimal


def initialize_flowsheet(m, tear_guesses={}):
    """
    Function to initialize process flowsheet
    Arguments:
        tear_guesses: dict
            Dictionary with guess values for tear stream
    """
    ### Initialization for H101
    m.fs.H101.initialize(outlvl=idaeslog.INFO_HIGH,
                    state_args=tear_guesses)
    propagate_state(m.fs.s02)
    m.fs.H101.report()

    ### Initialization for R101
    m.fs.R101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s03)
    m.fs.R101.report()

    ### Initialization for H102
    m.fs.H102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s04)
    m.fs.H102.report()

    ### Initialization for S101
    m.fs.S101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s05)
    propagate_state(m.fs.s06)
    m.fs.S101.report()

    ### Initialization for T101
    m.fs.T101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s07)
    m.fs.T101.report()

    ### Initialization for E101
    m.fs.E101.initialize(outlvl=idaeslog.INFO_HIGH)
    m.fs.E101.report()

    ### Initialization for H103
    m.fs.H103.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s08)
    m.fs.H103.report()

    ### Initialization for R102
    m.fs.R102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s09)
    m.fs.R102.report()

    ### Initialization for H104
    m.fs.H104.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s10)
    m.fs.H104.report()

    ### Initialization for T102
    m.fs.T102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s11)
    m.fs.T102.report()

    ### Initialization for H105
    m.fs.H105.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s12)
    m.fs.H105.report()

    ### Initialization for F101
    m.fs.F101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s13)
    propagate_state(m.fs.s15)
    m.fs.F101.report()

    ### Initialization for H106
    m.fs.H106.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s14)
    m.fs.H106.report()

    ### Initialization for F102
    m.fs.F102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s16)
    m.fs.F102.report()

    ### Initialization for T103
    m.fs.T103.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s17)
    m.fs.T103.report()

    ### Initialization for T104
    m.fs.T104.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s18)
    m.fs.T104.report()

    ### Initialization for C101
    m.fs.C101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s19)
    m.fs.C101.report()

    ### Initialization for T105
    m.fs.T105.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s20)
    m.fs.T105.report()

    ### Initialization for M102
    m.fs.M102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s21)
    m.fs.M102.report()

    ### Initialization for S102
    m.fs.S102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s22)
    m.fs.S102.report()

    ### Initialization for C102
    m.fs.C102.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s23)
    m.fs.C102.report()

    ### Initialization for M101
    m.fs.M101.initialize(outlvl=idaeslog.INFO_HIGH)
    propagate_state(m.fs.s01)
    m.fs.M101.report()


def set_scaling_factors(m, flow_mol_scaling_factor = 1e-3,inlet_composition_dict={}):
    """
    Function to set scaling for model components
    """
    def scale_variables(m, flow_mol_scaling_factor):
        for var in m.fs.component_data_objects(Var, descend_into=True):
            if 'flow_mol' in var.name:
                iscale.set_scaling_factor(var, flow_mol_scaling_factor)
            if 'temperature' in var.name:
                iscale.set_scaling_factor(var, 1e-2)
                
            if 'pressure' in var.name:
                iscale.set_scaling_factor(var, 1e-6)
            if 'enth_mol' in var.name:
                iscale.set_scaling_factor(var, 1e-6)
            if 'mole_frac' in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if 'mole_frac_comp' in var.name:
                iscale.set_scaling_factor(var, 1e1)
            if 'entr_mol' in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if "_enthalpy_flow_term" in var.name:
                iscale.set_scaling_factor(var, 1e-6)
            if 'rate_reaction_extent' in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if 'heat' in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if 'work' in var.name:
                iscale.set_scaling_factor(var, 1e-2)
            if 'flow_mol_phase_comp' in var.name:
                iscale.set_scaling_factor(var, flow_mol_scaling_factor)
            # if 'mol_frac_phase_comp' in var.name:
            #     iscale.set_scaling_factor(var, flow_mol_scaling_factor)
            if '_enthalpy_flow' in var.name:
                iscale.set_scaling_factor(var, 1e-6)
        
        # manual scaling
        for unit in ('M101','H101', 'R101', 'H102', 'S101', 'T101','E101', 'H103', 'R102', 'H104', 'T102', 'H105',
                    'F101', 'T103', 'H106', 'F102', 'T104', 'C101','T105', 'M102', 'S102', 'C102'):
            block = getattr(m.fs, unit)
    #         if 'mixer' in unit:
            if unit == 'M101':
                iscale.set_scaling_factor(
                    block.feed_state[0.0].enth_mol_phase['Vap'], 1e-6)
                iscale.set_scaling_factor(
                    block.recycle_state[0.0].enth_mol_phase['Vap'], 1e-6)
                iscale.set_scaling_factor(
                    block.mixed_state[0.0].enth_mol_phase['Vap'], 1e-6)

            if unit == 'M102':
                iscale.set_scaling_factor(
                    block.vapor_recycle_F101_state[0.0].enth_mol_phase['Vap'], 1e-6)
                iscale.set_scaling_factor(
                    block.vapor_recycle_F102_state[0.0].enth_mol_phase['Vap'], 1e-6)
                iscale.set_scaling_factor(
                    block.mixed_state[0.0].enth_mol_phase['Vap'], 1e-6)
                    
            if unit in ['F101','F102']:
                for c in inlet_composition_dict.keys():
                    if c!='hydrogen' and c!= 'methane':
                        iscale.set_scaling_factor(
                            block.split._Liq_flow_mol_phase_comp_ref[0.0,'Liq',c], 1e2)
                        iscale.set_scaling_factor(
                            block.split._Liq_flow_mol_phase_comp_ref[0.0,'Vap',c], 1e2)
                        iscale.set_scaling_factor(
                            block.split._Vap_flow_mol_phase_comp_ref[0.0,'Liq',c], 1e2)
                        iscale.set_scaling_factor(
                            block.split._Vap_flow_mol_phase_comp_ref[0.0,'Vap',c], 1e2)
                    else:
                        iscale.set_scaling_factor(
                            block.split._Liq_flow_mol_phase_comp_ref[0.0,'Vap',c], 1e2)
                        
                        iscale.set_scaling_factor(
                            block.split._Vap_flow_mol_phase_comp_ref[0.0,'Vap',c], 1e2)
            
            if 'splitter' in unit:
                iscale.set_scaling_factor(block.mixed_state[0.0].mole_frac_comp,
                                        1e1)
                iscale.set_scaling_factor(block.outlet_1_state[0.0].mole_frac_comp,
                                        1e1)
                iscale.set_scaling_factor(block.outlet_2_state[0.0].mole_frac_comp,
                                        1e1)
                iscale.set_scaling_factor(block.mixed_state[0.0].enth_mol_phase,
                                        1e-6)
                iscale.set_scaling_factor(block.outlet_1_state[0.0].enth_mol_phase,
                                        1e-6)
                iscale.set_scaling_factor(block.outlet_2_state[0.0].enth_mol_phase,
                                        1e-6)
            
            if hasattr(block, "control_volume"):
                if unit != 'R102':
                    iscale.set_scaling_factor(
                        block.control_volume.properties_in[0.0].mole_frac_comp, 1e1)
                    iscale.set_scaling_factor(
                        block.control_volume.properties_out[0.0].mole_frac_comp, 1e1)
                if unit == 'R102':
                    for i in block.control_volume.length_domain:
                        iscale.set_scaling_factor(block.control_volume.properties[0.0,i].enth_mol_phase, 1e-3)
                else:
                    iscale.set_scaling_factor(
                        block.control_volume.properties_in[0.0].enth_mol_phase, 1e-3)
                    iscale.set_scaling_factor(
                        block.control_volume.properties_out[0.0].enth_mol_phase, 1e-3)
                if hasattr(block.control_volume, "rate_reaction_extent"):
                    iscale.set_scaling_factor(block.control_volume
                                            .rate_reaction_extent, 1)
                if hasattr(block.control_volume, "heat"):
                    iscale.set_scaling_factor(block.control_volume.heat, 1e-3)
                if hasattr(block.control_volume, "work"):
                    iscale.set_scaling_factor(block.control_volume.work, 1e-3)

            if hasattr(block, "properties_isentropic"):
                iscale.set_scaling_factor(
                    block.properties_isentropic[0.0].mole_frac_comp, 1e1)
                iscale.set_scaling_factor(
                    block.properties_isentropic[0.0].enth_mol_phase, 1e-3)
            if hasattr(block, "properties"):
                iscale.set_scaling_factor(
                    block.properties[0.0].mole_frac_comp, 1e1)
                iscale.set_scaling_factor(
                    block.properties[0.0].log_mole_frac_phase_comp, 1e-3)
        
        iscale.set_scaling_factor(m.fs.R102.control_volume.area, 1e1)
        
        iscale.calculate_scaling_factors(m)
        
    def scale_constraints(m):
        # set scaling for unit constraints
        for name in ('M101','H101', 'R101', 'H102', 'S101', 'T101','E101', 'H103', 'R102', 'H104', 'T102', 'H105',
                    'F101', 'T103', 'H106', 'F102', 'T104', 'C101', 'T105', 'M102', 'S102', 'C102'):
            unit = getattr(m.fs, name)
            
            # flash constraints
            if hasattr(unit, 'material_splitting_eqn'):
                for (t, o, j), c in unit.material_splitting_eqn.items():
                    iscale.constraint_scaling_transform(c, 1e3, overwrite=False)
            if hasattr(unit, 'temperature_equality_eqn'):
                for (t, o), c in unit.temperature_equality_eqn.items():
                    iscale.constraint_scaling_transform(c, 1e3, overwrite=False)
            if hasattr(unit, 'molar_enthalpy_equality_eqn'):
                for (t, o), c in unit.molar_enthalpy_equality_eqn.items():
                    iscale.constraint_scaling_transform(c, 1e-3, overwrite=False)
            if hasattr(unit, 'molar_enthalpy_splitting_eqn'):
                for (t, o), c in unit.molar_enthalpy_splitting_eqn.items():
                    iscale.constraint_scaling_transform(c, 1e-3, overwrite=False)
            if hasattr(unit, 'pressure_equality_eqn'):
                for (t, o), c in unit.pressure_equality_eqn.items():
                    iscale.constraint_scaling_transform(c, 1e-6, overwrite=False)
            if hasattr(unit, 'sum_split_frac'):
                for t, c in unit.sum_split_frac.items():
                    iscale.constraint_scaling_transform(c, 1e1, overwrite=False)
            if hasattr(unit, 'split_fraction_eq'):
                for (t, o), c in unit.split_fraction_eq.items():
                    iscale.constraint_scaling_transform(c, 1e1, overwrite=False)

            # heater and reactor only add 0D control volume constraints
            if hasattr(unit, 'material_holdup_calculation'):
                for (t, p, j), c in unit.material_holdup_calculation.items():
                    iscale.constraint_scaling_transform(c, 1e3, overwrite=False)
            if hasattr(unit, 'rate_reaction_stoichiometry_constraint'):
                for (t, p, j), c in (
                        unit.rate_reaction_stoichiometry_constraint.items()):
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)
            if hasattr(unit, 'equilibrium_reaction_stoichiometry_constraint'):
                for (t, p, j), c in (
                        unit.equilibrium_reaction_stoichiometry_constraint
                        .items()):
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)
            if hasattr(unit, 'inherent_reaction_stoichiometry_constraint'):
                for (t, p, j), c in (
                        unit.inherent_reaction_stoichiometry_constraint.items()):
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)
            if hasattr(unit, 'material_balances'):
                for (t, p, j), c in unit.material_balances.items():
                    iscale.constraint_scaling_transform(c, 1e3, overwrite=False)
            if hasattr(unit, 'element_balances'):
                for (t, e), c in unit.element_balances.items():
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)
            if hasattr(unit, 'elemental_holdup_calculation'):
                for (t, e), c in unit.elemental_holdup_calculation.items():
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)
            if hasattr(unit, 'enthalpy_balances'):
                for t, c in unit.enthalpy_balances.items():
                    iscale.constraint_scaling_transform(c, 1e-3, overwrite=False)
            if hasattr(unit, 'energy_holdup_calculation'):
                for (t, p), c in unit.energy_holdup_calculation.items():
                    iscale.constraint_scaling_transform(c, 1e-3, overwrite=False)
            if hasattr(unit, 'pressure_balance'):
                for t, c in unit.pressure_balance.items():
                    iscale.constraint_scaling_transform(c, 1e-5, overwrite=False)
            if hasattr(unit, 'sum_of_phase_fractions'):
                for t, c in unit.sum_of_phase_fractions.items():
                    iscale.constraint_scaling_transform(c, 1e1, overwrite=False)
            if hasattr(unit, "material_accumulation_disc_eq"):
                for (t, p, j), c in unit.material_accumulation_disc_eq.items():
                    iscale.constraint_scaling_transform(c, 1e3, overwrite=False)
            if hasattr(unit, "_teq_constraint_Vap_Liq"):
                for t, c in unit._teq_constraint_Vap_Liq.items():
                    iscale.constraint_scaling_transform(c, 1e1, overwrite=False)

            if hasattr(unit, "energy_accumulation_disc_eq"):
                for (t, p), c in unit.energy_accumulation_disc_eq.items():
                    iscale.constraint_scaling_transform(c, 1e-3, overwrite=False)

            if hasattr(unit, "element_accumulation_disc_eq"):
                for (t, e), c in unit.element_accumulation_disc_eq.items():
                    iscale.constraint_scaling_transform(c, 1, overwrite=False)

        # equality constraints between ports at Arc sources and destinations
        for arc in m.fs.component_data_objects(Arc, descend_into=True):
            for c in arc.component_data_objects(Constraint, descend_into=True):
                if hasattr(unit, "enth_mol_equality"):
                    for t, c in unit.enth_mol_equality.items():
                        iscale.constraint_scaling_transform(c, 1e-6,
                                                            overwrite=False)
                if hasattr(unit, "flow_mol_equality"):
                    for t, c in unit.flow_mol_equality.items():
                        iscale.constraint_scaling_transform(c, 1e3,
                                                            overwrite=False)
                if hasattr(unit, "mole_frac_comp_equality"):
                    for (t, j), c in unit.mole_frac_comp_equality.items():
                        iscale.constraint_scaling_transform(c, 1e1,
                                                            overwrite=False)
                if hasattr(unit, "pressure_equality"):
                    for t, c in unit.pressure_equality.items():
                        iscale.constraint_scaling_transform(c, 1e-5,
                                                            overwrite=False)

        iscale.calculate_scaling_factors(m)
    

    scale_variables(m, flow_mol_scaling_factor)
    scale_constraints(m)