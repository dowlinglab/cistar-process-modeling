##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Alkene phase equilibrium package using Prng-Robinson liquid and vapor.

Example property package using the Generic Property Package Framework.
"""
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.core.base.phases import PhaseType as PT
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import smooth_VLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import \
        LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity

from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure import RPP4 as RPP

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for a Prng Robinson alkene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
#     Converted to J/mol.K, mol/m^3
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019
# [4] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid

configuration = {
    # Specifying components
    "components": {
        'hydrogen': {"type": Component,
                    "elemental_composition": {'H': 2, 'C': 0},
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "valid_phase_types": PT.vaporPhase,
                    "parameter_data": {
                        "mw": (2.016E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (12.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (33.2, pyunits.K),  # [1]
                        "omega": -0.218,
                        "cp_mol_ig_comp_coeff": {
                            'A': (2.714e1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                            'B': (9.274e-3, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.981e-5, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (7.645e-9, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            0.0, pyunits.J/pyunits.mol)}},

        'methane': {"type": Component,
                    "elemental_composition": {'H': 4, 'C': 1},
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "valid_phase_types": PT.vaporPhase,
                    "parameter_data": {
                        "mw": (16.043E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (46e5, pyunits.Pa),  # [1]
                        "temperature_crit": (190.4, pyunits.K),  # [1]
                        "omega": 0.011,
                        "cp_mol_ig_comp_coeff": {
                            'A': (1.925e1, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (5.213e-2, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (1.197e-5, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-1.132e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -74600, pyunits.J/pyunits.mol)}},

        'ethane': {"type": Component,
                    "elemental_composition": {'H': 6, 'C': 2},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (30.070E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (48.8e5, pyunits.Pa),  # [1]
                        "temperature_crit": (305.4, pyunits.K),  # [1]
                        "omega": 0.099,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.9122, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27937, None),
                            '3': (305.32, pyunits.K),
                            '4': (0.29187, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (5.409e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (1.781e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-6.938e-5, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (8.713e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (90.9, pyunits.J/pyunits.mol/pyunits.K)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            -83800, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.95405, None),  # [4]
                                                    'B': (663.720, pyunits.K),
                                                    'C': (256.681, pyunits.K)}}},
        
        'ethylene': {"type": Component,
                    "elemental_composition": {'H': 4, 'C': 2},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (28.054E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (50.5e5, pyunits.Pa),  # [1]
                        "temperature_crit": (282.4, pyunits.K),  # [1]
                        "omega": 0.089,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (2.0961, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27657, None),
                            '3': (282.34, pyunits.K),
                            '4': (0.29147, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (3.806e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (1.566e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-8.348e-5, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (1.755e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (2.4739e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-4.4280e3, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (4.0936e1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-1.6970e-1, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            52283.264, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            219.5, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.67374, None),  # [4]
                                                    'B': (528.6700, pyunits.K),
                                                    'C': (228.790, pyunits.K)}}}},

    # Specifying phases
    "phases":  {'Vap': {"type": VaporPhase,
                        "equation_of_state": Ideal}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 5000, pyunits.mol/pyunits.s),
                     "temperature": (273.15, 300, 1500, pyunits.K),
                     "pressure": (5e4, 1e5, 1e10, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),


    # Defining phase equilibria
    "parameter_data": {"PR_kappa": {("hydrogen", "hydrogen"): 0.000,
                                ("hydrogen", "methane"): 0.000,
                                ("hydrogen", "ethane"): 0.000,
                                ("hydrogen", "ethylene"): 0.000,
                                ("methane", "hydrogen"): 0.000,
                                ("methane", "methane"): 0.000,
                                ("methane", "ethane"): 0.000,
                                ("methane", "ethylene"): 0.000,
                                ("ethane", "hydrogen"): 0.000,
                                ("ethane", "methane"): 0.000,
                                ("ethane", "ethane"): 0.000,
                                ("ethane", "ethylene"): 0.000,
                                ("ethylene", "hydrogen"): 0.000,
                                ("ethylene", "methane"): 0.000,
                                ("ethylene", "ethane"): 0.000,
                                ("ethylene", "ethylene"): 0.000,
                                }}}
