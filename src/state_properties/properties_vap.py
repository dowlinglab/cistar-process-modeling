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
# Configuration dictionary for a Peng Robinson alkene system

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
        'propane': {"type": Component,
                    "elemental_composition": {'H': 8, 'C': 3},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (44.094E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (42.5e5, pyunits.Pa),  # [1]
                        "temperature_crit": (369.8, pyunits.K),  # [1]
                        "omega": 0.144,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.3757, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.274253, None),
                            '3': (369.83, pyunits.K),
                            '4': (0.29147, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-4.224e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.063e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.586e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (3.215e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (112.71, pyunits.J/pyunits.mol/pyunits.K)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_liq_comp_ref": (
                            0.0, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            -104700, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.92828, None),  # [4]
                                                    'B': (803.9970, pyunits.K),
                                                    'C': (247.040, pyunits.K)}}},
        'nbutane': {"type": Component,
                    "elemental_composition": {'H': 10, 'C': 4},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (58.124E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (38.0e5, pyunits.Pa),  # [1]
                        "temperature_crit": (425.2, pyunits.K),  # [1]
                        "omega": 0.199,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.0677, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27188, None),
                            '3': (425.12, pyunits.K),
                            '4': (0.28688, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (9.487e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.313e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.108e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (-2.822e-9, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (132.42, pyunits.J/pyunits.mol/pyunits.K)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            -125600, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.93266, None),  # [4]
                                                    'B': (935.7730, pyunits.K),
                                                    'C': (238.789, pyunits.K)}}},
        'ibutane': {"type": Component,
                    "elemental_composition": {'H': 10, 'C': 4},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (58.124E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (36.5e5, pyunits.Pa),  # [1]
                        "temperature_crit": (408.2, pyunits.K),  # [1]
                        "omega": 0.183,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.0463, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27294, None),
                            '3': (408.14, pyunits.K),
                            '4': (0.27301, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.390e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.847e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.846e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.895e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.7237e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-1.7839e3, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (1.4759e1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-4.7909e-2, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (5.8050e-5, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -135600, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            310.1, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.00272, None),  # [4]
                                                    'B': (947.5400, pyunits.K),
                                                    'C': (248.870, pyunits.K)}}},
        'pentane': {"type": Component,
                    "elemental_composition": {'H': 12, 'C': 5},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (72.149E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (33.6e5, pyunits.Pa),  # [1]
                        "temperature_crit": (469.8, pyunits.K),  # [1]
                        "omega": 0.252,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.84947, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26726, None),
                            '3': (469.7, pyunits.K),
                            '4': (0.27789, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.390e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.847e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.846e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.895e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.5908e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-2.7050e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (9.9537e-1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -146900, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            349.0, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.01432, None),  # [4]
                                                    'B': (1075.78, pyunits.K),
                                                    'C': (232.546, pyunits.K)}}},
        'hexane': {"type": Component,
                    "elemental_composition": {'H': 14, 'C': 6},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (86.178E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (30.31e5, pyunits.Pa),  # [1]
                        "temperature_crit": (507.44, pyunits.K),  # [1]
                        "omega": 0.300,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.70824, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26411, None),
                            '3': (507.6, pyunits.K),
                            '4': (0.27537, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.390e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.847e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.846e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.895e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.7212e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-1.8378e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (8.8734e-1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -166900, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            388.0, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.02105, None),  # [4]
                                                    'B': (1168.72, pyunits.K),
                                                    'C': (224.216, pyunits.K)}}},
        'heptane': {"type": Component,
                    "elemental_composition": {'H': 16, 'C': 7},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (100.198E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (27.36e5, pyunits.Pa),  # [1]
                        "temperature_crit": (540.61, pyunits.K),  # [1]
                        "omega": 0.350,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.61259, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26211, None),
                            '3': (540.2, pyunits.K),
                            '4': (0.28141, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.390e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.847e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.846e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.895e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (224.7, pyunits.J/pyunits.mol/pyunits.K)},
                        "entr_mol_form_vap_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            -187600, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.03755, None),  # [4]
                                                    'B': (1264.37, pyunits.K),
                                                    'C': (216.636, pyunits.K)}}},
        'octane': {"type": Component,
                    "elemental_composition": {'H': 18, 'C': 8},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (114.224E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (24.86e5, pyunits.Pa),  # [1]
                        "temperature_crit": (568.8, pyunits.K),  # [1]
                        "omega": 0.399,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.53731, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26115, None),
                            '3': (568.7, pyunits.K),
                            '4': (0.28034, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.390e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.847e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.846e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.895e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (2.2483e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-1.8663e3, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (9.5891e-1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0.0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -208500, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            467.0, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.04132, None),  # [4]
                                                    'B': (1349.82, pyunits.K),
                                                    'C': (209.385, pyunits.K)}}},
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
                                                    'C': (228.790, pyunits.K)}}},
        'propene': {"type": Component,
                    "elemental_composition": {'H': 6, 'C': 3},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (42.081E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (46.2e5, pyunits.Pa),  # [1]
                        "temperature_crit": (365.0, pyunits.K),  # [1]
                        "omega": 0.144,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.4094, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26465, None),
                            '3': (365.57, pyunits.K),
                            '4': (0.295, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (3.710e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (2.345e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.160e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (2.205e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (7.1720e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-3.8632e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (1.2348e0, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            20413.736, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            266.9, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.95606, None),  # [4]
                                                    'B': (789.6200, pyunits.K),
                                                    'C': (247.580, pyunits.K)}}},
        'butene': {"type": Component,
                    "elemental_composition": {'H': 8, 'C': 4},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (56.104E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (40.2e5, pyunits.Pa),  # [1]
                        "temperature_crit": (419.3, pyunits.K),  # [1]
                        "omega": 0.191,
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (1.0972, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.2649, None),
                            '3': (419.95, pyunits.K),
                            '4': (0.29043, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-2.994e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (3.532e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-1.990e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (4.463e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.3589e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-4.7739e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (2.1835e0, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-2.2230e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            1171.52, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            305.6, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.96640, None),  # [4]
                                                    'B': (927.2100, pyunits.K),
                                                    'C': (238.630, pyunits.K)}}},
        'pentene': {"type": Component,
                    "elemental_composition": {'H': 10, 'C': 5},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (70.135E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (40.5e5, pyunits.Pa),  # [1]
                        "temperature_crit": (464.7, pyunits.K),  # [1]
                        "omega": 0.233,  # [1]
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.9038, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26648, None),
                            '3': (464.78, pyunits.K),
                            '4': (0.2905, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.340e-1, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (4.329e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-2.317e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (4.681e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.5467e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-4.2600e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (1.964, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-1.8034e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -20920.00, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            347.82, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.96914, None),  # [4]
                                                    'B': (1044.010, pyunits.K),
                                                    'C': (233.450, pyunits.K)}}},
        'hexene': {"type": Component,
                    "elemental_composition": {'H': 12, 'C': 6},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (84.162E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (31.7e5, pyunits.Pa),  # [1]
                        "temperature_crit": (504.0, pyunits.K),  # [1]
                        "omega": 0.285,  # [1]
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.7389, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26147, None),
                            '3': (504.03, pyunits.K),
                            '4': (0.2902, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-1.746e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (5.309e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-2.903e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (6.054e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.9263e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-5.7116e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (2.4004e0, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-1.9758e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -41672.64, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            385.97, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (3.98260, None),  # [4]
                                                    'B': (1148.620, pyunits.K),
                                                    'C': (225.340, pyunits.K)}}},
        'heptene': {"type": Component,
                    "elemental_composition": {'H': 14, 'C': 7},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (98.189E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (25.4e5, pyunits.Pa),  # [1]
                        "temperature_crit": (537.2, pyunits.K),  # [1]
                        "omega": 0.358,  # [1]
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.63734, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26319, None),
                            '3': (537.29, pyunits.K),
                            '4': (0.27375, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-3.303e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (6.297e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-3.512e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (7.607E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.8997e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-1.5670e2, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (3.4300e-1, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (1.5222e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -62299.76, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            424, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.02677, None),  # [4]
                                                    'B': (1258.340, pyunits.K),
                                                    'C': (219.300, pyunits.K)}}},
        'octene': {"type": Component,
                    "elemental_composition": {'H': 16, 'C': 8},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (112.216E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (26.2e5, pyunits.Pa),  # [1]
                        "temperature_crit": (566.6, pyunits.K),  # [1]
                        "omega": 0.386,  # [1]
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.5871, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27005, None),
                            '3': (566.65, pyunits.K),
                            '4': (0.27187, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-4.099e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (7.239e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-4.036e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (8.675e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (3.7930e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-2.1175e3, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (8.2362e0, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-9.0093e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -82926.88, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            454.36, pyunits.J/pyunits.mol/pyunits.K),  # [3]
                        "pressure_sat_comp_coeff": {'A': (4.05985, None),  # [4]
                                                    'B': (1355.460, pyunits.K),
                                                    'C': (213.050, pyunits.K)}}},

        'nonene': {"type": Component,
                    "elemental_composition": {'H': 18, 'C': 9},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "entr_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (112.216E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (26.2e5, pyunits.Pa),  # [1]
                        "temperature_crit": (566.6, pyunits.K),  # [1]
                        "omega": 0.386,  # [1]
                        "dens_mol_liq_comp_coeff": {"eqn_type": 1,
                            '1': (0.5871, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.27005, None),
                            '3': (566.65, pyunits.K),
                            '4': (0.27187, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-4.099e0, pyunits.J/pyunits.mol/pyunits.K),
                            'B': (7.239e-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-4.036e-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (8.675e-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (3.7930e5, pyunits.J*pyunits.kmol**-1*pyunits.K**-1),  # [2]
                            '2': (-2.1175e3, pyunits.J*pyunits.kmol**-1*pyunits.K**-2),
                            '3': (8.2362e0, pyunits.J*pyunits.kmol**-1*pyunits.K**-3),
                            '4': (-9.0093e-3, pyunits.J*pyunits.kmol**-1*pyunits.K**-4),
                            '5': (0, pyunits.J*pyunits.kmol**-1*pyunits.K**-5)},
                        "enth_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_liq_comp_ref": (
                            0, pyunits.J/pyunits.mol/pyunits.K),
                        "enth_mol_form_vap_comp_ref": (
                            -82926.88, pyunits.J/pyunits.mol),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            454.36, pyunits.J/pyunits.mol/pyunits.K)}}},  # [3]
                        # "pressure_sat_comp_coeff": {'A': (4.05985, None),  # [4]
                        #                             'B': (1355.460, pyunits.K),
                        #                             'C': (213.050, pyunits.K)}}}},



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
                                ("hydrogen", "propane"): 0.000,
                                ("hydrogen", "nbutane"): 0.000,
                                ("hydrogen", "ibutane"): 0.000,
                                ("hydrogen", "pentane"): 0.000,
                                ("hydrogen", "hexane"): 0.000,
                                ("hydrogen", "heptane"): 0.000,
                                ("hydrogen", "octane"): 0.000,
                                ("hydrogen", "ethylene"): 0.000,
                                ("hydrogen", "propene"): 0.000,
                                ("hydrogen", "butene"): 0.000,
                                ("hydrogen", "pentene"): 0.000,
                                ("hydrogen", "hexene"): 0.000,
                                ("hydrogen", "heptene"): 0.000,
                                ("hydrogen", "octene"): 0.000,
                                ("hydrogen", "nonene"): 0.000,
                                ("methane", "hydrogen"): 0.000,
                                ("methane", "methane"): 0.000,
                                ("methane", "ethane"): 0.000,
                                ("methane", "propane"): 0.000,
                                ("methane", "nbutane"): 0.000,
                                ("methane", "ibutane"): 0.000,
                                ("methane", "pentane"): 0.000,
                                ("methane", "hexane"): 0.000,
                                ("methane", "heptane"): 0.000,
                                ("methane", "octane"): 0.000,
                                ("methane", "ethylene"): 0.000,
                                ("methane", "propene"): 0.000,
                                ("methane", "butene"): 0.000,
                                ("methane", "pentene"): 0.000,
                                ("methane", "hexene"): 0.000,
                                ("methane", "heptene"): 0.000,
                                ("methane", "octene"): 0.000,
                                ("methane", "nonene"): 0.000,
                                ("ethane", "hydrogen"): 0.000,
                                ("ethane", "methane"): 0.000,
                                ("ethane", "ethane"): 0.000,
                                ("ethane", "propane"): 0.000,
                                ("ethane", "nbutane"): 0.000,
                                ("ethane", "ibutane"): 0.000,
                                ("ethane", "pentane"): 0.000,
                                ("ethane", "hexane"): 0.000,
                                ("ethane", "heptane"): 0.000,
                                ("ethane", "octane"): 0.000,
                                ("ethane", "ethylene"): 0.000,
                                ("ethane", "propene"): 0.000,
                                ("ethane", "butene"): 0.000,
                                ("ethane", "pentene"): 0.000,
                                ("ethane", "hexene"): 0.000,
                                ("ethane", "heptene"): 0.000,
                                ("ethane", "octene"): 0.000,
                                ("ethane", "nonene"): 0.000,
                                ("propane", "hydrogen"): 0.000,
                                ("propane", "methane"): 0.000,
                                ("propane", "ethane"): 0.000,
                                ("propane", "propane"): 0.000,
                                ("propane", "nbutane"): 0.000,
                                ("propane", "ibutane"): 0.000,
                                ("propane", "pentane"): 0.000,
                                ("propane", "hexane"): 0.000,
                                ("propane", "heptane"): 0.000,
                                ("propane", "octane"): 0.000,
                                ("propane", "ethylene"): 0.000,
                                ("propane", "propene"): 0.000,
                                ("propane", "butene"): 0.000,
                                ("propane", "pentene"): 0.000,
                                ("propane", "hexene"): 0.000,
                                ("propane", "heptene"): 0.000,
                                ("propane", "octene"): 0.000,
                                ("propane", "nonene"): 0.000,
                                ("nbutane", "hydrogen"): 0.000,
                                ("nbutane", "methane"): 0.000,
                                ("nbutane", "ethane"): 0.000,
                                ("nbutane", "propane"): 0.000,
                                ("nbutane", "nbutane"): 0.000,
                                ("nbutane", "ibutane"): 0.000,
                                ("nbutane", "pentane"): 0.000,
                                ("nbutane", "hexane"): 0.000,
                                ("nbutane", "heptane"): 0.000,
                                ("nbutane", "octane"): 0.000,
                                ("nbutane", "ethylene"): 0.000,
                                ("nbutane", "propene"): 0.000,
                                ("nbutane", "butene"): 0.000,
                                ("nbutane", "pentene"): 0.000,
                                ("nbutane", "hexene"): 0.000,
                                ("nbutane", "heptene"): 0.000,
                                ("nbutane", "octene"): 0.000,
                                ("nbutane", "nonene"): 0.000,
                                ("ibutane", "hydrogen"): 0.000,
                                ("ibutane", "methane"): 0.000,
                                ("ibutane", "ethane"): 0.000,
                                ("ibutane", "propane"): 0.000,
                                ("ibutane", "nbutane"): 0.000,
                                ("ibutane", "ibutane"): 0.000,
                                ("ibutane", "pentane"): 0.000,
                                ("ibutane", "hexane"): 0.000,
                                ("ibutane", "heptane"): 0.000,
                                ("ibutane", "octane"): 0.000,
                                ("ibutane", "ethylene"): 0.000,
                                ("ibutane", "propene"): 0.000,
                                ("ibutane", "butene"): 0.000,
                                ("ibutane", "pentene"): 0.000,
                                ("ibutane", "hexene"): 0.000,
                                ("ibutane", "heptene"): 0.000,
                                ("ibutane", "octene"): 0.000,
                                ("ibutane", "nonene"): 0.000,
                                ("pentane", "hydrogen"): 0.000,
                                ("pentane", "methane"): 0.000,
                                ("pentane", "ethane"): 0.000,
                                ("pentane", "propane"): 0.000,
                                ("pentane", "nbutane"): 0.000,
                                ("pentane", "ibutane"): 0.000,
                                ("pentane", "pentane"): 0.000,
                                ("pentane", "hexane"): 0.000,
                                ("pentane", "heptane"): 0.000,
                                ("pentane", "octane"): 0.000,
                                ("pentane", "ethylene"): 0.000,
                                ("pentane", "propene"): 0.000,
                                ("pentane", "butene"): 0.000,
                                ("pentane", "pentene"): 0.000,
                                ("pentane", "hexene"): 0.000,
                                ("pentane", "heptene"): 0.000,
                                ("pentane", "octene"): 0.000,
                                ("pentane", "nonene"): 0.000,
                                ("hexane", "hydrogen"): 0.000,
                                ("hexane", "methane"): 0.000,
                                ("hexane", "ethane"): 0.000,
                                ("hexane", "propane"): 0.000,
                                ("hexane", "nbutane"): 0.000,
                                ("hexane", "ibutane"): 0.000,
                                ("hexane", "pentane"): 0.000,
                                ("hexane", "hexane"): 0.000,
                                ("hexane", "heptane"): 0.000,
                                ("hexane", "octane"): 0.000,
                                ("hexane", "ethylene"): 0.000,
                                ("hexane", "propene"): 0.000,
                                ("hexane", "butene"): 0.000,
                                ("hexane", "pentene"): 0.000,
                                ("hexane", "hexene"): 0.000,
                                ("hexane", "heptene"): 0.000,
                                ("hexane", "octene"): 0.000,
                                ("hexane", "nonene"): 0.000,
                                ("heptane", "hydrogen"): 0.000,
                                ("heptane", "methane"): 0.000,
                                ("heptane", "ethane"): 0.000,
                                ("heptane", "propane"): 0.000,
                                ("heptane", "nbutane"): 0.000,
                                ("heptane", "ibutane"): 0.000,
                                ("heptane", "pentane"): 0.000,
                                ("heptane", "hexane"): 0.000,
                                ("heptane", "heptane"): 0.000,
                                ("heptane", "octane"): 0.000,
                                ("heptane", "ethylene"): 0.000,
                                ("heptane", "propene"): 0.000,
                                ("heptane", "butene"): 0.000,
                                ("heptane", "pentene"): 0.000,
                                ("heptane", "hexene"): 0.000,
                                ("heptane", "heptene"): 0.000,
                                ("heptane", "octene"): 0.000,
                                ("heptane", "nonene"): 0.000,
                                ("octane", "hydrogen"): 0.000,
                                ("octane", "methane"): 0.000,
                                ("octane", "ethane"): 0.000,
                                ("octane", "propane"): 0.000,
                                ("octane", "nbutane"): 0.000,
                                ("octane", "ibutane"): 0.000,
                                ("octane", "pentane"): 0.000,
                                ("octane", "hexane"): 0.000,
                                ("octane", "heptane"): 0.000,
                                ("octane", "octane"): 0.000,
                                ("octane", "ethylene"): 0.000,
                                ("octane", "propene"): 0.000,
                                ("octane", "butene"): 0.000,
                                ("octane", "pentene"): 0.000,
                                ("octane", "hexene"): 0.000,
                                ("octane", "heptene"): 0.000,
                                ("octane", "octene"): 0.000,
                                ("octane", "nonene"): 0.000,
                                ("ethylene", "hydrogen"): 0.000,
                                ("ethylene", "methane"): 0.000,
                                ("ethylene", "ethane"): 0.000,
                                ("ethylene", "propane"): 0.000,
                                ("ethylene", "nbutane"): 0.000,
                                ("ethylene", "ibutane"): 0.000,
                                ("ethylene", "pentane"): 0.000,
                                ("ethylene", "hexane"): 0.000,
                                ("ethylene", "heptane"): 0.000,
                                ("ethylene", "octane"): 0.000,
                                ("ethylene", "ethylene"): 0.000,
                                ("ethylene", "propene"): 0.000,
                                ("ethylene", "butene"): 0.000,
                                ("ethylene", "pentene"): 0.000,
                                ("ethylene", "hexene"): 0.000,
                                ("ethylene", "heptene"): 0.000,
                                ("ethylene", "octene"): 0.000,
                                ("ethylene", "nonene"): 0.000,
                                ("propene", "hydrogen"): 0.000,
                                ("propene", "methane"): 0.000,
                                ("propene", "ethane"): 0.000,
                                ("propene", "propane"): 0.000,
                                ("propene", "nbutane"): 0.000,
                                ("propene", "ibutane"): 0.000,
                                ("propene", "pentane"): 0.000,
                                ("propene", "hexane"): 0.000,
                                ("propene", "heptane"): 0.000,
                                ("propene", "octane"): 0.000,
                                ("propene", "ethylene"): 0.000,
                                ("propene", "propene"): 0.000,
                                ("propene", "butene"): 0.000,
                                ("propene", "pentene"): 0.000,
                                ("propene", "hexene"): 0.000,
                                ("propene", "heptene"): 0.000,
                                ("propene", "octene"): 0.000,
                                ("propene", "nonene"): 0.000,
                                ("butene", "hydrogen"): 0.000,
                                ("butene", "methane"): 0.000,
                                ("butene", "ethane"): 0.000,
                                ("butene", "propane"): 0.000,
                                ("butene", "nbutane"): 0.000,
                                ("butene", "ibutane"): 0.000,
                                ("butene", "pentane"): 0.000,
                                ("butene", "hexane"): 0.000,
                                ("butene", "heptane"): 0.000,
                                ("butene", "octane"): 0.000,
                                ("butene", "ethylene"): 0.000,
                                ("butene", "propene"): 0.000,
                                ("butene", "butene"): 0.000,
                                ("butene", "pentene"): 0.000,
                                ("butene", "hexene"): 0.000,
                                ("butene", "heptene"): 0.000,
                                ("butene", "octene"): 0.000,
                                ("butene", "nonene"): 0.000,
                                ("pentene", "hydrogen"): 0.000,
                                ("pentene", "methane"): 0.000,
                                ("pentene", "ethane"): 0.000,
                                ("pentene", "propane"): 0.000,
                                ("pentene", "nbutane"): 0.000,
                                ("pentene", "ibutane"): 0.000,
                                ("pentene", "pentane"): 0.000,
                                ("pentene", "hexane"): 0.000,
                                ("pentene", "heptane"): 0.000,
                                ("pentene", "octane"): 0.000,
                                ("pentene", "ethylene"): 0.000,
                                ("pentene", "propene"): 0.000,
                                ("pentene", "butene"): 0.000,
                                ("pentene", "pentene"): 0.000,
                                ("pentene", "hexene"): 0.000,
                                ("pentene", "heptene"): 0.000,
                                ("pentene", "octene"): 0.000,
                                ("pentene", "nonene"): 0.000,
                                ("hexene", "hydrogen"): 0.000,
                                ("hexene", "methane"): 0.000,
                                ("hexene", "ethane"): 0.000,
                                ("hexene", "propane"): 0.000,
                                ("hexene", "nbutane"): 0.000,
                                ("hexene", "ibutane"): 0.000,
                                ("hexene", "pentane"): 0.000,
                                ("hexene", "hexane"): 0.000,
                                ("hexene", "heptane"): 0.000,
                                ("hexene", "octane"): 0.000,
                                ("hexene", "ethylene"): 0.000,
                                ("hexene", "propene"): 0.000,
                                ("hexene", "butene"): 0.000,
                                ("hexene", "pentene"): 0.000,
                                ("hexene", "hexene"): 0.000,
                                ("hexene", "heptene"): 0.000,
                                ("hexene", "octene"): 0.000,
                                ("hexene", "nonene"): 0.000,
                                ("heptene", "hydrogen"): 0.000,
                                ("heptene", "methane"): 0.000,
                                ("heptene", "ethane"): 0.000,
                                ("heptene", "propane"): 0.000,
                                ("heptene", "nbutane"): 0.000,
                                ("heptene", "ibutane"): 0.000,
                                ("heptene", "pentane"): 0.000,
                                ("heptene", "hexane"): 0.000,
                                ("heptene", "heptane"): 0.000,
                                ("heptene", "octane"): 0.000,
                                ("heptene", "ethylene"): 0.000,
                                ("heptene", "propene"): 0.000,
                                ("heptene", "butene"): 0.000,
                                ("heptene", "pentene"): 0.000,
                                ("heptene", "hexene"): 0.000,
                                ("heptene", "heptene"): 0.000,
                                ("heptene", "octene"): 0.000,
                                ("heptene", "nonene"): 0.000,
                                ("octene", "hydrogen"): 0.000,
                                ("octene", "methane"): 0.000,
                                ("octene", "ethane"): 0.000,
                                ("octene", "propane"): 0.000,
                                ("octene", "nbutane"): 0.000,
                                ("octene", "ibutane"): 0.000,
                                ("octene", "pentane"): 0.000,
                                ("octene", "hexane"): 0.000,
                                ("octene", "heptane"): 0.000,
                                ("octene", "octane"): 0.000,
                                ("octene", "ethylene"): 0.000,
                                ("octene", "propene"): 0.000,
                                ("octene", "butene"): 0.000,
                                ("octene", "pentene"): 0.000,
                                ("octene", "hexene"): 0.000,
                                ("octene", "heptene"): 0.000,
                                ("octene", "octene"): 0.000,
                                ("octene", "nonene"): 0.000,
                                ("nonene", "hydrogen"): 0.000,
                                ("nonene", "methane"): 0.000,
                                ("nonene", "ethane"): 0.000,
                                ("nonene", "propane"): 0.000,
                                ("nonene", "nbutane"): 0.000,
                                ("nonene", "ibutane"): 0.000,
                                ("nonene", "pentane"): 0.000,
                                ("nonene", "hexane"): 0.000,
                                ("nonene", "heptane"): 0.000,
                                ("nonene", "octane"): 0.000,
                                ("nonene", "ethylene"): 0.000,
                                ("nonene", "propene"): 0.000,
                                ("nonene", "butene"): 0.000,
                                ("nonene", "pentene"): 0.000,
                                ("nonene", "hexene"): 0.000,
                                ("nonene", "heptene"): 0.000,
                                ("nonene", "octene"): 0.000,
                                ("nonene", "nonene"): 0.000}}}