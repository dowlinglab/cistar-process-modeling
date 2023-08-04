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
Example of defining a reaction package for dehydrogenation reaction system:

Author: Alejandro Garciadiego
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

from idaes.core import VaporPhase, Component

from idaes.models.properties.modular_properties.base.generic_reaction import (
        ConcentrationForm)
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import \
    arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import \
    power_law_rate
from idaes.models.properties.modular_properties.reactions.equilibrium_constant  import \
    van_t_hoff
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import \
    power_law_equil
from idaes.core import MaterialFlowBasis

# Create the reaction property definition which describes the system on
# reactions to be modeled.
rxn_configuration = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "rate_reactions": {
        "R1": {"stoichiometry": {("Vap", "hydrogen"): 1,
                                ("Vap", "methane"): 0,
                                ("Vap", "ethane"): -1,
                                ("Vap", "propane"): 0,
                                ("Vap", "nbutane"): 0,
                                ("Vap", "ibutane"): 0,
                                ("Vap", "ethylene"): 1,
                                ("Vap", "propene"): 0,
                                ("Vap", "butene"): 0,
                                ("Vap", "pentene"): 0,
                                ("Vap", "hexene"): 0,
                                ("Vap", "heptene"): 0,
                                ("Vap", "octene"): 0},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (1.24e5, pyunits.J/pyunits.mol),
                   "arrhenius_const": 1,
                   "energy_activation": (1000, pyunits.J/pyunits.mol)}},
        "R2": {"stoichiometry": {("Vap", "hydrogen"): 1,
                                ("Vap", "methane"): 0,
                                ("Vap", "ethane"): 0,
                                ("Vap", "propane"): -1,
                                ("Vap", "nbutane"): 0,
                                ("Vap", "ibutane"): 0,
                                ("Vap", "ethylene"): 0,
                                ("Vap", "propene"): 1,
                                ("Vap", "butene"): 0,
                                ("Vap", "pentene"): 0,
                                ("Vap", "hexene"): 0,
                                ("Vap", "heptene"): 0,
                                ("Vap", "octene"): 0},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (1.24e5, pyunits.J/pyunits.mol),
                   "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation": (1000, pyunits.J/pyunits.mol)}},
        "R3": {"stoichiometry": {("Vap", "hydrogen"): 1,
                                ("Vap", "methane"): 0,
                                ("Vap", "ethane"): 0,
                                ("Vap", "propane"): 0,
                                ("Vap", "nbutane"): -1,
                                ("Vap", "ibutane"): 0,
                                ("Vap", "ethylene"): 0,
                                ("Vap", "propene"): 0,
                                ("Vap", "butene"): 1,
                                ("Vap", "pentene"): 0,
                                ("Vap", "hexene"): 0,
                                ("Vap", "heptene"): 0,
                                ("Vap", "octene"): 0},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (1.24e5, pyunits.J/pyunits.mol),
                   "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation": (1000, pyunits.J/pyunits.mol)}},
                   }}
