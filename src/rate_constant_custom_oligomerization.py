#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Methods for calculating rate constants
"""
from pyomo.environ import exp, Var, Param, log, value, units as pyunits

from idaes.core import MaterialFlowBasis
from idaes.models.properties.modular_properties.base.utility import \
    ConcentrationForm
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant dh_rxn
class CoupledOligomerizationCrackingReactionClass():

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        rbasis = parent.config.reaction_basis
        if rbasis == MaterialFlowBasis.molar:
            r_base = units["amount"]
        elif rbasis == MaterialFlowBasis.mass:
            r_base = units["mass"]
        else:
            raise BurntToast(
                "{} for unexpected reaction basis {}. This should not happen "
                "so please contact the IDAES developers with this bug."
                .format(rblock.name, rbasis))

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction or
              c_form == ConcentrationForm.activity):
            r_units = r_base*units["volume"]**-1*units["time"]**-1
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += -rblock.reaction_order[p, j].value

            if c_form == ConcentrationForm.molarity:
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                c_units = units["amount"]*units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units["pressure"]
            else:
                raise BurntToast(
                    "{} received unrecognised ConcentrationForm ({}). "
                    "This should not happen - please contact the IDAES "
                    "developers with this bug."
                    .format(rblock.name, c_form))
            
            r_units = (r_base *
                       units["volume"]**-1 *
                       units["time"]**-1 *
                       c_units**order)
        
        rblock.C_n = Var(
            doc="Carbon number of adsorbed reactant",
            units=None)
        set_param_from_config(rblock, param="C_n",config=config)
            
        rblock.C_m = Var(
            doc="Carbon number of gas-phase reactant reactant",
            units=None)
        set_param_from_config(rblock, param="C_m",config=config)

        rblock.model_code = config.parameter_data["model_code"]
        
        # Retroactively adding packed bed reactor characteristics for catalyst to GenericReactionParameterBlock (parent)
        # so that PFR model can be used to mimic PBR characteristics.
        if hasattr(parent, 'catalyst_mass'):
            pass
        else:
            parent.catalyst_mass = Param(initialize=config.parameter_data['catalyst_mass'], 
                                         mutable=True,units=units["mass"])
        
        if hasattr(parent, 'catalyst_density'):
            pass
        else:
            parent.catalyst_density = Param(initialize=config.parameter_data['catalyst_density'], 
                                            mutable=True,units=units["mass"]*units["volume"]**-1)
        
        if hasattr(parent, 'bulk_density'):
            pass
        else:
            parent.bulk_density = Param(initialize=config.parameter_data['bulk_density'], 
                                        mutable=True,units=units["mass"]*units["volume"]**-1)


        # Retroactively adding ROK model parameters to GenericReactionParameterBlock (parent) so that
        # only 1 copy of the parameters exist in the reactor model. This facilitates the use of
        # uncertainty propagation functions.
        if hasattr(parent,'alpha_olig'):
            pass
        else:
            parent.alpha_olig = Var(initialize=config.parameter_data['alpha_olig'],
                                bounds=(config.parameter_data['alpha_olig'],config.parameter_data['alpha_olig']))
        
        if hasattr(parent,'alpha_crack'):
            pass
        else:
            parent.alpha_crack = Var(initialize=config.parameter_data['alpha_crack'],
                                bounds=(config.parameter_data['alpha_crack'],config.parameter_data['alpha_crack']))
        
        if hasattr(parent,'alpha_ads'):
            pass
        else:
            parent.alpha_ads = Var(initialize=config.parameter_data['alpha_ads'],
                                bounds=(config.parameter_data['alpha_ads'],config.parameter_data['alpha_ads']))
        
        if hasattr(parent,'beta_ads'):
            pass
        else:
            parent.beta_ads = Var(initialize=config.parameter_data['beta_ads'],
                            bounds=(config.parameter_data['beta_ads'],config.parameter_data['beta_ads']))
        if hasattr(parent,'A_phys'):
            pass
        else:
            parent.A_phys = Var(initialize=config.parameter_data['A_phys'],
                                bounds=(config.parameter_data['A_phys'],config.parameter_data['A_phys']))
        
        if hasattr(parent,'A_des'):
            pass
        else:
            parent.A_des = Var(initialize=config.parameter_data['A_des'],
                                bounds=(config.parameter_data['A_des'],config.parameter_data['A_des']))


        if rblock.model_code == 2 or rblock.model_code == 3:
            if hasattr(parent,'beta_olig'):
                pass
            else:
                parent.beta_olig = Var(initialize=config.parameter_data['beta_olig'],
                                bounds=(config.parameter_data['beta_olig'],config.parameter_data['beta_olig']))
            if hasattr(parent,'beta_crack'):
                pass
            else:
                parent.beta_crack = Var(initialize=config.parameter_data['beta_crack'],
                                bounds=(config.parameter_data['beta_crack'],config.parameter_data['beta_crack']))
            if hasattr(parent,'gamma_olig'):
                pass
            else:
                parent.gamma_olig = Var(initialize=config.parameter_data['gamma_olig'],
                                bounds=(config.parameter_data['gamma_olig'],config.parameter_data['gamma_olig']))
            if hasattr(parent,'gamma_crack'):
                pass
            else:
                parent.gamma_crack = Var(initialize=config.parameter_data['gamma_crack'],
                                bounds=(config.parameter_data['gamma_crack'],config.parameter_data['gamma_crack']))
            if hasattr(parent,'phi_crack'):
                pass
            else:
                parent.phi_crack = Var(initialize=config.parameter_data['phi_crack'],
                                bounds=(config.parameter_data['phi_crack'],config.parameter_data['phi_crack']))
            if hasattr(parent,'Ea_olig'):
                pass
            else:
                parent.Ea_olig = Var(initialize=config.parameter_data['Ea_olig'],
                                bounds=(config.parameter_data['Ea_olig'],config.parameter_data['Ea_olig']))
            if hasattr(parent,'Ea_crack'):
                pass
            else:
                parent.Ea_crack = Var(initialize=config.parameter_data['Ea_crack'],
                                bounds=(config.parameter_data['Ea_crack'],config.parameter_data['Ea_crack']))

        elif rblock.model_code == 4 or rblock.model_code == 5:
            if hasattr(parent,'gamma'):
                pass
            else:
                parent.gamma = Var(initialize=config.parameter_data['gamma'],
                                bounds=(config.parameter_data['gamma'],config.parameter_data['gamma']))
            if hasattr(parent,'delta'):
                pass
            else:
                parent.delta = Var(initialize=config.parameter_data['delta'],
                                bounds=(config.parameter_data['delta'],config.parameter_data['delta']))
            if hasattr(parent,'E0'):
                pass
            else:
                parent.E0 = Var(initialize=config.parameter_data['E0'],
                                bounds=(config.parameter_data['E0'],config.parameter_data['E0']))
            if hasattr(parent,'kappa_olig'):
                pass
            else:
                parent.kappa_olig = Var(initialize=config.parameter_data['kappa_olig'],
                                bounds=(config.parameter_data['kappa_olig'],config.parameter_data['kappa_olig']))
            if hasattr(parent,'kappa_crack'):
                pass
            else:
                parent.kappa_crack = Var(initialize=config.parameter_data['kappa_crack'],
                                bounds=(config.parameter_data['kappa_crack'],config.parameter_data['kappa_crack']))
        

        # string to indicate reaction type: 'oligomerization', 'cracking'
        # TODO: add validation for r_type
        rblock.r_type = config.parameter_data["r_type"]
        # set_param_from_config(rblock, param="r_type",config=config)
        
    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        parent = rblock.parent_block()
        units = rblock.parent_block().get_metadata().derived_units
            
        lamda = log(parent.A_phys*101325.0/parent.A_des)
        
        Ea_olig_nm = None
        Ea_crack_nm = None
        delH_formation_n = None
        delH_formation_m = None
        delH_formation_nm = None


        if rblock.model_code == 2 or rblock.model_code == 3:
            Ea_olig_nm = parent.Ea_olig/((exp(parent.gamma_olig*rblock.C_m) - exp(-parent.gamma_olig*rblock.C_m))/(
                        exp(parent.gamma_olig*rblock.C_m)+exp(-parent.gamma_olig*rblock.C_m)))
            Ea_crack_nm = parent.Ea_crack*(1+parent.phi_crack*abs(rblock.C_n-rblock.C_m))/(
                        (exp(parent.gamma_crack*(rblock.C_m+rblock.C_n)) - exp(-parent.gamma_crack*(rblock.C_m+rblock.C_n)))/(
                        exp(parent.gamma_crack*(rblock.C_m+rblock.C_n))+exp(-parent.gamma_crack*(rblock.C_m+rblock.C_n))))
        elif rblock.model_code == 4 or rblock.model_code == 5:
            delH_formation_n = parent.gamma * rblock.C_n + parent.delta
            delH_formation_m = parent.gamma * rblock.C_m + parent.delta
            delH_formation_nm = parent.gamma * (rblock.C_n + rblock.C_m) + parent.delta
        
        # delH_reaction = None
        
        Ea = None
        
        log_k_nm = None
        
        
        if rblock.r_type == 'oligomerization':
            delH_ads = parent.alpha_ads + rblock.C_m * parent.beta_ads
            if rblock.model_code == 4 or rblock.model_code == 5:
                delH_reaction = delH_formation_nm - delH_formation_m - delH_formation_n
                if value(delH_reaction) <= 0.0:
                    Ea = parent.E0 + parent.kappa_olig * delH_reaction
                else:
                    Ea = parent.E0 + (1.0 - parent.kappa_olig) * delH_reaction
            
            if rblock.model_code == 2:
                # alpha multiplier: used to convert units from 1/atm/min/g to 1/Pa/s/kg (1000 (g to kg)/60 (s to min) * 101325 (atm to Pa))
                log_k_nm = parent.alpha_olig + log(1000.0/60.0/101325.0) + lamda - log(101325.0) + Ea_olig_nm * parent.beta_olig  + (delH_ads-Ea_olig_nm)/(
                                                                                                                                    pyunits.convert(c.gas_constant,
                                                                                                                                    to_units=units["gas_constant"])*T)
            elif rblock.model_code == 3:
                # alpha multiplier: used to convert units from 1/atm/min/g to 1/Pa/s/kg (1000 (g to kg)/60 (s to min) * 101325 (atm to Pa))
                log_k_nm = parent.alpha_olig + log(1000.0/60.0/101325.0) + log(T/298.15/units["temperature"]) + lamda - log(101325.0) + Ea_olig_nm * parent.beta_olig  + (delH_ads-Ea_olig_nm)/(
                                                                                                                                                                        pyunits.convert(c.gas_constant,
                                                                                                                                                                        to_units=units["gas_constant"])*T)
            elif rblock.model_code == 4:
                # alpha multiplier: used to convert units from 1/atm/min/g to 1/Pa/s/kg (1000 (g to kg)/60 (s to min) * 101325 (atm to Pa))
                log_k_nm = parent.alpha_olig + log(1000.0/60.0/101325.0) + lamda - log(101325.0) + (delH_ads-Ea)/(
                                                                                    pyunits.convert(c.gas_constant,
                                                                                    to_units=units["gas_constant"])*T)
            elif rblock.model_code == 5:
                # alpha multiplier: used to convert units from 1/atm/min/g to 1/Pa/s/kg (1000 (g to kg)/60 (s to min) * 101325 (atm to Pa))
                log_k_nm = parent.alpha_olig + log(1000.0/60.0/101325.0) + log(T/298.15/units["temperature"]) + lamda - log(101325.0) + (delH_ads-Ea)/(
                                                                                                                                        pyunits.convert(c.gas_constant,
                                                                                                                                        to_units=units["gas_constant"])*T)
            
        elif rblock.r_type == 'cracking':
            delH_ads = parent.alpha_ads + (rblock.C_n+rblock.C_m) * parent.beta_ads
            if rblock.model_code == 4 or rblock.model_code == 5:
                delH_reaction = delH_formation_m + delH_formation_n - delH_formation_nm
                if value(delH_reaction) <= 0.0:
                    Ea = parent.E0 + parent.kappa_crack * delH_reaction
                else:
                    Ea = parent.E0 + (1.0 - parent.kappa_crack) * delH_reaction
            
            if rblock.model_code == 2:
                # alpha multiplier: used to convert units from 1/min/g to 1/s/kg (1000 (g to kg)/60 (s to min))
                log_k_nm = parent.alpha_crack + log(1000.0/60.0) + lamda - log(101325.0) + Ea_crack_nm * parent.beta_crack + (delH_ads-Ea_crack_nm)/(
                                                                                                                            pyunits.convert(c.gas_constant,
                                                                                                                            to_units=units["gas_constant"])*T)
            elif rblock.model_code == 3:
                # alpha multiplier: used to convert units from 1/min/g to 1/s/kg (1000 (g to kg)/60 (s to min))
                log_k_nm = parent.alpha_crack + log(1000.0/60.0) + log(T/298.15/units["temperature"]) + lamda - log(101325.0) + Ea_crack_nm * parent.beta_crack + (delH_ads-Ea_crack_nm)/(
                                                                                                                                                                pyunits.convert(c.gas_constant,
                                                                                                                                                                to_units=units["gas_constant"])*T)
            elif rblock.model_code == 4:
                # alpha multiplier: used to convert units from 1/min/g to 1/s/kg (1000 (g to kg)/60 (s to min))
                log_k_nm = parent.alpha_crack + log(1000.0/60.0) + lamda - log(101325.0) + (delH_ads-Ea)/(
                                                                                            pyunits.convert(c.gas_constant,
                                                                                            to_units=units["gas_constant"])*T)
            elif rblock.model_code == 5:
                # alpha multiplier: used to convert units from 1/min/g to 1/s/kg (1000 (g to kg)/60 (s to min))
                log_k_nm = parent.alpha_crack + log(1000.0/60.0) + log(T/298.15/units["temperature"]) + lamda - log(101325.0) + (delH_ads-Ea)/(
                                                                                                                                pyunits.convert(c.gas_constant,
                                                                                                                                to_units=units["gas_constant"])*T)
        
        return parent.bulk_density*exp(log_k_nm)
