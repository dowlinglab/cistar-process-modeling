import pandas as pd
from idaes.models.properties.modular_properties.reactions.dh_rxn import \
    constant_dh_rxn

from idaes.models.properties.modular_properties.base.utility import (
        ConcentrationForm)

from src.rate_constant_custom_oligomerization import \
    CoupledOligomerizationCrackingReactionClass

from src.rate_forms_oligo import \
    power_law_rate

def return_reaction_network(model_code = 2,
                            catalyst_mass = 1167.003367, 
                            catalyst_density=1.2*10**3,
                            bulk_density=7.2*10**2
                            ):
    """
    Function to generate and return reaction network dictionary
    Arguments:
        model_code: Int
            Integer to denote ROK model to use
        catalyst_mass: float
            Amount of catalyst loading in PFR, kg
        catalyst_density: float
            Catalyst density with void, kg/m^3
        bulk_density: float
            Catalyst density without void, kg/m^3
    """
    # Define maximum number of C-atoms in a molecule in reaction network
    nC = 9

    # list of defined species
    reacting_species = ["ethylene","propene","butene","pentene","hexene","heptene",
                        "octene","nonene"]

    c_number_list = [n for n in range(2,nC+1)]

    # Define counter for number of reactions in network
    counter = 0

    # Dictionary to store all reactions in network
    rate_reactions = {}

    # Read-in kinetic parameters
    kinetic_model_df = pd.read_csv('src/../data/fitted_parameters_M{}.csv'.format(model_code))
    kinetic_model_params_dict = {}
    for i, row in kinetic_model_df.iterrows():
        kinetic_model_params_dict[row[0]] = row[1]
    
    # Create the reaction property definition which describes the system on
    # reactions to be modeled.
    for i, Ci in enumerate(c_number_list):
        for j, Cj in enumerate(c_number_list):
            if Ci + Cj <= nC:
                ''' Forward oligomerization reaction: Ci + Cj --> Ci+j'''
                # dictionary to keep track of reactions by reactant/product index
                i_j_dict = {}
                # dictionary to keep track of reaction stoichiometries
                stoichiometry_dict = {}
                stoich_dict = {key:0 for key in reacting_species} # initialize
                for k, key in enumerate(stoich_dict):
                    if k + 2 == Ci:
                        if ("Vap", key) in stoichiometry_dict.keys():
                            stoichiometry_dict[("Vap", key)] = -2 # Eg.: C2+C2 -> C4
                        else:
                            stoichiometry_dict[("Vap", key)] = -1 # -1 for reactant
                    if k + 2 == Cj:
                        if ("Vap", key) in stoichiometry_dict.keys():
                            stoichiometry_dict[("Vap", key)] = -2 # Eg.: C2+C2 -> C4
                        else:
                            stoichiometry_dict[("Vap", key)] = -1 # -1 for reactant
                    if k + 2 == Ci + Cj:
                        stoichiometry_dict[("Vap", key)] = 1 # +1 for product
                i_j_dict["stoichiometry"] = stoichiometry_dict
                i_j_dict["heat_of_reaction"] = constant_dh_rxn
                i_j_dict["rate_constant"] = CoupledOligomerizationCrackingReactionClass
                i_j_dict["rate_form"] = power_law_rate
                i_j_dict["concentration_form"]= ConcentrationForm.partialPressure
                if model_code == 2 or model_code == 3:
                    i_j_dict["parameter_data"] = {
                                                "model_code":model_code,
                                                "dh_rxn_ref": (1.24e5),
                                                "alpha_olig": (kinetic_model_params_dict['log_alpha_olig']),
                                                "alpha_crack": (kinetic_model_params_dict['log_alpha_crack']),
                                                "beta_olig": (kinetic_model_params_dict['beta_olig']),
                                                "beta_crack": (kinetic_model_params_dict['beta_crack']),
                                                "gamma_olig": (kinetic_model_params_dict['gamma_olig']),
                                                "gamma_crack": (kinetic_model_params_dict['gamma_crack']),
                                                "phi_crack": (kinetic_model_params_dict['phi_crack']),
                                                "alpha_ads": (kinetic_model_params_dict['alpha_ads']),
                                                "beta_ads": (kinetic_model_params_dict['beta_ads']),
                                                "Ea_olig": (kinetic_model_params_dict['Ea_olig']),
                                                "Ea_crack": (kinetic_model_params_dict['Ea_crack']),
                                                "A_phys": (kinetic_model_params_dict['A_phys']),
                                                "A_des": (kinetic_model_params_dict['A_des']),
                                                "C_n": (Ci),
                                                "C_m": (Cj),
                                                "r_type": ('oligomerization'),
                                                "catalyst_mass":(catalyst_mass),
                                                "catalyst_density":(catalyst_density),
                                                "bulk_density":(bulk_density)}

                elif model_code == 4 or model_code == 5:
                    delH_formation_n = kinetic_model_params_dict['gamma'] * Ci + kinetic_model_params_dict['delta']
                    delH_formation_m = kinetic_model_params_dict['gamma'] * Cj + kinetic_model_params_dict['delta']
                    delH_formation_nm = kinetic_model_params_dict['gamma'] * (Ci + Cj) + kinetic_model_params_dict['delta']
                    i_j_dict["parameter_data"] = {
                                                "model_code":model_code,
                                                "dh_rxn_ref": (delH_formation_nm - delH_formation_m - delH_formation_n),
                                                "alpha_olig": (kinetic_model_params_dict['log_alpha_olig']),
                                                "alpha_crack": (kinetic_model_params_dict['log_alpha_crack']),
                                                "gamma": (kinetic_model_params_dict['gamma']),
                                                "delta": (kinetic_model_params_dict['delta']),
                                                "E0": (kinetic_model_params_dict['E0']),
                                                "alpha_ads": (kinetic_model_params_dict['alpha_ads']),
                                                "beta_ads": (kinetic_model_params_dict['beta_ads']),
                                                "kappa_olig": (kinetic_model_params_dict['kappa_olig']),
                                                "kappa_crack": (kinetic_model_params_dict['kappa_crack']),
                                                "A_phys": (kinetic_model_params_dict['A_phys']),
                                                "A_des": (kinetic_model_params_dict['A_des']),
                                                "C_n": (Ci),
                                                "C_m": (Cj),
                                                "r_type": ('oligomerization'),
                                                "catalyst_mass":(catalyst_mass),
                                                "catalyst_density":(catalyst_density),
                                                "bulk_density":(bulk_density)}
                else:
                    pass
                # assign reaction details to reaction dictionary
    #             if Ci == 3 and Cj == 3:
                rate_reactions['R{}'.format(counter)] = i_j_dict

                counter += 1
                
    for i, Ci in enumerate(c_number_list):
        for j, Cj in enumerate(c_number_list):
            if Ci + Cj <= nC:
                ''' Backward cracking reaction: Ci+j --> Ci + Cj'''
                # dictionary to keep track of reactions by reactant/product index
                i_j_dict = {}
                # dictionary to keep track of reaction stoichiometries
                stoichiometry_dict = {}
                stoich_dict = {key:0 for key in reacting_species} # initialize
                for k, key in enumerate(stoich_dict):
                    if k + 2 == Ci:
                        if ("Vap", key) in stoichiometry_dict.keys():
                            stoichiometry_dict[("Vap", key)] = 2 # +2 for 2 same products
                        else:
                            stoichiometry_dict[("Vap", key)] = 1 # +1 for product
                    if k + 2 == Cj:
                        if ("Vap", key) in stoichiometry_dict.keys():
                            stoichiometry_dict[("Vap", key)] = 2 # +2 for 2 same products
                        else:
                            stoichiometry_dict[("Vap", key)] = 1 # +1 for product
                    if k + 2 == Ci + Cj:
                        stoichiometry_dict[("Vap", key)] = -1 # -1 for reactant
                i_j_dict["stoichiometry"] = stoichiometry_dict
                i_j_dict["heat_of_reaction"] = constant_dh_rxn
                i_j_dict["rate_constant"] = CoupledOligomerizationCrackingReactionClass
                i_j_dict["rate_form"] = power_law_rate
                i_j_dict["concentration_form"]= ConcentrationForm.partialPressure
                if model_code == 2 or model_code == 3:
                    i_j_dict["parameter_data"] = {
                                                "model_code":model_code,
                                                "dh_rxn_ref": (1.24e5),
                                                "alpha_olig": (kinetic_model_params_dict['log_alpha_olig']),
                                                "alpha_crack": (kinetic_model_params_dict['log_alpha_crack']),
                                                "beta_olig": (kinetic_model_params_dict['beta_olig']),
                                                "beta_crack": (kinetic_model_params_dict['beta_crack']),
                                                "gamma_olig": (kinetic_model_params_dict['gamma_olig']),
                                                "gamma_crack": (kinetic_model_params_dict['gamma_crack']),
                                                "phi_crack": (kinetic_model_params_dict['phi_crack']),
                                                "alpha_ads": (kinetic_model_params_dict['alpha_ads']),
                                                "beta_ads": (kinetic_model_params_dict['beta_ads']),
                                                "Ea_olig": (kinetic_model_params_dict['Ea_olig']),
                                                "Ea_crack": (kinetic_model_params_dict['Ea_crack']),
                                                "A_phys": (kinetic_model_params_dict['A_phys']),
                                                "A_des": (kinetic_model_params_dict['A_des']),
                                                "C_n": (Ci),
                                                "C_m": (Cj),
                                                "r_type": ('cracking'),
                                                "catalyst_mass":(catalyst_mass),
                                                "catalyst_density":(catalyst_density),
                                                "bulk_density":(bulk_density)}

                elif model_code == 4 or model_code == 5:
                    delH_formation_n = kinetic_model_params_dict['gamma'] * Ci + kinetic_model_params_dict['delta']
                    delH_formation_m = kinetic_model_params_dict['gamma'] * Cj + kinetic_model_params_dict['delta']
                    delH_formation_nm = kinetic_model_params_dict['gamma'] * (Ci + Cj) + kinetic_model_params_dict['delta']
                    i_j_dict["parameter_data"] = {
                                                "model_code":model_code,
                                                "dh_rxn_ref": (delH_formation_m + delH_formation_n - delH_formation_nm),
                                                "alpha_olig": (kinetic_model_params_dict['log_alpha_olig']),
                                                "alpha_crack": (kinetic_model_params_dict['log_alpha_crack']),
                                                "gamma": (kinetic_model_params_dict['gamma']),
                                                "delta": (kinetic_model_params_dict['delta']),
                                                "E0": (kinetic_model_params_dict['E0']),
                                                "alpha_ads": (kinetic_model_params_dict['alpha_ads']),
                                                "beta_ads": (kinetic_model_params_dict['beta_ads']),
                                                "kappa_olig": (kinetic_model_params_dict['kappa_olig']),
                                                "kappa_crack": (kinetic_model_params_dict['kappa_crack']),
                                                "A_phys": (kinetic_model_params_dict['A_phys']),
                                                "A_des": (kinetic_model_params_dict['A_des']),
                                                "C_n": (Ci),
                                                "C_m": (Cj),
                                                "r_type": ('cracking'),
                                                "catalyst_mass":(catalyst_mass),
                                                "catalyst_density":(catalyst_density),
                                                "bulk_density":(bulk_density)
                                                }
                
                else:
                    pass
                # assign reaction details to reaction dictionary
    #             if Ci == 3 and Cj == 3:
                rate_reactions['R{}'.format(counter)]= i_j_dict

                counter += 1
    return rate_reactions