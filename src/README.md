## Folders
### state_properties
This folder contains python files that contain properties of liquids and gases used in this work:
* Vapor-phase properties: properties_vap.py
* Vapor-phase properties for H2-rich gas mixture: properties_vap_H2_permeate.py
* Vapor-phase properties without H2 and CH4: properties_vap_post_flash_2.py
* Vapor-liquid phase properties: properties_VLE_FpcTP.py

## Python files

### Declaring oligomerization kinetics: rate_constant_custom_oligomerization.py, rate_forms_oligo.py, and reaction_network_generator.py

### Declaring simple dehydrogenation kinetics: dehydro_reactions.py

#### Declaring process and unit models: unit_initialization.py
This file includes all functions used to define and initialize process unit models in this implementation.

#### Calculating GHG emissions: emissions_calculations.py
This file includes functions to calculate GHG emissions associated with the process

#### Calculating process-related costs: costing_function.py
This file includes functions to calculate process costs such as CAPEX, OPEX, and MSP

#### Plotting functions: plotting_functions.py
This file includes all functions used to generate plots using matplotlib

All functions take the flowsheet model (m) as the first argument and return `None` by default unless specified otherwise.

#### Functions defined in unit_initialization.py

##### create_flowsheet()
    """
    Function to create and return ConcreteModel flowsheet
    """

##### create_properties()
    """
    Function to create stream properties
    Arguments:
	m: int
            Integer to identify ROK model. Possible values: 2,3,4,5
    """

##### define_models()
    """
    Function to define unit models and properties
    Arguments:
        catalyst_mass: float
            Mass of catalyst loading in oligomerization reactor, kg
    """

##### define_arcs()
    """
    Function to define connections using arcs
    """

##### set_unit_model_variables()
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

##### replace_heater_heat_duty_constraint_with_bounds()
    """
    Helper function to check and remove inequality constraints on heater/cooler heat_duty and replace with bounds using setlb() and setub()
    """

##### update_model_for_optimization()
    """
    Helper function to:
        - Add/update constraints and update DOF bounds for optimization
        - Add objective function for optimization problem
    Arguments:
        obj_reset_flag: boolean, default=False
            True: Delete and re-define objective function to ensure previous solve does not affect current solve
            False: Use existing objective function value
    """

##### unfix_DOFs_pre_optimization()
    """
    Helper function to unfix optimization degrees of freedom before optimization
    """

##### fix_DOFs_post_optimization()
    """
    Helper function to fix optimization degrees of freedom post-optimization
    """

##### vapor_only_to_vapor_liquid_reformulate()
    """
    Function to check if translator block with vapor to vapor-liquid flow has zero flow in vapor or liquid phase. If either phase flow is zero, alter model equations to remove degenerate zero-flow phase composition equations and fix values to 1e-8
    Arguments:
        blk: idaes Translator block model object
             to investigate and alter
        error_check: boolean
            True (default): ensure if all liquid or vapor flowrates are zero
            False: skip check, no need to alter model equations
        zero_liquid_outlet: boolean
            True (default): when translator outlet is supposed to have only liquid flow
            False: If translator outlet has zero vapor flow (provisional for future)
    """

##### H106_inlet_vapor_reformulate()
    """
    Function to check if unit block with vapor-liquid flow has zero flow in vapor phase. If vapor phase flow is zero, alter model equations to remove degenerate zero-flow phase composition equations and fix values to 1e-8
    Arguments:
        blk: idaes unit block model object
             block to investigate and alter
        error_check: boolean
            True (default): ensure if all vapor flowrates are zero
            False: skip check, no need to alter model equations
    """

##### update_model_after_initialization()
    """
    Function to add additional constraints to model and unfix corresponding variables after initialization
    """

##### create_model_and_return_optimization_results()
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

##### initialize_flowsheet()
    """
    Function to initialize process flowsheet
    Arguments:
        tear_guesses: dict
            Dictionary with guess values for tear stream
    """

##### set_scaling_factors()
    """
    Function to set scaling for model components
    """


#### Functions defined in costing_function.py

##### calculate_H102_cooler_cap_cost()
    """
    Function to calculate capital cost of a cooler. Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of Process and Product Design Principles: Synthesis, Analysis, and Evaluation. Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons Chapter 22. Cost Accounting and Capital Cost Estimation 22.2 Cost Indexes and Capital Investment Cooling fluid is water entering at 16 C and leaving at 40 C. Heat transfer coefficient estimates obtained from: https://www.engineersedge.com/heat_transfer/overall_heat_transfer_coefficients_13827.htm#:~:text=Engineering\%20Physics&text=Note\%20that\%20the\%20overall\%20heat,exchangers\%20that\%20involve\%20phase\%20changes.
    """

##### calculate_H106_cooler_cap_cost()
    """
    Function to calculate capital cost of a cooler. Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of Process and Product Design Principles: Synthesis, Analysis, and Evaluation. Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons Chapter 22. Cost Accounting and Capital Cost Estimation 22.2 Cost Indexes and Capital Investment Cooling fluid is water entering at 16 C and leaving at 40 C
    """

##### calculate_H104_H105_cooler_cap_cost()
    """
    Function to calculate capital cost of a cooler. Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of Process and Product Design Principles: Synthesis, Analysis, and Evaluation. Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons Chapter 22. Cost Accounting and Capital Cost Estimation 22.2 Cost Indexes and Capital Investment.
    Special case: If 2 heater models are used to replicate 1 Heat Exchanger, this function is used to calculate the capital cost of the combined heaters and the capital cost is assigned to first Heat Exchanger in flowsheet.
    """

##### add_costing()
    """
    Main function to calculate unit model CAPEX
    """

##### calculate_costs_for_objective()
    """
    Function to calculate utility and CAPEX, OPEX costs for MSP objective function calculation
    Arguments:
        c_tax_flag: bool
            Flag to check if CO2 emissions need to be taxed
            True (default): Emissions taxed
            False: Emission not taxed
        c_tax_val: float
            Emisisons tax rate, USD/kg CO2e
        include_HI: bool
            Flag to denote if heat integration calculations are included in model formulation
            True (default): Heat integration calculations included in problem formulation
            False: Heat integration calculations not included in problem formulation
    """

#### Functions defined in emissions_calculation.py

##### calc_lhv_values()
    """
    Function to calculate LHV values
    Arguments:
        case_name: str
            Name of region (case) whose NGL composition is considered
        lhv_values_file: str
            Filename for .xlsx file containing data for LHV value calculations
        ngl_compositions_file: str
            Filename for .csv file containing NGL compositions for all cases considered
        ngl_fractions_file: str
            Filename for .csv file containing data for fraction of shale gas that is NGL for all cases considered
    """

##### calculate_stream_energies()
    """
    Function to calculate energy content (LHV) of each inlet and outlet stream of the process
    """

##### calculate_emissions()
    """
    Function to calculate upstream and downstream emissions
    Arguments:
        case_name: str
            Name of region (case) whose NGL composition is considered
        emissions_factor_file: str
            Filename for .csv file with emissions factors for upstream emissions calculation for all cases considered
        include_HI: bool
            Flag to denote if heat integration calculations are included in model formulation
            True (default): Heat integration calculations included in problem formulation
            False: Heat integration calculations not included in problem formulation
    """

##### create_ghg_objective()
    """
    Function to calculate GHG objective expression
    """
    
##### delete_region_specific_components()
    """
    Function to delete region-specific emisisons data in order to update emissions calculations from a single initialization file (single region)
    """

#### Functions defined in plotting_functions.py

##### generate_and_save_composite_curves()
    """
    Function for plotting composite curves
    Arguments:
        CD: Class containing curve data
        model_code: ROK model code
        C_tax_rate: CO2 emissions tax rate
        region: shale region
        optimal_solution: bool
            True: data is for optimal solution
            False: data is for initialized flowsheet
    """
    
##### plot_outlet_flowrate_by_ROK_model_horizontal_stacked_bars()
    """
    Function to create horizontal stacked plot of liquid hydrocarbon composition w.r.t. ROK model, mol/s
    Arguments:
        liquid_hc_component_flow_dict: dict
            Dictionary of liquid outlet flow rate per component, mol/s
        ROK_model_list: list
            List of ROK model code strings
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_outlet_flowrate_by_region_horizontal_stacked_bars()
    """
    Function to create horizontal stacked plot of liquid hydrocarbon composition w.r.t. shale region, mol/s
    Arguments:
        liquid_hc_component_flow_dict: dict
            Dictionary of liquid outlet flow rate per component, mol/s
        region_list: list
            List of shale gas region strings
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_C4_C9_olefin_outlet_flowrate_by_ROK_model_vertical_bars()
    """
    Function to create horizontal stacked plot of liquid hydrocarbon composition in mol/s w.r.t. ROK model
    Arguments:
        liquid_C4_C9_percent_dict: dict
                Dictionary of liquid outlet flow rate of C4+ alkenes, mol/s
        liquid_C4_C9_percent_dict: dict
            Dictionary of liquid outlet flow rate of C4+ alkenes, mol/s
        ROK_models_list: list
            List of ROK model code strings
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """
    
##### plot_outlet_flow_rate_by_component_vs_ROK_model_lines()
    """
    Function to create line plot of liquid hydrocarbon flow rates by component vs ROK models, mol/s
    Arguments:
        liquid_hc_component_flow_dict: dict
            Dictionary of liquid outlet flow rate per component, mol/s
        ROK_model_list: list
            List of ROK model code strings
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_outlet_flow_rate_by_component_vs_region_lines()
    """
    Function to create line plot of liquid hydrocarbon flow rates by component vs shale region, mol/s
    Arguments:
        liquid_hc_component_flow_dict: dict
            Dictionary of liquid outlet flow rate per component, mol/s
        region_list: list
            List of shale gas region strings
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_emissions_by_oligo_model_vertical_stacked_bars()
    """
    Function to create stacked bar plot of upstream and downstream emissions in g CO2e/MJ fuel w.r.t. model instance
    Arguments:
        upstream_emissions_dict: dict
            Dictionary of upstream emissions by region, g CO2e/MJ fuel
        downstream_emissions_dict: dict
            Dictionary of downstream emissions by region, g CO2e/MJ fuel
        region_list: list
            List of regions
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_emissions_by_region_vertical_stacked_bars()
    """
    Function to create stacked bar plot of upstream and downstream emissions in g CO2e/MJ fuel w.r.t. shale region
    Arguments:
        upstream_emissions_dict: dict
            Dictionary of upstream emissions by region, g CO2e/MJ fuel
        downstream_emissions_dict: dict
            Dictionary of downstream emissions by region, g CO2e/MJ fuel
        region_list: list
            List of regions
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_LHV_by_region_stacked_vertical_bars()
    """
    Function to create stacked bar plot of LHV by component vs region, MJ/kmol fuel w.r.t. model instance
    Arguments:
        LHV_dict: dict
            Dictionary of LHV by component for each region, MJ/kmol fuel
        region_list: list
            List of regions
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """

##### plot_MSP_vs_c_tax_rate()
    """
    Function to create stacked bar plot of LHV by component vs region, MJ/kmol fuel w.r.t. model instance
    Arguments:
        optimal_MSP_dict: dict
            Dictionary of MSP for each C-tax rate, USD/MJ fuel
        c_tax_rates: list
            List of C-tax rates considered
        xlabel_string: str
            String with x-axis label for plot
        ylabel_string: str
            String with y-axis label for plot
        file_name_string: str
            String of name to save pdf of plot generated
    """