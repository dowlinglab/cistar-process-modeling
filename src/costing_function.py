from idaes.core import UnitModelBlock, UnitModelCostingBlock
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
from idaes.models.properties import iapws95

from idaes.models.costing.SSLW import (
    SSLWCosting,
    SSLWCostingData,
    HXMaterial,
    HXTubeLength,
    HXType,
    VesselMaterial,
    TrayType,
    TrayMaterial,
    HeaterMaterial,
    HeaterSource,
    PumpType,
    PumpMaterial,
    PumpMotorType,
    CompressorDriveType,
    CompressorMaterial,
    CompressorType,
)

from pyomo.environ import Param, Expression, Set, Objective, Var, Constraint, exp, log

from pyomo.environ import units as pyunits

def calculate_H102_cooler_cap_cost(blk):
    """
    Function to calculate capital cost of a cooler
    Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of 
    Process and Product Design Principles: Synthesis, Analysis, and Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
    Cooling fluid is water entering at 16 C and leaving at 40 C
    Heat transfer coefficient estimates obtained from:
    https://www.engineersedge.com/heat_transfer/overall_heat_transfer_coefficients_13827.htm#:~:text=Engineering%20Physics&text=Note%20that%20the%20overall%20heat,exchangers%20that%20involve%20phase%20changes.
    """
    blk.LMTD = Param(initialize=10.0, mutable=True)
    blk.heat_transfer_coeff = Param(initialize=100.0, mutable=True) # W/m^2/K, assumption

    blk.area = Expression(expr=-blk.heat_duty[0]/(blk.heat_transfer_coeff*blk.LMTD)) # m^2

    if blk.area()*10.7639 >= 15000.0: # Air fin-fan
        blk.capital_cost = Expression(expr=2500 * (blk.area*10.7639)**0.40, doc="Air fin-fan") #USD
    elif blk.area()*10.7639 >= 2000.0: # Plate-and-frame, 150 - 15000 ft^2
        blk.capital_cost = Expression(expr=8880 * (blk.area*10.7639)**0.42, doc="Plate-and-frame") #USD
    elif blk.area()*10.7639 >= 500.0: # Spiral plate, 20 - 2000 ft^2
        blk.capital_cost = Expression(expr=6200 * (blk.area*10.7639)**0.42, doc="spiral plate") #USD
    else: # Spiral tube
        blk.capital_cost = Expression(expr=exp(8.0757 + 0.4343* log(blk.area*10.7639) + 0.03812*(log(blk.area*10.7639)**2)), doc="Spiral tube") #USD
    return 

def calculate_H106_cooler_cap_cost(blk):
    """
    Function to calculate capital cost of a cooler
    Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of 
    Process and Product Design Principles: Synthesis, Analysis, and Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
    Cooling fluid is water entering at 16 C and leaving at 40 C
    """
    blk.LMTD = Param(initialize=10.0, mutable=True)
    blk.heat_transfer_coeff = Param(initialize=100.0, mutable=True) # W/m^2/K, assumption

    blk.area = Expression(expr=-blk.heat_duty[0]/(blk.heat_transfer_coeff*blk.LMTD)) # m^2
    if blk.area()*10.7639 >= 15000.0: # Air fin-fan
        blk.capital_cost = Expression(expr=2500 * (blk.area*10.7639)**0.40, doc="Air fin-fan") #USD
    elif blk.area()*10.7639 >= 2000.0: # Plate-and-frame, 150 - 15000 ft^2
        blk.capital_cost = Expression(expr=8880 * (blk.area*10.7639)**0.42, doc="Plate-and-frame") #USD
    elif blk.area()*10.7639 >= 500.0: # Spiral plate, 20 - 2000 ft^2
        blk.capital_cost = Expression(expr=6200 * (blk.area*10.7639)**0.42, doc="spiral plate") #USD
    else: # Spiral tube
        blk.capital_cost = Expression(expr=exp(8.0757 + 0.4343* log(blk.area*10.7639) + 0.03812*(log(blk.area*10.7639)**2)), doc="Spiral tube") #USD
    return 

def calculate_H104_H105_cooler_cap_cost(blk1,blk2):
    """
    Function to calculate capital cost of a cooler
    Formulae adapted from table 22.32 (P. 592, heat exchangers, other) of 
    Process and Product Design Principles: Synthesis, Analysis, and Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment

    If 2 heater models are used to replicate 1 Heat Exchanger, this function is used to calculate
    the capital cost of the combined heaters and the capital cost is assigned to blk1
    """
    blk1.LMTD = Param(initialize=10.0, mutable=True)
    blk1.heat_transfer_coeff = Param(initialize=100.0, mutable=True) # W/m^2/K, assumption
    blk1.area = Expression(expr=-(blk1.heat_duty[0] + blk2.heat_duty[0])/(blk1.heat_transfer_coeff*blk1.LMTD)) # m^2

    if blk1.area()*10.7639 >= 15000.0: # Air fin-fan
        blk1.capital_cost = Expression(expr=2500 * (blk1.area*10.7639)**0.40, doc="Air fin-fan") #USD
    elif blk1.area()*10.7639 >= 2000.0: # Plate-and-frame
        blk1.capital_cost = Expression(expr=8880 * (blk1.area*10.7639)**0.42, doc="Plate-and-frame") #USD
    elif blk1.area()*10.7639 >= 500.0: # Spiral plate
        blk1.capital_cost = Expression(expr=6200 * (blk1.area*10.7639)**0.42, doc="spiral plate") #USD
    else: # Spiral tube
        blk1.capital_cost = Expression(expr=exp(8.0757 + 0.4343* log(blk1.area*10.7639) + 0.03812*(log(blk1.area*10.7639)**2)), doc="Spiral tube") #USD
    return

def add_costing(m):
    """
    Main function to calculate unit model CAPEX
    """

#     assert degrees_of_freedom(m) == 0

    m.fs.costing = SSLWCosting()

    # Computing heater capital costs
    # loop over heater units
    for unit in [m.fs.H101, m.fs.H103]:
        unit.number_of_units = 1
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.parent_block().costing,
            costing_method=SSLWCostingData.cost_fired_heater,
            costing_method_arguments={
                    "material_type": HeaterMaterial.CarbonSteel,
                    "heat_source": HeaterSource.Fuel,
                }
        )
    # m.fs.H101.cost_heater = Expression(
    #     expr=0.036158*m.fs.H101.heat_duty[0] + 63931.475*pyunits.W,
    #     doc='capital cost of heater in $')
    
    # Computing reactor capital cost
    # bulk density of dehydro catalyst is similar to oligo catalyst
    # https://www.researchgate.net/profile/Maryam-Nikoo-2/publication/269872080_3_rd_Determination_of_Operational_Condition_for_Propane_Dehydrogenation_over_a_Commercial_Pt-_Sn-KAl_2_O_3_Catalyst/links/549824920cf2eeefc30f76cc/3-rd-Determination-of-Operational-Condition-for-Propane-Dehydrogenation-over-a-Commercial-Pt-Sn-K-Al-2-O-3-Catalyst.pdf
    # capital cost is 2.5 of oligo -> volume 2.5 times
    # Diameter -> 1.58 times
    m.fs.R101.diameter = (2.5*m.fs.R102.area/3.14)**0.5
    
    m.fs.R101.length = 1 * pyunits.m
    m.fs.R101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing)
    
    # Cooler H102 capital cost
    calculate_H102_cooler_cap_cost(m.fs.H102)
    
    # Expander capital cost
    
    m.fs.E101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
        )

    # Computing reactor capital cost
    m.fs.R102.diameter = (m.fs.R102.area/3.14)**0.5 # get area from model
    m.fs.R102.length = 1 * pyunits.m
    m.fs.R102.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing)
    
    
    # Cooler H104 and H105 costs assigned to H104
    calculate_H104_H105_cooler_cap_cost(m.fs.H104,m.fs.H105)
    
    # Computing flash capital cost
    m.fs.F101.diameter = Param(initialize=3.5, units=pyunits.m)
    m.fs.F101.length = Param(initialize=14, units=pyunits.m)
#     m.fs.F101.diameter = 3.5 * pyunits.m
#     m.fs.F101.length = 14 * pyunits.m
    m.fs.F101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing)

    # Cooler H106 capital cost
    calculate_H106_cooler_cap_cost(m.fs.H106)

    # Computing flash capital cost
    m.fs.F102.diameter = Param(initialize=3.5, units=pyunits.m)
    m.fs.F102.length = Param(initialize=14, units=pyunits.m)
#     m.fs.F102.diameter = 3.5 * pyunits.m
#     m.fs.F102.length = 14 * pyunits.m
    m.fs.F102.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing)

    # compressor 1
    m.fs.C101.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

    # compressor 2
    m.fs.C102.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing
    )

def calculate_costs_for_objective(m,c_tax_flag = True, c_tax_val = 0.0, include_HI=True):
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
    # Annualizing capital cost to same scale as operating costs (per year)
    m.fs.capex = Expression(
        expr=(m.fs.H101.costing.capital_cost
              + m.fs.R101.costing.capital_cost
              + m.fs.H102.capital_cost
              + m.fs.E101.costing.capital_cost
              + m.fs.H103.costing.capital_cost
              + m.fs.R102.costing.capital_cost
              + m.fs.H104.capital_cost
              + m.fs.F101.costing.capital_cost
              + m.fs.H106.capital_cost
              + m.fs.F102.costing.capital_cost
              + m.fs.C101.costing.capital_cost
              + m.fs.C102.costing.capital_cost))
    
    # utilities
    m.fs.cost_ngl = Param(initialize = 8.53, mutable=True) # USD/GJ NGL cost, 9 $/MMBtu (average price from 2021), 1 MJ = 0.000947817 MMBtu
    m.fs.cost_heating = Param(initialize = 5.34, mutable=True) # USD/GJ, natural gas for heating, 5.63 USD/MMBtu (average price for 2023)
    m.fs.cost_cooling = Param(initialize = 0.354, mutable=True) # USD/GJ
    m.fs.cost_electricity = Param(initialize = 42.6, mutable=True) # USD/GJ, 15.32 c/kWh (average across the US 2023)
    m.fs.c_tax_rate = Param(initialize = c_tax_val,mutable=True) ## [USD/kg CO2_e], $35/tonne CO2e [0, 17, 45, 190, 410] USD/ton CO2_e; 190, 410 - EPA proposal; 17, 45 - paper
    
    m.fs.RMC = Expression(expr = 365 * 24 * m.fs.feed_energy*m.fs.cost_ngl) # Yearly raw material cost [GJ/h]*[USD/GJ]

    if include_HI:
        m.fs.cooling_cost = Expression(
            expr=m.fs.cost_cooling * m.fs.Qw/3600) # heating demand in USD/s
        m.fs.heating_cost = Expression(
            expr=m.fs.cost_heating * (m.fs.Qs - m.fs.purge_heat)/3600) # heating demand in USD/s
        m.fs.heating_cost_constraint = Constraint(expr = m.fs.Qs >= m.fs.purge_heat) # constraint to make sure heating cost is non-negative
    else:
        m.fs.cooling_cost = Expression(
            expr=-m.fs.cost_cooling * (m.fs.H102.heat_duty[0] + m.fs.H104.heat_duty[0] +
                                      m.fs.H105.heat_duty[0] + m.fs.H106.heat_duty[0] +
                                      sum(m.fs.R102.control_volume.heat[0.0,i] for 
                                      i in m.fs.R102.control_volume.length_domain if i != 0)*0.05)*
                                      0.001*0.001*0.001) # heating demand in USD/s
        m.fs.heating_cost = Expression(
            expr=m.fs.cost_heating * ((m.fs.H101.heat_duty[0] + m.fs.H103.heat_duty[0] + 
                                       m.fs.H107.heat_duty[0] + m.fs.R101.heat_duty[0])*
                                       0.001*0.001*0.001 - m.fs.purge_heat/3600)) # heating demand in USD/s

    # Expression to compute the total electricity cost (utilities - credit), [USD/s]
    m.fs.electricity_cost = Expression(
        expr=m.fs.cost_electricity * m.fs.total_electricity/3600)
    
    # H2-membrane cost
    m.fs.H2_recovery_rate = Param(initialize=0.105, mutable=True) # kmol/h/m^2
    m.fs.membrane_unit_cost = Param(initialize = 50, mutable=True) # $/m^2
    m.fs.H2_membrane_cost = Expression(expr = m.fs.S101.H2_purge.flow_mol[0]*m.fs.S101.H2_purge.mole_frac_comp[0,'hydrogen'] * m.fs.membrane_unit_cost/m.fs.H2_recovery_rate) # area calculated from 0.105 kmol/h/m^2 recovery of H2
    m.fs.H2_purge_flow_kg = Expression(expr= m.fs.S101.H2_purge.flow_mol[0]*m.fs.S101.H2_purge.mole_frac_comp[0,'hydrogen']*m.fs.Vap_props.hydrogen.mw+
                                             m.fs.S101.H2_purge.flow_mol[0]*m.fs.S101.H2_purge.mole_frac_comp[0,'methane']*m.fs.Vap_props.methane.mw+
                                             m.fs.S101.H2_purge.flow_mol[0]*m.fs.S101.H2_purge.mole_frac_comp[0,'ethane']*m.fs.Vap_props.ethane.mw+
                                             m.fs.S101.H2_purge.flow_mol[0]*m.fs.S101.H2_purge.mole_frac_comp[0,'ethylene']*m.fs.Vap_props.ethylene.mw) # kg/s
    m.fs.H2_purge_unit_sell_price = Param(initialize = 1.0, mutable=True) # $/kg, general rate around $2/kg, discounted for lower purity

    # Catalyst cost from https://www.acsmaterial.com/zsm-5-catalyst.html
    m.fs.olig_catalyst_unit_cost = Param(initialize = 336, mutable=True) # $/kg
    
    m.fs.dehydro_catalyst_cost = Expression(expr = m.fs.olig_catalyst_unit_cost * 2.5 * m.fs.reactions.catalyst_mass) # dehydro reactor uses 2.5 times the catalyst used by the oligo reactor
    m.fs.oligo_catalyst_cost = Expression(expr = m.fs.olig_catalyst_unit_cost * m.fs.reactions.catalyst_mass)
    
    # annual carbon emissions tax
    m.fs.CO2_tax = Expression(expr = 24 * 365 * (m.fs.em_upstream + m.fs.em_process)*m.fs.c_tax_rate) # USD/year

    # annual H2 purge sell price
    m.fs.H2_purge_sell_price = Expression(expr = m.fs.H2_purge_unit_sell_price*m.fs.H2_purge_flow_kg*3600*24*365) # USD/year
    
    # Expression to compute the total operating cost
    m.fs.opex = Expression(
        expr=(
            3600 * 24 * 365 * (m.fs.heating_cost + m.fs.cooling_cost + m.fs.electricity_cost)
            + m.fs.H2_membrane_cost
            + m.fs.dehydro_catalyst_cost
            + m.fs.oligo_catalyst_cost))
    
    m.fs.capex_factor = Param(initialize = 0.18, mutable=True) # construction, auxilliary pipelines, etc.
    m.fs.opex_factor = Param(initialize = 1.23, mutable=True) # maintenance, personnel, etc.
    
    m.fs.OPEX_tot = Expression(expr = m.fs.capex_factor*m.fs.capex + m.fs.opex_factor*(m.fs.RMC + m.fs.opex))## $/year

    m.fs.dr = Param(initialize = 0.04,mutable=True) ## interest rate per year
    m.fs.n_year = Param(initialize = 20.0, mutable=True) ## lifetime
    m.fs.AF = Expression(expr = (m.fs.dr*(1.0+m.fs.dr)**m.fs.n_year)/(((1.0+m.fs.dr)**m.fs.n_year)-1.0)) ## Annualizing factor

    if c_tax_flag:
        m.fs.TAC = Expression(expr = m.fs.OPEX_tot + m.fs.AF*m.fs.capex + m.fs.CO2_tax - m.fs.H2_purge_sell_price) # [USD/year]
    else:
        m.fs.TAC = Expression(expr = m.fs.OPEX_tot + m.fs.AF*m.fs.capex - m.fs.H2_purge_sell_price) # [USD/year]
    # economic objective
    m.fs.min_sell_price = Var(initialize = 10.0, bounds=(1e-4, None))
    m.fs.min_sell_price_constraint = Constraint(expr = m.fs.min_sell_price * (m.fs.fuel_energy*24*365) == m.fs.TAC) #[USD/GJ] ## OBJECTIVE FUNCTION
    # m.fs.min_sell_price = Expression(expr = m.fs.TAC/(m.fs.fuel_energy*24*365)) #[USD/GJ] ## OBJECTIVE FUNCTION
