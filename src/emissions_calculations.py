import pandas as pd
from pyomo.environ import Param, Expression, Set, Var, Constraint
import numpy as np

def calc_lhv_values(m,case_name,lhv_values_file,ngl_compositions_file,ngl_fractions_file):
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
    
    ### LHV calc
    file_data = lhv_values_file
    LHV_in = pd.read_excel(file_data,sheet_name=0) 

    LHV_header = ['deltaH','H20','C','H']
    LHV_values = LHV_in[LHV_header]
    LHV_values = LHV_values.to_numpy() # [kJ/mol], molecules, atoms, atoms.

    chem_sp = LHV_in['Name']
    chem_sp = chem_sp.to_numpy()

    nsp = len(chem_sp) # number of chem sp.

    vap_h2o = 44.002 # kJ/mol

    LHV = np.zeros(nsp)
    MW = np.zeros(nsp)

    for i in range(nsp):
        LHV[i] = (LHV_values[i,0] - vap_h2o*LHV_values[i,1])/1000 # MJ/mol
        MW[i] = LHV_values[i,2]*12.0107 + LHV_values[i,3]*1.00784

    data = np.stack((LHV, MW), axis=-1)
    data = data.T
    df_LHV_MW = pd.DataFrame(data,columns=chem_sp)
    df_LHV_MW

    m.fs.component_LHV_values = Param(m.fs.vapor,initialize=df_LHV_MW.loc[0,:]) # MJ/mol
    m.fs.component_MW_values = Param(m.fs.vapor,initialize=df_LHV_MW.loc[1,:])

    ngl_inlet_df = pd.read_csv(ngl_compositions_file) # feed flow to process
    # well_comp_frac_df = pd.read_csv('shale_compositions.csv')
    well_comp_frac_df = pd.read_csv(ngl_fractions_file) # fraction of shale flow that is feed to process

    def get_composition_dict(comp_df,case_name):
        comp_dict = {}
        for i, r in comp_df.iterrows():
            if r[case_name] == 0.0:
                comp_dict[r['Species']] = 1e-6
            else:
                comp_dict[r['Species']] = r[case_name]
        return comp_dict

    def get_shale_and_NG_compositions(comp_df,frac_df):
        shale_comp_dict = {}
        shale_flow_by_comp = {}
        ng_comp_dict = {}
        ng_flow_by_comp = {}
        for i,r in comp_df.iterrows():
            if r[case_name] == 0.0:
                ng_flow_by_comp[r['Species']] = 1e-6
                shale_flow_by_comp[r['Species']] = 1e-6
            else:
                shale_flow_by_comp[r['Species']] = m.fs.M101.feed.flow_mol[0]()*r[case_name]/frac_df[case_name].iloc[0]
                ng_flow_by_comp[r['Species']] = m.fs.M101.feed.flow_mol[0]()*r[case_name]*(1/frac_df[case_name].iloc[0] - 1)
        total_shale_flow = 1342.778 # mol/s raw shale feed, data from Ridha et al. 2018
        total_ng_flow = total_shale_flow * (1.0 - frac_df[case_name].iloc[0])
        for k,v in shale_flow_by_comp.items():
            ng_comp_dict[k] = ng_flow_by_comp[k]/total_ng_flow
            shale_comp_dict[k] = v/total_shale_flow
        
        return shale_flow_by_comp, shale_comp_dict, total_shale_flow,\
                ng_flow_by_comp, ng_comp_dict, total_ng_flow

    shale_flow_by_comp, shale_comp_dict, total_shale_flow,\
    NG_flow_by_comp, NG_comp_dict, total_NG_flow = get_shale_and_NG_compositions(ngl_inlet_df,well_comp_frac_df)

    m.fs.shale_feed_comp = Param(m.fs.vapor, initialize = shale_comp_dict)
    m.fs.shale_flowrate = Param(initialize = total_shale_flow) # mol/s
    m.fs.shale_flow_by_comp = Param(m.fs.vapor,initialize = shale_flow_by_comp) # mol/s

    m.fs.NG_feed_comp = Param(m.fs.vapor, initialize = NG_comp_dict)
    m.fs.NG_flowrate = Param(initialize = total_NG_flow) # mol/s
    m.fs.NG_flow_by_comp = Param(m.fs.vapor,initialize = NG_flow_by_comp) # mol/s

    # m.fs.OutletStreams = Set(initialize = ['feed','liq_outlet','purge'])
    m.fs.Streams = Set(initialize = ['shale','NG','feed','H2_vent','liq_outlet','purge'])


    def return_LHV_per_component_per_stream(blk,j,i):
        if j == 'shale':
            return blk.shale_feed_comp[i]*blk.component_LHV_values[i]
        elif j == 'NG':
            return blk.NG_feed_comp[i]*blk.component_LHV_values[i]
        elif j == 'feed':
            return blk.M101.feed.mole_frac_comp[0,i]*blk.component_LHV_values[i]
        elif j == 'H2_vent':
            return blk.S101.H2_purge.mole_frac_comp[0,i]*blk.component_LHV_values[i]
        elif j == 'liq_outlet':
            if i in m.fs.liquid:
                return blk.F102.liq_outlet.flow_mol_phase_comp[0.0, 'Liq', i]*blk.component_LHV_values[i]/sum(blk.F102.liq_outlet.flow_mol_phase_comp[0.0, 'Liq', i] for i in blk.liquid)
            else:
                return 0.0
        elif j == 'purge':
            return blk.S102.purge.mole_frac_comp[0,i]*blk.component_LHV_values[i]

    m.fs.LHV_per_component_per_stream = Expression(m.fs.Streams,m.fs.vapor,rule=return_LHV_per_component_per_stream) # MJ/mol

    def return_LHV_per_unit_stream(blk,j):
        return sum(blk.LHV_per_component_per_stream[j,i] for i in blk.vapor)

    m.fs.LHV_per_unit_stream = Expression(m.fs.Streams,rule=return_LHV_per_unit_stream) # MJ/mol

    def return_LHV_per_stream(blk,j):
        if j == 'shale':
            return blk.LHV_per_unit_stream[j]*blk.shale_flowrate*3600/1000
        elif j == 'NG':
            return blk.LHV_per_unit_stream[j]*blk.NG_flowrate*3600/1000
        elif j == 'feed':
            return blk.LHV_per_unit_stream[j]*blk.M101.feed.flow_mol[0]*3600/1000
        elif j == 'H2_vent':
            return blk.LHV_per_unit_stream[j]*blk.S101.H2_purge.flow_mol[0]*3600/1000
        elif j == 'liq_outlet':
            return blk.LHV_per_unit_stream[j]*sum(blk.F102.liq_outlet.flow_mol_phase_comp[0.0, 'Liq', i] for i in blk.liquid)*3600/1000
        elif j == 'purge':
            return blk.LHV_per_unit_stream[j]*blk.S102.purge.flow_mol[0]*3600/1000

    m.fs.LHV_per_stream = Expression(m.fs.Streams,rule=return_LHV_per_stream) # GJ/h

    def return_MW_per_component_per_stream(blk,j,i):
        if j == 'shale':
            return blk.shale_feed_comp[i]*blk.component_MW_values[i]
        elif j == 'NG':
            return blk.NG_feed_comp[i]*blk.component_MW_values[i]
        elif j == 'feed':
            return blk.M101.feed.mole_frac_comp[0,i]*blk.component_MW_values[i]
        elif j == 'H2_vent':
            return blk.S101.H2_purge.mole_frac_comp[0,i]*blk.component_MW_values[i]
        elif j == 'liq_outlet':
            if i in m.fs.liquid:
                return blk.F102.liq_outlet.flow_mol_phase_comp[0.0, 'Liq', i]*blk.component_MW_values[i]
            else:
                return 0.0
        elif j == 'purge':
            return blk.S102.purge.mole_frac_comp[0,i]*blk.component_MW_values[i]

    m.fs.MW_upstream = Expression(m.fs.Streams,m.fs.vapor,rule=return_MW_per_component_per_stream) # g/mol

    def return_MW_totals(blk,j):
        return sum(blk.MW_upstream[j,i] for i in m.fs.vapor)

    m.fs.MW_upstream_total = Expression(m.fs.Streams,rule=return_MW_totals) #g/mol


def calculate_stream_energies(m):
    """
    Function to calculate energy content (LHV) of each inlet and outlet stream of the process
    """
    ### LHV from IDAES flow streams ## include H2
    m.fs.feed_energy = Expression(expr = m.fs.LHV_per_stream['feed']) # GJ NGL/h
    m.fs.NG_energy = Expression(expr = m.fs.LHV_per_stream['NG']) # GJ NG/h
    m.fs.fuel_energy = Expression(expr = m.fs.LHV_per_stream['liq_outlet']) # GJ Liquid HC/h
    m.fs.h2_energy = Expression(expr = m.fs.LHV_per_stream['H2_vent']) # GJ H2-gas/h
    m.fs.gas_energy = Expression(expr = m.fs.LHV_per_stream['purge']) # GJ purge/h

def calculate_emissions(m,case_name,emissions_factor_file, include_HI=True):
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
    upstream_emissions_factor_df = pd.read_csv(emissions_factor_file)
    ## Upstream production
    # Emissions factor to account for upstream methane and NGL processing emissions
    m.fs.upstream_emission_factor = Param(initialize=upstream_emissions_factor_df[case_name].iloc[0], mutable=True) ## kg CO2e/GJ of NGL (process feed) factor from GREET 2022

    m.fs.em_upstream = Expression(expr = m.fs.upstream_emission_factor*m.fs.feed_energy) ## kg co2/h


    ### Compressors
    m.fs.total_electricity = Expression(expr = (m.fs.E101.work_mechanical[0] + m.fs.C101.work_mechanical[0] + 
                                                m.fs.C102.work_mechanical[0])*0.001*0.001*0.001*3600) ## GJ/h

    ### Emission factors
    ### Electricity
    m.fs.heating_factor = Param(initialize = 82.8862946, mutable=True) ## kg CO2e/GJ
    m.fs.co2_electricity_factor = Param(initialize = 134.7489276, mutable=True) ## kg CO2e/GJ

    m.fs.purge_heat_efficiency = Param(initialize=0.6) # efficiency of burning purge stream for heat
    m.fs.purge_heat = Expression(expr = m.fs.purge_heat_efficiency*m.fs.gas_energy) # GJ/h

    if include_HI:
        m.fs.co2_em_heating = Expression(expr = m.fs.heating_factor*m.fs.Qs) ## kg CO2e/h
    else:
        m.fs.co2_em_heating = Expression(expr = m.fs.heating_factor*(m.fs.H101.heat_duty[0] + 
                                                m.fs.H103.heat_duty[0] + m.fs.H107.heat_duty[0] + 
                                                m.fs.R101.heat_duty[0])*0.001*0.001*0.001*3600) ## kg CO2e/h

    m.fs.co2_em_electricity = Expression(expr = m.fs.co2_electricity_factor*m.fs.total_electricity) ## kg CO2e/h

    ### Process emissions
    m.fs.em_process = Expression(expr = m.fs.co2_em_heating + m.fs.co2_em_electricity) ## kg CO2e/h

    m.fs.upstream_emissions = Var(initialize = 1.0, bounds = (1e-6, None)) # upstream, kg CO2e/GJ fuel
    m.fs.downstream_emissions = Var(initialize = 1.0, bounds = (1e-6, None)) # downstream, kg CO2e/GJ fuel

    m.fs.emissions_total_constraint = Constraint(expr = m.fs.upstream_emissions * (m.fs.fuel_energy+m.fs.gas_energy+m.fs.h2_energy) == m.fs.em_upstream)
    m.fs.emissions_process_constraint = Constraint(expr = m.fs.downstream_emissions * (m.fs.fuel_energy+m.fs.gas_energy+m.fs.h2_energy) == m.fs.em_process)

def create_ghg_objective(m):
    """
    Function to calculate GHG objective expression
    """
    ## Aggregated emissions for Liquid Fuels, kg CO2e/GJ fuel
    m.fs.GHG_obj_2 = Expression(expr = m.fs.upstream_emissions + m.fs.downstream_emissions)


def delete_region_specific_components(m):
    """
    Function to delete region-specific emisisons data in order to update emissions 
    calculations from a single initialization file (single region)
    """
    m.fs.del_component(m.fs.component_LHV_values)
    m.fs.del_component(m.fs.component_MW_values)
    m.fs.del_component(m.fs.shale_feed_comp)
    m.fs.del_component(m.fs.shale_flowrate)
    m.fs.del_component(m.fs.shale_flow_by_comp)
    m.fs.del_component(m.fs.NG_feed_comp)
    m.fs.del_component(m.fs.NG_flowrate)
    m.fs.del_component(m.fs.NG_flow_by_comp)
    m.fs.del_component(m.fs.Streams)
    m.fs.del_component(m.fs.LHV_per_component_per_stream)
    m.fs.del_component(m.fs.LHV_per_component_per_stream_index)
    m.fs.del_component(m.fs.LHV_per_unit_stream)
    m.fs.del_component(m.fs.LHV_per_stream)
    m.fs.del_component(m.fs.MW_upstream)
    m.fs.del_component(m.fs.MW_upstream_index)
    m.fs.del_component(m.fs.MW_upstream_total)