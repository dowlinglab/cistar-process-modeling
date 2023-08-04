import matplotlib.pyplot as plt
import numpy as np
from pyomo.environ import value
from utility_minimization_1d import gen_curves
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

colors_dict = {'hydrogen':'lightgray',
                'methane':'gray',
                'ethane':'rosybrown',
                'propane':'yellowgreen',
                'nbutane':'tan',
                'ibutane':'orange',
                'pentane':'lime',
                'hexane':'paleturquoise',
                'heptane':'lightsteelblue',
                'octane':'violet',
                'ethylene':'maroon',
                'propene':'darkolivegreen',
                'butene':'darkgoldenrod',
                'pentene':'darkgreen',
                'hexene':'lightseagreen',
                'heptene':'royalblue',
                'octene':'darkviolet',
                'nonene':'crimson'
                }

def generate_and_save_composite_curves(CD,model_code,C_tax_rate,region,optimal_solution):
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
    # Transform information from CD class into arrays for plotting function

    # Transforms cooling equipment data
    CTin = CD.Cooling_Tin
    CTout = CD.Cooling_Tout
    CQ = CD.Cooling_Q
    CTin, CTout, CQ = np.array((CTin, CTout, CQ))

    # Transforms heating equipment data
    HTin = CD.Heating_Tin
    HTout = CD.Heating_Tout
    HQ = CD.Heating_Q
    HTin, HTout, HQ = np.array((HTin, HTout, HQ))

    Thot, Qhot = gen_curves(CTin, CTout, CQ)
    Tcold, Qcold = gen_curves(HTin, HTout, -HQ)
    # Brings Qw scalar into scope out of the class
    Qw = CD.Qw
    

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.rcParams.update({'font.size': 20})
    # Plot values for Hot streams
    # Negative sign corrects for Q < 0 for hot streams
    p_hot = ax.plot([-Qhot_i*1000/3600 for Qhot_i in Qhot], Thot, color="r", label="Hot Streams", linewidth=3)

    # Need to shift curve to account for the cooling utility
    Qcold = Qcold + sum(HQ)
    # Plot values for cold streams
    plt.plot(
        [(Qcold_i+value(Qw))*1000/3600 for Qcold_i in Qcold],
        Tcold,
        color="b",
        label="Cold Streams", linewidth=3
    )
    ax.set_xlabel(
        "Cumulative process-wide \n heat exchange [MW]", fontsize=24, weight='bold'
    )
    ax.set_ylabel("Temperature [" + str(CD.T_unit) + "]", fontsize=24, weight='bold')
    
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    # plt.title("Composite Curves for minimum utilities", size=20)
    plt.legend(loc='best', fontsize=20)
    plt.grid()
    if optimal_solution:
        plt.savefig('./plots/'+'composite_curve_M{}_C-tax_{}_region_{}_optimal.pdf'.format(model_code,C_tax_rate,region),bbox_inches='tight',dpi=200)
    else:
        plt.savefig('./plots/'+'composite_curve_M{}_C-tax_{}_region_{}_initialized.pdf'.format(model_code,C_tax_rate,region),bbox_inches='tight',dpi=200)
    plt.show()
    return


def plot_outlet_flowrate_by_ROK_model_horizontal_stacked_bars(liquid_hc_component_flow_dict,ROK_model_list,
                                                            xlabel_string, ylabel_string, file_name_string,
                                                            ):
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
    # ROK_models_list = ['M2','M3','M4','M5']

    assert (len(liquid_hc_component_flow_dict.keys()) == len(ROK_model_list))

    width = 0.5

    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})
    left = np.zeros(len(ROK_model_list))

    for c, comp in liquid_hc_component_flow_dict['M2'].items():
    #     colors_dict[c] = 
        p = ax.barh(ROK_model_list, [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()], width,label=c, left=left,color=colors_dict[c])
        
        left += [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()] # add flowrates already plotted for stacking

    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')
    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')

    ax.tick_params(axis='y', labelsize= 16)
    ax.tick_params(axis='x', labelsize= 16)#, rotation=45)
    ax.bar_label(p, label_type='edge',fontsize=16,weight='bold', fmt='%.1f')

    ax.set_xlim(0,210)

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.2),ncol=4,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()

def plot_outlet_flowrate_by_region_horizontal_stacked_bars(liquid_hc_component_flow_dict,region_list,
                                                            xlabel_string, ylabel_string, file_name_string,
                                                            ):
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
    assert (len(liquid_hc_component_flow_dict.keys()) == len(region_list))

    width = 0.5

    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})
    left = np.zeros(len(region_list))

    for c, comp in liquid_hc_component_flow_dict['EF-1'].items():
    #     colors_dict[c] = 
        p = ax.barh(region_list, [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()], width,label=c, left=left,color=colors_dict[c])
        
        left += [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()] # add flowrates already plotted for stacking

    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')
    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')

    ax.tick_params(axis='y', labelsize= 16)
    ax.tick_params(axis='x', labelsize= 16)#, rotation=45)
    ax.bar_label(p, label_type='edge',fontsize=16,weight='bold', fmt='%.1f')

    ax.set_xlim(0,300)

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.2),ncol=4,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()


def plot_C4_C9_olefin_outlet_flowrate_by_ROK_model_vertical_bars(liquid_C4_C9_percent_dict,ROK_models_list,
                                                                xlabel_string, ylabel_string, file_name_string
                                                                ):
    """
    Function to create horizontal stacked plot of liquid hydrocarbon composition in mol/s
    w.r.t. ROK model
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
    # ROK_models_list = ['M2','M3','M4','M5']
    xlocs = np.arange(len(ROK_models_list))
    fig, ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})
    # for model_code,m in model_data_dict.items():
    bar = ax.bar(ROK_models_list,[v for k,v in liquid_C4_C9_percent_dict.items()], align='center')

    ax.set_ylabel(ylabel_string,fontsize=24, weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24, weight='bold')
    ax.set_xticks(xlocs)
    ax.set_ylim(0,100)

    ax.tick_params(axis='y', labelsize= 20)
    ax.bar_label(bar, label_type='edge',weight='bold',fontsize=16, fmt='%.2f')

    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()


def plot_outlet_flow_rate_by_component_vs_ROK_model_lines(liquid_hc_component_flow_dict, ROK_model_list,
                                                            xlabel_string, ylabel_string,file_name_string
                                                            ):
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
    # ROK_models_list = ['M2','M3','M4','M5']
    assert (len(liquid_hc_component_flow_dict.keys()) == len(ROK_model_list))
    width = 0.5

    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})

    for c, comp in liquid_hc_component_flow_dict['M2'].items():
    #     colors_dict[c] = 
        p = ax.plot(ROK_model_list, [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()], width,label=c, color=colors_dict[c])
        

    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')

    ax.tick_params(axis='y', labelsize= 16)
    ax.tick_params(axis='x', labelsize= 16)#, rotation=45)

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.2),ncol=4,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()

def plot_outlet_flow_rate_by_component_vs_region_lines(liquid_hc_component_flow_dict, region_list,
                                                            xlabel_string, ylabel_string,file_name_string
                                                            ):
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
    # ROK_models_list = ['M2','M3','M4','M5']
    assert (len(liquid_hc_component_flow_dict.keys()) == len(region_list))
    width = 0.5

    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})

    for c, comp in liquid_hc_component_flow_dict['EF-1'].items():
    #     colors_dict[c] = 
        p = ax.plot(region_list, [comp_dict[c] for r,comp_dict in liquid_hc_component_flow_dict.items()], width,label=c, color=colors_dict[c])
        

    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')

    ax.tick_params(axis='y', labelsize= 16)
    ax.tick_params(axis='x', labelsize= 16)#, rotation=45)

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.2),ncol=4,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()


def plot_emissions_by_oligo_model_vertical_stacked_bars(upstream_emissions_dict,downstream_emissions_dict,
                                                    oligo_model_list,xlabel_string, ylabel_string,
                                                    file_name_string
                                                    ):
    """
    Function to create stacked bar plot of upstream and downstream emissions in g CO2e/MJ fuel
    w.r.t. model instance
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
    assert (len(upstream_emissions_dict.keys())+1 == len(oligo_model_list) and len(downstream_emissions_dict.keys())+1 == len(oligo_model_list)) # +1 to add-in literature data inline
    width = 0.5

    xlocs = np.arange(len(oligo_model_list))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})

    p1 = ax.bar(xlocs, [9.26]+[v for k,v in upstream_emissions_dict.items()], width=0.75,label='Upstream')
    p2 = ax.bar(xlocs, [14.2]+[v for k,v in downstream_emissions_dict.items()], width=0.75,bottom = [9.26]+[v for k,v in upstream_emissions_dict.items()],label='Process')

    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')
    ax.set_xticks(xlocs)
    ax.set_xticklabels(oligo_model_list ,fontsize=16,rotation=45)
    ax.tick_params(axis='y', labelsize= 16)

    ax.bar_label(p1, label_type='center',fontsize=14,weight='bold',rotation=90, fmt='%.2f')
    ax.bar_label(p2, label_type='center',fontsize=14,weight='bold',rotation=90, fmt='%.2f')

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.34),ncol=2,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()

def plot_emissions_by_region_vertical_stacked_bars(upstream_emissions_dict,downstream_emissions_dict,
                                                    region_list,xlabel_string, ylabel_string,
                                                    file_name_string
                                                    ):
    """
    Function to create stacked bar plot of upstream and downstream emissions in g CO2e/MJ fuel
    w.r.t. shale region
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
    assert (len(upstream_emissions_dict.keys()) == len(region_list) and len(downstream_emissions_dict.keys()) == len(region_list))
    width = 0.5

    xlocs = np.arange(len(region_list))
    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})

    p1 = ax.bar(xlocs, [v for k,v in upstream_emissions_dict.items()], width=0.75,label='Upstream')
    p2 = ax.bar(xlocs, [v for k,v in downstream_emissions_dict.items()], width=0.75,bottom = [v for k,v in upstream_emissions_dict.items()],label='Process')

    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')
    ax.set_xticks(xlocs)
    ax.set_xticklabels(region_list ,fontsize=16,rotation=45)
    ax.tick_params(axis='y', labelsize= 16)

    ax.bar_label(p1, label_type='center',fontsize=14,weight='bold',rotation=90, fmt='%.2f')
    ax.bar_label(p2, label_type='center',fontsize=14,weight='bold',rotation=90, fmt='%.2f')

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.34),ncol=2,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()


def plot_LHV_by_region_stacked_vertical_bars(LHV_dict,regions_list,xlabel_string, ylabel_string,
                                            file_name_string):
    """
    Function to create stacked bar plot of LHV by component vs region, MJ/kmol fuel
    w.r.t. model instance
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
    width = 0.5

    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})
    bottom = np.zeros(len(regions_list))

    for c, comp in LHV_dict['EF-1'].items():
        p = ax.bar(regions_list, [comp_dict[c] for r,comp_dict in LHV_dict.items()], width,label=c, bottom=bottom,color=colors_dict[c])
        bottom += [comp_dict[c] for r,comp_dict in LHV_dict.items()]

    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')

    ax.tick_params(axis='y', labelsize= 16)
    ax.tick_params(axis='x', labelsize= 16, rotation=45)
    ax.bar_label(p, label_type='edge',fontsize=16,weight='bold', fmt='%.f',rotation=90)
    ax.set_ylim(0,950)

    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.3),ncol=4,fontsize=20)
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
    plt.show()

def plot_MSP_vs_c_tax_rate(optimal_MSP_dict, c_tax_rates, xlabel_string,
                           ylabel_string, file_name_string):
    """
    Function to create stacked bar plot of LHV by component vs region, MJ/kmol fuel
    w.r.t. model instance
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
    fig,ax = plt.subplots(figsize=(8,6))
    plt.rcParams.update({'font.size': 20})
    p1 = ax.semilogx([c_tax*1000 for c_tax in c_tax_rates],[v/1000 for k,v in optimal_MSP_dict.items()], marker='o',markersize=10, linewidth=3) # convert MSP to $/MJ, tax to $/tonne CO2e
    ax.set_xlabel(xlabel_string,fontsize=24,weight='bold')
    ax.set_ylabel(ylabel_string,fontsize=24,weight='bold')
    ax.grid(which='major')

    from matplotlib import ticker
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax.yaxis.set_major_formatter(formatter)

    ax.tick_params(axis='y', labelsize= 16, length=6)
    ax.tick_params(axis='x', labelsize= 16, length=6, which='both')
    plt.savefig('./plots/'+file_name_string,bbox_inches='tight',dpi=200)
