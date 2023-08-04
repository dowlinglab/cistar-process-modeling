# -*- coding: UTF-8 -*-
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
General purpose heat integration block for IDAES models
"""

# Import Python libraries
import logging

# Import Plotting and Numpy libraries for composite curves
import matplotlib.pyplot as plt
import numpy as np

# Import Pyomo libraries
from pyomo.environ import Constraint, sqrt, Expression, value, Var

from pyomo.environ import units as pyunits

import pandas as pd

__author__ = "Alejandro Garciadiego, Kanishka Ghosh, Alexander Dowling"

# Temperature Epsilon to add or subsract 1 degree to avoid division over
# 0 in equipment not changing temperature
EpsT = 1


def min_utility(blk, heating, cooling, DTmin, eps=1e-6, DG_units=pyunits.Gjoule/pyunits.hour):
    """
    Function for Duran-Grossmann for minimization of utilities
    This function adds variables and constraints to the flowsheet to
    minimize utilities

    Args:
        blk: flowsheet
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
        DTmin: HRAT (Heat Recovery Approximation Temperature)
        eps: Epsilon for smoothing operation

    Returns:
        Constraint and variable object representing the calculation
        for hot utility and cold utility
    """

    # Generate lists out of strings from Args
    exchanger_list = heating + cooling
    pinch_streams = exchanger_list

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    Q = {}

    for i in exchangerdict:
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            Q[i] = pyunits.convert(sum(pinch_streamsdict[i].heat_duty[0.0,j] for j in pinch_streamsdict[i].control_volume.length_domain if j != 0) * 0.05 * pyunits.m, to_units=DG_units)
        else:
            Q[i] = pyunits.convert(pinch_streamsdict[i].heat_duty[0], to_units=DG_units)

    # Run function to fill exchanger data for pinch calculations for
    # initialization of variables
    exchangeData = heat_data(blk, heating, cooling, DG_units)

    # Call pinch calculation for initialization of variables in PD class
    PD = pinch_calc(heating, cooling, exchangeData, DTmin, eps)

    # Define dictionary for addition of DTmin to heating equipment
    dT = {}

    for i, el in enumerate(cooling):
        dT[str(el)] = 0

    for i, el in enumerate(heating):
        dT[str(el)] = DTmin

    def T_in(blk, i):
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            return pinch_streamsdict[i].control_volume.properties[0,0.0].temperature
        else:
            return pinch_streamsdict[i].control_volume.properties_in[0].temperature
        
    blk.Tin = Expression(
        pinch_streamsdict.keys(), rule=T_in, doc="Inlet temperature in exchangers"
    )

    def T_out(blk, i):
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            return pinch_streamsdict[i].control_volume.properties[0,1.0].temperature
        else:
            return pinch_streamsdict[i].control_volume.properties_out[0].temperature

    blk.Tout = Expression(
        pinch_streamsdict.keys(), rule=T_out, doc="Outlet temperature in exchangers"
    )

    blk.Theta = Var(pinch_streamsdict.keys(),initialize = 1.0,bounds=(1e-8,None))
    # Expression for cp of equimpent with heat exchange
    def Theta_(blk, i):
        if i in heatingdict.keys():
            return blk.Theta[i] * (
                0.5
                * (
                    (blk.Tout[i] - blk.Tin[i] + EpsT)
                    + sqrt((blk.Tout[i] - blk.Tin[i] + EpsT) ** 2 + eps)
                )
            )== Q[i]
        else:
            return blk.Theta[i] *(
                -0.5
                * (
                    (blk.Tin[i] - blk.Tout[i] + EpsT)
                    + sqrt((blk.Tin[i] - blk.Tout[i] + EpsT) ** 2 + eps)
                )
            )== Q[i]

    blk.Theta_constraint = Constraint(
        pinch_streamsdict.keys(), rule=Theta_, doc="FCp in exchangers"
    )

    # Define expression for pinch candidate temperature
    def T_(blk, i):
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            return pinch_streamsdict[i].control_volume.properties[0,0.0].temperature + dT[i]
        else:
            return pinch_streamsdict[i].control_volume.properties_in[0].temperature + dT[i]

    blk.T_ = Expression(
        pinch_streamsdict.keys(), rule=T_, doc="Pinch candidate temperature"
    )

    # Define variable for heat content above the pinch point
    blk.QAh = Var(
        pinch_streamsdict.keys(),
        initialize=PD.initQAh,
        bounds=(1e-8, None),
        doc="Heat content above pinch",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to calculate the varaible QAh
    def rule_heat_above_pinch(blk, p):
        return blk.QAh[p] == sum(
            blk.Theta[i]
            * (
                0.5
                * ((blk.Tin[i] - blk.T_[p] + EpsT) + sqrt((blk.Tin[i] - blk.T_[p] + EpsT) ** 2 + eps))
                - 0.5
                * (
                    (blk.Tout[i] - blk.T_[p])+ sqrt((blk.Tout[i] - blk.T_[p]) ** 2 + eps))
            )
            for i in coolingdict.keys()
        )

    blk.heat_above_pinch = Constraint(
        pinch_streamsdict.keys(), rule=rule_heat_above_pinch
    )

    # Define variable for heat content below the pinch point
    blk.QAc = Var(
        pinch_streamsdict.keys(),
        initialize=PD.initQAc,
        bounds=(1e-8, None),
        doc="Heat content below pinch",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to calculate the varaible QAc
    def rule_heat_below_pinch(blk, p):
        return blk.QAc[p] == sum(
            blk.Theta[i]
            * (
                0.5
                * (
                    (blk.Tout[i] - blk.T_[p] + DTmin + EpsT)
                    + sqrt((blk.Tout[i] - blk.T_[p] + DTmin + EpsT) ** 2 + eps)
                )
                - 0.5
                * (
                    (blk.Tin[i] - blk.T_[p] + DTmin)
                    + sqrt((blk.Tin[i] - blk.T_[p] + DTmin) ** 2 + eps)
                )
            )
            for i in heatingdict.keys()
        )

    blk.heat_below_pinch = Constraint(
        pinch_streamsdict.keys(), rule=rule_heat_below_pinch
    )

    # Define variable for Heat of hot utility
    blk.Qs = Var(
        initialize=PD.initQs,
        bounds=(1e-8, None),
        doc="Heating utilities",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to solve for Qs
    # Where 1E-6 is added to both sides of the constraint as a scaling factor
    def rule_heating_utility(blk, p):
        return blk.Qs >= (blk.QAc[p] - blk.QAh[p])

    blk.heating_utility = Constraint(
        pinch_streamsdict.keys(), rule=rule_heating_utility
    )

    # Define variable for Heat of cold utility
    blk.Qw = Var(
        initialize=PD.initQw,
        bounds=(1e-8, None),
        doc="Cooling utilities",
        units=pyunits.get_units(Q[str(pinch_streams[0])]),
    )

    # Define a constraint to solve for Qw
    # Where 1E-6 is added to both sides of the constraint as a scaling factor
    def rule_cooling_utility(blk):
        return blk.Qw == -sum(Q[i] for i in exchangerdict.keys()) + blk.Qs

    blk.cooling_utility = Constraint(rule=rule_cooling_utility)


def heat_data(blk, heating, cooling, DG_units=pyunits.Gjoule/pyunits.hour):
    """
    Function for generating necesary data for Duran-Grossmann initialization
    Allows the calculation of reactors
    Args:
        blk: flowsheet
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
    Returns:
        Dictionary for the equipment exchanging heat containing:
            Tin
            Tout
            FCp
            Q
    """
    # Generate lists out of strings from Args
    exchanger_list = heating + cooling

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate dictionaries for information for heat data
    Q = {}
    T_in = {}
    T_out = {}
    FCp_ = {}

    # Defining inlet temperature for inlet from equipment control volume
    for i in heatingdict.keys():
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            T_in[i] = value(
                pinch_streamsdict[i].control_volume.properties[0,0.0].temperature
            )
            T_out[i] = (
                value(pinch_streamsdict[i].control_volume.properties[0,1.0].temperature)
                + EpsT
            )
        else:
            T_in[i] = value(
                pinch_streamsdict[i].control_volume.properties_in[0].temperature
            )
            T_out[i] = (
                value(pinch_streamsdict[i].control_volume.properties_out[0].temperature)
                + EpsT
            )

    # Defining inlet temperature for utlet from equipment control volume
    for i in coolingdict.keys():
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            T_in[i] = value(
                pinch_streamsdict[i].control_volume.properties[0,0.0].temperature
                + EpsT
            )
            T_out[i] = value(
                pinch_streamsdict[i].control_volume.properties[0,1.0].temperature
            )
                
            
        else:
            T_in[i] = (
                value(pinch_streamsdict[i].control_volume.properties_in[0].temperature)
                + EpsT
            )
            T_out[i] = value(
                pinch_streamsdict[i].control_volume.properties_out[0].temperature
            )

    # Calculating FCp out of heat in control volume
    # Obtaines Q from equipment's control volume heat
    for i in pinch_streamsdict.keys():
        if hasattr(pinch_streamsdict[i].control_volume, 'length_domain'):
            Q[i] = value(
                pyunits.convert(sum(pinch_streamsdict[i].heat_duty[0.0,j] for j in pinch_streamsdict[i].control_volume.length_domain if j != 0) * 0.05 * pyunits.m, to_units=DG_units)
            )
        else:
            Q[i] = value(
                pyunits.convert(pinch_streamsdict[i].heat_duty[0], to_units=DG_units)
            )
        FCp_[i] = Q[i] / (T_out[i] - T_in[i])

    # Generate a large dictioary containing all the data obtained
    # from the equipment
    exchangeData = {}
    for i in pinch_streamsdict.keys():
        exchangeData[i] = {
            "T_in": T_in[i],
            "T_out": T_out[i],
            "FCp_": FCp_[i],
            "Q": Q[i],
        }

    return exchangeData


def pinch_calc(heating, cooling, exchangeData, DTmin, eps):
    """
    Function for calculating heat data for Duran-Grossmann initialization

    Args:
        heating: Equipment that requires Heat
        cooling: Equipment from wich heat is being removed
        exchangeData: Dictionary containing Tin, Tout, FCp and Q
        for each equipment in lists
        DTmin: HRAT (Heat Recovery Approximation Temperature)
        eps: Epsilon for smoothing operation
    Returns:
        PD (Dictionary containint initialized values for QAh, QAc, Qs and Qw)
    """
    # Generate lists out of strings from Args
    exchanger_list = heating + cooling

    # Generate dictionaries out of strings to convert to pyomo objects
    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate dictionaries to contain initialized data
    T_ = {}
    initQAh = {}
    initQAc = {}
    b = []

    # Define dictionary for addition of DTmin to heating equipment
    dT = {}
    for i, el in enumerate(cooling):
        dT[str(el)] = 0

    for i, el in enumerate(heating):
        dT[str(el)] = DTmin

    # Calculate pinch temperature candidate
    # Calculate QAh and QAc
    for i in pinch_streamsdict.keys():
        T_[i] = exchangeData[i]["T_in"] + dT[i]
        initQAh[i] = sum(
            exchangeData[j]["FCp_"]
            * (
                0.5
                * (
                    (exchangeData[j]["T_in"] - T_[i])
                    + sqrt((exchangeData[j]["T_in"] - T_[i]) ** 2 + eps)
                )
                - 0.5
                * (
                    (exchangeData[j]["T_out"] - T_[i])
                    + sqrt((exchangeData[j]["T_out"] - T_[i]) ** 2 + eps)
                )
            )
            for j in coolingdict.keys()
        )
        initQAc[i] = sum(
            exchangeData[j]["FCp_"]
            * (
                0.5
                * (
                    (exchangeData[j]["T_out"] - T_[i] + DTmin)
                    + sqrt((exchangeData[j]["T_out"] - T_[i] + DTmin) ** 2 + eps)
                )
                - 0.5
                * (
                    (exchangeData[j]["T_in"] - T_[i] + DTmin)
                    + sqrt((exchangeData[j]["T_in"] - T_[i] + DTmin) ** 2 + eps)
                )
            )
            for j in heatingdict.keys()
        )

    # Generate array with all possible QS
    for i in exchangerdict.keys():
        b.append(initQAc[i] - initQAh[i])

    # Define largest value of QS
    c = max([max(b), 0.0])
    initQs = c

    # Calculate Qw from largest value of Qs
    initQw = -sum(value(exchangeData[i]["Q"]) for i in exchangerdict.keys()) + initQs
    initQw = max([initQw, 0.0])

    # Fill Class with all the data to initialize Duran-Grossmann variables
    PD = PinchDataClass(initQs, initQw)
    PD.initQAh = initQAh
    PD.initQAc = initQAc
    PD.T_ = T_

    return PD


def generate_curves(CD):
    """
    Function for plotting composite curves

    Args:
        CD: Class containing curve data
    Returns:
        Composite curves
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
    

    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)

    # Plot values for Hot streams
    # Negative sign corrects for Q < 0 for hot streams
    plt.plot(-Qhot, Thot, color="r", label="Hot Streams / Streams that are Cooled")

    # Need to shift curve to account for the cooling utility
    Qcold = Qcold + sum(HQ)
    # Plot values for cold streams
    plt.plot(
        (Qcold + value(Qw)),
        Tcold,
        color="b",
        label="Cold Streams / Streams that are Heated",
    )
    plt.xlabel(
        "Cumulative Process-Wide Heat Exchange [" + str(CD.Q_unit) + "]", fontsize=18
    )
    plt.ylabel("Temperature [" + str(CD.T_unit) + "]", fontsize=18)
    plt.title("Composite Curves for minimum utilities", size=20)
    plt.legend(loc='upper center',bbox_to_anchor=(0.5, -0.3), fontsize=18)
    plt.grid()
    plt.show()
    return

def heat_ex_data(blk, heating, cooling):
    """
    Function for turning IDAES heat exchanging equipment into a class for
    use in plotting

    Args:
        blk: Flowsheet
        heating: Equipment that heats streams
        cooling: Equipment that cools streams
    Returns:
        CD: Class with heating and cooling equipment as arrays
    """
    # Generate dictionaries out of strings to convert to pyomo objects
    exchanger_list = heating + cooling
    pinch_streams = exchanger_list

    pinch_streamsdict = {}
    coolingdict = {}
    heatingdict = {}
    exchangerdict = {}

    T_units = pyunits.get_units(
        pinch_streams[0].control_volume.properties_in[0].temperature
    )
    DG_units = pyunits.get_units(blk.Qw)

    for i, el in enumerate(exchanger_list):
        pinch_streamsdict[str(el)] = el

    for i, el in enumerate(cooling):
        coolingdict[str(el)] = el

    for i, el in enumerate(heating):
        heatingdict[str(el)] = el

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # Generate zero arrays from length of the cooling list
    CHX = len(coolingdict)
    CTin = np.zeros(CHX)
    CTout = np.zeros(CHX)
    CQ = np.zeros(CHX)

    # Convert Pyomo model values into arrays
    j = 0
    for i in coolingdict.keys():
        if hasattr(coolingdict[i].control_volume, 'length_domain'):
            CTin[j] = value(coolingdict[i].control_volume.properties[0,0.0].temperature) + 1
            CTout[j] = value(coolingdict[i].control_volume.properties[0,1.0].temperature)
            CQ[j] = value(
                pyunits.convert(sum(coolingdict[i].control_volume.heat[0,k] for k in coolingdict[i].control_volume.length_domain if k != 0) * 0.05 * pyunits.m, to_units=DG_units)
            )
        else:
            CTin[j] = value(coolingdict[i].control_volume.properties_in[0].temperature) + 1
            CTout[j] = value(coolingdict[i].control_volume.properties_out[0].temperature)
            CQ[j] = value(
                pyunits.convert(coolingdict[i].control_volume.heat[0], to_units=DG_units)
            )
        j += 1

    # Generate zero arrays from length of the heating list
    HHX = len(heatingdict)
    HTin = np.zeros(HHX)
    HTout = np.zeros(HHX)
    HQ = np.zeros(HHX)

    # Convert Pyomo model values into arrays
    j = 0
    for i in heatingdict.keys():
        if hasattr(heatingdict[i].control_volume, 'length_domain'):
            HTin[j] = value(heatingdict[i].control_volume.properties[0,0.0].temperature)
            HTout[j] = (
                value(heatingdict[i].control_volume.properties[0,1.0].temperature) + 1
            )
            HQ[j] = value(
                pyunits.convert(sum(heatingdict[i].control_volume.heat[0,k] for k in heatingdict[i].control_volume.length_domain if k != 0) * 0.05 * pyunits.m, to_units=DG_units)
            )
        else:
            HTin[j] = value(heatingdict[i].control_volume.properties_in[0].temperature)
            HTout[j] = (
                value(heatingdict[i].control_volume.properties_out[0].temperature) + 1
            )
            HQ[j] = value(
                pyunits.convert(heatingdict[i].control_volume.heat[0], to_units=DG_units)
            )
        j += 1

        # Fill class with values and arrays
    Qw = value(blk.Qw)
    CD = CurveData(Qw, T_units, DG_units)
    CD.Cooling_Tin = CTin
    CD.Cooling_Tout = CTout
    CD.Cooling_Q = CQ
    CD.Heating_Tin = HTin
    CD.Heating_Tout = HTout
    CD.Heating_Q = HQ
    return CD


def gen_curves(Tin, Tout, Q):
    """
    Function to do add cumulative heat arrays
    Args:
        Tin: Inlet temperatures of cooling/heating equipment
        Toout: Outlet temperatures of cooling/heating equipment
        Q: Heat of cooling/heating equipment
    Returns:
        Tstar: Cumulative Temperature array
        NQstar: Cumulative heat array
    """
    # Q < 0 for hot streams = heat removing = cooling units
    # Q > 0 for cold streams = heat added = heating units

    # Ignoring edge cases for phase changes
    # Shaping Temperature array
    # Calls unique function to avoid  repeating
    Tstart = unique(np.vstack([Tin, Tout]).reshape((1, -1)))
    Tstar = np.sort(Tstart[0])

    # Generate vector of Tin size and Tunique
    nTin = len(Tin)
    nTstar = len(Tstar)
    Qmat = np.zeros((nTin, nTstar))

    # Generate cumulative arays for temperature and heat
    for i in range(nTin):
        for j in range(nTstar):
            Qmat[i, j] = linear_interpolation([Tin[i], Tout[i]], [0, Q[i]], Tstar[j])
    Qstar = sum(Qmat, 0)
    NQstar = Qstar.reshape(nTstar)
    return Tstar, NQstar


def linear_interpolation(x, y, t):
    """
    Function to do Linear interpolation with nearest neighbor extrapolation
    Args:
        x: Inlet and Outlet Temperature values
        y: 0 and Heat value
        t: Unique inlet and outlet temperatures
    Returns:
        Qstar: New array of Q Values
    """
    # Set the upper or lower value
    if x[0] < x[1]:
        lo = 0
        hi = 1
    else:
        lo = 1
        hi = 0
    # do Linear interpolation with nearest neighbor extrapolation
    alpha = (x[hi] - t) / (x[hi] - x[lo])
    alpha = max([0, alpha])
    alpha = min([alpha, 1])
    Qsta = alpha * (y[hi] - y[lo]) + y[lo]
    return Qsta


def print_HX_results(blk, exchanger_list):
    """
    Function to print results of Heat Exchangers
    Args:
        blk: flowsheet of the equipment
        exchanger_list: List of equipment to print data
    Returns:
        Printed List
    """
    # Initialize dictionary and fill with exchanger list
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # initialize null dictionaries for data to be printed
    Tin_ = {}
    Tout_ = {}
    f_ = {}
    Q_ = {}
    DG_units_ = {}

    # Loop over heat exchangers
    for i in exchangerdict.keys():
        if hasattr(exchangerdict[i].control_volume, 'length_domain'):
            Tin_[i] = value(exchangerdict[i].control_volume.properties[0,0.0].temperature)
            Tout_[i] = value(exchangerdict[i].control_volume.properties[0,1.0].temperature)
            f_[i] = value(exchangerdict[i].control_volume.properties[0,1.0].flow_mol)
            Q_[i] = value(sum(exchangerdict[i].control_volume.heat[0,j] for j in exchangerdict[i].control_volume.length_domain if j != 0))*0.05
            T_units = pyunits.get_units(
                exchangerdict[i].control_volume.properties[0,0.0].temperature
            )
        else:
            Tin_[i] = value(exchangerdict[i].control_volume.properties_in[0].temperature)
            Tout_[i] = value(exchangerdict[i].control_volume.properties_out[0].temperature)
            f_[i] = value(exchangerdict[i].control_volume.properties_out[0].flow_mol)
            Q_[i] = value(exchangerdict[i].control_volume.heat[0])
            T_units = pyunits.get_units(
                exchangerdict[i].control_volume.properties_in[0].temperature
            )
        if hasattr(exchangerdict[i].control_volume, 'length_domain'):
            DG_units_[i] = pyunits.get_units(exchangerdict[i].heat_duty[0,1.0]) * pyunits.m
        else:
            DG_units_[i] = pyunits.get_units(exchangerdict[i].heat_duty[0])

    # Print the header
    print("Heat Exchanger Summary: ")

    # Print Inlet Temperature, Outlet Temperature and Heat
    for i in exchangerdict.keys():
        print("Heat exchanger: ", exchangerdict[i])
        print(f'Inlet T: {" "*3} {Tin_[i] : 0.3f} {T_units}')
        print(f'Outlet T: {" "*2} {Tout_[i] : 0.3f} {T_units}')
        print(f'Q : {" "*9} {Q_[i]: 0.3f} {DG_units_[i]}')

def return_HX_results(blk, exchanger_list):
    """
    Function to return results of Heat Exchangers
    Args:
        blk: flowsheet of the equipment
        exchanger_list: List of equipment to print data
    Returns:
        heat_ex_df: pandas dataframe with results from utility minimization
    """
    # Initialize dictionary and fill with exchanger list
    exchangerdict = {}

    for i, el in enumerate(exchanger_list):
        exchangerdict[str(el)] = el

    # initialize null dictionaries for data to be printed
    Tin_ = {}
    Tout_ = {}
    f_ = {}
    Q_ = {}
    DG_units_ = {}

    # Loop over heat exchangers
    for i in exchangerdict.keys():
        if hasattr(exchangerdict[i].control_volume, 'length_domain'):
            Tin_[i] = value(exchangerdict[i].control_volume.properties[0,0.0].temperature)
            Tout_[i] = value(exchangerdict[i].control_volume.properties[0,1.0].temperature)
            f_[i] = value(exchangerdict[i].control_volume.properties[0,1.0].flow_mol)
            Q_[i] = value(sum(exchangerdict[i].control_volume.heat[0,j] for j in exchangerdict[i].control_volume.length_domain if j != 0))*0.05
            T_units = pyunits.get_units(
                exchangerdict[i].control_volume.properties[0,0.0].temperature
            )
        else:
            Tin_[i] = value(exchangerdict[i].control_volume.properties_in[0].temperature)
            Tout_[i] = value(exchangerdict[i].control_volume.properties_out[0].temperature)
            f_[i] = value(exchangerdict[i].control_volume.properties_out[0].flow_mol)
            Q_[i] = value(exchangerdict[i].control_volume.heat[0])
            T_units = pyunits.get_units(
                exchangerdict[i].control_volume.properties_in[0].temperature
            )
        if hasattr(exchangerdict[i].control_volume, 'length_domain'):
            DG_units_[i] = pyunits.get_units(exchangerdict[i].heat_duty[0,1.0]) * pyunits.m
        else:
            DG_units_[i] = pyunits.get_units(exchangerdict[i].heat_duty[0])

    # pandas dataframe
    heat_ex_df = pd.DataFrame(columns=['Quantity','Unit']+[str(i) for i in exchangerdict.keys()])
    T_inlet_dict = {'Quantity':'T inlet',
                    'Unit':str(T_units)
                    }
    T_outlet_dict = {'Quantity':'T outlet',
                    'Unit':str(T_units)
                    }
    Q_dict = {'Quantity':'Heat duty'
                    }
    for i in exchangerdict.keys():
        T_inlet_dict[str(i)] = Tin_[i]
        T_outlet_dict[str(i)] = Tout_[i]
        Q_dict[str(i)] = Q_[i]
        Q_dict['Unit'] = str(DG_units_[i])
    heat_ex_df = pd.concat([heat_ex_df,pd.DataFrame([T_inlet_dict,T_outlet_dict,Q_dict])], ignore_index=True)
    return heat_ex_df
    


def unique(list1):
    """
    Function to remove not unique elements of a list
    Args:
        list1: List from where to obtain only unique values
    Returns:
        unique_list
    """
    # intilize a null list
    unique_list = []

    # traverse for all elements
    for i in list1:
        # check if exists in unique_list or not
        if i not in unique_list:
            unique_list.append(i)

    return unique_list


class PinchDataClass:
    """
    Class containing all the initialized values for all heat Variables
    """

    def __init__(self, initQs, initQw):
        """
        Args:
            initQs: Initial value of heat of hot utility
            initQw: Initial value of heat to be removed by cold utility
        Returns:
            None
        """
        self.initQs = initQs
        self.initQw = initQw
        self.initQAh = {}
        self.initQAc = {}
        self.T_ = {}

    def HeatAbove(self, data_dic):
        """
        Args:
            data_dic: Initial value of heat above pinch
        Returns:
            None
        """
        self.initQAh = data_dic

    def HeatBellow(self, data_dic2):
        """
        Args:
            data_dic2: Initial value of heat below pinch
        Returns:
            None
        """
        self.initQAc = data_dic2


class CurveData:
    """
    Class containing necessary data to generate composite curves
    """

    def __init__(self, Qw, T_unit, Q_unit):
        self.Qw = Qw
        self.T_unit = T_unit
        self.Q_unit = Q_unit

    def Cooling_Tin(self, list):
        Cooling_Tin = list

    def Cooling_Tout(self, list):
        Cooling_Tout = list

    def Cooling_Q(self, list):
        Cooling_Q = list

    def Heating_Tin(self, list):
        Heating_Tin = list

    def Heating_Tout(self, list):
        Heating_Tout = list

    def Heating_Q(self, list):
        Heating_Q = list


def return_data(blk,heating, cooling, hx_name_tuple, dTmin=10, eps=1):
    df = pd.DataFrame(columns= ['Unit','T_ (smoothed)','T_ (max_op)' , 'T_in (smoothed)','T_in (max_op)','T_out (smoothed)','T_out (max_op)','FCp/Theta (smoothed)', 'FCp/Theta (max_op)','QAh (smoothed)','QAh (max_op)','QAc (smoothed)','QAc (max_op)'])
    hx_data = heat_data(blk, heating, cooling)
    pd_obj = pinch_calc(heating,cooling,hx_data, dTmin,eps)
    heater_list = hx_name_tuple

    for hx in heater_list:
        unit_vals = {}
        unit = getattr(blk,hx)
        unit_vals['Unit'] = hx
        unit_vals['T_ (smoothed)'] = blk.T_[str(unit)]()
        unit_vals['T_ (max_op)'] = pd_obj.T_[str(unit)]
        unit_vals['T_in (smoothed)'] = blk.Tin[str(unit)]()
        unit_vals['T_in (max_op)'] = hx_data[str(unit)]['T_in']
        unit_vals['T_out (smoothed)'] = blk.Tout[str(unit)]()
        unit_vals['T_out (max_op)'] = hx_data[str(unit)]['T_out']
        unit_vals['FCp/Theta (smoothed)'] = blk.Theta[str(unit)]()
        unit_vals['FCp/Theta (max_op)'] = hx_data[str(unit)]['FCp_']
        unit_vals['QAh (smoothed)'] = blk.QAh[str(unit)]()
        unit_vals['QAh (max_op)'] = pd_obj.initQAh[str(unit)]
        unit_vals['QAc (smoothed)'] = blk.QAc[str(unit)]()
        unit_vals['QAc (max_op)'] = pd_obj.initQAc[str(unit)]

        df = pd.concat([df,pd.DataFrame([unit_vals])],ignore_index=True)
    
    return df