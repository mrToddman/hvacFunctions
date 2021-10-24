# FunctionsLibrary
This is a function library for work in the HVAC industry. Included are psychrometric state point and process functions, evaporation, heat transfer, etc.

The phychchrometric function: psych.h

P is the barometric pressure in PSI or Pa.
Tdb is the dry bulb in F or C
inValue is another parameter of choice (Wet bulb, Dew point, RH, Humidity Ratio, or Enthalpy)
inType is the number that corresponds to your choice of InV's parameter (1 through 4 or 7 respectively)
outType is the value requested.  It should be an integer between 1 and 10 excluding 8.  See below
SIq is the unit selector.  0 is IP, 1 is SI


The choices for inType and outType are:

1 Web Bulb Temp            F or C                              Valid for Input
2 Dew point                F or C                              Valid for input
3 RH                       between 0 and 1                     Valid for input
4 Humidity Ratio           Mass Water/ Mass Dry Air            Valid for input
5 Water Vapor Pressure     PSI or Pa
6 Degree of Saturation     between 0 and 1
7 Enthalpy                 BTU/lb dry air or kJ/kg dry air     Valid for input
    Warning 0 state for IP is ~0F, 0% RH ,and  1 ATM, 0 state for SI is 0C, 0%RH and 1 ATM
8 NOT VALID, Should be entropy
9 Specific Volume          ft^3/lbm or m^3/kg dry air
10 Moist Air Density       lb/ft^3 or m^3/kg

The evaporation function evap.h may be a fork to pool water evaporation and/or hydronic snowmelt calculator

Duct fittings and pressure losses may be a fork to Munual Q method
