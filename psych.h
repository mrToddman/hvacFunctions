/*
 * psych.h
 *
 *  Created on: Dec 26, 2016
 *      Author: Todd Miller
 *
 */



#ifndef	PSYCH_H
#define PSYCH_H
#include <math.h>



double part_press( double P, double W )
/*
 * Function to compute partial vapor pressure in [kPa]
 * From page 6.9 equation 38 in ASHRAE Fundamentals handbook (2005)
 * P = ambient pressure [kPa]
 * W = humidity ratio [kg/kg dry air]
*/

{
	return P * W / (0.62198 + W);
}


double sat_press( double Tdb)
/*
 * Function to compute saturation vapor pressure in [kPa]
 * ASHRAE Fundamentals handbook (2005) p 6.2, equation 5 and 6
 * Tdb = Dry bulb temperature [degC]
 * Valid from -100C to 200 C
*/

{
	double
	TK = 173.15,
	C1 = -5674.5359,
	C2 = 6.3925247,
    C3 = -0.009677843,
    C4 = 0.00000062215701,
    C5 = 2.0747825E-09,
    C6 = -9.484024E-13,
    C7 = 4.1635019,
    C8 = -5800.2206,
    C9 = 1.3914993,
    C10 = -0.048640239,
    C11 = 0.000041764768,
    C12 = -0.000000014452093,
    C13 = 6.5459673;

    TK = Tdb + 273.15; //Converts from degC to degK

    if(TK <= 273.15)
    	{
    	return exp(C1 / TK + C2 + C3 * TK + C4 * pow(TK, 2) + C5 * pow(TK, 3) + C6 * pow(TK, 4) + C7 * log(TK)) / 1000;
    	}
	else
		{
		return exp(C8 / TK + C9 + C10 * TK + C11 * pow(TK, 2) + C12 * pow(TK, 3) + C13 * log(TK)) / 1000;
		}
}


double hum_rat(double Tdb, double Twb, double P)
/*
 * Function to calculate humidity ratio [kg H2O/kg air]
 * Given dry bulb and wet bulb temperature inputs [degC]
 * ASHRAE Fundamentals handbook (2005)
 * Tdb = Dry bulb temperature [degC]
 * Twb = Wet bulb temperature [degC]
 * P = Ambient Pressure [kPa]
 */

{
	double Pws = sat_press(Twb);
	double Ws = 0.62198 * Pws / (P - Pws);	// Equation 23, p6.8
	if(Tdb >= 0)
	{
		// Equation 35, p6.9
		return ((2501 - 2.326 * Twb) * Ws - 1.006 * (Tdb - Twb)) / (2501 + 1.86 * Tdb - 4.186 * Twb);
	}
	else
	{
		// Equation 37, p6.9
		return ((2830 - 0.24 * Twb) * Ws - 1.006 * (Tdb - Twb)) / (2830 + 1.86 * Tdb - 2.1 * Twb);
	}
}


double hum_rat2(double Tdb, double RH, double P)
/*
 * Function to calculate humidity ratio [kg H2O/kg air]
 * Given dry bulb and wet bulb temperature inputs [degC]
 * ASHRAE Fundamentals handbook (2005)
 * Tdb = Dry bulb temperature [degC]
 * RH = Relative Humidity [Fraction or %/100]
 * P = Ambient Pressure [kPa]
 */
{
	double Pws = sat_press(Tdb);
	return 0.62198 * RH * Pws / (P - RH * Pws); // Equation 22, 24, p6.8
}


double rel_hum(double Tdb, double Twb, double P)
/*
 * Calculates relative humidity ratio
 * ASHRAE Fundamentals handbood (2005)
 * Tdb = Dry bulb temperature [degC]
 * Twb = Wet bulb temperature [degC]
 * P = Ambient Pressure [kPa]
 */
{
	double W = hum_rat(Tdb, Twb, P);
	return part_press(P, W) / sat_press(Tdb); // Equation 24, p6.8
}


double rel_hum2(double Tdb, double W, double P)
/*
 * Calculates the relative humidity
 * Tdb = Dry bulb temperature [degC]
 * W = humidity ratio [kg/kg dry air]
 * P = ambient pressure [kPa]
 */
{
	return part_press(P, W) / sat_press(Tdb);
}
double wet_bulb(double Tdb, double RH, double P)
/*
 * Calculates the Wet Bulb temperature [degC]
 * Uses Newton-Rhapson iteration to converge quickly
 * Tdb = Dry bulb temperature [degC]
 * RH = Relative humidity ratio [Fraction or %]
 * P = Ambient Pressure [kPa]
 */
{
	double W_normal = hum_rat2(Tdb, RH, P);

	// Solve to within 0.001% accuracy using Newton-Rhapson
	double Wet_bulb = Tdb; // initialize at saturation
	double W_new = hum_rat(Tdb, Wet_bulb, P);

	do
		{
			double W_new2 = hum_rat(Tdb, Wet_bulb - 0.001, P);
			double dw_dtwb = (W_new - W_new2) / 0.001;
			double Wet_bulb = Wet_bulb - (W_new - W_normal) / dw_dtwb;
			W_new = hum_rat(Tdb, Wet_bulb, P);
		}
		while (abs((W_new - W_normal) / W_normal) > 0.00001);
	return Wet_bulb;

}

double enthalpy_air_h2o(double Tdb, double W)
/*
 * Calculates enthalpy in [kJ/kg dry air]
 * From 2005 ASHRAE Handbook - Fundamentals - SI P6.9 eqn 32
 * Tdb = Dry bulb temperature [degC]
 * W = Humidity Ratio [kg/kg dry air]
 */
{
	return 1.006 * Tdb + W * (2501 + 1.86 * Tdb);
}


double dew_point(double P, double W)
/*
 * Calculates dew point temperature [deg C]
 * From page 6.9 equation 39 and 40 in ASHRAE Fundamentals handbook (2005)
 * P = ambient pressure [kPa]
 * W = humidity ratio [kg/kg dry air]
 * Valid for Dew Points less than 93 C
 */
{
	double
    C14 = 6.54,
    C15 = 14.526,
    C16 = 0.7389,
    C17 = 0.09486,
    C18 = 0.4569;

	double Pw = part_press(P, W);
	double alpha = log(Pw);
	double Tdp1 = C14 + C15 * alpha + C16 * pow(alpha, 2) + C17 * pow(alpha,  3) + C18 * pow(Pw, 0.1984);
	double Tdp2 = 6.09 + 12.608 * alpha + 0.4959 * pow(alpha, 2);

	if (Tdp1 >= 0)
	{
		return Tdp1;
	}
	else
	{
		return Tdp2;
	}
}


double dry_air_density(double P, double Tdb, double W)
/*
 * Calculates dry air density [kg_dry_air/m^3]
 * From page 6.8 equation 28 ASHRAE Fundamentals handbook (2005)
 * P = pressure [kPa]
 * Tdb = Dry bulb temperature [degC]
 * W = humidity ratio [kg/kg dry air]
 *
 * Note that total density of air-h2o mixture is:
 * rho_air_h2o = rho_dry_air * (1 + W)
 */
{
	double R_da = 287.055; // gas constant for dry air
	return 1000 * P / (R_da * (273.15 + Tdb) * (1 + 1.6078 * W));
}


/*
 * Use these functions below to calculate atmospheric pressure
 * Try the MPL3115A2, BMP180, or T5403 pressure sensor from Sparkfun.com
 * to compute real-time pressure readings, or better if measured in an air
 * duct use two pressure sensors with a Dwyer Instruments 160E pitot tube.
 * Readings from the pitot tube will give static pressure and tip pressure.
 * See https://www.grc.nasa.gov/WWW/K-12/airplane/pitot.html to solve for
 * air speed.
 */


double STD_press(double elevation)
/*
 * Calculates the standard pressure [kPa]
 * elevation = height relative to sea level [m]
 * ASHRAE Fundamentals 2005 - chap 6, eqn 3
 * Valid from -5000m to 11000m
 */
{
	return 101.325 * pow(1 - 0.0000225577 * elevation, 5.2559);
}


double STD_temp(double elevation)
/*
 * Calculates the standard temperature [degC] at given elevation [m]
 * ASHRAE Fundamentals 2005 - chap 6, eqn 4
 * Valid from -5000m to 11000m
 */
{
	return 15 - 0.0065 * elevation;
}

double psych(double P, double Tdb, int inValue, int inType, int outType, int SIq)
{
/*
 * P is the barometric pressure in PSI or Pa.
 * Tdb is the dry bulb in F or C
 * inValue is another parameter of choice (Wet bulb, Dew point, RH, Humidity Ratio, or Enthalpy)
 * inType is the number that corresponds to your choice of InV's parameter (1 through 4 or 7 respectively)
 * outType is the value requested.  It should be an integer between 1 and 10 excluding 8.  See below
 * SIq is the unit selector.  0 is IP, 1 is SI


 * The choices for inType and outType are:

 * 1 Web Bulb Temp            F or C                              Valid for Input
 * 2 Dew point                F or C                              Valid for input
 * 3 RH                       between 0 and 1                     Valid for input
 * 4 Humidity Ratio           Mass Water/ Mass Dry Air            Valid for input
 * 5 Water Vapor Pressure     PSI or Pa
 * 6 Degree of Saturation     between 0 and 1
 * 7 Enthalpy                 BTU/lb dry air or kJ/kg dry air     Valid for input
 *     Warning 0 state for IP is ~0F, 0% RH ,and  1 ATM, 0 state for SI is 0C, 0%RH and 1 ATM
 * 8 NOT VALID, Should be entropy
 * 9 Specific Volume          ft^3/lbm or m^3/kg dry air
 * 10 Moist Air Density       lb/ft^3 or m^3/kg
 */

	double Twb, Dew, RH, W, h, out;

	if(SIq == 1)
	{
	    P = P / 1000;  // Turns Pa to kPA
	    switch(inType)
	    {
	    case 1:
	        Twb = inValue;
	        break;

	    case 2:
	    	Dew = inValue;
	    	break;

	    case 3:
	    	RH = inValue;
	    	break;

	    case 4:
	    	W = inValue;
	    	break;

	    case 7:
	    	h = inValue;
	    	break;
	    }
	}
	else  // This section turns US Customary Units to SI units
	{
		Tdb = (Tdb - 32) / 1.8;
		P = P * 4.4482216152605 / pow(0.0254, 2) / 1000;    // PSI to kPa  Conversion factor exact
		switch(inType)
		{
		case 1:
			Twb = (inValue- 32) / 1.8;			// F to C
			break;

		case 2:
			Dew = (inValue- 32) / 1.8;			// F to C
			break;

		case 3:
		   	RH = inValue;						// no need to change
		    break;

		case 4:
		    W = inValue;						// no need to change
		    break;

		case 7:
		    h = inValue * 1.055056 / 0.45359237 - 17.884444444;
		    // 1.055056 kJ/(ISO_BTU)  .45359237 kg/lb
		    // 17.884444 kJ/kg 0 pt difference [Dry air at 0C and  dry air at 0F are both 0 enthalpy in their respective units]
		    break;
			}

	    }

	if(outType == 3 || outType == 1)			// Find RH
	    switch(inType)
	    {
	    case 1:									// given Twb
	        RH = rel_hum(Tdb, Twb, P);
	        break;
	    case 2:									// given Dew
	        RH = sat_press(Dew) / sat_press(Tdb);
	        break;
	    case 3:									// given RH
	        break;
	        // RH already Set
	    case 4:									// given W
	        RH = part_press(P, W) / sat_press(Tdb);
	        break;
	    case 7:
	        W = (1.006 * Tdb - h) / (-(2501 + 1.86 * Tdb));
	        // Algebra from 2005 ASHRAE Handbook - Fundamentals - SI P6.9 eqn 32
	        RH = part_press(P, W) / sat_press(Tdb);
	        break;
	    }
	else										// find W
	    switch(inType)
	    {
	    case 1:									// Given Twb
	    	W = hum_rat(Tdb, Twb, P);
	    	break;
	    case 2:									// Given Dew
	        W = 0.621945 * sat_press(Dew) / (P - sat_press(Dew));
	        break;
	        // Equation taken from eq 20 of 2009 Fundamentals chapter 1
	    case 3:									// Given RH
	        W = hum_rat2(Tdb, RH, P);
	        break;
	    case 4:									// Given W
	        // W already known
	    case 7:									// Given h
	        W = (1.006 * Tdb - h) / (-(2501 + 1.86 * Tdb));
	        // Algebra from 2005 ASHRAE Handbook - Fundamentals - SI P6.9 eqn 32
	        break;
		}

		// P, Tdb, and W are now available
		switch(outType)
		{
		case 1:									// requesting Twb
			out = wet_bulb(Tdb, RH, P);
			break;
		case 2:									// requesting Dew
			out = dew_point(P, W);
			break;
		case 3:									// Request RH
			out = RH;
			break;
		case 4:									// Request W
			out = W;
			break;
		case 5:									// Request Pw
			out = part_press(P, W) * 1000;
			break;
		case 6:									// Request deg of sat
			out = W / hum_rat2(Tdb, 1, P);
			// the middle arg of Hum_rat2 is suppose to be RH.  RH is suppose to be 100%
			break;
		case 7:									// Request enthalpy
			out = enthalpy_air_h2o(Tdb, W);
			break;
		case 8:									// Request entropy
			out = -9999;
			// don't have equation for Entropy
			break;
		case 9:									// Request specific volume
	    	out = 1 / (dry_air_density(P, Tdb, W));
	    	break;
		case 10:								// Request density
	    	out = dry_air_density(P, Tdb, W) * (1 + W);
	    	break;
		}

		if(SIq == 0)							// Convert to IP
			switch(outType)
			{
			case 1:								// Temperature
				out = (1.8 * out) + 32;
				break;
			case 2:								// Temperature
				out = (1.8 * out) + 32;
				break;
			// OutNum 3 and 4 (RH and W) are unitless
			case 5:								// Pressure
				out = out * pow(0.0254, 2) / 4.448230531;
				break;
			case 7:								// Enthalpy
				out = (out + 17.88444444444) * 0.45359237 / 1.055056;
				// Warning, 0 convention changes.  Be careful with units.
				break;
			case 9:								// Specific Volume
	        	out = out * 0.45359265 / (pow(12 * 0.0254, 3));
	        	break;
			case 10:							// Density
				out = out * pow(12 * 0.0254, 3) / 0.45359265;
				break;
			}
	return out;
}


#endif
