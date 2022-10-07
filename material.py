import math

class material:
    pass

# stainless steel 304L
# References: Choong S. Kim - Thermophysical Properties of Stainless Steels
class SS304L(material):
    def __init__(self):
        self.name = "Stainless Steel 304L"

    def get_name(self):
        return self.name

    def get_melting_point(self, unit="K"):
        if unit == "C":
            return 1400
        else:
            return 1673

    def get_thermal_conductivity(self, temp):
        # takes temperature in K
        # returns thermal conductivity in (W m-1 K-1)
        if temp < 1673:
            return (0.08116 + 0.0001618 * temp) * 100
        else:
            #print("304L melting!")
            return (0.1229 + (3.279 * 10**(-5)) * temp) * 100

    def get_thermal_diffusivity(self, temp):
        # takes temperature in K
        # returns thermal diffusivity in (m2 s-1)
        if temp < 1673:
            return (0.02276 + 3.285*10**(-5) * temp + 2.762*10**(-9) * (temp**2)) / 10000
        else:
            #print("304L melting!")
            return (0.02514 + 1.996*10**(-7) * temp + 2.386*10**(-9) * (temp**2)) / 10000

    def get_specific_heat(self, temp):
        # takes temperature in K
        # returns specific heat in J kg-1 K-1

        # initially calculated in cal g-1 K-1
        cgk = (0.1122 + 3.222*10**(-5) * temp)

        # do conversion when returning final value
        return cgk * 4.184 * 1000

    def get_density(self):
        # returns density in kg m-3
        return 8050

# copper chromium zirconium
# References: G. Pintsuk - Interlaboratory Test on Thermophysical Properties of the
#                          ITER Grade Heat Sink Material Copper–Chromium–Zirconium
#
# (Only data available was between 25 - 500 degrees C)
class CuCrZr(material):
    def __init__(self):
        self.name = "Copper-Chromium-Zirconium"

    def get_name(self):
        return self.name

    def get_melting_point(self, unit="K"):
        if unit == "C":
            return 1020
        else:
            return 1293

    def get_thermal_conductivity(self, temp):
        # takes temperature in K
        # returns thermal conductivity in W m-1 K-1
        return 358.07 * temp ** (-0.005)
        #return 353

    def get_specific_heat(self, temp):
        # takes temperature in K
        # returns specific heat in J kg-1 K-1
        return 0.0948 * temp + 367.97

    def get_density(self):
        # returns density in kg m-3
        return 8.75 * 1000

# Jet A-1
# References: James T. Edwards - Jet Fuel Properties
#                                AFRL-RQ-WP-TR-2020-0017
#                                January 2020 Interim Report
#
#             Tim Edwards - Reference Jet Fuels for Combustion Testing
class Jet_A1(material):
    def __init__(self):
        self.name = "Kerosene (Jet A1)"

    def get_name(self):
        return self.name

    def get_avg_MW(self):
        return 160

    def get_flash_point(self):
        return 273 + 50 # K

    def get_freeze_point(self):
        return 273 - 52 # K

    def get_density(self, temp):
        # takes temperature in K
        # returns specific heat in kg m-3
        # inter/extra-polated from Figure 9 of AFRL-RQ-WP-TR-2020-0017
        #return 0.805 * 1000 # kg m-3

        # initially calculate in g cm-3
        gcm = -0.0007 * temp + 1.0054

        # convert to kg m-3
        return gcm * 1000

    def get_heat_of_combustion(self):
        return 42.8 * 1000 # kJ kg-1

    def get_specific_heat(self, temp):
        # takes temperature in K
        # returns specific heat in J kg-1 K-1
        # inter/extra-polated from Figure 71 of AFRL-RQ-WP-TR-2020-0017

        # Figure 70
        return (0.0038 * temp + 0.832) * 1000

    def get_specific_gravity(self, temp):
        # Jet A-1 density divided by water density
        return self.get_density(temp)/997.77

    def get_viscosity(self, temp):
        # takes temperature in K
        # returns viscosity in Pa s

        # this function is based on a logarithmic excel fit
        # of Figure 10 from Edwards_AIAA-2017-0146_Reference_Jet_Fuels.pdf

        # this calculates viscosity in centistokes
        A = 19.8506
        B = 3.5317
        visc_cst = math.exp(math.exp(A-B*math.log(temp))) + 0.7

        # centipoise = centistokes * specific_gravity
        # now do the conversion
        visc_cP = visc_cst * self.get_specific_gravity(temp)

        # now convert to Pascal second
        # 1 centipoise = 0.001 Pascal second
        visc = visc_cP * 0.001

        return visc

    def get_thermal_conductivity(self, temp):
        # takes temperature in K
        # returns thermal conductivity in W m-1 K-1
        # based on Figure 82 of AFRL-RQ-WP-TR-2020-0017

        return -0.0002 * temp + 0.1553

    def get_heat_of_vaporization(self, temp):
        # taken from Edwards_AIAA-2017-0146_Reference_Jet_Fuels page 14
        # it is not a graph, just a simple text statement
        # should be good enough to roughly approximate things
        return 300000 # J kg-1

class water(material):
    def __init__(self):
        self.name = "water"

    def get_name(self):
        return self.name

    def get_density(self):
        return 997.77 # kg m-3
