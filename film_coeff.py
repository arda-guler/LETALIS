import math
pi = math.pi
from material import SS304L, CuCrZr, Jet_A1

SS304L = SS304L()
CuCrZr = CuCrZr()
Jet_A1 = Jet_A1()

# Bartz Equation

# D. R. Bartz - A Simple Equation for Rapid Estimation of Rocket
#               Nozzle Convective Heat Transfer Coefficients

# NASA - Experimental Studies of the Heat Transfer to RBCC Rocket
#        Nozzles for CFD Application to Design Methodologies

def get_convection_coeff(D_star, vis, Cp, Pr, Pc, c_star, r_c, A_star, A, gamma, M, Tw, T0):

    # takes:
    # D_star = throat diameter in m
    # vis = viscosity in millipoise
    # Cp = specific heat at const. press. kJ kg-1 K-1
    # Pr = prandtl number
    # Pc = stagnation press in Pa
    # r_c = throat radius of curvature in m
    # c_star = characteristic velocity in m/s
    # A_star = throat area in m2
    # A = local area in m2
    # M = local mach number

    # gives:
    # h_g = convection coefficient in W m-2 K-1

    vis = vis * 0.0001 # convert millipoise to Pa.s
    Cp = Cp * 1000 # convert kJ kg-1 K-1 to J kg-1 K-1

    # sigma = correction factor for boundary layer
    sigma = 1 / ((0.5 * (Tw/T0) * (1 + (gamma-1)/2 * M**2) + 0.5)**(0.68)) * ((1 + (gamma-1)/2 * M**2)**(0.12))

    h_g = (0.026 / D_star**0.2) * (vis**0.2 * Cp)/Pr**0.6 * (Pc/c_star)**0.8 * (D_star/r_c)**0.1 * (A_star/A)**0.9 * sigma
    
    return h_g

# Lebedinsky E.V., Kalmykov G.P., et al. Working processes in liquid-propellant rocket
# engine and their simulation. Moscow, Mashinostroenie, 2008
def get_h_clt_kerosene(mtl_clt, T_clt, mdot_clt, D_hydro, cy):

    Pr = mtl_clt.get_specific_heat(T_clt) * mtl_clt.get_viscosity(T_clt) / mtl_clt.get_thermal_conductivity(T_clt)
    
    # compute Reynold's number
    # https://en.wikipedia.org/wiki/Hydraulic_diameter
    Reynolds_num = (mdot_clt * D_hydro) / (mtl_clt.get_viscosity(T_clt) * cy.A_cochan_flow)

    # compute Nusselt number
    Nusselt_num = 0.021 * (Reynolds_num**(0.8)) * (Pr**(0.4)) * (0.64 + 0.36 * (T_clt/cy.T))

    h_liq = Nusselt_num * mtl_clt.get_thermal_conductivity(T_clt) / D_hydro
    return h_liq

# Dittus-Boelter
def get_h_clt_dittus_boelter(mtl_clt, T_clt, mdot_clt, D_hydro, cy):

    clt_density = mtl_clt.get_density(T_clt)
    clt_thermoConduct = mtl_clt.get_thermal_conductivity(T_clt)
    clt_visc = mtl_clt.get_viscosity(T_clt)
    Prandtl = mtl_clt.get_specific_heat(T_clt) * mtl_clt.get_viscosity(T_clt) / mtl_clt.get_thermal_conductivity(T_clt)
    clt_vel = mdot_clt / (clt_density * cy.A_cochan_flow)
    
    Reynold = clt_density * clt_vel * D_hydro / clt_visc
    Nusselt = 0.023 * (Reynold**0.8) * (Prandtl**0.4)

    return Nusselt * clt_thermoConduct / D_hydro

def film_cooling_length(cy, x_inject, mdot_inject, n_cochan, mtl_clt, T_clt_inject, Pr, Cp, visc_gas, avgMolecularMass, T_gas):
    
    # Cp = gas spec heat = Cpg
    T_clt_saturation = 660 # K
    T_clt_film = T_clt_inject # K (assuming latent heat of vaporization does not change significantly, this is a constant)
    temp_diff = T_gas - T_clt_saturation
    
    L_circumference = 2 * pi * cy.r_in # m
    D = cy.r_in * 2 # m
    A = pi * cy.r_in**2 # m2
    mass_flow_per_area = (mdot_inject * n_cochan)/A # kg m-2 s-1
    mass_flow_per_circumference = (mdot_inject * n_cochan)/L_circumference # kg m-1 s-1

    MW_g = avgMolecularMass
    MW_c = mtl_clt.get_avg_MW()

    # lmbda = latent heat of vaporization of coolant
    lmbda = mtl_clt.get_heat_of_vaporization(T_clt_film)
    lmbda_star = lmbda + mtl_clt.get_specific_heat(T_clt_film) * (T_clt_saturation - T_clt_film)

    if MW_g > MW_c:
        a = 0.6
    else:
        a = 0.35
    Km = (MW_g / MW_c) ** a
    
    bigH = (Cp * temp_diff / lmbda_star) * Km
    
    # h_ratio = h/h0
    h_ratio = math.log(1+bigH)/bigH

    L_lqdFilm = (61.62 * mtl_clt.get_viscosity(T_clt_film))/(mass_flow_per_area) # unit becomes: Pa s m2 s kg-1
    L_lqdFilm *= ( (lmbda_star * mass_flow_per_circumference) / (Cp * temp_diff * visc_gas * h_ratio) )**(1.25)
    L_lqdFilm *= (Pr ** 0.78)

    return L_lqdFilm

