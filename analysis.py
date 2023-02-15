# - - - - - - - - - - - - - - - - - - - -
# LIQUID PROPELLANT ROCKET ENGINE THERMAL
# ANALYSIS TOOL
# - - - - - - - - - - - - - - - - - - - -
# Program to compute transient heat
# transfers in an LRE.
# - - - - - - - - - - - - - - - - - - - -
# Authors:
# H. Arda Güler
# - - - - - - - - - - - - - - - - - - - -

import math
pi = math.pi
euler = math.e
import re
import time

from film_coeff import *
from geometry import *
from material import SS304L, CuCrZr, Jet_A1
from mach import *
from plot import *
from ui import *

SS = SS304L()
CCZ = CuCrZr()
JetA1 = Jet_A1()

materials = [SS, CCZ, JetA1]

def get_material_by_name(mtlname):
    global materials
    
    if mtlname == "SS":
        return materials[0]
    elif mtlname == "CCZ":
        return materials[1]
    elif mtlname == "Jet_A1":
        return materials[2]

def get_cylinder_index_at(x, L_engine, fineness_vertical):
    return int(x * fineness_vertical / L_engine)

def perform(params, config_filename=None, getchar=True):

    # - - - ENGINE GEOMETRY - - -
    L_engine = params[0] # m
    D_chm = params[1] # m
    D_thrt = params[2] # m
    D_exit = params[3] # m
    A_star = pi * (D_thrt/2)**2 # m2

    a_chmContract = params[4] # deg
    ROC_chm = params[5] # m

    type_nozzle = params[6]

    a_nzlExp = params[7] # deg
    if type_nozzle == "conic":
        ROC_thrtDn = params[8] # m
        ROC_thrtUp = params[9] # m
    else:
        ROC_thrtDn = (D_thrt/2) * 0.382
        ROC_thrtUp = (D_thrt/2) * 1.5

    percentLength_nzl = params[10]
    theta_n_nzl = params[11]
    theta_e_nzl = params[12]

    if type_nozzle == "bell":
        if not percentLength_nzl:
            percentLength_nzl = 80

        R_throat = D_thrt * 0.5
        R_exit = D_exit * 0.5
        expansion_ratio = (R_exit**2) / (R_throat**2)

        # theta_n not given, get it from Rao's graph
        if not theta_n_nzl:
            if percentLength_nzl <= 70:
                if expansion_ratio < 10:
                    theta_n_nzl = 30
                else:
                    theta_n_nzl = 35

            elif percentLength_nzl <= 85:
                if expansion_ratio < 30:
                    theta_n_nzl = 25
                else:
                    theta_n_nzl = 30

            else:
                if expansion_ratio <= 15:
                    theta_n_nzl = 20
                else:
                    theta_n_nzl = 25

        # theta_e not given, get it from Rao's graph
        if not theta_e_nzl:
            if percentLength_nzl <= 70:
                if expansion_ratio < 15:
                    theta_e_nzl = 20
                else:
                    theta_e_nzl = 15

            elif percentLength_nzl <= 85:
                if expansion_ratio < 10:
                    theta_e_nzl = 15
                elif expansion_ratio < 40:
                    theta_e_nzl = 10
                else:
                    theta_e_nzl = 5

            else:
                if expansion_ratio < 6:
                    theta_e_nzl = 10
                else:
                    theta_e_nzl = 5

    n_cochan = params[13] # number of coolant channels
    L_cochanInnerWallDist = params[14] # m
    L_cochanTangentialWidth = params[15] # m
    L_cochanDepth = params[16] # m

    L_filmInject1 = params[17] # m
    mdot_filmInject1 = params[18] # m

    L_filmInject2 = params[19] # m
    mdot_filmInject2 = params[20] # m

    # - - - COMBUSTION / CEA - - -
    D_star = D_thrt # m
    mdot_chamber = params[21] # kg s-1
    P_c = params[22] # Pa
    r_c = ROC_thrtDn # m
    T_c = params[23] # K
    c_star = params[24] # m/s, CEA
    gasConductivity = params[25] # W m-1 K-1, CEA
    avgMolecularMass = params[26] # g mol-1

    T_w = params[27] # K, wall temp

    # - - - COMBUSTION CHAMBER INPUTS - - -
    visc_chm = params[28] # millipoise, CEA
    gamma_chm = params[29] # CEA

    # - - - THROAT INPUTS - - -
    visc_thrt = params[30] # millipoise, CEA
    gamma_thrt = params[31] # CEA

    # - - - MATERIALS - - -
    mtl_innerWall = get_material_by_name(params[32])
    mtl_outerShell = get_material_by_name(params[33])

    # - - - COOLANT - - -
    mtl_clt = get_material_by_name(params[34])
    mdot_clt = params[35]/n_cochan # kg s-1 (per channel)
    T_clt = params[36] # manifold coolant temp. K
    P_clt = params[37] # manifold coolant press. Pa

    # - - - ANALYSIS - - -
    fineness_vertical = params[38]
    time_end = params[39] # s
    time_step = params[40] # s
    n_steps = int(time_end/time_step)

    # calculate engine geometry
    if type_nozzle == "conic":
        geom_x, geom_y, x_step, engine_lengths = calculate_geometry(L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm,
                                                                    a_nzlExp, ROC_thrtDn, ROC_thrtUp, fineness_vertical)
    else:
        geom_x, geom_y, x_step, engine_lengths = calculate_geometry_bell(L_engine, D_chm, D_thrt, D_exit, ROC_chm, a_chmContract,
                                                                    fineness_vertical, percentLength_nzl, theta_n_nzl, theta_e_nzl, 20, 1000)

    # generate 3D object
    print("\nGenerating 3D model...")
    vis_model = generate_3D_blade(geom_x, geom_y, n_cochan, L_cochanInnerWallDist, L_cochanTangentialWidth, L_cochanDepth,
                                  mdot_filmInject1, L_filmInject1, mdot_filmInject2, L_filmInject2)

    # calculate Mach distribution
    print("Calculating Mach distribution...")
    if type_nozzle == "conic":
        subsonic_x, subsonic_M, supersonic_x, supersonic_M = calc_mach_num(L_engine, engine_lengths[4], T_c, gamma_thrt, avgMolecularMass, fineness_vertical,
                                                                           L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)
    else:
        subsonic_x, subsonic_M, supersonic_x, supersonic_M = calc_mach_num_bell(L_engine, engine_lengths[4], T_c, gamma_thrt, avgMolecularMass, fineness_vertical,
                                                                                L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, percentLength_nzl, theta_n_nzl, theta_e_nzl)

    # generate cylinders
    print("Generating segments...")
    m_engine = 0
    r_prev = None
    cylinders = []

    if type_nozzle == "conic":
        for i in range(fineness_vertical):
            x = i * x_step
            r_in = get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)
            r_clt = r_in + L_cochanInnerWallDist
            r_out = r_clt + L_cochanDepth
            a_clt = L_cochanTangentialWidth
            b_clt = L_cochanDepth
            Mach = get_mach_num_at(x, subsonic_M, subsonic_x, supersonic_M, supersonic_x, engine_lengths)
            new_cylinder = cylinder(x, r_in, r_out, x_step, n_cochan, r_clt, a_clt, b_clt, mtl_innerWall, T_w, Mach, r_prev)
            cylinders.append(new_cylinder)
            m_engine += new_cylinder.get_m()

            r_prev = r_in
    else:
        for i in range(fineness_vertical):
            x = i * x_step
            r_in = get_inner_radius_at_bell(x, L_engine, D_chm, D_thrt, D_exit, ROC_chm, a_chmContract, percentLength_nzl, theta_n_nzl, theta_e_nzl)
            r_clt = r_in + L_cochanInnerWallDist
            r_out = r_clt + L_cochanDepth
            a_clt = L_cochanTangentialWidth
            b_clt = L_cochanDepth
            Mach = get_mach_num_at(x, subsonic_M, subsonic_x, supersonic_M, supersonic_x, engine_lengths)
            new_cylinder = cylinder(x, r_in, r_out, x_step, n_cochan, r_clt, a_clt, b_clt, mtl_innerWall, T_w, Mach, r_prev)
            cylinders.append(new_cylinder)
            m_engine += new_cylinder.get_m()

            r_prev = r_in

    # get important coolant channel widths
    L_skirt_chan_width = (2*pi*cylinders[-1].r_clt) * (cylinders[-1].a_clt/360)
    
    L_min_chan_width = None
    L_max_chan_width = None
    
    for cylin in cylinders:
        if not L_min_chan_width or (2*pi*cylin.r_clt) * (cylin.a_clt/360) < L_min_chan_width:
            L_min_chan_width = (2*pi*cylin.r_clt) * (cylin.a_clt/360)

        if not L_max_chan_width or (2*pi*cylin.r_clt) * (cylin.a_clt/360) > L_max_chan_width:
            L_max_chan_width = (2*pi*cylin.r_clt) * (cylin.a_clt/360)

    L_chamber_chan_width = (2*pi*cylinders[0].r_clt) * (cylinders[0].a_clt/360)

    # calculate Cp and Pr
    Cp_chm = (gamma_chm/(gamma_chm-1)) * uni_gas_const / avgMolecularMass # kJ kg-1 K-1, CEA
    Pr_chm = (4*gamma_chm) / (9*gamma_chm - 5) # unitless

    Cp_thrt = (gamma_thrt/(gamma_thrt-1)) * uni_gas_const / avgMolecularMass # kJ kg-1 K-1, CEA
    Pr_thrt = (4*gamma_thrt) / (9*gamma_thrt - 5) # unitless

    Q_ins = []
    Q_in_per_areas = []
    Q_outs = []
    xs = []
    cylinder_temps = []
    cylinder_temps_out = []
    cylinder_temps_in = []
    coolant_temps = []
    coolant_presses = []
    Reynolds = []
    Nusselts = []
    T_gases = []
    h_gs = []
    h_ls = []
    clt_vels = []
    Q_in_fulls = []
    Q_out_fulls = []
    flow_areas = []
    wet_perimeters = []
    D_hydros = []
    mdot_clts = []
    T_films = []
    rT_layers_plot = []
    T_effectives = []
    coolant_press_drops = []
    total_clt_press_drops = []
        
    time = 0
    j = 0

    T_film = 350 # initial guess for the first cycle
    
    for t_step in range(n_steps):

        if t_step % 100 == 0:
            Q_ins.append([])
            Q_in_per_areas.append([])
            Q_outs.append([])
            cylinder_temps.append([])
            cylinder_temps_out.append([])
            cylinder_temps_in.append([])
            coolant_temps.append([])
            coolant_presses.append([])
            Reynolds.append([])
            Nusselts.append([])
            h_gs.append([])
            h_ls.append([])
            T_films.append([])
            clt_vels.append([])
            rT_layers_plot.append([])
            T_effectives.append([])
            coolant_press_drops.append([])
            Q_in_full = 0
            Q_out_full = 0
            total_clt_press_drop = 0

        T_clt_current = T_clt # revert to manifold temperature
        P_clt_current = P_clt # revert to manifold pressure
        mdot_clt_current = mdot_clt # revert to manifold mass flow

        film_exists1 = False
        film_exists2 = False
        cylinder_film_exists = [False] * len(cylinders)
        rT_layers = [1] * len(cylinders) # T_film/T_gas ratio
        T_aw = 0 # adiabatic wall temp.

        i_cylinder_film = -1
        mdot_film_current = mdot_filmInject1
        # loop forwards when computing film cooling
        for cy in cylinders:
            i_cylinder_film += 1

            if cy.x < engine_lengths[4]: # before throat
                vis = visc_chm
                gamma = gamma_chm
                Cp = Cp_chm
                Pr = Pr_chm
            else: # after throat
                vis = visc_thrt
                gamma = gamma_thrt
                Cp = Cp_thrt
                Pr = Pr_thrt

            M = cy.get_Mach()
            A = pi * cy.r_in**2

            xd1 = cy.x - L_filmInject1
            xd2 = cy.x - L_filmInject2
            Hs = 0.025 * cy.r_in
            Kt = 0.05 * 10**(-2)

            # calculate flow temperature at point
            T_gas = (1 + (gamma-1)/2 * M**2)**(-1) * T_c
                
            # calculate heat transfer coeff
            h_g = get_convection_coeff(D_star, vis, Cp, Pr, P_c, c_star, r_c, A_star, A, gamma, M, T_film, T_c)

            if mdot_filmInject1:
                # no film yet
                if cy.x < L_filmInject1:
                    film_exists1 = False
                    film_exists2 = False

                # first injection here
                if not film_exists1 and cy.x >= L_filmInject1:
                    film_exists1 = True
                    mdot_clt_current = mdot_clt - (mdot_filmInject1 / n_cochan)
                    mdot_film_current = mdot_filmInject1

                # second injection here
                if not film_exists2 and cy.x >= L_filmInject2:
                    film_exists2 = True
                    mdot_clt_current -= (mdot_filmInject2 / n_cochan)
                    mdot_film_current += mdot_filmInject2

                if cy.x >= L_filmInject1 and mdot_film_current >= 0:
                    
                    mdot_totalFilm = mdot_filmInject1
                        
                    if cy.x >= L_filmInject2:
                        mdot_totalFilm += mdot_filmInject2
                        
                    stability_coeff = 0.6
                    dT_film = ( h_g * (T_gas - T_film) * cy.get_A_chm() )
                    dT_film *= (stability_coeff * mdot_totalFilm * mtl_clt.get_specific_heat(T_film))**(-1)

                    if T_film + dT_film < 600: # TODO: add this to material properties 
                        T_film += dT_film # increase film temperature
                        cylinder_film_exists[i_cylinder_film] = True

                    else: # check vaporization
                        dmdot_film = ( h_g * (T_gas - T_film) * cy.get_A_chm() ) / mtl_clt.get_heat_of_vaporization(T_film)

                        if mdot_film_current - dmdot_film > 0: # still not completely vaporized
                            mdot_film_current -= dmdot_film
                            cylinder_film_exists[i_cylinder_film] = True

                        else: # liquid film has completely vaporized (and is now in gas form)
                            cylinder_film_exists[i_cylinder_film] = False
                            mdot_film_current = 0

                            if cy.x >= L_filmInject1 and cy.x < L_filmInject2:
                                mbar_f = mdot_filmInject1 / mdot_chamber
                                A_surface_layer = pi * cy.r_in**2 - pi * (cy.r_in - Hs)**2
                                mdot_surface_layer = (A_surface_layer/(pi * cy.r_in**2)) * mdot_chamber
                                mbar_s = mdot_surface_layer / mdot_chamber
                                x_squared = xd1/Hs
                                bigM = Kt * (mbar_s/mbar_f)
                                xeta = 1 - euler**(-x_squared * bigM)
                                T_layer = (T_film + xeta * (T_c - T_film)) 
                                rT_layers[i_cylinder_film] = xeta

                            elif cy.x >= L_filmInject2 and cy.x >= L_filmInject2:
                                mbar_f = (mdot_filmInject1 + mdot_filmInject2) / mdot_chamber
                                A_surface_layer = pi * cy.r_in**2 - pi * (cy.r_in - Hs)**2
                                mdot_surface_layer = (A_surface_layer/(pi * cy.r_in**2)) * mdot_chamber
                                mbar_s = mdot_surface_layer / mdot_chamber
                                x_squared = xd2/Hs
                                bigM = Kt * (mbar_s/mbar_f)
                                xeta = 1 - euler**(-x_squared * bigM)
                                T_layer = (T_film + xeta * (T_c - T_film)) 
                                rT_layers[i_cylinder_film] = xeta

            if t_step % 100 == 0:
                T_films[j].append(T_film)

        film_inject_point1 = False
        film_inject_point2 = False
        i_cylinder_regen = len(cylinders)
        # loop backwards when computing cylinder heat transfer (go from manifold to injector face)
        for cy in cylinders[::-1]:
            i_cylinder_regen -= 1
            
            if cy.x < engine_lengths[4]:
                vis = visc_chm
                gamma = gamma_chm
                Cp = Cp_chm
                Pr = Pr_chm
            else:
                vis = visc_thrt
                gamma = gamma_thrt
                Cp = Cp_thrt
                Pr = Pr_thrt

            M = cy.get_Mach()
            A = pi * cy.r_in**2

            # calculate flow temperature at point
            T_gas = (1 + (gamma-1)/2 * M**2)**(-1) * T_c
                
            # calculate heat transfer
            h_g = get_convection_coeff(D_star, vis, Cp, Pr, P_c, c_star, r_c, A_star, A, gamma, M, cy.T, T_c)

            # get film cooling injection point temperature
            if not film_inject_point1 and cy.x <= L_filmInject1:
                film_inject_point1 = True
                T_film = T_clt_current
                mdot_clt_current -= mdot_filmInject1/n_cochan

            if not film_inject_point2 and cy.x <= L_filmInject2:
                film_inject_point2 = True
                T_film = T_clt_current
                mdot_clt_current -= mdot_filmInject2/n_cochan

            if not cylinder_film_exists[i_cylinder_regen]: # no film cooling on this cylinder

                T_effective = min(T_film + rT_layers[i_cylinder_regen] * (T_gas - T_film) * 3, T_gas) # Archvile!
                Q_in = h_g * (T_effective - (cy.T + cy.T_diff/2)) * cy.get_A_chm() * time_step
                Q_in_per_area = h_g * (T_effective - (cy.T + cy.T_diff/2)) # W per m2 (this is only for plotting)
                Q_in_full += Q_in

            else: # there is liquid film cooling on this cylinder
                Q_in = 0
                Q_in_per_area = 0

            # compute Reynold's number
            # https://en.wikipedia.org/wiki/Hydraulic_diameter
            wet_perimeter = cy.a_clt + 2 * cy.b_clt
            flow_area = cy.A_cochan_flow
            D_hydro = 4 * (flow_area / wet_perimeter)
            Reynolds_num = (mdot_clt_current * D_hydro) / (mtl_clt.get_viscosity(T_clt_current) * cy.A_cochan_flow)

            if not mdot_clt == 0: 
                # compute heat absorption into regen cooling channels
                #h_l = get_h_clt_kerosene(mtl_clt, T_clt_current, mdot_clt, D_hydro, cy)
                h_l = get_h_clt_dittus_boelter(mtl_clt, T_clt, mdot_clt_current, D_hydro, cy)
                Q_out = h_l * ((cy.T - cy.T_diff/2) - T_clt_current) * cy.get_A_clt() * time_step
                Q_out_full += Q_out

                # increase cylinder temp
                Q_net = Q_in - Q_out
                dT = Q_net/cy.get_heat_capacity()
                cy.T += dT
                cy.T_diff = Q_net * cy.get_thermal_resistance() / time_step

                # increase coolant fluid temp.
                clt_vel = mdot_clt_current / (mtl_clt.get_density(T_clt_current) * cy.A_cochan_flow)
                dT_clt = (Q_out/(n_cochan * time_step)) / (mdot_clt_current * mtl_clt.get_specific_heat(T_clt_current))
                T_clt_current += dT_clt

                # compute Nusselt number (Dittus Boelter)
                Pr_clt = mtl_clt.get_specific_heat(T_clt) * mtl_clt.get_viscosity(T_clt) / mtl_clt.get_thermal_conductivity(T_clt)
                Nusselt_num = 0.023 * Reynolds_num**0.8 * Pr_clt**0.3

                # compute coolant pressure drop and update pressures
                epsilon_f = 1
                
                if Reynolds_num <= 2320:
                    friction_loss_coeff = (64/Reynolds_num) * epsilon_f
                elif Reynolds_num < 10E5:
                    friction_loss_coeff = (0.3164/(Reynolds_num**(1/4))) * epsilon_f
                else:
                    friction_loss_coeff = (0.0032 + (0.221/(Reynolds_num**(0.237)))) * epsilon_f
                    
                coolant_press_drop = friction_loss_coeff * (cy.h/D_hydro) * mtl_clt.get_density(T_clt_current) * ((clt_vel**2)/2)
                P_clt_current -= coolant_press_drop
                total_clt_press_drop += coolant_press_drop
                    
            else:
                h_l = 0
                Q_out = 0
                Q_net = Q_in

                # increase cylinder temp (no cooling)
                dT = Q_net/cy.get_heat_capacity()
                cy.T += dT
                
                clt_vel = 0
                m_flow = 0
                dT_clt = 0
                T_clt_current += 0
                Pr_clt = 0
                Nusselt_num = 0
                coolant_press_drop = 0

            # axial conduction of heat (positive direction is from nozzle exit towards injector face)
            if cylinders.index(cy) - 1 >= 0:
                cy_upper = cylinders[cylinders.index(cy) - 1]
                A_axial = get_area_of_sector(cy.r_in, cy.r_out, 360) - (cy.A_cochan_flow * n_cochan)
                k_axial = cy.mtl.get_thermal_conductivity(cy.T)
                dx_axial = cy.h
                dT_axial = cy.T - cy_upper.T
                #print(A_axial, k_axial, dx_axial, dT_axial)
                Q_axial = (k_axial * A_axial * dT_axial / dx_axial) * time_step
                cy.T -= Q_axial/cy.get_heat_capacity()
                cy_upper.T += Q_axial/cy.get_heat_capacity()
                
            # record data for plotting
            if t_step == 0:
                xs.insert(0, cy.x)
                T_gases.insert(0, T_gas)
                flow_areas.insert(0, flow_area)
                wet_perimeters.insert(0, wet_perimeter)
                D_hydros.insert(0, D_hydro)
                mdot_clts.insert(0, mdot_clt_current)
            if t_step % 100 == 0:
                Q_ins[j].insert(0, Q_in/time_step) # convert to W
                Q_in_per_areas[j].insert(0, Q_in_per_area)
                Q_outs[j].insert(0, Q_out/time_step) # convert to W
                cylinder_temps[j].insert(0, cy.T - 273) # convert to celcius
                cylinder_temps_out[j].insert(0, cy.T - cy.T_diff/2 - 273) # convert to celcius
                cylinder_temps_in[j].insert(0, cy.T + cy.T_diff/2 - 273) # convert to celcius
                coolant_temps[j].insert(0, T_clt_current - 273) # convert to celcius
                coolant_presses[j].insert(0, P_clt_current)
                Reynolds[j].insert(0, Reynolds_num)
                Nusselts[j].insert(0, Nusselt_num)
                h_gs[j].insert(0, h_g)
                h_ls[j].insert(0, h_l)
                clt_vels[j].insert(0, clt_vel)
                rT_layers_plot[j].insert(0, rT_layers[i_cylinder_regen])
                T_effectives[j].insert(0, T_effective)
                coolant_press_drops[j].insert(0, coolant_press_drop)

        if t_step % 100 == 0:
            Q_in_fulls.append(Q_in_full)
            Q_out_fulls.append(Q_out_full)
            total_clt_press_drops.append(total_clt_press_drop)

        # proceed to next time step
        time += time_step
        if t_step % 100 == 0:
            j+=1

            clear_cmd_terminal()
            print("")
            print("= = = SINGLE THERMAL ANALYSIS = = =")
            print("")
            print("Current analysis:")
            print(generate_progress_bar((t_step/n_steps) * 100))

    plot_data(time_step, xs, cylinder_temps, cylinder_temps_out, cylinder_temps_in, coolant_temps, coolant_presses, Q_ins,
              Q_in_per_areas, Q_outs, Reynolds, Nusselts, T_gases, h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls, geom_x, geom_y,
              flow_areas, wet_perimeters, D_hydros, m_engine, L_skirt_chan_width, L_chamber_chan_width, L_min_chan_width,
              L_max_chan_width, engine_lengths, mdot_clts, T_films, rT_layers_plot, T_effectives, coolant_press_drops,
              total_clt_press_drops, vis_model, config_filename)

    if getchar:
        qc = input("Press Enter to move on...")
        clear_cmd_terminal()
        
    
