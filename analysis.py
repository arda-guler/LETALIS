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

    a_nzlExp = params[6] # deg
    ROC_thrtDn = params[7] # m
    ROC_thrtUp = params[8] # m

    n_cochan = params[9] # number of coolant channels
    L_cochanInnerWallDist = params[10] # m
    L_cochanSideWall = params[11] # m
    L_cochanDepth = params[12] # m

    L_filmInject = params[13] # m
    mdot_filmInject = params[14] # m

    # - - - COMBUSTION / CEA - - -
    D_star = D_thrt # m
    mdot_chamber = params[15] # kg s-1
    P_c = params[16] # Pa
    r_c = ROC_thrtDn # m
    T_c = params[17] # K
    c_star = params[18] # m/s, CEA
    gasConductivity = params[19] # W m-1 K-1, CEA
    avgMolecularMass = params[20] # g mol-1

    T_w = params[21] # K, wall temp

    # - - - COMBUSTION CHAMBER INPUTS - - -
    visc_chm = params[22] # millipoise, CEA
    gamma_chm = params[23] # CEA

    # - - - THROAT INPUTS - - -
    visc_thrt = params[24] # millipoise, CEA
    gamma_thrt = params[25] # CEA

    # - - - MATERIALS - - -
    mtl_innerWall = get_material_by_name(params[26])
    mtl_outerShell = get_material_by_name(params[27])

    # - - - COOLANT - - -
    mtl_clt = get_material_by_name(params[28])
    mdot_clt = params[29]/n_cochan # kg s-1 (per channel)
    T_clt = params[30] # manifold coolant temp. K
    P_clt = params[31] # manifold coolant press. Pa

    # - - - ANALYSIS - - -
    fineness_vertical = params[32]
    time_end = params[33] # s
    time_step = params[34] # s
    n_steps = int(time_end/time_step)

    # calculate engine geometry
    geom_x, geom_y, x_step, engine_lengths = calculate_geometry(L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm,
                                                                a_nzlExp, ROC_thrtDn, ROC_thrtUp, fineness_vertical)

    # generate 3D object
    print("\nGenerating 3D model...")
    vis_model = generate_3D(geom_x, geom_y, n_cochan, L_cochanInnerWallDist, L_cochanSideWall, L_cochanDepth)

    # calculate Mach distribution
    print("Calculating Mach distribution...")
    subsonic_x, subsonic_M, supersonic_x, supersonic_M = calc_mach_num(L_engine, engine_lengths[4], T_c, gamma_thrt, avgMolecularMass, fineness_vertical,
                                                                       L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)

    # generate cylinders
    print("Generating segments...")
    m_engine = 0
    r_prev = None
    cylinders = []
    for i in range(fineness_vertical):
        x = i * x_step
        r_in = get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)
        r_clt = r_in + L_cochanInnerWallDist
        r_out = r_clt + L_cochanDepth
        a_clt = (360/n_cochan) - (360 * L_cochanSideWall)/(2 * pi * r_clt)
        Mach = get_mach_num_at(x, subsonic_M, subsonic_x, supersonic_M, supersonic_x, engine_lengths)
        new_cylinder = cylinder(x, r_in, r_out, x_step, n_cochan, r_clt, a_clt, mtl_innerWall, T_w, Mach, r_prev)
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
        film_exists = False
        mdot_clt_current = mdot_clt # revert to manifold mass flow

        cylinder_film_exists = [False] * len(cylinders)
        rT_layers = [1] * len(cylinders)
        T_aw = 0

        i_cylinder_film = -1
        mdot_film_current = mdot_filmInject
        # loop forwards when computing film cooling
        for cy in cylinders:
            i_cylinder_film += 1

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

            xd = cy.x - L_filmInject
            Hs = 0.025
            Kt = 0.05 * 10**(-2)
            r_full = cy.r_in
            r_surf = cy.r_in - cy.r_in * Hs
            
            mf0 = mdot_filmInject # surface layer film flow
            ms0 = mdot_chamber * (pi * (r_full**2 - r_surf**2)/A) # surface layer gas flow

            # calculate flow temperature at point
            T_gas = (1 + (gamma-1)/2 * M**2)**(-1) * T_c
                
            # calculate heat transfer coeff
            h_g = get_convection_coeff(D_star, vis, Cp, Pr, P_c, c_star, r_c, A_star, A, gamma, M, T_film, T_c)

            # film cooling starts here
            if not film_exists and cy.x >= L_filmInject:
                film_exists = True
                mdot_clt_current = mdot_clt - mdot_filmInject/n_cochan
                mdot_film_current = mdot_filmInject
                ms = ms0
                mf = mf0
                
            # there is film cooling!
            if mdot_filmInject > 0 and cy.x >= L_filmInject and mdot_film_current >= 0: # TODO: add this to material properties 
                stability_coeff = 0.6 # how do you even get this
                dT_film = ( h_g * (T_gas - T_film) * cy.get_A_chm() )
                dT_film *= (stability_coeff * mdot_filmInject * mtl_clt.get_specific_heat(T_film))**(-1)

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

                        bigM = Kt * ms / mf
                        xbarsquared = x / Hs
                        xeta = 1 - euler**(-bigM * xbarsquared)
                        ms = ms0 * (1-xeta/2) + mf0 * xeta/2
                        mf = ms0 * xeta/2 + mf0 * (1-xeta/2)
                        T_layer = T_film + xeta * (T_c - T_film)
                        rT_layers[i_cylinder_film] = xeta * 20 # There is an archvile!

            if t_step % 100 == 0:
                T_films[j].append(T_film)

        film_inject_point = False  
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
            if not film_inject_point and cy.x <= L_filmInject:
                film_inject_point = True
                T_film = T_clt_current
                mdot_clt_current -= mdot_filmInject/n_cochan

            if not cylinder_film_exists[i_cylinder_regen]: # no film cooling on this cylinder

                T_effective = T_film + rT_layers[i_cylinder_regen] * (T_gas - T_film)
                Q_in = h_g * (T_effective - cy.T) * cy.get_A_chm() * time_step
                Q_in_per_area = h_g * (T_effective - cy.T) # W per m2 (this is only for plotting)
                Q_in_full += Q_in

            else: # there is liquid film cooling on this cylinder
                Q_in = 0
                Q_in_per_area = 0

            # compute Reynold's number
            # https://en.wikipedia.org/wiki/Hydraulic_diameter
            wet_perimeter = (2 * pi * cy.r_clt) * (cy.a_clt/360) + (2 * pi * cy.r_out) * (cy.a_clt/360) + 2 * L_cochanDepth
            flow_area = cy.A_cochan_flow
            D_hydro = 4 * (flow_area / wet_perimeter)
            Reynolds_num = (mdot_clt_current * D_hydro) / (mtl_clt.get_viscosity(T_clt_current) * cy.A_cochan_flow)

            if not mdot_clt == 0: 
                # compute heat absorption into regen cooling channels
                #h_l = get_h_clt_kerosene(mtl_clt, T_clt_current, mdot_clt, D_hydro, cy)
                h_l = get_h_clt_dittus_boelter(mtl_clt, T_clt, mdot_clt_current, D_hydro, cy)
                Q_out = h_l * (cy.T - T_clt_current) * cy.get_A_clt() * time_step
                Q_out_full += Q_out

                # increase cylinder temp
                Q_net = Q_in - Q_out
                dT = Q_net/cy.get_heat_capacity()
                cy.T += dT

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

    plot_data(time_step, xs, cylinder_temps, coolant_temps, coolant_presses, Q_ins, Q_in_per_areas, Q_outs, Reynolds, Nusselts,
              T_gases, h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls, geom_x, geom_y,
              flow_areas, wet_perimeters, D_hydros, m_engine, L_skirt_chan_width, L_chamber_chan_width, L_min_chan_width,
              L_max_chan_width, engine_lengths, vis_model, mdot_clts, T_films, rT_layers_plot, T_effectives, coolant_press_drops,
              total_clt_press_drops, config_filename)

    if getchar:
        qc = input("Press Enter to move on...")
        clear_cmd_terminal()
        
    
