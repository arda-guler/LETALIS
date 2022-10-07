import math
pi = math.pi
uni_gas_const = uni_gas_const = 8.314472 # m2 kg s-2 K-1 mol-1

from geometry import *

pseudo_infinity = 10e6

# this function calculates the distribution of the combustion gas flow mach number
# across the engine
# see https://www.grc.nasa.gov/WWW/k-12/airplane/nozzled.html for more
# also especially see https://www.grc.nasa.gov/WWW/K-12/airplane/astar.html
# (the R here isn't the universal constant)
def calc_mach_num(x_end, x_thrt, Tc, gamma_thrt, avg_molecular_mass, fineness,
                  L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp):
    global uni_gas_const, pi

    subsonic_step_size = x_thrt/fineness

    subsonic_M = []
    subsonic_x = []
    
    # calculate for subsonic region (go backwards from throat)
    # m_dot = density * velocity * area
    # density and m_dot are const
    # so velocity changes with area
    for i in range(fineness):
        x = x_thrt - (i * subsonic_step_size)

        A_star = pi * get_inner_radius_at(x_thrt, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)**2
        A = pi * get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)**2
        A_ratio = A/A_star

        # Reference (1)
        # - - - - - - BEGIN QUOTE - - - - - -
        gp1 = gamma_thrt + 1
        gm1 = gamma_thrt - 1

        arat = A_ratio
        aro = 2
        macho = 0.30 #

        fac1 = gp1/(2*gm1)
        machn = macho + 0.05

        infinity_fuse = 0
        while abs(A_ratio - aro) > 0.0001:
            fac = 1 + 0.5 * gm1 * machn**2
            arn = 1/(machn * fac**(-fac1) * (gp1/2)**fac1)
            deriv = (arn-aro)/(machn-macho)
            aro = arn
            macho = machn
            machn = macho + (arat - aro)/deriv

        # - - - - - - END QUOTE - - - - - -

            infinity_fuse += 1

            if infinity_fuse > pseudo_infinity:
                print("Mach number calculator might have entered an infinite loop, because it hasn't converged for", pseudo_infinity, "iterations. (A)bort or (C)ontinue for another", pseudo_infinity, "iterations?")
                fuse_replacement = input(" > ")
                if fuse_replacement.lower() == "c":
                    infinity_fuse = 0
                elif fuse_replacement.lower() == "a":
                    print("Analysis aborted.")
                    input("Press Enter to quit...")
                    quit()
                else:
                    print("Invalid choice!")

        # required failsafe
        # NASA's code fails spectacularly if
        # mach number is too low or something
        if macho >= 1 or macho < 0.2:
            macho = 0
            
        subsonic_M.insert(0, macho)
        subsonic_x.insert(0, x)

    # calculate for supersonic region
    # density is variable
    # m_dot is const
    supersonic_M = []
    supersonic_x = []

    supersonic_step_size = (x_end - x_thrt)/fineness
    for i in range(fineness):
        x = x_thrt + (i * supersonic_step_size)
        
        A_star = pi * get_inner_radius_at(x_thrt, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)**2
        A = pi * get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)**2
        A_ratio = A/A_star

        # Reference (1)
        # - - - - - - BEGIN QUOTE - - - - - -
        gp1 = gamma_thrt + 1
        gm1 = gamma_thrt - 1

        arat = A_ratio
        aro = 2
        macho = 2.2

        fac1 = gp1/(2*gm1)
        machn = macho + 0.05

        infinity_fuse = 0
        while abs(A_ratio - aro) > .0001:
            fac = 1 + 0.5 * gm1 * machn**2
            arn = 1/(machn * fac**(-fac1) * (gp1/2)**fac1)
            deriv = (arn-aro)/(machn-macho)
            aro = arn
            macho = machn
            machn = macho + (arat - aro)/deriv
        # - - - - - - END QUOTE - - - - - -

            infinity_fuse += 1

            if infinity_fuse > pseudo_infinity:
                print("Mach number calculator might have entered an infinite loop, because it hasn't converged for", pseudo_infinity, "iterations. (A)bort or (C)ontinue for another", pseudo_infinity, "iterations?")
                fuse_replacement = input(" > ")
                if fuse_replacement.lower() == "c":
                    infinity_fuse = 0
                elif fuse_replacement.lower() == "a":
                    print("Analysis aborted.")
                    input("Press Enter to quit...")
                    quit()
                else:
                    print("Invalid choice!")
        
        supersonic_M.insert(0, macho)
        supersonic_x.insert(0, x)

    return subsonic_x, subsonic_M, supersonic_x, supersonic_M

# this one below isn't any rocket science
def get_index_of_closest_num_in_list(x, lst):
    min_diff = None
    closest_index = None
    for i in range(len(lst)):
        if not min_diff or abs(x - lst[i]) < min_diff:
            min_diff = abs(x - lst[i])
            closest_index = i
        
    return closest_index

# this function just reads values from the already-calculated mach number distribution
def get_mach_num_at(x, subsonic_mach, subsonic_x, supersonic_mach, supersonic_x, engine_lengths):

    # subsonic region
    if x < engine_lengths[4]:
        index = get_index_of_closest_num_in_list(x, subsonic_x)
        return subsonic_mach[index]

    # supersonic region
    elif engine_lengths[4] <= x <= engine_lengths[6]:
        index = get_index_of_closest_num_in_list(x, supersonic_x)
        return supersonic_mach[index]
