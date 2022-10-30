# PROGRAM: RAO-STYLE THRUST OPTIMIZED BELL NOZZLE DESIGNER
# FOR LIQUID PROPELLANT ROCKET ENGINES
#
# Author: H. A. Güler (arda-guler @ Github)

import matplotlib.pyplot as plt
import math
import csv

from plot import *

def get_parabola_point(Nx, Ny, Qx, Qy, Ex, Ey, t):
    x = ((1-t)**2) * Nx + 2*(1-t)*t*Qx + (t**2) * Ex
    y = ((1-t)**2) * Ny + 2*(1-t)*t*Qy + (t**2) * Ey

    return x, y

def compute_bell_geometry(D_throat, D_exit, length_percent=80, theta_ch=30, theta_n=None, theta_e=None, x_fine=None, t_fine=None):

    xs = []
    ys = []

    if not length_percent:
        length_percent = 80

    if not x_fine:
        x_fine = 100

    if not t_fine:
        t_fine = 20

    if length_percent < 60:
        length_percent = 60
    
    R_throat = D_throat * 0.5
    R_exit = D_exit * 0.5
    expansion_ratio = (R_exit**2) / (R_throat**2)

    # theta_n not given, get it from Rao's graph
    if not theta_n:
        if length_percent <= 70:
            if expansion_ratio < 10:
                theta_n = 30
            else:
                theta_n = 35

        elif length_percent <= 85:
            if expansion_ratio < 30:
                theta_n = 25
            else:
                theta_n = 30

        else:
            if expansion_ratio <= 15:
                theta_n = 20
            else:
                theta_n = 25

    # theta_e not given, get it from Rao's graph
    if not theta_e:
        if length_percent <= 70:
            if expansion_ratio < 15:
                theta_e = 20
            else:
                theta_e = 15

        elif length_percent <= 85:
            if expansion_ratio < 10:
                theta_e = 15
            elif expansion_ratio < 40:
                theta_e = 10
            else:
                theta_e = 5

        else:
            if expansion_ratio < 6:
                theta_e = 10
            else:
                theta_e = 5

    # compute where the parabola starts
    # print("Parabola start angle:", theta_n)
    # print("Parabola end angle:", theta_e)
    x_throat = 0
    x_parabola = 0.382*R_throat*math.sin(math.radians(theta_n))
    x_exit = (length_percent/100) * (((expansion_ratio**0.5) - 1) * R_throat)/math.tan(math.radians(15))

    Nx = x_parabola
    Ny = R_throat + (R_throat * 0.382) - (((R_throat * 0.382)**2) - (x_parabola**2))**(0.5)

    Ex = x_exit
    Ey = R_exit

    m1 = math.tan(math.radians(theta_n))
    m2 = math.tan(math.radians(theta_e))
    C1 = Ny - m1*Nx
    C2 = Ey - m2*Ex

    Qx = (C2-C1)/(m1-m2)
    Qy = (m1*C2 - m2*C1)/(m1-m2)

    x = -math.sin(math.radians(theta_ch)) * (1.5 * R_throat)
    ## x = -R_throat * 1.5 / 2
    throat_dx = (x_parabola - x)/x_fine

    # throat downstream and upstream arcs
    while x < x_parabola:
        if x < 0:
            R_arc = R_throat * 1.5
        else:
            R_arc = R_throat * 0.382

        xs.append(x)
        ys.append(R_throat + R_arc - ((R_arc**2) - (x**2))**(0.5))

        x += throat_dx

    # parabola
    t = 0
    dt = 1/t_fine

    while t <= 1:
        x, y = get_parabola_point(Nx, Ny, Qx, Qy, Ex, Ey, t)
        xs.append(x)
        ys.append(y)
        t += dt

    # if the program somehow skips the very last point, compute it
    if t > 1 and t < 1 + dt:
        t = 1
        x, y = get_parabola_point(Nx, Ny, Qx, Qy, Ex, Ey, t)
        xs.append(x)
        ys.append(y)

    return xs, ys

##def design_and_analyze(params):
##    D_throat = params["Throat Diameter"]
##    D_exit = params["Exit Diameter"]
##    length_percent = params["% Length to Equivalent Cone (optional)"]
##    theta_n = params["Parabola Start Angle (optional)"]
##    theta_e = params["Exit Angle (optional)"]
##
##    x_fine = params["Throat Plot Fineness (optional)"]
##    t_fine = params["Parabola Fineness (optional)"]
##
##    gamma = params["Combustion Gases Gamma (optional)"]
##    
##    xs, ys = compute_geometry(D_throat, D_exit, length_percent, theta_n, theta_e,
##                              x_fine, t_fine)
##
##    mach_data = None
##
##    if gamma:
##        try:
##            subsonic_x, subsonic_M, supersonic_x, supersonic_M = calc_mach_num(xs, ys, gamma)
##            mach_data = [subsonic_x, subsonic_M, supersonic_x, supersonic_M]
##        except:
##            print("WARNING: Can not calculate Mach profile, likely because the program does some subtractive cancellation. Try reducing fineness.")
##
##    plot_all(xs, ys, mach_data)
##    return
