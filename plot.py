import time
import datetime
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil

def plot_data(time_step, xs, cylinder_temps, cylinder_temps_out, cylinder_temps_in, coolant_temps, coolant_presses,
              Q_ins, Q_in_per_areas, Q_outs, Reynolds, Nusselts, T_gases, h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls,
              geom_x, geom_y, flow_areas, wet_perimeters, D_hydros, m_engine, L_skirt_chan_width, L_chamber_chan_width,
              L_min_chan_width, L_max_chan_width, engine_lengths, mdot_clts, T_films, rT_layers_plot, T_effectives,
              coolant_press_drops, total_clt_press_drops, vis_model, filename=None):

    # stop max. number of figures warning (because death to your RAM, that's why!)
    plt.rcParams.update({'figure.max_open_warning': 0})

    # PRINT TOTAL Q
    Q_in_total = 0
    for Q in Q_in_fulls:
        Q_in_total += Q

    Q_out_total = 0
    for Q in Q_out_fulls:
        Q_out_total += Q

##    print("\nTotal Q_in:", Q_in_total * 100)
##    print("Total Q_out:", Q_out_total * 100)
##    print("Net Q:", (Q_in_total - Q_out_total) * 100)

    def plot_engine_contour(ax):
        ax.set_aspect("equal")
        ax.set_ylim([-max(geom_y)*0.7, max(max(geom_x), max(geom_y))])
        ax.plot(geom_x, geom_y, color="black")
        ax.fill_between(geom_x, geom_y, where=[True]*len(geom_x), interpolate=True, color='black')
        ax.set_ylabel("Engine Contour (m)")
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
    
    num_frames = len(cylinder_temps)
    fig, ax = plt.subplots()
    plotnum = 0

    # WALL TEMP. PLOT
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()
    
    for i in range(0, num_frames, int(num_frames/10)):
        red = max(min(1, max(cylinder_temps[i])/600),0)
        blue = 1 - red
        ax2.plot(xs, cylinder_temps[i], color=(red, 0, blue))

    plt.grid()
    plt.title("Wall Temperature")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Temperature (C)")

    # WALL OUTER - INNER TEMP. PLOT (FOR THE LAST FRAME)
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()
    
    ax2.plot(xs, cylinder_temps_out[-1], color=(0,0,1))
    ax2.plot(xs, cylinder_temps_in[-1], color=(1,0,0))

    plt.grid()
    plt.title("Inner & Outer Wall Surface Temperatures")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Temperature (C)")
    
    # COOLANT TEMP. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        red = max(min(1, max(coolant_temps[i])/350),0)
        blue = 1 - red
        ax2.plot(xs, coolant_temps[i], color=(red, 0, blue))

    plt.grid()
    plt.title("Coolant Temperature")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Temperature (C)")

    # COOLANT PRESS. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, coolant_presses[i])

    plt.grid()
    plt.title("Coolant Pressure")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Pressure (Pa)")

    # HEAT PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, Q_ins[i], color=(1,0,0))
        ax2.plot(xs, Q_outs[i], color=(0,0,1))

    plt.grid()
    plt.title("Heat Transfer")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Heat Transferred (W)")

    # HEAT PER AREA PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, Q_in_per_areas[i], color=(1,0,0))

    plt.grid()
    plt.title("Heat Transfer Per Area")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Heat Transferred (W m-2)")

    # REYNOLDS NUMBER PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, Reynolds[i])

    plt.grid()
    plt.title("Reynolds Number")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Reynolds Number")

    # NUSSELT NUMBER PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, Nusselts[i])

    plt.grid()
    plt.title("Nusselts Number")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Nusselts Number")

    # COMBUSTION GAS TEMP. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    ax2.plot(xs, T_gases)

    plt.grid()
    plt.title("Gas Temperature (K)")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Gas Temperature")

    # GAS CONVECTION COEFF. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, h_gs[i])

    plt.grid()
    plt.title("hg")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("hg")

    # LIQUID FILM COEFF. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, h_ls[i])

    plt.grid()
    plt.title("hl")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("hl")

    # COOLANT VELOCITY PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, clt_vels[i])

    plt.grid()
    plt.title("Coolant Flow Velocity")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Velocity (m s-1)")

    # TOTAL Q IN PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    plt.plot(Q_in_fulls)

    plt.grid()
    plt.title("Total Q In")
    plt.xlabel("Time")
    plt.ylabel("Total Q In")

    # TOTAL Q OUT PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    plt.plot(Q_out_fulls)

    plt.grid()
    plt.title("Total Q Out")
    plt.xlabel("Time")
    plt.ylabel("Total Q Out")

##    # HEAT PLOT (3D)
##    plt.figure(12)
##    ax = plt.axes(projection='3d')
##    
##    times = []
##    for i in range(num_frames):
##        times.append(i*time_step)
##
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_ins[i], xs, times[i], color="red")
##
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_outs[i], xs, times[i], color="blue")
##
##    # HEAT DIFF. PLOT (3D)
##    plt.figure(13)
##    ax = plt.axes(projection='3d')
##
##    Q_nets = []
##    for i in range(len(Q_ins)):
##        Q_nets.append([])
##        for j in range(len(Q_ins[0])):
##            Q_nets[i].append(Q_ins[i][j] - Q_outs[i][j])
##    
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_nets[i], xs, times[i], color="green")

    # ENGINE GEOMETRY
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax.set_aspect('equal')
    
    geom_y_negative = []
    for y in geom_y:
        geom_y_negative.append(-y)
        
    ax.plot(geom_x, geom_y)
    ax.plot(geom_x, geom_y_negative)

    plt.grid()
    plt.title("Thrust Chamber Geometry")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # FLOW AREA
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)
    plt.plot(xs, flow_areas)
    plt.grid()
    plt.title("Coolant Flow Area")
    plt.xlabel("X")
    plt.ylabel("Area m2")

    # WET PERIMETER
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)
    plt.plot(xs, wet_perimeters)
    plt.grid()
    plt.title("Wet Perimeter")
    plt.xlabel("X")
    plt.ylabel("Perimeter m")

    # HYDRAULIC DIAMETERS
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)
    plt.plot(xs, D_hydros)
    plt.grid()
    plt.title("Hydraulic Diameter")
    plt.xlabel("X")
    plt.ylabel("Diameter m")

    # COOLANT MDOT PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    ax2.plot(xs, mdot_clts)

    plt.grid()
    plt.title("Cooling Channel Mass Flow")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Flow (kg s-1)")

    # FILM TEMPERATURE PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, T_films[i])

    plt.grid()
    plt.title("Film Temperature")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Film Temperature (K)")

    # LAYER RATIO PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, rT_layers_plot[i])

    plt.grid()
    plt.title("Layer Comp.")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("S/F Ratio")

    # EFFECTIVE TEMP. PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, T_effectives[i])

    plt.grid()
    plt.title("Effective Temp.")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Effective Temp. (K)")

    # COOLANT PRESS. DROP PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    ax2 = ax.twinx()
    plot_engine_contour(ax)
    ax2.set_aspect("auto")
    ax2.yaxis.set_label_position("left")
    ax2.yaxis.tick_left()

    for i in range(0, num_frames, int(num_frames/10)):
        ax2.plot(xs, coolant_press_drops[i])

    plt.grid()
    plt.title("Coolant Pressure Drop")
    ax.set_xlabel("Position (m)")
    ax2.set_ylabel("Pressure Drop")

    # TOTAL COOLANT PRESS. DROP PLOT
    _, ax = plt.subplots()
    plotnum += 1
    plt.figure(plotnum)

    plt.plot(total_clt_press_drops)

    plt.grid()
    plt.title("Total Coolant Pressure Drop")
    plt.xlabel("Time")
    plt.ylabel("Total Coolant Pressure Drop (Pa)")

    ### ANIMATED FIGURES
    anims = []
    anim_plotnum = 0
    
    # WALL TEMP ANIMATED
    fig, ax = plt.subplots()
    anim_plotnum += 1
    plt.figure(plotnum + anim_plotnum)
    plt.title("Wall Temp.")
    plt.xlabel("X Position (mm)")
    plt.ylabel("Temperature (C)")
    plt.grid()

    plt.xlim = (0, xs[-1])
    plt.ylim(0, max(cylinder_temps[-1]) + 50)
    wall_temp_line, = ax.plot(xs, cylinder_temps[-1])

    def animate_wall_temp(i):
        wall_temp_line.set_ydata(cylinder_temps[i])
        return wall_temp_line,

    anim_wall_temp = animation.FuncAnimation(
        fig, animate_wall_temp, blit=True, frames=list(range(len(cylinder_temps))))

    anims.append(anim_wall_temp)

    # COOLANT TEMP ANIMATED
    fig, ax = plt.subplots()
    anim_plotnum += 1
    plt.figure(plotnum + anim_plotnum)
    plt.title("Coolant Temp.")
    plt.xlabel("X Position (mm)")
    plt.ylabel("Coolant Temp. (C)")
    plt.grid()

    plt.xlim = (0, xs[-1])
    plt.ylim(0, max(coolant_temps[-1]) + 50)
    coolant_temp_line, = ax.plot(xs, coolant_temps[-1])

    def animate_coolant_temp(i):
        coolant_temp_line.set_ydata(coolant_temps[i])
        return coolant_temp_line,

    anim_coolant_temp = animation.FuncAnimation(
        fig, animate_coolant_temp, blit=True, frames=list(range(len(coolant_temps))))

    anims.append(anim_coolant_temp)

    # HEAT FLOWS TEMP ANIMATED
    fig, ax = plt.subplots()
    anim_plotnum += 1
    plt.figure(plotnum + anim_plotnum)
    plt.title("Heat Transferred")
    plt.xlabel("X Position (mm)")
    plt.ylabel("Heat (W)")
    plt.grid()

    plt.xlim = (0, xs[-1])
    plt.ylim(min(Q_outs[-1]), max(Q_ins[-1]) + 500)
    Q_ins_line, = ax.plot(xs, Q_ins[-1])
    Q_outs_line, = ax.plot(xs, Q_outs[-1])

    def animate_heat_transfer(i):
        Q_ins_line.set_ydata(Q_ins[i])
        Q_outs_line.set_ydata(Q_outs[i])

        return Q_ins_line, Q_outs_line,

    anim_heat_trans = animation.FuncAnimation(
        fig, animate_heat_transfer, blit=True, frames=list(range(len(Q_ins))))

    anims.append(anim_heat_trans)

    print("")

    if not filename:
        folder_name = "heat_analysis_" + datetime.datetime.now().strftime("%y%m%d%H%M%S")
    else:
        if "/" in filename:
            folder_name = filename.split("/")[2].split(".")[0]
        else:
            folder_name = filename.split(".")[0]

    print("Exporting data to folder: " + folder_name)

    try:
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    except:
        print("ERROR: Could not create folder. Try saving figures by hand.")
        plt.show()
        return

    try:
        with open(str(folder_name + "/geometry.txt"), "w") as f:
            f.write("xs=" + str(geom_x))
            f.write("\n\n")
            f.write("ys=" + str(geom_y))
            f.write("\n\n")
            f.write("Engine mass (kg): " + str(m_engine) + "\n\n")
            #f.write("Min. coolant channel width (m): " + str(L_min_chan_width) + "\n")
            #f.write("Max. coolant channel width (m): " + str(L_max_chan_width) + "\n")
            #f.write("Chamber coolant channel width (m): " + str(L_chamber_chan_width) + "\n")
            #f.write("Skirt coolant channel width (m): " + str(L_skirt_chan_width) + "\n\n")
            f.write("Engine contour change positions (m) = " + str(engine_lengths))
    except:
        print("WARNING: Could not export geometry data.")

    try:
        with open(str(folder_name + "/3d_model.txt"), "w") as f:
            for vertex in vis_model:
                f.write(str(vertex)+"\n")
    except:
        print("WARNING: Couldn't export 3D model!")

    try:
        shutil.copy(filename, folder_name)
    except:
        print("WARNING: Could not copy inputs file to analysis folder.")

    try:
        shutil.copy("modelviewer.py", folder_name)
    except:
        print("WARNING: Could not copy modelviewer utility to analysis folder.")

    try:
        for i in range(1, plotnum + 1):
            new_fig = plt.figure(i)
            save_str = folder_name + "/figure_" + str(i) + ".png"
            new_fig.savefig(save_str)
    except:
        print("ERROR: Could not save some or all of the figures. Try saving figures by hand.")
        plt.show()
        return

    try:
        for anim_idx in range(len(anims)):
            plt.figure(plotnum + anim_idx)
            anim = anims[anim_idx]
            save_str = folder_name + "/anim_" + str(anim_idx + 1) + ".gif"
            anim.save(save_str)
    except:
        print("ERROR: Could not save some or all of animated figures. Try saving figures by hand.")
        #plt.show()
        return

    print("Figures exported successfully!")
    
    print("Clearing figures from memory...")
    for i in range(1, plotnum + anim_plotnum + 1):
        new_fig = plt.figure(i)
        new_fig.clear()
        plt.close()

    print("Analysis done!")
