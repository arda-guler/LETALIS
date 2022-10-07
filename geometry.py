import math
pi = math.pi

from material import SS304L, CuCrZr

def get_area_of_sector(r_in, r_out, angle):
    return pi * (r_out**2 - r_in**2) * (angle/360)

def get_volume_of_sector(r_in, r_out, angle, height):
    return get_area_of_sector(r_in, r_out, angle) * height

def get_sector_face_area(r, angle, height):
    return (2*pi*r)*height * (angle/360)

def get_sector_face_area_slanted(r1, r2, angle, height):
    return (pi*r2 + pi*r1) * height * (angle/360)

class cylinder:
    def __init__(self, x, r_in, r_out, h, n_clt, r_clt, a_clt, mtl, T_init, M, r_prev=None):
        self.x = x
        self.r_in = r_in
        self.thickness = r_out - r_in
        self.r_out = r_out
        self.h = h # height

        # inner wall length
        if r_prev and not r_prev == r_in: # cylinder is on a curve
            self.a_axial = math.atan((r_in - r_prev)/h)
            self.L_in = self.h * math.cos(self.a_axial)
        else: # cylinder has constant radius
            self.a_axial = 0
            self.L_in = self.h

        self.n_clt = n_clt # number of coolant channels
        self.r_clt = r_clt # coolant channel inner radius
        self.a_clt = a_clt # coolant channel angle (degrees)
        
        self.mtl = mtl # material
        self.T = T_init # temperature
        self.Mach = M # flow mach number at cylinder position

        if r_prev:
            dr = r_clt - r_in
            r_clt2 = r_prev + dr
            self.A_chm = get_sector_face_area_slanted(r_in, r_prev, 360, h) # area facing chamber (m2)
            self.A_cochan_flow = get_area_of_sector(r_clt, r_out, a_clt) # flow area of single coolant channel (m2)
            self.A_clt = (get_sector_face_area_slanted(r_clt, r_clt2, a_clt, h) + 2*(r_out - r_clt)*h) * n_clt # area facing coolant (m2)

        else:
            self.A_chm = get_sector_face_area(r_in, 360, h) # area facing chamber (m2)
            self.A_cochan_flow = get_area_of_sector(r_clt, r_out, a_clt) # flow area of single coolant channel (m2)
            self.A_clt = (get_sector_face_area(r_clt, a_clt, h) + 2*(r_out - r_clt)*h) * n_clt # area facing coolant (m2)

        self.V = get_volume_of_sector(r_in, r_out, 360, h) - (n_clt * get_volume_of_sector(r_clt, r_out, a_clt, h)) # volume (m3)
        self.m = self.V * self.mtl.get_density() # mass

    def get_m(self):
        return self.m

    def get_Mach(self):
        return self.Mach

    def get_mtl(self):
        return self.mtl

    def get_T(self):
        return self.temp

    def get_A_chm(self):
        return self.A_chm

    def get_A_clt(self):
        return self.A_clt

    def get_A_cochan_flow(self):
        return self.A_cochan_flow

    def get_spec_heat(self):
        return self.mtl.get_specific_heat(self.T)

    def get_heat_capacity(self):
        return self.get_spec_heat() * self.get_m()

# calculate_geometry() calculates the whole geometry all at once and show it
# to the user so that they can see if there are any problems with the mathematical model.
# also it generates a list that can be used as a look-up table for radius at various
# x coordinates should a function ever need it
#
# but why don't I just use get_inner_radius_at() in a while loop to generate the function?
# that's because I don't want get_inner_radius_at() to return L1 - L7 values. maybe there
# is a better way to do that but I'm not a computer engineer and it does the job as well as
# any other method.
def calculate_geometry(L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm,
                       a_nzlExp, ROC_thrtDn, ROC_thrtUp, fineness):
    global pi

    x_step = L_engine/fineness

    def arc_diff(R, theta):
        
        if theta > 2*pi:
            print("Func: arc_diff() -- Make sure theta is in radians!")
            quit()
            
        # takes R in meters, theta in radians
        # returns radius difference in meters
        return R * (1 - math.cos(theta))
    
    # R and theta are used as temporary variables throughout this fucntion

    R_thrt = D_thrt/2 # m
    R_chm = D_chm/2 # m
    
    # CHAMBER CONTRACTION CURVE
    
    R = ROC_chm # temporary variable to keep the code short (m)
    theta = math.radians(a_chmContract) # temporary variable to keep the code short (radians)
    L_chamber_contract_curve = R * math.sin(theta) # m
    R_chamber_contract_cone_max = R_chm - arc_diff(R, theta) # m

    R = ROC_thrtUp # m
    theta = math.radians(a_chmContract) # radians
    R_chamber_contract_cone_min = R_thrt + arc_diff(R, theta) # m
    L_nzl_upstream_curve = R * math.sin(theta)

    L_chm_contract_cone = (R_chamber_contract_cone_max - R_chamber_contract_cone_min) * (1/math.tan(theta)) ## change
    L_chm_contract_total = L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve
    
    # NOZZLE DOWNSTREAM CURVE
    
    L_nzl_downstream_curve = ROC_thrtDn * math.tan(math.radians(a_nzlExp)) # m
    L_nzl_upstream_curve = ROC_thrtUp * math.sin(math.radians(a_chmContract)) # m

    R = ROC_thrtDn # temporary variable to keep the code short (m)
    theta = math.radians(a_nzlExp) # temporary variable to keep the code short (radians)
    R_nzl_downstream_curve_max = arc_diff(R, theta) + R_thrt # m
    L_nzl_downstream_curve = R * math.sin(theta) # m
    R_exit = D_exit/2 # m
    L_nzl_cone = (R_exit - R_nzl_downstream_curve_max) * (1/math.tan(theta)) # m
    L_nzl = L_nzl_downstream_curve + L_nzl_cone # m

    # COMBUSTION CHAMBER CYLINDER
    L_chm_cylinder = L_engine - (L_chm_contract_total + L_nzl)

    if L_chm_cylinder <= 0:
        print("Invalid geometry! No room for cylindrical combustion chamber segment!")
        quit()

    # show calculated engine geometry
    L1 = 0
    L2 = L_chm_cylinder
    L3 = L_chm_cylinder + L_chamber_contract_curve
    L4 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone
    L5 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve
    L6 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve + L_nzl_downstream_curve
    L7 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve + L_nzl_downstream_curve + L_nzl_cone

    x1 = []
    y1 = []

    x2 = []
    y2 = []

    x3 = []
    y3 = []

    x4 = []
    y4 = []

    x5 = []
    y5 = []

    x6 = []
    y6 = []

    i = 0
    while i < L_engine:
        x_current = i

        # cylindrical combustion chamber segment
        if L1 <= i <= L2:
            y_current = R_chm
            x1.append(x_current)
            y1.append(y_current)

        # combustion chamber curving inwards
        elif L2 < i <= L3:
            theta = math.asin((i - L2) / ROC_chm)
            y_current = R_chm - arc_diff(ROC_chm, theta)
            x2.append(x_current)
            y2.append(y_current)

        # combustion chamber contraction with constant angle
        elif L3 < i <= L4:
            y_current = (R_chamber_contract_cone_max - R_chamber_contract_cone_min) * ((L4 - i)/(L4 - L3)) + R_chamber_contract_cone_min
            x3.append(x_current)
            y3.append(y_current)

        # nozzle upstream curve
        elif L4 < i <= L5:
            theta = math.asin(((L5 - i)/ROC_thrtUp))
            y_current = R_thrt + arc_diff(ROC_thrtUp, theta)
            x4.append(x_current)
            y4.append(y_current)

        # nozzle downstream curve
        elif L5 < i <= L6:
            theta = math.asin((i - L5)/ROC_thrtDn)
            y_current = R_thrt + arc_diff(ROC_thrtDn, theta)
            x5.append(x_current)
            y5.append(y_current)

        # nozzle cone
        elif L6 < i <= L7:
            y_current = R_nzl_downstream_curve_max + (R_exit - R_nzl_downstream_curve_max) * ((i - L6)/L_nzl_cone)
            x6.append(x_current)
            y6.append(y_current)

        i += L_engine/1000

    geom_x = x1 + x2 + x3 + x4 + x5 + x6
    geom_y = y1 + y2 + y3 + y4 + y5 + y6

    ## This section can be uncommented if you need to debug geometry creation formulas or such
    ## plt.plot(x1, y1)
    ## plt.plot(x2, y2)
    ## plt.plot(x3, y3)
    ## plt.plot(x4, y4)
    ## plt.plot(x5, y5)
    ## plt.plot(x6, y6)
    ## plt.show()

    return geom_x, geom_y, x_step, [L1, L2, L3, L4, L5, L6, L7]

# The function below is basically a copy of the calculate_geometry() function but
# it only calculates the radius at one point, which is required for other mathematical
# functions to work
#
# you could technically generate engine geometry once and just read the radius values
# from the generated list when needed, but just having the mathematical function to get
# the radius at any arbitrary point is more desirable, because then you don't need to do
# interpolations if you choose a value that's in the middle of two values on the pre-generated
# plot for radius vs. x coordinates
#
# this is basically the more correct way of doing it
def get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp):
    global pi

    def arc_diff(R, theta):
        
        if theta > 2*pi:
            print("Func: arc_diff() -- Make sure theta is in radians!")
            quit()
            
        # takes R in meters, theta in radians
        # returns radius difference in meters
        return R * (1 - math.cos(theta))
    
    # R and theta are used as temporary variables throughout this fucntion

    R_thrt = D_thrt/2 # m
    R_chm = D_chm/2 # m
    
    # CHAMBER CONTRACTION CURVE
    
    R = ROC_chm # temporary variable to keep the code short (m)
    theta = math.radians(a_chmContract) # temporary variable to keep the code short (radians)
    L_chamber_contract_curve = R * math.sin(theta) # m
    R_chamber_contract_cone_max = R_chm - arc_diff(R, theta) # m

    R = ROC_thrtUp # m
    theta = math.radians(a_chmContract) # radians
    R_chamber_contract_cone_min = R_thrt + arc_diff(R, theta) # m
    L_nzl_upstream_curve = R * math.sin(theta)

    L_chm_contract_cone = (R_chamber_contract_cone_max - R_chamber_contract_cone_min) * (math.tan(theta))
    L_chm_contract_total = L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve
    
    # NOZZLE DOWNSTREAM CURVE
    
    L_nzl_downstream_curve = ROC_thrtDn * math.tan(math.radians(a_nzlExp)) # m
    L_nzl_upstream_curve = ROC_thrtUp * math.sin(math.radians(a_chmContract)) # m

    R = ROC_thrtDn # temporary variable to keep the code short (m)
    theta = math.radians(a_nzlExp) # temporary variable to keep the code short (radians)
    R_nzl_downstream_curve_max = arc_diff(R, theta) + R_thrt # m
    L_nzl_downstream_curve = R * math.sin(theta) # m
    R_exit = D_exit/2 # m
    L_nzl_cone = (R_exit - R_nzl_downstream_curve_max) * (1/math.tan(theta)) # m
    L_nzl = L_nzl_downstream_curve + L_nzl_cone # m

    # COMBUSTION CHAMBER CYLINDER
    L_chm_cylinder = L_engine - (L_chm_contract_total + L_nzl)

    if L_chm_cylinder <= 0:
        print("Invalid geometry! No room for cylindrical combustion chamber segment!")
        quit()

    # show calculated engine geometry
    L1 = 0 # injector mount
    L2 = L_chm_cylinder
    L3 = L_chm_cylinder + L_chamber_contract_curve
    L4 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone
    L5 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve # throat
    L6 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve + L_nzl_downstream_curve
    L7 = L_chm_cylinder + L_chamber_contract_curve + L_chm_contract_cone + L_nzl_upstream_curve + L_nzl_downstream_curve + L_nzl_cone # exit

    i = x

    # cylindrical combustion chamber segment
    if L1 <= i <= L2:
        y_current = R_chm

    # combustion chamber curving inwards
    elif L2 < i <= L3:
        theta = math.asin((i - L2) / ROC_chm)
        y_current = R_chm - arc_diff(ROC_chm, theta)

    # combustion chamber contraction with constant angle
    elif L3 < i <= L4:
        y_current = (R_chamber_contract_cone_max - R_chamber_contract_cone_min) * ((L4 - i)/(L4 - L3)) + R_chamber_contract_cone_min

    # nozzle upstream curve
    elif L4 < i <= L5:
        theta = math.asin(((L5 - i)/ROC_thrtUp))
        y_current = R_thrt + arc_diff(ROC_thrtUp, theta)

    # nozzle downstream curve
    elif L5 < i <= L6:
        theta = math.asin((i - L5)/ROC_thrtDn)
        y_current = R_thrt + arc_diff(ROC_thrtDn, theta)

    # nozzle cone
    elif L6 < i <= L7:
        y_current = R_nzl_downstream_curve_max + (R_exit - R_nzl_downstream_curve_max) * ((i - L6)/L_nzl_cone)

    else:
        print("ERROR:", x, "is out of bounds to get radius!")
        quit()

    return y_current # meters

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
def get_mach_num_at(x):
    global subsonic_mach, subsonic_x, supersonic_mach, supersonic_x, engine_lengths

    # subsonic region
    if x < engine_lengths[4]:
        index = get_index_of_closest_num_in_list(x, subsonic_x)
        return subsonic_mach[index]

    # supersonic region
    elif engine_lengths[4] <= x <= engine_lengths[6]:
        index = get_index_of_closest_num_in_list(x, supersonic_x)
        return supersonic_mach[index]

# generate a 3D model of the engine geometry
def generate_3D(geom_x, geom_y, n_cochan, L_cochanInnerWallDist, L_cochanSideWall, L_cochanDepth):
    
    vertices = []
    faces = []
    n_vertices_half_circle = max(int(n_cochan * 5), 90)

    x_index = 0
    for x in geom_x:

        r_in = geom_y[x_index]
        theta = -pi/2

        # build inner wall
        for i in range(n_vertices_half_circle):
            current_vertex = [x, r_in * math.sin(theta), r_in * math.cos(theta)]
            vertices.append(current_vertex)
            theta += pi/n_vertices_half_circle

        vertices.append("OUTERSHELL")
        
        theta = -pi/2
        n_halfCochan = n_cochan/2
        r_clt = r_in + L_cochanInnerWallDist
        r_out = r_clt + L_cochanDepth
        a_sideWall = L_cochanSideWall / r_out
        a_channel = (pi - (a_sideWall * n_halfCochan)) / n_halfCochan
        a_full = a_sideWall + a_channel

        # build cooling channels
        for i in range(n_vertices_half_circle):
            # it is a wall vertex
            if theta % a_full <= a_sideWall:
                current_vertex = [x, r_out * math.sin(theta), r_out * math.cos(theta)]
                
            # it is a channel vertex
            else:
                current_vertex = [x, r_clt * math.sin(theta), r_clt * math.cos(theta)]

            vertices.append(current_vertex)
            theta += pi/n_vertices_half_circle

        vertices.append("NEWX")

        x_index += 1

    return vertices
