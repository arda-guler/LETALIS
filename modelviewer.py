import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *
import glfw

import keyboard
from pyquaternion import Quaternion
import numpy
import math

def show_scale(origin, scale, end):
    current_x = origin[0]
    current_y = origin[1]
    current_z = origin[2]

    glPushMatrix()
    glColor(0,1,0)
    glBegin(GL_LINES)
    while current_x < end[0]:
        glVertex3f(current_x, current_y, current_z)
        current_x += scale

    glEnd()
    glPopMatrix()
    current_x = origin[0]

    glPushMatrix()
    glBegin(GL_LINES)
    while current_y < end[1]:
        glVertex3f(current_x, current_y, current_z)
        current_y += scale

    glEnd()
    glPopMatrix()
    current_y = origin[1]

    glPushMatrix()
    glBegin(GL_LINES)
    while current_z < end[2]:
        glVertex3f(current_x, current_y, current_z)
        current_z += scale

    glEnd()
    glPopMatrix()
    current_z = origin[2]

# rotate an orientation matrix
def rotate_matrix(orientation_matrix, rotation):
    # orientation matrix is a 3x3 matrix, rotation is a list of three angles in degrees
    orientation_matrix = numpy.array(orientation_matrix)
        
    if rotation[0]:
        rotator = Quaternion(axis=orientation_matrix[0], angle=math.radians(rotation[0]))
        orientation_matrix = (numpy.array([rotator.rotate(orientation_matrix[0]), rotator.rotate(orientation_matrix[1]), rotator.rotate(orientation_matrix[2])]))

    if rotation[1]:
        rotator = Quaternion(axis=orientation_matrix[1], angle=math.radians(rotation[1]))
        orientation_matrix = (numpy.array([rotator.rotate(orientation_matrix[0]), rotator.rotate(orientation_matrix[1]), rotator.rotate(orientation_matrix[2])]))

    if rotation[2]:
        rotator = Quaternion(axis=orientation_matrix[2], angle=math.radians(rotation[2]))
        orientation_matrix = (numpy.array([rotator.rotate(orientation_matrix[0]), rotator.rotate(orientation_matrix[1]), rotator.rotate(orientation_matrix[2])]))

    return orientation_matrix.tolist()

class camera:
    def __init__(self, pos, orient):
        self.pos = pos
        self.orient = orient

    def move(self, movement):
        glTranslate((movement[0] * self.orient[0][0]) + (movement[1] * self.orient[1][0]) + (movement[2] * self.orient[2][0]),
                    (movement[0] * self.orient[0][1]) + (movement[1] * self.orient[1][1]) + (movement[2] * self.orient[2][1]),
                    (movement[0] * self.orient[0][2]) + (movement[1] * self.orient[1][2]) + (movement[2] * self.orient[2][2]))

        self.pos = [self.pos[0] + (movement[0] * self.orient[0][0]) + (movement[1] * self.orient[1][0]) + (movement[2] * self.orient[2][0]),
                    self.pos[1] + (movement[0] * self.orient[0][1]) + (movement[1] * self.orient[1][1]) + (movement[2] * self.orient[2][1]),
                    self.pos[2] + (movement[0] * self.orient[0][2]) + (movement[1] * self.orient[1][2]) + (movement[2] * self.orient[2][2])]

    def rotate(self, rotation):
        about_pos = self.pos
        
        glTranslate(-about_pos[0], -about_pos[1], -about_pos[2])
        glRotate(-rotation[0], self.orient[0][0], self.orient[0][1], self.orient[0][2])
        glTranslate(about_pos[0], about_pos[1], about_pos[2])

        glTranslate(-about_pos[0], -about_pos[1], -about_pos[2])
        glRotate(-rotation[1], self.orient[1][0], self.orient[1][1], self.orient[1][2])
        glTranslate(about_pos[0], about_pos[1], about_pos[2])

        glTranslate(-about_pos[0], -about_pos[1], -about_pos[2])
        glRotate(-rotation[2], self.orient[2][0], self.orient[2][1], self.orient[2][2])
        glTranslate(about_pos[0], about_pos[1], about_pos[2])

        self.orient = rotate_matrix(self.orient, rotation)

def import_model():
    import_file = open("3d_model.txt", "r")
    import_lines = import_file.readlines()

    axials = [ [ [],[] ] ]
    x_index = 0
    outer = 0

    for line in import_lines:
        
        if line == "OUTERSHELL\n":
            outer = 1

        elif line == "NEWX\n":
            axials.append([[],[]])
            x_index += 1
            outer = 0
            
        else:
            line = line[1:-2]
            line = line.split(",")
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            axials[x_index][outer].append([x, y, z])

    return axials

def main():
    print("Reading 3D model...")
    model_data = import_model()

    window_x = 1000
    window_y = 600
    window_ratio = window_x/window_y
    fov = 70
    near_clip = 0.0005
    far_clip = 100

    glfw.init()

    window = glfw.create_window(window_x, window_y, "Engine Viewer", None, None)
    glfw.set_window_pos(window,100,100)
    glfw.make_context_current(window)
    
    gluPerspective(fov, window_ratio, near_clip, far_clip)
    glEnable(GL_CULL_FACE)
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT)

    main_cam = camera([0,0,0], [[1,0,0],[0,1,0],[0,0,1]])

    show_inner = True
    show_outer = True

    x_end = model_data[-2][0][0][0]
    x_mid = x_end/2

    y_end = model_data[0][0][0][1]
    z_end = y_end
    main_cam.move([-x_mid, 0, y_end*5])

    n_angular = len(model_data[0][0])
    
    scale_size = 0.01 #m
    print("Reference scale:", scale_size*1000, "mm")
    scale_origin = [0, y_end, z_end]
    scale_end = [scale_origin[0] + x_end, scale_origin[1]-2*y_end, scale_origin[2]-2*z_end]
    
    while not glfw.window_should_close(window):
        glfw.poll_events()

        # camera controls
        if keyboard.is_pressed("shift"):
            if keyboard.is_pressed("w"):
                main_cam.move([0, 0, 0.004])
            if keyboard.is_pressed("s"):
                main_cam.move([0, 0, -0.004])
            if keyboard.is_pressed("d"):
                main_cam.move([-0.004, 0, 0])
            if keyboard.is_pressed("a"):
                main_cam.move([0.004, 0, 0])
            if keyboard.is_pressed("f"):
                main_cam.move([0, 0.004, 0])
            if keyboard.is_pressed("r"):
                main_cam.move([0, -0.004, 0])

            if keyboard.is_pressed("i"):
                main_cam.rotate([5,0,0])
            if keyboard.is_pressed("k"):
                main_cam.rotate([-5,0,0])
            if keyboard.is_pressed("j"):
                main_cam.rotate([0,5,0])
            if keyboard.is_pressed("l"):
                main_cam.rotate([0,-5,0])
            if keyboard.is_pressed("u"):
                main_cam.rotate([0,0,5])
            if keyboard.is_pressed("o"):
                main_cam.rotate([0,0,-5])
        else:
            if keyboard.is_pressed("w"):
                main_cam.move([0, 0, 0.001])
            if keyboard.is_pressed("s"):
                main_cam.move([0, 0, -0.001])
            if keyboard.is_pressed("d"):
                main_cam.move([-0.001, 0, 0])
            if keyboard.is_pressed("a"):
                main_cam.move([0.001, 0, 0])
            if keyboard.is_pressed("f"):
                main_cam.move([0, 0.001, 0])
            if keyboard.is_pressed("r"):
                main_cam.move([0, -0.001, 0])

            if keyboard.is_pressed("i"):
                main_cam.rotate([2,0,0])
            if keyboard.is_pressed("k"):
                main_cam.rotate([-2,0,0])
            if keyboard.is_pressed("j"):
                main_cam.rotate([0,2,0])
            if keyboard.is_pressed("l"):
                main_cam.rotate([0,-2,0])
            if keyboard.is_pressed("u"):
                main_cam.rotate([0,0,2])
            if keyboard.is_pressed("o"):
                main_cam.rotate([0,0,-2])

        # scale controls
        if keyboard.is_pressed("1") and not scale_size == 0.001:
            scale_size = 0.001
            print("Reference scale:", scale_size*1000, "mm")
        elif keyboard.is_pressed("2") and not scale_size == 0.01:
            scale_size = 0.01
            print("Reference scale:", scale_size*1000, "mm")
        elif keyboard.is_pressed("3") and not scale_size == 0.1:
            scale_size = 0.1
            print("Reference scale:", scale_size*1000, "mm")

        if keyboard.is_pressed("shift"):
            if keyboard.is_pressed("g"):
                scale_origin[2] += 0.005
                scale_end[2] += 0.005
            if keyboard.is_pressed("t"):
                scale_origin[2] -= 0.005
                scale_end[2] -= 0.005
            if keyboard.is_pressed("n"):
                scale_origin[0] += 0.005
                scale_end[0] += 0.005
            if keyboard.is_pressed("b"):
                scale_origin[0] -= 0.005
                scale_end[0] -= 0.005
            if keyboard.is_pressed("y"):
                scale_origin[1] += 0.005
                scale_end[1] += 0.005
            if keyboard.is_pressed("h"):
                scale_origin[1] -= 0.005
                scale_end[1] -= 0.005
        else:
            if keyboard.is_pressed("g"):
                scale_origin[2] += 0.001
                scale_end[2] += 0.001
            if keyboard.is_pressed("t"):
                scale_origin[2] -= 0.001
                scale_end[2] -= 0.001
            if keyboard.is_pressed("n"):
                scale_origin[0] += 0.001
                scale_end[0] += 0.001
            if keyboard.is_pressed("b"):
                scale_origin[0] -= 0.001
                scale_end[0] -= 0.001
            if keyboard.is_pressed("y"):
                scale_origin[1] += 0.001
                scale_end[1] += 0.001
            if keyboard.is_pressed("h"):
                scale_origin[1] -= 0.001
                scale_end[1] -= 0.001

        if keyboard.is_pressed("z"):
            show_inner = False
        if keyboard.is_pressed("x"):
            show_inner = True
        if keyboard.is_pressed("c"):
            show_outer = False
        if keyboard.is_pressed("v"):
            show_outer = True
            
        # rendering
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        
        for axial in model_data[::10]:

            # radial draw
            if show_inner:
                glColor(0.8, 0.2, 0.2)
                # inner wall vertices
                glPushMatrix()
                glBegin(GL_LINE_STRIP)
                for ivertex in axial[0][::15]:
                    glVertex3f(ivertex[0], ivertex[1], ivertex[2])
                glEnd()
                glPopMatrix()

            if show_outer:
                glColor(0.2, 0.5, 0.8)
                # outer wall vertices
                glPushMatrix()
                glBegin(GL_LINE_STRIP)
                for overtex in axial[1]:
                    glVertex3d(overtex[0], overtex[1], overtex[2])
                glEnd()
                glPopMatrix()

        rad_index = 0
        
        # axial draw
        glPushMatrix()
        while rad_index < n_angular:

            if show_inner:
                glColor(0.8, 0.2, 0.2)
                glBegin(GL_LINE_STRIP)
                ax_index = 0
                
                while ax_index < len(model_data) - 1:
                    glVertex3f(model_data[ax_index][0][rad_index][0], model_data[ax_index][0][rad_index][1], model_data[ax_index][0][rad_index][2])
                    ax_index += int(len(model_data)/50)
                    
                glEnd()

            if show_outer:
                glColor(0.2, 0.5, 0.8)
                glBegin(GL_LINE_STRIP)
                ax_index = 0
                
                while ax_index < len(model_data) - 1:
                    glVertex3f(model_data[ax_index][1][rad_index][0], model_data[ax_index][1][rad_index][1], model_data[ax_index][1][rad_index][2])
                    ax_index += int(len(model_data)/50)
                    
                glEnd()
            
            rad_index += int(n_angular/50)
        
        glPopMatrix()

        show_scale(scale_origin, scale_size, scale_end)
            
        glfw.swap_buffers(window)

    glfw.terminate()

main()
