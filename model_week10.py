import espressomd
from espressomd import lb

from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions
from espressomd.shapes import Rhomboid

import object_in_fluid as oif
import os, glob, sys, shutil
import numpy as np
import argparse
import math

# lb_grid = 4
parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/sim"+str(args.sim_no)
os.makedirs(directory)

fluid_speed = 0.6

canal_width_Y = 80.0
canal_width_Z = 40.0
side_wall_width = 4.0 #lb_grid

wall_lenght = 280
height_in_turn = 160

boxX = 560.0  # vyratanie velkosti kvadru v ktorom bude simulacia
boxY = int(height_in_turn + (math.sqrt(3)/2) * wall_lenght) + 2*side_wall_width
boxZ = canal_width_Z + 2*side_wall_width

system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = 0.1

print("boxX: ", boxX)
print("boxY: ", boxY)
print("boxZ: ", boxZ)

# creating the template for RBCs
type = oif.OifCellType(nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat", # hodnoty na vyratanie krviniek
                        check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                        kag=0.9, kv=0.5, resize=[4.0, 4.0, 4.0], normal=False)
# ks = koeficient natiahnutia
# kb = koeficient ohybnosti
# kal = koeficient local area
# kag = koeficient global area
# kv = koeficient objemu

# creating the RBCs
# cell = oif.OifCell(cell_type=type, particle_type=0, origin=[5.0, 8.0, 6.0])
# cell1 = oif.OifCell(cell_type=type, particle_type=1, origin=[5.0, 12.0, 13.0], rotate=[0.0, np.pi/2, 0.0], particle_mass=0.5)

#print("Cells created.")

# cell-wall interactions
#system.non_bonded_inter[0, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)
#system.non_bonded_inter[1, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)

# cell-cell interactions
#system.non_bonded_inter[0, 1].membrane_collision.set_params(a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

#print("Interactions created.")

# fluid
lb_params = {'agrid': 1, 'dens': 1, 'visc': 1.5, 'tau': system.time_step, 'ext_force_density': [0.001, 0.0, 0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------

# left velocity wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width/2, height_in_turn + 2*side_wall_width, 0.0], a=[side_wall_width, 0.0, 0.0], b=[
                            0.0, -height_in_turn-2*side_wall_width, 0.0], c=[0.0, 0.0, canal_width_Z + 2*side_wall_width], direction=1)
system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape, velocity=[
                        0.0, -fluid_speed, 0.0]))  # toto velocity = [] je smer ktorym chem zacat
# todo: velocity=[fluid_speed, 0.0, 0.0]
oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/velocityWallLeft.vtk")

# right velocity wall
# todo: zrusit - right velocity wall
tmp_shape = shapes.Rhomboid(corner=[boxX-side_wall_width,height_in_turn + 2*side_wall_width,0.0], a=[side_wall_width, 0.0, 0.0], b=[
                                 0.0, -height_in_turn-2*side_wall_width, 0.0], c=[0.0, 0.0, canal_width_Z + 2*side_wall_width], direction=1)
system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape, velocity=[
                        0.0, -fluid_speed, 0.0]))  # toto velocity = [] je smer ktorym chem zacat
# todo: velocity=[fluid_speed, 0.0, 0.0]
oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/velocityWallRight.vtk")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------


# Konštanty pre kanál
zigzag_angle = np.radians(60)  # Uhol v radiánoch
height_change = wall_lenght * np.sin(zigzag_angle)  # Zmena výšky segmentu
horizontal_change = wall_lenght / 2 # Horizontálna zložka

# Funkcia na vytváranie segmentov
# Funkcia na vytváranie segmentov pre jednu stenu
def create_wall(start_corner, direction, upward_initial, wall_type):
    corners = [start_corner]
    upward = upward_initial

    for _ in range(4):  # Štyri segmenty na vytvorenie jednej steny
        if upward:
            offset = [horizontal_change, -height_change, 0.0]
        else:
            offset = [horizontal_change, height_change, 0.0]

        corners.append(corners[-1] + np.array(offset))
        upward = not upward  # Striedanie smeru

    if wall_type == "top_bottom":  # Bočné steny (ľavá/pravá) ... horne/zadne steny
        a_vectors = [corners[i + 1] - corners[i] for i in range(4)]
        b_vector = [0.0, side_wall_width, 0.0]
        c_vector = [0.0, 0.0, boxZ]
    elif wall_type == "front_back":  # Horné/spodné steny ... predna/zadna stena
        a_vectors = [corners[i + 1] - corners[i] for i in range(4)]
        b_vector = [0.0, height_in_turn + side_wall_width, 0.0]
        c_vector = [0.0, 0.0, side_wall_width]
    else:
        raise ValueError("Unknown wall type: must be 'top_bottom' or 'front_back'.")

    # Spojenie štyroch segmentov do jednej Rhomboid
    walls = []
    for i in range(4):
        walls.append(
            Rhomboid(corner=corners[i], a=a_vectors[i], b=b_vector, c=c_vector, direction=direction)
        )
    return walls


# Generovanie kanála
boundaries = {}

# Generovanie stien
current_corner_top = np.array([0.0, height_in_turn + side_wall_width, 0.0])  # Pravá stena začína pri pravom okraji ... horna stena začína pri hornom okraji
current_corner_bottom = np.array([0.0, 0.0, 0.0])  # Ľavá stena začína pri ľavom okraji ... zadna stena začína pri zadnom okraji
current_corner_front = np.array([0.0, 0.0, boxZ - side_wall_width])  # Vrchná stena začína pri vrchu ... predna stena začína pri predku
current_corner_back = np.array([0.0, 0.0, 0.0])  # Spodná stena začína pri spodku ... zadna stena začína pri zadku

# Bočné steny
boundaries["top_wall"] = create_wall(current_corner_top, direction=1, upward_initial=False, wall_type="top_bottom")
boundaries["bottom_wall"] = create_wall(current_corner_bottom, direction=1, upward_initial=False, wall_type="top_bottom")

# Vrchná a spodná stena
boundaries["front_wall"] = create_wall(current_corner_front, direction=1, upward_initial=False, wall_type="front_back")
boundaries["back_wall"] = create_wall(current_corner_back, direction=1, upward_initial=False, wall_type="front_back")

# Export stien do VTK a pridanie do systému
for wall_name, wall_segments in boundaries.items():
    for idx, segment in enumerate(wall_segments):
        vtk_filename = f"{directory}/{wall_name}_segment_{idx}.vtk"
        oif.output_vtk_rhomboid(rhom_shape=segment, out_file=vtk_filename)
        system.lbboundaries.add(lbboundaries.LBBoundary(shape=segment))
        system.constraints.add(shape=segment, particle_type=10, penetrable=False)

print("Boundaries created.")

maxCycle = 500
# main integration loop
#cell.output_vtk_pos_folded(file_name=directory+"/cell_0.vtk")
#cell1.output_vtk_pos_folded(file_name=directory+"/cell1_0.vtk")
for i in range(1, maxCycle):
    if (i % 25 == 0):
        lbf.print_vtk_velocity(directory + "/fluid.vtk")
    system.integrator.run(steps=500)
    #cell.output_vtk_pos_folded(file_name=directory+"/cell_" + str(i) + ".vtk")
    #cell1.output_vtk_pos_folded(file_name=directory+"/cell1_" + str(i) + ".vtk")
    print("time: ", str(i*system.time_step*500))
    #print(cell1.get_velocity()[0], cell.get_velocity()[0])
    #print("......... ")
print("Simulation completed.")
