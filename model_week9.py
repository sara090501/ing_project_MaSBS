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
import csv
import math
import random

#nova uloha week 10
# sirka steny - lb_grid * 2
# dlzka romboidu msi byt presne 280microm lebo je to rovnoramenny trojuholnik - 60°
# sirka v kolmom priereze (na na hranach) je 80
# sirka na hranach je 160
# if i == 100 - po 100 krokoch zaznamenavat kvapalinu do .vtk podla vzorca z word-u lebo je jej vela
parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/model_week9_" + str(args.sim_no)
os.makedirs(directory)

boxX = 560.0 # vyratanie velkosti kvadru v ktorom bude simulacia
boxY = 40.0
boxZ = 220.0
system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = 0.1

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
cell = oif.OifCell(cell_type=type, particle_type=0, origin=[5.0, 8.0, 6.0])
cell1 = oif.OifCell(cell_type=type, particle_type=1, origin=[5.0, 12.0, 13.0], rotate=[0.0, np.pi/2, 0.0], particle_mass=0.5)

print("Cells created.")

# cell-wall interactions
system.non_bonded_inter[0, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)

# cell-cell interactions
system.non_bonded_inter[0, 1].membrane_collision.set_params(a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

print("Interactions created.")

# fluid
lb_params = {'agrid': 1, 'dens': 1, 'visc': 1.5, 'tau': system.time_step, 'ext_force_density': [0.001, 0.0, 0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

# Konštanty pre kanál
segment_length = 198.0
zigzag_angle = np.radians(60)  # Uhol v radiánoch
sirka_kanala = 80
height_change = segment_length * np.sin(zigzag_angle)  # Zmena výšky segmentu
horizontal_change = segment_length * np.cos(zigzag_angle)  # Horizontálna zložka

# Funkcia na vytváranie segmentov
# Funkcia na vytváranie segmentov pre jednu stenu
def create_wall(start_corner, direction, upward_initial, wall_type):
    corners = [start_corner]
    upward = upward_initial

    for _ in range(4):  # Štyri segmenty na vytvorenie jednej steny
        if upward:
            offset = [horizontal_change, 0.0, height_change]
        else:
            offset = [horizontal_change, 0.0, -height_change]

        corners.append(corners[-1] + np.array(offset))
        upward = not upward  # Striedanie smeru

    if wall_type == "side":  # Bočné steny (ľavá/pravá)
        a_vectors = [corners[i + 1] - corners[i] for i in range(4)]
        b_vector = [0.0, 1.0, 0.0]
        c_vector = [0.0, 0.0, sirka_kanala]
    elif wall_type == "top_bottom":  # Horné/spodné steny
        a_vectors = [corners[i + 1] - corners[i] for i in range(4)]
        b_vector = [0.0, boxY, 0.0]
        c_vector = [0.0, 0.0, 1.0]
    else:
        raise ValueError("Unknown wall type: must be 'side' or 'top_bottom'.")

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
current_corner_right = np.array([0.0, boxY - 1.0, 0.0])  # Pravá stena začína pri pravom okraji
current_corner_left = np.array([0.0, 0.0, 0.0])  # Ľavá stena začína pri ľavom okraji
current_corner_top = np.array([0.0, 0.0, sirka_kanala - 1.0])  # Vrchná stena začína pri vrchu
current_corner_bottom = np.array([0.0, 0.0, 0.0])  # Spodná stena začína pri spodku

# Bočné steny
boundaries["right_wall"] = create_wall(current_corner_right, direction=1, upward_initial=False, wall_type="side")
boundaries["left_wall"] = create_wall(current_corner_left, direction=1, upward_initial=False, wall_type="side")

# Vrchná a spodná stena
boundaries["top_wall"] = create_wall(current_corner_top, direction=1, upward_initial=False, wall_type="top_bottom")
boundaries["bottom_wall"] = create_wall(current_corner_bottom, direction=1, upward_initial=False, wall_type="top_bottom")

# Export stien do VTK a pridanie do systému
for wall_name, wall_segments in boundaries.items():
    for idx, segment in enumerate(wall_segments):
        vtk_filename = f"{directory}/{wall_name}_segment_{idx}.vtk"
        oif.output_vtk_rhomboid(rhom_shape=segment, out_file=vtk_filename)
        system.lbboundaries.add(lbboundaries.LBBoundary(shape=segment))
        system.constraints.add(shape=segment, particle_type=10, penetrable=False)

print("Boundaries created.")

maxCycle = 50
# main integration loop
cell.output_vtk_pos_folded(file_name=directory+"/cell_0.vtk")
cell1.output_vtk_pos_folded(file_name=directory+"/cell1_0.vtk")
for i in range(1, maxCycle):
    system.integrator.run(steps=500)
    cell.output_vtk_pos_folded(file_name=directory+"/cell_" + str(i) + ".vtk")
    cell1.output_vtk_pos_folded(file_name=directory+"/cell1_" + str(i) + ".vtk")
    print("time: ", str(i*system.time_step*500))
    #print(cell1.get_velocity()[0], cell.get_velocity()[0])
    #print("......... ")
print("Simulation completed.")
