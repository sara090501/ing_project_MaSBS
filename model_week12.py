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
import random
import csv

parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/sim" + str(args.sim_no)
os.makedirs(directory)

LBgrid = 2

fluid_speed = 0.3

canal_width_Y = 80.0
canal_width_Z = 40.0
side_wall_width = 4.0

wall_lenght = 280
height_in_turn = 160

boxX = 560.0  # vyratanie velkosti kvadru v ktorom bude simulacia
boxY = int(height_in_turn + (math.sqrt(3) / 2) * wall_lenght) + 2 * side_wall_width - 2
boxZ = canal_width_Z + 2 * side_wall_width

system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = 0.1

print("boxX: ", boxX)
print("boxY: ", boxY)
print("boxZ: ", boxZ)


# Definícia priamok tvoriacich obdĺžnik
def line1(x, y):
    return 40 * x + 40 * math.sqrt(3) * y - 6400 * math.sqrt(3)


def line2(x, y):
    return (140 * math.sqrt(3) - 120) * x + (40 * math.sqrt(3) - 140) * y


def line3(x, y):
    return 40 * x + 40 * math.sqrt(3) * y - 22400


def line4(x, y):
    return (-120 + 140 * math.sqrt(3)) * x + (-140 + 40 * math.sqrt(3)) * y + 22400 - 6400 * math.sqrt(3)


# Kontrola, ci bod (x, y) je vo vnutri obdlznika
def inside_rectangle(x, y):
    return line1(x, y) >= 0 and line2(x, y) <= 0 and line3(x, y) <= 0 and line4(x,
                                                                                y) >= 0  # Martin 19.12.2024 bolo potrebne zmenit znamienka u kazdej podmienky

#Kontrola zapisu origin-u bunky do csv
def lineForOrigin(x, y):
    return 40 * x + 40 * math.sqrt(3) * y - 6400 * math.sqrt(3) - 16800

# Rozmery buniek
small_cell_size = 5  # mm
big_cell_size = 10  # mm
margin = 10  # mm od stien

# creating the template for RBCs
small_cell_type = oif.OifCellType(nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat",
                                  # hodnoty na vyratanie krviniek
                                  check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                                  kag=0.9, kv=0.5, resize=[2.5, 2.5, 2.5], normal=False)
big_cell_type = oif.OifCellType(nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat",
                                # hodnoty na vyratanie krviniek
                                check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                                kag=0.9, kv=0.5, resize=[5.0, 5.0, 5.0], normal=False)
# ks = koeficient natiahnutia
# kb = koeficient ohybnosti
# kal = koeficient local area
# kag = koeficient global area
# kv = koeficient objemu

# creating the RBCs
# cell = oif.OifCell(cell_type=small_cell_type, particle_type=0, origin=[5.0, 8.0, 6.0])
# cell1 = oif.OifCell(cell_type=big_cell_type, particle_type=1, origin=[5.0, 12.0, 13.0])

# Generovanie pozícií buniek
cells = []

# Definovanie gridu
# x_random = random.uniform(0, 140)
# y_random = random.uniform(40 * math.sqrt(3), 140 * math.sqrt(3))

cell = oif.OifCell(cell_type=small_cell_type, particle_type=0, origin=[0.0, 0.0, 0.0])
cells.append(cell)  # Martin 19.12.2024 pridanie buniek do pola cells
cell1 = oif.OifCell(cell_type=big_cell_type, particle_type=1, origin=[0.0, 0.0, 0.0])
cells.append(cell1)

while True:
    x_random = random.uniform(0,
                              140)  # Martin 19.12.2024 random suradnice sa vyberaju z obdlznika moznych suradnic https://ctrlv.link/SdK8 a potom sa este testuje ci su v kanali
    y_random = random.uniform(160 - (40 * math.sqrt(3)), 140 * math.sqrt(3))

    if inside_rectangle(x_random + small_cell_size + margin, y_random + small_cell_size + margin) \
            and inside_rectangle(x_random - small_cell_size - margin, y_random - small_cell_size - margin):
        cell.set_origin([x_random, y_random, 20.0])

        break  # Bod leží vo vnútri obdĺžnika, prerušíme generovanie

while True:
    x_random = random.uniform(0,
                              140)  # Martin 19.12.2024 random suradnice sa vyberaju z obdlznika moznych suradnic https://ctrlv.link/SdK8 a potom sa este testuje ci su v kanali
    y_random = random.uniform(160 - (40 * math.sqrt(3)), 140 * math.sqrt(3))

    if inside_rectangle(x_random + big_cell_size + margin, y_random + big_cell_size + margin) \
            and inside_rectangle(x_random - big_cell_size - margin, y_random - big_cell_size - margin):
        cell1.set_origin([x_random, y_random, 20.0])

        break  # Bod leží vo vnútri obdĺžnika, prerušíme generovanie

# for x in x_random :
#     for y in y_random:
#         if inside_rectangle(x + margin, y + margin):
#             # Cell 1 (5 mm)
#             cell.set_origin([x, y, 20.0])

#         # Generuj pozíciu cell2 s odstupom
#         x2 = x + big_cell_size + margin
#         y2 = y + big_cell_size + margin
#         if inside_rectangle(x2, y2):
#             cell1.set_origin([x2, y2, 20.0])

print("Cells created.")

# cell-wall interactions
system.non_bonded_inter[0, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)

# print("Interactions created.")

# fluid
lb_params = {'agrid': LBgrid, 'dens': 1, 'visc': 1.5, 'tau': system.time_step}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# left velocity wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, height_in_turn + 2 * side_wall_width, 0.0],
                            a=[side_wall_width, 0.0, 0.0], b=[0.0, -height_in_turn - 2 * side_wall_width, 0.0],
                            c=[0.0, 0.0, canal_width_Z + 2 * side_wall_width], direction=1)
system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape, velocity=[
    fluid_speed, 0.0, 0.0]))  # toto velocity = [] je smer ktorym chem zacat
oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/velocityWallLeft.vtk")

# front wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, canal_width_Z + side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0],
                            b=[0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
                            c=[0.0, 0.0, 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/frontWall.vtk")
# back wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, -side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0],
                            b=[0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
                            c=[0.0, 0.0, 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/backWall.vtk")
# -------------------------------------------------------------------------------------------------------------------------------------------------------------


# Konštanty pre kanál
zigzag_angle = np.radians(60)  # Uhol v radiánoch
height_change = wall_lenght * np.sin(zigzag_angle)  # Zmena výšky segmentu
horizontal_change = wall_lenght / 2  # Horizontálna zložka


# Funkcia na vytváranie segmentov
# Funkcia na vytváranie segmentov pre jednu stenu
def create_wall(direction, wall_type):
    corners = []

    for i in range(5):  # 5 4-uholnikov na vytvorenie zig-zag kanala  #zdvojnasobujem velkost 4 uholnikov
        if i < 3:
            offset = [2 * horizontal_change * i, height_in_turn + 4 * height_change, 0.0]
        elif i == 3:
            offset = [horizontal_change * (i - 2), height_change, 0.0]
        elif i == 4:
            offset = [horizontal_change * (i - 1), height_change, 0.0]
        corners.append(np.array(offset))

    if wall_type == "steny" and i < 3:  # 4-uholniky tvoriace kanal #zdvojnasobujem velkost 4 uholnikov
        a_vector = [wall_lenght, 0.0, 0.0]
        b_vector = [-horizontal_change * 2, int(boxY) + height_change, 0.0]
        c_vector = [0.0, 0.0, int(boxZ)]
    elif wall_type == "steny" and i >= 3:
        a_vector = [wall_lenght, -2 * height_change, 0.0]
        b_vector = [-horizontal_change * 2, -2 * height_change, 0.0]
        c_vector = [0.0, 0.0, int(boxZ)]

    # Spojenie piatich tvarov
    walls = []
    for i in range(5):
        walls.append(
            Rhomboid(corner=corners[i], a=a_vector, b=b_vector, c=c_vector, direction=direction)
        )
    return walls


# Generovanie kanála
boundaries = {}

# Bočné steny
boundaries["top_bottom_wall"] = create_wall(direction=1, wall_type="steny")

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
cell.output_vtk_pos_folded(file_name=directory + "/cell_0.vtk")
cell1.output_vtk_pos_folded(file_name=directory + "/cell1_0.vtk")
for i in range(1, maxCycle):
    # if (i % 2 == 0):
    #      lbf.print_vtk_velocity(directory + "/fluid.vtk")
    system.integrator.run(steps=500)
    cell.output_vtk_pos_folded(file_name=directory + "/cell_" + str(i) + ".vtk")
    cell1.output_vtk_pos_folded(file_name=directory + "/cell1_" + str(i) + ".vtk")
    for j in range(0, len(cells)):
        origin = cells[j].get_origin()
        x = origin[0]
        y = origin[1]
        #bunka presla za hranicu priamky a je v druhej polovici geometrie
        if lineForOrigin(x, y) >= 0 and x > boxX / 2:
            cells[j].output_vtk_pos_folded(file_name=directory+"/cell" + str(j) + "_" + str(i) + ".vtk")
            with open(directory + '/cell' + str(j) + 'positions.csv', 'a', newline='') as file:
                writer = csv.writer(file, delimiter=';')
                writer.writerow(origin)
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------Martin 19.12.2024 reseeding
    for ii in range(0, 2):
        origin = cells[ii].get_origin()
        if origin[0] > (1.5 * wall_lenght):
            n_x1 = origin[0] - wall_lenght
            cells[ii].set_origin([n_x1, origin[1], origin[2]])
            print("Reseeding cell: " + str(ii) + "\n")
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("Iteration: " + str(i))
    print("time: ", str(i * system.time_step * 500))
    print("......... ")
    # print(cell1.get_velocity()[0], cell.get_velocity()[0])
    # print("......... ")
print("Simulation completed.")