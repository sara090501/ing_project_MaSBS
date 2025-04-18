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

parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/sim" + str(args.sim_no)
os.makedirs(directory)

numberOfSmallCells = 5
numberOfBigCells = 5

LBgrid = 1

fluid_speed = 0.6
#0.3 - polovicna oproti prierezu

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

# creating the template for RBCs
# zmenit typy buniek - solid spheres
typeCellBig = oif.OifCellType(nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat",
                       # hodnoty na vyratanie krviniek
                       check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                       kag=0.9, kv=0.5, resize=[10.0, 10.0, 10.0], normal=False)
typeCellSmall = oif.OifCellType(nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat",
                       # hodnoty na vyratanie krviniek
                       check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                       kag=0.9, kv=0.5, resize=[5.0, 5.0, 5.0], normal=False)
# ks = koeficient natiahnutia
# kb = koeficient ohybnosti
# kal = koeficient local area
# kag = koeficient global area
# kv = koeficient objemu

# creating the RBCs
cell = oif.OifCell(cell_type=typeCellSmall, particle_type=0, origin=[(boxX/8)-5, boxY/2, boxZ/2])
cell1 = oif.OifCell(cell_type=typeCellBig, particle_type=1, origin=[(boxX/8)+15, boxY/2, boxZ/2])

print("Cells created.")

# cell-wall interactions
system.non_bonded_inter[0, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)

print("Interactions created.")

# fluid
lb_params = {'agrid': LBgrid, 'dens': 1, 'visc': 1.5, 'tau': system.time_step}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# left velocity wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, height_in_turn + 2 * side_wall_width, 0.0],
                            a=[side_wall_width, 0.0, 0.0], b=[
        0.0, -height_in_turn - 2 * side_wall_width, 0.0], c=[0.0, 0.0, canal_width_Z + 2 * side_wall_width],
                            direction=1)
system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape, velocity=[
    fluid_speed, 0.0, 0.0]))  # toto velocity = [] je smer ktorym chem zacat
oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/velocityWallLeft.vtk")

# top safe wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, height_in_turn + (math.sqrt(3) / 2) * wall_lenght, 0.0],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0], b=[
        0.0, 2 * side_wall_width, 0.0], c=[0.0, 0.0, canal_width_Z + 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/safeWallTop.vtk")
# bottom safe wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, 0.0, 0.0], a=[boxX + side_wall_width / 2, 0.0, 0.0], b=[
    0.0, -2 * side_wall_width, 0.0], c=[0.0, 0.0, canal_width_Z + 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/safeWallBottom.vtk")

# front wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, canal_width_Z + side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0], b=[
        0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
                            c=[0.0, 0.0, 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/frontWall.vtk")
# back wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, -side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0], b=[
        0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
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

    for i in range(5):  # 5 4-uholnikov na vytvorenie zig-zag kanala
        if i < 3:
            offset = [2 * horizontal_change * i, height_in_turn + 2 * height_change, 0.0]
        elif i == 3:
            offset = [horizontal_change * (i - 2), height_change, 0.0]
        elif i == 4:
            offset = [horizontal_change * (i - 1), height_change, 0.0]
        corners.append(np.array(offset))

    if wall_type == "steny" and i < 3:  # 4-uholniky tvoriace kanal
        a_vector = [wall_lenght / 2, 0.0, 0.0]
        b_vector = [-horizontal_change * i, int(boxY), 0.0]
        c_vector = [0.0, 0.0, int(boxZ)]
    elif wall_type == "steny" and i >= 3:
        a_vector = [wall_lenght / 2, -height_change, 0.0]
        b_vector = [-horizontal_change * (i - 3), -height_change, 0.0]
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

maxCycle = 25
# main integration loop
cell.output_vtk_pos_folded(file_name=directory+"/cell_0.vtk")
cell1.output_vtk_pos_folded(file_name=directory+"/cell1_0.vtk")
for i in range(1, maxCycle):
    # if (i % 2 == 0):
    #     lbf.print_vtk_velocity(directory + "/fluid.vtk")
    system.integrator.run(steps=500)
    cell.output_vtk_pos_folded(file_name=directory+"/cell_" + str(i) + ".vtk")
    cell1.output_vtk_pos_folded(file_name=directory+"/cell1_" + str(i) + ".vtk")
    print("Iteration: " + str(i))
    print("time: ", str(i * system.time_step * 500))
    print("......... ")
    # print(cell1.get_velocity()[0], cell.get_velocity()[0])
    # print("......... ")
print("Simulation completed.")
