import espressomd
from espressomd import lb

from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions

import object_in_fluid as oif
import os, glob, sys, shutil
import numpy as np
import argparse
import csv

#nova uloha week 8
# reseedovanie buniek - set_origin
# nahodne generovanie prekazok - stacia 4 - do strednej casti , 2/3 - kotrola prekrytia
# nahodne generovanie buniek - 8mikrometrov - 8 buniek
# dve polohy 10 a 230 - zmena y suradneci zo 240 na 10, x a z ostavaju
# vacsie prekazky, 1,5x vacsie
parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/model_week7" + str(args.sim_no)
os.makedirs(directory)

boxX = 240.0 #boxX = 240.0
boxY = 80.0 #boxY = 80.0
boxZ = 20.0 #boxZ = 20.0
system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = 0.1 # zmensi casovy krok - 0.05

# creating the template for RBCs
type = oif.OifCellType(nodes_file="input/rbc374nodes.dat", triangles_file="input/rbc374triangles.dat",
                        check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                        kag=0.9, kv=0.5, resize=[4.0, 4.0, 4.0], normal=True)

# CELLS ----------------------------------------------------------------------------------------------------------------
listOfCells = []
file = open("input/cells_position_week7.txt", "r")

particle_type = 0
for line in file:
    a_str, b_str, c_str = line.split(" ")

    a = float(a_str)
    b = float(b_str)
    c = float(c_str)

    cell = oif.OifCell(cell_type=type, particle_type=particle_type, origin=[a, b, c])
    listOfCells.append(cell)
    print(a, b, c)
    particle_type = particle_type + 1

print("Cells created.")

# INTERACTIONS ----------------------------------------------------------------------------------------------------------
# cell-wall interactions
for i in range(0, len(listOfCells)):
    system.non_bonded_inter[i, len(listOfCells) + 1].soft_sphere.set_params(a=0.02, n=1.5, cutoff=1.0, offset=0.0)
    # nove nastavenie
    # n=2, cutoff=1.5

# cell-cell interactions
for i in range(0, len(listOfCells)):
    for j in range(i + 1, len(listOfCells)):
        system.non_bonded_inter[i, j].membrane_collision.set_params(a=0.02, n=1.5, cutoff=1.0, offset=0.0)

print("Interactions created.")

# fluid
lb_params = {'agrid': 1, 'dens': 1, 'visc': 1.5, 'tau': system.time_step, 'ext_force_density': [0.001, 0.0, 0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)
    # gamma - ako velmi sa snazi bunka priblizit v rychlosti tekutiny

# WALLS ----------------------------------------------------------------------------------------------------------------
boundaries = []

# bottom of the channel
tmp_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0], b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0],
                          direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_rhomboid(rhom_shape=tmp_shape, out_file=directory+"/wallBottom.vtk")

# top of the channel
tmp_shape = shapes.Rhomboid(corner=[0.0, 0.0, boxZ-1], a=[boxX, 0.0, 0.0], b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0],
                            direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_rhomboid(rhom_shape=tmp_shape, out_file=directory+"/wallTop.vtk")

# front wall of the channel
tmp_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0], b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ],
                            direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_rhomboid(rhom_shape=tmp_shape, out_file=directory+"/wallFront.vtk")

# back wall of the channel
tmp_shape = shapes.Rhomboid(corner=[0.0, boxY-1.0, 0.0], a=[boxX, 0.0, 0.0], b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ],
                            direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_rhomboid(rhom_shape=tmp_shape, out_file=directory+"/wallBack.vtk")

# CYLINDERS ------------------------------------------------------------------------------------------------------------
# prvy cyklus pre 3 cylindre

# Nastavenie parametrov pre jednotlivé valce
y_positions = [15.0, 40.0, 65.0]  # y-pozície pre jednotlivé valce
cylinder_names_1 = ["A", "B", "C"]  # Názvy súborov

# Prechádzanie cez y-pozície a názvy súborov v cykle
for y, name in zip(y_positions, cylinder_names_1):
    # Vytvorenie valca s aktuálnymi parametrami
    tmp_shape = shapes.Cylinder(center=[160.0, y, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
    # Pridanie valca do zoznamu hraníc
    boundaries.append(tmp_shape)
    # Export valca do VTK formátu so správnym názvom
    oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=f"{directory}/cylinder{name}.vtk")

# druhy cyklus pre 2 cylindre
y_positions_2 = [27.5, 52.5]
cylinder_names_2 = ["D", "E"] # Názvy súborov

# Prechádzanie cez y-pozície a názvy súborov v cykle
for y, name in zip(y_positions_2, cylinder_names_2):
    # Vytvorenie valca s aktuálnymi parametrami
    tmp_shape = shapes.Cylinder(center=[185.0, y, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)

    # Pridanie valca do zoznamu hraníc
    boundaries.append(tmp_shape)

    # Export valca do VTK formátu so správnym názvom
    oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=f"{directory}/cylinder{name}.vtk")

# treti cyklus pre 3 cylindre
cylinder_names_3 = ["F", "G", "H"]  # Názvy súborov

# Prechádzanie cez y-pozície a názvy súborov v cykle
for y, name in zip(y_positions, cylinder_names_3):
    # Vytvorenie valca s aktuálnymi parametrami
    tmp_shape = shapes.Cylinder(center=[210.0, y, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
    # Pridanie valca do zoznamu hraníc
    boundaries.append(tmp_shape)
    # Export valca do VTK formátu so správnym názvom
    oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=f"{directory}/cylinder{name}.vtk")

for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=len(listOfCells) + 1, penetrable=False)

print("Boundaries created.")

# SIMULATION -----------------------------------------------------------------------------------------------------------
maxCycle = 2
for i in range(0, len(listOfCells)):
    listOfCells[i].output_vtk_pos_folded(file_name=directory+"/cell" + str(i) + "_0.vtk")

directoryForCellsPositions = "output/model_week7" + str(args.sim_no) + "/cellsPositions"
os.makedirs(directoryForCellsPositions)

for i in range(1, maxCycle):
    system.integrator.run(steps=100)
    print(str(lbf[100,40,10].velocity))
    for j in range(0, len(listOfCells)):
        listOfCells[j].output_vtk_pos_folded(file_name=directory+"/cell" + str(j) + "_" + str(i) + ".vtk")
        with open(directoryForCellsPositions + '/cell' + str(j) + 'positions.csv', 'a', newline='') as file:
            writer = csv.writer(file, delimiter=';')
            origin = listOfCells[j].get_origin()
            writer.writerow(origin)
    print("Iteration: " + str(i))
    print("time: ", str(i*system.time_step*500))
    print("----------------------------------------")
print("Simulation completed.")
