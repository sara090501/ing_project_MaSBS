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
import math
import random

#nova uloha week 9
# vytvorenie geometrie
# zistenie dlzky stien kanala - 60 stupnov medzi kanalmi, sirka kanalu cez cos vysky, tvar pismeno M
# dopocitat sirku pri vrcholoch
# 4x LBGRID - romnoid co pusta tekutinu
# nejake bunky na skusku
# rychlost 62,5 raynolds number zigzag kennel
parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/model_week8_" + str(args.sim_no)
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
def is_collision(new_pos, existing_positions, min_distance):
    for pos in existing_positions:
        distance = math.sqrt((new_pos[0] - pos[0])**2 + (new_pos[1] - pos[1])**2)
        if distance < min_distance:
            return True
    return False

listOfCells = []
existing_positions = []
min_distance = 5.0

numberOfCells = 8
particle_type = 0

for _ in range(numberOfCells):
    z_fixed = 10.0
    while True:
        random_x = random.uniform(10, 70)
        random_y = random.uniform(10, 70)
        random_pos = (random_x, random_y, z_fixed)

        if not is_collision(random_pos, existing_positions, min_distance):
            break

    existing_positions.append(random_pos)
    cell = oif.OifCell(cell_type=type, particle_type=particle_type, origin=list(random_pos))
    listOfCells.append(cell)

    particle_type += 1

print("Cells created.")

# INTERACTIONS ----------------------------------------------------------------------------------------------------------
# cell-wall interactions
for i in range(0, len(listOfCells)):
    system.non_bonded_inter[i, len(listOfCells) + 1].soft_sphere.set_params(a=0.02, n=1.5, cutoff=0.5, offset=0.0)
    # nove nastavenie
    # n=2, cutoff=1.5

# cell-cell interactions
for i in range(0, len(listOfCells)):
    for j in range(i + 1, len(listOfCells)):
        system.non_bonded_inter[i, j].membrane_collision.set_params(a=0.02, n=1.5, cutoff=0.5, offset=0.0)

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
# nahodne generovane 4 cylindre
positions = []
cylinder_names_4 = ["I", "J", "K", "L"]  # Názvy súborov pre nové valce
min_distance = 10.0  # Minimálna vzdialenosť medzi valcami, aby sa neprekrývali


# Funkcia na kontrolu, či sa nová pozícia neprekrýva s existujúcimi
def is_position_valid(new_pos, positions, min_distance):
    for pos in positions:
        distance = math.sqrt((new_pos[0] - pos[0]) ** 2 + (new_pos[1] - pos[1]) ** 2)
        if distance < min_distance:
            return False
    return True


# Generovanie 4 náhodných pozícií pre (x, y) v daných intervaloch bez prekrývania
while len(positions) < 4:
    # Náhodne vyberieme pozície na osiach x a y v daných rozmedziach
    new_x = random.uniform(80.0, 160.0)
    new_y = random.uniform(15.0, 65.0)

    # Overíme, či sa nová pozícia neprekrýva s existujúcimi
    if is_position_valid((new_x, new_y), positions, min_distance):
        positions.append((new_x, new_y))  # Ak je pozícia platná, pridáme ju do zoznamu

# Vytvorenie valcov s generovanými (x, y) pozíciami
for (x, y), name in zip(positions, cylinder_names_4):
    # Vytvorenie valca s aktuálnymi parametrami
    tmp_shape = shapes.Cylinder(center=[x, y, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=8.0, direction=1)

    # Pridanie valca do zoznamu hraníc
    boundaries.append(tmp_shape)

    # Export valca do VTK formátu so správnym názvom
    oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=f"{directory}/cylinder{name}.vtk")

for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=len(listOfCells) + 1, penetrable=False)

print("Boundaries created.")

# SIMULATION -----------------------------------------------------------------------------------------------------------
maxCycle = 100
for i in range(0, len(listOfCells)):
    listOfCells[i].output_vtk_pos_folded(file_name=directory+"/cell" + str(i) + "_0.vtk")

directoryForCellsPositions = "output/model_week8_" + str(args.sim_no) + "/cellsPositions"
os.makedirs(directoryForCellsPositions)

for i in range(1, maxCycle):
    system.integrator.run(steps=100)
    # print(str(lbf[100,40,10].velocity)) - rychlost tekutiny na danych suradniciach
    for j in range(0, len(listOfCells)):
        if (listOfCells[j].get_origin()[0] > 230): #reseeding of cells
            listOfCells[j].set_origin((10.0, listOfCells[j].get_origin()[1], listOfCells[j].get_origin()[2]))
        print("Cell " + str(j) + " position x: " + str(listOfCells[j].get_origin()[0]))
        listOfCells[j].output_vtk_pos_folded(file_name=directory+"/cell" + str(j) + "_" + str(i) + ".vtk")
        with open(directoryForCellsPositions + '/cell' + str(j) + 'positions.csv', 'a', newline='') as file:
            writer = csv.writer(file, delimiter=';')
            origin = listOfCells[j].get_origin()
            writer.writerow(origin)
    print("Iteration: " + str(i))
    print("time: ", str(i*system.time_step*500))
    print("----------------------------------------")
print("Simulation completed.")
