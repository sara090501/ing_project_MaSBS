import espressomd
from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions
from espressomd.shapes import Rhomboid

import object_in_fluid as oif
import os
import numpy as np
import argparse
import math
import random
import csv

parser = argparse.ArgumentParser()
parser.add_argument("cells_positions_rotations_number", metavar="cells_positions_rotations_number",
                    type=int, help="cells positions rotations number")
parser.add_argument("cell_type", metavar="cell_type", type=str, help="cell type")
parser.add_argument("small_cell_size", metavar="small_cell_size", type=float, help="small cell size")
parser.add_argument("big_cell_size", metavar="big_cell_size", type=float, help="big cell size")
parser.add_argument("fluid_speed", metavar="fluid_speed", type=float, help="fluid speed")
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

print("\nCells positions and rotations number: " + str(args.cells_positions_rotations_number) + "\n" +
      "Cell type: " + str(args.cell_type) + "\n" +
      "Small cell size: " + str(args.small_cell_size) + "\n" +
      "Big cell size: " + str(args.big_cell_size) + "\n" +
      "Fluid speed: " + str(args.fluid_speed) + "\n" +
      "Simulation identifier: " + str(args.sim_no) + "\n")

directory_position_rotations = "input/cells_positions_rotations_" + str(args.cells_positions_rotations_number)
fluid_speed = args.fluid_speed
cell_type = args.cell_type
small_cell_size = args.small_cell_size
big_cell_size = args.big_cell_size

directory = "output/sim" + str(args.sim_no)
os.makedirs(directory)


LBgrid = 2


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

margin = 10  # mm od stien


# Kontrola zapisu origin-u bunky do csv
def lineForOrigin(x, y):
    return 40 * x + 40 * math.sqrt(3) * y - 6400 * math.sqrt(3) - 16800


# creating the template for RBCs
if cell_type == "bi":
    small_cell_type = oif.OifCellType(
                        nodes_file="input/Biclustersolid284nodes.dat",
                        triangles_file="input/Biclustersolid284triangles.dat",
                        check_orientation=False, system=system, ks=1.0, kb=0.0, kal=0.0,
                        kag=0.0, kv=0.0, resize=[small_cell_size, small_cell_size, small_cell_size], normal=False)
    big_cell_type = oif.OifCellType(
                        nodes_file="input/Biclustersolid284nodes.dat",
                        triangles_file="input/Biclustersolid284triangles.dat",
                        check_orientation=False, system=system, ks=1.0, kb=0.0, kal=0.0,
                        kag=0.0, kv=0.0, resize=[big_cell_size, big_cell_size, big_cell_size], normal=False)
else:
    small_cell_type = oif.OifCellType(nodes_file="input/solidball126nodes.dat",
                                      triangles_file="input/solidball126triangles.dat",
                                      check_orientation=False, system=system, ks=1.0, kb=0.0, kal=0.0,
                                      kag=0.0, kv=0.0, resize=[1.25, 1.25, 1.25], normal=False)
    big_cell_type = oif.OifCellType(nodes_file="input/solidball126nodes.dat",
                                    triangles_file="input/solidball126triangles.dat",
                                    check_orientation=False, system=system, ks=1.0, kb=0.0, kal=0.0,
                                    kag=0.0, kv=0.0, resize=[2.5, 2.5, 2.5], normal=False)
# ks = koeficient natiahnutia
# kb = koeficient ohybnosti
# kal = koeficient local area
# kag = koeficient global area
# kv = koeficient objemu
cells = []

# Počet malých a veľkých buniek
num_small_cells = 0 # počet závisí od počtu pozícií v cvs sobore nižšie
num_big_cells = 0 # počet závisí od počtu pozícií v cvs sobore nižšie

# Nacitanie nahodnych rotacii
x_angles = []
y_angles = []
z_angles = []
# Otvorenie CSV súboru na čítanie
with open(f"{directory_position_rotations}/cellsRotations.csv", 'r', newline='') as file:
    reader = csv.reader(file, delimiter=';')
    for row in reader:
        # Predpokladáme, že každá riadok obsahuje tri hodnoty
        if row:  # kontrola, že riadok nie je prázdny
            # Premeníme reťazce na čísla (float) a priradíme do jednotlivých listov
            x, y, z = map(float, row)
            x_angles.append(x)
            y_angles.append(y)
            z_angles.append(z)

with open(f"{directory_position_rotations}/smallCellPositions.csv", newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=';')  # Použitie bodkočiarky ako oddeľovača

    for i, row in enumerate(reader):  # Prechádza cez riadky CSV
        try:
            # Konverzia hodnôt na float
            x_random, y_random, z_random = map(float, row[:3])

            # Vytvorenie bunky
            cell = oif.OifCell(
                cell_type=small_cell_type,
                particle_type=i,
                origin=[x_random, y_random, z_random],
                inner_particles=False,
                rotate=[x_angles[i], y_angles[i], z_angles[i]]
            )

            # Pridanie bunky do zoznamu
            cells.append(cell)

        except ValueError:
            print(f"Chyba pri spracovaní riadku {i + 1}: {row}")  # Debug výstup

num_small_cells = len(cells)
print("\nNumber of small cells: " + str(num_small_cells))

with open(f"{directory_position_rotations}/bigCellPositions.csv", newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=';')  # Definovanie oddeľovača

    for i, row in enumerate(reader):  # Prechádza cez riadky CSV
        try:
            # Odstránenie medzier a konverzia na float
            x_random, y_random, z_random = map(lambda v: float(v.strip()), row[:3])

            # Vytvorenie bunky
            cell = oif.OifCell(
                cell_type=big_cell_type,
                particle_type=i,
                origin=[x_random, y_random, z_random],
                inner_particles=False,
                rotate=[x_angles[num_small_cells + i], y_angles[num_small_cells + i], z_angles[num_small_cells + i]]
            )

            # Pridanie bunky do zoznamu
            cells.append(cell)

        except ValueError:
            print(f"Chyba pri spracovaní riadku {i + 1}: {row}")  # Debug výstup

num_big_cells = len(cells) - num_small_cells
print("Number of big cells: " + str(num_big_cells))

print("Number of cells: " + str(len(cells)))
print("\nCells created.")

# --------------------------------------------------------------------------------------------------------------------------------------------------------
wall = num_small_cells + num_big_cells + 10 # je to lepsie takto, aby nevznikla chyba

# cell-wall interactions

# Nastavenie interakcií cell-wall pre malé bunky
for i in range(num_small_cells):
    system.non_bonded_inter[i, wall].soft_sphere.set_params(a=0.02, n=1.5, cutoff=1, offset=0.0)

# Nastavenie interakcií cell-wall pre veľké bunky
for i in range(num_big_cells):
    system.non_bonded_inter[num_small_cells + i, wall].soft_sphere.set_params(a=0.02, n=1.5, cutoff=1, offset=0.0)


print("\nInteractions created.")

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

# right velocity wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2 + boxX, height_in_turn + 2 * side_wall_width, 0.0],
                            a=[side_wall_width, 0.0, 0.0], b=[0.0, -height_in_turn - 2 * side_wall_width, 0.0],
                            c=[0.0, 0.0, canal_width_Z + 2 * side_wall_width], direction=1)
system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape, velocity=[
    fluid_speed, 0.0, 0.0]))  # toto velocity = [] je smer ktorym chem zacat
oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/velocityWallRight.vtk")

# front wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, canal_width_Z + side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0],
                            b=[0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
                            c=[0.0, 0.0, 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/frontWall.vtk")

system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape))
system.constraints.add(shape=tmp_shape, particle_type=wall, penetrable=False)


# back wall
tmp_shape = shapes.Rhomboid(corner=[-side_wall_width / 2, -2 * side_wall_width, -side_wall_width],
                            a=[boxX + side_wall_width / 2, 0.0, 0.0],
                            b=[0.0, 6 * side_wall_width + (height_in_turn + (math.sqrt(3) / 2) * wall_lenght), 0.0],
                            c=[0.0, 0.0, 2 * side_wall_width], direction=1)

oif.output_vtk_rhomboid(rhom_shape=tmp_shape,
                        out_file=directory + "/backWall.vtk")

system.lbboundaries.add(lbboundaries.LBBoundary(shape=tmp_shape))
system.constraints.add(shape=tmp_shape, particle_type=wall, penetrable=False)
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
        system.constraints.add(shape=segment, particle_type=wall, penetrable=False)

print("\nBoundaries created.")

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Hlavná integračná slučka
maxCycle = 50000

# Množina na sledovanie, či bunka už bola zaznamenaná za priamkou
cells_recorded_once = set()

# Export počiatočných pozícií buniek
for i, cell in enumerate(cells):
    cell.output_vtk_pos_folded(file_name=f"{directory}/cell_{i}_0.vtk")

for cycle in range(1, maxCycle):
    system.integrator.run(steps=20)

    # Zapis kvapaliny na 100. kroku
    if (cycle == 100):
        lbf.print_vtk_velocity(directory + "/fluid.vtk")

    # Množina na sledovanie zapísaných buniek v aktuálnom cykle
    recorded_cells = set()

    # Export pozícií buniek po každom cykle
    # for i, cell in enumerate(cells):
    #     cell.output_vtk_pos_folded(file_name=f"{directory}/cell_{i}_{cycle}.vtk")

    # Kontrola buniek, ktoré prešli za hranicu priamky
    for j in range(0, len(cells)):
        if j in recorded_cells or j in cells_recorded_once:
            continue  # Preskočíme bunky, ktoré už boli zapísané v tomto cykle alebo už boli zaznamenané

        origin = cells[j].get_origin()
        x = origin[0]
        y = origin[1]

        # Ak bunka prešla za hranicu priamky a je v druhej polovici geometrie
        if lineForOrigin(x, y) >= 0 and x > boxX / 2:
            #cells[j].output_vtk_pos_folded(file_name=f"{directory}/cell{j}_{cycle}.vtk")
            with open(f"{directory}/cell{j}positions.csv", 'a', newline='') as file:
                writer = csv.writer(file, delimiter=';')
                writer.writerow(origin)

            # Pridať bunku do množiny zapísaných a do celkovej sledovanej množiny
            recorded_cells.add(j)
            cells_recorded_once.add(j)

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Reseed buniek, ak prekročili hranice (voliteľné)
    for ii, cell in enumerate(cells):
        origin = cell.get_origin()
        if origin[0] > (1.5 * wall_lenght):
            new_x = origin[0] - wall_lenght
            cell.set_origin([new_x, origin[1], origin[2]])
            print(f"Reseeding cell: {ii}")

            # Po reseede odstrániť bunku zo sledovanej množiny
            if ii in cells_recorded_once:
                cells_recorded_once.remove(ii)
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------

    print(f"Iteration: {cycle}")
    print(f"time: {cycle * system.time_step * 500}")
    print("......... ")

print("Simulation completed.")


