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
parser.add_argument("csv_no", metavar="csv_no", type=int, help="csv identifier")
parser.add_argument("small_cell_size", metavar="small_cell_size", type=float, help="number of small cells")
parser.add_argument("big_cell_size", metavar="big_cell_size", type=float, help="number of small cells")
parser.add_argument("num_small_cells", metavar="num_small_cells", type=int, help="number of small cells")
parser.add_argument("num_big_cells", metavar="num_big_cells", type=int, help="number of big cells")
parser.add_argument("cell_type", metavar="cell_type", type=str, help="cell type")
parser.add_argument("generate_positions", metavar="generate_positions", type=str, help="generate positions")
parser.add_argument("generate_rotations", metavar="generate_rotations", type=str, help="generate rotations")
args = parser.parse_args()

print("Number of csvs: " + str(args.csv_no) + "\n" +
      "Small cell size: " + str(args.small_cell_size) + "\n" +
      "Big cell size: " + str(args.big_cell_size) + "\n" +
      "Number of small cell size: " + str(args.num_small_cells) + "\n" +
      "Number of big cell size: " + str(args.num_big_cells) + "\n" +
      "Cell type: " + str(args.cell_type) + "\n" +
      "Generate positions: " + str(args.generate_positions) + "\n" +
      "Generate rotations: " + str(args.generate_rotations) + "\n")

directory = "input/cells_positions_rotations_" + str(args.csv_no)
os.makedirs(directory)

# Veľkosti malých a veľkých buniek
small_cell_size = args.num_small_cells
big_cell_size = args.num_small_cells

# Počet malých a veľkých buniek
num_small_cells = args.num_small_cells
num_big_cells = args.num_big_cells
num_of_cells = num_small_cells + num_big_cells

# Typ buniek
cell_type = args.cell_type

# generovanie
generate_positions = args.generate_positions
generate_rotations = args.generate_rotations

# Medzera od stien
margin = 10


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
                                                                                y) >= 0


# Kontrola zapisu origin-u bunky do csv
def lineForOrigin(x, y):
    return 40 * x + 40 * math.sqrt(3) * y - 6400 * math.sqrt(3) - 16800


def random_positions(angle):
    x_random = random.uniform(0, angle)
    x_int = int(x_random)
    y_random = random.uniform(0, angle)
    x_int = int(y_random)
    z_random = random.uniform(0, angle)
    x_int = int(z_random)
    return [x_random, y_random, z_random]


# -----------------------------------------------------------------------------------------------------------------------
# Funkcia na vytvorenie bunky s podobným seeding procesom ako v originálnom kóde
# Generovanie pozicii malých buniek
if generate_positions == "Y" or generate_positions == "y":
    if cell_type == "bi":
        small_cell_size = small_cell_size * 2
        big_cell_size = big_cell_size * 2

    for i in range(num_small_cells):
        while True:
            x_random = random.uniform(0, 140)
            y_random = random.uniform((160 - (40 * math.sqrt(3))) + 15, (140 * math.sqrt(3)) - 15)

            # Rozdelenie na dve skupiny podľa indexu
            if i < num_small_cells // 2:
                z_random = random.uniform(9, 17)  # Prvá polovica v Z rozsahu [9, 17]
            else:
                z_random = random.uniform(23, 31)  # Druhá polovica v Z rozsahu [23, 31]

            if inside_rectangle(x_random + small_cell_size + margin, y_random + small_cell_size + margin) and \
                    inside_rectangle(x_random - small_cell_size - margin, y_random - small_cell_size - margin):
                with open(f"{directory}/smallCellPositions.csv", 'a', newline='') as file:
                    writer = csv.writer(file, delimiter=';')
                    origin = [x_random, y_random, z_random]
                    writer.writerow(origin)
                break  # Vygenerovaná validná bunka, ukončenie cyklu

    # Generovanie veľkých buniek
    for i in range(num_big_cells):
        while True:
            x_random = random.uniform(0, 140)
            y_random = random.uniform((160 - (40 * math.sqrt(3))) + 15, (140 * math.sqrt(3)) - 15)

            # Rozdelenie na dve skupiny podľa indexu
            if i < num_big_cells // 2:
                z_random = random.uniform(9, 17)  # Prvá polovica v Z rozsahu [9, 17]
            else:
                z_random = random.uniform(23, 31)  # Druhá polovica v Z rozsahu [23, 31]

            if inside_rectangle(x_random + big_cell_size + margin, y_random + big_cell_size + margin) and \
                    inside_rectangle(x_random - big_cell_size - margin, y_random - big_cell_size - margin):
                with open(f"{directory}/bigCellPositions.csv", 'a', newline='') as file:
                    writer = csv.writer(file, delimiter=';')
                    origin = [x_random, y_random, z_random]
                    writer.writerow(origin)
                break  # Vygenerovaná validná bunka, ukončenie cyklu

    print("Cells positions created.")

if generate_rotations == "Y" or generate_rotations == "y":
    for i in range(num_of_cells):
        # Generovanie náhodného uhla v radiánoch medzi 0 a 2*pi
        random_angle_x = random.uniform(0, 2 * math.pi)
        random_angle_y = random.uniform(0, 2 * math.pi)
        random_angle_z = random.uniform(0, 2 * math.pi)

        with open(f"{directory}/cellsRotations.csv", 'a', newline='') as file:
            writer = csv.writer(file, delimiter=';')
            origin = [random_angle_x, random_angle_y, random_angle_z]
            writer.writerow(origin)

    print("Cells rotations created.")
