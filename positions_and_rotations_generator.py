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
parser.add_argument("side_wall_width", metavar="side_wall_width", type=int, help="side wall width")
parser.add_argument("canal_width_z", metavar="canal_width_z", type=int, help="canal width z")
args = parser.parse_args()

print("Number of csvs: " + str(args.csv_no) + "\n" +
      "Small cell size: " + str(args.small_cell_size) + "\n" +
      "Big cell size: " + str(args.big_cell_size) + "\n" +
      "Number of small cell size: " + str(args.num_small_cells) + "\n" +
      "Number of big cell size: " + str(args.num_big_cells) + "\n" +
      "Cell type: " + str(args.cell_type) + "\n" +
      "Generate positions: " + str(args.generate_positions) + "\n" +
      "Generate rotations: " + str(args.generate_rotations) + "\n" +
      "Side wall width: " + str(args.side_wall_width) + "\n" +
      "Canal width Z: " + str(args.canal_width_z) + "\n")

directory = "input/cells_positions_rotations_" + str(args.csv_no)
os.makedirs(directory)

# Veľkosti malých a veľkých buniek
small_cell_size = args.small_cell_size
big_cell_size = args.big_cell_size

# Počet malých a veľkých buniek
num_small_cells = args.num_small_cells
num_big_cells = args.num_big_cells
num_of_cells = num_small_cells + num_big_cells

# Typ buniek
cell_type = args.cell_type

# generovanie
generate_positions = args.generate_positions
generate_rotations = args.generate_rotations

# udaje o geometrii kanalu simulacie
side_wall_width = args.side_wall_width
canal_width_z = args.canal_width_z
center_of_canal_z = ((side_wall_width * 2) + canal_width_z) // 2

# Medzera od stien
margin = 10

# Medzera od dna a stropu
margin_top_and_bottom = 2


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
    new_small_cell_size = small_cell_size
    new_big_cell_size = big_cell_size

    if cell_type == "bi":
        new_small_cell_size = small_cell_size * 2
        new_big_cell_size = big_cell_size * 2

    for i in range(num_small_cells):
        while True:
            x_random = random.uniform(0, 140)
            y_random = random.uniform((160 - (40 * math.sqrt(3))) + 15, (140 * math.sqrt(3)) - 15)

            while True:
                z_random = random.uniform(
                    side_wall_width + margin_top_and_bottom + small_cell_size,
                    side_wall_width + canal_width_z - margin_top_and_bottom - small_cell_size)
                if z_random != center_of_canal_z:
                    break

            if inside_rectangle(x_random + new_small_cell_size + margin, y_random + new_small_cell_size + margin) and \
                    inside_rectangle(x_random - new_small_cell_size - margin, y_random - new_small_cell_size - margin):
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

            while True:
                z_random = random.uniform(
                    side_wall_width + margin_top_and_bottom + small_cell_size,
                    side_wall_width + canal_width_z - margin_top_and_bottom - small_cell_size)
                if z_random != center_of_canal_z:
                    break

            if inside_rectangle(x_random + new_big_cell_size + margin, y_random + new_big_cell_size + margin) and \
                    inside_rectangle(x_random - new_big_cell_size - margin, y_random - new_big_cell_size - margin):
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
