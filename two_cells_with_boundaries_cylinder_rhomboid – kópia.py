import espressomd
from espressomd import lb

from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions

import object_in_fluid as oif
import os, glob, sys, shutil
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sim_no", metavar="sim_no", type=int, help="simulation identifier")
args = parser.parse_args()

directory = "output/sim_du"+str(args.sim_no)
os.makedirs(directory)

boxX = 40.0 #boxX = 240.0
boxY = 20.0 #boxY = 80.0
boxZ = 20.0 #boxZ = 20.0
system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = 0.1

# creating the template for RBCs
type = oif.OifCellType(nodes_file="rbc374nodes.dat", triangles_file="rbc374triangles.dat",
                        check_orientation=False, system=system, ks=0.02, kb=0.016, kal=0.02,
                        kag=0.9, kv=0.5, resize=[4.0, 4.0, 4.0], normal=True)

# creating the RBCs
# subor polohy.txt: 2,3,12;3,12,18;...
# listOfCells b1, b2, b3,...
# particle_type - pre kazdu bunku iny
cell = oif.OifCell(cell_type=type, particle_type=0, origin=[5.0, 8.0, 6.0])
cell1 = oif.OifCell(cell_type=type, particle_type=1, origin=[5.0, 12.0, 13.0], rotate=[0.0, np.pi/2, 0.0], particle_mass=0.5)

print("Cells created.")

# cell-wall interactions
# cislo pre stenu listOfCells.size() + 1
system.non_bonded_inter[0, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(a=0.0002, n=1.2, cutoff=0.1, offset=0.0)

# cell-cell interactions
# dvojity cyklus
system.non_bonded_inter[0, 1].membrane_collision.set_params(a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

print("Interactions created.")

# fluid
lb_params = {'agrid': 1, 'dens': 1, 'visc': 1.5, 'tau': system.time_step, 'ext_force_density': [0.001, 0.0, 0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1.5)

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

# obstacle - cylinder A
tmp_shape = shapes.Cylinder(center=[21.0, 2.0, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=directory+"/cylinderA.vtk")

# obstacle - cylinder B
tmp_shape = shapes.Cylinder(center=[30.0, 8.0, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=directory+"/cylinderB.vtk")

# obstacle - cylinder C
tmp_shape = shapes.Cylinder(center=[21.0, 18.0, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=directory+"/cylinderC.vtk")

# obstacle - cylinder D
tmp_shape = shapes.Cylinder(center=[16.0, 8.0, 10.0], axis=[0.0, 0.0, 1.0], length=20.0, radius=3.0, direction=1)
boundaries.append(tmp_shape)
oif.output_vtk_cylinder(cyl_shape=tmp_shape, n=20, out_file=directory+"/cylinderD.vtk")

for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=10, penetrable=False)

print("Boundaries created.")

maxCycle = 300
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
