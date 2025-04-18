import espressomd
from espressomd import lb
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
import numpy as np
import os, sys, shutil
import object_in_fluid as oif

boxX = 20
boxY = 20
boxZ = 20
time_step = 0.1
system = espressomd.System(box_l=(boxX,boxY,boxZ))
system.time_step = time_step
system.cell_system.skin = 0.2

cell_type_rbc = oif.OifCellType(
    nodes_file="input/rbc482nodes.dat", triangles_file="input/rbc482triangles.dat",
    system = system, 
    ks=0.1, kb=0.1, kal=0.1, 
    kag=0.1, kv=0.1, 
    resize = [4.0, 4.0, 4.0],  normal=False)
cell_rbc = oif.OifCell(
    cell_type=cell_type_rbc, particle_type=0, 
    origin=[10.0,10.0,10.0])

lb_params = {'agrid': 1.0, 'dens': 1.0, 'visc':1.0, 'tau': system.time_step, 'ext_force_density':[0.01,0.0,0.0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

nSteps = 0
while nSteps < 50:
    print( "steps: " + str(nSteps))
    print("cell pos " + str(cell_rbc.get_origin()))
    system.integrator.run(steps = 10)
    nSteps += 1
print( "Simulation completed.")
