"""
particle_manybody, a program to simulate many particles in 3D

Maya Khela & Francesca Elverson
"""
import math
import numpy as np
import sys
import random
import MDUtilities as mdu
from Particle3D import Particle3D as par3d
import time

start_time = time.time()

#main code
def main():
    
    """
    Read sys arguments to define parameters
    **Particle data format**
    x_pos y_pos z_pos x_vel y_vel z_vel mass label
    
    Label should just be the element type
    
    **Param data format**
    numpar numstep dt rc T ρ ε σ
    
    **Traj Out format**
    numpar \\ point number \\ label[0] x_pos y_pos z_pos \\...\\ label[numpar] x_pos y_pos z_pos

    **Obsv Out format**
    KineticEnergy  PotentialEnergy  TotalEnergy \\ MeanSquaredDisplacement \\ RadialDistributionFunction

    **Sys argument format**
    parData paramData trajOut obsvOut
    """
    # Open input and output files, set parameters
    par_input = sys.argv[1] #don't open this file, par3d.from_file opens it in the method
    param_input = sys.argv[2]
    traj_output = sys.argv[3]
    obsv_output = sys.argv[4]
    
    with open(param_input, "r") as paramdata:
        line = paramdata.readline()
        param = line.split()
        numpar = int(param[0])
        numstep = int(param[1])
        dt = float(param[2])
        rc = float(param[3])
        T = float(param[4])
        rho = float(param[5])
        epsilon = float(param[6])
        sigma = float(param[7])
        
    # Initialise particles
    particles = []
    for i in range(numpar):
        particles.append(par3d.from_file(par_input))
        particles[i].label += str(i)
    
    
    # Set box length and initial conditions
    L = mdu.set_initial_positions(rho, particles)[0]
    mdu.set_initial_velocities(T, particles)

    # List of particles in their initial state as ref point
    particles_init = particles

    # Time integrator with outputs
    with open(traj_output, "w") as trajout , open(obsv_output, "w") as obsvout:
        len_par = len(particles)
        obsvout_counter = 0
        trajout_counter = 0
        for i in range(numstep):
        # Write trajectory data every 5 steps
            if trajout_counter != 5 and i != 0:
                trajout_counter += 1
            else:
                trajout.write("%s \n Point = %s \n"%(numpar,i+1))
                for j in particles:
                    trajout.write("%s \n"%(j))
                trajout_counter = 0
        
            # Write observable data every 20 steps and print step to terminal
            if obsvout_counter != 19 and i != 0:
                obsvout_counter += 1
            else:
                ke,pe,e_total = par3d.energy(particles, rc, L)
                msd = par3d.mean_squared_displacement(particles, particles_init, L)
                rdf = par3d.radial_distribution(particles, L)
                obsvout.write("%s  %s  %s \n %s \n %s"%(ke, pe, e_total, msd, rdf))
                print("SIMULATION STEP %s OF %d"%(i,numstep))
                obsvout_counter = 0

            # Velocity and position updates every step
            par3d.leap_pos(particles,dt,rc,L)
            par3d.leap_vel(particles,dt,rc,L)

main()

print("--- %s seconds ---" %(time.time() - start_time))
