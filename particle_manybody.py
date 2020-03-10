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
from matplotlib import pyplot as plt
import copy

start_time = time.time()

#main code
def main():
    """
    Read sys arguments to define parameters
    **Particle data format**
    x_pos y_pos z_pos x_vel y_vel z_vel mass label
    
    Label should just be the element type
    
    **Param data format**
    numpar numstep dt rc T ρ
    
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
        
    # Initialise particles
    particles = []
    for i in range(numpar):
        particles.append(par3d.from_file(par_input))
        particles[i].label += str(i)

    # Set box length and initial conditions
    L = mdu.set_initial_positions(rho, particles)[0]
    mdu.set_initial_velocities(T, particles)

    # List of particles in their initial state as ref point
    particles_init = copy.deepcopy(particles)

    # Time integrator with outputs
    with open(traj_output, "w") as trajout , open(obsv_output, "w") as obsvout:

        len_par = len(particles)
        progress_five_percent = int(numstep/20)
        obsv_data_points = 0
        
        #test plots init
        msd_x = []
        msd_y = []
        energy_x = []
        energy_y = []
        ke_y = []
        pe_y = []
        
        rdf, rdf_x_axis = par3d.radial_distribution(particles,L)
        
        # Initial force
        forces = par3d.lj_forces(particles,rc,L)

        for i in range(numstep):
            # Write trajectory data every 3 steps
            
            if i%3 == 0:
                trajout.write("%s \n Point = %s \n"%(numpar,i+1))
                for j in particles:
                    trajout.write("%s \n"%(j))
        
            # Write observable data every 10 steps and print progress to terminal
            
            if i ==0 or i%i == 0:
                ke,pe,e_total = par3d.energy(particles, rc, L)

                msd=(par3d.mean_squared_displacement(particles, particles_init, L))
                msd_y.append(msd)
                msd_x.append(i)
                
                energy_x.append(i)
                ke_y.append(ke)
                pe_y.append(pe)
                energy_y.append(e_total)

                rdf += par3d.radial_distribution(particles,L)[0]
                obsvout.write("%s  %s  %s \n %s \n %s \n \n"%(ke, pe, e_total, msd, rdf))
                obsv_data_points += 1
                
            # Output progress to terminal
            if i%progress_five_percent == 0:
                print("%s%%  Completed"%(int(5*i/progress_five_percent)))

            # Velocity and position updates every step
            particles = par3d.leap_pos(particles,dt,rc,L,forces)
            forces_new = par3d.lj_forces(particles,rc,L)
            particles = par3d.leap_vel(particles,dt,rc,L,(forces+forces_new)/2)
            forces = forces_new

        # Normalise rdf
        rdf_norm = par3d.rdf_normalized(rho,rdf,rdf_x_axis)

        # Plot obsv data
        plt.title('RDF vs data points')
        plt.plot(rdf_x_axis, rdf/obsv_data_points)
        plt.show()
        plt.title('MSD vs data points')
        plt.plot(msd_x,msd_y)
        plt.show()
        plt.title('Energies vs data points')
        plt.plot(energy_x,energy_y, label='Total Energy')
        plt.plot(energy_x,ke_y, label='KE')
        plt.plot(energy_x,pe_y, label='PE')
        plt.legend()
        plt.show()
    
main()

print("--- %s seconds ---" %(time.time() - start_time))

