"""
particle_manybody, a program to simulate many particles in 3D

Maya Khela & Francesca Elverson
"""
import math
import numpy as np
import sys
import random
import MDUtilites as mdu
import Particle3D as par3d

#calculate energies
def energy(particles, rc, L):
    KE = 0
    for i in range(len(particles)):
        KE += 0.5 * particles[i].mass * (np.linalg.norm(particles[i].velocity))^2
    
    PE = 0
    for i in range(len(particles)):
        for j in range(len(particles)):
            if j>i:
                PE += par3d.lj_potential_single(par[i], par[j], rc, L)
            else:
                continue
    energy = KE + PE
    return KE, PE, energy
#main code
def main():
    
    """
    Read sys arguments to define parameters
    **Particle data format**
    x_pos y_pos z_pos x_vel y_vel z_vel mass label
    
    **Param data format**
    numpar numstep dt rc T ρ ε σ
    
    **Traj Out format**
    numpar \\ point number \\ label[0] x_pos y_pos z_pos \\...\\ label[numpar] x_pos y_pos z_pos

    **Obsv Out format**
    !!!TO BE DECIDED!!!
    
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
    with open(traj_output, "w") as trajout and open(obsv_output, "w") as obsvout:
        for i in range(numstep):
        
