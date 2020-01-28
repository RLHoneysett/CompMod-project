"""
Program to model the trajectories of randomly generated particles
Values given in LJ reduced units
E* = E/epsilon
m* = 1
r* = r/sigma
T* = epsilon/ k_B
t* = sigma sqrt(m/epislon)
"""

import sys
import math
import numpy as np
import random
from Particle3D import Particle3D

def pe_lj(particle1,particle2,sigma):
    """
    LJ params for Argon
    epsilon = 119.8 k_B
    sigma = 3.405e-10
    m = 0.03994 kg/mol
    """
    r = np.linalg.norm(Particle3D.separation(particle1,particle2))/sigma
    pe = 4*((r)*10^(-12) - (r)*10^(-6))
    return pe

def force_lj(particle1,particle2,sigma):
    r = np.linalg.norm(Particle3D.separation(particle1,particle2))/sigma
    force = 48*(r*10^(-14) - 0.5*r*10^(-8))* Particle3D.separation(particle1,particle2)/sigma
    return force

def user_input():
    natoms = int(input("Number of atoms: "))
    rho = float(input("Rho: "))
    sigma = float(input("Sigma in units of 10^-10: "))
    return natoms,rho,sigma

def main():
    natoms,rho,sigma = user_input()
    force = np.array([[0,0,0],[0,0,0]])
    position = np.array([[],[]])
    particles = []
    for i in range(natoms):
        particles.append(Particle3D([0,0,0] , [random.random(),random.random(),random.random()] , 1 , i))
    for j in range(natoms):
        for i in range(natoms):
            force[[i],[j]]= force_lj(particles[j],particles[i],sigma)
    #output files
    outfile_name = sys.argv[1]
    outfile = open(outfile_name, "w")
    print (str(natoms) +" "+ str(rho)+" " + str(sigma/(10e10)))

    for k in range (40):
        force = 0
        for i in range(natoms):
            for j in range(natoms-1):
                force += force[[i][j]]
            particles[i].leap_pos2nd(0.01,force)
            position[[i],[k]] = particles[i].position
main()
