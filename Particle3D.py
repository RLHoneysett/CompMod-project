"""
Particle3D, a class to describe particles in 3D

Maya Khela & Francesca Elverson
"""
import math
import numpy as np
import sys
import random
import MDUtilities as MDU
import matplotlib.pyplot as plt
import time

class Particle3D(object):
    """
    Class to describe 3D particles
    
    Properties:
    position[float,float,float] - position in 3D space
    velocity[float,float,float] - velocity in 3D space
    mass(float) - particle mass
    label(string) - label for the particle
    
    Methods:
    """
    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle3D instance
        
        :param pos: position as array of floats
        :param vel: velocity as array of floats
        :param mass: mass as float
        """
        self.position = np.array(pos)
        self.velocity = np.array(vel)
        self.mass = float(mass)
        self.label = label

    def __str__(self):
        """
        Define output format
        For particle p = ([x,y,z],[vx,vy,vz],mass,label) this will print as
        "<label> x y z"
        This is the VMD compatible format
        """
        return '%s %s %s %s'%(self.label,self.position[0],self.position[1],self.position[2])

    def kinetic_energy(self):
        """
        Return kinetic energy as 1/2*mass*vel^2
        """
        return 0.5*self.mass*(self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)
        
    def leap_vel_single(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m
        
        :param dt: timestep as float
        :param force: force on particle as array of floats
        """
        force = np.array(force)
        self.velocity += dt * force / self.mass

    def leap_pos_single(self, dt, force, L):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)/m

        :param dt: timestep as float
        :param force: current force as array of floats
        :param L: side length of simulation box as float

        Includes Periodic Boundary Condition
        """
        force = np.array(force)
        self.position += (dt*self.velocity + 0.5 * dt**2 * force/self.mass)
        self.position = np.mod(self.position,L)    # Gives position of image of particle in original box

    @staticmethod
    def from_file(particle_input):
        """
        Read content from file
        Call Particle3D __init__ method

        Allows for creation of multiple different kinds of particle in a given ratio. Each time this method is 
        run a particle will be created for each line in the data file.
        """
        with open(particle_input,'r') as pardata:
            for line in pardata:
                data = line.split()
                for i in range(6):
                    data[i] = float(data[i])
                return Particle3D([data[0],data[1],data[2]] , [data[3],data[4],data[5]] , data[6] , data[7])

    @staticmethod
    def separation(par1, par2, L):
        """
        Gives the vector separation between two particles

        :param L: side length of simulation box as float

        Takes account of minimum image convention- Gives separation of closest image of particle2
        MIC algorithm shifts par2 L/2 up in all directions. If this moves it out of the simulation
        box, the algorithm gives the position of the image in the original box. It then shifts this
        image back to its original position. This is the location of the closest image.
        """
        shift = np.array([L/2, L/2, L/2])
        separation = -(np.mod((par2.position - par1.position + shift) , L) - shift)
        return separation

    @staticmethod
    def lj_force_single(par1, par2, rc, L):
        """
        Gives the Lennard_Jones force between two Particle3D objects as a vector. If particle2 is
        further than the cut-off radius rc then the force will be given as zero

        :param rc: cut-off distance
        :param L: box side length

        All values are in reduced units
        """
        r = np.linalg.norm(Particle3D.separation(par1, par2, L))
        if 0 < r < rc and par1 != par2 :
            lj_force = 48 * (1/(r**14) - 1/(2*r**8)) * Particle3D.separation(par1, par2, L)
        else:
            lj_force = 0

        return lj_force
        
    @staticmethod
    def lj_potential_single(par1, par2, rc, L):
        """
        Gives the Lennard-Jones potential energy between two interacting Particle3D objects.
        Gives the potential as zero if the particles have a separation greater than rc

        :param rc: cut-off distance

        All values are in reduced units
        """
        r = np.linalg.norm(Particle3D.separation(par1, par2, L))
        if 0 < r < rc:
            lj_pe = 4 * (r**-12 - r**-6)
        else:
            lj_pe = 0
        return lj_pe

    @staticmethod
    def lj_forces(particles,rc,L):
        """
        Takes a list of Particle3D objects
        Calculates the total force on each particle and returns a list of the forces

        :param rc: cut-off distance
        :L: box side length
        """
        start_time_forces=time.time()

        loop_range = len(particles)

        forces = np.zeros((loop_range,3))
        for i in range(loop_range):
            par1 = particles[i]
            for j in range(i):
                par2 = particles[j]
                if np.linalg.norm(Particle3D.separation(par1,par2,L)) < rc:
                    force_ij = Particle3D.lj_force_single(par1, par2,rc, L)
                    forces[i] += force_ij
                    forces[j] += -force_ij
                else:
                    continue
 
        #print("--- %s s for Forces to run ---"%(time.time()-start_time_forces))
        return forces
    
    @staticmethod
    def leap_pos(particles,dt,rc,L,forces):
        """
        Performs a second-order position update on every object in a list of
        Particle3D objects.

        :param dt: time-step
        :param rc: cut-off distance
        :param L: box side length
        """
        start_time_leappos = time.time()

        # Slower, may delete

        #forces = Particle3D.lj_forces(particles,rc,L)
        #loop_range = len(particles)
        #for i in range(loop_range):
            #force = forces[i]
            #particles[i].leap_pos_single(dt,force,L)

        def leap_pos_sing(particles,p,dt,forces,L):
            p.position += dt*p.velocity + 0.5*dt**2 * forces[particles.index(p)]/p.mass
            p.position = np.mod(p.position,L)
            return p
        particles = list(map(lambda p: leap_pos_sing(particles,p,dt,forces,L), particles))

        #print("--- %s s for leap_pos to run ---"%(time.time() - start_time_leappos))
        return particles
        
    @staticmethod
    def leap_vel(particles,dt,rc,L,forces):
        """
        Performs a first-order velocity update on every object in a list of
        Particle3D objects
        
        :param dt: time-step
        :param rc: cut-off distance
        :param L: box side length
        """
        start_time_leapvel = time.time()

        #Slower, may delete

        #forces = Particle3D.lj_forces(particles,rc,L)
        #loop_range = len(particles)
        #for i in range(loop_range):
            #force = forces[i]
            #particles[i].leap_vel_single(dt,force)

        def leap_vel_sing(particles,p,dt,forces):
            p.velocity += dt * forces[particles.index(p)] / p.mass
            return p
        particles = list(map(lambda p: leap_vel_sing(particles,p,dt,forces), particles))
            
        #print("---%s s for leap_vel to run ---"%(time.time()- start_time_leapvel))
        return particles

    @staticmethod
    def energy(particles, rc, L):
        """
        Calculates the kinetic, potential, and total energy of the system
        of Particle3D objects
        
        :param rc: cut-off distance
        :param L: box side length
        """
        KE = 0
        PE = 0
        for i in particles:
            KE += 0.5 * i.mass * (np.linalg.norm(i.velocity))**2
            
            for j in (particles):
                if particles.index(j) < particles.index(i):
                    PE += Particle3D.lj_potential_single(i, j, rc, L)

        energy = KE + PE
        return KE, PE, energy

    @staticmethod
    def mean_squared_displacement(particles, particles_init, L):
        """
        Compares the current and initial positions of every particle in the list to
        determine the average displacement of the system.

        :param L: box side length
        """
        displacement = 0
        loop_range = len(particles)
        for i in range(loop_range):
            displacement += (np.linalg.norm(Particle3D.separation(particles[i],particles_init[i],L)))**2
        msd = displacement / loop_range
        return msd

    @staticmethod
    def radial_distribution(particles, L):
        """
        Calculates the radial distribution function of the system by calculating the
        probability of finding a particle within a spherical slice dr, distance r away from a
        reference particle.
        
        :param L: box side length
        """
        distances = []
        for p in particles:
            for q in particles:
                if p != q:
                    distances.append(np.linalg.norm(Particle3D.separation(p,q,L)))
                else:
                    continue
        distances = np.around(np.sort(np.array(distances)), decimals=1 )
        
        #Loop through and add up how many particles are in each 0.1 bin
        rdf = []
        rdf_x_axis = np.around(np.arange(0,L,0.1),decimals=1)
        loop_range = len(rdf_x_axis)
        j = 0
        for i in range(loop_range):
            rdf.append(0)
            while j < len(distances) and rdf_x_axis[i] == distances[j]:
                rdf[i] += 1
                j+=1
            else:
                continue
        return np.array(rdf), rdf_x_axis
        
