"""
Particle3D, a class to describe particles in 3D

Maya Khela & Francesca Elverson
"""
import math
import numpy as np
import sys
import random


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
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.mass = mass
        self.label = label

    def __str__(self):
        """
        Define output format
        For particle p = ([x,y,z],[vx,vy,vz],mass,label) this will print as
        "<label> x y z"
        This is the VMD compatible format
        """
        return str(self.label)+' '+str(self.pos[0])+' '+str(self.pos[1])+' '+str(self.pos[2])

    def kinetic_energy(self):
        """
        Return kinetic energy as 1/2*mass*vel^2
        """
        return 0.5*self.mass*(self.vel[0]**2 + self.vel[1]**2 + self.vel[2]**2)
        
    def leap_vel_single(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m
        
        :param dt: timestep as float
        :param force: force on particle as array of floats
        """
        force = np.array(force)
        self.vel += dt * force / self.mass

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
        self.pos += dt*self.vel + 0.5 * dt**2 * force/self.mass
        self.pos = np.mod(self.pos,L)    # Gives position of image of particle in original box

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
                    data[i] = float(data[i])    #Label argument left as string. !!!!!MAY CAUSE ERROR?!!!!
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
        separation = -(np.mod((par2.pos - par1.pos + shift) , L) - shift)
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
        if r < rc:
            lj_force = 48 * (r**-(14) - 2*r**-(8)) * Particle3D.separation(par1, par2, L)
        else:
            lj_force = 0

        return lj_force


#Test function
def main():
    par1 = Particle3D([0.,0.,0.],[0.0,0.0,0.0],1.,'par1')
    par2 = Particle3D([1.1,1.1,1.1],[0.,0.,0.],1.,'par2')
    outfile_name = sys.argv[1]
    outfile = open(outfile_name, "w")
    for i in range(5000):
        outfile.write(str(2)+"\n")
        outfile.write("Point = "+str(i+1)+"\n")
        outfile.write(str(par1)+"\n")
        outfile.write(str(par2) + "\n")
        force = Particle3D.lj_force_single(par1,par2,3.5,10)
        par1.leap_vel_single(0.005,force)
        par1.leap_pos_single(0.005,force,10)
        par2.leap_vel_single(0.005,-force)
        par2.leap_pos_single(0.005,-force,10)
main()
    
