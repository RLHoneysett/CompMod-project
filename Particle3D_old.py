e"""
CompMod Ex2: Particle3D, a class to describe particles in 3D
"""
import math
import numpy as np
from matplotlib import pyplot as plt

class Particle3D(object):
    """
    Class to describe 3D particles
    
    Properties:
    position[float,float,float] - position in 3D space
    velocity[float,float,float] - velocity in 3D space
    mass(float) - particle mass
    label(string) - label for the particle
    
    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first and second order position updates
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
        self.mass = mass
        self.label = label

    def __str__(self):
        """
        Define output format
        For particle p = ([x,y,z],[vx,vy,vz],mass,label) this will print as
        "<label> x y z"
        """
        return str(self.label)+' '+str(self.position[0])+' '+str(self.position[1])+' '+str(self.position[2])

    def kinetic_energy(self):
        """
        Return kinetic energy as 1/2*mass*vel^2
        """
        return 0.5*self.mass*(self.velocity[0]**2 + self.velocity[1]**2 + self.velocity[2]**2)
        
    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)
        
        :param dt: timestep as float
        :param force: force on particle as array of floats
        """
        self.velocity += np.multiply(dt,np.true_divide(force,self.mass))
    
    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)
        
        :param dt: timestep as float
        """
        self.position += np.multiply(dt,self.velocity)

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as array of floats
        """
        self.position += np.multiply(dt,self.velocity) + np.multiply(0.5*(dt**2),np.true_divide(force,self.mass))

    @staticmethod
    def from_file(data_file):
        """
        Read content from file
        Call Particle3D __init__ method
        """
        line = data_file.readline()
        data = line.split()
        for i in range(7):
            data[i] = float(data[i])
        return Particle3D([data[0],data[1],data[2]],[data[3],data[4],data[5]],data[6],data[7])

    @staticmethod
    def separation(particle1,particle2):
        """
        Gives the vector separation between two particles
        """
        return np.subtract(particle2.position,particle1.position)



