"""
Particle3D, a class to describe particles in 3D

Maya Khela & Francesca Elverson
"""
import math
import numpy as np
import sys
import random
import MDUtilities as MDU


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
        self.mass = mass
        self.label = label

    def __str__(self):
        """
        Define output format
        For particle p = ([x,y,z],[vx,vy,vz],mass,label) this will print as
        "<label> x y z"
        This is the VMD compatible format
        """
        return str(self.label)+' '+str(self.position[0])+' '+str(self.position[1])+' '+str(self.position[2])

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
        self.position += dt*self.velocity + 0.5 * dt**2 * force/self.mass
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
        if r < rc and par1 != par2:
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
        if r < rc:
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
        forces = []
        for i in range(len(particles)):
            force_i = np.array([0,0,0])
            par1 = particles[i]
            for j in range(len(particles)):
                par2 = particles[j]
                r = np.linalg.norm(Particle3D.separation(par1,par2,L))
                if j != i and r < rc:
                    force_ij = Particle3D.lj_force_single(par1,par2,rc,L)
                    force_i = force_i + force_ij
                else:
                    continue
            forces.append(force_i.tolist())
        forces = np.array(forces)
        return forces

"""
I feel like we can add some functionality here so that once it gets far away enough it doesn't 
keep checking every particle. Will have to look at how the crystal forms to see at what point 
we can be sure we can stop counting.
"""
        


#Test function


def main():
#new test code for list of particles
    #Take user input
    numpar = int(input("Enter number of particles"))
    rc = float(input("Enter cut-off radius (suggested: 3.5)"))
    L = float(input("Enter box-size (suggested: 10)"))
    numstep = int(input("Enter number of steps"))
    dt = float(input("Enter timestep"))
    
    #Open input and output files and define parameters
    par_input = sys.argv[1]
    param_input = sys.argv[2]
    traj_output = sys.argv[3]
    
    with open(param_input, "r") as paramdata:
        line = paramdata.readline()
        param = line.split()
        epsilon = float(param[0])   #well depth
        sigma = float(param[1])     #van der Waals radius
        T = float(param[2])         #reduced temperature
        rho = float(param[3])       #reduced density
        

    trajout = open(traj_output, "w")

    #Initialise particles
    particles = []
    for i in range(numpar):
        particles.append(Particle3D.from_file(par_input))
    MDU.set_initial_positions(rho, particles)
    
    print (str(Particle3D.lj_forces(particles,rc,L)))
    
        
    



#does two particles
"""
    par1 = Particle3D([0.,0.,0.],[0.0,0.0,0.0],1.,'par1')
    par2 = Particle3D([1.1,1.1,1.1],[0.,0.,0.],1.,'par2')
    print(str(Particle3D.separation(par1,par2,10)))
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
"""
main()
    
