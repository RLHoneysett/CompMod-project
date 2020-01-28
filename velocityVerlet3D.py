"""
CMod Ex2: velocity Verlet time integration of
two particles interacting via the Morse Potential

Produces plots of the seperation of the particles
and the system's energy, both as functions of time. Saves
both to file.
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

def pe_morse(particle1,particle2,re,de,a):
        """
        Gives the Morse potential energy between particles interacting via the Morse potential
        
        :param re: equilibrium bond distance
        :param de: parameter controlling depth of the potential minimum
        :param a: parameter controlling curvature of the potential minimum (small a -> wide well)
        """
        r12 = np.linalg.norm(Particle3D.separation(particle1,particle2))
        exp = math.exp(-a * (r12 - re))
        pe = de* (((1 - exp)**2) - 1)
        return pe


def force_morse(particle1,particle2,re,de,a):
        """
        Gives the Morse force between particles interacting via the Morse potential
        
        :param re: equilibrium bond distance
        :param de: parameter controlling depth of the potential minimum
        :param a: parameter controlling curvature of the potential minimum (small a -> wide well)
        """
        r12 = np.linalg.norm(Particle3D.separation(particle1,particle2))
        exp = math.exp(-a * (r12 - re))
        unit_r = np.true_divide(Particle3D.separation(particle1,particle2),r12)
        force = 2 * a * de * (1-exp) * exp * unit_r

        return force


#Main code
def main():
    # Read name of input and output files from command line
    if len(sys.argv) != 4:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <input file> " + "<separation output file> " + "<energy output file>")
        quit()
    else:
        datafile_name = sys.argv[1]
        outfile_separation_name = sys.argv[2]
        outfile_energy_name = sys.argv[3]

    # Open files
    outfile_separation = open(outfile_separation_name, "w")
    outfile_energy = open(outfile_energy_name, "w")
    """
    data file should be formated
    x1 y1 z1 Vx1 Vy1 Vz1 mass1 label1
    x2 y2 z2 Vx2 Vy2 Vz2 mass2 label2
    dt time numstep re de a
    """
    datafile = open(datafile_name, "r")
    
    # Set up particles and simulation parameters
    particle1 = Particle3D.from_file(datafile)
    particle2 = Particle3D.from_file(datafile)

    line = datafile.readline()
    param_data = line.split()
    dt = float(param_data[0])
    time = float(param_data[1])
    numstep = int(param_data[2])
    re = float(param_data[3])
    de = float(param_data[4])
    a = float(param_data[5])
    
    # Write initial conditions
    separation = np.linalg.norm(Particle3D.separation(particle1,particle2))
    energy = particle1.kinetic_energy() + particle2.kinetic_energy() + pe_morse(particle1,particle2,re,de,a)
    """
    ::::::FIX FORMAT FOR VMD USE:::::::

    outfile_separation.write("{0:f} , {1:4.4f}\n".format(time,separation))
    outfile_energy.write("{0:f} , {1:4.4f}\n".format(time,energy))
    """
    
    #Initial force (use '-force' for force on particle2)
    force = force_morse(particle1,particle2,re,de,a)

    # Initialise data lists for plotting
    time_list = [time]
    sep_list = [separation]
    energy_list = [energy]


    # Time integration loop
    for i in range(numstep):
        print('ITERATION ' + str(i+1)+':')
        print(str(particle1) +'  '+ str(particle2))
        # Update positions
        particle1.leap_pos2nd(dt,force)
        particle2.leap_pos2nd(dt,-force)
        
        # Calculate separation
        separation = np.linalg.norm(Particle3D.separation(particle1,particle2))
        
        # Calculate new force
        force_new = force_morse(particle1,particle2,re,de,a)

        # Update velocities with an average of current and new forces
        avg_force = 0.5 * (force + force_new)
        particle1.leap_velocity(dt,avg_force)
        particle2.leap_velocity(dt,-avg_force)

        # Update force value
        force = force_new
        
        # Update energy
        energy = particle1.kinetic_energy() + particle2.kinetic_energy() + pe_morse(particle1,particle2,re,de,a)
        
        # Update time
        time += dt * 10.1805057107594
        
        # Output information :::::::CLEAN THIS UP:::::::::
        outfile_separation.write(str(2)+"\n"+"Point = "+str(i+1)+"\n" + "s1 "+str(particle1.position[0]) +" "+str(particle1.position[1])+" "+str(particle1.position[2]) +"\n" +"s2 " + str(particle2.position[0]) + str(particle2.position[1])+str(particle2.position[2])+"\n")
        outfile_energy.write("{0:f} , {1:4.4f}\n".format(time,energy))

        # Append information to data lists
        time_list.append(time)
        sep_list.append(separation)
        energy_list.append(energy)
        
        

    # Close output files
    outfile_separation.close()
    outfile_energy.close()


    # Find the error in the energy
    energy_error = abs((max(energy_list)-min(energy_list))/energy_list[0])
    print ('Error in energy = '+ str(energy_error * 100) + ' %')
    
    # Find separation peaks and calculate period, frequency and wavenumber
    peak_locations = []
    for i in range((len(sep_list)-1)):
        if sep_list[i] > sep_list[i+1] and sep_list[i] > sep_list[i-1]:
                peak_locations.append(time_list[i])
        else:
                None

    if len(peak_locations)>=2:    
        period = (peak_locations[len(peak_locations)-1] - peak_locations[0])/(len(peak_locations)-1)
        frequency = 1/(period * 1.0e-15)
        wavenumber = 1/((period*1.0e-15)*29979245800)
        print ('Period = '+str(period)+'fs | Frequency = '+str(frequency)+' Hz | Wavenumber = '+str(wavenumber)+'cm^-1')
    else:
        print ('INSUFFICIENT TIME PERIOD. Use a smaller timestep or a larger number of steps')

    # Plots
    pyplot.title('Velocity Verlet: separation vs time')
    pyplot.xlabel('Time /fs')
    pyplot.ylabel('Separation /Å')
    pyplot.grid(which='both',axis='both')
    pyplot.plot(time_list, sep_list, color = 'purple')
    pyplot.show()

    pyplot.title('Velocity Verlet: energy vs time')
    pyplot.xlabel('Time /fs')
    pyplot.ylabel('Energy /eV')
    pyplot.grid(which='both',axis='both')
    pyplot.plot(time_list, energy_list, color = 'g')
    pyplot.show()

main()
    
