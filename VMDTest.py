import sys
import math
import numpy
import random

def main():
    outfile_name = sys.argv[1]
    outfile = open(outfile_name,"w")
    numstep = int(input("Enter number of steps"))
    numpar = int(input("Enter number of particles"))
    for i in range(numstep):
        outfile.write(str(numpar)+"\n")
        outfile.write("Point = "+str((i+1))+"\n")
        for j in range(numpar):
            outfile.write("s"+str(j+1) +" "+str(random.random()*j+math.sin(i))+" "+str((2-j)+math.cos(i))+" "+str(j*random.random()/50 + 2*math.sin(i))+"\n")
main()
