import os
import numpy as np
import fileinput
import math

#Define environment
#The path needs to be altered by user
PATH = os.getcwd()
file = PATH+"/Parameters.dat"
pyfile = PATH+"/pytest"

para = open(file, 'r')
data = para.readlines()
exe = open(pyfile, 'w')
char = "./DGM -f Parameters"

min=input("Choose minimal level of refinement\n")
max=input("Choose maximal level of refinement\n")
time=input("Choose timestep for min refinement level\n")
time_total=input("Choose final time\n")

#Save new files
for i in range(min,max):
  
    #Manipulate data
    data[11] = str(i)+"\n"
    data[14] = str(i)+"\n"

    #Adjust the timestep length according to rough estimate/ half every second refinement
    exp = math.ceil((i-min)/2)+1
    timestep = time/(2**exp)
    data[29] = str(timestep)+"\n"
    #Last timestep
    data[56] =str(int(time_total/timestep))+"\n"

    new_file = open(file.replace("Parameters","Parameters"+str(i)),'w')

    for j in range(0,102):
        new_file.write(data[j])
   
    new_file.close()
    if i == (int(max)-1):
        exe.write(char+str(i)+".dat"+"\n")
    else:
        exe.write(char+str(i)+".dat &&"+"\n")  
    

para.close()

exe.close()
#Write file for execution of tests
