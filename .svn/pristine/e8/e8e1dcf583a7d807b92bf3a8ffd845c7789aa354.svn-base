"""

  Routine for generating track.dat given angle, speed and initial state of eye
  Nicole Beisiegel, August 2013

"""
import numpy as np
import math
import tkMessageBox

#-------------------------------------------------------------------------------
# Determine parameters for track data

# Initial vortex position[(m,m)], speed [m/s], angle[degree], final time[h]
#                      +x degree
#                    /
#                   /
# angle is given by --- 0 degree
#                   \
#                    \
#                      -x degree

eye_init = np.array([0.0, 0.0])
v        =  25.0
angle    =  0.0
final_time = 30.0

# Holland's model
pc = 950.0
pn = 1005.0
A  = 23.0
B  = 1.5
#-------------------------------------------------------------------------------
"""
   FUNCTION for preprocessing input data
"""

def final_position(v,eye_init,final_time,angle):
    # Compute directional vector
    vec = np.array([math.cos(angle), math.sin(angle)])
    # Compute final time in sec
    ftime=final_time*float(3600.0)
    # Final position of vortex
    pos_final=np.zeros([2])
    for i in range(2):
        pos_final[i] = eye_init[i]+ftime*v*vec[i]
    return pos_final, ftime

#------------------------------------------------------------------------------
"""
   FUNCTION for printing the track.dat
"""

def print_input(pos_final, ftime):
    input = open("track_storm.dat", "w")
    input.write(str(pn)+ " "+ str(pc) + " " +str(B) + " " +str(A) + " "+ str(0.0) + " " +str(ftime) + " " +str(ftime) + " " +str(2) +"\n")
    input.write(str(0.0)+" " +str(0.0)+"\n")
    input.write(str(pos_final[0]) +" "+ str(pos_final[1])+ "\n")
    input.close()
    tkMessageBox.showinfo("Print Message ",
                          "Track data file has been successfully written.")

#------------------------------------------------------------------------------
# Creation of file
#------------------------------------------------------------------------------

(pos_final,time)=final_position(v,eye_init,final_time,angle)
print_input(pos_final, time)
