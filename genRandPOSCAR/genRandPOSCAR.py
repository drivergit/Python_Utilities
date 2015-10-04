#!/usr/bin/env python
                                                                   
import sys
import numpy as np



def computeDistances(n,Nfatoms,distTolf,POSfarray):

    distCheck = int(1) #initialize to pass

    for j in range(0, Nfatoms):
        if n != j: #skip the current atom with which we are computing distances to all other atoms
            for kx in range (-1, 2): #periodic image loop (note cell length is 1.0)
                for ky in range (-1, 2): #periodic image loop
                    for kz in range (-1, 2): #periodic image loop
                        xdiff = POSfarray[n,0]-POSfarray[j,0]+kx*1.0 #+kx*1.0 because the cell size is 1.0 in direct coordinates
                        ydiff = POSfarray[n,1]-POSfarray[j,1]+ky*1.0
                        zdiff = POSfarray[n,2]-POSfarray[j,2]+kz*1.0
                        dist = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                        distance = np.absolute(dist)
                        if distance < distTolf:
                            #                print "dist too low"
                            distCheck = int(0) #fails
                            break #leave the for loop
                    else:
                        continue
                    break
                else:
                    continue
                break
            else:
                continue
            break

    return distCheck




#generate initial set of random atom positions
###############################################
Natoms = 128
distTol = 0.165 #x/alatt = distTol 
POSarray = np.zeros((Natoms,3))
for i in range(0, Natoms):
    for j in range(0, 3):
        POSarray[i,j]=np.random.uniform(0,1)

print "initial POSarray = ",POSarray


#enforce that the atoms are sufficiently seperated
##################################################
for i in range(0, Natoms):
    POSarraytemp = np.zeros((Natoms,3)) #temporary array to hold positons as distances are checked
    failcount = 0
    checkDist = int(0) # initialize checkDist
    checkDist = computeDistances(i,Natoms,distTol,POSarray) #get the pass/fail value of checkdist
#    print "first checkDist check: ",checkDist

    if checkDist == 1:
        print "checkDist passed on first gen for atom ",i
    
#    print "checking distance for atom ",i
    while checkDist == int(0):
        for j in range(0, 3):
            POSarraytemp[i,j]=np.random.uniform(0,1) #generate a new position
            POSarray[i,j]=POSarraytemp[i,j]  #update the position in the orginal array           

        checkDist = computeDistances(i,Natoms,distTol,POSarray) #check if the new position satisfies distance criteria
        failcount += 1

        if checkDist == 0:
#            print "dist check failed again",i
            if failcount == 50000:
                sys.exit("Failed to find a position in 50000 iterations")

        if checkDist == 1:
            print "checkDist now satisfied for atom",i, "with", failcount, "tries"



print POSarray




