#!/usr/bin/env python

'''
VASP utility
Compute g(r) for simple cubic cell from XDATCAR
'''

                                                                   
import sys
import numpy as np


#function to compute the number of atoms in a range between r and r+dr from Nth atom 
def countAtomsInDistRtoDR(nthatom,Nfatoms,fr,fdr,POSfarray):

    distTrueCount = int(0) #initialize atom count to 0 within tolerance distance

    for j in range(0, Nfatoms):
        if nthatom != j: #skip the current atom with which we are computing distances to all other atoms
            for kx in range (-1, 2): #periodic image loop (note cell length is 1.0)
                for ky in range (-1, 2): #periodic image loop
                    for kz in range (-1, 2): #periodic image loop
                        xdiff = POSfarray[nthatom,0]-POSfarray[j,0]+kx*1.0 #+kx*1.0 because the cell size is 1.0 in direct coordinates
                        ydiff = POSfarray[nthatom,1]-POSfarray[j,1]+ky*1.0
                        zdiff = POSfarray[nthatom,2]-POSfarray[j,2]+kz*1.0
                        dist = np.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
                        distance = np.absolute(dist)
                        if distance >= fr and distance <= fr+fdr:
                            distTrueCount += 1

    return distTrueCount



###### Open and read XDATCAR; compute volume and rho=N/V ######
with open('XDATCAR', "r") as infile, open('Output.txt', 'w') as outfile:
    infile.readline() # comment line
    infile.readline() # alatt assumed 1.0 in this code

    acellread1=infile.readline().split() # vec1
    acell1=float(acellread1[0])
    acellread2=infile.readline().split() # vec2
    acell2=float(acellread2[1])
    acellread3=infile.readline().split() # vec3
    acell3=float(acellread3[2])

    CellVolume = acell1*acell2*acell3 #simple cubic cell

    #print CellVolume
    infile.readline() # species or blank

    Natoms=int(infile.readline()) # Natoms
    #print Natoms
    rho=Natoms/CellVolume
    #print rho
    
    infile.readline() #read Direct configuration line


########## Store the XDATCAR coordinates #############
    Coordarray = np.zeros((Natoms,3)) # (atom, coord)
    for i in range(Natoms): #loops over 0 to Natoms-1
        coordString = infile.readline().split()
        for j in range(3): #x,y,z coords
            Coordarray[i,j] = float(coordString[j])
    #print Coordarray


############## Compute g(r) #################
nBins=int(20) #number of bins
binNumber=int(0) #initialize bin counter to 0
GofRarray = np.zeros((Natoms,nBins)) # (atom, coord)
#for r in np.arange(0.05, 0.5, 0.05):
for r in np.linspace(0.05, 1.0, num=nBins):
    #print r
    dr=0.01 #shell thickness
    shellVol = 4*np.pi*acell1*r*acell1*r*dr*acell1 #volume of shell
    for ithatom in range(Natoms):    
        NthDistCount = countAtomsInDistRtoDR(ithatom,Natoms,r,dr,Coordarray)
        #print NthDistCount
        GofRarray[ithatom,binNumber] = NthDistCount/(shellVol*rho) # dividing actual count by ideal count
    binNumber += 1

#print GofRarray

#average the GoR data for each atom:
GofRarrayAvg = np.zeros(binNumber) # (atom, coord)
maxbinNumber=binNumber
binNumber = int(0) #reset bin count to 0

for i in range(maxbinNumber): #loop over bins
    BinSum=float(0.0) #reset BinSum to zero for each bin
    for ithatom in range(Natoms): #sum over atoms for like-bins
        BinSum += GofRarray[ithatom,binNumber]
        #print GofRarray[ithatom,binNumber]
    AvgBinSum = BinSum/Natoms
    GofRarrayAvg[binNumber]=AvgBinSum
    binNumber += 1

print GofRarrayAvg
