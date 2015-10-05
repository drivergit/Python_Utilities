#!/usr/bin/env python

'''
VASP utility by Kevin Driver 10-3-2015
Compute g(r) for simple cubic cell from XDATCAR
-NsnapshotMax and nBins are values set by the user

Functionality:
-Loops over several snapshots and computes and average g(r)


TODO:
-further smooth the g(r) by averaging over bins

'''

                                                                   
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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



    ###### Determine total number of snapshots in XDATCAR
    wc = subprocess.Popen(["wc", "XDATCAR"], stdout=subprocess.PIPE)
    awk = subprocess.Popen(["awk", "{print $1}"], stdin=wc.stdout, stdout=subprocess.PIPE)
    #print awk                                                                  
    output = float(awk.communicate()[0])
    NsnapshotsTotal = int(np.ceil((output-8)/(Natoms+1))) #remove top 8 lines, divide by Natoms+1blankline
    #print NsnapshotsTotal

    #Let's skip halfway down XDATCAR to ensure we are using equilibrated snapshots
    HalfNsnapshotsLine=(Natoms+1)*int(np.ceil(3*NsnapshotsTotal/4)) #number of lines to get to 75% down XDATCAR
    #print HalfNsnapshotsLine
    infile.readline() #read first Direct configuration line or blank

    #Next line to be read is line number 9; the first coordinate line.
    for i in range(HalfNsnapshotsLine):
        infile.readline() #move ahead 1 line until we are halfway down

    #checkline1 = infile.readline()
    #print checkline1



    Nsnapshot=int(0) #snapshot counter
    NsnapshotMax=int(4) # I want to loop over NsnapshotMax snapshots picked from the last half of the XDATCAR file
    for Nsnapshot in range(NsnapshotMax): #Big loop over snapshots

        ########## Store the XDATCAR coordinates #############
        Coordarray = np.zeros((Natoms,3)) # (atom, coord)
        for i in range(Natoms): #loops over 0 to Natoms-1
            coordString = infile.readline().split()
            if i == 0:
                print "first coord of snapshot",Nsnapshot, "is", coordString
            for j in range(3): #x,y,z coords
                Coordarray[i,j] = float(coordString[j])
        #print Coordarray


        #go ahead and advance the line number by a tenth of the second half
        infile.readline() #read Direct configuration line or blank
        fracofthelines = (Natoms+1)*int(np.ceil(1*NsnapshotsTotal/4/NsnapshotMax))
        for i in range(fracofthelines):
            infile.readline() #move ahead 1 line until a tenth of the way down
    
        #checkline1 = infile.readline()
        #print checkline1

        ############## Compute g(r) #################
        nBins=int(200) #number of bins
        binNumber=int(0) #initialize bin counter to 0
        GofRarray = np.zeros((Natoms,nBins)) # (atom, coord)
        distArray = np.zeros(nBins) # (atom, coord)
        #for r in np.arange(0.05, 0.5, 0.05):
        for r in np.linspace(0.05, 0.5, num=nBins): #r goes to L/2
            distArray[binNumber]=r*acell1 #r in angstroms for plotting g(r) at the end
            #print r
            dr=0.005 #shell thickness
            shellVol = 4*np.pi*acell1*r*acell1*r*dr*acell1 #volume of shell
            for ithatom in range(Natoms):    
                NthDistCount = countAtomsInDistRtoDR(ithatom,Natoms,r,dr,Coordarray)
                #print NthDistCount
                GofRarray[ithatom,binNumber] = NthDistCount/(shellVol*rho) # dividing actual count by ideal count
            binNumber += 1

        #print GofRarray

        #Average the GoR data for each atom:
        if Nsnapshot == 0: #intitalize GofRarrayAvg to zero for the first snapshot only
            GofRarrayAvg = np.zeros((NsnapshotMax,binNumber))
        maxbinNumber=binNumber
        binNumber = int(0) #reset bin count to 0

        for i in range(maxbinNumber): #loop over bins
            BinSum=float(0.0) #reset BinSum to zero for each bin
            for ithatom in range(Natoms): #sum over atoms for like-bins
                BinSum += GofRarray[ithatom,binNumber]
                #print GofRarray[ithatom,binNumber]
            AvgBinSum = BinSum/Natoms
            GofRarrayAvg[Nsnapshot,binNumber]=AvgBinSum #average for ith snapshot
            binNumber += 1

        #print GofRarrayAvg
        Nsnapshot += 1
    ###End of snapshot loop

    #Average g(r) over snapshots
    GofR = np.zeros(binNumber)
    binNumber = int(0) #reset bin count to 0
    for i in range(maxbinNumber): #loop over bins
        BinSum=float(0.0) #reset BinSum to zero for each bin
        for ithsnapshot in range(NsnapshotMax): #sum over atoms for like-bins
            BinSum += GofRarrayAvg[ithsnapshot,binNumber]
        AvgBinSum = BinSum/NsnapshotMax
        GofR[binNumber]=AvgBinSum
        binNumber += 1

print "Average GofR and the corresponding distance array:"
print " "
print GofR
print " "
print distArray




#######
#plot
#######
fig_size = [13 ,10]

mpl.rcParams['backend'] = 'ps'
mpl.rcParams['font.size'] = 30

mpl.rcParams['axes.labelsize'] = 30
mpl.rcParams['axes.linewidth'] = 3.0

mpl.rcParams['xtick.labelsize'] = 30
mpl.rcParams['xtick.major.size'] = 15
mpl.rcParams['xtick.major.width'] = 2.2
mpl.rcParams['xtick.minor.size'] = 8
mpl.rcParams['xtick.minor.width'] = 2.2
mpl.rcParams['xtick.major.pad'] = 15

mpl.rcParams['ytick.labelsize'] = 30
mpl.rcParams['ytick.major.size'] = 15
mpl.rcParams['ytick.major.width'] = 2.2
mpl.rcParams['ytick.major.pad'] = 15
mpl.rcParams['ytick.minor.size'] = 8
mpl.rcParams['ytick.minor.width'] = 2.2

mpl.rcParams['figure.figsize'] = fig_size

mpl.rcParams['legend.fontsize'] = 15
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.handlelength'] = 2.5
mpl.rcParams['legend.fontsize'] = 30

fig = plt.figure() #defines an overall 'big' figure that contains the subfigures; allows common axis label


fig1 = plt.subplot(111)
plt.plot(distArray, GofR,'r-o',linewidth=2.0,label='g(r)')

plt.xlabel('r ($\AA$)')
plt.xlim(0,4.0)
plt.xticks(np.arange(0,4.0,1.0))

plt.ylabel('g$_\mathrm{N-N}$(r)')
plt.ylim(0,2.5)
plt.yticks(np.arange(0,2.5,0.5))
minorLocatorx   = plt.MultipleLocator(0.1)
fig1.xaxis.set_minor_locator(minorLocatorx)
minorLocatory   = plt.MultipleLocator(0.1)
fig1.yaxis.set_minor_locator(minorLocatory)

plt.legend(loc='lower right',prop={'size':30},handlelength=3.5,handletextpad=0.4)


#plt.show()

plt.savefig('GofR.pdf', bbox_inches='tight')
#plt.savefig('GofR.eps', bbox_inches='tight')
