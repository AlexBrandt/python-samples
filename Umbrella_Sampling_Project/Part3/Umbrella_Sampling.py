import numpy
import random
import math
import matplotlib.pyplot as plt
import time

CC = 20 # Coupling constant
ExField = 0
UmbrellaWindows = 8
EdgeBuffer = 1

#This method is essentially a collection of all of the steps of the Monte Carlo process
#Its arguments are:
#n_steps: The number of steps you'd like the monte carlo algorithm to iterate
#edge_length: The edge dimension of your lattice.  So, for a 20x20 lattice, this would equal 20
#ext_field: The strength of the external field
#temp: The temperature
#start_flag: An optional flag which, if set to 1 will read a starting configuration from an external file.
#config_file: Another optional flag which specifies the external file which initializes the lattice.

def combineHistograms(new, total):
    for k in new.keys():
        if k in total.keys(): total[k] += new[k]
        else: total[k] = new[k]
    return total

def stripZeros(histogram):
    for k in histogram.keys():
        if histogram[k] is 0:
            del histogram[k]
    return histogram

def normalizeAndLogHistogram(histogram, temp):
    total = 0
    for k in histogram.keys(): total += histogram[k]
    for k in histogram.keys(): histogram[k] /= float(total)
    for k in histogram.keys(): histogram[k] = -1 * math.log(histogram[k])
    return histogram

def raiseHistogram(height_offset,histogram):
    for k in histogram.keys():
        histogram[k] += height_offset
    return histogram


def stitchHistogram(TotalHistogram,Start,End,Interval):
    junctions = []
    offset = 0.0
    previous_value = 0.0
    switch = False
    # Bad code I use: junctions = [50, -50, 0]
    for i in range(0, UmbrellaWindows-2):
        junctions.append(Interval * i)
        junctions.append(-Interval * i )
    #print junctions
    # Go backwards...
    for m in range(Start,End+1)[::-1]:
        if m in TotalHistogram.keys():
            if switch:
                offset = (previous_value - TotalHistogram[m])
                switch = False
            TotalHistogram[m] = TotalHistogram[m] + offset
            previous_value = TotalHistogram[m]
            if m in junctions:
                switch = True
    return TotalHistogram


def PerformRun(n_steps,edge_length,ext_field,temp,start_flag=0,config_file=None):
    LatticeEdgeDimension = edge_length
    ExField = ext_field
    T = temp

    End = edge_length**2 #TODO: verify
    Start = -edge_length**2
    Interval = int(((End-Start)/UmbrellaWindows))
    TotalHistogram = {}
    tempAddHistogram = {}
    ProbabilityHistogram = {}
    for x in range(UmbrellaWindows):
        myStart = Start + (Interval * x)
        myEnd = min(Start + Interval * (x + 1), End)
        print myStart
        print myEnd
        for i in range(20):
            SpinLattice, energy, netmag = initialize(start_flag,LatticeEdgeDimension,ExField,myStart,myEnd)
            SpinLattice, energy, netmag,histogram = monte_carlo(n_steps,SpinLattice,T,energy,netmag,myStart,myEnd)
            ProbabilityHistogram = combineHistograms(histogram, ProbabilityHistogram)
            # Keeps track off all histograms samples from ONE umbrella region
            tempAddHistogram = combineHistograms(histogram, tempAddHistogram)

        # Strip zeroes to prevent math errors.  The rest is just calculating before
        tempAddHistogram = normalizeAndLogHistogram(stripZeros(tempAddHistogram), temp)
        TotalHistogram = combineHistograms(tempAddHistogram, TotalHistogram)

        tempAddHistogram = {}

    before = TotalHistogram.copy()
    # Combine each new histogram
    TotalHistogram = stitchHistogram(TotalHistogram,Start,End,Interval)

    return SpinLattice, energy, netmag, ProbabilityHistogram, before, TotalHistogram

#Initialize, if given a start_flag of 0, initializes the lattice to be entirely magnetized.
#If start_flag is 1, it reads from config_file for the lattice configuration

def initialize(start_flag,edge_length,field,myStart,myEnd):
    SpinLattice = [[0 for x in xrange(edge_length)] for y in xrange(edge_length)]
    netmag = 0
    energy = 0.
    ExField = field
    magTarget = 1
    #while magTarget % 2 != 0:
    magTarget = random.choice(range(myStart + EdgeBuffer,myEnd - EdgeBuffer))


    print "Mag target set..."
    print magTarget
    print "Generating magnetized configuration...\n"
    for x in xrange(edge_length):
        for y in xrange(edge_length):
            SpinLattice[x][y]=1
            netmag += 1

        energy = -2.0*(edge_length**2)

    flipsNeeded = (netmag - magTarget) / 2
    possibleSites = []
    for x in range(edge_length):
        for y in range(edge_length): possibleSites.append((x, y))
    random.shuffle(possibleSites)
    print len(possibleSites)
    print flipsNeeded
    for i in range(flipsNeeded):
        posTuple = possibleSites[i]
        SpinLattice[posTuple[0]][posTuple[1]] *= -1
        netmag += -2
        energy += 4.0

    print "Initial magnetization per spin = " + str(netmag/(edge_length**2)) + "\n"
    print "Initial energy per spin = " + str(energy/(edge_length**2)) + "\n"

    return SpinLattice,energy,netmag

#This method contains the bulk of the calculation.  Code that collects statistics, however, is missing.
def monte_carlo(n_steps,spin,T,energy,netmag,start,end):

    histogram = {}
    avmag=0.0
    aven=0.0
    edge_length = len(spin[0])

    # with open("trajectory.dat",'w') as f:

    #Each 'step' of the algorithm attempts edge_length**2 trial moves.
    print "Generating trajectory...\n"
    for step in xrange(n_steps):
        for j in xrange(edge_length**2):
            x = int(edge_length*numpy.random.rand())
            y = int(edge_length*numpy.random.rand())

            spin,energy,netmag = trial_move(x,y,spin,energy,netmag,T,start,end)

            #if netmag == whatever m state you want:
            #    createSpinPlot(spin)

            if netmag not in histogram.keys():
                histogram[netmag] = 1
            else:
                histogram[netmag] += 1
            avmag = (avmag * (step) + netmag) / (step + 1) # step = (number of moves) - 1
            aven  = (aven  * (step) + energy) / (step + 1)

        for val in range(start,end):
            if val not in histogram.keys():
                histogram[val] = 0 #This fills in missing values

        #/* Output averages */
        print "Average magnetization per spin = " + str(avmag/(n_steps*edge_length**2)) + "\n"
        print "Average energy per spin = " + str(aven/(n_steps*edge_length**2)) + "\n"


    return spin, energy, netmag, histogram

#This method actually attempts the trial move, and decides whether to accept or reject the trial.
def trial_move(x,y,spin,energy,netmag,T,start,end):
    neighbor_mag=0
    edge_length = len(spin[0])
    up,down,left,right = 0,0,0,0
    deltae = 0.0
    # if start == 200: end = 401
    magRange = range(start,end)

    if (x==0):
        left=spin[edge_length-1][y]
    else:
        left=spin[x-1][y]
    if (x==edge_length-1):
        right=spin[0][y]
    else:
        right=spin[x+1][y]
    if (y==edge_length-1):
        up=spin[x][0]
    else:
        up=spin[x][y+1]
    if (y==0):
        down=spin[x][edge_length-1]
    else:
        down=spin[x][y-1]

    #/* left, right, up, and down are the states */
    #/* of spins neighboring (x,y) */
    #/* compute change in energy if spin[x][y] were flipped, */
    my_old_spin = spin[x][y]
    neighbor_sum = up + down + left + right
    delta_E = 2 * ExField * my_old_spin + 2 * CC * my_old_spin * neighbor_sum
    #/*************************************************************/

    #/* accept according to Metropolis Monte Carlo rules */
    #/* update magnetization and energy if necessary */
    if ((netmag - 2 * my_old_spin) in magRange and min(1, math.exp(-1 * delta_E / T)) > random.random() ):
        # Flip the spin if min(1, e^(-beta E)) > random float (from [0, 1))
        energy += delta_E
        spin[x][y] = -1 * my_old_spin
        netmag -= (2 * my_old_spin)
    #/*************************************************************/
    return spin, energy, netmag


#Writes a configuration to file.  Spins are written as integers, tab delimited.
def write_config(filnam,spin):
    with open(filnam,'w') as f:
        for x in xrange(len(spin[0])):
            f.write(str("\t".join([str(s) for s in spin[x]])) + "\n")

#Reads from a configuration file.  Reads from the same format as the write_config method above--i.e. tab delimited integers, with an equal
#number of rows as columns.
def read_config(filnam):
    spin = []
    with open(filnam,'r') as f:
        for line in f:
            ThisRow = []
            for s in line.split("\t"):
                ThisRow.append(int(s))
            spin.append(ThisRow)

    return spin

def createSpinPlot(SpinLattice):
    edge_length = len(SpinLattice)
    SpinColor = [[0 for x in xrange(edge_length)] for y in xrange(edge_length)]
    for y in range(edge_length):
        for x in range(edge_length):
            if SpinLattice[x][y] == 1:
                SpinColor[x][y] = [1.0, 1.0, 1.0]
            else:
                SpinColor[x][y] = [0.0, 0.0, 1.0]
    SpinArray = numpy.array(SpinColor)
    plt.imshow(SpinArray, interpolation='nearest')
    plt.show()
    plt.clf()

#This is how you might call the program to perform 10000 steps in a 20x20 lattice at 0 applied field and at 1.0 K, where the lattice is initially fully magnetized

tempPoints = []
amagPoints = []
temp = 40
step = 1

SpinLattice, energy, netmag, ProbabilityHistogram, histogram, histogramAfter = PerformRun(10,20,0,temp)
HistogramPoints = []
MagPoints = []

ProbabilityHistogramPoints = []
PMagPoints = []

ProbabilityHistogram = stripZeros(ProbabilityHistogram)

for x in range(-400,401):
    if x in ProbabilityHistogram.keys():
        ProbabilityHistogramPoints.append(ProbabilityHistogram[x])
        PMagPoints.append(x)
        print "At %i, %f" % (x, ProbabilityHistogram[x])

plt.plot(PMagPoints,ProbabilityHistogramPoints,'ro')
plt.show()
plt.clf()

for x in range(-400,401):
    if x in histogram.keys():
        HistogramPoints.append(histogram[x])
        print "At %i, %f" % (x, histogram[x])
        MagPoints.append(x)

plt.plot(MagPoints,HistogramPoints,'ro')
plt.show()
plt.clf()

OtherMagPoints = []
OtherHistogramPoints = []

current_min = 90
current_min_m = 1000

for x in range(-400,401):
    if x in histogramAfter.keys():
        #current_min = min(current_min, histogramAfter[x])
        #current_min_m = x
        #Use this code to get the Beta * A(M) for the 2nd derivative.
        OtherHistogramPoints.append(histogramAfter[x])
        OtherMagPoints.append(x)
        print "At %i, %f" % (x, histogramAfter[x])
plt.plot(OtherMagPoints,OtherHistogramPoints,'ro')
plt.show()

#print current_min
#print current_min_m

#This is how you might call a program to perform the same, except on a user-specified input file:
#Note that when using a user-specified input file, the edge dimension argument is ignored (in this case, 20 is ignored, and the
#edge dimension is set to the dimensions of the input file lattice
#PerformRun(10000,20,0.0,1.0,1,"startconfig.in")
