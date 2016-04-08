import numpy
import random
import math
import matplotlib.pyplot as plt

CC = 1 # Coupling constant
ExField = 0.0
H = 3
#This method is essentially a collection of all of the steps of the Monte Carlo process
#Its arguments are:
#n_steps: The number of steps you'd like the monte carlo algorithm to iterate
#edge_length: The edge dimension of your lattice.  So, for a 20x20 lattice, this would equal 20
#ext_field: The strength of the external field
#temp: The temperature
#start_flag: An optional flag which, if set to 1 will read a starting configuration from an external file.
#config_file: Another optional flag which specifies the external file which initializes the lattice.

def ExternalFieldBias(x,y):
    if y < 1:
        if x < 10:
            return -1 * H
        else:
            return H
    else:
        return 0

def PerformRun(n_steps,edge_length,ext_field,temp,start_flag=0,config_file=None):
    LatticeEdgeDimension = edge_length
    ExternalField = ext_field
    T = temp

    SpinLattice, energy, netmag = initialize(start_flag,LatticeEdgeDimension,ExternalField,config_file)
    SpinLattice, energy, netmag, stderr, cij = monte_carlo(n_steps,SpinLattice,T,energy,netmag)
    write_config("FinalConfig.out",SpinLattice)
    return SpinLattice, energy, netmag, stderr, cij

#Initialize, if given a start_flag of 0, initializes the lattice to be entirely magnetized.
#If start_flag is 1, it reads from config_file for the lattice configuration
def initialize(start_flag,edge_length,field,config_file=None):
    SpinLattice = [[0 for x in xrange(edge_length)] for y in xrange(edge_length)]
    netmag = 0
    energy = 0.
    ExField = field

    if (start_flag==0):
        print "Generating magnetized configuration...\n"
        for x in xrange(edge_length):
            for y in xrange(edge_length):
                SpinLattice[x][y]=1
                netmag += 1

        energy = -2.0*(edge_length**2)

    else:
        print "Reading configuration from file: " + str(config_file)

        SpinLattice = read_config(config_file)
        edge_length = len(SpinLattice[0])
        for x in xrange(edge_length):
            for y in xrange(edge_length):
                netmag += SpinLattice[x][y]
                up = 0
                right = 0

                if (x==edge_length-1):
                    right=SpinLattice[0][y]
                else:
                    right=SpinLattice[x+1][y]

                if (y==edge_length-1):
                    up=SpinLattice[x][0]
                else:
                    up=SpinLattice[x][y+1]

                energy -= SpinLattice[x][y]*(up+right)

    energy -= field*netmag

    print "Initial magnetization per spin = " + str(netmag/(edge_length**2)) + "\n"
    print "Initial energy per spin = " + str(energy/(edge_length**2)) + "\n"

    return SpinLattice,energy,netmag

#This method contains the bulk of the calculation.  Code that collects statistics, however, is missing.
def monte_carlo(n_steps,spin,T,energy,netmag):

    avmag=0.0
    aven=0.0
    edge_length = len(spin[0])

    HOWFAR = 6

    avmagList = []
    xlen = edge_length
    ylen = edge_length
    cij = [-1 for i in range(HOWFAR)] # Assume xlen = ylen...
    productRij = [[] for i in range(HOWFAR)]
    jSpinRij = [[] for i in range(HOWFAR)]
    mySpinList = []

    with open("trajectory.dat",'w') as f:

        #Each 'step' of the algorithm attempts edge_length**2 trial moves.
        print "Generating trajectory...\n"
        for step in xrange(n_steps):
            for j in xrange(edge_length**2):
                x = int(edge_length*numpy.random.rand())
                y = int(edge_length*numpy.random.rand())

                spin,energy,netmag = trial_move(x,y,spin,energy,netmag,T)

                xc, yc = 9, 0
                mySpin = float(spin[xc][yc])
                mySpinList.append(mySpin)
                for delta_x in range(HOWFAR): #TODO: Fix hard coding
                    '''Sj1 = float(spin[(xc + delta_x) % 20][yc])
                    # Sj2 = float(spin[(xc - delta_x) % 20][yc])
                    productRij[delta_x].append(mySpin * Sj1)
                    jSpinRij[delta_x].append(Sj1)
                    # productRij[delta_x].append(mySpin * Sj2)
                    # jSpinRij[delta_x].append(Sj2)'''


                    # Code for going "into" bulk
                    Sj1 = float(spin[(xc)][(yc + delta_x) % 20])
                    # Sj2 = float(spin[(xc - delta_x) % 20][yc])
                    productRij[delta_x].append(mySpin * Sj1)
                    jSpinRij[delta_x].append(Sj1)
                    # productRij[delta_x].append(mySpin * Sj2)
                    # jSpinRij[delta_x].append(Sj2)

            avmag = (avmag * (step) + netmag) / (step + 1) # step = (number of moves) - 1
            aven  = (aven  * (step) + energy) / (step + 1)

            #if (step%10==0):
            #    f.write("\t".join([str(s) for s in [step,netmag/(edge_length**2),energy/(edge_length**2)]]) + "\n")


        for r in range(HOWFAR):
            print r
            print numpy.mean(productRij[r])
            print float(sum(productRij[r]))/len(productRij[r])
            print numpy.mean(jSpinRij[r])
            print numpy.mean(mySpinList)
            cij[r] = (float(sum(productRij[r]))/len(productRij[r]) - (float(sum(jSpinRij[r]))/len(jSpinRij[r]) * float(sum(mySpinList))/len(mySpinList)))
            # raw_input()

        #/* Output averages */
        print "Average magnetization per spin = " + str(avmag/(n_steps*edge_length**2)) + "\n"
        print "Average energy per spin = " + str(aven/(n_steps*edge_length**2)) + "\n"
    stderr = (numpy.std(avmagList)/(n_steps**.5))
    return spin, energy, netmag, stderr, cij

#This method actually attempts the trial move, and decides whether to accept or reject the trial.
def trial_move(x,y,spin,energy,netmag,T):
    neighbor_mag=0
    edge_length = len(spin[0])
    up,down,left,right = 0,0,0,0
    deltae = 0.0

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
    delta_E = 2 * ExternalFieldBias(x,y) * my_old_spin + 2 * CC * my_old_spin * neighbor_sum

    #/*************************************************************/

    #/* accept according to Metropolis Monte Carlo rules */
    #/* update magnetization and energy if necessary */
    if (min(1, math.exp(-1 * delta_E / T)) > random.random() ):
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

rPoints, cPoints = [], []
tempPoints = []
amagPoints = []
initialTemp = 1.0
step = 1

for T in numpy.linspace(0,9,100):
    print T
    SpinLattice, energy, netmag, stderr, cij = PerformRun(100,20,0.0,T)
    tempPoints.append(T)
    amagPoints.append(netmag)
    # if deltaT == 13 or deltaT == 49 or deltaT == 99:
    #    createSpinPlot(SpinLattice)
    #    createCriticalPlots(SpinLattice)

    if T == 1.0 or T == 2.0:
        createSpinPlot(SpinLattice)
        for r in range(1,6):
            rPoints.append(r)
            cPoints.append(cij[r])
        plt.plot(rPoints,cPoints,'bo')
        plt.xlabel("r_ij")
        plt.ylabel("C(r_ij)")
        plt.show()
        plt.clf()
        rPoints = []
        cPoints = []

plt.plot(tempPoints,amagPoints,'ro')
plt.show()

#This is how you might call a program to perform the same, except on a user-specified input file:
#Note that when using a user-specified input file, the edge dimension argument is ignored (in this case, 20 is ignored, and the
#edge dimension is set to the dimensions of the input file lattice
#PerformRun(10000,20,0.0,1.0,1,"startconfig.in")
