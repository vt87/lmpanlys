#-----------------------------------------------------------------------------------------------
# lmpanlys toolkit. Inspired from Pizza.py written by Steve Plimpton
# Vikram Thapar, vt87@cornell.edu, Cornell University
# integratemartini tool
# File Name: paircoeffs_martini.py
# The code will read martini.inp.
# martini.inp file contains the martini label for each type. 
# Example with type index followed by martini label
#  1     TC5  
#  2     TP1 
# This class also requires the values of cutoff and inner cut off (lj/gromacs)
# paircoeffs_martini in return will generate one files.
#     paircoeffs.txt : This file will contain the Martini field (LJ) pair coefficients
# Pair coefficients will be generated with the help of martini database file
# NOTE THAT CUTOFF UNITS ARE SPECIFIED IN "nm".
# MARTINI LJ PARAMETERS ARE ALSO IN "nm" and "kJ/mol".
# LAMMPS REAL UNITS ARE "angstrom" AND "kcal/mol". UNIT CONVERSION WILL BE DONE IN THIS CODE. 
#-----------------------------------------------------------------------------------------------

import numpy as np
import subprocess
import sys
import os

# Class definition
class paircoeffs_martini:
    def __init__(self,fname = "martini.inp",cutoff = 1.2,inner_cutoff = 0.9,fftype = "lj/gromacs", ofile = "paircoeffs.txt", database = "martini_v3.0.0.itp"):
        
        #Class variables
        self.cutoff = cutoff #martini lj cut off
        self.inner_cutoff = inner_cutoff #martini lj inner cut off (For lj/gromacs)
        self.fftype = fftype #Force field type
        self.inpdata = []  #data reading the martini.inp file
        self.inpdictswap = dict() #dictionary containing martini label as the key and type as value
        self.inpdict = dict() #swapping inpdict containing type as key and martini label as value
        self.paircoeffs = [] #pair coefficients for writing in lammps
        self.fname = fname # martini input filename
        self.pairall = "" #string containing the style and coefficients (will be written to file)
        self.ofile = ofile #ofile saving pair coefficients

        #Extract the pair coefficients from Martini database
        self.database = database

        #error check
        self.errorcheck()

        #CONSTANTS for unit conversion
        self.nmTOang = 10
        self.kjTOkcal = 1.0/4.184

        #Read martini input file
        self.read_martini_inp()

        #Allot pair coefficients
        self.allotpaircoeffs()
    
        #Do the unit conversion
        self.convertunits()

        #Prepare the string containing the pair style and its coefficients
        self.setpairstyle()

        #Write pair coefficients
        self.writepaircoeffs()


    #get parameters based on martini label, value1, value2
    def getparameters(self,value1,value2):
        #when martini labels are different
        if(value1 != value2):
            STRING = '/%s /&&/%s /'%(value1,value2) 
            COMMAND = "awk '%s' %s"%(STRING,self.database)
            proc = subprocess.Popen(COMMAND,shell = True, stdout = subprocess.PIPE)
            out = proc.communicate()[0].decode("utf-8")
            if(len(out) == 0):
                print("EXITING : No match for %s %s pair coeffs"%(value1,value2))
                sys.exit(-1)
            out = out.split("\n")
            out = out[:-1]
            matches = []
            for i in range(len(out)):
                outtmp = out[i].split()
                if( (outtmp[0] == value1 and outtmp[1] == value2) or (outtmp[1] == value1 and outtmp[0] == value2)):
                    matches.append(outtmp)
            if(len(matches) == 0):
                print("EXITING : No match for %s %s pair coeffs"%(value1,value2))
                sys.exit(-1)
            elif(len(matches) > 1):
                print("EXITING : Multiple matches for %s %s pair coeffs"%(value1,value2))
                sys.exit(-1)
            else:
                outtmp = matches[0]
                print("MATCH FOUND %s"%outtmp)
                outstr = "%lf %lf"%(float(outtmp[4]),float(outtmp[3]))
        #when martini labels are same
        else:
            STRING = '/%s/'%(value1) 
            COMMAND = "awk '%s' %s"%(STRING,self.database)
            proc = subprocess.Popen(COMMAND,shell = True, stdout = subprocess.PIPE)
            out = proc.communicate()[0].decode("utf-8")
            out = out.split("\n")
            STRING = "%s "%(value1)
            matches = []
            for i in range(len(out)-1):
                if(out[i].count(STRING) == 2):
                    outtmp = out[i].split()
                    if( (outtmp[0] == value1 and outtmp[1] == value1)):
                        matches.append(outtmp)
            if(len(matches) == 0):
                print("EXITING : No match for %s %s pair coeffs"%(value1,value2))
                sys.exit(-1)
            elif(len(matches) > 1):
                print("EXITING : Multiple matches for %s %s pair coeffs"%(value1,value2))
                sys.exit(-1)
            else:
                outtmp = matches[0]
                print("MATCH FOUND %s"%outtmp)
                outstr = "%lf %lf"%(float(outtmp[4]),float(outtmp[3]))
        return outstr

    #Read the input martini file starting with fname
    def read_martini_inp(self):
        f = open(self.fname,'r')
        for line in f:
            li = line.strip()
            if not li.startswith("#") and not li.startswith("$") and not li.startswith("@"):
                self.inpdata.append(line.split())
        f.close()
        inparray = np.asarray(self.inpdata)
        self.inpdict = dict(zip(inparray[:,0],inparray[:,1]))
        self.inpdictswap = dict(zip(inparray[:,1],inparray[:,0]))
        for key1,value1 in self.inpdict.items():
            print(key1,value1)

    #Allot the pair coefficients using database file of martini
    def allotpaircoeffs(self):
        for key1, value1 in self.inpdict.items():
            for key2,value2 in self.inpdict.items():
                #print(key1,value1,key2,value2)
                #cross interactions
                if(int(key1) < int(key2)):
                    outstr = self.getparameters(value1,value2)
                    outstr = "%s %s %s"%(key1,key2,outstr)
                    self.paircoeffs.append(outstr)
                #self interactions
                if(int(key1) == int(key2)):  
                    outstr = self.getparameters(value1,value2)
                    outstr = "%s %s %s"%(key1,key2,outstr)
                    self.paircoeffs.append(outstr)
                        
        print("PARAMETERS IN kJ/mol and nm")
        for i in range(len(self.paircoeffs)):
            print(self.paircoeffs[i])
    
    #Convert the units to lammps unit
    def convertunits(self):
        self.cutoff = self.cutoff*self.nmTOang
        self.inner_cutoff = self.inner_cutoff*self.nmTOang
        tmp = []
        print("PARAMETERS IN kcal/mol and Angstrom")
        for i in range(len(self.paircoeffs)):
            data = self.paircoeffs[i].split()
            data[2] = float(data[2])*self.kjTOkcal
            data[3] = float(data[3])*self.nmTOang
            self.paircoeffs[i] = "%s %s %.6f %.6f"%(data[0],data[1],data[2],data[3])
            print(self.paircoeffs[i])

    #Prepare the string to be written to pair coefficients file.
    def setpairstyle(self):

        #Set the pair style
        if(self.fftype == "lj/gromacs"):
            initstring = "pair_style lj/gromacs %s %s" %(self.inner_cutoff,self.cutoff)
        elif(self.fftype == "lj/cut"):
            initstring = "pair_style lj/cut %s" %(self.cutoff)
        else:
            print("Exiting : LAMMPS pair style is invalid")
            sys.exit(-1)
        self.pairall += "%s\n"%initstring

        #Set the coeffficients
        for i in range(len(self.paircoeffs)):
            pc = self.paircoeffs[i]
            if(self.fftype == "lj/gromacs"):
                string = "pair_coeff %s %s %s"%(pc,self.inner_cutoff,self.cutoff)
            else:
                string = "pair_coeff %s %s"%(pc,self.cutoff)
            self.pairall += "%s\n"%string

        #Set the shifting. Only for lj/cut
        if(self.fftype == "lj/cut"):
            string = "pair_modify shift yes"
            self.pairall += "%s\n"%string            
            
    #Write pair coefficients
    def writepaircoeffs(self):
        f = open(self.ofile,'w')
        f.write("# MARTINI force field \n")
        f.write("# The bead id and its martini label are as follows\n")
        for key,value in self.inpdict.items():
            f.write("# %s %s\n"%(key,value))
        f.write(self.pairall)
        f.close()
        
    #Error check
    def errorcheck(self):
        if not os.path.exists(self.fname):
            raise Exception("ERROR : Martini Input file named %s does not exist"%self.fname)
        if not os.path.exists(self.database):
            raise Exception("ERROR : Martini database file in GROMACS format, named %s does not exist"%self.database)
        if self.fftype == "lj/gromacs" or self.fftype == "lj/cut":
            pass
        else:
            raise Exception("ERROR : Invalid string, %s for type of forcefeld. Only lj/gromacs or lj/cut are allowed"%(self.fftype))
