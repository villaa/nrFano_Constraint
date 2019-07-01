#!/usr/bin/env python

import numpy as np
import pandas as pd
import re
import os
#===============to suppress h5py warning see:
#https://github.com/h5py/h5py/issues/961
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py
warnings.resetwarnings()
import argparse
from argparse import ArgumentParser, ArgumentTypeError
# Author: A. Villano
# The purpose here is to convert a Geant4 output text file
# to an h5 file with a workable output structure.
#

######################functions############################
def listFiles(path='./',regex=re.compile(r'(.*?)')):

        f = []
        for (dirpath, dirnames, filenames) in os.walk(path):
            passfiles = []
            for i,file in enumerate(filenames):
              if(regex.search(file) is not None):
                passfiles.extend([file])
            f.extend(passfiles)
            break

        print(f)
        return dirpath,f 

def readFile(filename):

        f = open(filename)

        #assume the first line is the tags
        tags = f.readline().split()
        data = [x.split() for x in f.readlines()] 
        
        #convert to numpy array
        data = np.asarray(data,dtype=np.float64)

        f.close()
        return data,tags 

def readFiles(flist,dirpath='./'):

        d = []
        data,tags = readFile(dirpath+flist[0]) #FIXME bug in simcode only puts correct header for file 0 -- eventually need to check if n colums is same for all files and return error if not
        for n,f in enumerate(flist):
          print('{} ({} out of {})'.format(dirpath+f,n,np.shape(flist)[0]))
          data,t = readFile(dirpath+f)
          d.extend(data)

        #convert to numpy array
        d = np.asarray(d,dtype=np.float64)

        return d,tags 

def saveh5(ofile='data.h5',path='./',regex=re.compile(r'(.*?)')):

        #get the data
        dirpath, f = listFiles(path,regex)
        d, tags = readFiles(f,dirpath)

        #open and write file
        of = h5py.File(ofile, 'w')

        #hits dataset
        dset_hits = of.create_dataset("geant4/hits", np.shape(d), dtype=np.dtype('float64').type, compression="gzip", compression_opts=9)
        dset_hits[...] = d      

        of.close()
        return

def nrFanoCondense(infile='datain.h5',outfile='data.h5',path='./'):

        f = h5py.File(path+'/'+infile,"r")

        data = f['geant4/hits']
        #using the examples above let's get byzantine and make a bad-ass data frame.

        #first do some cuts:
        #first some hit-level cuts
        cHVDet = np.zeros(np.shape(data)[0],dtype=bool)
        cZeroEdep = np.zeros(np.shape(data)[0],dtype=bool)
        cNeutron = np.zeros(np.shape(data)[0],dtype=bool)
        cGamma = np.zeros(np.shape(data)[0],dtype=bool)
        cNR = np.zeros(np.shape(data)[0],dtype=bool)
        
        cHVDet[data[:,1]==1] = True
        cZeroEdep[data[:,6]==0] = True
        cNeutron[data[:,4]==2112] = True
        cGamma[data[:,4]==22] = True
        cNR[data[:,4]>3000] = True
        
        #now make a dataframe with the restricted data
        nr_data = data[:,[0,4,5,6,21]]
        nr_data = nr_data[cHVDet&~cZeroEdep&cNR,:]
        nr_dataframe = pd.DataFrame(data=nr_data)
        
        #print(nr_dataframe)


        groupbyvec=[0]
        #print(np.max(nr_dataframe.groupby([0,1],axis=0).size()))
        max_vec = np.max(nr_dataframe.groupby(groupbyvec,axis=0).size())
        
        evec = np.zeros((0,max_vec))
        nhit = np.zeros((0,1))
        
        for i in nr_dataframe.groupby(groupbyvec,axis=0)[3].apply(list):
            #print(i)
            #print(np.shape(i))
            vector = np.zeros((1,max_vec))
            #print(np.shape(vector[0,0:np.shape(i)[0]]))
            vector[0,0:np.shape(i)[0]] = np.transpose(np.asarray(i))
            evec = np.append(evec,vector,0)
            nhit = np.append(nhit,np.shape(i)[0])
            
        print(np.shape(evec))
        print(np.shape(nhit))
        #print(max_vec)
        #print(np.sum(evec,1))
        print(np.sum(nhit[nhit==1]))
        

        #open and write file
        print(outfile)
        of = h5py.File(outfile, 'w')
        
        d = evec 
        #hits dataset
        dset_hits = of.create_dataset("nr_Fano/nr_energies", np.shape(d), dtype=np.dtype('float64').type, compression="gzip", compression_opts=9)
        dset_hits[...] = d
        
        d = nhit
        dset_hits = of.create_dataset("nr_Fano/nr_hits", np.shape(d), dtype=np.dtype('float64').type, compression="gzip", compression_opts=9)
        dset_hits[...] = d
        
        of.close()

        return
        
#the stuff below is so this functionality can be used as a script
########################################################################
if __name__ == "__main__":

        #make a parser for the input
        parser = argparse.ArgumentParser(description='Input processing specifications')
        parser.add_argument('-i','--infile', type=str, dest='infile', default=None, help='file to pre-process like ms_simulation_yield.ipynb')
        parser.add_argument('-d','--filedir', type=str, dest='filedir', default='./', help='directory to look for files')
        parser.add_argument('-x','--regex', type=str, dest='regstr', default=r'(.*?)', help='regex for picking up files')
        parser.add_argument('-o','--outfile', type=str, dest='outfile', default='data.h5', help='output file for data')
        #parser.set_defaults(filedir='./');

        args = parser.parse_args()

        try:
          #v = getPosition(args.iteration[0],args.n,pit) 

          print(args.filedir)
          print(args.regstr)
          print(args.outfile)
          if(args.infile==None):
            saveh5(args.outfile,args.filedir,re.compile(args.regstr))
          else:
            nrFanoCondense(infile=args.infile,outfile=args.outfile,path=args.filedir) 
        except KeyboardInterrupt:
          print('Shutdown requested .... exiting')
        except Exception:
          traceback.print_exc(file=sys.stderr)


