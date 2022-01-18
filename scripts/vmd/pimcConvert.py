import numpy as np
import scipy as sp
import pandas as pd
import h5py
import argparse
import os
import re
import tqdm

def _pdbAtomRow( index,iParticle, x,y,z,group,iTime ):
    _atoms=["N","O"]
    template="ATOM  {:5d} {}    HIS A{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f} 00.00           {}"
    return template.format(index,_atoms[group],iTime,x,y,z,iParticle,_atoms[group])

class pdbParser:
    def __init__(self):
        self.index=0
        self.text=""
    
    def addAtom(self,particle, x, set, time):
        self.text+= _pdbAtomRow(index=self.index,iParticle=particle,x=x[0],y=x[1],z=x[2],group=set,iTime=time) + "\n"
        self.index+=1
                
    
    def __str__(self):
        return self.text


def pbc( x , lBox ):
    lBoxInv=1./lBox
    return x - np.floor(x*lBoxInv + 0.5)*lBox


def convertHDF5ToPdb( fileNameHDF5, fileNamePDB ,lBox=None ):
    file=h5py.File(fileNameHDF5,"r")

    parser=pdbParser()
    iSet=0
    for iSet,setName in enumerate(file):
        set=file[setName]
        particleData=set["particles"]
        nexts=set["nexts"]
        prevs=set["prevs"]
        heads=set["heads"]
        tails=set["tails"]

        _maxX=[9999,9999,9999]

        def isMasked(i,t):
            if (t>heads[i]) or (t<=tails[i]):
                    return True 


        dimensions=3
        
        X=np.array(particleData[:,:,:])

        if ( lBox is not None ):
            for d in range(dimensions):
                X[:,d,:]=pbc(X[:,d,:],lBox[d])

        for i in range(particleData.shape[2]):
            for t in range(particleData.shape[0]):
                x=[ X[t,d,i] for d in range(dimensions) ]


                if isMasked(i,t):
                    x=_maxX
                parser.addAtom(particle=i,x=x,time=t,set=iSet)

        
    with open(fileNamePDB,"w+") as f:
        f.write(str(parser))

    file.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert HDF5 output file to pdb file ')
    parser.add_argument('inputFiles', type=str, nargs='+',
                    help='File name of the HDF5 file to convert')
    parser.add_argument('--out', type=str ,
                        help='File name of the convert pdb file',default=None, nargs='*')
        
    parser.add_argument('--lBox', type=float ,
                        help='Size of the box',default=None, nargs=3 )
    
    args = parser.parse_args()

    inputFiles=args.inputFiles
    outputFiles=args.out
    lBox=args.lBox


    # generate output file names if --out was not given
    if outputFiles is None:
       outputFiles=[ re.sub("\.[a-zA-Z0-9]+$",".pdb" , inputFile )  for inputFile in inputFiles ]
    
    for inputFile,outputFile in tqdm.tqdm(zip(inputFiles,outputFiles) , total=len(outputFiles) ):
        
        convertHDF5ToPdb(inputFile,outputFile,lBox=args.lBox)

    