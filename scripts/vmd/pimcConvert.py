import numpy as np
import scipy as sp
import pandas as pd
import h5py
import argparse

def _pdbAtomRow( index,iParticle, x,y,z,group,iTime ):
    _atoms=["N","O"]
    template="ATOM  {:5d} {}    HIS A{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f} 00.00           {}"
    return template.format(index,_atoms[group],iTime,x,y,z,iParticle,_atoms[group])

class pdbParser:
    def __init__(self):
        self.index=0
        self.text=""
    def addAtoms(self, atomData,group=0):
        for i in range(atomData.shape[2]):
            for t in range(atomData.shape[0]):
                self.text+= _pdbAtomRow(index=self.index,iParticle=i,x=atomData[t,0,i],y=atomData[t,1,i],z=atomData[t,2,i],group=group,iTime=t) + "\n"
                self.index+=1
    
    
    def __str__(self):
        return self.text
            



def convertHDf5ToPdb(hdf5FileName,pdbFileName):
    fromFile=h5py.File(source,"r")
    with open(pdbFileName,"w+") as f:
        f.write("")
    fromFile.close()
    



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert HDF5 output file to pdb file ')
    parser.add_argument('fileNameHDF5', type=str,
                    help='File name of the HDF5 file to convert')
    parser.add_argument('fileNamePDB', type=str,
                    help='File name of the convert pdb file')
    
    args = parser.parse_args()

    
    fileNameHDF5=args.fileNameHDF5
    fileNamePDB=args.fileNamePDB
    file=h5py.File(fileNameHDF5,"r")
    
    parser=pdbParser()
    
    data=np.array(file["set0"]["particles"])
    parser.addAtoms(data)

    with open(fileNamePDB,"w+") as f:
        f.write(str(parser))

    file.close()
