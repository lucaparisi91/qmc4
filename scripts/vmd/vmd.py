import re
import argparse
import os
import os.path


_templateLocation="/home/luca/source/qmc4/PPA/qmc4/scripts/vmd/loadConfigurationTemplate.tcl"


def generateVMDScript( pdbFiles, M, lBox=None):

    with open(_templateLocation ) as fTemplate:
        tclScript=fTemplate.read()

    tclScript=re.sub("__FIRSTFILE__",pdbFiles[0],tclScript)

    tclScript=re.sub("__M__",str(M),tclScript)


    fileString=" ".join(pdbFiles[1:] )

    tclScript=re.sub("__PDBFILES__",fileString,tclScript)

    pbcScript=""
    if lBox is not None:
        pbcScript='''
pbc set { __LBOX__ } -all
pbc wrap -all
pbc join connected -all
pbc box_draw
    '''

        lBoxString= " ".join( [ str(l) for l in lBox]  )
        pbcScript=re.sub("__LBOX__",lBoxString,pbcScript)
    
    tclScript=re.sub("__PBCSCRIPT__",pbcScript,tclScript)

    return tclScript



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate a TCL script for vmd visualization of the molecule')
    parser.add_argument('inputFiles', type=str, nargs='+',
                    help='File names of the pdb files to convert')
    parser.add_argument('--M', type=int,
                    help='Number of beads in each simulation')
    
    parser.add_argument('--out', type=str,
                    help='Number of beads in each simulation',default="vmdTclScript.tcl")
    parser.add_argument('--box', type=float ,
                        help='Size of the box',default=None, nargs=3 )
    
    
    args = parser.parse_args()

    inputFiles= [ os.path.abspath(file) for file in args.inputFiles ]
    
    script=generateVMDScript(  inputFiles , args.M ,lBox=args.box)

    with open(args.out,"w") as f:
        f.write(script)
    




