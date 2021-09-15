import argparse
import subprocess
import re
import sys
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import tqdm
import seaborn as sns

sns.set_style("whitegrid")

def timePerf(command,n=10):

    labels=[]
    timings=[]
    for i in tqdm.tqdm(range(n)):

        output= subprocess.check_output([command],shell=True).decode(sys.stdout.encoding)
    
        expr = re.compile( r"Time\t(.*)\t(\d+\.\d+)") 
        ms=expr.findall(output)
        labels+=( [label for label,_ in ms ])
        timings+=([float(time) for _,time in ms ])

    data=pd.DataFrame({"label":labels , "time":timings})


    return data

def analPerformance(datas,bins=None):
    for label,data in datas.groupby("label"):
       
        time = np.array(data["time"])
        currentBins=bins
        if currentBins is None:
            currentBins=int( np.max( [ np.sqrt(len(time)) , 10 ] ) )
        
        hist=np.histogram(time,bins=currentBins,density=True)

        print ("{}\t{}\t{}".format(label,np.average(time),np.sqrt(np.var(time)) ) )


        y=np.array(hist[0])
        x=np.array(hist[1][0:len(hist[0])])
        dx=x[1]-x[0]

        
        #plt.plot(x,y,"--",label=label,linewidth=1)

        plt.bar(x,y,align="edge",label=label,alpha=0.9,width=dx)

        


        



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Time profiling of a command')

    parser.add_argument('command', help='command to time')
    
    args  = parser.parse_args()
    data=timePerf(args.command,n=100)
    #print(data)
    analPerformance(data,bins=10)
    plt.legend()
    plt.show()

