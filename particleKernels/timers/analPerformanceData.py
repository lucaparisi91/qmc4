import pandas as pd 
import matplotlib.pylab as plt
import seaborn as sns
from scipy.optimize import curve_fit
import numpy as np


if __name__ == "__main__":
    
    datas=pd.read_csv("2bDistancesTimings_2.dat",sep=" ",index_col=False)
    print(datas)
    
    cutOff=0;
    for hue,data in datas.groupby("method"):

        plt.plot(data["N"] , data["time"], "o",label=hue )
        
        f=lambda x,a,b : a + b*x
        x=data["N"]
        y=data["time"]
        y=y[x>cutOff]
        x=x[x>cutOff]
    
        params,sigma=curve_fit(f,np.log(x),np.log(y))
        errors=np.sqrt(np.diag(sigma))

        x=np.linspace( np.min(x),np.max(x) , num=1000 )
        plt.plot(x,np.exp(f(np.log(x),*params) ) , "--")
        print("{} {}".format(hue,params))

        





    plt.xscale("log")
    plt.yscale("log")
    
    plt.legend()
    plt.show()
    