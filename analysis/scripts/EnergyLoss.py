import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

"""
There are 6 data sets in the output of lise++'s dE/dX tables. With the energy values in between. 

0-[He-base] Hubert 
1-[H-base] Ziegler (low energy)
2-ATIMA 1.2 LS-theory (high energy)
3-ATIMA 1.2 w/o LS-correction
4-ATIMA 1.4 improved mean charge formula
5-Ziegler electric component
6-Ziegler nuclear component

Make sure the tables are formated with the columns in the above order.
Another annoying part of the 
"""
def LoadLayer(inFile,u):
    df = pd.read_csv(inFile, delimiter = "\t",header=1)
    df.columns = ["e0","He-base","e1","H-base","e2","atima-1.2-LS","e3","atima-1.2","e4","atima-1.4","e5","ziegler-el","e6","ziegler-nuc","end"]
    f = interp1d(df['e0']*u, df['H-base'], kind='cubic')
    return f

def main():
    cm2um=10000
    l=16*cm2um # in
    u=25
    
    # read data, y-axis is in MeV/micrometer, x-axis is in MeV/u
    gas=LoadLayer("../input/Al25_in_Butane.txt",25)
    kapton=LoadLayer("../input/Al25_in_Kapton.txt",25)
    target=LoadLayer("../input/Al25_in_CD2.txt",25)

    beam=np.linspace(60,102,100)
    
    y1=beam-target(beam)*0.516/2
    y2=y1-kapton(y1)*7
    y3=y2-gas(y2)*l

    #print("Energy loss = ",y3)
    plt.plot(y3,beam)

    coef=np.polyfit(y3,beam,2)
    poly=np.poly1d(coef)
    xfit=np.linspace(y3[0],y3[-1])
    yfit=poly(xfit)
    plt.plot(xfit,yfit)
    plt.show()

    print(coef)
    

if __name__ == "__main__":
    main()

