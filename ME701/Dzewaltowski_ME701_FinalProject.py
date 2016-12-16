#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:26:26 2016

@author: Alexdzewaltowski
"""

from __future__ import division 
import numpy as np
import glob
from scipy.integrate import dblquad
import matplotlib.pyplot as plt

Cp=[0,0,0,0,0,0]
Filenames=[0,0,0,0,0,0]
ListTorque=[0,0,0,0,0,0]

#user input
n=6
HubRadius=4/100
WindSpeed=11
B=3
theta=70
R=26/100
MaxChord=6/100
MinChord=3/100
rho=0.001225
#Reads lines from Xfoil files
for filename in glob.iglob('/Users/Alexdzewaltowski/Desktop/ME701/ProfileData/*.txt', recursive=True):
    f=open(filename, 'r')
    file=f.readlines()[12:]
    for i in range(len(file)):
        lines=file[i].split()
        phi=(eval(lines[0])+theta)*np.pi/180
        CL=eval(lines[1])
        CD=eval(lines[2])
        incrementChord=(MaxChord-MinChord)/n
        
        if CL>0 and CD>0:
            for k in range(0,n):
                    #Bounds
                    rLower=R-(R-HubRadius)*k/n-0.01
                    rUpper=R-(R-HubRadius)*(k+1)/n
                    cLower=MinChord+incrementChord*(1+k)
                    cUpper=MinChord+incrementChord*(k)
                    
                    #Calculation for maximum Tip Speed Ratio
                    sigma=B*MinChord/(2*np.pi*R)
                    X=4*np.sin(phi)**2/(sigma*CL*np.cos(phi))
                    alpha=1/(1+X)
                    X=4*np.cos(phi)/(sigma*CL)
                    alphaprime=1/(X-1)
                    LTSR=np.tan(phi)*(1-alpha)/(1+alphaprime)
                    

                    
                    r1,r2= rLower, rUpper
                    c1, c2= cLower, cUpper
                    
                    #had to integrat over the equation in terms of r and c.....
                    if 6<LTSR<8 and phi<90*np.pi/180:
                        def CPfunc(r,c,B,phi,CL,CD,R,LTSR):
                            return (8/LTSR**2)*((LTSR*r/R)**3)*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))\
                                                      *1/(4*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))\
                                                      *(np.cos(phi)/(CL*B*c/(2*np.pi*r)))-1)\
                                                      *(1-(1/(1+(4*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))\
                                                      *(np.sin(phi)**2)/((B*c*CL*np.cos(phi))/(2*np.pi*r))))))\
                                                      *(1-(CD*(np.tan(phi)/CL)))
        
                        def torque(r,c,B,WindSpeed,phi,CL,R,rho) :
                           return -1*4*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))*np.pi*rho*WindSpeed**2\
                                   *(1/(1+(4*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))\
                                   *(np.sin(phi)**2)/((B*c*CL*np.cos(phi))/(2*np.pi*r)))))*(1-\
                                   (1/(1+(4*((2/np.pi)*np.arccos(np.exp((-B/2)*(1-(r/R))/(r*np.sin(phi)/R))))\
                                   *(np.sin(phi)**2)/((B*c*CL*np.cos(phi))/(2*np.pi*r))))))*r
                        
                        Torque= dblquad(torque, c1, c2, lambda c:r1, lambda c: r2, args=(B, WindSpeed, phi, CL, R, rho))                                  
                        Torque=Torque[0]                
                                          
                        CP1=dblquad(CPfunc, c1, c2, lambda c: r1, lambda c: r2, args=(B,phi,CL,CD,R,LTSR))
                        CP=CP1[0]

                        if Cp[k]<CP:
                            Cp[k]=CP
                            Filenames[k]=filename
                            ListTorque[k]=Torque 

Power=[]
x=[]
for i in range (0,20):
    x.append(i)
    Power.append(sum(Cp)*np.pi*R**2*rho*i**3)     
                    
plt.plot(x,Power)
plt.xlabel("Wind Speed (m/s)")
plt.ylabel("Watts W")
plt.show()









