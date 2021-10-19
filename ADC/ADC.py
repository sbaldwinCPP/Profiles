# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 16:21:43 2021

Atmospheric Dispersion Comparability (ADC)

Started from scratch to generate approach profile plots for ADC report

@author: sbaldwin

3 defined approaches:
R1: T1=100, T2=0, S1=725, S2=0 ,    M=10,   B=10     ID: 300 SF @ Zo 0.1
R2: T1=300, T2=0, S1=0,   S2=922 ,  M=75,   B=75     ID: 180 SF @ Zo 1.0
R3: T1=400, T2=0, S1=0,   S2=715 ,  M=50,   B=50     ID: 180 SF @ Zo 0.5


"""

#%% Import
print('Importing modules...')

#default packages
import os

#non-default packages
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

print('Modules imported!')

#%% Get files & data
def ingest_data():
    
    inidir=os.getcwd()  #start at location of this script
    
    E_path= inidir + '\EAST_ADC.xlsx'
    W_path= inidir + '\WEST_ADC.xlsx'
    
    print('Reading data from WTE...')
    E_raw=pd.read_excel(E_path, header=None,skiprows=4)
    print('Reading data from WTW...')
    W_raw=pd.read_excel(W_path, header=None,skiprows=4)
    
    
    #inital filtering 
    #inculde any filters that rely on a column not needed for rest of script
    #other filters can be applied later
    print('Pre-processing raw data...')
    
    #remove calcheck runs, reference comments col
    E_f=E_raw[pd.isnull(E_raw[4])]
    W_f=W_raw[pd.isnull(W_raw[4])] 
    
    #copy select data
    #West and East files are not 1:1
    E=E_f[[10,16,17,18,190,191,230,231,233,235,374,382]].copy()
    W=W_f[[10,16,17,18,187,188,227,228,230,232,337,345]].copy()
    
    columns=['uref',
             'x',
             'y',
             'z',
             'uhf',
             'varhf',
             'Uomni',
             'Uvar',
             'Vvar',
             'Wvar',
             'trip',
             'blocks'
             ]
    
    E.columns=columns
    W.columns=columns
    
    #make all values numeric, some read in as text
    E=E.apply(pd.to_numeric)
    W=W.apply(pd.to_numeric)
    
    #remove low WS data (zeroes)
    E=E[E.uref>5] 
    W=W[W.uref>5]
    
    E=E.round({"trip":0,"blocks":0})
    W=W.round({"trip":0,"blocks":0})
    
    E['hfnorm']=E.uhf/E.uref
    E['omninorm']=E.Uomni/E.uref
    
    W['hfnorm']=W.uhf/W.uref
    W['omninorm']=W.Uomni/W.uref
    
    #remove bad data (extremely high u/urefs)
    E=E[(E.hfnorm < 1.5) & (E.omninorm < 1.5)]
    W=W[(W.hfnorm < 1.5) & (W.omninorm < 1.5)]
    
    #calc TI
    E['hfti']=(E.varhf**(1/2))/E.uhf
    E['Uti']=(E.Uvar**(1/2))/E.Uomni
    E['Vti']=(E.Vvar**(1/2))/E.Uomni
    E['Wti']=(E.Wvar**(1/2))/E.Uomni
    
    
    W['hfti']=(W.varhf**(1/2))/W.uhf
    W['Uti']=(W.Uvar**(1/2))/W.Uomni
    W['Vti']=(W.Vvar**(1/2))/W.Uomni
    W['Wti']=(W.Wvar**(1/2))/W.Uomni
    
    #normalize Z height
    E['znorm']=E.z/1000
    W['znorm']=W.z/1000
    
    print('Initial data ingested.')
    return E,W

E,W=ingest_data()
    
#%% functions from CPP AQ "Bible" and ADC report

#   velocity equations
#   from Counihan (1975)
def nFromZo(Zo):    
    return 0.096*np.log10(Zo)+0.016*(np.log10(Zo)**2)+.24

def pwr(z,zref,Zo):    #   u/Uref = (z/Zo)**n
    n=nFromZo(Zo)
    return (z/zref)**n

#   new log law from wikipedia
#   Holmes JD. Wind Loading of Structures. 3rd ed. Boca Raton, Florida: CRC Press; 2015.
def log(z,zref,Zo):
    d=0
    #d=  # zero-plane displacement:
         # height in meters above the ground at which zero wind speed is 
         # achieved as a result of flow obstacles such as trees or buildings
    return np.log((z-d)/Zo)/np.log((zref-d)/Zo)

def Ustar(u,z,Zo):
    k=.4
    d=0
    return u*k/(np.log(z-d)/Zo)

#   turbulence equation
#   EPA (1981)
#   valid between 5 and 100 m
#   above 100m, decreases linearly to 0.01 at 600 m
def TI(z,Zo):
    n=nFromZo(Zo)
    return n*(np.log(30/Zo))/(np.log(z/Zo))


#%% plot functions

def ADC(tun,R,speed):
    #   for debugging #  tun,R,speed= 'E','R1','fast'
    
    #   determine WT, SF , Zo target
    if tun=='E': df=E
    elif tun=='W': df=W
    else: raise ValueError('Tunnel designation not recognized, use "E" or "W"')
    # filter for tunnel speed
    if speed=='fast':   df=df[df.uref>8]    
    elif speed=='slow': df=df[df.uref<8]
    else: raise ValueError('Speed designation not recognized, use "fast" or "slow"')
    #filter for roughness    
    if R=='R1':     
        df=df[df.trip==100]
        SF=300 
        Zo=0.1
    elif R=='R2':   
        df=df[df.trip==300]
        SF=180 
        Zo=1.0
    elif R=='R3':   
        df=df[df.trip==400]
        SF=180 
        Zo=0.5
    else: raise ValueError('Roughness designation not recognized, use "R1", "R2", or "R3"')
     
    center= df[(df.y>-10) & (df.y <10)]

    X0=center[(center.x>-10) & (center.x<10)]
    X1=center[(center.x>780) & (center.x<800)]
    X2=center[(center.x>1580) & (center.x<1600)]
    X3=center[(center.x>2380) & (center.x<2400)]
    
    #   full-scale x positions    
    x1=round(790*SF/1000)
    x2=round(1590*SF/1000)
    x3=round(2390*SF/1000)
    
    #   "Fig B-2" #change all hfnorm to omninorm to plot 12hp data
    plt.scatter(X0.hfnorm, X0.znorm, marker='o', facecolors='none', edgecolors='k', label='x=0 m')
    plt.scatter(X1.hfnorm, X1.znorm, marker='s', facecolors='none', edgecolors='k', label='x={} m'.format(x1))
    plt.scatter(X2.hfnorm, X2.znorm, marker='^', facecolors='none', edgecolors='k', label='x={} m'.format(x2))
    plt.scatter(X3.hfnorm, X3.znorm, marker='x', color='k', label='x={} m'.format(x3))
    
    plt.legend()
    plt.xlabel('U/Uref')
    plt.ylabel('Z/Zref')
    
    title='MEAN VELOCITY PROFILES\n{} Tunnel ADC, 1:{}, Zo={}, y=0 m'.format(tun,SF,Zo)
    plt.suptitle(title)
    
    savepath=os.path.join(os.getcwd(),'Plots')
    if not os.path.exists(savepath): os.makedirs(savepath)
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B2'))
    plt.close()
    
    
    #   "Fig B-3a"
    fig,ax = plt.subplots(1,2, sharey=True)
    
    #    hotfilm data
    hf=ax[0].scatter(X0.hfnorm, X0.znorm, marker='s', facecolors='none', edgecolors='k', label='HF')
    #   omni data
    omni=ax[0].scatter(X0.omninorm, X0.znorm, marker='o', facecolors='none', edgecolors='k', label='12HP')
    #   target Zo power law line
    z_pwr=np.linspace(.01,1.2,num=1000)
    ax[0].plot(pwr(z_pwr,1,Zo), z_pwr, c='k', label='Target U (pwr, n={})'.format(round(nFromZo(Zo),2)))
    
    #### log law - has to be corrected for FS calculations
    u=X0.uref.mean()
    ustar=Ustar(u,SF,Zo)
    ax[0].plot(log(z_pwr*SF,SF,Zo), z_pwr, c='k', linestyle='dashed', label='Target U (log, u*={})'.format(round(ustar,3)))
    ####
    
    ax[0].set_title('MEAN PROFILE PLOT')
    ax[0].legend()
    
    #   turbulence
    ax[1].scatter(X0.Uti, X0.znorm, marker='s', facecolors='none', edgecolors='k', label="U'/U")
    ax[1].scatter(X0.Vti, X0.znorm, marker='^', facecolors='none', edgecolors='k', label="V'/U")
    ax[1].scatter(X0.Wti, X0.znorm, marker='o', facecolors='none', edgecolors='k', label="W'/U")
    ax[1].scatter(X0.hfti, X0.znorm, marker='x', color='k', label="Uhf'/U")
    
    
    z_ti=np.linspace(5,100,num=1000)
    ax[1].plot(TI(z_ti,Zo),z_ti/SF,'--',c='k', label='Target T.I.')
    
    
    ax[1].legend()
    ax[1].set_title('TURBULENCE PROFILE PLOT')
    
    ax[0].set_ylabel('Height, Z/Zref')
    ax[0].set_xlabel('Mean Velocity, U/Uref')
    ax[1].set_xlabel("Local Turbulent Intensity, U'/U")
    
    fig.suptitle('x=0 m')
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B3_a'))
    plt.close()
    
    #   "Fig B-3b"
    fig,ax = plt.subplots(1,2, sharey=True)
    
    #    hotfilm data
    hf=ax[0].scatter(X1.hfnorm, X1.znorm, marker='s', facecolors='none', edgecolors='k')
    #   omni data
    omni=ax[0].scatter(X1.omninorm, X1.znorm, marker='o', facecolors='none', edgecolors='k')
    #   target Zo power law line
    z_pwr=np.linspace(.01,1.2,num=1000)
    ax[0].plot(pwr(z_pwr,1,Zo), z_pwr, c='k')
    #fake artist for legend
    targ=mpl.lines.Line2D([],[],color='k')
    
    ax[0].legend([hf,omni,targ],['HF','12HP','Power law, n={}'.format(round(nFromZo(Zo),3))])
    ax[0].set_title('MEAN PROFILE PLOT')
    
    #   turbulence
    ax[1].scatter(X1.Uti, X1.znorm, marker='s', facecolors='none', edgecolors='k', label="U'/U")
    ax[1].scatter(X1.Vti, X1.znorm, marker='^', facecolors='none', edgecolors='k', label="V'/U")
    ax[1].scatter(X1.Wti, X1.znorm, marker='o', facecolors='none', edgecolors='k', label="W'/U")
    ax[1].scatter(X1.hfti, X1.znorm, marker='x', color='k', label="Uhf'/U")
    
    
    z_ti=np.linspace(5,100,num=1000)
    ax[1].plot(TI(z_ti,Zo),z_ti/SF,'--',c='k', label='Target T.I.')
    
    
    ax[1].legend()
    ax[1].set_title('TURBULENCE PROFILE PLOT')
    
    ax[0].set_ylabel('Height, Z/Zref')
    ax[0].set_xlabel('Mean Velocity, U/Uref')
    ax[1].set_xlabel("Local Turbulent Intensity, U'/U")
    
    fig.suptitle('x={} m'.format(x1))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B3_b'))
    plt.close()
    
    #   "Fig B-3c"
    fig,ax = plt.subplots(1,2, sharey=True)
    
    #    hotfilm data
    hf=ax[0].scatter(X2.hfnorm, X2.znorm, marker='s', facecolors='none', edgecolors='k')
    #   omni data
    omni=ax[0].scatter(X2.omninorm, X2.znorm, marker='o', facecolors='none', edgecolors='k')
    #   target Zo power law line
    z_pwr=np.linspace(.01,1.2,num=1000)
    ax[0].plot(pwr(z_pwr,1,Zo), z_pwr, c='k')
    #fake artist for legend
    targ=mpl.lines.Line2D([],[],color='k')
    
    ax[0].legend([hf,omni,targ],['HF','12HP','Power law, n={}'.format(round(nFromZo(Zo),3))])
    ax[0].set_title('MEAN PROFILE PLOT')
    
    #   turbulence
    ax[1].scatter(X2.Uti, X2.znorm, marker='s', facecolors='none', edgecolors='k', label="U'/U")
    ax[1].scatter(X2.Vti, X2.znorm, marker='^', facecolors='none', edgecolors='k', label="V'/U")
    ax[1].scatter(X2.Wti, X2.znorm, marker='o', facecolors='none', edgecolors='k', label="W'/U")
    ax[1].scatter(X2.hfti, X2.znorm, marker='x', color='k', label="Uhf'/U")
    
    
    z_ti=np.linspace(5,100,num=1000)
    ax[1].plot(TI(z_ti,Zo),z_ti/SF,'--',c='k', label='Target T.I.')
    
    
    ax[1].legend()
    ax[1].set_title('TURBULENCE PROFILE PLOT')
    
    ax[0].set_ylabel('Height, Z/Zref')
    ax[0].set_xlabel('Mean Velocity, U/Uref')
    ax[1].set_xlabel("Local Turbulent Intensity, U'/U")
    
    fig.suptitle('x={} m'.format(x2))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B3_c'))
    plt.close()
    
    #   "Fig B-3d"
    fig,ax = plt.subplots(1,2, sharey=True)
    
    #    hotfilm data
    hf=ax[0].scatter(X3.hfnorm, X3.znorm, marker='s', facecolors='none', edgecolors='k')
    #   omni data
    omni=ax[0].scatter(X3.omninorm, X3.znorm, marker='o', facecolors='none', edgecolors='k')
    #   target Zo power law line
    z_pwr=np.linspace(.01,1.2,num=1000)
    ax[0].plot(pwr(z_pwr,1,Zo), z_pwr, c='k')
    #   fake artist for legend
    targ=mpl.lines.Line2D([],[],color='k')
    
    ax[0].legend([hf,omni,targ],['HF','12HP','Power law, n={}'.format(round(nFromZo(Zo),3))])
    ax[0].set_title('MEAN PROFILE PLOT')
    
    #   turbulence
    ax[1].scatter(X3.Uti, X3.znorm, marker='s', facecolors='none', edgecolors='k', label="U'/U")
    ax[1].scatter(X3.Vti, X3.znorm, marker='^', facecolors='none', edgecolors='k', label="V'/U")
    ax[1].scatter(X3.Wti, X3.znorm, marker='o', facecolors='none', edgecolors='k', label="W'/U")
    ax[1].scatter(X3.hfti, X3.znorm, marker='x', color='k', label="Uhf'/U")
    
    
    z_ti=np.linspace(5,100,num=1000)
    ax[1].plot(TI(z_ti,Zo),z_ti/SF,'--',c='k', label='Target T.I.')
    
    
    ax[1].legend()
    ax[1].set_title('TURBULENCE PROFILE PLOT')
    
    ax[0].set_ylabel('Height, Z/Zref')
    ax[0].set_xlabel('Mean Velocity, U/Uref')
    ax[1].set_xlabel("Local Turbulent Intensity, U'/U")
    
    fig.suptitle('x={} m'.format(x3))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B3_d'))
    plt.close()
    
    #   lateral profiles
    lat0=df[(df.x>-10) & (df.x<10)]
    lat1=df[(df.x>780) & (df.x<800)]
    lat2=df[(df.x>1580) & (df.x<1600)]
    lat3=df[(df.x>2380) & (df.x<2400)]
    
    dh=5
    
    h1=.75*30/SF*1000
    h2=30/SF*1000
    h3=1.25*30/SF*1000
    
    lat0_h1=lat0[(lat0.z>h1-dh)&(lat0.z<h1+dh)]
    lat0_h2=lat0[(lat0.z>h2-dh)&(lat0.z<h2+dh)]
    lat0_h3=lat0[(lat0.z>h3-dh)&(lat0.z<h3+dh)]
    
    lat1_h1=lat1[(lat1.z>h1-dh)&(lat1.z<h1+dh)]
    lat1_h2=lat1[(lat1.z>h2-dh)&(lat1.z<h2+dh)]
    lat1_h3=lat1[(lat1.z>h3-dh)&(lat1.z<h3+dh)]
    
    lat2_h1=lat2[(lat2.z>h1-dh)&(lat2.z<h1+dh)]
    lat2_h2=lat2[(lat2.z>h2-dh)&(lat2.z<h2+dh)]
    lat2_h3=lat2[(lat2.z>h3-dh)&(lat2.z<h3+dh)]
    
    lat3_h1=lat3[(lat3.z>h1-dh)&(lat3.z<h1+dh)]
    lat3_h2=lat3[(lat3.z>h2-dh)&(lat3.z<h2+dh)]
    lat3_h3=lat3[(lat3.z>h3-dh)&(lat3.z<h3+dh)]
    
    #   lateral velocity plots
    #   "Fig B-4a"
    plt.scatter(lat0_h1.y*SF/1000, lat0_h1.hfnorm, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat0_h2.y*SF/1000, lat0_h2.hfnorm, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat0_h3.y*SF/1000, lat0_h3.hfnorm, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,1)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel('U/Uref')
    plt.suptitle('LATERAL VELOCITY PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,0))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B4_a'))
    plt.close()
    
    #   "Fig B-4b"
    plt.scatter(lat1_h1.y*SF/1000, lat1_h1.hfnorm, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat1_h2.y*SF/1000, lat1_h2.hfnorm, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat1_h3.y*SF/1000, lat1_h3.hfnorm, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,1)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel('U/Uref')
    plt.suptitle('LATERAL VELOCITY PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x1))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B4_b'))
    plt.close()
    

    #   "Fig B-4c"
    plt.scatter(lat2_h1.y*SF/1000, lat2_h1.hfnorm, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat2_h2.y*SF/1000, lat2_h2.hfnorm, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat2_h3.y*SF/1000, lat2_h3.hfnorm, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,1)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel('U/Uref')
    plt.suptitle('LATERAL VELOCITY PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x2))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B4_c'))
    plt.close()
    
    #   "Fig B-4d"
    plt.scatter(lat3_h1.y*SF/1000, lat3_h1.hfnorm, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat3_h2.y*SF/1000, lat3_h2.hfnorm, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat3_h3.y*SF/1000, lat3_h3.hfnorm, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,1)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel('U/Uref')
    plt.suptitle('LATERAL VELOCITY PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x3))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B4_d'))
    plt.close()    
    
    
    #   lateral turbulance plots
    #   "Fig B-5a"
    plt.scatter(lat0_h1.y*SF/1000, lat0_h1.hfti, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat0_h2.y*SF/1000, lat0_h2.hfti, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat0_h3.y*SF/1000, lat0_h3.hfti, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,.6)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel("U'/U")
    plt.suptitle('LATERAL TURBULENCE PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,0))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B5_a'))
    plt.close()
    
    #   "Fig B-5b"
    plt.scatter(lat1_h1.y*SF/1000, lat1_h1.hfti, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat1_h2.y*SF/1000, lat1_h2.hfti, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat1_h3.y*SF/1000, lat1_h3.hfti, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,.6)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel("U'/U")
    plt.suptitle('LATERAL TURBULENCE PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x1))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B5_b'))
    plt.close()
    
    #   "Fig B-5c"
    plt.scatter(lat2_h1.y*SF/1000, lat2_h1.hfti, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat2_h2.y*SF/1000, lat2_h2.hfti, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat2_h3.y*SF/1000, lat2_h3.hfti, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,.6)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel("U'/U")
    plt.suptitle('LATERAL TURBULENCE PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x2))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B5_c'))
    plt.close()
    
    #   "Fig B-5d"
    plt.scatter(lat3_h1.y*SF/1000, lat3_h1.hfti, marker='s', facecolors='none', edgecolors='k', label='0.75 h')
    plt.scatter(lat3_h2.y*SF/1000, lat3_h2.hfti, marker='o', facecolors='none', edgecolors='k', label='h')
    plt.scatter(lat3_h3.y*SF/1000, lat3_h3.hfti, marker='^', facecolors='none', edgecolors='k', label='1.25 h')
    
    plt.ylim(0,.6)
    plt.legend()
    plt.xlabel('y (m), Lateral Distance')
    plt.ylabel("U'/U")
    plt.suptitle('LATERAL TURBULENCE PROFILES \n{} Tunnel ADC, 1:{}, Zo={}, h=30 m, x={} m'.format(tun,SF,Zo,x3))
    
    plt.savefig('{}\{}-{}-{}-{}.png'.format(savepath,tun,R,speed,'Fig_B5_d'))
    plt.close()
    
    #END
    
#%% call function to generate plots

# =============================================================================
# #set plot style 
# mpl.rcdefaults()            #reset to defaults
# styles=plt.style.available  #show plot styles
# #plt.style.use(styles[9])
# =============================================================================

ADC('E','R1','fast')
ADC('W','R1','fast')

ADC('E','R2','fast')
ADC('W','R2','fast')

ADC('E','R3','fast')
ADC('W','R3','fast')

# =============================================================================
# ADC('E','R1','slow')
# ADC('W','R1','slow')
# 
# ADC('E','R2','slow')
# ADC('W','R2','slow')
# 
# ADC('E','R3','slow')
# ADC('W','R3','slow')
# =============================================================================
