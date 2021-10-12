# -*- coding: utf-8 -*-
"""
Created on Tue Oct 5 14:16:23 2021

@author: sbaldwin
"""

#%% Import
print('Importing modules...')

try:
    #default packages
    import sys
    import os
    import datetime
    #non default packages
    import numpy as np
    import pandas as pd
    import easygui
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    
    print('All modules imported')

except ModuleNotFoundError as err:
    input(str(err) +' \n Press enter to exit')
    sys.exit()

#%% Functions

# files
def Get_Data():
    print('Use GUI to select data file...')
    inidir = os.getcwd() #starts file selection at working directory
    txt='Select data file...'
    ftyp= '*.xlsx'
    dflt=inidir + "\\" + ftyp
    filepath = easygui.fileopenbox(default=dflt,msg=txt,filetypes=ftyp,multiple=False)
    if filepath is None:
        easygui.msgbox(msg='Nothing selected, nothing to do. Goodbye!')
        sys.exit()
        
    return filepath

def Inputs():
    print('Use GUI to select inputs...')
    msg = "Enter Settings"
    title = "Enter Settings"
    fieldNames = ["Tunnel:", "Scale:", "Target Zo:"]
    fieldValues = ['GDW', 180, 0.5]  #Defaults
    inputs=easygui.multenterbox(msg, title, fieldNames, fieldValues)
    inputs= inputs[0], int(inputs[1]), float(inputs[2])
    return inputs

# og dataframe filters
def F_f(df,fetch):
    return df[df.Fetch==fetch].copy()

def F_t(df,T):
    """
    'GDW' or 'GDE'
    """
    return df[df.Tunnel==T].copy()

def F_d(df,D):
    """
    'Vertical Profile' or 'Lateral Profile'
    """
    return df[df.Direction==D].copy()


#%% functions from CPP AQ "Bible" and ADC report

#   velocity equations
#   from Counihan (1975)
def nFromZo(Zo):    
    return 0.096*np.log10(Zo)+0.016*(np.log10(Zo)**2)+.24

def pwr(z,zref,Zo):    #   u/Uref = (z/Zo)**n
    n=nFromZo(Zo)
    return (z/zref)**n

#   turbulence equation
#   EPA (1981)
#   valid between 5 and 100 m
#   above 100m, decreases linearly to 0.01 at 600 m
def TI(z,Zo):
    n=nFromZo(Zo)
    return n*(np.log(30/Zo))/(np.log(z/Zo))*100

#%% Log functions from RP's matlab code
# not in use (yet)
def log(z,Zo):    #   u/Uinf = 2.5*(Ustar0/Uinf)*ln(z/Zo)
    Fr=F_ratio(Zo)
    return 2.5*Fr*np.log(z/Zo)

def F_ratio(Zo):
    return np.sqrt(2.75e-3+6e-4*np.log10(Zo))

#%% plotting function
# mostly for QA
def Profile_plot(fetch, tun):
    
    data=F_f(df_og,fetch)
    
    main=data['Main(mm)'].mean()
    boost=data['Boost(mm)'].mean()
    
    block_ht=max(main,boost)*SF/1000
    
    data=F_t(data,tun)
    
    data=F_d(data,'Vertical Profile')
    
    data=data[data.Location.isin([
        'Ref',
        '-3r/2',
        '-r/2',
        'TTC',
        '+r/2',        
            ])]
    
    data['Z_fs']=data['Z(mm)']*SF/1000
    
    Z = np.linspace(block_ht,Zmax,num=100)
    
    #actual plotting, all above is just data manipulation/prep
    
    figure,ax = plt.subplots(1,2, sharey=True)
    
    # velocity
    data.plot(kind='scatter', x='Uhf/Uref', y='Z_fs',ax=ax[0], marker='x', c='X(mm)', cmap=cm1, colorbar=False, label='HF')
    data.plot(kind='scatter', x='Uomni/Uref', y='Z_fs',ax=ax[0], marker='+', c='X(mm)', cmap=cm2, colorbar=False, label='12HP')
    
    # turbulence
    data.plot(kind='scatter', x='TIUhf(%)', y='Z_fs',ax=ax[1], marker='x', c='X(mm)', cmap=cm1, colorbar=False, label='HF')
    data.plot(kind='scatter', x='TIUomni(%)', y='Z_fs',ax=ax[1], marker='+', c='X(mm)', cmap=cm2, colorbar=True, label='12HP')
    
    # targets
    ax[0].plot(pwr(Z,SF,Zo_targ), Z, c='k', label='Target U')
    ax[1].plot(TI(Z,Zo_targ), Z, c='k', label='Target TI')
    
    # horizontal lines
    ax[1].axhline(y=30, color='g', linestyle='-.',label='30m FS')
    ax[1].axhline(y=100, color='grey', linestyle='-.',label='100m FS')
    ax[1].axhline(y=block_ht, color='b', linestyle='-.',label='Block Height')
    
    ax[0].axhline(y=block_ht, color='b', linestyle='-.',label='Block Height')
    
    
    ax[0].legend()
    ax[1].legend()
    
    
    ax[0].set_xlabel('U/Uref')
    ax[1].set_xlabel('T.I.  (%)')
    ax[0].set_ylabel('Z F.S. (m)')
    
    ax[0].set_xlim(.2, 1.2)
    ax[1].set_xlim(0, 50)
    
    plt.suptitle(tun+' - '+targets+' -'+fetch)
    
    #plt.show() #moved to after loop so that all lpots show at same time

#%% error calculations

def Calc_Error(fetch,tun):    
    #filter for fetch and WT
    df=F_t(F_f(df_og,fetch),tun)
    
    if not df.empty:
        #trim down data set to "good" window
        df=F_d(df,'Vertical Profile')
        df=df[df.Location.isin([
            #'Ref',
            #'-3r/2',
            '-r/2',
            'TTC',
            #'+r/2',        
                ])]
        
        df['Z_fs']=df['Z(mm)']*SF/1000
        
        main=df['Main(mm)'].max()
        boost=df['Boost(mm)'].max()
        block_ht=max(main,boost)*SF/1000
        
        # limit to "good" height region
        df=df[df['Z_fs'].between(block_ht*2,Zmax)]
        
        Z=df.Z_fs
        U_targ=pwr(Z,SF,Zo_targ)
        TI_targ=TI(Z,Zo_targ)
        
        #least squares method (i think)
        U_err=sum((U_targ-df['Uhf/Uref'])**2 + (U_targ-df['Uomni/Uref'])**2)**.5
        TI_err=sum((TI_targ-df['TIUhf(%)'])**2 + (TI_targ-df['TIUomni(%)'])**2)**.5
        
        return U_err, TI_err
    
    else: return 999, 999

#%% Get file & initialize

#hard-coded filepath of OG data (Cross-network)
#filepath='//gd-fs3/CentralDataCache/Projects/G153E/Results_Post_Processed/BLhunts/2020-06-GoldPlateOutput1/GoldPlateBrief.xlsx'

#local copy - set working directory to location of this script first
filepath='./GP_Data\GoldPlateBrief.xlsx'

#path of pickled (saved) dataframe
pklpath='./GP.pkl'

if os.path.isfile(pklpath):         #check if .pkl exists already   
    if easygui.ynbox("Saved data found, do you want to load the saved data?"):
        print('Loading saved dataframe...')
        df_og=pd.read_pickle(pklpath)
    else:
        if os.path.isfile(filepath) is False:   #open file dialogue if .xlsx path not found
            print('Default file not found')
            filepath=Get_Data()
        print('Reading data file...')
        df_og = pd.read_excel(filepath, header=0)
        print('Saving dataframe...')
        df_og.to_pickle(pklpath)            #save data to .pkl file
else:                                       #if no .pkl, try to load data
    if os.path.isfile(filepath) is False:   #open file dialogue if .xlsx path not found
        print('Default file not found')
        filepath=Get_Data()
    print('Reading data file...')
    df_og = pd.read_excel(filepath, header=0)
    print('Saving dataframe...')
    df_og.to_pickle(pklpath)
    

# get inputs
# add any extra gui interfaces here (tunnel, profile code, etc.)
tun,SF,Zo_targ=Inputs()
targets='{}-{}'.format(SF,Zo_targ)

#%% Begin real script
#t0 = datetime.datetime.now()    #start process timer

#set plot style 
mpl.rcdefaults()            #reset to defaults
styles=plt.style.available  #show plot styles
plt.style.use(styles[14])
cm1='rainbow'
cm2=cm1

Zmax=100

results=pd.DataFrame()
results['FetchID']=df_og['Fetch'].dropna().drop_duplicates().reset_index(drop=True)

c1='{}-{}-{}'.format(tun,targets,'U')
c2='{}-{}-{}'.format(tun,targets,'TI')

results[c1]=''
results[c2]=''

print('Calculating results...')

for i in results.index:
    fetch=results.FetchID[i]
    results[c1].loc[i],results[c2].loc[i]=Calc_Error(fetch,tun)
    
results['Score']=(results[c1]**2+results[c2]**2)**.5

best_score=results.Score.min()
best_index=results.index[results.Score==best_score][0]
best_fetch=results.FetchID[best_index]

ranked_U=results.sort_values(by=[c1]).copy().reset_index(drop=True)
top_U=ranked_U.FetchID[0:3]

print('Generating plots...')

for i in top_U:
    Profile_plot(i, tun)
    
ranked_TI=results.sort_values(by=[c2]).copy().reset_index(drop=True)
top_TI=ranked_TI.FetchID[0:3]

for i in top_TI:
    Profile_plot(i, tun)
    
ranked_score=results.sort_values(by=['Score']).copy().reset_index(drop=True)
top=ranked_score.FetchID[0:3]

for i in top:
    Profile_plot(i, tun)

#t1=datetime.datetime.now()

print('Done! Close all plot windows to exit')
plt.show()

#%% Done
# =============================================================================
# print('Done!')
# #t1=datetime.datetime.now()
# dt= t1-t0
# dt=dt.seconds
# easygui.msgbox(msg="Done!\n Process took: {} seconds".format(dt))  
# =============================================================================
