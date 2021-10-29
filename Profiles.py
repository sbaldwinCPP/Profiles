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
    fieldValues = ['GDW', 240, 0.35]  #Defaults
    inputs=easygui.multenterbox(msg, title, fieldNames, fieldValues)
    if inputs==None: sys.exit()
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


#%% functions from CPP AQ Handbook and ADC report

#   Power Law
#   ADC report - equation B.3
#   Also equation 21 in CPP Handbook
def pwr(z,zref,Zo):    #   u/Uref = (z/Zo)**n
    n=nFromZo(Zo)
    return (z/zref)**n

#   Exponent calculation from Zo
#   ADC report - equation B.4 (from Counihan, 1975)
#   Also equation 22 in CPP handbook
def nFromZo(Zo):    
    return 0.24 + 0.096*np.log10(Zo) + 0.016*(np.log10(Zo)**2)

#   turbulence equation
#   CPP handbook - equation 23 (from EPA, 1981)
#   valid between 5 and 100 m
def TI(z,Zo):
    n=nFromZo(Zo)
    return n*(np.log(30/Zo))/(np.log(z/Zo))*100

#   Log Law
#   ADC report - equation B.1
def log_ADC(z,Zo):      #   U/Uinf = 2.5*(UstarZero/Uinf)*ln(z/Zo)
    FR=F_ratio(Zo)      #   (UstarZero/Uinf)
    return 2.5*FR*np.log(z/Zo)

#   Friction Velocity Ratio
#   ADC report - equation B.2
def F_ratio(Zo):
    return np.sqrt(2.75e-3+6e-4*np.log10(Zo))

#%%
#   Log Law from excel example
def Ustar_xl(Zo,Zref,Uref):
    #=SQRT(0.00275+0.0006*LOG(O44))*(600/F38)^O41*VLOOKUP($F$37,$A$16:$I$36,9)
    n=nFromZo(Zo)
    return np.sqrt(0.00275+0.0006*np.log10(Zo))*(600/Zref)**n*Uref

def Log_xl(z,Zo,Uref,Ustar):
    #=IF(E16<>"",($O$42/0.4*(LN(D16)-LN($O$44)))/VLOOKUP($F$37,$A$16:$I$36,9),NA())
    return (Ustar/0.4*(np.log(z)-np.log(Zo)))/Uref
    
def TI_xl(z,Zo):
    #=IF(E16<>"",IF(D16<100,O$41*LN(30/(O$44))/(LN((E16/100*$E$3)/O$44)),($O$41*LN(30/$O$44)/LN(100/$O$44)+(D16-100)/500*(0.01-$O$41*LN(30/$O$44)/LN(100/$O$44)))),NA())
    n=n=nFromZo(Zo)
    try:
        t1=n*np.log(30/Zo)/np.log(z[z<100]/Zo)
        t2=n*(np.log(30/Zo))/np.log(100/Zo)+(z[z>=100]-100)/500*(0.01-n*np.log(30/Zo)/np.log(100/Zo))
        t=np.append(t1,t2)
    except:
        if z < 100:
            t=n*np.log(30/Zo)/np.log(z/Zo)
        else: 
            t=n*(np.log(30/Zo))/np.log(100/Zo)+(z-100)/500*(0.01-n*np.log(30/Zo)/np.log(100/Zo))
    return t*100
    
#%% log law from wikipedia

# =============================================================================
# #   Holmes JD. Wind Loading of Structures. 3rd ed. Boca Raton, Florida: CRC Press; 2015.
# def log_w(z,zref,Zo):
#     #d=  # zero-plane displacement:
#          # height in meters above the ground at which zero wind speed is 
#          # achieved as a result of flow obstacles such as trees or buildings
#     d=0
#     return np.log((z-d)/Zo)/np.log((zref-d)/Zo)
# 
# def Ustar_over_Uref(Zref,Zo): 
#     #adapted from U(z)=(u*/k)*ln((z-d)/Zo)
#     k=0.41
#     d=0
#     return k/np.log((Zref-d)/Zo)
# 
# def Ustar(u1,u2,z1,z2):
#     k=0.41
#     d=0
#     return k*(u2-u1)/np.log((z2-d)/(z1-d))
# =============================================================================


#%%   Log functions from RP's matlab code

# =============================================================================
# 
# #   not in use, does not provide a good fit to any data
# def old_log(z,Zo):    #   u/Uinf = 2.5*(Ustar0/Uinf)*ln(z/Zo)
#     Fr=F_ratio(Zo)
#     return 2.5*Fr*np.log(z/Zo)
# 
# def F_ratio(Zo):
#     return np.sqrt(2.75e-3+6e-4*np.log10(Zo))
# 
# # =============================================================================
# #   ###     Also found this in the .txt copy of RP code
# #   #LogLaw=@(zo,z)((sqrt(2.75*10^-3+6*10^-4*log10(zo))*...
# #           #(600/SF)^(0.24+0.096*log10(zo)+0.016*log10(zo)^2)*Uref)/0.41*log(z/zo));
# # =============================================================================
# 
# def sqrt(x):
#     return x**.5
# 
# def log_law(z,zo,Uref):
#     return sqrt(2.75e-3+6e-4*np.log10(zo))*(600/SF)**(0.24+0.096*np.log10(zo)+0.016*np.log10(zo)**2)*Uref/0.41*np.log(z/zo)
# =============================================================================


#%% plotting function
#   mostly for QA
def Profile_plot(fetch, tun, label=''):
    
    data=F_f(df_og,fetch)
    
    data['Z_fs']=data['Z(mm)']*SF/1000
    
    main=data['Main(mm)'].mean()
    boost=data['Boost(mm)'].mean()
    t1=data['T1(mm)'].mean()
    t2=data['T2(mm)'].mean()
    s1=data['S1(mm)'].mean()
    s2=data['S2(mm)'].mean()
    
    block_ht=max(main,boost)*SF/1000
    trip_ht=max(t1,t2)*SF/1000
    spire_ht=max(s1,s2)*SF/1000
    
    data=F_t(data,tun)
    
    data=data[data.Location.isin([
        'Ref',
        '-3r/2',
        '-r/2',
        'TTC',
        '+r/2',        
            ])]
    
    #sort by abs distance from TT center, show most centered points on top of plot
    data['abs_x']=abs(data['X(mm)'])
    data=data.sort_values('abs_x', ascending=False)
    
    horiz=F_d(data,'Lateral Profile')
    data=F_d(data,'Vertical Profile')

    z = np.linspace(1,SF*1.2,num=1000)

    #### actual plotting, all above is just data manipulation/prep
    
    figure,ax = plt.subplots(1,2, sharey=True, figsize=(11,7))
    
    # target U pwr - ADC
    n=nFromZo(Zo_targ)
    ax[0].plot(pwr(z,SF,Zo_targ), z , c='grey', label='Target Power, n={})'.format(round(n,2)))
    
    # target U log - EXCEL Example
    Uref = data['Uref(m/s)'].mean()
    Us = Ustar_xl(Zo_targ,SF,Uref)
    ax[0].plot(Log_xl(z,Zo_targ,Uref,Us), z, c='grey', linestyle='-.', label='Target Log, U*={}'.format(round(Us,2)))
    
    # target TI - EXCEL Example
    ax[1].plot(TI_xl(z,Zo_targ), z, c='grey', linestyle='--', label='Target TI')
    
    # horizontal profile setup
    horiz['Z_fsr']=round(horiz['Z_fs'],-1)
    hts=horiz['Z_fsr'].drop_duplicates()
    horiz['hf_avg']=''
    horiz['omni_avg']=''
    ps=SF*.75 # scalar to magnify percent variation
    
###   single horizontal profile
    h = min(hts, key=lambda x:abs(x-100)) #finds value in hts closest to 100
    mean_hf=horiz[horiz['Z_fsr']==h]['Uhf/Uref'].mean()
    mean_omni=horiz[horiz['Z_fsr']==h]['Uomni/Uref'].mean()
    mean_combo=(mean_omni+mean_hf)/2
    horiz.loc[horiz['Z_fsr']==h, 'hf_avg']=mean_hf
    horiz.loc[horiz['Z_fsr']==h, 'omni_avg']=mean_omni
    horiz.loc[horiz['Z_fsr']==h, 'combo_avg']=mean_combo
    
    ax[0].axhline(y=(h+ps*.05), color='grey', linestyle='--', xmax=0.4)
    ax[0].axhline(y=(h-ps*.05), color='grey', linestyle='--', xmax=0.4)
###

    #enable to plot all horizontal profiles, messy for some data sets
# =============================================================================
#     for h in hts:
#         mean_hf=horiz[horiz['Z_fsr']==h]['Uhf/Uref'].mean()
#         mean_omni=horiz[horiz['Z_fsr']==h]['Uomni/Uref'].mean()
#         mean_combo=(mean_omni+mean_hf)/2
#         horiz.loc[horiz['Z_fsr']==h, 'hf_avg']=mean_hf
#         horiz.loc[horiz['Z_fsr']==h, 'omni_avg']=mean_omni
#         horiz.loc[horiz['Z_fsr']==h, 'combo_avg']=mean_combo
#         
#         ax[0].axhline(y=(h+ps*.05), color='grey', linestyle='--', xmax=0.4)
#         ax[0].axhline(y=(h-ps*.05), color='grey', linestyle='--', xmax=0.4)
# =============================================================================
    
    # just for single label, plots way off screen
    ax[0].axhline(y=3*SF, color='grey', linestyle='--', label='+- 5% bounds')
    
    # horizontal lines for FS heights and tunnel setup
    ax[1].axhline(y=30, color='g', linestyle='-.',label='30m FS')
    ax[1].axhline(y=100, color='grey', linestyle='-.',label='100m FS')
    ax[1].axhline(y=block_ht, color='b', linestyle='--',label='Block, Trip, Spire Heights')
    ax[1].axhline(y=trip_ht, color='b', linestyle='--')
    ax[1].axhline(y=spire_ht, color='b', linestyle='--')
    
    # velocity data
    data.plot(kind='scatter', x='Uhf/Uref', y='Z_fs',ax=ax[0], marker='x', c='X(mm)', cmap=cm, colorbar=False, label='HF')
    data.plot(kind='scatter', x='Uomni/Uref', y='Z_fs',ax=ax[0], marker='+', c='X(mm)', cmap=cm, colorbar=False, label='12HP')
    
    # turbulence data
    data.plot(kind='scatter', x='TIUhf(%)', y='Z_fs',ax=ax[1], marker='x', c='X(mm)', cmap=cm, colorbar=False)
    data.plot(kind='scatter', x='TIUomni(%)', y='Z_fs',ax=ax[1], marker='+', c='X(mm)', cmap=cm, colorbar=True)

    # horizontal profile plots
    horiz['HFxvar']=0.3*horiz['Yhf(mm)']/horiz['Yhf(mm)'].max()+0.3    # normalize y position to go from 0 to .6
    horiz['Omnixvar']=0.3*horiz['Yomni(mm)']/horiz['Yomni(mm)'].max()+0.3
    
    horiz['HFyvar']=horiz['Z_fsr']+ps*(horiz['Uhf/Uref']-horiz['combo_avg'])/horiz['combo_avg'] 
    horiz['Omniyvar']=horiz['Z_fsr']+ps*(horiz['Uomni/Uref']-horiz['combo_avg'])/horiz['combo_avg']  
    
    horiz.plot(kind='scatter', x='HFxvar', y='HFyvar',ax=ax[0], marker='x', c='X(mm)', cmap=cm, colorbar=False)
    horiz.plot(kind='scatter', x='Omnixvar', y='Omniyvar',ax=ax[0], marker='+', c='X(mm)', cmap=cm, colorbar=False)
    
    ax[0].legend()
    ax[1].legend()
    
    ax[0].set_xlabel('U/Uref')
    ax[1].set_xlabel('T.I.  (%)')
    ax[0].set_ylabel('Z F.S. (m)')
    
    ax[0].set_xlim(.2, 1.2)
    ax[1].set_xlim(0, 50)
    
    ax[0].set_ylim(0, 1.5*SF)
    
    ID=results[results.FetchID==fetch].index[0]
    
    plt.suptitle('{}:{} \n ID#{}:{} \n {}'.format(tun, targets, ID, fetch, label))
    

    Plots_Folder=os.path.join(inidir,'Plots')
    if not os.path.exists(Plots_Folder): os.makedirs(Plots_Folder)
    Plot_Path=os.path.join(Plots_Folder,'{}_{}_{}.png'.format(tun,targets,label))
    plt.savefig(Plot_Path)
    plt.close()

#%% error calculations

def Calc_Error(fetch,tun):    
    #filter for fetch and WT
    df=F_t(F_f(df_og,fetch),tun)
    
    #trim down data set to "good" window
    df=F_d(df,'Vertical Profile')
    df=df[df.Location.isin([
                            #'Ref',
                            #'-3r/2',
                            '-r/2',
                            'TTC',
                            #'+r/2',        
                                ])]
    
    if not df.empty:
        df['Z_fs']=df['Z(mm)']*SF/1000
        
        # limit to "good" height region for each law
        # TI uses the pwr range
        # need to check if its worth the effort/complexity
        df_pwr=df[df['Z_fs'].between(20,SF)]
        df_log=df[df['Z_fs'].between(5,100)]
        
        Z_pwr=df_pwr.Z_fs
        Z_log=df_log.Z_fs
        n_pwr=len(Z_pwr)
        n_log=len(Z_log)

        U_pwr=pwr(Z_pwr,SF,Zo_targ)
        
        Uref = df_log['Uref(m/s)'].mean()
        Us = Ustar_xl(Zo_targ,SF,Uref)        
        U_log=Log_xl(Z_log,Zo_targ,Uref,Us)
        
        TI_targ=TI_xl(Z_pwr,Zo_targ)
        
        uhf_pwr=df_pwr['Uhf/Uref']
        uomni_pwr=df_pwr['Uomni/Uref']
        Uavg_pwr=(uhf_pwr+uomni_pwr)/2
        
        uhf_log=df_log['Uhf/Uref']
        uomni_log=df_log['Uomni/Uref']
        Uavg_log=(uhf_log+uomni_log)/2
        
        TIhf=df_pwr['TIUhf(%)']
        TIomni=df_pwr['TIUomni(%)']
        TIavg=(TIhf+TIomni)/2
        
        #calculate average error     
        U_pwr_err=abs(U_pwr-Uavg_pwr)
        U_log_err=abs(U_log-Uavg_log)
        TI_err=abs(TI_targ-TIavg)
        
        U_pwr_err_avg=U_pwr_err.mean()
        U_log_err_avg=U_log_err.mean()
        TI_err_avg=TI_err.mean()
        
        return U_pwr_err_avg, U_log_err_avg, TI_err_avg
    
    else: return 999, 999, 999

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

# get inputs, in future may turn into an array of values to loop thru
tun,SF,Zo_targ=Inputs()
targets='{}-{}'.format(SF,Zo_targ)

#%% Begin real script
t0 = datetime.datetime.now()    #start process timer

results=pd.DataFrame()
results['FetchID']=df_og['Fetch'].dropna().drop_duplicates().reset_index(drop=True)

c1='{}-{}-{}'.format(tun,targets,'U_pwr')
c2='{}-{}-{}'.format(tun,targets,'U_log')
c3='{}-{}-{}'.format(tun,targets,'TI')

results[c1]=''
results[c2]=''
results[c3]=''

print('Calculating results...')

for i in results.index:
    fetch=results.FetchID[i]
    results[c1].loc[i],results[c2].loc[i],results[c3].loc[i]=Calc_Error(fetch,tun)
    
results['Score']=(results[c1]*results[c2]*results[c3]**2)

inidir=os.getcwd()
Results_Folder=os.path.join(inidir,'Results')
Results_Path=os.path.join(Results_Folder,'{}_{}_results.csv'.format(tun,targets))
if not os.path.exists(Results_Folder): os.makedirs(Results_Folder)

results.to_csv(Results_Path)

#%% plot results
#fetch=results.FetchID[10]  #for manual fetch selection

#set plot style 
mpl.rcdefaults()            #reset to defaults
styles=plt.style.available  #save all plot styles to  list
plt.style.use(styles[14])   #set style
cm='twilight'               #set colormap
#cm='twilight_shifted'      #if plot style has dark background

print('Generating plots...')

n=3 #number of plots per set

ranked_pwr=results.sort_values(by=[c1]).copy().reset_index(drop=True)
top_pwr=ranked_pwr.FetchID[0:n]
for i in range(n):
    f=top_pwr[i]
    label='Best Power Law fit #{}'.format(i+1)
    Profile_plot(f, tun, label)
    
ranked_log=results.sort_values(by=[c2]).copy().reset_index(drop=True)
top_log=ranked_log.FetchID[0:n]
for i in range(n):
    f=top_log[i]
    label='Best Log Law fit #{}'.format(i+1)
    Profile_plot(f, tun, label)
    
ranked_TI=results.sort_values(by=[c3]).copy().reset_index(drop=True)
top_TI=ranked_TI.FetchID[0:n]
for i in range(n):
    f=top_TI[i]
    label='Best TI fit #{}'.format(i+1)
    Profile_plot(f, tun, label)
    
ranked_score=results.sort_values(by=['Score']).copy().reset_index(drop=True)
top_score=ranked_score.FetchID[0:n]
for i in range(n):
    f=top_score[i]
    label='Best overall #{}'.format(i+1)
    Profile_plot(f, tun, label)

t1=datetime.datetime.now()
print('Close all plot windows to continue...')
#plt.show()          #shows all plots at once, holds until closed

#%% Done
print('Done!')

dt= t1-t0
dt=dt.seconds
easygui.msgbox(msg="Done!\n Process took: {} seconds".format(dt))  
