#%% GP Plotting
# SAB 6/22/20
# python 3.7 in spyder

#%% import modules
# =============================================================================
# import tkinter
# from tkinter import filedialog, messagebox as mb
# =============================================================================

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
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    
    print('All modules imported')

except ModuleNotFoundError as err:
    input(str(err) +' \n Press enter to exit')

#%% Functions

# files
def Get_Data():
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
    msg = "Enter Settings"
    title = "Enter Settings"
    fieldNames = ["Scale:","Target Zo:"]
    fieldValues = [240,0.35]  #Defaults
    inputs=easygui.multenterbox(msg, title, fieldNames, fieldValues)
    inputs= int(inputs[0]),float(inputs[1])
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
    return n*(np.log(30/Zo))/(np.log(z/Zo))

#%% Log functions from RP's matlab code
def log(z,Zo):    #   u/Uinf = 2.5*(Ustar0/Uinf)*ln(z/Zo)
    Fr=F_ratio(Zo)
    return 2.5*Fr*np.log(z/Zo)

def F_ratio(Zo):
    return np.sqrt(2.75e-3+6e-4*np.log10(Zo))

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
            print('Default file not found, use GUI to select data file...')
            filepath=Get_Data()
        print('Reading data file...')
        df_og = pd.read_excel(filepath, header=0)
        print('Saving dataframe...')
        df_og.to_pickle(pklpath)                #save data to .pkl file
else:                               #if no .pkl, try to load data
    if os.path.isfile(filepath) is False:   #open file dialogue if .xlsx path not found
        print('Default file not found, use GUI to select data file...')
        filepath=Get_Data()
    print('Reading data file...')
    df_og = pd.read_excel(filepath, header=0)
    print('Saving dataframe...')
    df_og.to_pickle(pklpath)
    

# get inputs
# add any extra gui interfaces here (tunnel, profile code, etc.)
SF,Zo_targ=Inputs()

#%% Begin real script
t0 = datetime.datetime.now()    #start process timer

# build dataframe of approaches and coef's(added later)
results=pd.DataFrame()
results['FetchID']=df_og['Fetch'].dropna().drop_duplicates().reset_index(drop=True)
results['Zo']=''
#data['{}-Zo'.format(SF)]=''   #include SF in column name


#%% workspace for plot function

#set plot style 
mpl.rcdefaults()            #reset to defaults
styles=plt.style.available  #show plot styles
plt.style.use(styles[17])



fetch=results['FetchID'][11]  # this will become input arg for function

data=F_f(df_og,fetch)

main=data['Main(mm)'].mean()
boost=data['Boost(mm)'].mean()

block_ht=max(main,boost)*SF/1000

data=data[data.Location=='TTC']
data=F_d(data,'Vertical Profile')

data['Z_fs']=data['Z(mm)']*SF/1000

Z = np.linspace(block_ht,SF,num=100)

#actual plotting, all above is just data manipulation/prep

figure,ax = plt.subplots(1,2, sharey=True)

# velocity
F_t(data,'GDE').plot(kind='scatter', x='Uhf/Uref', y='Z_fs',ax=ax[0], marker='x', c='r',label='GDE HF')
F_t(data,'GDE').plot(kind='scatter', x='Uomni/Uref', y='Z_fs',ax=ax[0], marker='x', c='m',label='GDE 12HP')

F_t(data,'GDW').plot(kind='scatter', x='Uhf/Uref', y='Z_fs',ax=ax[0], marker='+', c='b',label='GDW HF')
F_t(data,'GDW').plot(kind='scatter', x='Uomni/Uref', y='Z_fs',ax=ax[0], marker='+', c='g',label='GDW 12HP')

# turbulence
F_t(data,'GDE').plot(kind='scatter', x='TIUhf(%)', y='Z_fs',ax=ax[1], marker='x', c='r',label='GDE HF')
F_t(data,'GDE').plot(kind='scatter', x='TIUomni(%)', y='Z_fs',ax=ax[1], marker='x', c='m',label='GDE 12HP')

F_t(data,'GDW').plot(kind='scatter', x='TIUhf(%)', y='Z_fs',ax=ax[1], marker='+', c='b',label='GDW HF')
F_t(data,'GDW').plot(kind='scatter', x='TIUomni(%)', y='Z_fs',ax=ax[1], marker='+', c='g',label='GDW 12HP')

# targets
ax[0].plot(pwr(Z,SF,Zo_targ), Z, c='c')

ax[1].axhline(y=30, color='g', linestyle='-.',label='30m FS')
ax[1].axhline(y=100, color='k', linestyle='-.',label='100m FS')
ax[1].axhline(y=block_ht, color='b', linestyle='-.',label='Block Height')

ax[0].axhline(y=block_ht, color='b', linestyle='-.',label='Block Height')


#ax[0].legend()
ax[1].legend()


ax[0].set_xlabel('U/Uref')
ax[1].set_xlabel('T.I.  (%)')
ax[0].set_ylabel('Z F.S. (m)')

plt.suptitle(fetch)

plt.show()


#%% error calcs

def Calc_error(fetch,tun):
    #fetch=results['FetchID'][10]  # this will become input arg for function
    #tun='GDW'
    
    #trim down data set to "good" window
    df=F_t(F_f(df_og,fetch),tun)
    df=F_d(df,'Vertical Profile')
    df=df[df.Location=='TTC']
    
    df['Z_fs']=df['Z(mm)']*SF/1000
    
    main=df['Main(mm)'].max()
    boost=df['Boost(mm)'].max()
    block_ht=max(main,boost)*SF/1000
    
    # limit to "good" height region
    df=df[df['Z_fs'].between(block_ht*1.5,100)]
    
    
    Z_pwr=df.Z_fs
    U_pwr=pwr(Z_pwr,SF,Zo_targ)
    
    err=sum((U_pwr-df['Uhf/Uref'])**2)**.5
    
    return err


#df.plot(kind='scatter', x='Uhf/Uref', y='Z_fs',ax=ax[0], marker='x', c='r',label='GDE HF')


#%%
for i in results.index:
    fetch=results.FetchID[i]
    tun='GDE'
    results.Zo.loc[i]=Calc_error(fetch,tun)



#%% old code for reference
#%%
#Fit functions
    
# =============================================================================
# def fit(df,func,independantkw,dependantkw):
#     f=pd.DataFrame()
#     f['xval']=df[independantkw]*SF      #Z(mm) MS 
#     f['yval']=df[dependantkw]           #U/Uref
#     #f=f.sort_values(by='xval')
#     #f=f.sort_values(by='yval')
#     coef, pcov = opt.curve_fit(func, f['xval'], f['yval']);
#     err = np.sqrt(np.diag(pcov))
#     return f, coef, err
# 
# def pwr(z,zref,Zo):    #   u/Uref = (z/Zo)**n
#     n=nFromZo(Zo)
#     return (z/zref)**n
# 
# def log(z,Zo):    #   u/Uinf = 2.5*(Ustar0/Uinf)*ln(z/Zo)
#     #u=u.astype(float) # do this outside of function
#     Fratio=F_ratio(Zo)
#     return 2.5*Fratio*np.log(z/Zo)
# 
# def free_pwr(z,Zo,n):    #   u/Uinf = (z/Zo)**n
#     return (z/Zo)**n
# 
# def free_log(z,Fratio,Zo):    #   u/Uinf = 2.5*(Ustar0/Uinf)*ln(z/Zo)
#     #u=u.astype(float) # do this outside of function
#     return 2.5*Fratio*np.log(z/Zo)
# 
# def F_ratio(Zo):
#     return np.sqrt(2.75e-3+6e-4*np.log10(Zo))
#     
# def nFromZo(Zo):
#     return 0.096*np.log10(Zo)+0.016*(np.log10(Zo)**2)+.24
# =============================================================================

#%% Iterate
    
# =============================================================================
# for i in data['Fetch ID']:      #iterate through list of fetch IDs (placeholder)
#     print (i)
# =============================================================================


#%% single plot
# turn this into a function after it works    
    
# =============================================================================
# mpl.rcdefaults()           #reset to defaults
# styles=plt.style.available  #show plot styles
# plt.style.use(styles[5])
# =============================================================================
# =============================================================================
# 
# fetch=data['Fetch ID'][11]  # this will become input arg for function
# 
# 
# #read strings
# # =============================================================================
# # def left(s, amount):
# #     return s[:amount]
# # 
# # def right(s, amount):
# #     return s[-amount:]
# # =============================================================================
# 
# def mid(s, offset, amount):
#     return s[offset:offset+amount]
# 
# #FS/MS calculations
# 
# #block height    
# Mindex=fetch.find('M-')
# Bindex=fetch.find('B-')
# TTindex=fetch.find('TT')
# 
# main=int(mid(fetch,Mindex+2,Bindex-Mindex-2))
# boost=int(mid(fetch,Bindex+2,TTindex-Bindex-2))
# 
# ##
# f=fetchfilter(df_og,fetch)
# 
# f1=pd.DataFrame()
# f1['Z_n']=f['Z(mm)']/1000
# f=f.join(f1)
# 
# v=directionfilter(f,'Vertical Profile')
# #v=f
# vE=tunnelfilter(v,'GDE')
# vW=tunnelfilter(v,'GDW')
# 
# h=directionfilter(f,'Lateral Profile')
# hE=tunnelfilter(h,'GDE')
# hW=tunnelfilter(h,'GDW')
# =============================================================================

#fits

# =============================================================================
# Uref=np.mean(vE['Uref(m/s)'])
# SF=180
# =============================================================================

# =============================================================================
# f0,coef0,err0=fit(vE,free_pwr,'Z(mm)','Uhf/Uref')
# f01,coef01,err01=fit(vE,free_log,'Z(mm)','Uhf/Uref')
# 
# f02,coef02,err02=fit(vE,pwr,'Z(mm)','Uhf/Uref')
# f03,coef03,err03=fit(vE,log,'Z(mm)','Uhf/Uref')
# =============================================================================


#actual plotting, all above is just data manipulation
# =============================================================================
# 
# figure,ax = plt.subplots(1,2, sharey=True)
# 
# 
# # velocity
# vE.plot(kind='scatter', x='Uhf/Uref', y='Z_n',ax=ax[0], marker='.', c='b',label='GDE HF')
# vW.plot(kind='scatter', x='Uhf/Uref', y='Z_n',ax=ax[0], marker='.', c='r',label='GDW HF')
# 
# 
# # TI
# vE.plot(kind='scatter', x='TIUhf(%)', y='Z_n', ax=ax[1], marker='.', c='b',label='GDE HF')
# vW.plot(kind='scatter', x='TIUhf(%)', y='Z_n', ax=ax[1], marker='.', c='r',label='GDW HF')
# 
# 
# #fits
# 
# 
# z_pwr=np.linspace(.01,1.2,num=1000)
# 
# ax[0].plot(pwr(z_pwr,1,Zo_targ), z_pwr,  c='m', label='target pwr')
# 
# #refspeed=log(SF,Zo_targ)
# 
# #ax[0].plot(log(z_in,Zo_targ), z_in/SF,  c='c', label='target log')
# 
# # =============================================================================
# # ax[0].plot(log(z_in,Zo_targ*2), z_in/SF,  c='b', label='target log*2')
# # ax[0].plot(log(z_in,Zo_targ*.5), z_in/SF,  c='m', label='target log/2')
# # =============================================================================
# 
# # =============================================================================
# # ax[0].plot(free_pwr(z_in,*coef0), z_in,  c='g', label='free pwr')
# # ax[0].plot(free_log(z_in,*coef01), z_in, c='y', label='free log')
# # 
# # ax[0].plot(pwr(z_in,*coef02), z_in,  c='g', label='pwr')
# # ax[0].plot(log(z_in,*coef03), z_in, c='y', label='log')
# # =============================================================================
# 
# 
# ax[0].legend()
# ax[1].legend()
# 
# 
# ax[0].set_xlabel('U/Uref')
# ax[1].set_xlabel('T.I.')
# ax[0].set_ylabel('Z/Zref')
# 
# plt.suptitle(fetch)
# 
# plt.show()
# =============================================================================

#%% Done
print('Done!')
t1=datetime.datetime.now()
dt= t1-t0
dt=dt.seconds
easygui.msgbox(msg="Done!\n Process took: {} seconds".format(dt))  
