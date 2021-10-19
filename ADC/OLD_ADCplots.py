#%% ADC Plotting
# SAB 4/18/20
# python 3.7 in spyder

# =============================================================================
# D:\Conda\python.exe 
# D:\Projects\ADC & Benchmark\ADC/NEW_ADCplots.py
# =============================================================================


#%% Get_file
print('Importing modules...')
import datetime
import tkinter
from tkinter import filedialog
from tkinter import messagebox as mb
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import sys

print('Get File...')
#get file to use
root = tkinter.Tk()
root.withdraw()                 #use to hide tkinter window
root.overrideredirect(True)
root.geometry('0x0+0+0')
root.deiconify()
root.lift()
root.focus_force()              #force key focus to window

inidir = os.path.dirname(os.path.abspath(__file__)) #starts file selection at location of this script 
  
filepath = filedialog.askopenfilename(parent=root, initialdir=inidir, 
                                      title='Please select a File',
                                      filetypes=[('.csv','.csv xlsx')])

# =============================================================================
# filepath = filedialog.askopenfilename(parent=root, initialdir=inidir, 
#                                       title='Please select a File',
#                                       filetypes=[('.xls','.xlsx')])
# =============================================================================

if filepath == '':
    print('No file')
    mb.showerror('Bye Felica','Bye Felica')
    sys.exit()

#start process timer
print('On your marks, get set, go!')
t0 = datetime.datetime.now()

#%% Import_data_&_Process_headers

df_og = pd.read_excel(filepath, header=None)
#df_og = pd.read_csv(filepath, header=None)

#concatenate top 4 rows into  1 header
#probably a better way to to do this with a loop, brute force instead
def conc(df):
    a= df.iloc[0,:].astype(str)
    b= df.iloc[1,:].astype(str)
    c= df.iloc[2,:].astype(str)
    d= df.iloc[3,:].astype(str)
    df.columns= [i + j + k + l for i, j, k, l in zip(a, b, c, d)]
    df = df.iloc[4:].reset_index(drop=True)     # reset indexes
    return df

df=conc(df_og)

#remove special characters, make all lowercase
def rep(s):
    s1=''.join(e for e in s if e.isalnum()).lower()
    return s1

df.columns =[rep(column) for column in df.columns]

#remove nan
df.columns =[column.replace("nan", "") for column in df.columns]

#build dictionary from dataframe 
#so far unused
#dct = df.to_dict()           

#search_utility

def srch(kw):
    s=[col for col in df.columns if kw in col]
    return s

#%% Build data sets for plotting
    
#use search utl to define & store keywords for useful columns
#not the most robust method, requires some testing, should be improved
    
data=pd.DataFrame()

data['T1']= df[srch('trip1a')[0]] # a for actual
data['T2']= df[srch('trip2a')[0]] # a for actual

data['S1']= df[srch('spires1a')[0]] # a for actual
data['S2']= df[srch('spires2a')[0]] # a for actual

data['M1']= df[srch('main1a')[0]] # a for actual
data['B1']= df[srch('boost1a')[0]] # a for actual

data['uref']= df[srch('urefmean')[0]]
data['uref_var']=df[srch('urefvar')[0]]
data['uhf']=df[srch('uhfmean')[0]]
data['uhf_var']=df[srch('uhfvar')[0]]


data['X']=df[srch('uhfxpos')[0]].astype(float)

#data['Y']=df[srch('uhfypos')[0]].astype(float)
data['Y']=df[srch('12hp1ypos')[0]].astype(float)

data['Z']=df[srch('uhfzpos')[0]].astype(float)

# =============================================================================
# omniX=srch('12xpos')[0]
# omniY=srch('12ypos')[0]
# omniZ=srch('12zpos')[0]
# use hf for all xyz locations
# =============================================================================

#select data set for omni 1:raw 2:afgoodv 3:afgoodvp
ds='afgoodv'   
data['uvw']=df[srch(ds+'overlineuvw')]
data['U']=df[srch(ds+'overlineu')[0]]
data['V']=df[srch(ds+'overlinev')[0]]
data['W']=df[srch(ds+'overlinew')[0]]

data['U_var']=df[srch(ds+'sigma2u')[0]]
data['V_var']=df[srch(ds+'sigma2v')[0]]
data['W_var']=df[srch(ds+'sigma2w')[0]]
data['uvw_var']=df[srch(ds+'sigma2uvw')[0]]

#%% make computations
data['uhf_n']=data['uhf']/data['uref']
data['U_n']=data['U']/data['uref']
data['z_n']=data['Z']/1000

#remove calcheck and pressure zero runs (not effective??) filters moved to plots
#df=df.query('usercomments != "refcalchkhs"', inplace = False)
#df=df.query(uref+' > 3',  inplace = False)

#Good HF
a=data.query('uhf_n < 2', inplace= False)
max(a)


# TI%
data['uhf_ti']= data['uhf_var']**(1/2)/data['uhf']
data['U_ti']= data['U_var']**(1/2)/data['U']
data['V_ti']= data['V_var']**(1/2)/data['U']
data['W_ti']= data['W_var']**(1/2)/data['U']

#filter by trip 1 to sort roughnii
R1=data.query('97< T1 < 103', inplace = False) #300-Zo0.1
R2=data.query('297< T1 < 303', inplace = False) #180-Zo1.0
R3=data.query('397< T1 < 403', inplace = False) #180-Zo0.5
#%%Plots

mpl.rcdefaults()           #reset to defaults
styles=plt.style.available  #show plot styles
plt.style.use(styles[2])

def pwr(x,a,b):
    return x**a+b   

def ipwr(x,a,b):
    return x**-a+b

def exp(u,z0,a):
    u=u.astype(float)
    return z0*np.exp(u*a)

def fit(df,func,independantkw,dependantkw):
    f=pd.DataFrame()
    f['xval']=df[independantkw]
    f['yval']=df[dependantkw]
    f=f.sort_values(by='xval')
    coef, pcov = opt.curve_fit(func, f['xval'], f['yval']);
    err = np.sqrt(np.diag(pcov))
    return f, coef, err

def vert(dat,x,y,uref):
    #filter for desired data, adjust tolerance for each search parameter here
    v=dat.query(str(x-5)+'<X<'+str(x+5),inplace=False)
    v=v.query(str(y-20)+'<Y<'+str(y+20),inplace=False)           
    v=v.query(str(uref-2)+'<uref<'+str(uref+2),inplace=False) 
    
# =============================================================================
#     f0,coef0,err0=fit(v,exp,'U_n','z_n')
#     ustar0=0.4/coef0[1]
#     
#     f1,coef1,err1=fit(v,pwr,'U_n','z_n')
#     f2,coef2,err2=fit(v,ipwr,'U_ti','z_n')
# =============================================================================
    
    figure,ax = plt.subplots(1,2, sharey=True)
    
#velocity
    #data
    v.plot(kind='scatter', x='uhf_n', y='z_n',ax=ax[0], marker='D', c='b',label='U HF')
    v.plot(kind='scatter', x='U_n', y='z_n',ax=ax[0], marker='D', c='r',label='U omni')
    
# =============================================================================
#     #fits
#     ax[0].plot(f0['xval'], exp(f0['xval'], *coef0), c='g', 
#     label="Z0={:.2f}, U*={:.2f} err={:.2e}".format(coef0[0],ustar0,err0[0]))
#     
#     ax[0].plot(f1['xval'], pwr(f1['xval'], *coef1), c='r',
#     label="x^{:.1f}+{:.1f}, err={:.2f}".format(coef1[0],coef1[1],err1[0]));
#     
# =============================================================================
    
#turbulence
    #data
    v.plot(kind='scatter', x='uhf_ti', y='z_n', ax=ax[1], marker='o', c='b',label='HF TI')
    v.plot(kind='scatter', x='U_ti', y='z_n', ax=ax[1], marker='o', c='r',label='omni TI')

# =============================================================================
#     #fits
#     ax[1].plot(f2['xval'], ipwr(f2['xval'], *coef2), c='r', 
#     label="x^-{:.1f}+{:.1f}, err={:.2f}".format(coef2[0],coef2[1],err2[0]));
#     #ax[1].plot(f2['xval'], ipwr(f2['xval'], *coef2), c='r', label="fit2")
#     
# =============================================================================
    ax[0].set_xlabel('U/Uref')
    ax[1].set_xlabel('T.I.')
    ax[0].set_ylabel('Z/Zref')
    #ax[0].set_ylim([0,1.25])
    #ax[0].set_xlim([0,1.25])
    ax[0].legend()
    ax[1].legend()
    plt.title(label='X={}, Y={}, uref={}'.format(x,y,uref), loc='left')
    plt.show()                     #turn on to plot new fig each call
    
# =============================================================================
#     coef=pd.DataFrame()
#     coef[0,1,2]=coef0,coef1,coef2
#     
# =============================================================================
    #plt.savefig('demo.png')
#    return coef
    
#v1=vert(R1,790,500,12)
#v2=vert(R2,790,500,12)
#v3=vert(R3,0,500,12)

#array of inputs to feed plot functions
R=[R1,R2,R3]
xloc=[0,790,1590,2390]
yloc=[0,500]
refspeed=[12]
for dat in R:
    for x in xloc:
        for y in yloc:
            for uref in refspeed:
                vert(dat,x,y,uref)
            

#plt.show()                         #flip-flop w/ above to plot all, 1 fig
#%% Done
t1 = datetime.datetime.now()
dt= t1-t0
dt=dt.seconds

mb.showinfo('Task Complete',"Total elapsed time: {} seconds".format(dt), parent=root) 
root.destroy()
