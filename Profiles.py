# -*- coding: utf-8 -*-
"""
Created on Tue Oct 5 14:16:23 2021

@author: sbaldwin
New branch to plot east and west data together


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
    from matplotlib import cm, colors
    import matplotlib as mpl
    from scipy  import stats
    
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
    fieldNames = ["Tunnel:", "Scale:", "Target Zo:", "Include TT roughness?:"]
    fieldValues = ['ALL', 240, 0.3, False]  #Defaults
    inputs=easygui.multenterbox(msg, title, fieldNames, fieldValues)
    if inputs==None: sys.exit()
    #inputs= inputs[0], int(inputs[1]), float(inputs[2]), inputs[3]
    if inputs[3].lower() in ['false','f','no']:
        inputs= inputs[0], int(inputs[1]), float(inputs[2]), False
    elif inputs[3].lower() in ['true','t','yes']:
        inputs= inputs[0], int(inputs[1]), float(inputs[2]), True
    return inputs

# og dataframe filters - optimize later to be 1 func w/ mult. cases
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
def Pwr_Law(z,zref,Zo):    #   u/Uref = (z/Zo)**n
    n=nFromZo(Zo)
    return (z/zref)**n

#   Exponent calculation from Zo
#   ADC report - equation B.4 (from Counihan, 1975)
#   Also equation 22 in CPP handbook
def nFromZo(Zo):    
    return 0.24 + 0.096*np.log10(Zo) + 0.016*(np.log10(Zo)**2)

#%% EXCEL-derived formulas
#   U* - theory
def Ustar_xl(Zo,Zref,Uref):
    #=SQRT(0.00275+0.0006*LOG(O44))*(600/F38)^O41*VLOOKUP($F$37,$A$16:$I$36,9)
    n=nFromZo(Zo)
    return np.sqrt(0.00275+0.0006*np.log10(Zo))*(600/Zref)**n*Uref

#   Log Law 
def Log_xl(z,Zo,Uref,Ustar):
    #=IF(E16<>"",($O$42/0.4*(LN(D16)-LN($O$44)))/VLOOKUP($F$37,$A$16:$I$36,9),NA())
    l=pd.Series((Ustar/0.4*(np.log(z)-np.log(Zo)))/Uref)
    l=l.rename('U/Ur')
    return l
    
#   Turbulence
def TI_xl(z,Zo):
    #=IF(E16<>"",IF(D16<100,O$41*LN(30/(O$44))/(LN((E16/100*$E$3)/O$44)),($O$41*LN(30/$O$44)/LN(100/$O$44)+(D16-100)/500*(0.01-$O$41*LN(30/$O$44)/LN(100/$O$44)))),NA())
    n=nFromZo(Zo)
    # operaate on a series of z
    try:
        z=pd.Series(z)
        t1=n*np.log(30/Zo)/np.log(z[z<100]/Zo)
        t2=n*(np.log(30/Zo))/np.log(100/Zo)+(z[z>=100]-100)/500*(0.01-n*np.log(30/Zo)/np.log(100/Zo))
        t=t1.append(t2)
    # operate on a single z
    except:
        if z < 100:
            t=n*np.log(30/Zo)/np.log(z/Zo)
        else: 
            t=n*(np.log(30/Zo))/np.log(100/Zo)+(z-100)/500*(0.01-n*np.log(30/Zo)/np.log(100/Zo))
    return t*100
    
#%% plotting functions

#colorbar setup func
def cbScale(bounds):
    #bounds as a (low, high) tuple
    s=np.arange(bounds[0],bounds[1]+1,1)
    norm = colors.Normalize()
    norm.autoscale(s)
    sm = cm.ScalarMappable(cmap=cmap,norm=norm)
    sm.set_array([])
    return sm

#   mostly for QA
def Profile_plot(fetch, label, ID, results='', tun='ALL'):
    #for debugging
    #FetchID=df_og['Fetch'].dropna().drop_duplicates().reset_index(drop=True)
    #fetch=FetchID[27]  #for manual fetch selection
    
    #parse results info for given fetch
    try:
        R2pwr=results.iloc[ID,4].round(3)
        R2log=results.iloc[ID,5].round(3)
        R2ti=results.iloc[ID,6].round(3)
    except AttributeError:
        R2pwr=''
        R2log=''
        R2ti=''
    
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
    
    if tun=='GDE' or tun=='GDW':
        data=F_t(data,tun)
        
# =============================================================================
#     #limit data to specific locations
#     data=data[data.Location.isin([
#         'Ref',
#         '-3r/2',
#         '-r/2',
#         'TTC',
#         '+r/2',        
#             ])]
# =============================================================================
    
    # sort by abs distance from TT center, show most centered points at front
    data['abs_x']=abs(data['X(mm)']+1000)
    data=data.sort_values('abs_x', ascending=False)
    
    # separate horizontal search from vertical
    horiz=F_d(data,'Lateral Profile')
    data=F_d(data,'Vertical Profile')

    # z heights for ideal curves
    z = np.linspace(1,SF*1.2,num=1000)

    #### actual plotting, all above is just data manipulation/prep
    
    fig,ax = plt.subplots(1,2, sharey=True, figsize=(11,7)) #custom fig size
    #fig,ax = plt.subplots(1,2, sharey=True)
    
    # target U pwr - ADC
    n=nFromZo(Zo_targ)
    ax[0].plot(Pwr_Law(z,SF,Zo_targ), z , c='grey', label='Target Power, n={}, r^2={})'.format(round(n,2),R2pwr))
    
    # target U log - EXCEL Example
    Uref = data['Uref(m/s)'].mean()
    Us = Ustar_xl(Zo_targ,SF,Uref)
    ax[0].plot(Log_xl(z,Zo_targ,Uref,Us), z, c='grey', linestyle='-.', label='Target Log, U*={}, r^2={}'.format(round(Us,2),R2log))
    
    # target TI - EXCEL Example
    ax[1].plot(TI_xl(z,Zo_targ), z, c='grey', linestyle='--', label='Target TI, r^2={}'.format(R2ti))
    
    # horizontal profile setup
    horiz['Z_fsr']=round(horiz['Z_fs'],-1)
    hts=horiz['Z_fsr'].drop_duplicates()
    horiz['hf_avg']=''
    horiz['omni_avg']=''
    ps=SF*.75 # scalar to magnify percent variation
    
###vvv      single horizontal profile
    h = min(hts, key=lambda x:abs(x-100)) #finds value in hts closest to 100
    mean_hf=horiz[horiz['Z_fsr']==h]['Uhf/Uref'].mean()
    mean_omni=horiz[horiz['Z_fsr']==h]['Uomni/Uref'].mean()
    mean_combo=(mean_omni+mean_hf)/2
    horiz.loc[horiz['Z_fsr']==h, 'hf_avg']=mean_hf
    horiz.loc[horiz['Z_fsr']==h, 'omni_avg']=mean_omni
    horiz.loc[horiz['Z_fsr']==h, 'combo_avg']=mean_combo
    
    ax[0].axhline(y=(h+ps*.05), color='grey', linestyle='--', xmax=0.45)
    ax[0].axhline(y=(h-ps*.05), color='grey', linestyle='--', xmax=0.45)
###^^^      single horizontal profile

###vvv      enable to plot all horizontal profiles, messy for many data sets
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
###^^^      enable to plot all horizontal profiles, messy for some data sets

    # Horizontal profile data manipulation
    horiz['HFxvar']=0.2*horiz['Yhf(mm)']/horiz['Yhf(mm)'].max()+0.425           # normalize y position to go from 0 to .6
    horiz['Omnixvar']=0.2*horiz['Yomni(mm)']/horiz['Yomni(mm)'].max()+0.425
    horiz['HFyvar']=horiz['Z_fsr']+ps*(horiz['Uhf/Uref']-horiz['combo_avg'])/horiz['combo_avg'] 
    horiz['Omniyvar']=horiz['Z_fsr']+ps*(horiz['Uomni/Uref']-horiz['combo_avg'])/horiz['combo_avg']
    
    # just for single label, plots way off screen
    ax[0].axhline(y=-SF, color='grey', linestyle='--', label='+- 5% bounds')
    
    # horizontal lines for FS heights and tunnel setup
    ax[1].axhline(y=30, color='g', linestyle='-.',label='30m FS')
    ax[1].axhline(y=100, color='grey', linestyle='-.',label='100m FS')
    ax[1].axhline(y=block_ht, color='b', linestyle='--',label='Block, Trip, Spire Heights')
    ax[1].axhline(y=trip_ht, color='b', linestyle='--')
    ax[1].axhline(y=spire_ht, color='b', linestyle='--')
    
    # Manual colorbar setup
    bounds = (-3500,1200)
    sm=cbScale(bounds)
    cb=plt.colorbar(sm)
    cb.ax.set_title('Xpos (mm)', 
                    fontsize=9,
                    ha='left')
    
    # fake data for legend
    
    # East
    if tun=='GDE' or tun=='ALL':
        ax[0].scatter(x=-1, y=-SF, marker='x', c='k', label='GDE HF')
        ax[0].scatter(x=-1, y=-SF, marker='+', c='k', label='GDE 12HP')
    # West
    if tun=='GDW' or tun=='ALL':
        ax[0].scatter(x=-1, y=-SF, marker='s', edgecolors='k', facecolors='none', label='GDW HF')
        ax[0].scatter(x=-1, y=-SF, marker='^', edgecolors='k', facecolors='none', label='GDW 12HP')
    
    
    ### actual tunnel data
    # East
    if tun=='GDE' or tun=='ALL':
        E=F_t(data,'GDE')
        horiz_e=F_t(horiz,'GDE')
        
        color=sm.cmap(sm.norm(E['X(mm)']))
        color_h=sm.cmap(sm.norm(horiz_e['X(mm)']))
        
        # velocity data
        ax[0].scatter(x=E['Uhf/Uref'], y=E['Z_fs'], marker='x', c=color)
        ax[0].scatter(x=E['Uomni/Uref'], y=E['Z_fs'], marker='+', c=color)
        
        # turbulence data
        ax[1].scatter(x=E['TIUhf(%)'], y=E['Z_fs'], marker='x', c=color)
        ax[1].scatter(x=E['TIUomni(%)'], y=E['Z_fs'], marker='+', c=color)
        
        # lateral profile plots
        ax[0].scatter(x=horiz_e['HFxvar'], y=horiz_e['HFyvar'], marker='x', c=color_h)
        ax[0].scatter(x=horiz_e['Omnixvar'], y=horiz_e['Omniyvar'], marker='', c=color_h)
    
    
    # West
    if tun=='GDW' or tun=='ALL':
        W=F_t(data,'GDW')
        horiz_w=F_t(horiz,'GDW')
        color=sm.cmap(sm.norm(W['X(mm)']))
        color_h=sm.cmap(sm.norm(horiz_w['X(mm)']))
        
        # velocity data
        ax[0].scatter(x=W['Uhf/Uref'], y=W['Z_fs'], marker='s', edgecolors=color, facecolors='none')
        ax[0].scatter(x=W['Uomni/Uref'], y=W['Z_fs'], marker='^', edgecolors=color, facecolors='none')
        
        # turbulence data
        ax[1].scatter(x=W['TIUhf(%)'], y=W['Z_fs'], marker='s', edgecolors=color, facecolors='none')
        ax[1].scatter(x=W['TIUomni(%)'], y=W['Z_fs'], marker='^', edgecolors=color, facecolors='none')
        
        # lateral profile plots
        ax[0].scatter(x=horiz_w['HFxvar'], y=horiz_w['HFyvar'], marker='s', edgecolors=color_h, facecolors='none')
        ax[0].scatter(x=horiz_w['Omnixvar'], y=horiz_w['Omniyvar'], marker='^', edgecolors=color_h, facecolors='none')

    #legend, axis bounds, and title setup
    ax[0].legend()
    ax[1].legend()
    
    ax[0].set_xlabel('U/Uref')
    ax[1].set_xlabel('T.I.  (%)')
    ax[0].set_ylabel('Z F.S. (m)')
    
    ax[0].set_xlim(.2, 1.2)
    ax[1].set_xlim(0, 50)
    
    ax[0].set_ylim(0, 1.5*SF)
    
    plt.suptitle('{}:{} \n ID#{}:{} \n {}'.format(tun, targets, ID, fetch, label))
    
    # save the plot
    Plot_Path=os.path.join(os.getcwd(),'Plots','Include_TT_{}'.format(TT),tun,targets,'{}_{}_{}.png'.format(tun,targets,label))
    Plots_Folder=os.path.dirname(Plot_Path)
    if not os.path.exists(Plots_Folder): os.makedirs(Plots_Folder)
    plt.savefig(Plot_Path)
    
    #plt.close()

#%% error calculations

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2


def Calc_Error(fetch, tun):    
    #for debugging:
    #   FetchID=df_og['Fetch'].dropna().drop_duplicates().reset_index(drop=True)
    #   fetch=FetchID[10] 
    
    #filter for fetch and WT
    df=F_f(df_og,fetch)
    
    if tun=='GDW' or tun=='GDE':
        df=F_t(df,tun)
    
    #trim down data set to "good" window
    df=F_d(df,'Vertical Profile')
    df=df[df.Location.isin([
                            #'Ref',
                            '-3r/2',    #-2600 mm
                            '-r/2',     #-900 mm
                            #'TTC',
                            #'+r/2',        
                                ])]
    
    if not df.empty:
        df['Z_fs']=df['Z(mm)']*SF/1000
        
        #new column names for columns with annoying special characters
        new=['Ur',
            'Uh',
            'Uo',
            'Th',
            'To']
        
        old=['Uref(m/s)',
            'Uhf/Uref',
            'Uomni/Uref',
            'TIUhf(%)',
            'TIUomni(%)']
        
        df[new] = df[old]
        
        cols=['Tunnel',
            'Fetch',
            'Ur',
            'Location',
            'Z(mm)',
            'Z_fs',
            'Uh',
            'Uo',
            'Th',
            'To']

        # limit to "good" height region for each law    
        pwr=df[df['Z_fs'].between(20,SF)][cols].copy()
        log=df[df['Z_fs'].between(5,100)][cols].copy()
        ti=df[df['Z_fs'].between(5,SF)][cols].copy()
        
        # calc target values
        pwr['Targ']=Pwr_Law(pwr.Z_fs, SF, Zo_targ)
        log['Us']=Ustar_xl(Zo_targ,SF,log.Ur)
        log['Targ']=Log_xl(log.Z_fs, Zo_targ, log.Ur, log.Us)
        ti['Targ']=TI_xl(ti.Z_fs,Zo_targ)
        

        # calc errors HF only
        U_pwr_ErAvg=pd.DataFrame.mean(abs(pwr.Uh-pwr.Targ))
        U_log_ErAvg=pd.DataFrame.mean(abs(log.Uh-log.Targ))
        TI_ErAvg=pd.DataFrame.mean(abs(ti.Th-ti.Targ))
        # calc r^2 HF only
        U_pwr_R2=rsquared(pwr.Uh,pwr.Targ)
        U_log_R2=rsquared(log.Uh,log.Targ)
        TI_R2=rsquared(ti.Th,ti.Targ)

# =============================================================================
#         # calc errors Omni only
#         U_pwr_ErAvg=pd.DataFrame.mean(pwr.Uo-pwr.Targ)
#         U_log_ErAvg=pd.DataFrame.mean(log.Uo-log.Targ)
#         TI_ErAvg=pd.DataFrame.mean(ti.To-ti.Targ)
#         # calc r^2 OMNI only
#         U_pwr_R2=rsquared(pwr.Uo,pwr.Targ)
#         U_log_R2=rsquared(log.Uo,log.Targ)
#         TI_R2=rsquared(ti.To,ti.Targ)
# =============================================================================
        
# =============================================================================
#         # Avg of HF and Omni
#         pwr['Ua']=(pwr.Uh+pwr.Uo)/2
#         log['Ua']=(log.Uh+log.Uo)/2
#         ti['Ta']=(ti.Th+ti.To)/2
#         # calc errors AVG
#         U_pwr_ErAvg=pd.DataFrame.mean(pwr.Ua-pwr.Targ)
#         U_log_ErAvg=pd.DataFrame.mean(log.Ua-log.Targ)
#         TI_ErAvg=pd.DataFrame.mean(ti.Ta-ti.Targ)
#         # calc r^2 AVG
#         U_pwr_R2=rsquared(pwr.Ua,pwr.Targ)
#         U_log_R2=rsquared(log.Ua,log.Targ)
#         TI_R2=rsquared(ti.Ta,ti.Targ)
# =============================================================================

        out=[
            U_pwr_ErAvg, U_log_ErAvg, TI_ErAvg,
            U_pwr_R2, U_log_R2, TI_R2
            ]
        return out
    
    else: 
        out=[
            None, None, None,
            None, None, None
             ]
        return out
    
    
def Results(df_og, tun, targets):
    df=df_og.copy()
    results=pd.DataFrame()
    results['FetchID']=df['Fetch'].dropna().drop_duplicates().reset_index(drop=True)
    #df=df[df['Notes'].isna()] #drop any trials with notes, assume data is bad 
    # i think this causes problems
    
    #set column names as variables, easy to change later
    col=[   'U_pwr_ERR',
            'U_log_ERR',
            'TI_ERR',
            'U_pwr_r^2',
            'U_log_r^2',
             'TI_r^2']
    
    results[col]=None
    
    print('Calculating results...')
    
    for i in results.index:
        try:
            fetch=results.FetchID[i]
            results.loc[i,col]=Calc_Error(fetch,tun)
        except:
            print('Something went wrong with FetchID:{}'.format(i))
            pass
        
    # save results
    print('Saving results...')
    Results_Path=os.path.join(os.getcwd(),'Results','Include_TT_{}'.format(TT),tun,'{}_results.csv'.format(targets))
    Results_Folder=os.path.dirname(Results_Path)
    if not os.path.exists(Results_Folder): os.makedirs(Results_Folder)
    results.to_csv(Results_Path)
    
    return results
    
def Ranked_Plots(results):
    #columns=results.columns
    #col=columns[~columns.isin(['FetchID','Score'])]
    
    print('Generating plots...')
    
    n=3 #number of plots per set
    
    ###vvv USE THIS FOR MEAN ERROR RANKING
    results['Score']=results.U_pwr_ERR*results.U_log_ERR*(results.TI_ERR**2)
    
    ranked_pwr=results.sort_values(by='U_pwr_ERR',ascending=True).copy().reset_index(drop=True)
    top_pwr=ranked_pwr.FetchID[0:n]
    
    ranked_log=results.sort_values(by='U_log_ERR',ascending=True).copy().reset_index(drop=True)
    top_log=ranked_log.FetchID[0:n]    
    
    ranked_TI=results.sort_values(by='TI_ERR',ascending=True).copy().reset_index(drop=True)
    top_TI=ranked_TI.FetchID[0:n]
    
    ranked_score=results.sort_values(by=['Score'],ascending=True).copy().reset_index(drop=True)
    top_score=ranked_score.FetchID[0:n]
    
    for i in range(n):
        f_pwr=top_pwr[i]
        label='Best Power Law fit #{}'.format(i+1)
        ID=results[results.FetchID==f_pwr].index[0]
        Profile_plot(f_pwr, label, ID, results, tun)
        
        f_log=top_log[i]
        label='Best Log Law fit #{}'.format(i+1)
        ID=results[results.FetchID==f_log].index[0]
        Profile_plot(f_log, label, ID, results, tun)
        
        f_ti=top_TI[i]
        label='Best TI fit #{}'.format(i+1)
        ID=results[results.FetchID==f_ti].index[0]
        Profile_plot(f_ti, label, ID, results, tun)
        
        f_s=top_score[i]
        label='Best overall #{}'.format(i+1)
        ID=results[results.FetchID==f_s].index[0]
        Profile_plot(f_s, label, ID, results, tun)        
    ###^^^ USE THIS FOR ERROR RANKING

# this needs fixing if used
# =============================================================================
#     ###vvv USE THIS FOR R^2 RANKING
#     results['Score']=results['U_pwr_r^2']*results['U_log_r^2']*(results['TI_r^2']**2)
#     ranked_pwr=results.sort_values(by='U_pwr_r^2',ascending=False).copy().reset_index(drop=True)
#     top_pwr=ranked_pwr.FetchID[0:n]
#     for i in range(n):
#         f=top_pwr[i]
#         label='Best Power Law fit #{}'.format(i+1)
#         ID=results[results.FetchID==f].index[0]
#         Profile_plot(f, tun, label, ID, results)
#         
#     ranked_log=results.sort_values(by='U_log_r^2',ascending=False).copy().reset_index(drop=True)
#     top_log=ranked_log.FetchID[0:n]
#     for i in range(n):
#         f=top_log[i]
#         label='Best Log Law fit #{}'.format(i+1)
#         ID=results[results.FetchID==f].index[0]
#         Profile_plot(f, tun, label, ID, results)
#         
#     ranked_TI=results.sort_values(by='TI_r^2',ascending=False).copy().reset_index(drop=True)
#     top_TI=ranked_TI.FetchID[0:n]
#     for i in range(n):
#         f=top_TI[i]
#         label='Best TI fit #{}'.format(i+1)
#         ID=results[results.FetchID==f].index[0]
#         Profile_plot(f, tun, label, ID, results)
#         
#     ranked_score=results.sort_values(by=['Score'],ascending=False).copy().reset_index(drop=True)
#     top_score=ranked_score.FetchID[0:n]
#     for i in range(n):
#         f=top_score[i]
#         label='Best overall #{}'.format(i+1)
#         ID=results[results.FetchID==f].index[0]
#         Profile_plot(f, tun, label, ID, results)
#     ###^^^ USE THIS FOR R^2 RANKING
# =============================================================================


#%% Get file & initialize - should make this a function later

#hard-coded filepath of OG data (Cross-network)
#filepath='//gd-fs3/CentralDataCache/Projects/G153E/Results_Post_Processed/BLhunts/2020-06-GoldPlateOutput1/GoldPlateBrief.xlsx'

#local copy - set working directory to location of this script first
filepath='./GP_Data\GoldPlateBrief.xlsx'

#path of pickled (saved) dataframe
pklpath='./GP.pkl'

if os.path.isfile(pklpath):         #check if .pkl exists already   
    LoadSaved=easygui.ynbox("Saved data found, do you want to load the saved data?")
    if LoadSaved:
        print('Loading saved dataframe...')
        df_og=pd.read_pickle(pklpath)
    elif LoadSaved==None: sys.exit()    
    elif LoadSaved==False:
        if os.path.isfile(filepath) is False:   #open file dialogue if .xlsx path not found
            print('Default file not found')
            filepath=Get_Data()
        print('Reading data file...')
        df_og = pd.read_excel(filepath, header=0)
        print('Saving dataframe...')
        df_og.to_pickle(pklpath)                #save data to .pkl file
else:                                       #if no .pkl, try to load data
    if os.path.isfile(filepath) is False:   #open file dialogue if .xlsx path not found
        print('Default file not found')
        filepath=Get_Data()
    print('Reading data file...')
    df_og = pd.read_excel(filepath, header=0)
    print('Saving dataframe...')
    df_og.to_pickle(pklpath)


#TT=False
        
### single roughness                
# get inputs, in future should turn into an array of values to loop thru
tun,SF,Zo_targ,TT=Inputs()
targets='{}-{}'.format(SF,Zo_targ)            
###

### pre-process fetch ID to ignore TT roughness
if not TT:
    df_og['fetch_og']=df_og.Fetch.copy()
    for i in df_og.index:
        #print(i)
        s=df_og['fetch_og'][i]
        try:
            df_og.loc[i,'Fetch']=s[:s.find('TT')].strip()
        except AttributeError:
            pass
###

#%% Begin real script
t0 = datetime.datetime.now()    #start process timer
inidir=os.getcwd()


#set plot style 
mpl.rcdefaults()            #reset to defaults
styles=plt.style.available  #save all plot styles to  list
plt.style.use(styles[14])   #set style
cmap='twilight'             #set colormap
#cmap='twilight_shifted'    #if plot style has dark background

### single roughness
results=Results(df_og, tun , targets)
if easygui.ynbox('Generate Plots?'):
    Ranked_Plots(results)
###
  

# =============================================================================
# ####vvv         enable to loop through array of targets
# 
# tunnel=['GDW','GDE']
# #tunnel= "ALL"
# Scale=[180,240,300]
# Zo=[0.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
# 
# for tun in tunnel:
#     for SF in Scale:
#         #print(SF)
#         for Zo_targ in Zo:
#             #print(Zo_targ)
#             targets='{}-{}'.format(SF,Zo_targ)
#             print(targets)
#             results=Results(df_og, tun , targets)
#             Ranked_Plots(results)
# 
# ####^^^         enable to loop through arrays of tun, scale, Zo
# =============================================================================

#%% Done
print('Done!')
t1=datetime.datetime.now()
dt= t1-t0
dt=dt.seconds
easygui.msgbox(msg="Done!\n Process took: {} seconds\nPress ok to view plots".format(dt))  
plt.show()
