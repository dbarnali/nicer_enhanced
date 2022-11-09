### End-to-end Handling of NICER Data (EnHaNCED)
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
from scipy import interpolate
from scipy.optimize import curve_fit
from xspec import *
import glob
import spectral_fitting_modular as spec_fit
import nicer_functions as ni_func
from input_enhanced import *
# Module that has the ephemeric information in it.


'''
PLEASE KEEP THIS CODE AND THE ACCOMPANYING CODES IN THE SAME DIRECTORY THAT CONTAINS THE OBSID FOLDERS
You must enable heasoft before running this code
It assumes the default NICER data structure
e.g. there are event files at obsid/xti/event_cl
It assumes there exists a filename containing the list of ObsIDs
Right now it is hard-coded to use the 3C50 background model

OUTPUT: The diagnostic plots will be saved in the directory given by diagnostic_plot_path
        The lc, pha, and other kinds of outputs (mcmc fits file, mcmc plot, spectra, norm values, flux tables) will be saved in path_to_code/cl/full
'''

mpl.rcParams['font.size']           =12
mpl.rcParams['axes.linewidth']      =1.5
mpl.rcParams['xtick.major.size']    =6
mpl.rcParams['xtick.minor.size']    =3.5
mpl.rcParams['xtick.major.width']   =1.2
mpl.rcParams['xtick.minor.width']   =1.1
mpl.rcParams['xtick.direction']     ='in'
mpl.rcParams['ytick.major.size']    =6
mpl.rcParams['ytick.major.width']   =1.2
mpl.rcParams['ytick.minor.width']   =1.1
mpl.rcParams['ytick.direction']     ='in'
mpl.rcParams['xtick.top']           =True
mpl.rcParams['ytick.right']         =True

bkg_type        ='3C50'     ### MEANT TO BE HARD-CODED, DO NOT CHANGE IT

obsid = np.loadtxt(obsid_f, dtype=str)
if len(np.shape(obsid))==0: #if there is just one observation ID
    obsid   =np.array([obsid])

if input_custom_gti==True:
    folder          =obs_dir
    ni_func.nimaketime_command_create(folder,obsid, include=include_gti, run_f='run_nimaketime.sh')
    os.chdir(folder)
    os.system('bash run_nimaketime.sh')
os.chdir(primary_dir)

if run_nicerl2==True:
    ni_func.nicerl2_command_create(obsid, output=obs_dir+'/nicerl2_command.sh',custom_gti=input_custom_gti)  
    os.environ["TYPE"] = clean_type.split("_",1)[-1]
    os.environ["FLAGS"] = inp_flags
    os.system('bash nicerl2_command.sh')




if extract_cl_lc_pha==True:
    is_make =os.path.isdir(obs_dir+'/cl')
    if is_make==False:
        os.mkdir(obs_dir+'/cl')
    extract_mid_energy_lc_pha_cl   =True
    extract_low_energy_lc_pha_cl   =False#True
    extract_high_energy_lc_pha_cl  =True
    extract_full_energy_lc_pha_cl  =True
else:
    extract_mid_energy_lc_pha_cl   =False
    extract_low_energy_lc_pha_cl   =False
    extract_high_energy_lc_pha_cl  =False
    extract_full_energy_lc_pha_cl  =False

if extract_ufa_lc_pha==True:
    is_make =os.path.isdir(obs_dir+'/ufa')
    if is_make==False:
        os.mkdir(obs_dir+'/ufa')
    extract_mid_energy_lc_pha_ufa   =True
    extract_low_energy_lc_pha_ufa   =False#True
    extract_high_energy_lc_pha_ufa  =True
    extract_full_energy_lc_pha_ufa  =True
else:
    extract_mid_energy_lc_pha_ufa   =False
    extract_low_energy_lc_pha_ufa   =False
    extract_high_energy_lc_pha_ufa  =False
    extract_full_energy_lc_pha_ufa  =False

low_en_chan1                =0
low_en_chan2                =30
mid_en_chan1                =40
mid_en_chan2                =200
high_en_chan1               =1200
high_en_chan2               =1500

chan1_arr       =[low_en_chan1,mid_en_chan1,high_en_chan1]
chan2_arr       =[low_en_chan2,mid_en_chan2,high_en_chan2]

ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_low_energy_lc_pha_cl,'cl','low',low_en_chan1,low_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_mid_energy_lc_pha_cl,'cl','mid',mid_en_chan1,mid_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_high_energy_lc_pha_cl,'cl','high',high_en_chan1,high_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_full_energy_lc_pha_cl,'cl','full',low_en_chan1,high_en_chan2)

ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_low_energy_lc_pha_ufa,'ufa','low',low_en_chan1,low_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_mid_energy_lc_pha_ufa,'ufa','mid',mid_en_chan1,mid_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_high_energy_lc_pha_ufa,'ufa','high',high_en_chan1,high_en_chan2)
ni_func.run_xselect_custom_energy_range(obsid,obs_dir,clean_type,time_bin_xselect,extract_full_energy_lc_pha_ufa,'ufa','full',low_en_chan1,high_en_chan2)

os.chdir(primary_dir)


###########################################3

if do_flag==True or make_diagnostic_plot==True:
    ufa_lc_type_trigger=[False,True,True,True]#[extract_low_energy_lc_pha_ufa,extract_mid_energy_lc_pha_ufa,extract_high_energy_lc_pha_ufa,extract_full_energy_lc_pha_ufa]
    cl_lc_type_trigger=[False,True,True,True]#[extract_low_energy_lc_pha_cl,extract_mid_energy_lc_pha_cl,extract_high_energy_lc_pha_cl,extract_full_energy_lc_pha_cl]
    ni_func.make_diagnostic_lc(obsid,obs_dir,'.mkf',clean_type,ufa_lc_type_trigger,cl_lc_type_trigger,obs_dir,chan1_arr,chan2_arr,merge_gap=300)
    if do_flag==True:
        print('\nPlease check the diagnostic plots and let us know if you want to do flagging or not.\n')
        ans =input('Do you want to flag the data? (press 1 for yes and 0 for no) ')
    else:
        ans='0'
else:
    ans='0'

if ans=='1': 
    underonly_flag  =input('Do you want to flag based on underonly count lightcurves? (press 1 for yes and 0 for no, recommended is 0) ')
    if underonly_flag=='1':
        max_count   =float(input('Enter the maximum underonly count above which to clip (recommended is 15) '))
        grad_max_count=float(input('Enter the maximum gradient of the underonly count above which to clip (recommended is 1) '))
    high_en_flag    =input('Do you want to flag based on the high energy (12-15 keV) lightcurves? (press 1 for yes and 0 for no) ')  
else:
    do_flag=False

while do_flag==True:
    new_clean_type=''
    if high_en_flag=='1':
        new_clean_type+='_flagged_high_energy'
        lc_path =obs_dir+'/cl/high'
        ni_func.identify_bad_gti_lc(obs_dir,obsid,lc_path,clean_type)
        folder=obs_dir
        ni_func.nimaketime_command_create(folder,obsid,False,'run_nimaketime_high_en.sh')
        os.chdir(folder)
        os.system('bash run_nimaketime_high_en.sh')

        os.chdir(primary_dir)        
        ni_func.nicerl2_command_create(obsid, output=obs_dir+'/nicerl2_command.sh',custom_gti=True)
        os.environ["TYPE"] = new_clean_type.split('_',1)[-1]
        os.environ["FLAGS"] = inp_flags
        os.system('bash nicerl2_command.sh')
        
        
    if underonly_flag=='1': #!!!!! YET TO BE TESTED !!!!!
        new_clean_type+='_flagged_underonly'
        identify_bad_gti_underonly(obspath=obs_dir,obsid_file=obsid_f,mkf_ext='mkf',max_count=15,grad_max_count=1)
        folder=obs_dir
        ni.nimaketime_command_create(folder,obsid_f,False,'run_nimaketime_high_en.sh')
        os.chdir(folder)
        os.system('bash run_nimaketime_high_en.sh')

        os.chdir(primary_dir)
        ni_func.nicerl2_command_create(obsid, output=obs_dir+'/nicerl2_command.sh',custom_gti=True)
        os.environ["TYPE"] = new_clean_type.split('_',1)[-1]
        os.environ["FLAGS"] = inp_flags
        os.system('bash nicerl2_command.sh')
        

    
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'cl','low',low_en_chan1,low_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'cl','mid',mid_en_chan1,mid_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'cl','high',high_en_chan1,high_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'cl','full',low_en_chan1,high_en_chan2)

    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'ufa','low',low_en_chan1,low_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'ufa','mid',mid_en_chan1,mid_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'ufa','high',high_en_chan1,high_en_chan2)
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'ufa','full',low_en_chan1,high_en_chan2)

    ufa_lc_type_trigger=[False,True,True,True]#[extract_low_energy_lc_pha_ufa,extract_mid_energy_lc_pha_ufa,extract_high_energy_lc_pha_ufa,extract_full_energy_lc_pha_ufa]
    cl_lc_type_trigger=[False,True,True,True]#[extract_low_energy_lc_pha_cl,extract_mid_energy_lc_pha_cl,extract_high_energy_lc_pha_cl,extract_full_energy_lc_pha_cl]

    ni_func.make_diagnostic_lc(obsid,obs_dir,'.mkf',new_clean_type,ufa_lc_type_trigger,cl_lc_type_trigger,obs_dir,chan1_arr,chan2_arr,merge_gap=300)
    os.chdir(primary_dir)
    print('\nPlease check the diagnostic plots and let us know if you want to do flagging or not.\n')
    ans =input('Do you want to flag the data? (press 1 for yes and 0 for no) ')

    if ans=='1':
        do_flag=True
    else:
        do_flag=False
    clean_type=new_clean_type

new_clean_type  =clean_type   

####### MERGE some EVENTS

if merge_events==True:
    ni_func.merge_events(merge_obs_arr,out_ID,new_clean_type,obs_dir,filter_type_arr)
######### UPDATE the obsid arr
if update_obs_list==True:
    obsid   =list(obsid)
    for i in range(len(merge_obs_arr)):
        obsid.remove(merge_obs_arr[i])
    obsid.append(out_ID)
obsid   =np.array(obsid)

####### Now ready for spectral analysis
####### Will be done for the full spectra
######################    
pha_dir ='cl/full'
out_dir =obs_dir+'/'+pha_dir
if extract_bkg==True:
    ni_func.get_bkg_3c50(obsid,obs_dir,bkg_dir,out_dir,new_clean_type)

if extract_arf_rmf==True:
    ni_func.get_arf_rmf_3c50_tot_spec(obsid,obs_dir,pha_dir,new_clean_type)
#pha_dir should be inside the obs_dir
    ni_func.run_xselect_custom_energy_range(obsid,obs_dir,new_clean_type,time_bin_xselect,True,'cl','full',low_en_chan1,high_en_chan2)

os.chdir(primary_dir)
if grp_pha==True:
    ni_func.perform_optimal_binning(obsid,new_clean_type,out_dir,min_count,min_chan_num,max_chan_num)
os.chdir(primary_dir)


model_type_arr  =model_type.split('_')
if do_spectral_analysis==True: 
    if multiresponse==False:
        if 'powerlaw' in model_type_arr:
            spec_fit.spectral_analysis_apec_pow(obsid,new_clean_type,out_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,chain_length,burn_length,input_wt,out_dir,min_count,save_fig=True,write_norm=False)
            spec_fit.calculate_flux_apec_pow(ism_correct,obsid,new_clean_type,out_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,input_wt,out_dir,min_count)
        else:
            spec_fit.spectral_analysis_apec(obsid,new_clean_type,out_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,chain_length,burn_length,input_wt,out_dir,min_count,save_fig=True,write_norm=True)
            spec_fit.calculate_flux_apec(ism_correct,obsid,new_clean_type,out_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,input_wt,out_dir,min_count)
    else:
        ni_func.perform_optimal_binning_for_bkg(obsid,new_clean_type,out_dir,min_count,min_chan_num,max_chan_num)
        spec_fit.spectral_analysis_apec_multiresponse(obsid,new_clean_type,out_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,chain_length,burn_length,input_wt,out_dir,min_count,save_fig=True,write_norm=True)
        spec_fit.calculate_flux_apec_multiresponse(ism_correct,obsid,new_clean_type,out_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE,endE,bkg_type,input_wt,out_dir,min_count)
    


os.chdir(primary_dir)

