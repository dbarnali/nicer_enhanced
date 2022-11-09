### End-to-end Handling of NICER Data (EnHaNCED)
import os
import numpy as np
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


#### INPUTS that you must provide
'''
The following three inputs must be provided accurately
for the code to run
'''
obsid_f             ='ObsID_small.dat'    # The name of the file that contains the list of ObsIDs
bkg_dir             ='/home/dbarnali/heasoft/3c50_bgmodel/bg_models_3C50'   ### The directory with the 3c50 model info
clean_type          ='_default'#     ###### VERY IMPORTANT INPUT ######
                        #If you opt to do manual flagging, you can either do underonly based flagging (not well tested yet), or 12-15 keV based flagging, or both
                        #Since, I don't want to overwrite any files, I will use a new 'clean_type' to add to the out names for these different schemes
                        #For default, clean_type='_default', for underonly based flagging, clean_type='_flagged_underonly', for high energy based flagging, clean_type='_flagged_high_energy'
                        #For both, clean_type='_flagged_high_energy_flagged_underonly'


#### Where should I start?
'''
Carefully go through these inputs, they tell the code what it is supposed to do
The default should work in the first run
'''
run_nicerl2         =True   #Run the nicerl2 pipeline
input_custom_gti    =False  #if True, I will expect a file named obsid_time.txt inside each obsid folder
include_gti         =False  #This is used only if input_custom_gti is set True, 
                            #it says whether the selected GTIs in obsid_time.txt is to be excluded (False) 
                            #or those are the only GTIs that are to be included (True)                                    
extract_cl_lc_pha   =True   #Extract the lc and spectrum for the cleaned files   
extract_ufa_lc_pha  =True   #Extract the lc and spectrum for the ufa files
make_diagnostic_plot=True   #if set to True, it will create some diagnostic plots
do_flag             =False   #Manual flagging, if set to True, it will create some diagnostic plots,
                                    #and then ask you to decide whether you indeed want manual flagging based on those plots
merge_events        =False  #If you want to merge some ObsIDs, 
                                    #if True, you will have to provide the list of ObsIDs that you want to merge (merge_obs_arr),
                                    #and also the name of the output ObsID for the merged event (out_ID)
update_obs_list     =False  #Set to True if you have merged some of the ObsIDs in the original ObsID file given in obsid_f,
                            #and want to use the merged event instead of the individual ObsIDs that have been merged
extract_bkg         =True   #Extract the background spectra
extract_arf_rmf     =True   #Extract the arf and rmf
grp_pha             =True   #Group the spectrum so as to avoid having too few counts in any energy bin
do_spectral_analysis=True     

######### INPUTS for xselect########
time_bin_xselect    =25
min_chan_xselect    =0
max_chan_xselect    =1500
######## Relevant only if merge_events=True or update_obs_list   =True
merge_obs_arr       =['3627010401','3627010402']#['3627010101','3627010301','3627010501','3627010901']#['3627010101','3627010301','3627010501','3627010701','3627010901']#
out_ID              ='401_402'#'periastron_no_07'#'merged_periastron'#
filter_type_arr     =['mkf']
######### INPUTS for grouping########
min_count           =15
min_chan_num        =26
max_chan_num        =1000
####################################
# INPUTS for SPECTRAL ANALYSIS
'''
This code can only handle models like tbabs(apec+apec+....).
Do not use this code if you want to use different models.
Alernatively, you can edit the functions in spectral_fitting_modular.py.
It will perform MCMC to get constraints on the fitting parameters.
'''
model_type      ='2apec'                #a name of your choice to add to the output files
kT_arr          =np.array([0.1,0.8])    #Guess value of kTs
freeze_kT_arr   =np.array([True,True])#Do you want to freeze any of the kTs?
startE,endE     =0.3, 10.0              #Energy range to be used in spectral analysis       
freeze_nH       =True                  #Freeze nH in the tbabs model?
nH0             =0.01                   #Guess value for the nH parameter of the tbabs model
ism_correct     =True                   #Do you want the ISM absorption corrected flux values
startE_arr      =[0.5,2.0]              #Energy ranges over which you want to calculate the fluxes
endE_arr        =[2.0,10.0]
multiresponse   =False                  #This is not a general input, better set it to False 
input_wt        ='standard'             
chain_length    =200000                 #Inputs for the MCMC analysis, if the spectral part is taking too much time
                                        #reduce the number in chain_length
burn_length     =2000                   
####################################################
# Miscellanous inuts, go through these, but defaults should be okay
inp_flags       ='overonly_range=0-1 niprefilter2_coltypes=base,3c50 nicersaafilt=YES saafilt=NO' #HARD-CODED for 3C50 model
obs_dir             =os.getcwd()    #where all the obs_id folders are present
diagnostic_plot_path=obs_dir #This is where the diagnostics plots will be saved
primary_dir         =os.getcwd()        #the directory where the code is located
if multiresponse==True:     
    model_type+='_multiresponse_bkg_pow_gaussian_fixed_mu_sigma'
############################


