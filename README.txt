### End-to-end Handling of NICER Data (EnHaNCED)
This set of codes is supposed to make NICER data analysis more convenient.
You need to provide a whole bunch of inputs in 'inputs_enhanced.py'.
Defaults should work for all but the first three inputs.
The default inputs assume that the codes are in the same directory that contains the ObsID folders.
But you can change the inputs appropriately in case you don't want to keep the codes and the ObsID folders in the same path.

****** Important ***********************
You must enable heasoft before running the code, and have the 3C50 background model installed.

*****************

************* INPUT explanation (the self-explanatory inputs are not included) **********
obsid_f		: The name of the file that contains the ObsIDs. 
         	e.g. obsid_f             ='my_path/ObsID.dat'
bkg_dir		: The director containing the 3C50 background model informatiom.
         	see https://heasarc.gsfc.nasa.gov/docs/nicer/tools/README_nibackgen3C50_v7b.txt
clean_type	: Whether you want to do the default filtering (clean_type='_default'), or default+manual filtering.
            	If you opt to do manual flagging, you can either do underonly based flagging (not well tested yet, clean_type='_flagged_underonly'), or 
            	12-15 keV based flagging (clean_type='_flagged_high_energy', NICER response is ideally 0 over 12-15 keV), or 
            	both (clean_type='_flagged_high_energy_flagged_underonly').
            	Note that specifying clean_type will not automatically lead to the selection of a given mode of filtering.
merge_events	: The need to merge some of the ObsIDs may arise under certain circumstances (for example, to merge all the ObsIDs to obtain an average event).
              	If you set it to True, you must also provide inputs for merge_obs_arr (e.g. merge_obs_arr       =['3627010401','3627010402']), and
	      	out_ID, which is the name of the output ObsID (e.g. out_ID='merged_event').
              	Note that, right now, you can form only one 'merged event'at one go.
update_obs_list	: This is when you want to replace the ObsIDs which you merged to form the out_ID, with the out_ID. For example, let's say your obs_id=['1','2','3','4'], 
                 and you merged '3' and '4' to '34'. update_obs_list =True will change obs_id to obs_id=['1','2','34']
filter_type_arr	: Do not change it unless you know what you are doing
model_type	:This can be anything, but try to give a name that reflects the model that you will be using.
multiresponse	: This option is not very general one. It allows one to model the total spectrum by modeling both the target and background spectra together.
inp_flags	: This allows you to pass custom inputs to nicerl2. The inputs for nicerl2 can be found in 
		https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerl2.html
obs_dir		: Path containing the ObsID folders
primary_dir	: Path containing the codes


************* OUTPUTS *******************
1) Cleaned event files	: can be found in the obsid/xti/event_cl. The cleaned event file will have names like niobsID_0mpu7_cl_default.evt, 
			where '_default' is the clean_type
2) Lightcurve and spectra: Outputs of xselect. The ufa (unfiltered) and cleaned outputs will be in obs_dir/ufa and obs_dir/cl respectively. 
			  Inside each of them, there will be four directories (at max): full (0-15 keV), low (0.3-0.4 keV), mid (0.4-2.0 keV) and high (12-15 keV).
                          The most important outputs will be in obs_dir/cl/full
3) Background spectra	: These will be in obs_dir/cl/full
4) Responses (arf and rmf): These will be in obs_dir/cl/full
5) Grouped spectra	: These will be in obs_dir/cl/full
6) MCMC outputs		: These will be in obs_dir/cl/full
7) Text file containing norms, fluxes:  These will be in obs_dir/cl/full



