import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import os
from nicergof.bkg import bkg_estimator as be
import matplotlib as mpl
from astropy.time import Time
from astropy.time import TimeDelta
from scipy import interpolate
from scipy.optimize import curve_fit
from xspec import *
import glob
import time
import corner
from functools import reduce
import operator
from scipy.stats import median_abs_deviation
from PyAstronomy import pyasl


def nicerl2_command_create(obsid, output='nicerl2_command.sh',custom_gti=False):
    """
    Creates a bash script that contains command to run Nicerl2 for each ObsID in the list (read from filename='ObsID.dat'.
    
    The options for the Nicerl2 are listed as ${FLAGS} in the script. By defining this environement variable, the same script can be used for different processing options.
    
    The cleaned event list are renamed to cl_${TYPE}, where the environement variable ${TYPE} contains the desired suffix. This is so that the Nicerl2 can be used with different options, for comparison (the default behavior of Nicerl2 is to clobber).
    
    See the nicerl2_source.sh template to see how to set the environement variables.
    
    :param filename: ('ObsID.dat'), the filename that contains the list of ObsID.
    :param output: ('nicerl2_command.sh'), the filename of the output bash script
    :rtype: None. Creates a bash script called output='nicerl2_command.sh'
    
    """
    
    f = open( output  ,'w')
    
    for obs in obsid:
        if custom_gti==True and os.path.isfile(obs+'/ni'+obs+'.gti')==True:
            f.write( 'nicerl2 indir={} clobber=YES gtifiles={}/ni{}.gti ${{FLAGS}}\n'.format(obs,obs,obs)  )
        else:
            f.write( 'nicerl2 indir={} clobber=YES ${{FLAGS}}\n'.format(obs)  )
        f.write( 'mv {}/xti/event_cl/ni{}_0mpu7_cl.evt {}/xti/event_cl/ni{}_0mpu7_cl_${{TYPE}}.evt\n\n'.format(obs, obs, obs, obs) )

            
    f.close()

def xselect_command_create(folder,obsid, template_f='template.dat', run_f='run.sh'):
    '''
    This function reads in a xselect template and substitute the ObsID numbers
    
    In the template, the obsID is replaced with {}.
    The template is expected to be located in the specified folder.
    A xselect .xco script will be created for each ObsID.
    Also, a run file (default 'run.sh') will be created, to execute all of the .xco script at once, with a `bash run.sh`
    
    :param folder: the folder containing the template, and that will contain the created scripts
    :param filename: ('ObsID.dat') name of the file containing the list of ObsID
    :param template_f: ('template.dat') name of the file containing the template to use (expected location in `folder`)
    :param run_f: ('run.sh') name of the file containing all of the run commands (created in `folder`)
    '''

    template_file = open(folder + '/' + template_f, 'r')
    template = template_file.readlines()
    template_file.close()

    run = open( folder + '/'+ run_f, 'w')

    for obs in obsid:

        f = open( '{}/ni{}.xco'.format(folder,obs)  ,'w')
        run.write('xselect @ni{}.xco\n'.format(obs))
        
        for line in template:
        
            # This will place the obs in any {} in the template.
            # It there isn't any, it just returns the string.
            # Also here, I do not need a /n to break the line,
            # as it is already included in the string from reading the template file
            f.write( line.format(obs, obs, obs, obs, obs, obs, obs, obs) )
            # MAKE SURE that there are more "obs" in there than the mak # of {} in one line.
            
        f.close()
        
    run.close()


def create_xselect_template_file(folder,file_type,template_name,obs_dir,clean_type,time_bin_xselect,min_chan_xselect,max_chan_xselect):
    f=open(folder+'/'+template_name,'w')
    f.write('{}\nread events\n')
    f.write(obs_dir+'/{}/xti/event_cl\n')
    if file_type=='cl':
        f.write('ni{}_0mpu7_cl'+clean_type+'.evt\nyes\n')
    if file_type=='ufa':
        f.write('ni{}_0mpu7_ufa.evt\nyes\n')
    f.write('set bin '+str(time_bin_xselect)+'\n')
    f.write('set phaname PI\n')
    f.write('filter pha_cut\n')
    f.write(str(min_chan_xselect)+'\n'+str(max_chan_xselect)+'\nextract curve\nsave curve\nni{}'+clean_type+'\n')
    f.write('extract spectrum\nsave spectrum\nni{}'+clean_type+'\nexit\nno')
    f.close()


def rename_lc_pha(obsid,clean_type,file_type):
    for obs in obsid:
        filename    ='ni'+obs+clean_type+'.'+file_type
        if os.path.isfile(filename)==True:
            if os.path.isfile(filename.split('.')[0]+'_previous.'+file_type)==True:
                os.system('rm '+filename.split('.')[0]+'_previous.'+file_type)
            os.system('rename s/.'+file_type+'/_previous.'+file_type+'/ '+filename)

def run_xselect_custom_energy_range(obsid,obs_dir,clean_type='_default',time_bin_xselect=25,trigger=True,file_type='cl',en_type='mid',chan1=40,chan2=200):
    '''
    xselect will give an error if there already exists a file with the same name that it wants to create.
    Any *.lc and *.pha files in the destination directory will hence be renamed by as original_name_previous.lc/original_name_previous.pha
    '''
    if trigger==True:
        template_f    ='template_lc_pha_'+en_type+'.dat'
        my_dir  =obs_dir+'/'+file_type+'/'+en_type
        is_make =os.path.isdir(my_dir)
        if is_make==True:
            os.chdir(my_dir)
            rename_lc_pha(obsid,clean_type,'lc')
            rename_lc_pha(obsid,clean_type,'pha')
            
        if is_make==False:
            os.mkdir(my_dir)
        create_xselect_template_file(my_dir,file_type,template_f,obs_dir,clean_type,time_bin_xselect,chan1,chan2)
        xselect_command_create(my_dir,obsid, template_f=template_f, run_f='run.sh')
        os.chdir(my_dir)
        os.system('bash run.sh')
        return
    else:
        return

def load_curve(filename, T0=0):
    """
    Load a light curve from Xselect.
    
    The light curve starts from the first GTI.
    
    :param filename: The name of the file containing the Xselect LC.
    :param T0: (0) Difference between the start of the GTI in the Xselect LC, and that of the desired "zero" time stamp.
    :rtype time_LC: = Xselect Time stamp + start of Xselect first GTI - T0.
    :rtype LC: The light curve
    
    """
    # Open the light curve
    hdul = fits.open(filename)
    LC = hdul[1].data
    gti_LC = hdul[2].data # get the gti
    hdul.close()
    # The light curve starts from the first GTI.
    # Need to add it back then reset to our T0
    time_LC = LC['TIME']-T0+gti_LC[0][0]
    
    return(time_LC, LC)

def gti_info(gti, verbose=False):
    """
    Returns the duration of each GTI and the gap in between each GTI
    
    :param gti: a GTI object read from the GTI fits extension.
    :param verbose: (False) flag to print out the GTI information
    :rtype duration: array with the duration of each GTI (in sec)
    :rtype gap: array with the gap (in sec) between each GTI
    
    """

    n = gti.shape[0]
    duration = gti[:,1] - gti[:,0]
    gap = gti[1:,0] - gti[:-1,1]
    gap = np.append(0,gap)
    
    if verbose:
        print('{:10} {:10} {:10} {:10}'.format('start', 'stop', 'duration', 'gap'))
        for i in range(0,n):
            print('{:10.1f} {:10.1f} {:10.1f} {:10.1f}'.format(gti[i,0], gti[i,1], duration[i], gap[i] ))

    return(duration, gap)

def make_diagnostic_lc(obsid,obspath,mkf_ext='.mkf',clean_type='_default',ufa_lc_type_trigger=[False,True,True,True],cl_lc_type_trigger=[False,True,True,True],PDF_path='.',chan1_arr=[0,40,1200],chan2_arr=[30,200,1500],merge_gap=300):
    '''
    ufa_lc_type_trigger is an array of four booleans. 
    ufa[0]=True will trigger the plotting for the low
    ufa[1]=True will trigger the plotting for the mid
    ufa[2]=True will trigger the plotting for the high
    ufa[3]=True will trigger the plotting for the full energy range

    The same is the case with cl_lc_type_trigger
    '''

    col_arr =['c','b','g','r']

    
    for N in range(len(obsid)): 
        obs=obsid[N]
          
        # Open the ufa event file to get the original GTI information
        os.chdir(obspath)
        filename = obs+'/xti/event_cl/ni'+obs+'_0mpu7_ufa.evt'
        hdul = fits.open(filename)
        gti_ufa = hdul['GTI'].data
        hdul.close()
        
        # get the GTI into a (x, 2) array to facilitate manipulation later on.
        n_gti_ufa = len(gti_ufa)
        gti_ufa = np.array([np.array(x) for x in gti_ufa])
        # defining a T0 as the start of the first ufa GTI.
        # The graph will use that T0 as the zero point of the time axis.
        T0 = gti_ufa[0,0]
        gti_ufa = gti_ufa - T0
        
        
        # get the duration of each GTI, and the gap in between.
        gti_ufa_duration, gti_ufa_gap = gti_info(gti_ufa)

        filename = obs+'/xti/event_cl/ni'+obs+'_0mpu7_cl'+clean_type+'.evt'
        hdul = fits.open(filename)
        gti_clean = hdul['GTI'].data
        hdul.close()

        # get the GTI into a (x, 2) array to facilitate manipulation later on.
        n_gti_clean = len(gti_clean)
        gti_clean = np.array([np.array(x) for x in gti_clean])
        gti_clean = gti_clean - T0 # The zero point is the start of the UFA ****
        gti_clean_duration, gti_clean_gap = gti_info(gti_clean)
        

        #Get the mkf data
        filename = obs+'/auxil/ni'+obs+mkf_ext
        hdul = fits.open(filename)
        head_mk = hdul[1].header
        data = hdul[1].data
        hdul.close()
        # Make t=0 the start of the first UFA GTI.
        MKtime = data['TIME']-T0

        # This piece of code merge together some GTIs that have very small gap
        # so that we can make a better use of a page size.
        # If the gap is rather large, then the next GTI will be printed on a new page.
        gti_plot = np.array( (1,2) )
        gti_plot.shape = (1,2)
        gti_plot[0,0] = gti_ufa[0,0] # set the start of the plit gti to that of the first GTI.
        k=0 # iteration variable for changes
        for i in range(0,n_gti_ufa):
            if gti_ufa_gap[i] > merge_gap: # if the gap with the last seg is more than 5 minutes
                gti_plot[k,1] = gti_ufa[i-1,1] #Set the stop of the k GTI to that of the previous segment
                gti_plot = np.vstack( [gti_plot, np.array( [gti_ufa[i,0], 0] ) ] ) # append a new plot_GTI
                                # and set the start time to the current GTI.
                k = k+1
        gti_plot[-1,1]=gti_ufa[-1,1] # set the stop of last gti
        n_gti_plot = k+1
        
        gti_plot_duration, gti_plot_gap = gti_info(gti_plot)
        
        num_lc_panel=0

        cl_ax_num   =2
        if any(ufa_lc_type_trigger)==True:
            cl_ax_num=3
            num_lc_panel+=1
        if any(cl_lc_type_trigger)==True:
            num_lc_panel+=1

        #Open the lightcurves
        max_val_ufa =1.
        if ufa_lc_type_trigger[0]==True:
            os.chdir(obspath+'/ufa/low')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_ufa_low, LC_ufa_low    =load_curve(lc_name, T0)
            max_val_ufa =max(max_val_ufa,max(LC_ufa_low['RATE']))

        if ufa_lc_type_trigger[1]==True:
            os.chdir(obspath+'/ufa/mid')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_ufa_mid, LC_ufa_mid    =load_curve(lc_name, T0)
            max_val_ufa =max(max_val_ufa,max(LC_ufa_mid['RATE']))


        if ufa_lc_type_trigger[2]==True:
            os.chdir(obspath+'/ufa/high')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_ufa_hi, LC_ufa_hi    =load_curve(lc_name, T0) 
            max_val_ufa =max(max_val_ufa,max(LC_ufa_hi['RATE']))



        if ufa_lc_type_trigger[3]==True:
            os.chdir(obspath+'/ufa/full')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_ufa_full, LC_ufa_full    =load_curve(lc_name, T0)
            max_val_ufa =max(max_val_ufa,max(LC_ufa_full['RATE']))


        max_val_cl  =1.0
        if cl_lc_type_trigger[0]==True:
            os.chdir(obspath+'/cl/low')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_cl_low, LC_cl_low    =load_curve(lc_name, T0)
            max_val_cl =max(max_val_cl,max(LC_cl_low['RATE']))
        
        if cl_lc_type_trigger[1]==True:
            os.chdir(obspath+'/cl/mid')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_cl_mid, LC_cl_mid    =load_curve(lc_name, T0)
            max_val_cl =max(max_val_cl,max(LC_cl_mid['RATE']))
        

        if cl_lc_type_trigger[2]==True:
            os.chdir(obspath+'/cl/high')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_cl_hi, LC_cl_hi    =load_curve(lc_name, T0) 
            max_val_cl =max(max_val_cl,max(LC_cl_hi['RATE']))
        

        if cl_lc_type_trigger[3]==True:
            os.chdir(obspath+'/cl/full')
            lc_name     ='ni'+obs+clean_type+'.lc'
            time_cl_full, LC_cl_full    =load_curve(lc_name, T0)
            max_val_cl =max(max_val_cl,max(LC_cl_full['RATE']))
        

        low_en_chan1,mid_en_chan1,high_en_chan1=chan1_arr[0],chan1_arr[1],chan1_arr[2]
        low_en_chan2,mid_en_chan2,high_en_chan2=chan2_arr[0],chan2_arr[1],chan2_arr[2]

        low_label   =str(low_en_chan1/100)+'-'+str(low_en_chan2/100)+' keV'
        mid_label   =str(mid_en_chan1/100)+'-'+str(mid_en_chan2/100)+' keV'
        high_label  =str(high_en_chan1/100)+'-'+str(high_en_chan2/100)+' keV'
        full_label  ='all keV'
        
        with PdfPages('{}/ni{}_mk{}.pdf'.format(PDF_path, obs,clean_type)) as pdf:
            for i in range(0,n_gti_plot):
            #for i in range(1,2):
                tmin = gti_plot[i,0]
                tmax = gti_plot[i,1]

                n = np.where( np.logical_and( MKtime >= tmin, MKtime <= tmax ) )
                fig, ax = plt.subplots(4,1, figsize=(8,10))
                fig, ax = plt.subplots(2+num_lc_panel,1, figsize=(8,10))
                #-----------------------
                # Overonly
                ax[0].scatter(MKtime[n], data['FPM_OVERONLY_COUNT'][n], s=3, zorder=500)
                ax[0].axhline(y=1.0, c='k', ls='--', label='abs cutoff')
                ax[0].plot( MKtime[n], 1.52*data['COR_SAX'][n]**(-0.633), c='orchid', lw=2, zorder=1000, label='overonly_exp' )
                ax[0].set_ylim(0,3)
                ax[0].set_ylabel('overonly')
                ax[0].set_xlabel('Time since T0 (s)')
                ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[0].legend(loc=0)
                #-----------------------
                # Underonly
                ax[1].scatter(MKtime[n], data['FPM_UNDERONLY_COUNT'][n], s=3, zorder=500)
                ax[1].set_ylim(0,300)
                ax[1].set_ylabel('underonly')
                ax[1].axhline(y=200, c='k', ls='--', label='def cutoff')
                ax[1].set_xlabel('Time since T0 (s)')
                #ax[0].set_title('{}, segment {}'.format(obs, i))
                ax[1].legend(loc=0)
                #-----------------------
                # The light curves from the ufa
                
                if ufa_lc_type_trigger[0]==True:
                    n2 = np.where( np.logical_and(  time_ufa_low >= tmin, time_ufa_low <= tmax ) )
                    ax[2].errorbar(time_ufa_low[n2], LC_ufa_low['RATE'][n2], yerr=LC_ufa_low['ERROR'][n2], fmt='.', ms=2, color=col_arr[0], label=low_label, zorder=500 )
                if ufa_lc_type_trigger[1]==True:
                    n2 = np.where( np.logical_and(  time_ufa_mid >= tmin, time_ufa_mid <= tmax ) )
                    ax[2].errorbar(time_ufa_mid[n2], LC_ufa_mid['RATE'][n2], yerr=LC_ufa_mid['ERROR'][n2], fmt='.', ms=2, color=col_arr[1], label=mid_label, zorder=500 )
                if ufa_lc_type_trigger[2]==True:
                    n2 = np.where( np.logical_and(  time_ufa_hi >= tmin, time_ufa_hi <= tmax ) )
                    ax[2].errorbar(time_ufa_hi[n2], LC_ufa_hi['RATE'][n2], yerr=LC_ufa_hi['ERROR'][n2], fmt='.', ms=2, color=col_arr[2], label=high_label, zorder=500 )
                if ufa_lc_type_trigger[3]==True:
                    n2 = np.where( np.logical_and(  time_ufa_full >= tmin, time_ufa_full <= tmax ) )
                    ax[2].errorbar(time_ufa_full[n2], LC_ufa_full['RATE'][n2], yerr=LC_ufa_full['ERROR'][n2], fmt='.', ms=2, color=col_arr[3], label=full_label, zorder=500 )
                if cl_ax_num==3:
                    ax[2].plot(MKtime[n], data['SAA'][n], c='k', label='SAA (1:yes, 0:no)' )
                    ax[2].set_ylabel('UFA Count Rate')
                    ax[2].set_ylim(0,max_val_ufa)
                    #ax[2].axhline(y=1.0, ls='--', c='0.5')
                    #ax[2].axhline(y=2.0, ls='--', c='0.5')
                    ax[2].legend(loc=0)
                #-----------------------
     


                # The light curves from the current filtering
                if cl_lc_type_trigger[0]==True:
                    n2 = np.where( np.logical_and(  time_cl_low >= tmin, time_cl_low <= tmax ) )
                    ax[cl_ax_num].errorbar(time_cl_low[n2], LC_cl_low['RATE'][n2], yerr=LC_cl_low['ERROR'][n2], fmt='.', ms=2, color=col_arr[0], label=low_label, zorder=500 )
                if cl_lc_type_trigger[1]==True:
                    n2 = np.where( np.logical_and(  time_cl_mid >= tmin, time_cl_mid <= tmax ) )
                    ax[cl_ax_num].errorbar(time_cl_mid[n2], LC_cl_mid['RATE'][n2], yerr=LC_cl_mid['ERROR'][n2], fmt='.', ms=2, color=col_arr[1], label=mid_label, zorder=500 )
                if cl_lc_type_trigger[2]==True:
                    n2 = np.where( np.logical_and(  time_cl_hi >= tmin, time_cl_hi <= tmax ) )
                    ax[cl_ax_num].errorbar(time_cl_hi[n2], LC_cl_hi['RATE'][n2], yerr=LC_cl_hi['ERROR'][n2], fmt='.', ms=2, color=col_arr[2], label=high_label, zorder=500 )
                if cl_lc_type_trigger[3]==True:
                    n2 = np.where( np.logical_and(  time_cl_full >= tmin, time_cl_full <= tmax ) )
                    ax[cl_ax_num].errorbar(time_cl_full[n2], LC_cl_full['RATE'][n2], yerr=LC_cl_full['ERROR'][n2], fmt='.', ms=2, color=col_arr[3], label=full_label, zorder=500 )
                if any(cl_lc_type_trigger)==True:
                    ax[cl_ax_num].set_ylabel('Filtered Count Rate')
                    ax[cl_ax_num].set_ylim(0,max_val_cl)
                    #ax[cl_ax_num].axhline(y=1.0, ls='--', c='0.5')
                    #ax[cl_ax_num].axhline(y=2.0, ls='--', c='0.5')
                    ax[cl_ax_num].legend(loc=0)
                #-----------------------
                # Overplotting the clean GTIs
                for j in range(0,n_gti_clean):
                    ttmin = gti_clean[j,0]
                    ttmax = gti_clean[j,1]
                    if np.logical_and( ttmin >= tmin, ttmin <= tmax ):
                        ax[0].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        ax[1].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        if any(ufa_lc_type_trigger)==True:
                            ax[2].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)
                        if any(cl_lc_type_trigger)==True:
                            ax[cl_ax_num].axvspan(ttmin, ttmax, alpha=0.1, color='green', zorder=100)

                for item in ax:
                    item.set_xlim(tmin, tmax)
                plt.tight_layout()
                pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
                plt.close()
                
    return

def nimaketime_expr(infile,outfile,include=False):
    # include determines whether the time intervals given in infile are to be included or excluded in the gtifile
   
    data=np.loadtxt(infile)
    if len(data)==0:
        return False,''
    else:
        if len(np.shape(data[0]))>0:
            t_start,t_end=data[:,0],data[:,1]
        else:
            t_start,t_end   =np.array([data[0]]),np.array([data[1]])

        f=open(outfile,'w')
        n   =len(t_start)
        x=''
        for i in range(n):
            if i==0:
                if include==True:
                    f.write('(time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
                    x+='(time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'
                else:
                    f.write('!((time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
                    x+='!((time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'
            else:
                f.write(' || (time>='+str(t_start[i])+' && time<='+str(t_end[i])+')')
                x+=' || (time>='+str(t_start[i])+' && time<='+str(t_end[i])+')'

        if include==False:
            f.write(')')
            x+=')'
            
        f.close()
        return True,x


def nimaketime_command_create(folder,obsid,include=False, run_f='run.sh'):

    '''
    This function creates the command for the nimaketime task.
    The commands are written in run.sh inside the folder (see below) directory.
    After running bash run.sh from the terminal, the gti file with name niObsID.gti will be created inside each of the ObsID directory.
    Inside folder/ObsID/, there should be text file with name obsid_time.txt, 
    the obsid_time.txt is expected to have a list of time intervals which are to be included (include=True) or excluded (include=False) while running nicerl2.

    In addition, this function also writes down the expression used in nimaketime inside the ObsID directory under the name of ObsID_nimaketime_expr.txt.
   
    :param folder: folder is the directory containing directories with names given in filename. 
    :param filename: ('ObsID.dat') name of the file containing the list of ObsID
    :param include: whether to include the gtis, or exclude them
    :param run_f: ('run.sh') name of the file containing all of the run commands (created in `folder`)

    '''

    os.chdir(folder)
    run = open(run_f, 'w')

    for obs in obsid:
        infile  =obs+'/'+obs+'_time.txt'
        outfile =obs+'/'+obs+'_nimaketime_expr.txt'

        truth_val,expr    =nimaketime_expr(infile,outfile,include)
        if truth_val==True:
            run.write('nimaketime infile=\'{}/{}/auxil/ni{}.mkf\' outfile=\'{}/{}/ni{}.gti\' expr=\'{}\' clobber=YES\n'.format(folder,obs,obs,folder,obs,obs,expr))
        
                
    run.close()


def get_bad_gti(bad_time_arr,f_mode,obsid,outdir,time_res=1):
    f   =open(outdir+str(obsid)+'_time.txt',f_mode)
    N   =len(bad_time_arr)
    
    diff_time   =np.diff(bad_time_arr)
    
    pos=np.where(diff_time>time_res)[0]
    if time_res>1:
        time_res=1*time_res                                     #@@@@@@@@@@@@@@@@@@@@@@@CAUTION**********************************
    if (len(pos)>0):
        if f_mode=='w':
            f.write(str(bad_time_arr[0]-0.5*time_res)+'\t'+str(bad_time_arr[pos[0]]+0.5*time_res))
        else:    
            f.write('\n'+str(bad_time_arr[0]-0.5*time_res)+'\t'+str(bad_time_arr[pos[0]]+0.5*time_res))
        for i in range(1,len(pos)):
            f.write('\n'+str(bad_time_arr[pos[i-1]+1]-0.5*time_res)+'\t'+str(bad_time_arr[pos[i]]+0.5*time_res))
        f.write('\n'+str(bad_time_arr[pos[-1]+1]-0.5*time_res)+'\t'+str(bad_time_arr[N-1]+0.5*time_res))
    else:
        if N>0:
            f.write('\n'+str(min(bad_time_arr)-0.5*time_res)+'\t'+str(max(bad_time_arr)+0.5*time_res))
    f.close()

def identify_bad_gti_lc(obspath,obsid,lc_path,clean_type,T_ref=0,time_res=50,fmode='w'):
    #print(time_res) 
    time,rate,rate_err  =[],[],[]
    os.chdir(lc_path)
    for obs in obsid:        
        filename    ='ni'+obs+clean_type+'.lc'
        hdul = fits.open(filename)
        LC = hdul[1].data
        LC_head=hdul[1].header
        hdul.close()
 
        lc_type=LC_head['TTYPE2']
        if lc_type=='COUNTS':
            x   =LC[lc_type]/LC_head['EXPOSURE']
        if lc_type=='RATE':
            x   =LC[lc_type]
    
        T0=LC_head['TIMEZERO']
        time_LC =LC['TIME']+T0-T_ref
        
        time.append(list(time_LC))
        rate.append(list(x))
        rate_err.append(list(LC['ERROR']))
        
    my_time=np.array(reduce(operator.concat, time))
    my_rate    =np.array(reduce(operator.concat,rate))
    my_rate_err=np.array(reduce(operator.concat,rate_err))
    pos =np.where(my_rate>0.3)[0]
    median_rate=np.median(np.delete(my_rate,pos)) 
    MAD =median_abs_deviation(np.delete(my_rate,pos),nan_policy='omit')
    
    factor  =1#1.48
    new_pos =np.where(((my_rate-my_rate_err>median_rate+(2/factor)*MAD) | ((my_rate+my_rate_err>median_rate+(3/factor)*MAD))))[0]
    tol =1.2
    for i in range(len(obsid)):
        time_arr,rate_arr,rate_err_arr=np.array(time[i]),np.array(rate[i]),np.array(rate_err[i])
        pos =np.where(((rate_arr-rate_err_arr>median_rate+(2/factor)*MAD) | ((rate_arr+rate_err_arr>median_rate+(3/factor)*MAD))))[0]
        bad_time=time_arr[pos]
        
        del_pos_index   =[]
        for j in range(len(pos)):
            my_pos=pos[j]
            if my_pos>0 and my_pos<len(time_arr)-1:
                if np.any(bad_time==time_arr[my_pos-1])==False and np.any(bad_time==time_arr[my_pos+1])==False:
                    rate_val    =rate_arr[my_pos]
                    rate_err_val    =rate_err_arr[my_pos]
                    if rate_val-rate_err_val<median_rate+tol*(2/factor)*MAD and rate_val+rate_err_val<median_rate+tol*(3/factor)*MAD:
                        del_pos_index.append(j)
                     
        if len(del_pos_index)>0:
            pos =np.delete(pos,np.array(del_pos_index))
           
        os.chdir(obspath)
        
        get_bad_gti(time_arr[pos],fmode,obsid[i],outdir=obspath+'/'+str(obsid[i])+'/',time_res=time_res)



def get_arf_rmf_3c50_tot_spec(obsid,obs_dir='.',pha_dir='.',clean_type='_default'):
    #pha_dir should be inside the obs_dir
    

    os.chdir(obs_dir)   

    for i in range(len(obsid)):
        print(obsid[i])
        obs =obsid[i]
        mkf_file    =obs+'/auxil/ni'+obs+'.mkf'
        hdul    =fits.open(mkf_file)
        RA,dec  =str(hdul[1].header["RA_NOM"]) , str(hdul[1].header["DEC_NOM"])
        hdul.close()
        event_file  =obs+'/xti/event_cl/ni'+obs+'_0mpu7_cl'+clean_type+'.evt'
        spec_file   =pha_dir+'/nibackgen3C50_'+obs+clean_type+'_tot.pi'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'      
        outwtfile   ='nibackgen3C50_'+obs+'_tot_wt.lis'
        print(os.getcwd())
        a=os.system('nicerarf '+spec_file+' '+RA+' '+dec+' '+mkf_file+' '+event_file+' '+arf_file+' outwtfile='+outwtfile+' updatepha=YES clobber=YES')
        if a!=0:
            os.system('nicerarf '+spec_file+' '+RA+' '+dec+' '+mkf_file+' '+mkf_file+' '+arf_file+' outwtfile='+outwtfile+' updatepha=YES clobber=YES')
        os.system('nicerrmf '+spec_file+' '+mkf_file+' '+rmf_file+' detlist=@'+outwtfile+' clobber=YES')
        os.system('mv '+arf_file+' '+obs_dir+'/'+pha_dir)
        os.system('mv '+rmf_file+' '+obs_dir+'/'+pha_dir)
        os.system('mv '+outwtfile+' '+obs_dir+'/'+pha_dir)



def get_bkg_3c50(obsid,obs_dir='.',bkg_dir='/home/dbarnali/heasoft/3c50_bgmodel/bg_models_3C50',out_dir='.',clean_type='_flagged_v2'):
    #out_dir must end with /
    ini_path    =os.getcwd()
    

    for obs in obsid:
        os.chdir(obs_dir+'/'+obs+'/xti/event_cl')
        ufa_file    ='ni'+obs+'_0mpu7_ufa.evt'
        cl_file     ='ni'+obs+'_0mpu7_cl'+clean_type+'.evt'
        command_str ='nibackgen3C50 rootdir=NONE obsid=NONE bkgidxdir='+bkg_dir+' bkglibdir='+bkg_dir+' gainepoch=2019 calevtdir=NONE ufafile='+ufa_file+' clfile='+cl_file+' clobber=YES'
        os.system(command_str)
        os.system('rename s/_/_{}{}_/ nibackgen3C50*.pi'.format(obs,clean_type))
        os.system('mv nibackgen3C50*.pi '+out_dir)

    os.chdir(ini_path)


def perform_optimal_binning(obsid,clean_type,spec_dir,min_count=15,min_chan_num=26,max_chan_num=1000):
    os.chdir(spec_dir)
    for obs in obsid:
        in_spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot.pi'
        bkg_file        ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        out_spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        rmf_file        ='nibackgen3C50_'+obs+'_tot.rmf'
        #rmf_file        =in_spec_file.split('.')[0]+'.rmf'
        command_str ='ftgrouppha infile='+in_spec_file+' backfile='+bkg_file+' outfile='+ out_spec_file+' grouptype=optmin groupscale='+str(min_count)+' minchannel='+str(min_chan_num)+' maxchannel='+str(max_chan_num)+' respfile='+rmf_file
        os.system(command_str)
    
    return

def perform_optimal_binning_for_bkg(obsid,clean_type,spec_dir,min_count=15,min_chan_num=26,max_chan_num=1000):
    os.chdir(spec_dir)
    for obs in obsid:
        bkg_file        ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        in_spec_file    =bkg_file
        out_spec_file   ='nibackgen3C50_'+obs+clean_type+'_bkg_grp'+str(min_count)+'.pi'
        rmf_file        ='nibackgen3C50_'+obs+'_tot.rmf'
        #rmf_file        =in_spec_file.split('.')[0]+'.rmf'
        command_str ='ftgrouppha infile='+in_spec_file+' outfile='+ out_spec_file+' grouptype=optmin groupscale='+str(min_count)+' minchannel='+str(min_chan_num)+' maxchannel='+str(max_chan_num)+' respfile='+rmf_file
        os.system(command_str)
    
    return


def merge_events(obsid,out_ID,clean_type,obs_dir,filter_type_arr=['mkf']):
    os.chdir(obs_dir)
    if os.path.isdir(out_ID)==False:
        os.mkdir(out_ID)
        os.chdir(out_ID)
        os.mkdir('xti')
        os.mkdir('auxil')
        os.chdir('xti')
        os.mkdir('event_cl')
        out_dir =obs_dir+'/'+out_ID+'/xti/event_cl'
        out_mkfdir=obs_dir+'/'+out_ID+'/auxil'
        os.chdir(obs_dir)
    else:
        if os.path.isdir(out_ID+'/xti')==False:
            os.chdir(out_ID)
            os.mkdir('xti')
            os.chdir('xti')
            os.mkdir('event_cl')
            out_dir =obs_dir+'/'+out_ID+'/xti/event_cl'
            os.chdir(obs_dir)
        else:
            if os.path.isdir(out_ID+'xti/event_cl')==False:
                os.chdir(out_ID+'/xti')
                os.mkdir('event_cl')
                out_dir =obs_dir+'/'+out_ID+'/xti/event_cl'
                os.chdir(obs_dir)
            else:
                out_dir =obs_dir+'/'+out_ID+'/xti/event_cl'
        
        if os.path.isdir(out_ID+'/auxil')==False:
            os.chdir(out_ID)
            os.mkdir('auxil')
            out_mkfdir=obs_dir+'/'+out_ID+'/auxil'
            os.chdir(obs_dir)
        else:
            out_mkfdir=obs_dir+'/'+out_ID+'/auxil'
    #########For ufa event files ###########
    event_type  ='ufa'
    out_evt ='ni'+out_ID+'_0mpu7_'+event_type+'.evt'
    event_listfile  ='evtfiles.lis'

    f   =open(event_listfile,'w')
    for obs in obsid:
        f.write('{}/{}/xti/event_cl/ni{}_0mpu7_{}.evt\n'.format(obs_dir,obs,obs,event_type))
    f.close()
    command_str ='nimpumerge infiles=@evtfiles.lis outfile='+out_dir+'/'+out_evt+' mpulist=7 clobber=yes'
    os.system(command_str)

    #########For cl event files ###########
    event_type  ='cl'+clean_type
    out_evt     ='ni'+out_ID+'_0mpu7_'+event_type+'.evt'
    event_listfile  ='evtfiles.lis'

    f   =open(event_listfile,'w')
    for obs in obsid:
        f.write('{}/{}/xti/event_cl/ni{}_0mpu7_{}.evt\n'.format(obs_dir,obs,obs,event_type))
    f.close()
    command_str ='nimpumerge infiles=@evtfiles.lis outfile='+out_dir+'/'+out_evt+' mpulist=7, clobber=yes'
    os.system(command_str)
    
    ####### For filter files ###########
    for filter_type in filter_type_arr:
        out_evt ='ni'+out_ID+'.'+filter_type 
        event_listfile  ='mkffiles.lis'
        f   =open(event_listfile,'w')
        for obs in obsid:
            f.write('{}/{}/auxil/ni{}.{}\n'.format(obs_dir,obs,obs,filter_type))
        f.close()
        command_str ='nimkfmerge infiles=@mkffiles.lis outfile='+out_mkfdir+'/'+out_evt+' clobber=yes'
        os.system(command_str)


'''

def underonly_statistics(underonly_arr,max_count=15,savefig=False,figname=None,showfig=False):
    #plt.plot(underonly_arr,'o')
    #plt.plot(np.gradient(underonly_arr))
    #plt.axhline(y=0.2)
    #plt.axhline(y=-0.2)
    #plt.show()
    
    N   =len(underonly_arr)
    pos0    =np.arange(0,N,1)
    
    cont=1
    my_arr  =underonly_arr
    my_pos  =pos0
    bad_pos=[]

    pos =np.where(abs(my_arr)>max_count)[0]
    if (len(pos)>0):
        my_arr =np.delete(my_arr,pos)
        bad_pos.append(list(my_pos[pos]))
        my_pos  =np.delete(my_pos,pos)

    if savefig==True or showfig==True:
        plt.figure()
        plt.plot(my_pos,my_arr,'^')   
    while cont==1:
        med =np.nanmedian(my_arr)
        MAD =median_abs_deviation(my_arr,nan_policy='omit')
        dev =abs(my_arr-med)
        pos =np.where(dev>2.8*MAD)[0]
        if len(pos)==0:
            cont=0
        else:
            #print(my_arr)
            my_arr =np.delete(my_arr,pos)
            bad_pos.append(list(my_pos[pos]))
            my_pos  =np.delete(my_pos,pos)
        
    #plt.figure()
    if savefig==True or showfig==True:
        plt.plot(my_pos,my_arr,'d',alpha=0.5,markersize=10,markerfacecolor='w',markeredgecolor='k',markeredgewidth=2)
        plt.plot(pos0,underonly_arr,'r-',alpha=1,lw=2)
        plt.ylabel('Underonly count')
        plt.ylim(0,3*max_count)
        plt.tight_layout()
        if showfig==True:
            plt.show()
        else:    
            plt.savefig(figname)
            plt.close()
    
    if len(bad_pos)>0:
        bad_pos=reduce(operator.concat, bad_pos)
        sort_pos    =np.argsort(bad_pos)
        bad_pos  =np.array(bad_pos)[sort_pos]       
    return  bad_pos   
    #print(med,MAD,len(underonly_arr),len(r_mean))

def underonly_grad_statistics(underonly_arr,max_count=9,grad_max_count=1,savefig=False,figname=None,showfig=False,n_smooth=3):
    maxpos =np.where(underonly_arr>max_count)[0]
    N   =len(underonly_arr)
    pos0    =np.arange(0,N,1)
    #n_smooth    =3
    grad_under_arr=np.gradient(underonly_arr)
    sm_grad_under_arr   =np.convolve(grad_under_arr,np.ones(n_smooth)/n_smooth,mode='valid')[::n_smooth]
    
    bad_pos =underonly_statistics(sm_grad_under_arr,grad_max_count,False,None,False)
    if n_smooth>1:
        new_bad_pos=[]
        for i in range(len(bad_pos)):
            val =bad_pos[i]
            new_bad_pos.append(list(np.arange(n_smooth*val,n_smooth*(val+1),1)))
        if len(new_bad_pos)>0:    
            new_bad_pos=reduce(operator.concat, new_bad_pos)
        new_bad_pos =np.array(new_bad_pos)    
    else:
        new_bad_pos =bad_pos
        
    if len(maxpos)>0:
        new_bad_pos=np.concatenate((new_bad_pos,maxpos))
    new_bad_pos =new_bad_pos[np.argsort(new_bad_pos)]    
    
    if savefig==True or showfig==True:
        plt.figure()
        if len(new_bad_pos)>0:
            plt.plot(pos0[new_bad_pos],underonly_arr[new_bad_pos],'o',alpha=0.5,markersize=10,markerfacecolor='w',markeredgecolor='k',markeredgewidth=2)
        plt.plot(pos0,underonly_arr,'r-',alpha=1,lw=2)
        plt.ylabel('Underonly count')
        #plt.ylim(-3*max_count,3*max_count)
        plt.tight_layout()
        if showfig==True:
            plt.show()
        else:    
            plt.savefig(figname)
            plt.close()
    return  new_bad_pos   



def identify_bad_gti_underonly(obspath=obs_dir,obsid_file=obsid_f,mkf_ext='mkf',max_count=15,grad_max_count=1):
    os.chdir(obspath)
    obsid = np.loadtxt(filename, dtype=str)
    if len(np.shape(obsid))==0: #if there is just one observation ID
        obsid   =np.array([obsid])

    for N in range(len(obsid)): 
        obs=obsid[N]
        filename = obs+'/xti/event_cl/ni'+obs+'_0mpu7_ufa.evt'
        hdul = fits.open(filename)
        gti_ufa = hdul['GTI'].data
        hdul.close()
        
        # get the GTI into a (x, 2) array to facilitate manipulation later on.
        n_gti_ufa = len(gti_ufa)
        gti_ufa = np.array([np.array(x) for x in gti_ufa])
        # defining a T0 as the start of the first ufa GTI.
        # The graph will use that T0 as the zero point of the time axis.
        T0 = gti_ufa[0,0]
        gti_ufa = gti_ufa - T0
        
        
        # get the duration of each GTI, and the gap in between.
        gti_ufa_duration, gti_ufa_gap = gti_info(gti_ufa)


        mkf_file    =obs+'/auxil/ni'+obs+mkf_ext
        hdul = fits.open(filename)
        head_mk = hdul[1].header
        data = hdul[1].data
        hdul.close()
        # Make t=0 the start of the first UFA GTI.
        MKtime = data['TIME']-T0

        gti_plot = np.array( (1,2) )
        gti_plot.shape = (1,2)
        gti_plot[0,0] = gti_ufa[0,0] # set the start of the plit gti to that of the first GTI.
        k=0 # iteration variable for changes
        for i in range(0,n_gti_ufa):
            if gti_ufa_gap[i] > merge_gap: # if the gap with the last seg is more than 5 minutes
                gti_plot[k,1] = gti_ufa[i-1,1] #Set the stop of the k GTI to that of the previous segment
                gti_plot = np.vstack( [gti_plot, np.array( [gti_ufa[i,0], 0] ) ] ) # append a new plot_GTI
                                # and set the start time to the current GTI.
                k = k+1
        gti_plot[-1,1]=gti_ufa[-1,1] # set the stop of last gti
        n_gti_plot = k+1
        
        gti_plot_duration, gti_plot_gap = gti_info(gti_plot)
        for i in range(0,n_gti_plot):
            #for i in range(1,2):
            tmin = gti_plot[i,0]
            tmax = gti_plot[i,1]

            n = np.where( np.logical_and( MKtime >= tmin, MKtime <= tmax ) )
            time_arr=MKtime[n]+T0
            obsid_arr.append(obs)
            # Underonly
            underonly_arr=data['FPM_UNDERONLY_COUNT'][n]
            my_bad_pos=underonly_grad_statistics(underonly_arr,max_count, grad_max_count)
            if len(my_bad_pos)>0:
                if i==0:
                    f_mode='w'
                else:
                    f_mode='a'
                get_bad_gti(time_arr[my_bad_pos],f_mode,obs,outdir=str(obs)+'/')





'''
