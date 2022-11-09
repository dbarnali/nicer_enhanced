import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from xspec import *
import os
from astropy.io import fits
import time
import corner


def spectral_analysis_apec(obsid,clean_type,spec_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE=0.4,endE=10.0,bkg_type='3c50',chain_length=200000,burn_length=2000,input_wt='standard',out_dir='.',min_count=15,save_fig=False,write_norm=False):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_varying_nH_mcmc_with_nH_bound_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    

    for N in range(len(obsid)):
        obs =obsid[N]
        out_name    =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        cont=os.path.isfile(out_name)
        if cont==True:
            continue
        else:
            os.chdir(spec_dir)
            bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
            rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
            arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
            s           =Spectrum(spec_file)
            s.background    =bkg_file
            s.response      =rmf_file
            s.response.arf  =arf_file
            AllData.ignore("bad")
            s.ignore(f"**-{startE} {endE}-**")
            ############################# Part to edit
            model_str='tbabs('
            for i in range(len(kT_arr)):
                if i==0:
                    model_str+='apec'
                else:
                    model_str+='+apec'
            model_str+=')'        
            m1=Model(model_str)
            m1.TBabs.nH =nH0
            m1.TBabs.nH.frozen  =freeze_nH
            if freeze_nH==False:
                m1.TBabs.nH =[nH0,0.01,0.0,0.0,10*nH0,10*nH0]
            for i in range(len(kT_arr)):
                m1(2+4*i).values=kT_arr[i]
                m1(2+4*i).frozen=freeze_kT_arr[i]
                if freeze_kT_arr[i]==False:
                    m1(2+4*i).values=[kT_arr[i],0.01,kT_arr[i]/3.,kT_arr[i]/3.,3*kT_arr[i],3*kT_arr[i]] # HARD-CODED
            #########################
            Fit.query = "yes"
            Fit.nIterations=100
            #Fit.statMethod = "cstat"
            Fit.weight  =input_wt   
            Fit.renorm()
            Fit.perform()
         
            AllChains.defBurn = burn_length
            AllChains.defLength = chain_length
            print('I am going to run chain')
    
            c1 = Chain(out_name,rand=False)
            
            if save_fig==True:
                data=fits.getdata(out_name,1)
                head=fits.getheader(out_name,1)
                data1=np.array([np.array(x) for x in data])
                median_vals =np.median(data1,axis=0)
                labels=[head['TTYPE'+str(i)] for i in range(1,len(data1[0])+1)]
                corner.corner(data1[:,:], labels=labels[:],truths=median_vals[:],quantiles=[0.16, 0.5, 0.84],
                       show_titles=False, title_kwargs={"fontsize": 12})
                plt.savefig(out_dir+'/'+obs+'_'+model_type+out_mcmc_name)
                plt.close()

                if freeze_nH==False:
                    start_index=1
                    m1.TBabs.nH     =median_vals[0]
                else:
                    start_index=0

                j=0
                for i in range(len(kT_arr)):
                    if freeze_kT_arr[i]==True:
                        m1(5+4*i).values=median_vals[start_index+j]
                    else:
                        m1(2+4*i).values=median_vals[start_index+j]
                        j+=1
                        m1(5+4*i).values=median_vals[start_index+j]
                    j+=1    

                 ###############
                Plot.xAxis='keV'
                Plot("data resid")
                energies = Plot.x()
                edeltas = Plot.xErr()
                rates = Plot.y(1,1)
                errors = Plot.yErr(1,1)
                foldedmodel = Plot.model()
                dataLabels = Plot.labels(1)
                residLabels = Plot.labels(2)
                # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
                # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
                nE = len(energies)
                stepenergies = list()
                for j in range(nE):
                    stepenergies.append(energies[j] - edeltas[j])
                stepenergies.append(energies[-1]+edeltas[-1])
                foldedmodel.append(foldedmodel[-1])
                resid = Plot.y(1,2)
                residerr = Plot.yErr(1,2)
                #print('The current ObsID is',obs)
                plt.subplot(211)
                plt.xscale('log')
                #plt.yscale('log')
                plt.ylabel(dataLabels[1])
                plt.title(dataLabels[2])
                plt.errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.')
                plt.step(stepenergies,foldedmodel,where='post')
                plt.subplot(212)
                plt.xscale('log')
                plt.xlabel(residLabels[0])
                plt.ylabel(residLabels[1])
                plt.errorbar(energies,resid,xerr=edeltas,yerr=residerr,fmt='.')
                plt.hlines(0.0,stepenergies[0],stepenergies[-1],linestyles='dashed')
                plt.savefig(out_dir+'/'+obs+'_'+model_type+out_spec_name)
                plt.close()
            if write_norm==True:
                f=open(out_dir+'/'+obs+'_'+model_type+norm_file_name,'w')
                f.write('### '+obs+'\n\n')
                
                if save_fig==False:
                    data=fits.getdata(out_name,1)
                    head=fits.getheader(out_name,1)
                    data1=np.array([np.array(x) for x in data])
                    median_vals =np.median(data1,axis=0)
                u_vals_1sig  =np.percentile(data1,84.0,axis=0)
                d_vals_1sig  =np.percentile(data1,16.0,axis=0)
                u_vals_2sig  =np.percentile(data1,97.5,axis=0)
                d_vals_2sig  =np.percentile(data1,2.5,axis=0)  
                u_vals_3sig  =np.percentile(data1,99.85,axis=0)
                d_vals_3sig  =np.percentile(data1,0.15,axis=0)
                my_str=''
                if freeze_nH==False:
                    my_str+='## nH='+str(median_vals[0])+'\n'
                    my_str+='## 1sigma nH interval='+str(u_vals_1sig[0])+'-'+str(d_vals_1sig[0])+'\n'
                    my_str+='## 2sigma nH interval='+str(u_vals_2sig[0])+'-'+str(d_vals_2sig[0])+'\n'
                    my_str+='## 3sigma nH interval='+str(u_vals_3sig[0])+'-'+str(d_vals_3sig[0])+'\n'
                    start_index=1
                else:
                    start_index=0
                my_str+='## Chi-sq='+str(median_vals[-1])+'\n'
                my_str+='## 1sigma Chi-sq interval='+str(u_vals_1sig[-1])+'-'+str(d_vals_1sig[-1])+'\n'
                my_str+='## 2sigma Chi-sq interval='+str(u_vals_2sig[-1])+'-'+str(d_vals_2sig[-1])+'\n'
                my_str+='## 3sigma chi-sq interval='+str(u_vals_3sig[-1])+'-'+str(d_vals_3sig[-1])+'\n'
                    
                my_str+='\n#kT\tkT_2u\tkT_2d\tNorm\t\tNorm_1u\t\tNorm_1d\t\tNorm_2u\t\tNorm_2d\t\tNorm_3u\t\tNorm_3d\n'            
                f.write(my_str)
                j=0
    
                for i in range(len(kT_arr)):
                    if freeze_kT_arr[i]==True:
                        my_str=str(kT_arr[i])+'\t0\t0\t'
            
                        my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                    else:
                        my_str=str(median_vals[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'
                        j+=1
                        my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                    j+=1
                    f.write(my_str)

                f.close()

            AllChains.clear()
            AllData.clear()
            AllModels.clear()
    return


def calculate_flux_apec(ism_correct,obsid,clean_type,spec_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE0=0.4,endE0=10.0,bkg_type='3c50',input_wt='standard',out_dir='.',min_count=15):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'  
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'

    for N in range(len(obsid)):
        os.chdir(spec_dir)
        obs     =obsid[N]
        chain_name =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        obs0        =obs#obs.split('_')[0]
        data=fits.getdata(chain_name,1)   
        data=np.array([np.array(x) for x in data])
        ##############################
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            
        s           =Spectrum(spec_file)
        s.background    =bkg_file
        s.response      =rmf_file
        s.response.arf  =arf_file
        AllData.ignore("bad")
        s.ignore(f"**-{startE0} {endE0}-**")
        ############################# Part to edit
        model_str='tbabs('
        for i in range(len(kT_arr)):
            if i==0:
                model_str+='apec'
            else:
                model_str+='+apec'
        model_str+=')'        
        m1=Model(model_str)
        for j in range(0,len(data),1000):
            if freeze_nH==False:
                if ism_correct==False:
                    m1.TBabs.nH=data[j,0]
                else:
                    m1.TBabs.nH=0.0
                start_index =1
            else:
                if ism_correct==False:
                    m1.TBabs.nH =nH0
                else:
                    m1.TBabs.nH=0.0    
                start_index =0
            
            l=0
            for i in range(len(kT_arr)):
                if freeze_kT_arr[i]==True:
                    m1(2+4*i).values=kT_arr[i]
                    m1(5+4*i).values=data[j,start_index+l]
                else:
                    m1(2+4*i).values=data[j,start_index+l]
                    l+=1
                    m1(5+4*i).values=data[j,start_index+l]
                l+=1  
            m1.show()
            for E1, E2 in zip(startE_arr, endE_arr):
                if freeze_nH==True:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH_ISM_corrected.txt'
                else:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH_ISM_corrected.txt'
                if j==0:
                    f=open(out_file,'w')
                    f.write('#ObsID: '+obs+', startE='+str(startE0)+'keV, endE='+str(endE0)+'\n')
                    f.write('#targetFlux\n')
                else:
                    f=open(out_file,'a')
                
                AllModels.calcFlux(str(E1)+' '+str(E2)+' err')
                f.write(str(s.flux[0])+'\n')
                f.close()
                
        AllData.clear()
        AllModels.clear()
        AllChains.clear()


def spectral_analysis_apec_multiresponse(obsid,clean_type,spec_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE=0.4,endE=10.0,bkg_type='3c50',chain_length=200000,burn_length=2000,input_wt='standard',out_dir='.',min_count=15,save_fig=False,write_norm=False):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_varying_nH_mcmc_with_nH_bound_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    

    for N in range(len(obsid)):
        obs =obsid[N]
        out_name    =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        cont=os.path.isfile(out_name)

        os.chdir(spec_dir)
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg_grp'+str(min_count)+'.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        my_str      ='1:1 '+spec_file+' 2:2 '+bkg_file
        AllData(my_str)
        s1  =AllData(1)
        s2  =AllData(2)
        s1.background   ='none'
        s1.response     =rmf_file
        s1.response.arf =arf_file
        s2.response     =rmf_file
        s2.response.arf ='none'
        AllData.ignore("bad")
        s1.ignore(f"**-{startE} {endE}-**")
        s2.ignore(f"**-{startE} {endE}-**")
        ############################# Part to edit
        model_str='tbabs('
        for i in range(len(kT_arr)):
            if i==0:
                model_str+='apec'
            else:
                model_str+='+apec'
        model_str+=')'        
        Model(model_str)
        m1 = AllModels(1)
        m2 = AllModels(2)
        m1.TBabs.nH =nH0
        m1.TBabs.nH.frozen  =freeze_nH
        if freeze_nH==False:
            m1.TBabs.nH =[nH0,0.01,0.0,0.0,10*nH0,10*nH0]
        for i in range(len(kT_arr)):
            m1(2+4*i).values=kT_arr[i]
            m1(2+4*i).frozen=freeze_kT_arr[i]
            m2(5+4*i).values=0.0
            m2(5+4*i).frozen    =True

        s1.multiresponse[1] = rmf_file
        s2.multiresponse[1] = rmf_file
        m1_2 = Model("pow+gaussian",'myback',2)
        m1_2.gaussian.LineE=0.56#[0.56,0.01,0.54,0.55,0.57,0.58]
        m1_2.gaussian.Sigma=0.01#[0.01,0.005,0,0,0.1,0.1]
        m1_2.gaussian.LineE.frozen=True
        m1_2.gaussian.Sigma.frozen=True

        if cont==False:
            #########################
            Fit.query = "yes"
            Fit.nIterations=100
            #Fit.statMethod = "cstat"
            Fit.weight  =input_wt   
            Fit.renorm()
            Fit.perform()
         
            AllChains.defBurn = burn_length
            AllChains.defLength = chain_length
            print('I am going to run chain')
    
            c1 = Chain(out_name,rand=False)
            
        if save_fig==True:
            data=fits.getdata(out_name,1)
            head=fits.getheader(out_name,1)
            data1=np.array([np.array(x) for x in data])
            median_vals =np.median(data1,axis=0)
            labels=[head['TTYPE'+str(i)] for i in range(1,len(data1[0])+1)]
            corner.corner(data1[:,:], labels=labels[:],truths=median_vals[:],quantiles=[0.16, 0.5, 0.84],
            show_titles=False, title_kwargs={"fontsize": 12})
            plt.savefig(out_dir+'/'+obs+'_'+model_type+out_mcmc_name)
            plt.close()

            if freeze_nH==False:
                start_index=1
                m1.TBabs.nH     =median_vals[0]
            else:
                start_index=0

            j=0
            for i in range(len(kT_arr)):
                if freeze_kT_arr[i]==True:
                    m1(5+4*i).values=median_vals[start_index+j]
                else:
                    m1(2+4*i).values=median_vals[start_index+j]
                    j+=1
                    m1(5+4*i).values=median_vals[start_index+j]
                j+=1    
                
            m1_2.powerlaw.PhoIndex  =median_vals[start_index+j]
            m1_2.powerlaw.norm      =median_vals[start_index+j+1]
            m1_2.gaussian.norm      =median_vals[start_index+j+2]

            #################
            Plot.xAxis='keV'
            Plot("data resid")
            energies = Plot.x()
            edeltas = Plot.xErr()
            rates = Plot.y(1,1)
            errors = Plot.yErr(1,1)
            data_max    =max(rates)+max(errors)

            foldedmodel = Plot.model()
            dataLabels = Plot.labels(1)
            residLabels = Plot.labels(2)
            # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
            # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
            nE = len(energies)
            stepenergies = list()
            for j in range(nE):
                stepenergies.append(energies[j] - edeltas[j])
            stepenergies.append(energies[-1]+edeltas[-1])
            foldedmodel.append(foldedmodel[-1])
            resid = Plot.y(1,2)
            residerr = Plot.yErr(1,2)
            max_resid   =max(resid)+max(residerr)
            min_resid   =min(resid)-max(residerr)
    
            #print('The current ObsID is',obs)
            plt.subplot(221)
            plt.xscale('log')
            plt.yscale('log')
            plt.ylabel(dataLabels[1])
            plt.title(dataLabels[2])
            plt.errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.')
            plt.step(stepenergies,foldedmodel,where='post')
            plt.ylim(0.01,data_max)
            plt.subplot(223)
            plt.xscale('log')
            plt.xlabel(residLabels[0])
            plt.ylabel(residLabels[1])
            plt.errorbar(energies,resid,xerr=edeltas,yerr=residerr,fmt='.')
            plt.hlines(0.0,stepenergies[0],stepenergies[-1],linestyles='dashed')
            plt.ylim(min_resid,max_resid)

                #plt.figure()
            Plot.xAxis='keV'
            Plot("data resid")
            energies = Plot.x(2)
            edeltas = Plot.xErr(2)
            rates = Plot.y(2,1)
            errors = Plot.yErr(2,1)
            foldedmodel = Plot.model(2)
            dataLabels = Plot.labels(1)
            residLabels = Plot.labels(2)
            # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
            # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
            nE = len(energies)
            stepenergies = list()
            for i in range(nE):
                stepenergies.append(energies[i] - edeltas[i])
            stepenergies.append(energies[-1]+edeltas[-1])
            foldedmodel.append(foldedmodel[-1])
            resid = Plot.y(2,2)
            residerr = Plot.yErr(2,2)
            plt.subplot(222)
            plt.xscale('log')
            plt.yscale('log')
            #plt.ylabel(dataLabels[1])
            plt.title(dataLabels[2])
            plt.errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.')
            plt.step(stepenergies,foldedmodel,where='post')
            plt.ylim(0.01,data_max)
            plt.subplot(224)
            plt.xscale('log')
            plt.xlabel(residLabels[0])
            #plt.ylabel(residLabels[1])
            plt.errorbar(energies,resid,xerr=edeltas,yerr=residerr,fmt='.')
            plt.hlines(0.0,stepenergies[0],stepenergies[-1],linestyles='dashed')   
            plt.ylim(min_resid,max_resid)

            plt.savefig(out_dir+'/'+obs+'_'+model_type+out_spec_name)
            plt.close()

        if write_norm==True:
            f=open(out_dir+'/'+obs+'_'+model_type+norm_file_name,'w')
            f.write('### '+obs+'\n\n')
                
            if save_fig==False:
                data=fits.getdata(out_name,1)
                head=fits.getheader(out_name,1)
                data1=np.array([np.array(x) for x in data])
                median_vals =np.median(data1,axis=0)
            u_vals_1sig  =np.percentile(data1,84.0,axis=0)
            d_vals_1sig  =np.percentile(data1,16.0,axis=0)
            u_vals_2sig  =np.percentile(data1,97.5,axis=0)
            d_vals_2sig  =np.percentile(data1,2.5,axis=0)  
            u_vals_3sig  =np.percentile(data1,99.85,axis=0)
            d_vals_3sig  =np.percentile(data1,0.15,axis=0)
            my_str=''
            if freeze_nH==False:
                my_str+='## nH='+str(median_vals[0])+'\n'
                my_str+='## 1sigma nH interval='+str(u_vals_1sig[0])+'-'+str(d_vals_1sig[0])+'\n'
                my_str+='## 2sigma nH interval='+str(u_vals_2sig[0])+'-'+str(d_vals_2sig[0])+'\n'
                my_str+='## 3sigma nH interval='+str(u_vals_3sig[0])+'-'+str(d_vals_3sig[0])+'\n'
                start_index=1
            else:
                start_index=0
            my_str+='## Chi-sq='+str(median_vals[-1])+'\n'
            my_str+='## 1sigma Chi-sq interval='+str(u_vals_1sig[-1])+'-'+str(d_vals_1sig[-1])+'\n'
            my_str+='## 2sigma Chi-sq interval='+str(u_vals_2sig[-1])+'-'+str(d_vals_2sig[-1])+'\n'
            my_str+='## 3sigma chi-sq interval='+str(u_vals_3sig[-1])+'-'+str(d_vals_3sig[-1])+'\n'
                    
            my_str+='\n#kT\tkT_u\tkT_d\tNorm\t\tNorm_1u\t\tNorm_1d\t\tNorm_2u\t\tNorm_2d\t\tNorm_3u\t\tNorm_3d\n'            
            f.write(my_str)
            j=0
    
            for i in range(len(kT_arr)):
                if freeze_kT_arr[i]==True:
                    my_str=str(kT_arr[i])+'\t0\t0\t'
            
                    my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                else:
                    my_str=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'
                    j+=1
                    my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                j+=1
                f.write(my_str)

            f.close()

            AllChains.clear()
            AllData.clear()
            AllModels.clear()
    return

def calculate_flux_apec_multiresponse(ism_correct,obsid,clean_type,spec_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE0=0.4,endE0=10.0,bkg_type='3c50',input_wt='standard',out_dir='.',min_count=15):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'  
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'

    for N in range(len(obsid)):
        os.chdir(spec_dir)
        obs     =obsid[N]
        chain_name =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        obs0        =obs#obs.split('_')[0]
        data=fits.getdata(chain_name,1)   
        data=np.array([np.array(x) for x in data])
        ##############################
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg_grp'+str(min_count)+'.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            
        my_str      ='1:1 '+spec_file+' 2:2 '+bkg_file
        AllData(my_str)
        s1  =AllData(1)
        s2  =AllData(2)
        s1.background   ='none'
        s1.response     =rmf_file
        s1.response.arf =arf_file
        s2.response     =rmf_file
        s2.response.arf ='none'
        AllData.ignore("bad")
        s1.ignore(f"**-{startE0} {endE0}-**")
        s2.ignore(f"**-{startE0} {endE0}-**")
        ############################# Part to edit
        model_str='tbabs('
        for i in range(len(kT_arr)):
            if i==0:
                model_str+='apec'
            else:
                model_str+='+apec'
        model_str+=')' 

        Model(model_str)
        m1 = AllModels(1)
        m2 = AllModels(2)
        
        if freeze_nH==True:
            if ism_correct==True:
                m1.TBabs.nH=0.0
            else:    
                m1.TBabs.nH =nH0
            start_index =0

        for i in range(len(kT_arr)):
            m2(5+4*i).values=0.0

        s1.multiresponse[1] = rmf_file
        s2.multiresponse[1] = rmf_file
        m1_2 = Model("pow+gaussian",'myback',2)
        m1_2.gaussian.LineE=0.56
        m1_2.gaussian.Sigma=0.01

        for j in range(0,len(data),1000):
            if freeze_nH==False:
                if ism_correct==True:
                    m1.TBabs.nH=0.0
                else:
                    m1.TBabs.nH=data[j,0]
                start_index =1
                        
            l=0
            for i in range(len(kT_arr)):
                if freeze_kT_arr[i]==True:
                    m1(2+4*i).values=kT_arr[i]
                    m1(5+4*i).values=data[j,start_index+l]
                else:
                    m1(2+4*i).values=data[j,start_index+l]
                    l+=1
                    m1(5+4*i).values=data[j,start_index+l]
                l+=1  
            m1_2.powerlaw.PhoIndex  =data[j,start_index+l]
            m1_2.powerlaw.norm      =data[j,start_index+l+1]
            m1_2.gaussian.norm      =data[j,start_index+l+2]


            for E1, E2 in zip(startE_arr, endE_arr):
                if freeze_nH==True:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH_ISM_corrected.txt'
                else:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH_ISM_corrected.txt'
                if j==0:
                    f=open(out_file,'w')
                    f.write('#ObsID: '+obs+', startE='+str(startE0)+'keV, endE='+str(endE0)+'\n')
                    f.write('#targetFlux\tbkgFlux\n')
                else:
                    f=open(out_file,'a')
                
                AllModels.calcFlux(str(E1)+' '+str(E2)+' err')
                f.write(str(s1.flux[0])+'\t'+str(s1.flux[6])+'\n')
                f.close()
                
        AllData.clear()
        AllModels.clear()
        AllChains.clear()



def spectral_analysis_apec_pow(obsid,clean_type,spec_dir,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE=0.4,endE=10.0,bkg_type='3c50',chain_length=200000,burn_length=2000,input_wt='standard',out_dir='.',min_count=15,save_fig=False,write_norm=False):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_varying_nH_mcmc_with_nH_bound_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'
        out_mcmc_name       =clean_type+'_mcmc_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        out_spec_name       =clean_type+'_spectrum_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.png'
        norm_file_name      =clean_type+'_norm_values_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.txt'
    

    for N in range(len(obsid)):
        obs =obsid[N]
        out_name    =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        cont=os.path.isfile(out_name)
        if cont==True:
            continue
        else:
            os.chdir(spec_dir)
            bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
            rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
            arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
            s           =Spectrum(spec_file)
            s.background    =bkg_file
            s.response      =rmf_file
            s.response.arf  =arf_file
            AllData.ignore("bad")
            s.ignore(f"**-{startE} {endE}-**")
            ############################# Part to edit
            model_str='tbabs('
            for i in range(len(kT_arr)):
                if i==0:
                    model_str+='apec'
                else:
                    model_str+='+apec'
            model_str+='+powerlaw)'        
            m1=Model(model_str)
            m1.TBabs.nH =nH0
            m1.TBabs.nH.frozen  =freeze_nH
            if freeze_nH==False:
                m1.TBabs.nH =[nH0,0.01,0.0,0.0,10*nH0,10*nH0]
            for i in range(len(kT_arr)):
                m1(2+4*i).values=kT_arr[i]
                m1(2+4*i).frozen=freeze_kT_arr[i]
                if freeze_kT_arr[i]==False:
                    m1(2+4*i).values=[kT_arr[i],0.01,kT_arr[i]/3.,kT_arr[i]/3.,3*kT_arr[i],3*kT_arr[i]] # HARD-CODED
            #########################
            Fit.query = "yes"
            Fit.nIterations=100
            #Fit.statMethod = "cstat"
            Fit.weight  =input_wt   
            Fit.renorm()
            Fit.perform()
         
            AllChains.defBurn = burn_length
            AllChains.defLength = chain_length
            print('I am going to run chain')
    
            c1 = Chain(out_name,rand=False)
            
            if save_fig==True:
                data=fits.getdata(out_name,1)
                head=fits.getheader(out_name,1)
                data1=np.array([np.array(x) for x in data])
                median_vals =np.median(data1,axis=0)
                labels=[head['TTYPE'+str(i)] for i in range(1,len(data1[0])+1)]
                corner.corner(data1[:,:], labels=labels[:],truths=median_vals[:],quantiles=[0.16, 0.5, 0.84],
                       show_titles=False, title_kwargs={"fontsize": 12})
                plt.savefig(out_dir+'/'+obs+'_'+model_type+out_mcmc_name)
                plt.close()

                if freeze_nH==False:
                    start_index=1
                    m1.TBabs.nH     =median_vals[0]
                else:
                    start_index=0

                j=0
                for i in range(len(kT_arr)):
                    if freeze_kT_arr[i]==True:
                        m1(5+4*i).values=median_vals[start_index+j]
                    else:
                        m1(2+4*i).values=median_vals[start_index+j]
                        j+=1
                        m1(5+4*i).values=median_vals[start_index+j]
                    j+=1    
                m1.powerlaw.PhoIndex    =median_vals[start_index+j]
                m1.powerlaw.norm        =median_vals[start_index+j+1]
                ###############
                Plot.xAxis='keV'
                Plot("data resid")
                energies = Plot.x()
                edeltas = Plot.xErr()
                rates = Plot.y(1,1)
                errors = Plot.yErr(1,1)
                foldedmodel = Plot.model()
                dataLabels = Plot.labels(1)
                residLabels = Plot.labels(2)
                # note that for matplotlib step plots we need an x-axis array which includes the start and end value for each
                # bin and the y-axis has to be the same size with an extra value added equal to the value of the last bin
                nE = len(energies)
                stepenergies = list()
                for j in range(nE):
                    stepenergies.append(energies[j] - edeltas[j])
                stepenergies.append(energies[-1]+edeltas[-1])
                foldedmodel.append(foldedmodel[-1])
                resid = Plot.y(1,2)
                residerr = Plot.yErr(1,2)
                #print('The current ObsID is',obs)
                plt.subplot(211)
                plt.xscale('log')
                #plt.yscale('log')
                plt.ylabel(dataLabels[1])
                plt.title(dataLabels[2])
                plt.errorbar(energies,rates,xerr=edeltas,yerr=errors,fmt='.')
                plt.step(stepenergies,foldedmodel,where='post')
                plt.subplot(212)
                plt.xscale('log')
                plt.xlabel(residLabels[0])
                plt.ylabel(residLabels[1])
                plt.errorbar(energies,resid,xerr=edeltas,yerr=residerr,fmt='.')
                plt.hlines(0.0,stepenergies[0],stepenergies[-1],linestyles='dashed')
                plt.savefig(out_dir+'/'+obs+'_'+model_type+out_spec_name)
                plt.close()
            if write_norm==True:
                f=open(out_dir+'/'+obs+'_'+model_type+norm_file_name,'w')
                f.write('### '+obs+'\n\n')
                
                if save_fig==False:
                    data=fits.getdata(out_name,1)
                    head=fits.getheader(out_name,1)
                    data1=np.array([np.array(x) for x in data])
                    median_vals =np.median(data1,axis=0)
                u_vals_1sig  =np.percentile(data1,84.0,axis=0)
                d_vals_1sig  =np.percentile(data1,16.0,axis=0)
                u_vals_2sig  =np.percentile(data1,97.5,axis=0)
                d_vals_2sig  =np.percentile(data1,2.5,axis=0)  
                u_vals_3sig  =np.percentile(data1,99.85,axis=0)
                d_vals_3sig  =np.percentile(data1,0.15,axis=0)
                my_str=''
                if freeze_nH==False:
                    my_str+='## nH='+str(median_vals[0])+'\n'
                    my_str+='## 1sigma nH interval='+str(u_vals_1sig[0])+'-'+str(d_vals_1sig[0])+'\n'
                    my_str+='## 2sigma nH interval='+str(u_vals_2sig[0])+'-'+str(d_vals_2sig[0])+'\n'
                    my_str+='## 3sigma nH interval='+str(u_vals_3sig[0])+'-'+str(d_vals_3sig[0])+'\n'
                    start_index=1
                else:
                    start_index=0
                my_str+='## Chi-sq='+str(median_vals[-1])+'\n'
                my_str+='## 1sigma Chi-sq interval='+str(u_vals_1sig[-1])+'-'+str(d_vals_1sig[-1])+'\n'
                my_str+='## 2sigma Chi-sq interval='+str(u_vals_2sig[-1])+'-'+str(d_vals_2sig[-1])+'\n'
                my_str+='## 3sigma chi-sq interval='+str(u_vals_3sig[-1])+'-'+str(d_vals_3sig[-1])+'\n'
                    
                my_str+='\n#kT\tkT_2u\tkT_2d\tNorm\t\tNorm_1u\t\tNorm_1d\t\tNorm_2u\t\tNorm_2d\t\tNorm_3u\t\tNorm_3d\n'            
                f.write(my_str)
                j=0
    
                for i in range(len(kT_arr)):
                    if freeze_kT_arr[i]==True:
                        my_str=str(kT_arr[i])+'\t0\t0\t'
            
                        my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                    else:
                        my_str=str(median_vals[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'
                        j+=1
                        my_str+=str(median_vals[start_index+j])+'\t'+str(u_vals_1sig[start_index+j])+'\t'+str(d_vals_1sig[start_index+j])+'\t'+str(u_vals_2sig[start_index+j])+'\t'+str(d_vals_2sig[start_index+j])+'\t'+str(u_vals_3sig[start_index+j])+'\t'+str(d_vals_3sig[start_index+j])+'\n'
                    j+=1
                    f.write(my_str)

                f.close()

            AllChains.clear()
            AllData.clear()
            AllModels.clear()
    return


def calculate_flux_apec_pow(ism_correct,obsid,clean_type,spec_dir,startE_arr,endE_arr,model_type,kT_arr,freeze_kT_arr,freeze_nH,nH0,startE0=0.4,endE0=10.0,bkg_type='3c50',input_wt='standard',out_dir='.',min_count=15):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'  
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE0)+'_'+str(endE0)+'keV.fits'

    for N in range(len(obsid)):
        os.chdir(spec_dir)
        obs     =obsid[N]
        chain_name =out_dir+'/chain_'+obs+'_'+model_type+out_chain_name
        obs0        =obs#obs.split('_')[0]
        data=fits.getdata(chain_name,1)   
        data=np.array([np.array(x) for x in data])
        ##############################
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            
        s           =Spectrum(spec_file)
        s.background    =bkg_file
        s.response      =rmf_file
        s.response.arf  =arf_file
        AllData.ignore("bad")
        s.ignore(f"**-{startE0} {endE0}-**")
        ############################# Part to edit
        model_str='tbabs('
        for i in range(len(kT_arr)):
            if i==0:
                model_str+='apec'
            else:
                model_str+='+apec'
        model_str+='+powerlaw)'        
        m1=Model(model_str)
        for j in range(0,len(data),1000):
            if freeze_nH==False:
                if ism_correct==True:
                    m1.TBabs.nH=0.0
                else:    
                    m1.TBabs.nH=data[j,0]
                start_index =1
            else:
                if ism_correct==True:
                    m1.TBabs.nH=0.0
                else:
                    m1.TBabs.nH =nH0
                start_index =0
            
            l=0
            for i in range(len(kT_arr)):
                if freeze_kT_arr[i]==True:
                    m1(2+4*i).values=kT_arr[i]
                    m1(5+4*i).values=data[j,start_index+l]
                else:
                    m1(2+4*i).values=data[j,start_index+l]
                    l+=1
                    m1(5+4*i).values=data[j,start_index+l]
                l+=1  
            m1.powerlaw.PhoIndex    =data[j,start_index+l]
            m1.powerlaw.norm        =data[j,start_index+l+1]

            for E1, E2 in zip(startE_arr, endE_arr):
                if freeze_nH==True:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_fixed_nH_ISM_corrected.txt'
                else:
                    if ism_correct==False:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH.txt'
                    else:
                        out_file    =out_dir+'/'+obs+clean_type+'_'+model_type+'_flux_'+str(E1)+'_'+str(E2)+'keV'+'_varying_nH_ISM_corrected.txt'
                if j==0:
                    f=open(out_file,'w')
                    f.write('#ObsID: '+obs+', startE='+str(startE0)+'keV, endE='+str(endE0)+'\n')
                    f.write('#targetFlux\n')
                else:
                    f=open(out_file,'a')
                
                AllModels.calcFlux(str(E1)+' '+str(E2)+' err')
                f.write(str(s.flux[0])+'\n')
                f.close()
                
        AllData.clear()
        AllModels.clear()
        AllChains.clear()














def get_energy_bin(obsid,spec_dir,clean_type,min_count,startE,endE): ### only for background subtracted spectra, 3C50 model
    num_bin=np.zeros(len(obsid))
    for N in range(len(obsid)):
        os.chdir(spec_dir)
        obs     =obsid[N]
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            
        s           =Spectrum(spec_file)
        s.background    =bkg_file
        s.response      =rmf_file
        s.response.arf  =arf_file
        AllData.ignore("bad")
        s.ignore(f"**-{startE} {endE}-**")
        Plot('data')
        energies=Plot.x()
        num_bin[N]=len(energies)
        AllData.clear()
    return num_bin    
    



def check_chi_sq(obsid,spec_dir,clean_type='_flagged_high_energy',model_type='4apec_naze2014',kT_arr=[0.2,0.6,1.0,4.0],freeze_kT_arr=[True,True,True,True],freeze_nH=True,startE=0.3,endE=10.0,startE0=0.5,endE0=2.0,nH0=0.03,bkg_type='3c50',input_wt='standard',confidence=1,min_count=15,model_pow=False):
    if freeze_nH==False:
        out_chain_name      =clean_type+'_varying_nH_with_nH_bound_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'  
    else:
        out_chain_name      =clean_type+'_'+input_wt+'_'+str(startE)+'_'+str(endE)+'keV.fits'

    for N in range(len(obsid)):
        os.chdir(spec_dir)
        obs     =obsid[N]
        chain_name ='chain_'+obs+'_'+model_type+out_chain_name
        obs0        =obs#obs.split('_')[0]
        data=fits.getdata(chain_name,1)   
        data1=np.array([np.array(x) for x in data])
        median_vals =np.median(data1,axis=0)
        if confidence==1:
            u_vals  =np.percentile(data1,84.0,axis=0)
            d_vals  =np.percentile(data1,16.0,axis=0)
        if confidence==2:    
            u_vals  =np.percentile(data1,97.5,axis=0)
            d_vals  =np.percentile(data1,2.5,axis=0)  
        if confidence==3:
            u_vals  =np.percentile(data1,99.85,axis=0)
            d_vals  =np.percentile(data1,0.15,axis=0)
        ##############################
        spec_file   ='nibackgen3C50_'+obs+clean_type+'_tot_grp'+str(min_count)+'.pi'
        bkg_file    ='nibackgen3C50_'+obs+clean_type+'_bkg.pi'
        rmf_file    ='nibackgen3C50_'+obs+'_tot.rmf'
        arf_file    ='nibackgen3C50_'+obs+'_tot.arf'
            
        s           =Spectrum(spec_file)
        s.background    =bkg_file
        s.response      =rmf_file
        s.response.arf  =arf_file
        AllData.ignore("bad")
        s.ignore(f"**-{startE0} {endE0}-**")
        ############################# Part to edit
        model_str='tbabs('
        for i in range(len(kT_arr)):
            if i==0:
                model_str+='apec'
            else:
                model_str+='+apec'
        if model_pow==False:        
            model_str+=')'        
        else:
            model_str+='+powerlaw)'
        m1=Model(model_str)

        ### Median chi-sq
        
        if freeze_nH==False:
            m1.TBabs.nH=median_vals[0]
            start_index =1
        else:
            m1.TBabs.nH =nH0
            start_index =0
        l=0    
        for i in range(len(kT_arr)):
            if freeze_kT_arr[i]==True:
                m1(2+4*i).values=kT_arr[i]
                m1(5+4*i).values=median_vals[start_index+l]
            else:
                m1(2+4*i).values=median_vals[start_index+l]
                l+=1
                m1(5+4*i).values=median_vals[start_index+l]
            l+=1
        if model_pow==True:
            m1.powerlaw.PhoIndex=median_vals[start_index+l]
            m1.powerlaw.norm    =median_vals[start_index+l+1]
        m1.show()
        ans=input('Have you noted the chi-sq for median values? (press any to continue)')
        '''
        l=0
        if freeze_nH==False:
            m1.TBabs.nH=d_vals[0]
            start_index =1
        else:
            m1.TBabs.nH =nH0
            start_index =0
        for i in range(len(kT_arr)):
            if freeze_kT_arr[i]==True:
                m1(2+4*i).values=kT_arr[i]
                m1(5+4*i).values=d_vals[start_index+l]
            else:
                m1(2+4*i).values=d_vals[start_index+l]
                l+=1
                m1(5+4*i).values=d_vals[start_index+l]
            l+=1
        m1.show()
        ans=input('Have you noted the chi-sq for CI below for obsid? (press any to continue)')

        l=0
        if freeze_nH==False:
            m1.TBabs.nH=u_vals[0]
            start_index =1
        else:
            m1.TBabs.nH =nH0
            start_index =0
        for i in range(len(kT_arr)):
            if freeze_kT_arr[i]==True:
                m1(2+4*i).values=kT_arr[i]
                m1(5+4*i).values=u_vals[start_index+l]
            else:
                m1(2+4*i).values=u_vals[start_index+l]
                l+=1
                m1(5+4*i).values=u_vals[start_index+l]
            l+=1
        m1.show()
        ans=input('Have you noted the chi-sq for sigma CI above for obsid? (press any to continue)')
        '''
        ans=input('Do you want to continue? (press 1 for yes and 0 for no)? ')
        AllData.clear()
        AllModels.clear()
        AllChains.clear()
        if ans=='0':
            return

def get_flux_med_CI(spec_dir,obsid,clean_type,model_type,startE,endE,confidence=1,freeze_nH=True,ism_correct=False):
    os.chdir(spec_dir)
    flux,flux_u,flux_d  =np.zeros(len(obsid)),np.zeros(len(obsid)),np.zeros(len(obsid))

    for N in range(len(obsid)):
        obs =obsid[N]
        if freeze_nH==True:
            if ism_correct==False:
                out_file    =obs+clean_type+'_'+model_type+'_flux_'+str(startE)+'_'+str(endE)+'keV'+'_fixed_nH.txt'
            else:
                out_file    =obs+clean_type+'_'+model_type+'_flux_'+str(startE)+'_'+str(endE)+'keV'+'_fixed_nH_ISM_corrected.txt'
        else:
            if ism_correct==False:
                out_file    =obs+clean_type+'_'+model_type+'_flux_'+str(startE)+'_'+str(endE)+'keV'+'_varying_nH.txt'
            else:
                out_file    =obs+clean_type+'_'+model_type+'_flux_'+str(startE)+'_'+str(endE)+'keV'+'_varying_nH_ISM_corrected.txt'

        

        data    =np.loadtxt(out_file)
        flux[N]=np.median(data)
        if confidence==1:
            flux_u[N]  =np.percentile(data,84.0,axis=0)
            flux_d[N]  =np.percentile(data,16.0,axis=0)
        if confidence==2:    
            flux_u[N]  =np.percentile(data,97.5,axis=0)
            flux_d[N]  =np.percentile(data,2.5,axis=0)  
        if confidence==3:
            flux_u[N]  =np.percentile(data,99.85,axis=0)
            flux_d[N]  =np.percentile(data,0.15,axis=0)
    return flux,flux_u,flux_d
                
