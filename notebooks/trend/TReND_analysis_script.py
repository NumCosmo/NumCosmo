import sys
sys.path.insert(0, "/global/homes/c/cinlima/NumCosmo/notebooks/richness_proxy/")
sys.path.insert(0, "/global/homes/c/cinlima/NumCosmo/notebooks/richness_proxy/ESMCMC")

from esmcmc_rm_relation_script import catalog_fit, esmcmc

#NumCosmo
from numcosmo_py import Ncm, Nc, GObject
Ncm.cfg_init()
Ncm.cfg_set_log_handler(lambda msg: sys.stdout.write(msg) and sys.stdout.flush())

#Useful packages
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib as mpl
from getdist.mcsamples import  MCSamples

#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

class TReNDAnalysis:
    
    def __init__(self, MASS_CUT_list, RICH_CUT_list):
        self.RICH_CUT_list = RICH_CUT_list
        self.MASS_CUT_list = MASS_CUT_list



#----------------------------------------------------------------------------------------------------------------#
    def mcmc_results(self, model_type):
    
        mcmc_list = []
        
        for mcut in self.MASS_CUT_list:
              
            mcmc_i = pd.DataFrame(data = {'Min_Mass':[], 'Min_Richness': [],  
                                                        'mu0':[], 'mu0_sup': [], 'mu0_inf': [], 
                                                        'mu1':[], 'mu1_sup': [], 'mu1_inf': [], 
                                                        'mu2':[], 'mu2_sup': [], 'mu2_inf': [],
                                                        'sigma0':[], 'sigma0_sup': [], 'sigma0_inf': [], 
                                                        'sigma1':[], 'sigma1_sup': [], 'sigma1_inf': [], 
                                                        'sigma2':[], 'sigma2_sup': [], 'sigma2_inf': []})
             
            for rcut in self.RICH_CUT_list:
                    
                
                #------ NumCosmo's MCMC catalog  ----------------------------------------------------------------------------#
                N_WALKERS = 1200
                N_RUN = 300
                    
                RICH_CUT = rcut
                MASS_CUT = mcut 
    
                FILE_NAME = "/global/homes/c/cinlima/ESMCMC/"+model_type+"/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(int(MASS_CUT))+".fits"
 
                burnin_cat = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME, 0.0).peek_e_mean_stats().estimate_const_break(0)  
                
                mcat = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME, N_WALKERS * burnin_cat) 
                
                #------ Best fit  ----------------------------------------------------------------------------#
                mu0 = mcat.get_bestfit_row().dup_array()[1]
                mu1 = mcat.get_bestfit_row().dup_array()[2]
                mu2 = mcat.get_bestfit_row().dup_array()[3]        
                sigma0 = mcat.get_bestfit_row().dup_array()[4]            
                sigma1 = mcat.get_bestfit_row().dup_array()[5]            
                sigma2 = mcat.get_bestfit_row().dup_array()[6]                        


                #------ GetDist's MCMC Samples  ----------------------------------------------------------------------------#
                data_fit_full = pd.DataFrame(fits.open(FILE_NAME)[1].data).iloc[:, 1:7].T
                data_fit_void = np.array(data_fit_full)
                data_fit = []
                for item in data_fit_void:
                    arr= np.array(item)
                    data_fit.append(np.asarray(arr.tolist()))
                    
                names = ['1', '2', '3', '4', '5', '6']
                labels=[r"\mu_{0}", r"\mu_{1}", r"\mu_{2}", r"\sigma_{0}", r"\sigma_{1}", r"\sigma_{2}"]
                
                settings = {
                    "mult_bias_correction_order": 0,
                    "smooth_scale_2D": 3,
                    "smooth_scale_1D": 3,
                    "boundary_correction_order": 0,
                }
                
                samples1 = MCSamples(samples=data_fit, names=names, labels=labels, settings=settings)
                samples1.removeBurn(0.3)

                #------ 2 sigma interval ---------------------------------------------------------------------#
                mu0_sup = samples1.get1DDensity(0).getLimits(0.95)[0]
                mu0_inf = samples1.get1DDensity(0).getLimits(0.95)[1]
            
                mu1_sup = samples1.get1DDensity(1).getLimits(0.95)[0]
                mu1_inf = samples1.get1DDensity(1).getLimits(0.95)[1]
    
                mu2_inf = samples1.get1DDensity(2).getLimits(0.95)[0]
                mu2_sup = samples1.get1DDensity(2).getLimits(0.95)[1]
            
                sigma0_inf = samples1.get1DDensity(3).getLimits(0.95)[0]
                sigma0_sup = samples1.get1DDensity(3).getLimits(0.95)[1]
            
                sigma1_inf = samples1.get1DDensity(4).getLimits(0.95)[0]
                sigma1_sup = samples1.get1DDensity(4).getLimits(0.95)[1]
            
                sigma2_inf = samples1.get1DDensity(5).getLimits(0.95)[0]
                sigma2_sup = samples1.get1DDensity(5).getLimits(0.95)[1]

                
                #------ MCMC dataframe ----------------------------------------------------------------------#    
                mcmc_i = pd.concat([mcmc_i, pd.DataFrame([{'Min_Mass':mcut, 'Min_Richness': rcut,   
                                                        'mu0':mu0, 'mu0_sup': mu0_sup, 'mu0_inf': mu0_inf, 
                                                        'mu1':mu1, 'mu1_sup': mu1_sup, 'mu1_inf': mu1_inf, 
                                                        'mu2':mu2, 'mu2_sup': mu2_sup, 'mu2_inf': mu2_inf,
                                                        'sigma0':sigma0, 'sigma0_sup': sigma0_sup, 'sigma0_inf': sigma0_inf, 
                                                        'sigma1':sigma1, 'sigma1_sup': sigma1_sup, 'sigma1_inf': sigma1_inf, 
                                                        'sigma2':sigma2, 'sigma2_sup': sigma2_sup, 'sigma2_inf': sigma2_inf}])], ignore_index=True)
                    
            mcmc_list.append(mcmc_i)
         
        
        return pd.concat(mcmc_list).reset_index(drop=True)
    
    

#----------------------------------------------------------------------------------------------------------------#
    
    def bayes_factor(self):
        bef_list = []
            
        for mcut in self.MASS_CUT_list:
             
            bef = pd.DataFrame(data = {'Min_Mass':[], 'Min_Richness': [],
                                       'BEQ': [], 'BEQ Err': [], 'BEL': [], 
                                       'BEL Err':[], 'ln(BF)': []})
                   
            for rcut in self.RICH_CUT_list:
                    
                N_WALKERS = 1200
                N_RUN = 300
                    
                RICH_CUT = rcut
                MASS_CUT = mcut 
        
                FILE_NAME1 = "/global/homes/c/cinlima/ESMCMC/with_correction/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(int(MASS_CUT))+".fits"   
                FILE_NAME2 = "/global/homes/c/cinlima/ESMCMC/without_correction/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(int(MASS_CUT))+".fits"
           
                burnin_cat1 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME1, 0.0).peek_e_mean_stats().estimate_const_break(0)   
                mcat1 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME1, N_WALKERS * burnin_cat1)
                be1, post_lnnorm_sd1 = mcat1.get_post_lnnorm()
            
                burnin_cat2 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME2, 0.0).peek_e_mean_stats().estimate_const_break(0)  
                mcat2 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME2, N_WALKERS * burnin_cat2)
                be2, post_lnnorm_sd2 = mcat2.get_post_lnnorm()
    
                # ln (BF) - BF is the Bayes Factor.
                lnbf = be1 - be2 
                
                bef = pd.concat([bef, pd.DataFrame([{'Min_Mass':mcut, 'Min_Richness': rcut,   
                                                        'BEQ': be1, 'BEQ Err': post_lnnorm_sd1, 'BEL': be2, 
                                                        'BEL Err':post_lnnorm_sd2, 'ln(BF)': lnbf}])], ignore_index=True)
    
            bef_list.append(bef)
    
        return pd.concat(bef_list).reset_index(drop=True)

        
#----------------------------------------------------------------------------------------------------------------#
    
    def latex_notation(self, x, precision=2):
    
        if x == 0:
            return f"{0:.{precision}f}"
    
        sci_str = f"{x:.{precision}e}"
        base, exponent = sci_str.split('e')
        exponent = int(exponent)
    
        return f"${base} \\times 10^{{{exponent}}}$"

#----------------------------------------------------------------------------------------------------------------#
    
    def mcmc_results_table(self, mcmc_table):
         
        for i in range(0, 3):
            
            mcmc_table['mu'+str(i)+'_sup'] = mcmc_table['mu'+str(i)+'_sup'] - mcmc_table['mu'+str(i)+'']
            mcmc_table['mu'+str(i)+'_inf'] = mcmc_table['mu'+str(i)+'_inf'] - mcmc_table['mu'+str(i)+'']
            
            mcmc_table['mu'+str(i)+''] = mcmc_table.apply(lambda row: f"${row['mu'+str(i)+'']:.2f}^{{+{row['mu'+str(i)+'_sup']:.2f}}}_{{{row['mu'+str(i)+'_inf']:.2f}}}$", axis=1)
            mcmc_table['sigma'+str(i)+'_sup'] = mcmc_table['sigma'+str(i)+'_sup'] - mcmc_table['sigma'+str(i)+'']
            mcmc_table['sigma'+str(i)+'_inf'] = mcmc_table['sigma'+str(i)+'_inf'] - mcmc_table['sigma'+str(i)+'']
             
            mcmc_table['sigma'+str(i)+''] = mcmc_table.apply(lambda row: f"${row['sigma'+str(i)+'']:.2f}^{{+{row['sigma'+str(i)+'_sup']:.2f}}}_{{{row['sigma'+str(i)+'_inf']:.2f}}}$", axis=1)

                
        mcmc_table['Min_Mass'] = mcmc_table['Min_Mass'].apply(self.latex_notation)
        
        
        # pd.options.display.float_format = '{:.2f}'.format
        
        
        mcmc_table_df = mcmc_table.rename(columns={'Min_Mass': '$M_c$', 'Min_Richness': '$\lambda_c$', 
                                                     'mu0': '$\mu_0$', 'mu1': '$\mu_1$','mu2': '$\mu_2$', 
                                                     'sigma0': '$\sigma_0$', 'sigma1': '$\sigma_1$','sigma2': '$\sigma_2$' })
        
        mcmc_table_df = mcmc_table_df.drop(['mu0_sup', 'mu1_sup', 'mu2_sup', 
                                            'mu0_inf', 'mu1_inf', 'mu2_inf', 
                                            'sigma0_sup', 'sigma1_sup', 'sigma2_sup', 
                                            'sigma0_inf', 'sigma1_inf', 'sigma2_inf'], axis=1)
        
        
        return mcmc_table_df.reset_index(drop=True)
    
    
    
#----------------------------------------------------------------------------------------------------------------#
    
    def get_global_ylim(self, mcmc_list1, mcmc_list2, key_inf, key_sup):
        all_vals = []
        for bef in mcmc_list1 + mcmc_list2:
            all_vals.extend(bef[key_inf])
            all_vals.extend(bef[key_sup])
        return np.min(all_vals), np.max(all_vals)
    

  
#----------------------------------------------------------------------------------------------------------------#
    
    def comparison_plot(self, calib1, calib2):    
        fig, axs = plt.subplots(3, 2, figsize=(10,8), sharex=True, sharey='row')
        fig.subplots_adjust(wspace=0.30, hspace=0) 
    
        mcmc_list1 = [calib1.iloc[i:i+7] for i in range(0, len(calib1), 7)]
        mcmc_list2 = [calib2.iloc[i:i+7] for i in range(0, len(calib2), 7)]
        
        ylim_mu0 = self.get_global_ylim(mcmc_list1, mcmc_list2, "mu0_inf", "mu0_sup")
        ylim_mu1 = self.get_global_ylim(mcmc_list1, mcmc_list2, "mu1_inf", "mu1_sup")
        ylim_mu2 = self.get_global_ylim(mcmc_list1, mcmc_list2, "mu2_inf", "mu2_sup")
    
        
        for i in range(5):
            axs[0, 0].plot(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu0"], label=f'{self.MASS_CUT_list[i]:.2}', ls='-.', lw=0.5, marker='o', ms=3)
            axs[0, 0].fill_between(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu0_inf"], mcmc_list1[i]["mu0_sup"], alpha=0.2)
        
            axs[1, 0].plot(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu1"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[1, 0].fill_between(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu1_inf"], mcmc_list1[i]["mu1_sup"], alpha=0.2)
        
            axs[2, 0].plot(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu2"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[2, 0].fill_between(mcmc_list1[i]["Min_Richness"], mcmc_list1[i]["mu2_inf"], mcmc_list1[i]["mu2_sup"], alpha=0.2)
        
        for i in range(5):
            axs[0, 1].plot(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu0"], label=f'{self.MASS_CUT_list[i]:.2}', ls='-.', lw=0.5, marker='o', ms=3)
            axs[0, 1].fill_between(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu0_inf"], mcmc_list2[i]["mu0_sup"], alpha=0.2)
        
            axs[1, 1].plot(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu1"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[1, 1].fill_between(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu1_inf"], mcmc_list2[i]["mu1_sup"], alpha=0.2)
        
            axs[2, 1].plot(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu2"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[2, 1].fill_between(mcmc_list2[i]["Min_Richness"], mcmc_list2[i]["mu2_inf"], mcmc_list2[i]["mu2_sup"], alpha=0.2)
        
        axs[0, 0].set_ylim(ylim_mu0)
        axs[0, 1].set_ylim(ylim_mu0)
        axs[1, 0].set_ylim(ylim_mu1)
        axs[1, 1].set_ylim(ylim_mu1)
        axs[2, 0].set_ylim(ylim_mu2)
        axs[2, 1].set_ylim(ylim_mu2)
        
        for row in range(3):
            axs[row, 0].set_ylabel(f'$\mu{row}$', fontsize=12)
            for col in range(2):
                axs[row, col].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f'{x:.2f}'))
                axs[row, col].tick_params(axis='y', labelleft=True)  # força números do eixo y
                axs[row, col].tick_params(axis='x')
        
        axs[2, 0].set_xlabel('$\lambda_{c}$', fontsize=12)
        axs[2, 1].set_xlabel('$\lambda_{c}$', fontsize=12)
        
        axs[0, 0].set_title("Simple log-normal")
        axs[0, 1].set_title("TReND")
        
        handles, labels = axs[0, 0].get_legend_handles_labels()
        lgd = fig.legend(handles, labels,  bbox_to_anchor=(0.45, 0.75), framealpha=1.0, title='Mass cut:')
        loc = "center"
    
        plt.show()
    
    

#----------------------------------------------------------------------------------------------------------------#    
    def bayes_factor_plot(self):
        plt.figure(figsize=(10,6))
    
        bf_table = self.bayes_factor()
        
        bef_list = [bf_table.iloc[i:i+7] for i in range(0, len(bf_table), 7)]
        
        for i in range(0,5): 
            
            plt.plot(bef_list[i]["Min_Richness"], bef_list[i]["ln(BF)"], label = f'{self.MASS_CUT_list[i]:.2}', ls = '-.', linewidth = 1.0, marker = 'o')
            plt.yscale("log")
            
        linewidth = 1.0            
        plt.axline((0, np.log(10)), (50, np.log(10)), c = 'darkcyan', ls= '--', label = '$\ln (10)$')
        plt.axline((0, np.log(100)), (50, np.log(100)), c = 'r', ls= '--', label = '$\ln (100)$')
        
        lgd = plt.legend(fontsize=12, bbox_to_anchor=(1.0, 1.0)) 
        
        plt.ylabel('$\ln BF$', fontsize=12)
        plt.xlabel(r'$\lambda_{c}$', fontsize=12)
        # plt.title('Bayes factor' )
        
        plt.grid()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        plt.show()
        