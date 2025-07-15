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

class ThresholdAnalysis:
    
    def __init__(self, MASS_CUT_list, RICH_CUT_list):
        self.RICH_CUT_list = RICH_CUT_list
        self.MASS_CUT_list = MASS_CUT_list



#----------------------------------------------------------------------------------------------------------------#

    def calibration_analysis(self, model_type):
    
        bef_list = []
        
        for mcut in self.MASS_CUT_list:
              
            bef = pd.DataFrame(data = {'Min_Mass':[], 'Min_Richness': [],  
                                                        'mu0':[], 'mu0_sup': [], 'mu0_inf': [], 
                                                        'mu1':[], 'mu1_sup': [], 'mu1_inf': [], 
                                                        'mu2':[], 'mu2_sup': [], 'mu2_inf': [],
                                                        'sigma0':[], 'sigma0_sup': [], 'sigma0_inf': [], 
                                                        'sigma1':[], 'sigma1_sup': [], 'sigma1_inf': [], 
                                                        'sigma2':[], 'sigma2_sup': [], 'sigma2_inf': []})
             
            for rcut in self.RICH_CUT_list:
                    
                N_WALKERS = 1200
                N_RUN = 300
                    
                RICH_CUT = rcut
                MASS_CUT = mcut 
    
                FILE_NAME = "/global/homes/c/cinlima/ESMCMC/"+model_type+"/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(MASS_CUT)+".fits"

                
                burnin_cat = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME, 0.0).peek_e_mean_stats().estimate_const_break(0) + 10  
                mcat = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME, N_WALKERS * burnin_cat)
                be2, post_lnnorm_sd2 = mcat.get_post_lnnorm()
                # lnevol2, glnvol2 = mcat2.get_post_lnvol(0.6827)
            
                mu0 = mcat.get_bestfit_row().dup_array()[1]
                # sd_mu0 = mcat2.peek_pstats().get_sd(1)
            
                mu1 = mcat.get_bestfit_row().dup_array()[2]
                # sd_mu1 = mcat2.peek_pstats().get_sd(2)
            
                mu2 = mcat.get_bestfit_row().dup_array()[3]
                # sd_mu2 = mcat2.peek_pstats().get_sd(3)
        
                sigma0 = mcat.get_bestfit_row().dup_array()[4]
                # sd_mu0 = mcat2.peek_pstats().get_sd(1)
            
                sigma1 = mcat.get_bestfit_row().dup_array()[5]
                # sd_mu1 = mcat2.peek_pstats().get_sd(2)
            
                sigma2 = mcat.get_bestfit_row().dup_array()[6]
                # sd_mu2 = mcat2.peek_pstats().get_sd(3)
                        
                    
                data_fit_full = pd.DataFrame(fits.open(FILE_NAME)[1].data).iloc[:, 1:7].T
                data_fit_void = np.array(data_fit_full)
                data_fit = []
                for item in data_fit_void:
                    arr= np.array(item)
                    data_fit.append(np.asarray(arr.tolist()))
                    
                names = [
                    '1',
                    '2',
                    '3',
                    '4',
                    '5',
                    '6',
                ]
                labels=[r"\mu_{0}", r"\mu_{1}", r"\mu_{2}", r"\sigma_{0}", r"\sigma_{1}", r"\sigma_{2}"]
                settings = {
                    "mult_bias_correction_order": 0,
                    "smooth_scale_2D": 3,
                    "smooth_scale_1D": 3,
                    "boundary_correction_order": 0,
                }
                samples1 = MCSamples(samples=data_fit, names=names, labels=labels, settings=settings)
                samples1.removeBurn(0.3)
                    
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
        
                    
                            
                bef = pd.concat([bef, pd.DataFrame([{'Min_Mass':mcut, 'Min_Richness': rcut,   
                                                        'mu0':mu0, 'mu0_sup': mu0_sup, 'mu0_inf': mu0_inf, 
                                                        'mu1':mu1, 'mu1_sup': mu1_sup, 'mu1_inf': mu1_inf, 
                                                        'mu2':mu2, 'mu2_sup': mu2_sup, 'mu2_inf': mu2_inf,
                                                        'sigma0':sigma0, 'sigma0_sup': sigma0_sup, 'sigma0_inf': sigma0_inf, 
                                                        'sigma1':sigma1, 'sigma1_sup': sigma1_sup, 'sigma1_inf': sigma1_inf, 
                                                        'sigma2':sigma2, 'sigma2_sup': sigma2_sup, 'sigma2_inf': sigma2_inf}])], ignore_index=True)
                    
            bef_list.append(bef)
        bf_datafame = pd.concat(bef_list)
        
        return bf_datafame.reset_index(drop=True)
    
    

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
        
                FILE_NAME1 = "/global/homes/c/cinlima/ESMCMC/with_correction/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(MASS_CUT)+".fits"   
                FILE_NAME2 = "/global/homes/c/cinlima/ESMCMC/without_correction/asc_rmin_"+str(RICH_CUT)+"_mmin_"+str(MASS_CUT)+".fits"
           
                burnin_cat1 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME1, 0.0).peek_e_mean_stats().estimate_const_break(0) + 10  
                mcat1 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME1, N_WALKERS * burnin_cat1)
                be1, post_lnnorm_sd1 = mcat1.get_post_lnnorm()
            
                burnin_cat2 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME2, 0.0).peek_e_mean_stats().estimate_const_break(0) + 10  
                mcat2 = Ncm.MSetCatalog.new_from_file_ro(FILE_NAME2, N_WALKERS * burnin_cat2)
                be2, post_lnnorm_sd2 = mcat2.get_post_lnnorm()
    
                # ln (BF)
                lnbf = be1 - be2 
                
                bef = pd.concat([bef, pd.DataFrame([{'Min_Mass':mcut, 'Min_Richness': rcut,   
                                                        'BEQ': be1, 'BEQ Err': post_lnnorm_sd1, 'BEL': be2, 
                                                        'BEL Err':post_lnnorm_sd2, 'ln(BF)': lnbf}])], ignore_index=True)
    
            bef_list.append(bef)
        
        bf_datafame = pd.concat(bef_list)
    
        return bf_datafame.reset_index(drop=True)
    
    
        
#----------------------------------------------------------------------------------------------------------------#
    
    def float_to_latex_sci(self, x, precision=2):
    
        if x == 0:
            return f"{0:.{precision}f}"
    
        sci_str = f"{x:.{precision}e}"
        base, exponent = sci_str.split('e')
        exponent = int(exponent)
    
        return f"${base} \\times 10^{{{exponent}}}$"

    

#----------------------------------------------------------------------------------------------------------------#
    
    def calibration_analysis_table(self, paper_trend):
         
        paper_trend['mu0_sup'] = paper_trend['mu0_sup'] - paper_trend['mu0']
        paper_trend['mu1_sup'] = paper_trend['mu1_sup'] - paper_trend['mu1']
        paper_trend['mu2_sup'] = paper_trend['mu2_sup'] - paper_trend['mu2']
        
        paper_trend['mu0_inf'] = paper_trend['mu0_inf'] - paper_trend['mu0']
        paper_trend['mu1_inf'] = paper_trend['mu1_inf'] - paper_trend['mu1']
        paper_trend['mu2_inf'] = paper_trend['mu2_inf'] - paper_trend['mu2']
        
        paper_trend['mu0'] = paper_trend.apply(lambda row: f"${row['mu0']:.2f}^{{+{row['mu0_sup']:.2f}}}_{{{row['mu0_inf']:.2f}}}$", axis=1)
        paper_trend['mu1'] = paper_trend.apply(lambda row: f"${row['mu1']:.2f}^{{+{row['mu1_sup']:.2f}}}_{{{row['mu1_inf']:.2f}}}$", axis=1)
        paper_trend['mu2'] = paper_trend.apply(lambda row: f"${row['mu2']:.2f}^{{+{row['mu2_sup']:.2f}}}_{{{row['mu2_inf']:.2f}}}$", axis=1)
        
        paper_trend['sigma0_sup'] = paper_trend['sigma0_sup'] - paper_trend['sigma0']
        paper_trend['sigma1_sup'] = paper_trend['sigma1_sup'] - paper_trend['sigma1']
        paper_trend['sigma2_sup'] = paper_trend['sigma2_sup'] - paper_trend['sigma2']
        
        paper_trend['sigma0_inf'] = paper_trend['sigma0_inf'] - paper_trend['sigma0']
        paper_trend['sigma1_inf'] = paper_trend['sigma1_inf'] - paper_trend['sigma1']
        paper_trend['sigma2_inf'] = paper_trend['sigma2_inf'] - paper_trend['sigma2']
        
        paper_trend['sigma0'] = paper_trend.apply(lambda row: f"${row['sigma0']:.2f}^{{+{row['sigma0_sup']:.2f}}}_{{{row['sigma0_inf']:.2f}}}$", axis=1)
        paper_trend['sigma1'] = paper_trend.apply(lambda row: f"${row['sigma1']:.2f}^{{+{row['sigma1_sup']:.2f}}}_{{{row['sigma1_inf']:.2f}}}$", axis=1)
        paper_trend['sigma2'] = paper_trend.apply(lambda row: f"${row['sigma2']:.2f}^{{+{row['sigma2_sup']:.2f}}}_{{{row['sigma2_inf']:.2f}}}$", axis=1)
        
        # Aplicar ao DataFrame
        paper_trend['Min_Mass'] = paper_trend['Min_Mass'].apply(self.float_to_latex_sci)
        
        
        # pd.options.display.float_format = '{:.2f}'.format
        
        
        paper_trend_df = paper_trend.rename(columns={'Min_Mass': '$M_c$', 'Min_Richness': '$\lambda_c$', 
                                                     'mu0': '$\mu_0$', 'mu1': '$\mu_1$','mu2': '$\mu_2$', 
                                                     'sigma0': '$\sigma_0$', 'sigma1': '$\sigma_1$','sigma2': '$\sigma_2$' })
        
        paper_trend_df = paper_trend_df.drop(['mu0_sup', 'mu1_sup', 'mu2_sup', 'mu0_inf', 'mu1_inf', 'mu2_inf', 
                                              'sigma0_sup', 'sigma1_sup', 'sigma2_sup', 'sigma0_inf', 'sigma1_inf', 'sigma2_inf'], axis=1)
        
        
        return paper_trend_df.reset_index(drop=True)
    
    
    
#----------------------------------------------------------------------------------------------------------------#
    
    def get_global_ylim(self, bef_list1, bef_list2, key_inf, key_sup):
        all_vals = []
        for bef in bef_list1 + bef_list2:
            all_vals.extend(bef[key_inf])
            all_vals.extend(bef[key_sup])
        return np.min(all_vals), np.max(all_vals)
    

  
#----------------------------------------------------------------------------------------------------------------#
    
    def comparison_plot(self, calib1, calib2):    
        fig, axs = plt.subplots(3, 2, sharex=True, sharey='row')
        fig.subplots_adjust(wspace=0.30, hspace=0) 
    
        bef_list1 = [calib1.iloc[i:i+7] for i in range(0, len(calib1), 7)]
        bef_list2 = [calib2.iloc[i:i+7] for i in range(0, len(calib2), 7)]
        
        # Calcular os y-limits para cada linha
        ylim_mu0 = self.get_global_ylim(bef_list1, bef_list2, "mu0_inf", "mu0_sup")
        ylim_mu1 = self.get_global_ylim(bef_list1, bef_list2, "mu1_inf", "mu1_sup")
        ylim_mu2 = self.get_global_ylim(bef_list1, bef_list2, "mu2_inf", "mu2_sup")
    
        
        # Plot - Primeira Coluna (bef_list1)
        for i in range(5):
            axs[0, 0].plot(bef_list1[i]["Min_Richness"], bef_list1[i]["mu0"], label=f'{self.MASS_CUT_list[i]:.2}', ls='-.', lw=0.5, marker='o', ms=3)
            axs[0, 0].fill_between(bef_list1[i]["Min_Richness"], bef_list1[i]["mu0_inf"], bef_list1[i]["mu0_sup"], alpha=0.2)
        
            axs[1, 0].plot(bef_list1[i]["Min_Richness"], bef_list1[i]["mu1"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[1, 0].fill_between(bef_list1[i]["Min_Richness"], bef_list1[i]["mu1_inf"], bef_list1[i]["mu1_sup"], alpha=0.2)
        
            axs[2, 0].plot(bef_list1[i]["Min_Richness"], bef_list1[i]["mu2"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[2, 0].fill_between(bef_list1[i]["Min_Richness"], bef_list1[i]["mu2_inf"], bef_list1[i]["mu2_sup"], alpha=0.2)
        
        # Plot - Segunda Coluna (bef_list2)
        for i in range(5):
            axs[0, 1].plot(bef_list2[i]["Min_Richness"], bef_list2[i]["mu0"], label=f'{self.MASS_CUT_list[i]:.2}', ls='-.', lw=0.5, marker='o', ms=3)
            axs[0, 1].fill_between(bef_list2[i]["Min_Richness"], bef_list2[i]["mu0_inf"], bef_list2[i]["mu0_sup"], alpha=0.2)
        
            axs[1, 1].plot(bef_list2[i]["Min_Richness"], bef_list2[i]["mu1"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[1, 1].fill_between(bef_list2[i]["Min_Richness"], bef_list2[i]["mu1_inf"], bef_list2[i]["mu1_sup"], alpha=0.2)
        
            axs[2, 1].plot(bef_list2[i]["Min_Richness"], bef_list2[i]["mu2"], ls='-.', lw=0.5, marker='o', ms=3)
            axs[2, 1].fill_between(bef_list2[i]["Min_Richness"], bef_list2[i]["mu2_inf"], bef_list2[i]["mu2_sup"], alpha=0.2)
        
        # Aplicar os mesmos y-limits por linha
        axs[0, 0].set_ylim(ylim_mu0)
        axs[0, 1].set_ylim(ylim_mu0)
        axs[1, 0].set_ylim(ylim_mu1)
        axs[1, 1].set_ylim(ylim_mu1)
        axs[2, 0].set_ylim(ylim_mu2)
        axs[2, 1].set_ylim(ylim_mu2)
        
        # Rótulos dos eixos
        for row in range(3):
            axs[row, 0].set_ylabel(f'$\mu{row}$')
            for col in range(2):
                axs[row, col].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f'{x:.2f}'))
                axs[row, col].tick_params(axis='y', labelleft=True)  # força números do eixo y
                axs[row, col].tick_params(axis='x')
        
        # Rótulo do eixo x na base
        axs[2, 0].set_xlabel('$\lambda_{c}$')
        axs[2, 1].set_xlabel('$\lambda_{c}$')
        
        axs[0, 0].set_title("Simple log-normal")
        axs[0, 1].set_title("TReND")
        
        # Legenda fora
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
        
        plt.ylabel('$\ln BF$', fontsize=16)
        plt.xlabel(r'$\lambda_{c}$', fontsize=16)
        # plt.title('Bayes factor' )
        
        plt.grid()
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        
        plt.show()
        