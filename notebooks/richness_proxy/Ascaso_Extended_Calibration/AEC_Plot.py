import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def cplot(model, model_name, data, bins_mean, bins_std, option: str="mean"):
    
    fig = plt.figure(figsize=(12,4))
    # ax1 = fig.add_subplot(1,3,1, projection='3d')
    # ax2 = fig.add_subplot(1,3,2)
    # ax3 = fig.add_subplot(1,3,3)

    ax2 = fig.add_subplot(1,2, 1)
    ax3 = fig.add_subplot(1,2, 2)

    plt.subplots_adjust(wspace=0.3)
    
    xs = data["redshift"]
    ys = np.log10(data["mass"])
    zs = model

    x2 = bins_mean["redshift"]
    y2 = np.log10(bins_mean["mass"])
    
    lnR_mean = np.log(bins_mean["richness"])
    std_mean = bins_std

    
    #Mu plot config --------------------------------------------------------------#
    
    if option == "mean":
        #lnR model, z_mean, y_mean
        
        lb_model = '<$\ln\lambda$| M, z>' #'<$Ln\lambda$| M, z>'
        lb_binned =  '$\hat{\mu_j}$' #'Mean $\ln \lambda_i$'
        m = lnR_mean
        
        if model_name == 'lnR_mean_ascaso': 
            fig.suptitle('Mean') #$<\ln\lambda| M, z>_{LINEAR}$
            
        elif model_name == 'lnR_mean_ascaso_c': 
            fig.suptitle('Mean (with correction)') #$<\ln\lambda| M, z>_{LINEAR}$
            
        elif model_name == 'lnR_mean_ext_ln1pz': 
            fig.suptitle('$<\ln\lambda| M, z>_{QUADRATIC}$')
            
        elif model_name == 'lnR_mean_ext_ln1pz_c': 
            fig.suptitle('$<\ln\lambda| M, z>_{QUADRATIC} (with correction)$')
            
        else: 
            pass
    
   
    
    #Std plot config -------------------------------------------------------------#
    
    elif option == "std":
        
        lb_model = '$\sigma_{\ln \lambda}$'
        lb_binned = '$\hat{\sigma_j}$'
        m = std_mean
        
        if model_name == 'lnR_std_ascaso': 
            fig.suptitle('Standard deviation')
            
        elif model_name == 'lnR_std_ascaso_c': 
            fig.suptitle('$\sigma$ Linear Model')
            
        elif model_name == 'lnR_std_ext_ln1pz': 
            fig.suptitle('$\sigma$ Quadratic Model')
            
        elif model_name == 'lnR_sd_ext_ln1pz_c': 
            fig.suptitle('$\sigma$ Quadratic Model ')
            # (with correction)
        else: 
            pass

    #-----------------------------------------------------------------------------#
    
    else:
        print('option: "mean" or "std"')

    #-----------------------------------------------------------------------------#
    
    # p1 =ax1.scatter(xs, ys, zs, c=zs, cmap='viridis')
    # ax1.set_xlabel('z')
    # ax1.set_ylabel('$log_{10}M$')
    # fig.colorbar(p1, ax=ax1, label= lb_model)

    #-----------------------------------------------------------------------------#

    ax2.scatter(y2, m, c='silver', s=2.0, label=lb_binned)
    p2 = ax2.scatter(ys, zs , c= xs, s=2.0, cmap='viridis')
    
    ax2.set_xlabel('$\log_{10}M$')
    ax2.set_ylabel(lb_model)
    fig.colorbar(p2, ax=ax2, label='z')
    ax2.legend()
    
    
    #-----------------------------------------------------------------------------#
    
    ax3.scatter(x2, m, c='silver', s=2.0, label=lb_binned)
    p3 = ax3.scatter(xs, zs , c=ys, s=2.0, cmap='viridis') 
    
    ax3.set_xlabel('z')
    ax3.set_ylabel(lb_model)
    fig.colorbar(p3, ax=ax3, label='lnM')
    ax3.legend()
    
    
    #-----------------------------------------------------------------------------#
    

