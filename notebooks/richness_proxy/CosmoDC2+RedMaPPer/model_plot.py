import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as gridspec
from matplotlib.gridspec import GridSpec
import numpy as np


def Plot_model(model, data_set,  mean):
    
    fig = plt.figure(figsize=(16,4))
#     gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[1, 1], height_ratios=[2.5, 1], wspace=0.4)
    
#     ax1 = fig.add_subplot(gs[0, :], projection='3d')
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax3 = fig.add_subplot(gs[1, 1])
    ax1 = fig.add_subplot(1,3,1, projection='3d')
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)

    #Mu plot config ---------------------------------------
    if mean == True:
        #lnR model, z_mean, y_mean
        halos_mean = bd.bin_meanf(b_z,b_m)[0]
        lnR_mean = np.log(halos_mean["richness"])

        xs = halos_mean["redshift_true"]
        ys = np.log(halos_mean["m200c"])
        zs = model

        lb_model = '<$Ln\lambda$| M, z>'
        lb_binned = '$\ln \lambda_i$ médio'
        m = lnR_mean
        
        if f'{model}' == f'{lnR_mean_ascaso}': 
            fig.suptitle('Modelo de Ascaso')
        elif f'{model}' == f'{lnR_mean_ext_ln1pz}': 
            fig.suptitle('Modelo Estendido (ln(1 + z))')
        elif f'{model}' == f'{lnR_mean_ext_z}': 
            fig.suptitle('Modelo Estendido (z)')
        else: 
            pass
    
    #Std plot config ---------------------------------------
    else:
        halos_mean = bd.bin_meanf(b_z, b_m)[0]
        std_mean = bd.bin_meanf(b_z, b_m)[4]
    
        xs = halos_mean["redshift_true"]
        ys = np.log(halos_mean["m200c"])
        zs = model 
        
        lb_model = '$\sigma_{\ln \lambda}$'
        lb_binned = '$\sigma^{i}$'
        m = std_mean
        
        if f'{model}' == f'{lnR_std_ascaso}': 
            fig.suptitle('Modelo de Ascaso', size=16)
        elif f'{model}' == f'{lnR_std_ext_ln1pz}': 
            fig.suptitle('Modelo Estendido (ln(1 + z))', size=16)
        elif f'{model}' == f'{lnR_std_ext_z}': 
            fig.suptitle('Modelo Estendido (z)', size=16)
        else: 
            pass
        
    #-------------------------------------------------------
    p1 =ax1.scatter(xs, ys, zs, c=zs, cmap='cool')
    ax1.set_xlabel('z')
    ax1.set_ylabel('lnM')
    fig.colorbar(p1, ax=ax1, label= lb_model)

    ax2.scatter(ys, m, c='k', s=2.0, label=lb_binned)
    p2 = ax2.scatter(ys, zs , c= xs, s=2.0, cmap='cool')
    ax2.set_xlabel('lnM')
    ax2.set_ylabel(lb_model)
    fig.colorbar(p2, ax=ax2, label='z')
    ax2.legend()

    ax3.scatter(xs, m, c='k', s=2.0, label=lb_binned)
    p3 = ax3.scatter(xs, zs , c=ys, s=2.0, cmap='cool') 
    ax3.set_xlabel('z')
    ax3.set_ylabel(lb_model)
    fig.colorbar(p3, ax=ax3, label='lnM')
    ax3.legend()         
    plt.show()
