import numpy as np
import matplotlib as mpl

from scipy.interpolate import interp1d as ip_1d
from scipy.interpolate import CubicSpline as c_spline
from scipy.optimize import fsolve

from matplotlib import pyplot as plt


# enable latex
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "times"
})
        

lw_arrow = 1
lw = 1
ms_arrow = 10
ms = 4

label_size = 25


def flow_curve(ax, N, sigma_xx_max, n, arrow_loc = 1):
    """Flow curve around sigma_xy = n."""    
    sigma_xy = np.pi * np.linspace(n-0.5, n+0.5, N)    
    sigma_xx = np.sin(sigma_xy) ** 2
    
    separatrix_color = "black"
    stable_fp_color = "blue"
    instable_fp_color = "red"
    
    ax.plot(sigma_xy, sigma_xx, color = separatrix_color, linewidth = lw)  
    
    arrow_inds = [int(arrow_loc* N/4), int(3 * arrow_loc* N/4)]
    
    for j in range(2):        
        ind = arrow_inds[j]
        ax.annotate('', xy=(sigma_xy[ind+1-j], sigma_xx[ind+1-j]), xytext=(sigma_xy[ind+j], sigma_xx[ind+j]), arrowprops=dict(arrowstyle='-|>', linewidth = lw_arrow, color = separatrix_color, mutation_scale = ms_arrow), zorder = 0)
    
    epsilon = 0.1 * sigma_xx_max
    
    colors_fp = [instable_fp_color, stable_fp_color, instable_fp_color]
    for j in range(-1,2):
        sigma_xx_j = sigma_xx_max * abs(j)
        sigma_xy_j = np.pi * (n + 0.5 * j) 
        ax.plot(sigma_xy_j, sigma_xx_j, color = colors_fp[j + 1], marker ="o", zorder = 3, ms = ms) 
        
        if abs(j) == 0:        
            ax.vlines(sigma_xy_j, color = separatrix_color, ymin = 0, ymax = 2 * sigma_xx_max, linewidth = lw)
            
            ax.annotate('', xy=(sigma_xy_j, (1 + 1/5) * sigma_xx_max), xytext=(sigma_xy_j, (1 + 1/5) * sigma_xx_max + epsilon), arrowprops=dict(arrowstyle='-|>', linewidth = lw_arrow, color = separatrix_color, mutation_scale = ms_arrow), zorder = 0)
            
        else:
            ax.vlines(sigma_xy_j, color = separatrix_color, ymin = 0, ymax = 2 * sigma_xx_max, linewidth = lw, linestyle = "dotted")
            
            ax.annotate('', xy=(sigma_xy_j, sigma_xx_max/2 + abs(j) * epsilon), xytext=(sigma_xy_j, sigma_xx_max/2 + (1 - abs(j)) * epsilon), arrowprops=dict(arrowstyle='-|>', linewidth = lw_arrow, 
                       color = separatrix_color, mutation_scale = ms_arrow), zorder = 0)
            ax.annotate('', xy=(sigma_xy_j, 3 * sigma_xx_max/2), xytext=(sigma_xy_j, 3 * sigma_xx_max/2 + epsilon), arrowprops=dict(arrowstyle='-|>', linewidth = lw_arrow, color = separatrix_color, mutation_scale = ms_arrow), zorder = 0)
            
    return 


def Fig_two_parameter_scaling(sigma_xx_max, sigma_xy_min, sigma_xy_max, N = 1000, sigma_xy_cuttoff = 0.8):
    """Figure to illustrate two parameter scaling theory."""
    
    fig = plt.figure(figsize=(12, 4), layout = "tight")
    a1 =  fig.add_subplot(1,1,1)
    
    x_upper = 2.5 * sigma_xx_max
            
    for n in range(-3,4):
        flow_curve(a1, N, sigma_xx_max, n)
                            
    a1.set_xlabel(r"$\frac{h}{e^2} \sigma_{xy}$", fontsize=label_size)
    a1.set_ylabel(r"$\frac{h}{e^2} \sigma_{xx}$", fontsize=label_size, rotation = 0)
    
    a1.xaxis.set_label_coords(0.975, -0.025)
    a1.yaxis.set_label_coords(0.55, 0.875)
    
    a1.tick_params(direction='out', length=0, width=2, labelsize = label_size, pad = 5)    
      
    # Move left y-axis and bottom x-axis to centre, passing through (0,0) and set linewidth
    a1.spines['left'].set_position("zero")
    a1.spines['bottom'].set_position('zero')
    a1.spines['left'].set_linewidth(lw)
    a1.spines['left'].set_bounds(low = 0, high = x_upper)
    a1.spines['bottom'].set_linewidth(lw)
    
    # Eliminate upper and right axes
    a1.spines['right'].set_color('none')
    a1.spines['top'].set_color('none')      
    
    a1.set_yticks([])
    #a1.set_yticklabels([r"$-1$", r"$1$"])
    
    n_min = int(sigma_xy_min/ np.pi)
    n_max = int(sigma_xy_max/ np.pi)
    
    x_ticks = np.pi * np.array(range(n_min, n_max))
    
    a1.set_xticks(x_ticks)
    a1.set_xticklabels(range(n_min, n_max))
    
    a1.set_xlim([sigma_xy_cuttoff * sigma_xy_min, sigma_xy_cuttoff * sigma_xy_max])
    a1.set_ylim([-0.1, x_upper])
        
    # make arrows
    a1.plot((1), (0), ls="", marker=">", ms=5, color="k", transform=a1.get_yaxis_transform(), clip_on=False)
    a1.plot((0), (x_upper), ls="", marker="^", ms=5, color="k", clip_on=False)
    
    fig.tight_layout()
    fig.savefig("2_param_scaling_DPG.png",  bbox_inches='tight', dpi = 300, transparent=True)

    plt.show()

def main():

    sigma_xx_max = 1
    sigma_xy_min = -4 * np.pi
    sigma_xy_max = 4 * np.pi
    
    N = 1000
    axis_cutoff = 0.85
    
    Fig_two_parameter_scaling(sigma_xx_max, sigma_xy_min, sigma_xy_max, N, axis_cutoff)

    
    
    
    
    
if __name__ == '__main__':
    main()
