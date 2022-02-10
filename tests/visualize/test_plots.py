
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from .settings import *
from ..utilize.auxiliary import *
from ..manage.cellattributes import CellAttributes

# Plots of different proerties:
def rdf_ideal_plotter(ax, bin_edges, nmon, dmon, rdf=True, dmon_large=0, *args, **kwargs):
    """
    plot the probability distribution function (pdf) of the end-to-end vector of an ideal linear chain
    
    inputs:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an ideal linear chainm
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    rdfFactor = 4*np.pi*bin_centers**2
    if rdf == False:
        rdfFactor = 1.0   
    pdf_values = rdfFactor*np.exp((-3.0*bin_centers**2)/(2.0*nmon*dmon**2))
    pdf_values = pdf_values/np.sum(pdf_values) # normalization
    bin_centers = bin_centers / dmon_large # normalized bins
    return ax.plot(bin_centers, pdf_values, **kwargs)

def pdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """
    pdf: probability distribution function of the end-to-end distance vector.
    pdf draws a step plot.

    inputs:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax
    return:
    a plot of (normalized) probability distribution function of the end-to-end distance vector.
    """
    histo_collections = histo_collections / np.sum(histo_collections) # normalization
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    rdfFactor = 4*np.pi*(bin_centers)**2
    histo_collections = np.divide(histo_collections,rdfFactor) 
    histo_collections = histo_collections/np.sum(histo_collections) # normalization
    bin_centers = bin_centers / dmon_large
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)

def rdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """
    rdf: radial distribution function or the probability distribution function of *the magnitude of the end-to-end vector* 
    (or *the end-to-end distance*) in between R and R+dR

    inputs:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax 

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    histo_collections = histo_collections / np.sum(histo_collections) # normalization
    bin_centers = bin_centers/dmon_large # normalized bins
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)

def rdf_real_plotter(ax, bin_edges, rFlory, nmon, dmon,rdf=True, dmon_large=0, *args, **kwargs):
    """
    plot the probability distribution function (pdf) of the end-to-end vector of a real linear chain
    
    inputs:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an ideal linear chainm
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    #rFlory = 1.12 * nmon ** 0.588 * dmon  # chain size or the root-mean-square end-to-end distancs or Flory radius
    rdfFactor = 4*np.pi*(bin_centers)**2
    if rdf == False:
        rdfFactor = 1.0   
    pdf_values = rdfFactor*(bin_centers/rFlory)**0.28*np.exp(-1.206*(bin_centers/rFlory)**2.43) # default coeff = 1.206
    pdf_values = pdf_values/np.sum(pdf_values) # normalized
    bin_centers = bin_centers/dmon_large # normalized bins
    return ax.plot(bin_centers, pdf_values, **kwargs)


def looping_plotter(plotting_df):
    with sns.plotting_context("paper"):
        font_name = "Times New Roman"

        fig, axes = plt.subplots(nrows=2, ncols=1,sharex=True, figsize=(16,9))
        plt.rcParams['font.sans-serif'] = font_name
        plt.rcParams['font.family'] = font_name
        plt.rcParams["font.serif"] = font_name
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'

        colors = sns.color_palette("husl", 3)
        
        nmon = plotting_df.loc[0,'nmon']
        dmon_large = plotting_df.loc[0,'dmon_large']
        dcrowd = plotting_df.loc[0,'dcrowd']

        style = {'marker':'s','markersize':10,'markerfacecolor':'none','linestyle':'-','linewidth':2}
        axes[0].plot('vfrc_crowd','looping_p',data=plotting_df,c=colors[0],**style)
        axes[0].set_title(r'$N={},a_M={},a_c={}$'.format(nmon,dmon_large,dcrowd),fontsize=20,fontname=font_name)
        axes[0].set_xlim(-0.002, 0.302)
        axes[0].set_ylim(-0.002, 0.16)
        axes[0].set_ylabel(r'$P_L(R_{min},R_{max})$',fontsize=18,fontname=font_name)
        axes[0].text(0.002,0.11,r'$\Delta=[0,a_M+a_m]$', fontsize=18, va='bottom',fontname=font_name)
        axes[0].xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[0].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[0].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[0].grid(which='major',ls=':',lw=1)
        axes[0].tick_params(axis='both', which='major', labelsize=16)
        axes[0].tick_params(axis='both', which='minor', labelsize=12)


 
        axes[1].set_xlim(-0.002, 0.302)
        axes[1].set_ylim(0.7, 1.05)
        axes[1].plot('vfrc_crowd','rflory_nor',data=plotting_df,c=colors[2],**style)
        axes[1].set_xlabel(r'$\phi_c$',fontsize=20,fontname=font_name)
        axes[1].set_ylabel(r'$R_F/R_0$',fontsize=18,fontname=font_name)
        axes[1].xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[1].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        axes[1].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
        axes[1].grid(which='major',ls=':',lw=1)
        axes[1].tick_params(axis='both', which='major', labelsize=16)
        axes[1].tick_params(axis='both', which='minor', labelsize=12)
        fig.align_ylabels()
        today = datetime.date.today().strftime('%Y%m%d')
        plt.savefig(today+'_loopP_and_R_vs_vfrc_c.png',dpi=300,format='png')  


def chainsize_plot(df, xcol, save_to=None, fontsize=20, ext='pdf'):
    """
    plot the different chain size types (furthermost distance, radius of gyration, and Flory radius (end-to-end distance)) for a given chain as function of the bulk volume fraction of crowders in varied systems. The difference in systems arises from the difference in the chain size (number of monomers), crowder size, and cylinder size. 
    
    Cuation: 
    This methid is indeed a fine-tuned wrapper to produce high-quality figure, so you probably need to change many of ints internal codes if you want to use in a very differetn style.
    
    Parameters:
    df (pandas dataframe): The dataframe that contains the statisitcal properties of chain in different systems.
    xcol (str): the x-axis variable.
    colors (a list of RGB tuples or matplotlib.colors.Colormap): the length of list choose by the chain property which is represented by hue keyword.
    fontsize (float): Font size for the x and y axis labels.
    save_to (str): the address to which the ensemble files of a group are saved.
    ext (str): extention of output image.
    
    Return:
    Save the plot to a file.
    
    Requirements:
    PipeLine, pandas, matplotlib, seaborn.
    """
    
    xlabels = {"phi_c_bulk":r"$\phi_c$","phi_c_bulk_normalized":r"${a\phi_c}/{a_c}$","phi_c_bulk_eff":r"$\phi_c^{eff}$",
               "phi_c_bulk_eff_normalized":r"${a\phi_c^{eff}}/{a_c}$"}
    xcols = list(xlabels.keys())
    
    ylabels = {"phi_c_bulk":[r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
               "phi_c_bulk_eff":[r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
              "phi_c_bulk_normalized":[r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$'],
              "phi_c_bulk_eff_normalized":[r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$']}
    
    if xcol not in xcols:
        raise ValueError(f"'{xcol}'"
                         " is not a valid x-axis label. Please select one of "
                         f"{xcols} x-axis labels.")
    
    sns.set_context('paper')
    sns.set_style("ticks")
    dcrowds = list(df.dcrowd.sort_values().drop_duplicates())
    dcrowds = list(map(str,dcrowds))
    dcyls = list(df.dcyl.sort_values().drop_duplicates())
    dcyls = list(map(str,dcyls))
    nmons = list(df.nmon.sort_values().drop_duplicates())
    nmons = list(map(str,nmons))
    colors = sns.color_palette(n_colors = len(dcrowds), palette="tab10", desat=1)
    legends_labels = [r'$\frac{a_c}{a}$ (Line color)']+dcrowds+[r'$N$ (Line size)']+nmons+[r'$\frac{D}{a}$ (Marker type)']+dcyls
    fig, axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(16,12))
    
    line1 = sns.lineplot(x=xcol, y="fsd_normalized", hue='dcrowd',style='dcyl', size='nmon', palette=colors, markers=True, markersize=8, data=df,ax=axes[0])
    line2 = sns.lineplot(x=xcol, y="gyr_normalized", hue='dcrowd',style='dcyl', size='nmon', palette=colors, markers=True, markersize=8, data=df, legend=False,ax=axes[1])
    line3 = sns.lineplot(x=xcol, y="rflory_normalized", hue='dcrowd',style='dcyl', size='nmon', palette=colors, markers=True, markersize=8, data=df, legend=False,ax=axes[2])

    for num, axis in enumerate(axes):
        axis.grid(True,ls=':',lw=1)
        axis.tick_params(axis ='both',direction='inout',width=1)
        axis.set_ylabel(ylabels[xcol][num],fontsize=fontsize)
        yticks(axis,(0, 1.0, 0.2), code=True, fontsize=14, decimals=3)
    
    xticks_dict = {"phi_c_bulk":(0.0, round_up_nearest(df['phi_c_bulk'].max(),0.05), 0.05),
                   "phi_c_bulk_normalized":(0.0, round_up_nearest(df['phi_c_bulk_normalized'].max(),0.05), 0.05),
                   "phi_c_bulk_eff":(0.0, round_up_nearest(df['phi_c_bulk_eff'].max(),0.05), 0.05),
                   "phi_c_bulk_eff_normalized":(0.0, round_up_nearest(df['phi_c_bulk_eff_normalized'].max(),0.05), 0.05)}
        
    xticks(axes[2], xticks_dict[xcol], code=True,fontsize=14, decimals=2)
    change_legend_name(line1,legends_labels)
    line3.set_xlabel(xlabels[xcol],fontsize=fontsize)
    if save_to != None:
        picname = save_to+"chainsize-"+xcol+"."+ext
    else :
        picname = "chainsize-"+xcol+"."+ext
    plt.savefig(picname, dpi=300, bbox_inches='tight')


def acf_plot_with_ci(acf, ens_names, group_name, attrs_dict, phi_crds, ylimits = (-0.3, 1.1, 0.2), fontsize = 20, nrows=4, ncols=3, lags= 7000, ext='pdf'):
    """
    plots the auto-correlation functions (AFCs) of a group ofphysical attributes (attrs_dict) with their associated confidence intervals (CIs) in one plot for all the ensembles in a simulation group.
    
    Issues:
    The function does not work if le(ens_names) != ncols*nrows - resloved, the user can set nrows and ncols.
    It does not guarantee that the a phi_crd is correctly associated with an ens_name - resolved, provided that the ens_names and phi_crds are given as ordered lists.
    
    
    Parameters:
    acf (pandas dataframe): the dataset of the ACFs.
    ens_names (list of str): the ordered names of ensembles in a simulations group.
    group_name (str): the name of simulation group
    attrs_dict (dict): the dictionary of attributes for them, the ACF graphs are plotted tegether.
    phi_crds (numpy array): the ordered bulk volume fraction of crowders.
    ylimits (tuple): the lower and upper limits, and the step for the y-ticks.
    fontsize (int): the maximum font size in the plot.
    lags (int): the maximum lag in th AFC plot.
    nrows (int): number of subplot rows
    ncols (int): number of subplot columns
    ext (str): foramt of the output file.
    
    Return:
    Save a plot in ext format into a file.
    
    Requirements:
    PipeLine, Matplotlib, Statsmodels
    """
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16,12), sharey=True, sharex=True)
    mpl.rcParams['font.family'] = "Times New Roman"
    mpl.rcParams['mathtext.default'] = "regular"
    
    def acf_plot(ax, acf_db, time, attr, color):
        
        acf_t = ax.plot(acf_db[attr+'-acf_only'], color = color)
        acf_t = ax.fill_between(time, acf_db[attr+'-acf_ci_lower'], acf_db[attr+'-acf_ci_upper'], color = color, alpha=0.25)
        return acf_t
    
    time = np.arange(lags+1)
    for idx, (ens_name, phi_c, ax) in enumerate(zip(ens_names, phi_crds, axes.flat)):
        ax.axhline(y=0, c='black', ls='--', lw=1)
        acf_ens = acf.loc[acf["ens_name"]==ens_name,:]
        acf_ens.reset_index(inplace=True)
        
        legend_colors = []
        legend_labels = []
        for attr, attr_dict in attrs_dict.items():
            acft_fsd = acf_plot(ax, acf_ens, time, attr, attr_dict['color'])
            legend_colors.append(attr_dict['color'])
            legend_colors.append(mpl.colors.to_rgba(attr_dict["color"],0.25))
            legend_labels.append(attr_dict['symbol'])
            legend_labels.append(attr_dict['symbol']+" CIs")
        
        ax.legend(title=fr'$\phi_c^{{(bulk)}}={phi_c}$', title_fontsize=fontsize-2, framealpha=None, frameon=False)
        yticks(ax, ylimits, code=True, fontsize=fontsize-6, decimals=3)
        xticks(ax, (0, lags, 1000), code=True, fontsize=fontsize-6, decimals=3)
        if idx % 3 == 0:
            ax.set_ylabel(r"$C(\hat{t})$", fontsize=fontsize-2)
        if idx >= 9:
            ax.set_xlabel(r"$\hat{t}=lag\times {\Delta t_{sampling}}/{\tau}$", fontsize=fontsize-2) 
    
    phi_c_patches = color_patcher(legend_colors)
    phi_c_legends = mpl.legend.Legend(axes[0,2], handles=phi_c_patches, labels=legend_labels, title='Physical property', title_fontsize=fontsize-2, fontsize=fontsize-4, framealpha=None, frameon=False, bbox_to_anchor=(1.02,1.02))
    axes[0,2].add_artist(phi_c_legends)
    cell_attrs = CellAttributes(ens_name, geometry='cylindrical', cell_type=None, warning=False)
    fig.suptitle(f"the ACFs of three physical properties of a chain in a system with $N={cell_attrs.nmon}$, $D={cell_attrs.dcyl}$, $a={cell_attrs.dcrowd}$.",fontsize=fontsize+2)
    fname = "acf-confidence_intervals-"+group_name+"."+ext
    fig.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    plt.close()
    

def acf_plot_group(acf, group_names, attr, attr_dict, phi_crds, phi_colors, ylimits = (-0.2, 1.0, 0.2), fontsize = 20, lags=7000, nrows=4, ncols=3, ext='pdf'):
    """
    plots the auto-correlation functions (AFCs) of the physical property "attr" for all the simulation groups.
    
    Issues:
    The function does not work if le(ens_names) != ncols*nrows - resloved, the user can set nrows and ncols.
    It does not guarantee that the a phi_crd is correctly associated with an ens_name - resolved, provided that the ens_names and phi_crds are given as ordered lists.
    
    Parameters:
    acf (pandas dataframe): the dataset of the ACFs.
    group_names (list of str): the ordered names of simulation groups.
    attr (str): the name of physical property for which the ACF is measured.
    phi_crds (numpy array): the ordered bulk volume fraction of crowders.
    phi_colors (list): a list of matplotlib colors
    ylimits (tuple): the lower and upper limits, and the step for the y-ticks.
    fontsize (int): the maximum font size in the plot.
    lags (int): the maximum lag in th AFC plot.
    nrows (int): number of subplot rows
    ncols (int): number of subplot columns
    ext (str): foramt of the output file.
   
    Return:
    Save a plot in ext format into a file.
    
    Requirements:
    PipeLine, Matplotlib, Statsmodels
    """  
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16,12), sharey=True)
    mpl.rcParams['font.family'] = "Times New Roman"
    mpl.rcParams['mathtext.default'] = "regular"
    time = np.arange(lags+1)
    for idx, (group_name, ax) in enumerate(zip(group_names, axes.flat)):
        acf_group = acf.loc[acf['group_name']==group_name,:]
        ens_names = list(acf_group['ens_name'].drop_duplicates())
        ens_names = sorted(ens_names, key = sort_by_int)
        ax.grid(True,ls=':',lw=0.75,c='black')
        ax.tick_params(axis ='both', direction='inout', width=1, labelsize=fontsize-4, color='black')
        ax.axhline(y=0, c='black', ls='--', lw=1)
        for ens_name, phi_c, color in zip(ens_names, phi_crds, phi_colors):
            acf_ens = acf_group.loc[acf_group['ens_name']==ens_name, : ]
            acf_ens.reset_index(inplace=True)
            ax.plot(time, acf_ens[attr+"-acf_only"], color= color)
        yticks(ax, ylimits, code=True, fontsize=fontsize-4, decimals=3)
        xticks(ax,(0, lags, 1000), code=True, fontsize=fontsize-4, decimals=1)
        if idx % 3 == 0:
            ax.set_ylabel(attr_dict['symbol'], fontsize=fontsize-2)
        ax.set_xlabel(r"$\hat{t}=lag\times {\Delta t_{sampling}}/{\tau}$", fontsize=fontsize-2)
        cell_attrs = CellAttributes(ens_name, geometry='cylindrical', cell_type=None, warning=False)
        ax.set_title(fr"({idx+1}) $N={cell_attrs.nmon}$, $D={cell_attrs.dcyl}$, $a={cell_attrs.dcrowd}$",fontsize=fontsize)
    phi_c_patches = color_patcher(phi_colors )
    phi_c_legends = mpl.legend.Legend(axes[0,2], handles=phi_c_patches, labels=list(phi_crds), title=r'$\phi_c^{(bulk)}$', title_fontsize=fontsize-2, fontsize=fontsize-4, framealpha=None, frameon=False, bbox_to_anchor=(1.02,1.02))
    axes[0,2].add_artist(phi_c_legends)
    fname = "acf"+"-"+attr+"."+ext
    fig.tight_layout()
    plt.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close()