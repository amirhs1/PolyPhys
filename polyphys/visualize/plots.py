from datetime import datetime
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from matplotlib import axes
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from polyphys.manage import organizer
from polyphys.visualize import tuning as ptune
from polyphys.manage.parser import SumRule


def rdf_ideal_plotter(
    ax,
    bin_edges,
    nmon,
    dmon,
    rdf=True,
    dmon_large=0,
    **kwargs
):
    """plots the probability distribution function (pdf) of the end-to-end
    vector of an ideal linear chain.

    Parameters:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an
    ideal linear chainm
    dmon_large: size of the pair of monomers at the two ends of the chain.
    This is used for scaling bin centers defined below.

    Return:
    a plot of (normalized) radial probability distribution function of
        the magnitude of the end-to-end distance vector.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  \
        # the centers of the bins being used for plotting
    rdfFactor = 4 * np.pi * bin_centers**2
    if rdf is False:
        rdfFactor = 1.0
    pdf_values = \
        rdfFactor * np.exp((-3.0 * bin_centers**2) / (2.0 * nmon * dmon**2))
    pdf_values = pdf_values / np.sum(pdf_values)
    bin_centers = bin_centers / dmon_large
    return ax.plot(bin_centers, pdf_values, **kwargs)


def pdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """ plots probability distribution function (pdf) of the end-to-end
    distance vector.

    Parameters:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest
    over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large = size of the pair of monomers at the two ends of the chain.
    This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax

    Return:
    a plot of (normalized) probability distribution function of the
    end-to-end distance vector.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    rdfFactor = 4 * np.pi * (bin_centers)**2
    histo_collections = np.divide(histo_collections, rdfFactor)
    histo_collections = histo_collections / np.sum(histo_collections)
    bin_centers = bin_centers / dmon_large
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)


def rdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """plots radial distribution function (rdf) or the probability
    distribution function of the magnitude of the end-to-end vector
    (or *the end-to-end distance*) between R and R+dR.

    Parameters:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest
    over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large: size of the pair of monomers at the two ends of the chain.
    This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax.

    return:
    a plot of (normalized) radial probability distribution function of
    the magnitude of the end-to-end distance vector.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  \
        # the centers of the bins being used for plotting
    histo_collections = histo_collections / np.sum(histo_collections)
    bin_centers = bin_centers / dmon_large
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)


def rdf_real_plotter(
                ax, bin_edges, rFlory, rdf=True, dmon_large=0, **kwargs):
    """plots the probability distribution function (pdf) of the
    end-to-end vector of a real linear chain.

    Parameters:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an ideal
    linear chainm
    dmon_large: size of the pair of monomers at the two ends of the chain.
    This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of
    the magnitude of the end-to-end distance vector.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # rFlory = 1.12 * nmon ** 0.588 * dmon
    # chain size or the root-mean-square end-to-end distancs or Flory radius
    rdfFactor = 4 * np.pi * bin_centers**2
    if rdf is False:
        rdfFactor = 1.0
    pdf_values = \
        (rdfFactor * (bin_centers / rFlory)**0.28
            * np.exp(-1.206 * (bin_centers / rFlory)**2.43))  \
        # default coeff = 1.206
    pdf_values = pdf_values / np.sum(pdf_values)
    bin_centers = bin_centers / dmon_large
    return ax.plot(bin_centers, pdf_values, **kwargs)


def looping_plotter(plotting_df):
    with sns.plotting_context('paper'):
        font_name = "Times New Roman"
        fig, axes = plt.subplots(
            nrows=2, ncols=1, sharex=True, figsize=(16, 9))
        plt.rcParams['font.sans-serif'] = font_name
        plt.rcParams['font.family'] = font_name
        plt.rcParams["font.serif"] = font_name
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'
        colors = sns.color_palette('husl', 3)
        nmon = plotting_df.loc[0, 'nmon']
        dmon_large = plotting_df.loc[0, 'dmon_large']
        dcrowd = plotting_df.loc[0, 'dcrowd']
        style = {
            'marker': 's',
            'markersize': 10,
            'markerfacecolor': 'none',
            'linestyle': '-',
            'linewidth': 2}
        axes[0].plot(
            'vfrc_crowd', 'looping_p', data=plotting_df, c=colors[0], **style)
        axes[0].set_title(r'$N={},a_M={},a_c={}$'.format(
            nmon, dmon_large, dcrowd), fontsize=20, fontname=font_name)
        axes[0].set_xlim(-0.002, 0.302)
        axes[0].set_ylim(-0.002, 0.16)
        axes[0].set_ylabel(
            r'$P_L(R_{min},R_{max})$', fontsize=18, fontname=font_name)
        axes[0].text(
            0.002, 0.11, r'$\Delta=[0,a_M+a_m]$', fontsize=18, va='bottom',
            fontname=font_name)
        axes[0].xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[0].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[0].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[0].grid(which='major', ls=':', lw=1)
        axes[0].tick_params(axis='both', which='major', labelsize=16)
        axes[0].tick_params(axis='both', which='minor', labelsize=12)
        axes[1].set_xlim(-0.002, 0.302)
        axes[1].set_ylim(0.7, 1.05)
        axes[1].plot(
            'vfrc_crowd', 'rflory_nor', data=plotting_df, c=colors[2], **style)
        axes[1].set_xlabel(r'$\phi_c$', fontsize=20, fontname=font_name)
        axes[1].set_ylabel(r'$R_F/R_0$', fontsize=18, fontname=font_name)
        axes[1].xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
        axes[1].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.025))
        axes[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        axes[1].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05))
        axes[1].grid(which='major', ls=':', lw=1)
        axes[1].tick_params(axis='both', which='major', labelsize=16)
        axes[1].tick_params(axis='both', which='minor', labelsize=12)
        fig.align_ylabels()
        today = datetime.date.today().strftime('%Y%m%d')
        plt.savefig(
            today + '_loopP_and_R_vs_vfrc_c.png', dpi=300, format='png')


def chainsize_plot(df, xcol, save_to=None, fontsize=20, ext='pdf'):
    """
    plot the different measures of chain size (furthermost distance, radius
    of gyration, and Flory radius (end-to-end distance)) as function of
    the bulk volume fraction of crowders in varied systems.

    The difference in systems arises from the difference in the chain
    size (number of monomers), crowder size, and cylinder size.

    Cuation:
    This methid is indeed a fine-tuned wrapper to produce high-quality
    figure, so you probably need to change many of ints internal codes
    if you want to use in a very differetn style.

    Parameters:
    df (pandas dataframe): The dataframe that contains the statisitcal
    properties of chain in different systems.
    xcol (str): the x-axis variable.
    colors (a list of RGB tuples or matplotlib.colors.Colormap): the length
    of list choose by the chain property which is represented by hue keyword.
    fontsize (float): Font size for the x and y axis labels.
    save_to (str): address to which the ensemble files of a group are saved.
    ext (str): extention of output image.

    Return:
    Save the plot to a file.
    """
    xlabels = {
        'phi_c_bulk': r"$\phi_c$",
        'phi_c_bulk_normalized': r"${a\phi_c}/{a_c}$",
        'phi_c_bulk_eff': r"$\phi_c^{eff}$",
        'phi_c_bulk_eff_normalized': r"${a\phi_c^{eff}}/{a_c}$"
        }
    xcols = list(xlabels.keys())
    ylabels = {
        'phi_c_bulk': [
            r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',
            r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',
            r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
        'phi_c_bulk_eff': [
            r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',
            r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',
            r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
        'phi_c_bulk_normalized': [
            r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',
            r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',
            r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$'],
        'phi_c_bulk_eff_normalized': [
            r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',
            r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',
            r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$']
    }

    if xcol not in xcols:
        raise ValueError(f"'{xcol}'"
                         " is not a valid x-axis label. Please select one of "
                         f"{xcols} x-axis labels.")
    sns.set_context('paper')
    sns.set_style('ticks')
    dcrowds = list(df.dcrowd.sort_values().drop_duplicates())
    dcrowds = list(map(str, dcrowds))
    dcyls = list(df.dcyl.sort_values().drop_duplicates())
    dcyls = list(map(str, dcyls))
    nmons = list(df.nmon.sort_values().drop_duplicates())
    nmons = list(map(str, nmons))
    colors = sns.color_palette(n_colors=len(dcrowds), palette='tab10', desat=1)
    legends_labels = [r'$\frac{a_c}{a}$ (Line color)'] + dcrowds
    + [r'$N$ (Line size)'] + nmons + [r'$\frac{D}{a}$ (Marker type)'] + dcyls
    _, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(16, 12))
    line1 = sns.lineplot(
        x=xcol, y='fsd_normalized', hue='dcrowd', style='dcyl', size='nmon',
        palette=colors, markers=True, markersize=8, data=df, ax=axes[0])
    _ = sns.lineplot(
        x=xcol, y='gyr_normalized', hue='dcrowd', style='dcyl', size='nmon',
        palette=colors, markers=True, markersize=8, data=df, legend=False,
        ax=axes[1])
    line3 = sns.lineplot(
        x=xcol, y='rflory_normalized', hue='dcrowd', style='dcyl', size='nmon',
        palette=colors, markers=True, markersize=8, data=df, legend=False,
        ax=axes[2])
    for num, axis in enumerate(axes):
        axis.grid(True, ls=':', lw=1)
        axis.tick_params(axis='both', direction='inout', width=1)
        axis.set_ylabel(ylabels[xcol][num], fontsize=fontsize)
        ptune.yticks(axis, (0, 1.0, 0.2), code=True, fontsize=14, decimals=3)
    xticks_dict = {
        'phi_c_bulk': (0.0, ptune.round_up_nearest(
            df['phi_c_bulk'].max(), 0.05),
                0.05),
        'phi_c_bulk_normalized': (0.0, ptune.round_up_nearest(
            df['phi_c_bulk_normalized'].max(), 0.05),
                0.05),
        'phi_c_bulk_eff': (0.0, ptune.round_up_nearest(
            df['phi_c_bulk_eff'].max(), 0.05),
                0.05),
        'phi_c_bulk_eff_normalized': (0.0, ptune.round_up_nearest(
            df['phi_c_bulk_eff_normalized'].max(), 0.05),
                0.05)
    }
    ptune.xticks(
        axes[2], xticks_dict[xcol], code=True, fontsize=14, decimals=2)
    ptune.change_legend_name(line1, legends_labels)
    line3.set_xlabel(xlabels[xcol], fontsize=fontsize)
    if save_to is not None:
        picname = save_to + "chainsize-" + xcol + "." + ext
    else:
        picname = "chainsize-"+xcol+"."+ext
    plt.savefig(picname, dpi=300, bbox_inches='tight')


def acf_plot(
    ax: axes.Axes,
    acf_db: pd.DataFrame,
    time: np.array,
    property_: str,
    color: str
) -> axes.Axes:
    """plot the auto-correlation function (AFCs) with its associated
    confidence intervals (CIs) for the physical `property_`.

    Parameters
    ----------
    ax: matplotlib.axes.Axes
        The Axes in which the data is plotted.
    acf_db: pd.DataFrame
        The dataframe in which the ACF and CIs of the `property_` are
        located.
    time: np.array
        The time data plotted along x axis.
    property_: str
        The physical property of interest.
    color: str
        A color as it recognised by matplotlib.

    Return
    ------
    acf_t: matplotlib.axes.Axes
        the plotted ACF with its assocaited CIs.
    """
    acf_t = ax.plot(acf_db[property_+'-acf-mean'], color=color)
    acf_t = ax.fill_between(
        time,
        acf_db[property_+'-acfLowerCi-mean'],
        acf_db[property_+'-acfUpperCi-mean'],
        color=color,
        alpha=0.25
    )
    return acf_t


def acf_plot_with_ci(
    acf: pd.DataFrame,
    ensembles: List[str],
    space: str,
    properties: Dict[str, Dict[str, str]],
    phi_crds: List[float],
    xlimits: Tuple[float, float, float] = (0, 7000, 1000),
    ylimits: Tuple[float, float, float] = (-0.3, 1.1, 0.2),
    fontsize: int = 20,
    nrows: int = 4,
    ncols: int = 3,
    lags: int = 7000,
    ext: str = 'pdf',
    save_to: str = './'
) -> None:
    """plot the auto-correlation function (AFC) with its associated
    confidence intervals (CIs) of a group of physical `properties` in one
    plot for all the `ensembles` in a simulation `space`, and save it in
    `ext` format into the `save_to` path.

    Parameters
    ----------
    acf: pd.DataFrame
        The dataset of the ACFs.
    ensembles: dict
        The ordered list of ensemble names in the `space`.
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    properties: dict of dict
        A dictionary in which the keys are the name of physical properties for
        which the ACFs are plotted and the values are dictionaries. In each
        internal dictionary, the keys are 'name', 'symbol', 'color', and
        similar characteristics and the values are the spcific values of these
        characteristics for a physical property.
    phi_crds: np.array
        An ordered list of bulk volume fraction of crowders.
    xlimits: tuple, default (0, 7000, 1000),
        The lower and upper limits, and the step for the x-ticks.
    ylimits: tuple, default (-0.3, 1.1, 0.2)
        The lower and upper limits, and the step for the y-ticks.
    fontsize: int, default 20
        The maximum font size in the plot.
    nrows: int, default 4
        The number of subplot rows
    ncols: int, default 3
        The number of subplot columns
    lags: int, default 7000
        The maximum lag in th AFC plot.
    ext: str, default 'pdf'
        The format of the output file.
    save_to : str, default './'
        An/a absolute/relative path of a directory to which outputs are saved.

    Requirements
    ------------
    Matplotlib, Statsmodels
    """
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 12),
        sharey=True,
        sharex=True
    )
    mpl.rcParams['font.family'] = "Times New Roman"
    mpl.rcParams['mathtext.default'] = 'regular'
    time = np.arange(lags+1)
    for idx, (ens_name, phi_c, ax) in \
            enumerate(zip(ensembles, phi_crds, axes.flat)):
        ax.axhline(y=0, c='black', ls='--', lw=1)
        acf_ens = acf.loc[acf['ensemble'] == ens_name, :]
        acf_ens.reset_index(inplace=True)
        legend_colors = []
        legend_labels = []
        for property_, prop_dict in properties.items():
            _ = acf_plot(ax, acf_ens, time, property_, prop_dict['color'])
            legend_colors.append(prop_dict['color'])
            legend_colors.append(mpl.colors.to_rgba(prop_dict['color'], 0.25))
            legend_labels.append(prop_dict['symbol'])
            legend_labels.append(prop_dict['symbol'] + " CIs")
        ax.legend(
            title=fr'$\phi_c^{{(bulk)}}={phi_c}$', title_fontsize=fontsize-2,
            framealpha=None, frameon=False)
        ptune.yticks(ax, ylimits, code=True, fontsize=fontsize-6, decimals=3)
        ptune.xticks(ax, xlimits, code=True, fontsize=fontsize-6, decimals=3)
        if idx % 3 == 0:
            ax.set_ylabel(r"$C(\hat{t})$", fontsize=fontsize-2)
        if idx >= 9:
            ax.set_xlabel(
                r"$\hat{t}=lag\times {\Delta t_{sampling}}/{\tau}$",
                fontsize=fontsize-2)
    phi_c_patches = ptune.color_patcher(legend_colors)
    phi_c_legends = mpl.legend.Legend(
        axes[0, 2],
        handles=phi_c_patches,
        labels=legend_labels,
        title='Physical property',
        title_fontsize=fontsize-2,
        fontsize=fontsize-4,
        framealpha=None,
        frameon=False,
        bbox_to_anchor=(1.02, 1.02)
    )
    axes[0, 2].add_artist(phi_c_legends)
    space_info = SumRule(
        space, geometry='biaxial',
        group='bug',
        lineage='space',
        ispath=False
    )
    fig.suptitle(
        f"the ACFs of three physical properties of a chain in a"
        f"system with $N={space_info.nmon}$, $D={space_info.dcyl}$,"
        f" $a={space_info.dcrowd}$.", fontsize=fontsize+2)
    output = save_to + "acf-confidence_intervals-" + space + "." + ext
    fig.tight_layout()
    plt.savefig(output, bbox_inches='tight')
    plt.close()


def acf_plot_group(
    acf: pd.DataFrame,
    spaces: List[str],
    property_: str,
    property_dict: Dict[str, str],
    phi_crds: List[float],
    phi_colors: sns.palettes._ColorPalette,
    xlimits: Tuple[float, float, float] = (0, 7000, 1000),
    ylimits: Tuple[float, float, float] = (-0.3, 1.1, 0.2),
    fontsize: int = 20,
    nrows: int = 4,
    ncols: int = 3,
    legend_anchor: Tuple[float, float] = (1.25, 1.02),
    lags: int = 7000,
    ext: str = 'pdf',
    save_to: str = './'
) -> None:
    """plot the auto-correlation functions (AFCs) of the physical
    `property_` for all the simulation `spaces`.

    Parameters
    ----------
    acf: pd.DataFrame
        The dataset of the ACFs.
    ensembles: dict
        The ordered list of ensemble names in the `space`.
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    property_: str
        Name of the physical property
    property_dict: dict of str
        A dictionary in which the keys are the name of physical properties for
        which the ACFs are plotted and the values are dictionaries. In each
        internal dictionary, the keys are 'name', 'symbol', 'color', and
        similar characteristics and the values are the spcific values of these
        characteristics for a physical property.
    phi_crds: np.array
        An ordered list of bulk volume fraction of crowders.
    phi_colors: sns.palettes._ColorPalette
        A palette of colors for `phi_crds`.
    xlimits: tuple, default (0, 7000, 1000),
        The lower and upper limits, and the step for the x-ticks.
    ylimits: tuple, default (-0.3, 1.1, 0.2)
        The lower and upper limits, and the step for the y-ticks.
    fontsize: int, default 20
        The maximum font size in the plot.
    nrows: int, default 4
        The number of subplot rows
    ncols: int, default 3
        The number of subplot columns
    legend_anchor: tuple of two floats
        Position at which the legend is placed.
    lags: int, default 7000
        The maximum lag in th AFC plot.
    ext: str, default 'pdf'
        The format of the output file.
    save_to : str, default './'
        An/a absolute/relative path of a directory to which outputs are saved.

    Requirements
    ------------
    Matplotlib
    """
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(16, 12), sharey=True)
    mpl.rcParams['font.family'] = "Times New Roman"
    mpl.rcParams['mathtext.default'] = 'regular'
    time = np.arange(lags+1)
    if nrows * ncols > 1:
        axes_ = axes.flat
    else:
        axes_ = [axes]
    for idx, (space, ax) in enumerate(zip(spaces, axes_)):
        acf_space = acf.loc[acf['space'] == space, :]
        ensembles = list(acf_space['ensemble'].drop_duplicates())
        ensembles = sorted(ensembles, key=organizer.sort_by_alphanumeric)
        ax.grid(True, ls=':', lw=0.75, c='black')
        ax.tick_params(
            axis='both',
            direction='inout',
            width=1,
            labelsize=fontsize-4,
            color='black'
        )
        ax.axhline(y=0, c='black', ls='--', lw=1)
        for ens, color in zip(ensembles, phi_colors):
            acf_ens = acf_space[acf_space['ensemble'] == ens]
            acf_ens.reset_index(inplace=True)
            ax.plot(time, acf_ens[property_ + "-acf-mean"], color=color)
        ptune.yticks(
            ax, ylimits, code=True, fontsize=fontsize-6, decimals=3
        )
        ptune.xticks(
            ax, xlimits, code=True, fontsize=fontsize-6, decimals=3
        )
        if idx % ncols == 0:
            ax.set_ylabel(
                property_dict['symbol'],
                fontsize=fontsize-2
                )
        ax.set_xlabel(
            r"$\hat{t}=lag\times {\Delta t_{sampling}}/{\tau}$",
            fontsize=fontsize - 2)
        space_info = SumRule(
            space, geometry='biaxial',
            group='bug',
            lineage='space',
            ispath=False
        )
        ax.set_title(
            fr"({idx+1}) $N={space_info.nmon}$, $D={space_info.dcyl}$"
            fr", $a={space_info.dcrowd}$",
            fontsize=fontsize
        )
    if nrows == 1 & ncols == 1:
        leg_ax = axes
    elif nrows == 1:
        leg_ax = axes[ncols-1]
    elif ncols == 1:
        leg_ax = axes[0]
    else:
        leg_ax = axes[0, ncols-1]
    phi_c_patches = ptune.color_patcher(phi_colors)
    phi_c_legends = mpl.legend.Legend(
        leg_ax,
        handles=phi_c_patches,
        labels=list(phi_crds),
        title=r'$\phi_c^{(bulk)}$',
        title_fontsize=fontsize-2,
        fontsize=fontsize-4,
        framealpha=None,
        frameon=False,
        bbox_to_anchor=(1.12, 1.02)
    )
    leg_ax.add_artist(phi_c_legends)
    output = save_to + 'acf' + "-" + property_ + "." + ext
    fig.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches='tight')
    plt.close()
