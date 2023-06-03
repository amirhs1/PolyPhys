from datetime import datetime
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)
import numpy as np
import pandas as pd
from matplotlib import axes
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from ..manage import organizer
from ..visualize import tuner as ptuner
from ..analyze import correlations
from ..analyze import measurer
from ..manage.parser import (
    SumRuleCyl, TransFociCyl, TransFociCub, HnsCub, HnsCyl)


PROJECT_DETAILS = {
    'SumRuleCyl': {
        'group': 'bug',
        'divisor': 0.025,
        'geometry': 'cylindrical',
        'geometry_name': 'cylindrical_confinement',
        'chain_name': 'homogeneous_linear',
        'topology': 'linear',
        'parser': SumRuleCyl,
        'space_pat': 'N*D*ac*',
        'hierarchy': 'N*',
        'species': ['Mon', 'Crd'],
        'directions': ['r', 'z'],
        'cross_section': ['xy', 'xz', 'yz'],
        'edge_directions': ['x', 'y', 'z'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon', 'dcyl',
                       'dcrowd', 'phi_c_bulk'],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'rfloryTMon', 'shapeTMon', 'transSizeMon'],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['space', 'ensemble_long', 'ensemble', 'nmon',
                             'dcyl', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round', 'size_ratio'],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'fsdMon-mean',
                             'fsdMon-var', 'fsdMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'rfloryMon-mean',
                             'rfloryMon-var', 'rfloryMon-sem',
                             'shapeMon-mean', 'shapeMon-var', 'shapeMon-sem',
                             'transSizeMon-mean', 'transSizeMon-var',
                             'transSizeMon-sem']
    },
    'TransFociCyl': {
        'group': 'bug',
        'divisor': 0.025,
        'geometry': 'cylindrical',
        'geometry_name': 'cylindrical_confinement',
        'chain_name': 'heterogeneous_ring',
        'topology': 'ring',
        'parser': TransFociCyl,
        'space_pat': 'ns*nl*al*D*ac*',
        'hierarchy': 'eps*',
        'species': ['Mon', 'Foci', 'Crd'],
        'directions': ['r', 'z'],
        'cross_section': ['xy', 'xz', 'yz'],
        'edge_directions': ['x', 'y', 'z'],
        'space_hierarchy': 'ns*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon_small',
                       'nmon_large', 'dmon_large', 'dcyl', 'dcrowd',
                       'phi_c_bulk'],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'shapeTMon'],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'dcyl',
                             'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'fsdMon-mean',
                             'fsdMon-var', 'fsdMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem',
                             'transSizeMon-mean', 'transSizeMon-var',
                             'transSizeMon-sem']
    },
    'TransFociCub': {
        'group': 'bug',
        'divisor': 0.025,
        'geometry': 'cubic',
        'geometry_name': 'free_space',
        'chain_name': 'heterogeneous_ring',
        'topology': 'ring',
        'parser': TransFociCub,
        'space_pat': 'ns*nl*al*ac*',
        'hierarchy': 'al*',
        'species': ['Mon', 'Foci', 'Crd'],
        'directions': ['r'],
        'cross_section': ['xy', 'xz', 'yz'],
        'edge_directions': ['x', 'y', 'z'],
        'space_hierarchy': 'ns*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon_small',
                       'nmon_large', 'dmon_large', 'dcrowd', 'phi_c_bulk'],
        'time_varying_props': ['asphericityTMon', 'gyrTMon', 'shapeTMon'],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space',
                             'dmon_large', 'nmon_large', 'nmon_small',
                             'dcrowd', 'phi_c_bulk', 'phi_c_bulk_round'],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem']
    },
    'HnsCub': {
        'group': 'nucleoid',
        'divisor': 0.025,
        'geometry': 'cubic',
        'geometry_name': 'free_space',
        'chain_name': 'semiflexible_ring',
        'topology': 'ring',
        'parser': HnsCub,
        'space_pat': 'N*epshm*nh*ac*',
        'hierarchy': 'N*',
        'species': ['Mon', 'Hns', 'Crd'],
        'directions': ['r'],
        'cross_section': ['xy', 'xz', 'yz'],
        'edge_directions': ['x', 'y', 'z'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon',
                       'eps_hm', 'nhns', 'dcrowd', 'phi_c_bulk',
                       'rho_hns_bulk'],
        'time_varying_props': ['asphericityTMon', 'gyrTMon', 'shapeTMon'],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'nmon',
                             'eps_hm', 'nhns', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round', 'size_ratio', 'rho_hns_bulk'],
        'equil_properties': ['asphericityMon-mean', 'gyrMon-mean',
                             'shapeMon-mean']
    },
    'HnsCyl': {
        'group': 'nucleoid',
        'divisor': 0.01,
        'geometry': 'cylindrical',
        'geometry_name': 'free_space',
        'chain_name': 'semiflexible_ring',
        'topology': 'ring',
        'parser': HnsCyl,
        'space_pat': 'N*D*nh*ac*',
        'hierarchy': 'N*',
        'species': ['Mon', 'Hns', 'Crd'],
        'directions': ['r', 'z'],
        'cross_section': ['xy', 'xz', 'yz'],
        'edge_directions': ['x', 'y', 'z'],
        'space_hierarchy': 'N*',
        'attributes': ['space', 'ensemble_long', 'ensemble', 'nmon', 'dcyl',
                       'nhns', 'dcrowd', 'phi_c_bulk', 'rho_hns_bulk'],
        'time_varying_props': ['asphericityTMon', 'fsdTMon', 'gyrTMon',
                               'shapeTMon'],
        'equil_measures': [np.mean, np.var, measurer.sem],
        'equil_attributes': ['ensemble_long', 'ensemble', 'space', 'nmon',
                             'dcyl', 'nhns', 'dcrowd', 'phi_c_bulk',
                             'phi_c_bulk_round', 'rho_hns_bulk'],
        'equil_properties': ['asphericityMon-mean', 'asphericityMon-var',
                             'asphericityMon-sem', 'fsdMon-mean',
                             'fsdMon-var', 'fsdMon-sem', 'gyrMon-mean',
                             'gyrMon-var', 'gyrMon-sem', 'shapeMon-mean',
                             'shapeMon-var', 'shapeMon-sem',
                             'transSizeMon-mean', 'transSizeMon-var',
                             'transSizeMon-sem']
    }
}

SIZE_MEASURES_LABELS = {
    'rfloryTMon': {
        'name': 'Flory radius',
        'symbol': r'$R_F/a_m$',
        'symbol-norm': r'$R_F/\langle R_F \rangle$',
        'pdf': r'$f(R_F)$',
        'acf': r'$c_{R_FR_F}$'
    },
    'gyrTMon': {
        'name': 'radius of gyration',
        'symbol': r'$R_g/a_m$',
        'symbol-norm': r'$R_g/\langle R_g \rangle$',
        'pdf': r'$f(R_g)$',
        'acf': r'$c_{R_gR_g}$'
    },
    'fsdTMon': {
        'name': 'furthermost distance',
        'symbol': r'$L/a_m$',
        'symbol-norm': r'$L/\langle L \rangle$',
        'pdf': r'$f(L)$',
        'acf': r'$c_{LL}$'
    },
    'transSizeTMon': {
        'name': 'furthermost distance',
        'symbol': r'$L_{\perp}/a_m$',
        'symbol-norm': r'$L_{\perp}/\langle L_{\perp} \rangle$',
        'pdf': r'$f(L_{\perp})$',
        'acf': r'$c_{L_{\perp}L_{\perp}}$'
    },
    'asphericityTMon': {
        'name': 'asphericity',
        'symbol': r'$\Delta$',
        'acf': r'$c_{\Delta\Delta}$'
    },
    'shapeTMon': {
        'name': 'shape parameter',
        'symbol': r'$S$',
        'acf': r'$c_{SS}$'
    },
    'gyrTMon-zscoreNorm': {
        'name': 'radius of gyration',
        'symbol': r'$(R_g-\langle R_g \rangle)/\sigma_{R_g}$',
    },
    'asphericityTMon-zscoreNorm': {
        'name': 'asphericity',
        'symbol': r'$(\Delta-\langle \Delta \rangle)/\sigma_{\Delta}$',
    },
    'shapeTMon-zscoreNorm': {
        'name': 'shape parameter',
        'symbol': r'$(S-\langle S \rangle)/\sigma_{S}$',
    }
}

ATTRIBUTE_LABELS = {
    'rfloryMon-norm': r'$\langle R_F\rangle/\langle R_{F,0}\rangle$',
    'gyrMon-norm': r'$\langle R_g\rangle/\langle R_{g,0}\rangle$',
    'fsdMon-norm': r'$\langle L\rangle/\langle L_{0}\rangle$',
    'transSizeMon-norm':
        r'$\langle L_{\perp}\rangle/\langle L_{\perp,0}\rangle$',
    'asphericityMon-norm':
        r'$\langle \Delta\rangle/\langle \Delta_{0}\rangle$',
    'shapeMon-norm': r'$\langle S\rangle/\langle S_{0}\rangle$',
    'gyrMon-mean': r'$\langle R_g\rangle$',
    'asphericityMon-mean': r'$\langle \Delta\rangle$',
    'shapeMon-mean': r'$\langle S\rangle$',
    'rPhi': r'$\phi(r)$',
    'zPhi': r'$\phi(|z|)$',
    'rRho': r'$\rho(r)$',
    'zRho': r'$\rho(|z|)$',
    'rPhi-norm': r'$\phi(r)$',
    'rPhi-norm-Mon': r'$\phi_m(r)$',
    'rPhi-norm-Foci': r'$\phi_M(r)$',
    'rPhi-norm-Hns': r'$\phi_n(r)$',
    'rPhi-norm-Crd': r'$\phi_c(r)$',
    'rPhi-norm-Sum': r'$\sum_i\phi_i(r)/a_i$',
    'rPhi-norm-Sum_constant': r'$\sum_c\phi_c(\infty)/a_c$',
    'zPhi-norm': r'$\phi(r)$',
    'zPhi-norm-Mon': r'$\phi_m(|z|)$',
    'zPhi-norm-Foci': r'$\phi_M(|z|)$',
    'zPhi-norm-Hns': r'$\phi_n(|z|)$',
    'zPhi-norm-Crd': r'$\phi_c(|z|)$',
    'zPhi-norm-Sum': r'$\sum_i\phi_i(|z|)/a_i$',
    'zPhi-norm-Sum_constant': r'$\sum_c\phi_c(\infty)/a_c$',
    'rRho-norm': r'$\rho(r)$',
    'rRho-norm-Mon': r'${{\rho_m(r)}}/{{\rho_m(0)}}$',
    'rRho-norm-Foci': r'${{\rho_M(r)}}/{{\rho_M(0)}}$',
    'rRho-norm-Hns': r'${{\rho_n(r)}}/{{\rho_n(0)}}$',
    'rRho-norm-Crd': r'${{\rho_c(r)}}/{{\rho_c(\infty)}}$',
    'rRho-norm-Sum':
        r'${(\sum_i\rho_i(r)a_i^2)}/{(\sum_i\rho_i(\infty)a_i^2)}$',
    'rRho-norm-Sum_constant':
        r'${(\sum_c\rho_c(\infty)a_c^2)}//{(\sum_i\rho_i(\infty)a_i^2)}$',
    'zRho-norm': r'$\rho(|z|)$',
    'zRho-norm-Mon': r'${{\rho_m(|z|)}}/{{\rho_m(0)}}$',
    'zRho-norm-Foci': r'${{\rho_M(|z|)}}/{{\rho_M(0)}}$',
    'zRho-norm-Crd': r'${{\rho_c(|z|)}}/{{\rho_c(\infty)}}$',
    'zRho-norm-Sum':
        r'${(\sum_i\rho_i(|z|)a_i^2)}/{(\sum_i\rho_i(\infty)a_i^2)}$',
    'zRho-norm-Sum_constant':
        r'${(\sum_c\rho_c(\infty)a_c^2)}//{(\sum_i\rho_i(\infty)a_i^2)}$',
    'rPhi-norm-old-Mon': r'${{\phi_m(r)}}/{{\phi_m(0)}}$',
    'rPhi-norm-old-Foci': r'${{\phi_M(r)}}/{{\phi_M(0)}}$',
    'rPhi-norm-old-Hns': r'${{\phi_n(r)}}/{{\phi_n(0)}}$',
    'rPhi-norm-old-Crd': r'${{\phi_c(r)}}/{{\phi_c(\infty)}}$',
    'rPhi-norm-old-Sum':
        r'${(\sum_i\phi_i(r)/a_i)}/{(\sum_i\phi_i(\infty)/a_i)}$',
    'rPhi-norm-old-Sum_constant':
        r'${(\sum_c\phi_c(\infty)/a_c)}/{(\sum_i\phi_i(\infty)/a_i)}$',
    'zPhi-norm-old-Mon': r'${{\phi_m(|z|)}}/{{\phi_m(0)}}$',
    'zPhi-norm-old-Foci': r'${{\phi_M(|z|)}}/{{\phi_M(0)}}$',
    'zPhi-norm-old-Crd': r'${{\phi_c(|z|)}}/{{\phi_c(\infty)}}$',
    'zPhi-norm-old-Sum':
        r'${(\sum_i\phi_i(|z|)/a_i)}/{(\sum_i\phi_i(\infty)/a_i)}$',
    'zPhi-norm-old-Sum_constant':
        r'${(\sum_c\phi_c(\infty)/a_c)}/{(\sum_i\phi_i(\infty)/a_i)}$',
    'rcut_bonding': r'$r_{{thres}}/a_m$',
    'rcut_bonding-norm': r'$r_{{thres}}/(a_M+a_c)$',
    'bondsHistFoci-bin_center': r'$x_d$',
    'bondsHistFoci-bin_center-average': r'$\langle x_d\rangle$',
    'bondsHistFoci-norm': r'$p$',
    'bondsHistFoci-norm-full_name': r'$p(x_d)$',
    'bondsHistFoci-norm-bonding_threshold': r'$p(r_{thres};x_d)$',
    'clustersHistFoci-bin_center': r'$x_c$',
    'clustersHistFoci-bin_center-average': r'$\langle x_c\rangle$',
    'clustersHistFoci-norm': r'$p$',
    'clustersHistFoci-norm-full_name': r'$p(x_c)$',
    'pairDistHistFoci-legend_title':
        r'$(n_i,n_j,\Delta n_{ij}/\langle l_m \rangle)$',
    'pairDistHistFoci-legend_title-index_diff':
        r'$(n_i,n_j,\Delta n_{ij})$',
    'pairDistHistFoci': r'$\mathcal{H}(n_i,n_j,\Delta n_{ij})$',
    'pairDistRdfFoci-legend_title':
        r'$(n_i,n_j,\Delta n_{ij}/\langle l_m \rangle)$',
    'pairDistRdfFoci-legend_title-index_diff':
        r'$(n_i,n_j,\Delta n_{ij})$',
    'pairDistRdfFoci': r'$f(r)$',
    'pairDistRdfFoci-full_name': r'$f(n_i,n_j,\Delta n_{ij};r)$',
    'pairDistRdfGenDistAvg': r'$f(r)$',
    'pairDistRdfGenDistAvg-full_name': r'$f(\Delta n;r)$',
    'pairDistTFoci': r'$r(t)$',
    "size_ratio": "",
    "space": "",
    "ensemble": "",
    "phi_c_bulk_round": r"$\phi_c$",
    "time": r"$t$",
    "t_idx_norm": r"$t/t^{{max}}$",
    "t_index": r"${t}/{\delta t}$",
    "lag_time": r"${t}_{lag}$",
    "lag_index": r"${t}_{lag}$",
    "dmon_large": "${{a_M}}/{{a_m}}$",
    "dcrowd": "${{a_c}}/{{a_m}}$",
    "dcyl": "${{D}}/{{a_m}}$",
    "lcyl": "${{L_{{cyl}}}}/{{a_m}}$",
    "nmon_small": "$N_m$",
    "nmon_large": "$N_M$",
    'phi_c_bulk_norm': r"${a\phi_c}/{a_c}$",
    'genomic_distance':  r"$\Delta n/\langle l\rangle$",
    'index_difference':  r"$\Delta n$",
    'loopLengthHistMon-norm': r'$f$',
    'loopLengthHistMon-norm-full_name': r'$f(\Delta n)$',
    'bin_center': '$r/a_m$',  # Pair distance between two monomers in their PDF
    'bin_center-norm-r-dmon_large': '${{r}}/{{a_M}}$',
    'bin_center-norm': '${{r}}/{{r_{max}}}$',
    'bin_center-r': '$r$',
    'bin_center-norm-r-dcyl': '${{2r}}/{{D}}$',
    'bin_center-fsd_mean-r': '${{2r}}/{{D}}$',
    'bin_center-recentered-norm-r-cyl': '${{(2r-a^{shift})}}/{{D}}$',
    'bin_center-norm-r-cub': '${{r}}/{{r_{max}}}$',
    'bin_center-dcrowd-r': '${{2r}}/{{a_c}}$',
    'bin_center-dcrowd-recentered-r': '${{(2r-a^{shift})}}/{{a_c}}$',
    'bin_center-z': '$z$',
    'bin_center-norm-z': '${{|z|}}/{{z_{max}}}$',
    'bin_center-dcrowd-z': '${{|z|}}/{{z_{max}}}$',
    'bin_center-fsd_mean-z': r'${{2|z|}}/{{\langle L(\phi_c)\rangle}}$',
    'bin_center-dcrowd-recentered-z': '${{|z|}}/{{z_{max}}}$',
    'bin_center-theta': r'$\theta$',
    'bin_center-dcrowd-theta': r'${{\theta}}/{{\pi}}$',
    'bin_center-norm-theta': r'${{\theta}}/{{\pi}}$',
    'geometry': 'System geometry',
    'cub': 'Free space',
    'cyl': 'Cylindrical Confinement',
    'nmon': '$N$',
    'nhns': r'$N_{{hns}}$',
    'c_hns_bulk': r'$c_{hns} (\mu M)$',
    'c_hns_bulk-old': r'$[HNS] (\mu M)$',
    'phi_c_rescaled': r'$\\frac{{a_m\phi_c}}{{a_c}}$',
    'confinement_rate': r'$\kappa=\\frac{{a_c}}{{D-a_c}}$',
    'confinement_rate_r': r'$\kappa=\\frac{{D-a_c}}{{a_c}}$',
    'dep_energy_max':
        r'$\mathcal{F}_{dep}=\phi_c^{(bulk)}[1+{{3a_m}}{{2a_c}}]$',
    'int_energy_max':
        r'$\mathcal{F}_{int}=\\frac{{Na\phi_c^{(bulk)}}}{{a_c}}$' +
        r'$[3a_ma_c + \\frac{{3}}{{2}}]$',
    'species': 'Type',
    'nBoundHnsPatch-mean': r'$\langle N_{hns}^{bound}\rangle$',
    'nFreeHnsPatch-mean': r'$\langle N_{hns}^{free}\rangle$',
    'nEngagedHnsPatch-mean': r'$\langle N_{hns}^{engaged}\rangle$',
    'nFreeHnsCore-mean': r'$\langle N_{hns}^{free}\rangle$',
    'nBridgeHnsCore-mean': r'$\langle N_{hns}^{bridge}\rangle$',
    'nDangleHnsCore-mean': r'$\langle N_{hns}^{dangle}\rangle$',
    'nBoundHnsCore-mean': r'$\langle N_{hns}^{bound}\rangle$',
    'nCisMon-mean': r'$\langle N^{cis}\rangle$',
    'nTransMon-mean': r'$\langle N^{trans}\rangle$',
    'bondLengthMon-mean': r'$\langle l \rangle$',
    'nBoundHnsPatch-norm':
        r'$\frac{{\langle N_{hns}^{bound}\rangle}}{{2N_h}}$',
    'nFreeHnsPatch-norm':
        r'$\frac{{\langle N_{hns}^{free}\rangle}}{{2N_h}}$',
    'nHnsPatch-total':
        r'$\frac{{\langle N_{hns}^{total}\rangle}}{{2N_h}}$',
    'nEngagedHnsPatch-norm':
        r'$\frac{{\langle N_{hns}^{engaged}\rangle}}{{2N_h}}$',
    'nFreeHnsCore-norm':
        r'$\frac{{\langle N_{hns}^{free}\rangle}}{{N_h}}$',
    'nBridgeHnsCore-norm-old':
        r'$\frac{{\langle N_{hns}^{bridge}\rangle}}{{N_h}}$',
    'nDangleHnsCore-norm-old':
        r'$\frac{{\langle N_{hns}^{dangle}\rangle}}{{N_h}}$',
    'nHnsCore-total':
        r'$\frac{{\langle N_{hns}^{total}\rangle}}{{N_h}}$',
    'nBoundHnsCore-norm-old':
        r'$\frac{{\langle N_{hns}^{bound}\rangle}}{{N_h}}$',
    'nCisHnsCore-norm-old':
        r'$\frac{{\langle N_{hns}^{cis}\rangle}}{{N_h}}$',
    'nTransHnsCore-norm-old':
        r'$\frac{{\langle N_{hns}^{trans}\rangle}}{{N_h}}$',
    'nEngagedHnsPatch-norm-old':
        r'$\frac{{\langle N_{hns}^{engaged}\rangle}}{{2N_h}}$',
    'nBridgeHnsCore-norm': r'$\theta_{bridge}$',
    'nBridgeHnsCore-norm-rel':
        r'${\theta_{bridge}}/{\theta_{bound}}$',
    'nDangleHnsCore-norm-': r'$\theta_{dangle}$',
    'nDangleHnsCore-norm-rel':
        r'${\theta_{dangle}}/{\theta_{bound}}$',
    'nBoundHnsCore-norm': r'$\theta_{bound}$',
    'nCisHnsCore-norm': r'$\theta_{cis}$',
    'nCisHnsCore-norm-rel':
        r'$\theta_{cis}/\theta_{bridge}$',
    'nTransHnsCore-norm': r'$\theta_{trans}$',
    'nTransHnsCore-norm-rel':
        r'$\theta_{trans}/\theta_{bridge}$',
    'patch_binding_style': 'Binding State',
    'rcut_binding': r'$r_{{thres}}/a_m$',
    'cis_threshold': r'$\Delta n_{{thres}}/\langle l\rangle$',
    'binding_probability_rcut': r'$\theta$',
    'binding_style': 'Binding state',
    'bridging_probability_rcut': r'$\theta$',
    'bridging_probability_nhns_norm':
        r'$\theta_i/\theta_{bridge}$',
    'bridging_style': 'Bridging state',
}

# https://www.heavy.ai/blog/12-color-palettes-for-telling-better-stories-with-your-data
DUTCH_FEILD_COLORS = ["#e60049", "#0bb4ff", "#50e991", "#9b19f5",
                      "#ffa300", "#dc0ab4", "#b3d4ff", "#e6d800", "#00bfa0"]
# https://sashamaps.net/docs/resources/20-colors/
ACCES_COLORS = ['#e6194B', '#3cb44b', '#ffe119', '#911eb4', '#4363d8',
                '#f58231', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
                '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000',
                '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9']
# https://en.wikipedia.org/wiki/Help:Distinguishable_colors
DIST_COLORS = ['#F0A3FF', '#0075DC', '#993F00', '#4C005C', '#191919',
               '#005C31', '#2BCE48', '#FFCC99', '#808080', '#94FFB5',
               '#8F7C00', '#9DCC00', '#C20088', '#003380', '#FFA405',
               '#FFA8BB', '#426600', '#FF0010', '#5EF1F2', '#00998F',
               '#E0FF66', '#740AFF', '#990000', '#FFFF80', '#FFE100',
               '#FF5005']
# min luminocity 1%, max luminocity: 80%: https://mokole.com/palette.html
OTHER_COLORS = ['#2f4f4f', '#8b4513', '#228b22', '#00008b', '#ff0000',
                '#ffff00', '#00ff00', '#00ffff', '#ff00ff', '#1e90ff',
                '#eee8aa', '#ff69b4']

AMIRHSI_COLORS_ORDER = ["#e60049", '#ffe119', '#eee8aa', "#9b19f5", "#ffa300",
                        '#ff00ff', '#8b4513', "#00008b", '#00ffff', '#228b22',
                        '#00ff00', '#2f4f4f', '#a9a9a9']
AMIRHSI_COLORS = ["#00008b", "#e60049", "#9b19f5", '#228b22', 'gray',
                  '#ffe119', '#ff00ff', '#00ffff', "#ffa300", '#00ff00',
                  '#8b4513', '#000000']


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
    ideal linear chain
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
    linear chain
    dmon_large: size of the pair of monomers at the two ends of the chain.
    This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of
    the magnitude of the end-to-end distance vector.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # rFlory = 1.12 * nmon ** 0.588 * dmon
    # chain size or the root-mean-square end-to-end distances or Flory radius
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
    with sns.set_theme('paper'):
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


def p_vexc(
    axis: axes.Axes,
    vexc_data: pd.DataFrame,
    xname: str,
    yname: str,
    hue: str,
    palette: dict,
    **kwargs,
):
    """plot the excluded-volume data `yname` as function of `xname` for
    different values of the crowder size 'dcrowd'.

    Parameters
    ----------
    axis: matplotlib.axes
        The axis over which the curves are plotted
    vexc_data: pd.DataFrame
        The excluded-volume data.
    xname: str
        The name of column used as x data.
    yname: str
        The name of column used as y data.
    hue: str:
        The name of column used to color the data.
    palette: dict
        A dictionary in which the keys are the crowder sizes and the values
        are their associated colors.
    kwargs: dict
       The **kwargs in the `matplotlib.plot` used to customize the plot.
    Requirements
    ------------
    matplotlib, seaborn, pandas
    """
    _legend_titles = {
        'dcrowd': r'${a_c}/{a}$',
        'vexc_model': 'Model'
    }
    mpl.rcParams['font.family'] = "Times New Roman"
    mpl.rcParams['mathtext.default'] = 'regular'
    plt.rcParams.update({"text.usetex": True})
    line_plot = sns.lineplot(
        x=xname,
        y=yname,
        ax=axis,
        hue=hue,
        palette=palette,
        data=vexc_data,
        legend='full',
        **kwargs
    )
    handlers, labels = line_plot.get_legend_handles_labels()
    line_plot.legend_.remove()
    line_plot.legend(
        handlers,
        labels,
        title=_legend_titles[hue],
        ncol=2,
        framealpha=0.4
    )
    axis.grid(True, ls=':', lw=1)
    axis.tick_params(axis='both', direction='inout', width=1.25)
    axis.axhline(0, color='black', lw=1, ls='--')


def p_vexc_models(
    vexc_df: pd.DataFrame,
    dcrowds: List[float],
    phi_c_limit: float = 0.4,
    color_map: str = 'coolwarm',
    save_to: str = './'
) -> None:
    """plot the excluded-volume of monomers (scaled by their value in the
    absence of crowders)  as function of the bulk volume fraction of crowders
    and rescaled bulk volume fraction of crowder for different excluded-volume
    models.

    Parameters
    ----------
    vexc_df: pd.DataFrame
        The excluded-volume data.
    dcrowds: list of float
        List of crowder sizes.
    phi_c_limit: float, default 0.4
        The upper limit on the bulk volume fraction of crowders in the
        `vexc_df`.
    color_map: str, default 'coolwarm'
        Name of the color map used for curves with different crowder sizes.
    save_to : str, default './'
        An/a absolute/relative path of a directory to which outputs are saved.
    """
    _xnames = ['phi_c_bulk', 'phi_c_bulk_scaled']
    _vexc_styles = {
        'AO': {
            'model_fullname': 'Asakura-Oosawa interactions',
            'line': {
                'ls': '-',
                'lw': 1.75,
                'alpha': None
            },
        },
        'LJ': {
            'model_fullname': 'Lennard-Jones interactions',
            'line': {
                'ls': '--',
                'lw': 1.75,
                'alpha': None
            },
        },
        'Edwards': {
            'model_fullname': 'Edwards-type interactions',
            'line': {
                'ls': '-.',
                'lw': 1.75,
                'alpha': None
            },
        }
    }
    _label_dict = {
        "phi_c_bulk": r"$\phi_c$",
        "phi_c_bulk_scaled": r"${a\phi_c}/{a_c}$",
        "vexc_scaled": r"$v_{exc}/v_{0}$",
    }
    _mpl_settings = {
        'legend.fontsize': 10,
        'legend.title_fontsize': 12,
        'font.family': 'Times New Roman',
        'font.weight': 'bold',
        'mathtext.default': 'regular',
        'text.usetex': True
    }
    plt.rcParams.update(_mpl_settings)
    vexc_df = vexc_df.loc[vexc_df['phi_c_bulk'] <= phi_c_limit, :]
    vexc_df = vexc_df.loc[vexc_df['dcrowd'].isin(dcrowds), :]
    dcrowds_color = sns.color_palette(color_map, len(dcrowds))
    dcrowds_dict = dict(zip(dcrowds, dcrowds_color))
    fig, axes = plt.subplots(
        nrows=len(_vexc_styles.keys()),
        ncols=2,
        sharey=True,
        sharex='col',
        figsize=(16, 9)
    )
    for row_idx, vexc_model in enumerate(_vexc_styles.keys()):
        row_axes = axes[row_idx, :]
        for xname, axis in zip(_xnames, row_axes):
            model_cond = vexc_df['vexc_model'] == vexc_model
            vexc_per_model = vexc_df.loc[model_cond, :]
            p_vexc(
                axis,
                vexc_per_model,
                xname,
                'vexc_scaled',
                'dcrowd',
                dcrowds_dict,
                **_vexc_styles[vexc_model]['line']
            )
            axis.set_title(
                _vexc_styles[vexc_model]['model_fullname'],
                fontsize=16
            )
            axis.set_xlabel(_label_dict[xname], fontsize=16)
            axis.set_ylabel(_label_dict['vexc_scaled'], fontsize=16)
    fig.tight_layout()
    output = save_to + f"vExc-Models-phicLimit{str(phi_c_limit)}.pdf"
    plt.savefig(output, dpi=200)


def p_vexc_dcrowds(
    vexc_df,
    dcrowds,
    phi_c_limit: float = 0.4,
    color_map: str = 'flare',
    fontsize: int = 13,
    save_to: str = './'
):
    vexc_df = vexc_df.loc[vexc_df['phi_c_bulk'] <= phi_c_limit, :]
    vexc_df = vexc_df.loc[vexc_df['dcrowd'].isin(dcrowds), :]
    dcrowds_color = sns.color_palette(color_map, len(dcrowds))
    dcrowds_dict = dict(zip(dcrowds, dcrowds_color))
    plt.rcParams.update({
        'legend.fontsize': fontsize-3,
        'legend.title_fontsize': fontsize-3,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'font.family': 'Times New Roman',
        'font.weight': 'bold',
        'mathtext.default': 'regular',
        'text.usetex': True
    })
    fig, axes = plt.subplots(
        nrows=len(dcrowds),
        ncols=1,
        sharex=True,
        figsize=(8, 6)
    )
    for row_idx, (dcrowd, axis) in enumerate(zip(dcrowds, axes)):
        dcrowd_cond = vexc_df['dcrowd'] == dcrowd
        vexc_per_dcrowd = vexc_df.loc[dcrowd_cond, :]
        line_plot = sns.lineplot(
            x='phi_c_bulk',
            y='vexc_scaled',
            style='vexc_model',
            ax=axis,
            data=vexc_per_dcrowd,
            color=dcrowds_dict[dcrowd],
            dashes=[(1, 0), (2, 2), (4, 1)]
        )
        handlers, labels = line_plot.get_legend_handles_labels()
        line_plot.legend_.remove()
        if row_idx == 0:
            line_plot.legend(
                handlers,
                labels,
                title='Model',
                ncol=3,
                framealpha=0.6,
                loc='lower right'
            )
        color_patches = ptuner.color_patcher([dcrowds_dict[dcrowd]])
        color_legend = mpl.legend.Legend(
            line_plot,
            handles=color_patches,
            labels=[rf'${{a_c}}/{{a}}={dcrowd}$'],
            title=None,
            loc='lower left',
            framealpha=0.6
        )
        line_plot.add_artist(color_legend)
        axis.grid(True, ls=':', lw=1)
        axis.tick_params(axis='both', direction='inout', width=1.25)
        axis.axhline(0, color='black', lw=1, ls='--')
        axis.set_xlabel(r"$\phi_c$", fontsize=fontsize)
        axis.set_ylabel(r"$v_{exc}/v_{0}$", fontsize=fontsize)
    fig.tight_layout()
    output = save_to + f"vExc-Dcrowds-phicLimit{str(phi_c_limit)}.pdf"
    plt.savefig(output, dpi=200)


def p_acf(
    ax: axes.Axes,
    acf_db: pd.DataFrame,
    time: np.array,
    property_: str,
    color: str,
    fbetween_label: str = None,
    **acf_kwargs,
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
    time: np.ndarray
        The time data plotted along the x-axis.
    property_: str
        The physical property of interest.
    color: str
        A color as it recognised by matplotlib.
    fbetween_label: str, default None
        Label for the fill_between curve (the confidence interval).
    **acf_kwargs :
        Keyword arguments passed to `ax.plot` method below.

    Return
    ------
    acf_t: matplotlib.axes.Axes
        the plotted ACF with its associated CIs.
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "mathtext.default": "regular"
    })
    acf_t = ax.plot(acf_db[property_+'-acf-mean'], color=color, **acf_kwargs)
    if fbetween_label is not None:
        acf_t = ax.fill_between(
            time,
            acf_db[property_+'-acfLowerCi-mean'],
            acf_db[property_+'-acfUpperCi-mean'],
            color=color,
            alpha=0.25,
            label='CI'
        )
    else:
        acf_t = ax.fill_between(
            time,
            acf_db[property_+'-acfLowerCi-mean'],
            acf_db[property_+'-acfUpperCi-mean'],
            color=color,
            alpha=0.25,
        )
    return acf_t


def p_acf_with_ci_space(
    space_acf: pd.DataFrame,
    space: str,
    space_title: str,
    project_name: str,
    properties: Dict[str, Dict[str, str]],
    xlimits: Tuple[float, float, float] = (0, 7000, 1000),
    ylimits: Tuple[float, float, float] = (-0.3, 1.1, 0.2),
    fontsize: int = 20,
    ncols: int = 3,
    lags: int = 7000,
    ext: str = 'pdf',
    save_to: str = './',
    **save_kwargs
):
    """Plots the auto-correlation function (AFC) with its associated
    confidence intervals (CIs) of a group of physical `properties` in one
    plot for all the `ensembles` in a simulation `space`, and save it in
    `ext` format into the `save_to` path.

    Issues:
    1. Currently, this plot works well if we have 12=4*3 values for each
    property. Look at seaborn col_wrap parameter.

    Parameters
    ----------
    space_acf: pd.DataFrame
        The dataset of the ACFs.
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    space_title: str
        The formatted space name for the plot title.
    project_name: str
        The name of a project.
    properties: dict of dict
        A dictionary in which the keys are the name of physical properties for
        which the ACFs are plotted and the values are dictionaries. In each
        internal dictionary, the keys are 'name', 'symbol', 'color', and
        the like, and the values are the specific values of these
        characteristics for a physical property.
    xlimits: tuple, default (0, 7000, 1000),
        The lower and upper limits, and the step for the x-ticks.
    ylimits: tuple, default (-0.3, 1.1, 0.2)
        The lower and upper limits, and the step for the y-ticks.
    fontsize: int, default 20
        The maximum font size in the plot.
    ncols: int, default 3
        The number of subplot columns
    lags: int, default 7000
        The maximum lag in th AFC plot; use to slice the data and change the
        range of plotted data.
    ext: str, default 'pdf'
        The format of the output file.
    save_to : str, default './'
        An/a absolute/relative path of a directory to which outputs are saved.
    save_kwrags:
        Keywords arguments pass to the plt.savefig method.

    Requirements
    ------------
    Matplotlib, Statsmodels
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "mathtext.default": "regular"
    })
    ensembles = list(space_acf['ensemble'].drop_duplicates())
    ensembles = sorted(ensembles, key=organizer.sort_by_alphanumeric)
    time = np.arange(lags+1)
    n_ens = len(ensembles)
    nrows = int(np.ceil(n_ens / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 12),
        sharey=True,
        sharex=False
    )
    for idx, (ens_name, ax) in enumerate(zip(ensembles, axes.flat[:n_ens])):
        ax.axhline(y=0, c='black', ls='--', lw=1)
        ens_acf = space_acf.loc[space_acf['ensemble'] == ens_name, :]
        ens_acf.reset_index(inplace=True)
        # all the rows have the same 'phi_c_bulk_round' value:
        phi_c = ens_acf.loc[0, 'phi_c_bulk_round']
        legend_colors = []
        legend_labels = []
        for property_, prop_dict in properties.items():
            _ = p_acf(
                ax,
                ens_acf.loc[:lags],
                time,
                property_,
                prop_dict['color']
            )
            legend_colors.append(prop_dict['color'])
            legend_colors.append(mpl.colors.to_rgba(prop_dict['color'], 0.25))
            legend_labels.append(prop_dict['symbol'])
            legend_labels.append(prop_dict['symbol'] + ": CIs")
        ax_title = fr'({idx+1}) $\phi_c={phi_c}$'
        ax.set_title(ax_title, fontsize=fontsize-5)
        ptuner.yticks(ax, ylimits, code=True, fontsize=fontsize-6, decimals=3)
        ptuner.xticks(ax, xlimits, code=True, fontsize=fontsize-6, decimals=3)
        if idx % ncols == 0:
            ax.set_ylabel(
                r"$C_{\mathcal{X}\mathcal{X}}(\hat{t}_{lags})$",
                fontsize=fontsize-3
            )
        ax.set_xlabel(
            r"$\hat{t}_{lags}$",
            fontsize=fontsize-2
        )
    phi_c_patches = ptuner.color_patcher(legend_colors)
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
    fig.suptitle(
        r"the ACFs of different physical properties " + space_title,
        fontsize=fontsize+2
    )
    if ncols * nrows > n_ens:
        axes_to_del = ncols * nrows - n_ens
        for col_idx in range(1, axes_to_del + 1):
            fig.delaxes(axes[nrows-1, -1 * col_idx])
    output = "-".join(["acf-confidence_intervals", project_name, space])
    output = output + '-lags' + str(lags) + "." + ext
    fig.tight_layout()
    fig.savefig(save_to + output, bbox_inches='tight', **save_kwargs)
    plt.close()


def p_acf_fit_curve_space(
    space_acf: pd.DataFrame,
    space: str,
    space_title: str,
    project_name: str,
    space_stamps: pd.DataFrame,
    func_name: str,
    property_: str,
    properties: Dict[str, str],
    x_type: str = 'index',
    round_to: float = 0.025,
    fontsize: int = 20,
    nrows: int = 4,
    ncols: int = 3,
    ext: str = 'pdf',
    save_to: str = './',
    **save_kwargs
) -> None:
    """plot the auto-correlation function (AFC) with its associated
    confidence intervals (CIs) of a physical `property_` with its associated
    fitting function 'func' in one plot for all the `ensembles` in a simulation
    `space`, and save it in `ext` format into the `save_to` path.

    Issues:
    1. Currently, this plot works well if we have 12=4*3 values for each
    property. Look at seaborn col_wrap parameter.

    Parameters
    ----------
    space_acf: pd.DataFrame
        The dataset of the ACFs.
    space: str
        The name of the simulation `space` to which `ensembles` belong.
    space_title: str
        The formatted space name for the plot title.
    project_name: str
        The name of a project.
    space_stamps: pd.DataFrame         that contains the physical attributes, equilibrium
        properties, and the fitting parameters of `func` for each
        'ensemble' in the 'space'
    func_name: str
        Function fitted to the ACF data of `property_`.
    property_: str
        The physical property of interest.
    properties: dict of dict
        A dictionary in which the keys are the name of physical properties for
        which the ACFs are plotted and the values are dictionaries. In each
        internal dictionary, the keys are 'name', 'symbol', 'color', and
        the like and the values are the specific values of these
        characteristics for a physical property.
    x_type: {'index', 'time'}, default 'index'
        Either use the 'index' of the data set as x variable or use the real
    round_to: float, default 0.025
        The float number to which the bulk volume fraction of crowders is
        rounded.
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
    ext: str, default 'pdf'
        The format of the output file.
    save_to : str, default './'
        An/a absolute/relative path of a directory to which outputs are saved.
    save_kwrags:
        Keywords arguments pass to the plt.savefig method.

    Requirements
    ------------
    Matplotlib, Statsmodels
    """
    fit_info = {
        'mono_unit_exp': {
            'func': correlations.mono_unit_exp,
            'name': r'$\exp[-\omega t^{\alpha}]$',
            'params': ['omega', 'alpha'],
            'p_labels': r'$(\omega,\alpha)$'
        },
        'mono_unit_exp_tau': {
            'func': correlations.mono_unit_exp_tau,
            'name': r'$\exp[-\frac{t}{\tau}^{\alpha}]$',
            'params': ['tau', 'alpha'],
            'p_labels': r'$(\omega,\alpha)$'
        },
        'mono_exp_tau_res': {
            'func': correlations.mono_exp_tau_res,
            'name': r'$A\exp[-\frac{t}{\tau}^{\alpha}]+B$',
            'params': ['tau', 'alpha', 'amp', 'residue'],
            'p_labels': r'$(\omega,\alpha,A,B)$'
        },
        'mono_exp_res': {
            'func': correlations.mono_exp_res,
            'name': r'$A\exp[-\omega t^{\alpha}]+B$',
            'params': ['omega', 'alpha', 'amp', 'residue'],
            'p_labels': r'$(\omega,\alpha,A,B)$'
        },
        'mono_exp': {
            'func': correlations.mono_exp,
            'name': r'$A\exp[-\omega t^{\alpha}]$',
            'params': ['omega', 'alpha', 'amp'],
            'p_labels': r'$(\omega,\alpha,A)$'
        },
        'mono_exp_tau': {
            'func': correlations.mono_exp_tau,
            'name': r'$A\exp[-\frac{t}{\tau}^{\alpha}]$',
            'params': ['tau', 'alpha', 'amp'],
            'p_labels': r'$(\omega,\alpha,A)$'
        }
    }
    func = fit_info[func_name]['func']
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 12),
        sharey=True,
        sharex=True
    )
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "mathtext.default": "regular"
    })
    cols = ['phi_c_bulk']
    cols.extend(fit_info[func_name]['params'])
    for idx, (ax, (_, row)) in enumerate(
        zip(
            axes.flat, space_stamps[cols].iterrows()
        )
    ):
        phi_c = row['phi_c_bulk_round']
        params = row[fit_info[func_name]['params']].values
        ens_df = space_acf.loc[space_acf['phi_c_bulk_round'] == phi_c]
        ens_df.reset_index(inplace=True)
        phi_c = np.round(np.floor(phi_c/round_to)*round_to, 3)
        if x_type == 'time':
            x = ens_df['time']
        else:
            x = ens_df.index.values
        y = func(x, *params)
        x_axis = ens_df.index.values
        ax.axhline(y=0, c='black', ls='--', lw=1)
        ax.axhline(
            y=1/np.e,
            c='black',
            ls='-',
            lw=1,
            alpha=0.5,
            label=r'$ACF(t)=e^{-1}$'
        )
        _ = p_acf(
            ax,
            ens_df,
            x_axis,
            property_,
            'steelblue',
            fbetween_label='CI',
            label='data'
        )
        params = np.around(params, 3)
        params = [str(param) for param in list(params)]
        params = ','.join(params)
        ax.plot(
            x_axis,
            y,
            ls='--',
            lw=1.5,
            color='firebrick',
            label='fit')
        if idx == (ncols - 1):
            ax.legend(fontsize=fontsize-7, frameon=False, ncol=2)
        ax.tick_params(
            axis='both',
            direction='inout',
            width=1,
            labelsize=fontsize-8,
            color='black'
        )
        ax_title = fr'({idx+1}) $\phi_c={phi_c}$: '
        ax_title = ax_title + fit_info[func_name]['p_labels']
        ax_title = ax_title + r'$=$' + fr'$({params})$'
        ax.set_title(ax_title, fontsize=fontsize-5)
        if idx % ncols == 0:
            ax.set_ylabel(
                r"$C_{xx}(\hat{t}_{lags})$",
                fontsize=fontsize-3
            )
        ax.set_xlabel(
            r"$\hat(t)_{lags}$",
            fontsize=fontsize-2
        )
    title = (
        'Exponential fit (' + fit_info[func_name]['name']
        + ') to the auto-correlation function (ACF) of the '
        + properties[property_]['name']
        + space_title
    )
    fig.suptitle(title, fontsize=fontsize - 2)
    fig.tight_layout()
    output = '-'.join(
        ['acf', 'fitCurve', project_name, func_name, property_, space]
    )
    output = output + '.' + ext
    fig.savefig(save_to + output, bbox_inches='tight', **save_kwargs)
    plt.close()


def p_acf_allInOne_project(
    acf: pd.DataFrame,
    project: str,
    space_titles: Dict[str, str],
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
    save_to: str = './',
    **save_kwrags
) -> None:
    """plot the auto-correlation functions (AFCs) of the physical `property_`
    for all the simulation `spaces` in a project, using `attr_values` and
    `attr_colors` for setting axes' colors.

    Parameters
    ----------
    acf: pd.DataFrame
        The dataset of the ACFs.
    project: str
        The name of a project.
    space_titles: dict of str and str
        The dict of space and space titles for sub-plots.
    property_: str
        Name of the physical property
    property_dict: dict of str
        A dictionary in which the keys are the name of physical properties for
        which the ACFs are plotted and the values are dictionaries. In each
        internal dictionary, the keys are 'name', 'symbol', 'color', and
        the like and the values are the specific values of these
        characteristics for a physical property.
    phi_crds: np.ndarray
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
    **save_kwrags:
        Keywords arguments pass to the plt.savefig method.

    Requirements
    ------------
    Matplotlib, Pandas
    """
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(16, 12), sharey=True)
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "mathtext.default": "regular"
    })
    time = np.arange(lags + 1)
    if nrows * ncols > 1:
        axes_ = axes.flat
    else:
        axes_ = [axes]
    for idx, ((space, space_title), ax) in enumerate(
            zip(space_titles.items(), axes_)):
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
        ptuner.yticks(
            ax, ylimits, code=True, fontsize=fontsize-6, decimals=3
        )
        ptuner.xticks(
            ax, xlimits, code=True, fontsize=fontsize-6, decimals=3
        )
        if idx % ncols == 0:
            ax.set_ylabel(
                property_dict['symbol'],
                fontsize=fontsize-2
                )
        ax.set_xlabel(
            r"$\hat{t}_{lags}$",
            fontsize=fontsize - 2)

        ax.set_title(
            fr"({idx+1})" + space_title,
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
    phi_c_patches = ptuner.color_patcher(phi_colors)
    phi_c_legends = mpl.legend.Legend(
        leg_ax,
        handles=phi_c_patches,
        labels=list(phi_crds),
        title=r'$\phi_c$',
        title_fontsize=fontsize-2,
        fontsize=fontsize-4,
        framealpha=None,
        frameon=False,
        bbox_to_anchor=legend_anchor
    )
    leg_ax.add_artist(phi_c_legends)
    output = "-".join(['acf', project, property_]) + "." + ext
    fig.tight_layout(pad=1)
    plt.savefig(save_to + output, bbox_inches='tight', **save_kwrags)
    plt.close()


def p_hist_allInOne_project(
    data: pd.DataFrame,
    project_title: str,
    x_prop,
    y_prop,
    hue_attr: str,
    col_attr: str,
    prop_labels: Dict[str, Dict[str, str]],
    attr_labels: Dict[str, str],
    kind: Optional[str] = 'point',
    col_wrap: Optional[int] = 3,
    height: Optional[float] = 3,
    aspect: Optional[float] = 1.618,  # golden ratio
    color_palette: Optional[str] = 'rocket_r',
    plot_context: Optional[str] = 'paper',
    font_scale: Optional[int] = 2,
    rc_params: Optional[Dict[str, Union[str, float, int]]] = {
        'mathtext.default': 'regular',
        'text.usetex': True
        },
    axes_style: Optional[str] = 'ticks',
    font_family: Optional[str] = 'Times New Roman',
    facet_kws: Optional[dict] = {
        'sharey': False, 'sharex': False, 'legend_out': True
        },
    fig_title_kws: Optional[dict] = {'x': 0.5, 'y': 1.0},
    loc: str = 'center right',
    legend_style: str = 'full',
    move_legend_kws: Optional[dict] = None,
    **kwargs
) -> sns.FacetGrid:
    """Plots the time series of a physical `property_` in a given `space` for
    all the values of a given `hue_attr` (different colors) and all the values
    of a given `col_attr` (different subplots).

    Parameters
    ----------
    data: pd.DataFrame
        The input dataset
    project_title: str
        The formatted name of the project.
    x_prop: str
        The time-series index-like property.
    y_prop: str
        The time-series-like property
    hue_attr: str
        The categorical-like attribute used for coloring curves.
    col_attr: str
        The categorical-like attribute used for setting number of sub-plots.
    kind: str, default 'point'.
        The kind of catplot.
    prop_labels: dict of str
        The formatted details of the `y_prop`.
    attr_labels:
        The formatted titles of the `prop_x`, `hue_attr`, and `col_attr`.
    col_wrap: int, default 2
        The number of sub-figure columns
    height: float, default 3
        The height of each sub-figure.
    aspect: float, default 1.618
        The ratio of width to height.
    color_palette: str, default 'rocket_r'
        The color palette`
    plot_context: str, default 'paper'
        Set seaborn plotting context
    font_scale: float, default 2
        The scaling factor for fonts.
    rc_params: dict, default {
        'mathtext.default': 'regular',
        'text.usetex': True
        }
        The rc parameters.
    axes_style: str, default 'ticks',
        The seaborn axes style.
    font_family: str, default 'Times New Roman'
        The seaborn font family
    figsize: tuple, default (12,6)
        The size of the figure.
    leg_ncol: int, default 1
        The number of columns in the legend.
    facet_kws: dict, default {
        'sharey': False, 'sharex': False, 'legend_out': True
        }
        kwargs passed to FaceGrid.
    fig_title_kws: dict, default {'x': 0.5, 'y': 1.0}
        The kwargs passed to `FacetGrid.fig.suptitle`
    loc: str, default center right
        Location of the legend
    legend_style: str
        One of `seaborn.catplot` legend keywords.
    move_legend_kws: dict, default None
        kwargs pass to seaborn.move_legend
    kwargs:
        kwargs passed to relplot.

    Return
    ------
    tseries_grid: sns.FacetGrid
        A FacetGrid plot.

    Requirements
    ------------
    Matplotlib, Seaborn, Pandas
    """
    sns.set_theme(
        context=plot_context,
        style=axes_style,
        palette=color_palette,
        font=font_family,
        font_scale=font_scale,
        rc=rc_params
    )
    tseries_grid = sns.catplot(
        data=data,
        x=x_prop,
        y=y_prop,
        col=col_attr,
        hue=hue_attr,
        kind=kind,
        ls="-",
        height=height,
        alpha=0.7,
        aspect=aspect,
        palette=color_palette,
        col_wrap=col_wrap,
        legend=legend_style,
        facet_kws=facet_kws,
        **kwargs
    )
    tseries_grid.set_axis_labels(
        attr_labels[x_prop],
        prop_labels[y_prop]['symbol']
    )
    tseries_grid.set_titles(attr_labels[col_attr] + r"$={col_name}$")
    if legend_style is not False:
        tseries_grid.legend.set_title(attr_labels[hue_attr])
        sns.move_legend(tseries_grid, loc=loc, **move_legend_kws)
    tseries_grid.fig.suptitle(project_title, **fig_title_kws)
    tseries_grid.tight_layout(w_pad=1)
    return tseries_grid


def p_pairDist_allInOne_project_lineStyle(
    data: pd.DataFrame,
    project_title: str,
    x_prop,
    y_prop,
    hue_attr: str,
    style_attr: str,
    style_order: list,
    col_attr: str,
    prop_labels: Dict[str, Dict[str, str]],
    attr_labels: Dict[str, str],
    legend_labels,
    col_wrap: int = 3,
    height: float = 3,
    aspect: float = 1.618,  # golden ratio
    color_palette: str = 'rocket_r',
    plot_context: Optional[str] = 'paper',
    font_scale: int = 2,
    rc_params: Optional[Dict[str, Union[str, float, int]]] = {
        'mathtext.default': 'regular',
        'text.usetex': True
        },
    axes_style: Optional[str] = 'ticks',
    font_family: Optional[str] = 'Times New Roman',
    facet_kws: dict = {
        'sharey': False, 'sharex': False, 'legend_out': True
        },
    fig_title_kws: dict = {'x': 0.5, 'y': 1.0},
    alpha: float = 1.0,
    loc: str = 'center right',
    move_legend_kws: Optional[dict] = None,
    **kwargs
) -> sns.FacetGrid:
    """Plots the time series of a physical `property_` in a given `space` for
    all the values of a given `hue_attr` (different colors) and all the values
    of a given `col_attr` (different subplots).

    Parameters
    ----------
    data: pd.DataFrame
        The input dataset
    project_title: str
        The formatted name of the project.
    x_prop: str
        The time-series index-like property.
    y_prop: str
        The time-series-like property
    hue_attr: str
        The categorical-like attribute used for coloring curves.
    style_attr: str
        The categorical-like attribute used for line size.
    style_attr: list
        The order of categorical-like attribute used for line size.
    col_attr: str
        The categorical-like attribute used for setting number of sub-plots.
    prop_labels: dict of str
        The formatted details of the `y_prop`.
    attr_labels:
        The formatted titles of the `prop_x`, `hue_attr`, and `col_attr`.
    col_wrap: int, default 2
        The number of sub-figure columns
    height: float, default 3
        The height of each sub-figure.
    aspect: float, default 1.618
        The ratio of width to height.
    color_palette: str, default 'rocket_r'
        The color palette`
    plot_context: str, default 'paper'
        Set seaborn plotting context
    font_scale: float, default 2
        The scaling factor for fonts.
    rc_params: dict, default {
                               'mathtext.default': 'regular',
                               'text.usetex': True
                               }
        The rc parameters.
    axes_style: str, default 'ticks',
        The seaborn axes style.
    font_family: str, default 'Times New Roman'
        The seaborn font family
    figsize: tuple, default (12,6)
        The size of the figure.
    leg_ncol: int, default 1
        The number of columns in the legend.
    facet_kws: dict, default {
        'sharey': False, 'sharex': False, 'legend_out': True
        }
        kwargs passed to FaceGrid.
    fig_title_kws: dict, default {'x': 0.5, 'y': 1.0}
        The kwargs passed to `FacetGrid.fig.suptitle`
    loc: str, default center right
        Location of the legend
    move_legend_kws: dict, default None
        kwargs pass to seaborn.move_legend
    alpha: float, default 1,
        The transparency of curves
    kwargs:
        kwargs passed to relplot.

    Return
    ------
    tseries_grid: sns.FacetGrid
        A FacetGrid plot.

    Requirements
    ------------
    Matplotlib, Seaborn, Pandas
    """
    sns.set_theme(
        context=plot_context,
        style=axes_style,
        palette=color_palette,
        font=font_family,
        font_scale=font_scale,
        rc=rc_params
    )
    tseries_grid = sns.relplot(
        data=data,
        x=x_prop,
        y=y_prop,
        col=col_attr,
        hue=hue_attr,
        style=style_attr,
        style_order=style_order,
        kind='line',
        height=height,
        aspect=aspect,
        palette=color_palette,
        col_wrap=col_wrap,
        facet_kws=facet_kws,
        alpha=alpha,
        errorbar=None,
        **kwargs
    )
    for idx, new_label in enumerate(legend_labels):
        tseries_grid._legend.texts[idx].set_text(new_label)
    tseries_grid.set_ylabels(prop_labels[hue_attr]['symbol'])
    tseries_grid.set_xlabels(attr_labels[x_prop])
    tseries_grid.set_titles(attr_labels[col_attr] + r"$={col_name}$")
    tseries_grid.fig.suptitle(project_title, **fig_title_kws)
    tseries_grid.tight_layout(w_pad=1)
    sns.move_legend(tseries_grid, loc=loc, **move_legend_kws)
    return tseries_grid


def p_pairDist_allInOne_project_colStyle(
    data: pd.DataFrame,
    project_title: str,
    x_prop,
    y_prop,
    hue_attr: str,
    col_attr: str,
    col_order: list,
    row_attr: str,
    row_order: list,
    prop_labels: Dict[str, Dict[str, str]],
    attr_labels: Dict[str, str],
    legend_labels,
    height: float = 3,
    aspect: float = 1.618,  # golden ratio
    color_palette: str = 'rocket_r',
    plot_context: Optional[str] = 'paper',
    font_scale: int = 2,
    rc_params: Optional[Dict[str, Union[str, float, int]]] = {
        'mathtext.default': 'regular',
        'text.usetex': True
        },
    axes_style: Optional[str] = 'ticks',
    font_family: Optional[str] = 'Times New Roman',
    facet_kws: dict = {
        'sharey': False, 'sharex': False, 'legend_out': True
        },
    loc: str = 'center right',
    move_legend_kws: Optional[dict] = None,
    alpha: float = 1.0,
    fig_title_kws: dict = {'x': 0.5, 'y': 1.0},
    **kwargs
) -> sns.FacetGrid:
    """Plots the time series of a physical `property_` in a given `space` for
    all the values of a given `hue_attr` (different colors) and all the values
    of a given `col_attr` (different subplots).

    Parameters
    ----------
    data: pd.DataFrame
        The input dataset
    project_title: str
        The formatted name of the project.
    x_prop: str
        The time-series index-like property.
    y_prop: str
        The time-series-like property
    hue_attr: str
        The categorical-like attribute used for coloring curves.
    col_attr: str
        The categorical-like attribute used for setting number of cols.
    col_order: list
        The order of categorical-like attribute used for setting number of
        cols.
    row_attr: str
        The categorical-like attribute used for setting number of rows.
    row_order: list
        The order of categorical-like attribute used for setting number of
        rows.
    prop_labels: dict of str
        The formatted details of the `y_prop`.
    attr_labels:
        The formatted titles of the `prop_x`, `hue_attr`, and `col_attr`.
    height: float, default 3
        The height of each sub-figure.
    aspect: float, default 1.618
        The ratio of width to height.
    color_palette: str, default 'rocket_r'
        The color palette`
    plot_context: str, default 'paper'
        Set seaborn plotting context
    font_scale: float, default 2
        The scaling factor for fonts.
    rc_params: dict, default {
                               'mathtext.default': 'regular',
                               'text.usetex': True
                               }
        The rc parameters.
    axes_style: str, default 'ticks',
        The seaborn axes style.
    font_family: str, default 'Times New Roman'
        The seaborn font family
    facet_kws: dict, default {
        'sharey': False, 'sharex': False, 'legend_out': True
        }
        kwargs passed to FaceGrid.
    loc: str, default 'center right
        Location of the legend
    fig_title_kws: dict, default {'x': 0.5, 'y': 1.0}
        The kwargs passed to `FacetGrid.fig.suptitle`
    move_legend_kws: dict, default None
        kwargs pass to seaborn.move_legend
    alpha: float, default 1,
        The transparency of curves
    kwargs:
        kwargs passed to relplot.

    Return
    ------
    tseries_grid: sns.FacetGrid
        A FacetGrid plot.
    """
    sns.set_theme(
        context=plot_context,
        style=axes_style,
        palette=color_palette,
        font=font_family,
        font_scale=font_scale,
        rc=rc_params
    )
    tseries_grid = sns.relplot(
        data=data,
        x=x_prop,
        y=y_prop,
        col=col_attr,
        col_order=col_order,
        row=row_attr,
        row_order=row_order,
        hue=hue_attr,
        kind='line',
        height=height,
        aspect=aspect,
        palette=color_palette,
        facet_kws=facet_kws,
        alpha=alpha,
        errorbar=None,
        **kwargs
    )
    for idx, new_label in enumerate(legend_labels[1:]):
        tseries_grid._legend.texts[idx].set_text(new_label)
    tseries_grid.set_ylabels(prop_labels[hue_attr]['symbol'])
    tseries_grid.set_xlabels(attr_labels[x_prop])
    tseries_grid.set_titles(
        attr_labels[row_attr] + r"$={row_name}$" +
        " | " + attr_labels[col_attr] + r"$={col_name}$"
        )
    tseries_grid.legend.set_title(legend_labels[0])
    tseries_grid.fig.suptitle(project_title, **fig_title_kws)
    tseries_grid.tight_layout(w_pad=1)
    sns.move_legend(tseries_grid, loc=loc, **move_legend_kws)
    return tseries_grid
