import pandas as pd
import seaborn as sns
from glob import glob
from PipeLine import *

properties_path = "../sumrule_data/properties-ens_avg-all_in_one.csv"
properties_ens_avg = pd.read_csv(properties_path,header=0)
fontsize=20
colors = sns.color_palette(n_colors=4, palette="Set2", desat=1)
leg_labels = [r'$\frac{D}{a}$ (Line color)', '15.0', '20.0', '25.0','30.0', r'$N$ (Line size)', 1000, 2000, r'$\frac{a_c}{a}$ (Marker type)', '1.0', '2.0', '4.0']
PipeLine.chainsize_plot(properties_ens_avg, "phi_c_eff", leg_labels, colors, fontsize=fontsize)
PipeLine.chainsize_plot(properties_ens_avg, "phi_c", leg_labels, colors, fontsize=fontsize)
PipeLine.chainsize_plot(properties_ens_avg, "phi_c_eff_normalized", leg_labels, colors, fontsize=fontsize)
PipeLine.chainsize_plot(properties_ens_avg, "phi_c_normalized", leg_labels, colors, fontsize=fontsize)