"""
###########################################################
#
# funregulation_tf_heatmap.py
#
# FunRegulation: TF Heatmap (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB) / Universidade de Caxias do Sul (UCS)
# https://alexandrelenz@github.com/alexandrelenz/funregulation.git
# Last update: 08.05.2020 (Lenz, A. R.)
#
# Gene regulatory networks (GRN) of 
# Penicillium ucsensis 2HH and Penicillium oxalicum 114-2 
# inferred by a computational biology approach
#
# We propose the inference of global GRNs for Penicillium ucsensis
# 2HH and Penicillium oxalicum 114-2, based on TF-TG orthology relationships of
# related species combined with TFBSs prediction. First, global GRNs of related
# species (A. nidulans, N. crassa and S. cerevisiae) afford the mapping of
# orthologous interactions. Further, the TFBSs prediction provides accuracy to
# TF-TG relationships.
#
# In addition to the 37 Fungal TF families, two more families (HAP3 and HAP5) 
# deposited in interpro database, were included.
# This Python script plots a heatmap in PNG format that shows the distribution of 
# the TF families in the studied fungi.
#
#   INPUT FILE:
#
#    i) input_tf_file (CSV format)
#    
#   OUTPUT:
#    
#    i) heatmap ploted in PNG format
#
###########################################################
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns;

#Initialize Input File Path
input_tf_file = '/Users/arlenz/NetBeansProjects/funregulation/GRN/input/tfs/All_tf_domains.csv'

"""
    This Python script plots a Heatmap of TF Domains in PNG format
""" 
def plot_heatmap():
    input_tf_domains = pd.read_csv(input_tf_file)

    input_tf_domains = input_tf_domains.pivot('Species','Families', 'Quantity')

    ####  Free input parameters ####
    cols = [5, 39]

    # number of cells along width
    cells_in_row = 39

    # start by specifying figure width and margins
    figwidth= 25. #inch
    figheight=5.
    marg_top = 0.4
    marg_bottom = 2.
    marg_left = 3.
    marg_right = 0

    # set figure size
    fig = plt.figure(figsize=(figwidth, figheight))

    ####  Automatic calculation ####
    # determine cell size in inch
    cellsize = (figwidth-marg_left-marg_right)/float(cells_in_row)

    # now loop over cols:
    for cells_in_column in cols:
        # adjust margins (relative numbers) according to absolute values
        fig.subplots_adjust(bottom =marg_bottom/figheight ,top=1.-marg_top/figheight,
                            left=marg_left/figwidth, right=1.-marg_right/figwidth)
        # calculate figure height in inches
        figheight = cellsize*cells_in_column+marg_top+marg_bottom

    yticks = input_tf_domains.index
    xticks = input_tf_domains.columns

    # Blues_r = dark / Blues = light
    sns.heatmap(input_tf_domains, xticklabels=xticks, yticklabels=yticks, fmt="d", annot=True, cmap = 'Blues_r')

    # This sets the xticks "upright" with 0, as opposed to sideways with 90.
    plt.xticks(rotation=90)

    # This sets the yticks "upright" with 0, as opposed to sideways with 90.
    plt.yticks(rotation=0)

    plt.show()

"""
    Main function of this program
""" 
if __name__ == '__main__':
    plot_heatmap()