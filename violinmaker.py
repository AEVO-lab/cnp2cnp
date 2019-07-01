import matplotlib.pyplot as plt
import matplotlib.axes as axes
import csv
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import os



indir = "results_rcg"

def set_axis_style(ax, labels, positions):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylim([0,1])
    #ax.set_xlim(0.25, len(labels) + 0.5)
    #ax.set_xlabel('Sample name')




def parse_csv(infile, outfile):
    heur = []
    euc = []
    flat = []
    zzs = []
    csvfile = open(infile,'r')
    plots = csv.reader(csvfile, delimiter=',')
    nr = 0
    for row in plots:
        if nr != 0:
            flat.append(float(row[7].replace(",",".")))
            heur.append(float(row[6].replace(",",".")))
            euc.append(float(row[8].replace(",",".")))
            zzs.append(float(row[9].replace(",",".")))
        nr+=1

    font = {'family' : 'arial', 'size'   : 16}

    plt.rc('font', **font)
    

    labels = ["heuristic", "flat", "ZZS", "Euclidean"]
    positions = [1,1.25,1.5, 1.75]
    fig, ax = plt.subplots()
    ax.violinplot([heur, flat, zzs, euc], widths=0.25, positions = positions, showmeans=True, showextrema=True, showmedians=True, bw_method=0.5)
    set_axis_style(ax, labels, positions)

    lbl = os.path.basename(infile).replace("_", " ").replace(".csv", "") 
    lbl = lbl.replace("leaves ", "l = ")
    lbl = lbl.replace("genes ", "n = ")

    ax.set_xlabel(lbl)

    plt.savefig(outfile)



for fname in os.listdir(indir):
    if fname.endswith(".csv"):
        fullfname = indir + "/" + fname
		#ff = open(fullfname)
        newfullfname = fullfname.replace(".csv", ".pdf")
        parse_csv(fullfname, newfullfname)
        



