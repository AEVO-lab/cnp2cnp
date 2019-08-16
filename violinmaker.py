import matplotlib.pyplot as plt
import matplotlib.axes as axes
import csv
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys


xlabel_err_rate = False
xlabel_forced = False
xlabel = "r = 0.1"
xforced_prefix = "R"
indir = "results3_s001l1x1"


def set_axis_style(ax, labels, positions):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylim([0,1])
    #ax.set_xlim(0.25, len(labels) + 0.5)
    #ax.set_xlabel('Sample name')




def parse_csv(infile, outfile, error_rate = "0"):
    heur = []
    euc = []
    flat = []
    zzs = []
    csvfile = open(infile,'r')
    allrows = csv.reader(csvfile, delimiter=',')
    plots = []
    for row in allrows:
        if row[10] == error_rate:
            plots.append(row)

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

    if xlabel_err_rate:
        ax.set_xlabel("error rate " + error_rate)
    else:
        ax.set_xlabel(lbl)

    if xlabel_forced:
        ax.set_xlabel(xlabel)

    plt.tight_layout()
    plt.savefig(outfile)


for error_rate in ["0", "0.1", "0.25", "0.5", "1"]:
    for fname in os.listdir(indir):
        if fname.endswith(".csv"):
            fullfname = indir + "/" + fname
            errstr = error_rate.replace(".", "")
            newfullfname = fullfname.replace("0.25.csv", "25.csv").replace("0.5.csv", "5.csv").replace("0.75.csv", "75.csv").replace(".csv", "_err" + errstr + ".pdf")
            if xlabel_err_rate:
                newfullfname = newfullfname.replace(".pdf", "_E.pdf")
            if xlabel_forced:
                newfullfname = newfullfname.replace(".pdf", "_" + xforced_prefix + ".pdf")
            parse_csv(fullfname, newfullfname, error_rate)
        



