# This script do the boxplot
import matplotlib.pyplot as plt
import pandas
from pandas import DataFrame
import numpy as np
# Set the default colors
import brewer2mpl
import matplotlib as mpl
#http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=5
# colors = ['empirical','model','sample','bound','test']
# colors = ['#fdae61','#ffffbf','#abd9e9','#2c7bb6','#d7191c']
# colors = ['#fdae61','#ffffbf','#abdda4','#2b83ba','#d7191c']
colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
lw = 1.5 # the line width of the dashed lines
s = 20 # the size of dot
widthBox = 0.1 # the width of the box

lenTotal = 25 # the length of total position
lenOne = 0.15 # the length of one postion

address = './plot/case33box.xlsx'

# err_gbCPS = pandas.read_excel(address, sheet_name='gbErrCPS', header=None)
# err_gbFirst = pandas.read_excel(address, sheet_name='gbErrFirst', header=None)
# err_gbSecond = pandas.read_excel(address, sheet_name='gbErrSecond', header=None)
# bound_gb = pandas.read_excel(address, sheet_name='gbBound', header=None)

err_gCPS = pandas.read_excel(address, sheet_name='gErrCPS', header=None)
err_bCPS = pandas.read_excel(address, sheet_name='bErrCPS', header=None)
err_gFirst = pandas.read_excel(address, sheet_name='gErrFirst', header=None)
err_bFirst = pandas.read_excel(address, sheet_name='bErrFirst', header=None)
err_gSecond = pandas.read_excel(address, sheet_name='gErrSecond', header=None)
err_bSecond = pandas.read_excel(address, sheet_name='bErrSecond', header=None)

bound_g = pandas.read_excel(address, sheet_name='gBound', header=None)
bound_b = pandas.read_excel(address, sheet_name='bBound', header=None)

# val_gEvalCPS = pandas.read_excel(address, sheet_name='gEvalCPS', header=None)
# val_bEvalCPS = pandas.read_excel(address, sheet_name='bEvalCPS', header=None)
# val_gEvalFirst = pandas.read_excel(address, sheet_name='gEvalFirst', header=None)
# val_bEvalFirst = pandas.read_excel(address, sheet_name='bEvalFirst', header=None)
# val_gEvalSecond = pandas.read_excel(address, sheet_name='gEvalSecond', header=None)
# val_bEValSecond = pandas.read_excel(address, sheet_name='bEvalSecond', header=None)

# val_gReal = pandas.read_excel(address, sheet_name='gReal', header=None)
# val_bReal = pandas.read_excel(address, sheet_name='bReal', header=None)

numBus = err_gCPS.shape[0]

# pos_err_gbCPS = np.linspace(1, lenTotal, numBus)
# pos_err_gbFirst = np.linspace(1+lenOne, lenTotal+lenOne, numBus)
# pos_err_gbSecond = np.linspace(1+2*lenOne, lenTotal+2*lenOne, numBus)
# pos_bound_gb = np.linspace(1+lenOne*3,lenTotal+lenOne*3,numBus)

pos_err_gFirst = np.linspace(1, lenTotal, numBus)
pos_err_gSecond = np.linspace(1+lenOne, lenTotal+lenOne, numBus)
pos_err_gCPS = np.linspace(1+lenOne*2, lenTotal+lenOne*2, numBus)
pos_bound_g = np.linspace(1+lenOne*3, lenTotal+lenOne*3, numBus)
pos_labels = np.linspace(1+lenOne*1.5, lenTotal+lenOne*1.5, numBus)

fig = plt.figure(figsize = (lenTotal,5))
ax1 = fig.add_subplot(111)

# c = colors[0]
# ax1.boxplot(err_gFirst,positions=pos_err_gFirst,widths=widthBox,
#     boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
#     flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

# c = colors[1]
# ax1.boxplot(err_gSecond,positions=pos_err_gSecond,widths=widthBox,
#     boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
#     flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

# c = colors[2]
# ax1.boxplot(err_gCPS,positions=pos_err_gCPS,widths=widthBox,
#     boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
#     flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

# for i in range(numBus):
#     plt.vlines(x=pos_bound_g[i], ymin=0, ymax=bound_g[0][i], lw=lw, colors=colors[3], linestyles = "dashed")
#     plt.scatter(x=pos_bound_g[i], y=bound_g[0][i], color=colors[3], s=s)

c = colors[0]
ax1.boxplot(err_bFirst,positions=pos_err_gFirst,widths=widthBox,
    boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
    flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

c = colors[1]
ax1.boxplot(err_bSecond,positions=pos_err_gSecond,widths=widthBox,
    boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
    flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

c = colors[2]
ax1.boxplot(err_bCPS,positions=pos_err_gCPS,widths=widthBox,
    boxprops=dict(color=c),capprops=dict(color=c),whiskerprops=dict(color=c),
    flierprops=dict(color=c, markeredgecolor=c, markersize=s/6),medianprops=dict(color=c))

for i in range(numBus):
    plt.vlines(x=pos_bound_g[i], ymin=0, ymax=bound_b[0][i], lw=lw, colors=colors[3], linestyles = "dashed")
    plt.scatter(x=pos_bound_g[i], y=bound_b[0][i], color=colors[3], s=s)


# draw temporary red and blue lines and use them to create a legend
h0, = plt.plot([0.25,0.1],color=colors[0])
h1, = plt.plot([0.25,0.1],color=colors[1])
h2, = plt.plot([0.25,0.1],color=colors[2])
h3, = plt.plot([0.25,0.1],color=colors[3],linestyle = "dashed")
plt.legend((h0, h1, h2, h3),('First-order', 'Second-order', 'CPS', 'Bound'),loc='lower right')
h0.set_visible(False)
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)

# plt.ylabel('Estimation Error of Line Conductance(p.u.)')
plt.ylabel('Estimation Error of Line Susceptance(p.u.)')
# plt.xlabel('Line Numbers')

name_labels = [str(i+1) for i in range(numBus)]
ax1.set_xticklabels(name_labels)
ax1.set_xticks(pos_labels)

plt.yscale('log')
plt.show()