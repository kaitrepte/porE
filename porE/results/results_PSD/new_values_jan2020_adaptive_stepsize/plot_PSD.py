import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline

structs = ['DUT8Ni_open_exp','UiO66',    'UiO67','UiO68','IRMOF10','MOF5','HKUST1_open_Cu_sites']
tags    = ['DUT-8(Ni) open', 'UiO-66',    'UiO-67','UiO-68','IRMOF-10','MOF-5',  'HKUST1']
colors  = ['red',           'darkblue','purple', 'cornflowerblue', 'gray',                'green',  'gold']

pores = []
distr = []

for a in range(len(structs)):
    ffile = open('output_'+structs[a],'r')
    lines = ffile.readlines()
    ffile.close()

    tmp_pores = []
    tmp_distr = []
    for b in range(len(lines)):
        splitt = lines[b].split()
        if len(splitt) > 1:
            if splitt[-1] == '(fractional)':
                i = 0
                j = 0
                while i == 0:
                    j = j + 1
                    splitt2 = lines[b+j].split()
                    if len(splitt2) == 0:
                        i = 1
                    else:
                       tmp_pores.append(float(splitt2[0])) 
                       tmp_distr.append(float(splitt2[1])) 

    # sort both temporary list together
    tmp_pores, tmp_distr = zip(*sorted(zip(tmp_pores,tmp_distr)))

    pores.append(tmp_pores)
    distr.append(tmp_distr)

#### PLOT
for a in range(len(structs)):
    plot_pores = []
    plot_distr = []
    for b in range(len(pores[a])):
        plot_pores.append(pores[a][b]-0.050)
        plot_pores.append(pores[a][b]-0.025)
        plot_pores.append(pores[a][b])
        plot_pores.append(pores[a][b]+0.025)
        plot_pores.append(pores[a][b]+0.050)
        plot_distr.append(0.0)
        plot_distr.append(distr[a][b])
        plot_distr.append(distr[a][b])
        plot_distr.append(distr[a][b])
        plot_distr.append(0.0)


#    plt.plot(plot_pores,plot_distr,'-',color=colors[a],linewidth=2.0,label=tags[a])
    plt.fill_between(plot_pores,plot_distr,color=colors[a],linewidth=0.01,label=tags[a])
###    plt.plot(xnew,ynew,'-o',color='green',label=structs[a])
    plt.xlabel(r'Pore diameter [$\mathrm{\AA}$]',fontsize=25)
    plt.ylabel(r'Distribution [%]',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.xlim([4.0,19])
    plt.ylim([-0.05,100.5])
    plt.legend(prop={'size':25},loc=2)
plt.show()
