import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

grids = [3,5,10,15,20]
structs = ['do','u6','u7','u8','ir','h1','ho','m5']
colors  = ['blue','green','orange','black','red','purple','brown','yellow']
names = ['DUT-8 open','UiO-66','UiO-67','UiO-68','IRMOF-10','HKUST-1,a','HKUST-1,b','MOF-5']

for a in range(len(structs)):
    grid = []
    grid_points = []
    time = []
    for b in range(len(grids)):
        ffile = open(structs[a]+'_'+str(grids[b])+'.out')
        lines = ffile.readlines()
        ffile.close()
        tmp_grid = []
        for c in range(len(lines)):
            splitt = lines[c].split()
            if len(splitt) > 2:
                if splitt[0] == 'Grid' and splitt[1] == 'which':
                    tmp_grid = [int(splitt[4]),int(splitt[6]),int(splitt[8])]
                    grid.append(tmp_grid)
                    grid_points.append(int(lines[c+1].split()[5]))
                if splitt[0] == 'Total' and splitt[1] == 'CPU':
                    time.append(float(splitt[3]))

    print('\nStructure '+names[a])
    for c in range(len(grid_points)):
        print(str('%4.0i' % grid[c][0])+' x '+str('%4.0i' % grid[c][1])+' x '+str('%4.0i' % grid[c][2])+' = '+str('%15.0i' % grid_points[c])+'\t Time = '+str('%10.4f' % time[c])+' s\t Grid points per second (grid_points/time) = '+str('%9.5f' % float(grid_points[c]/time[c])))
    for c in range(len(grids)-1):
        print('Ratio of grid points between grid '+str('%3.0i' % grids[c])+' and grid '+str('%3.0i' % grids[c+1])+' = '+str('%7.4f' % float(grid_points[c+1]/grid_points[c])))
        print('Ratio of time        between grid '+str('%3.0i' % grids[c])+' and grid '+str('%3.0i' % grids[c+1])+' = '+str('%7.4f' % float(time[c+1]/time[c])))

    # linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(grid_points,time)
    print('Regression slope = '+str('%12.6f' % slope)+'\tIntersection '+str('%12.6f' % intercept)+'\t R**2 = '+str('%12.6f' % (r_value**2)))
    for c in range(len(grid_points)):
        grid_points[c] = grid_points[c]/10**6   # for plotting, write grid points in Mio
    # PLOT results
    plt.plot(grid_points,time,'-o',color=colors[a],label=names[a])
    plt.legend(prop={'size':25},loc=4)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.xlim([-0.01,float(grid_points[-1])*1.1])
    plt.ylim([0,float(time[-1])*1.1])
    plt.xlabel('$10^{6}$ grid points $N$',fontsize=25)
    plt.ylabel('Time [s]',fontsize=25)
    
plt.show()
