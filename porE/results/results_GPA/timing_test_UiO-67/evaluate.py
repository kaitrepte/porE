import os
import time

grids = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0]

for u in range(len(grids)):
    ffile = open('uio67_'+str(grids[u]),'r')
    lines = ffile.readlines()
    ffile.close()

    for x in range(len(lines)):
        splitt = lines[x].split()
        if len(splitt) > 3:
            if splitt[0] == 'Total' and splitt[3] == 'grid':
                grid_points = int(splitt[5]) 

            if splitt[0] == 'Porosity' and splitt[1] == '(accessible):':
                porosity = float(splitt[2])

    time = float(lines[-1].split()[3])


    print('Grid density = '+str(grids[u])+'\tGrid points = '+str(grid_points)+'\t Phi_acc = '+str('%8.2f' % porosity)+' %\t Time = '+str('%10.3f' % time)+' s')

