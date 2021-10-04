import os

structures = ['do','vo','dc','vc','u6','u7','u8','m5','ir','h1','ho']
mc_steps   = [250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000,20000]

for i in range(len(structures)):
    pores = []
    distr = []
    time = []
    for j in range(len(mc_steps)):
        try:
            ffile = open(structures[i]+'_'+str(mc_steps[j])+'.out','r')
            lines = ffile.readlines()
            ffile.close()

            tmp_pores = []
            tmp_distr = []
            for b in range(len(lines)):
                splitt = lines[b].split()
                if len(splitt) > 1:
                    if splitt[-1] == '(fractional)':
                        k = 0
                        l = 0
                        while k == 0:
                            l = l + 1
                            splitt2 = lines[b+l].split()
                            if len(splitt2) == 0:
                                k = 1
                            else:
                                tmp_pores.append(float(splitt2[0]))
                                tmp_distr.append(float(splitt2[1]))
                    if splitt[0] == 'Total':
                        time.append(float(splitt[3]))

            pores.append(tmp_pores)
            distr.append(tmp_distr)


        except FileNotFoundError:
            continue



    # if any pores were found
    print('Structure '+structures[i])
    if pores != []:
        for t in range(len(pores)):
            super_string = ''
            for ff in range(len(pores[t])):
                super_string = super_string+'  '+str('%10.4f' % pores[t][ff])+'('+str('%4.1f' % distr[t][ff])+')'
            print(str('%5.0i' % mc_steps[t])+' MC steps\t time = '+str('%10.4f' % time[t])+'\t'+super_string)

