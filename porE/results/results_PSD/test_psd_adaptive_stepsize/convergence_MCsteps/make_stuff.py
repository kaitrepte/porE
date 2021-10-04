import os

def write_input(struct,steps):
    ffile = open('input_PSD','w')
    ffile.write("STRUCTURE\n")
    ffile.write(struct+"\n")
    ffile.write("bla\n")
    ffile.write("bla\n")
    ffile.write("\n")
    ffile.write("CALCULATION PARAMETERS\n")
    ffile.write("200                             Number of starting points.      Recommended > 100\n")
    ffile.write(str(steps)+"                      Number of Monte-Carlo steps.    Recommended > 10000    \n")

#structures = ['do','vo','dc','vc','u6','u7','u8','m5','ir','h1','ho']
structures = ['do','ho']
mc_steps   = [250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,7000,8000,9000,10000,20000]

for i in range(len(structures)):
    for j in range(len(mc_steps)):
        write_input(structures[i],mc_steps[j])
        os.system('/data/postdoc/porE/github/porE/src/a.out')
        os.system('mv output_PSD results/'+structures[i]+'_'+str(mc_steps[j])+'.out')
