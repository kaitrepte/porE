import os

grids = [3,5,10,15,20]
structs = ['do','u6','u7','u8','ir','h1','ho','m5']

for z in range(len(grids)):
    for y in range(len(structs)):
        ffile = open('input_porosity','w')
        ffile.write("STRUCTURE\n")
        ffile.write(structs[y]+"                      Structure string. Several available options are given at the end of this file\n")
        ffile.write("n/a                     If Structure = 'ud': the full path to the xyz file you want to use\n")
        ffile.write("n/a                     If Structure = 'ud': the name of the structure (for output purposes only)\n")
        ffile.write("\n")
        ffile.write("EVALUATION METHOD\n")
        ffile.write("2                       1 : Overlapping sphere approach (OSA), 2 : Grid point approach (GPA)\n")
        ffile.write("\n")
        ffile.write("FOR GPA ONLY\n")
        ffile.write("2.16                    Probe radius for the evaluation of the accessible porosity\n")
        ffile.write("2                       Grid point initialization. 1 : grid points per unit cell vector, 2 : grid points per A (uniform grid)\n")
        ffile.write(str(grids[z])+"                       If initialization = 1: 3 integers for the number of grid points along each cell vector\n")
        ffile.close()

        # to be adjusted to your executable path
        os.system('/data/postdoc/porE/github/porE/src/porE_subgrid')
        os.system('mv output_porosity results/'+structs[y]+'_'+str(grids[z])+'.out')
