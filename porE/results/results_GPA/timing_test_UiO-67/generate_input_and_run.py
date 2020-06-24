import os

grids = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0]

for z in range(len(grids)):
    ffile = open('input_porosity','w')
    ffile.write("STRUCTURE\n")
    ffile.write("u7                      Structure string. Several available options are given at the end of this file\n")
    ffile.write("n/a                     If Structure = 'ud': the full path to the xyz file you want to use\n")
    ffile.write("n/a                     If Structure = 'ud': the name of the structure (for output purposes only)\n")
    ffile.write("\n")
    ffile.write("EVALUATION METHOD\n")
    ffile.write("2                       1 : Overlapping sphere approach (OSA), 2 : Grid point approach (GPA)\n")
    ffile.write("\n")
    ffile.write("FOR GPA ONLY\n")
    ffile.write("1.20                    Probe radius for the evaluation of the accessible porosity\n")
    ffile.write("2                       Grid point initialization. 1 : grid points per unit cell vector, 2 : grid points per A (uniform grid)\n")
    ffile.write(str(grids[z])+"                       If initialization = 1: 3 integers for the number of grid points along each cell vector\n")
    ffile.close()

    os.system('/data/postdoc/porE/github/porE/src/porE_subgrid')

    os.system('mv output_porosity uio67_'+str(grids[z]))
