from porE.hea.HEA import HEA
from pore import porosity as p
from pore import psd 
from porE.io.ase2pore import *

def run_HEA(cif):
    '''
        Execute HEA, using a cif file

        Input
        -----
        cif ... cif file containing the structural information

        Output
        ------
        Phi_points   ... Porosity from inserted points, in %
        Phi_subcell  ... Porosity from subcell volumes, in %
        Phi_volume   ... Porosity from packing factor,  in %
        V_total      ... Total unit cell volume, in A^3
        V_subcell    ... Total subcell volume, in A^3
        V_insert     ... Volume of the inserted species, in A^3
        N_insert     ... Initial number of inserted species
        N_insert_red ... Reduced number of inserted species
        f_calc       ... calculated packing factor
        ropt         ... optimal radius
    '''

    print('')
    print('-------')
    print('Run HEA')
    print('-------')
    [Phi_points, Phi_subcell, Phi_volume, V_total, V_subcell, V_insert, N_insert, N_insert_red, f_calc, ropt] = HEA(cif)
    return Phi_points, Phi_subcell, Phi_volume, V_total, V_subcell, V_insert, N_insert, N_insert_red, f_calc, ropt


def run_OSA(xyz):
    '''
        Execute porosity evaluation, 
        using the overlapping sphere approach (OSA)

        Input
        -----
        xyz ... structural information in a xyz format compatible with porE

        Output
        ------
        Phi       ... Porosity, in %
        density   ... density of the structure, in kg/m^3
        poreV     ... pore volume density, in cm^3/g
        V_total   ... Total unit cell volume, in A^3
        V_vdwSum  ... Sum of the volumes of the vdW spheres of all atoms, in A^3
        V_overlap ... Overlap volume, in A^3
    '''

    print('')
    print('-------')
    print('Run OSA')
    print('-------')
    Phi, density, poreV, V_total, V_vdwSum, V_overlap = p.osa(xyz)
    return Phi, density, poreV, V_total, V_vdwSum, V_overlap

def run_GPA_fullgrid(xyz,probe_R,grid_a,grid_b,grid_c):
    '''
        Execute porosity evaluation, using the grid point apporach (GPA)
        Here: Provide an explicit size of the grid in directions a, b, c

        A probe radius for the accessible porosity needs to be provided

        Input
        -----
        xyz     ... structural information in a xyz format compatible with porE
        probe_R ... probe radius, in A
        grid_a  ... number of grid points in crystallographic direction a
        grid_b  ... number of grid points in crystallographic direction b
        grid_c  ... number of grid points in crystallographic direction c

        Output
        ------
        Phi_void   ... Void porosity, in %
        Phi_acc    ... Accessible porosity, in %
        density    ... density of the structure, in kg/m^3
        poreV_void ... pore volume density wrt void porosity, in cm^3/g
        poreV_acc  ... pore volume density wrt accessible porosity, in cm^3/g
    '''

    print('')
    print('-------------------------------')
    print('Run GPA: grid_a, grid_b, grid_c')
    print('-------------------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_fullgrid(xyz,probe_R,grid_a,grid_b,grid_c)
    return Phi_void, Phi_acc, density, poreV_void, poreV_acc


def run_GPA_gridPerA(xyz,probe_R,grid_density):
    '''
        Execute porosity evaluation, using the grid point apporach (GPA)
        Here: Provide a grid density along all cell vectors a, b, c
              porE will determine the amount of grid points automatically

        A probe radius for the accessible porosity needs to be provided

        Input
        -----
        xyz          ... structural information in a xyz format compatible with porE
        probe_R      ... probe radius, in A
        grid_density ... grid point density along all cell vectors, in points/A

        Output
        ------
        Phi_void   ... Void porosity, in %
        Phi_acc    ... Accessible porosity, in %
        density    ... density of the structure, in kg/m^3
        poreV_void ... pore volume density wrt void porosity, in cm^3/g
        poreV_acc  ... pore volume density wrt accessible porosity, in cm^3/g
    '''

    print('')
    print('---------------------')
    print('Run GPA: grid_density')
    print('---------------------')
    Phi_void, Phi_acc, density, poreV_void, poreV_acc = p.gpa_gridpera(xyz,probe_R,grid_density)
    return Phi_void, Phi_acc, density, poreV_void, poreV_acc

def run_PSD(xyz,N_points,N_steps):
    '''
        Execute an analyis of the pore size distribution (PSD)

        Input
        -----
        xyz      ... structural information in a xyz format compatible with porE
        N_points ... Number of starting points
        N_steps  ... Number of Monte-Carlo steps

        Output
        ------
        pores         ... all pore diameters, in A
        distr         ... the corresponding ditribution, in %
        pore_pos_cart ... cartesian  coordinates of pore centers, in A
        pore_pos_frac ... fractional coordinates of pore centers
    '''

    print('')
    print('-------')
    print('Run PSD')
    print('-------')
    no_pores,tmp1,tmp2,tmp3,tmp4 = psd.get_psd(xyz,N_points,N_steps)
    pores           = tmp1[0:no_pores]
    distr           = tmp2[0:no_pores]
    pore_pos_cart   = tmp3[0:no_pores]
    pore_pos_frac   = tmp4[0:no_pores]
    return pores, distr, pore_pos_cart, pore_pos_frac




def main(cif,do_HEA,do_OSA,do_GPA,do_PSD):
    '''
        Execute any evaluation you want
    '''
    name = cif.split('/')[-1].replace('.cif','')
    # start from a cif file used in the HEA;
    # convert cif to porE xyz (i.e., pypore.xyz) 
    # To be used for OSA, PSD and GPA
    ase2pore(cif)
    xyz = name+'.xyz'

    ##############
    # PARAMETERS #
    ##############
    # for PSD: N_points, N_steps
    # Recommendation: N_points >= 200, N_steps >= 1000
    N_points = 200
    N_steps  = 2000
    # For GPA: probe radius, grid density
    # Recommendation: grid_density >= 5
    probe_R = 1.20
    grid_density = 5

    ### Run HEA. Only cif file needed
    if do_HEA: 
        Phi_points, Phi_subcell, Phi_volume, V_total, V_subcell, V_insert, N_insert, N_insert_red, f_calc, ropt = run_HEA(cif)
    ### Run OSA. Only xyz file needed
    if do_OSA:
        Phi_OSA, density, poreV, V_total, V_vdwSum, V_overlap = run_OSA(xyz)
    ### Run PSD. Needs xyz, N_points and N_steps 
    if do_PSD:
        pores, distr, pore_pos_cart, pore_pos_frac = run_PSD(xyz,N_points,N_steps)
    ## Run GPA. Needs xyz, probe radius and grid
    if do_GPA:
        Phi_void, Phi_acc, density, poreV_void, poreV_acc = run_GPA_gridPerA(xyz,probe_R,grid_density)

    # Write a little summary
    output = open(name+'.dat','w') 
    print('\n### porE SUMMARY ###')
    output.write('### porE SUMMARY ###\n')
    if do_HEA:
        print('Phi_HEA(points)  = {:10.5f} %'.format(Phi_points))
        print('Phi_HEA(subcell) = {:10.5f} %'.format(Phi_subcell))
        print('Phi_HEA(volume)  = {:10.5f} %'.format(Phi_volume))
        output.write('Phi_HEA(points)   = {:10.5f} %'.format(Phi_points)+'\n')
        output.write('Phi_HEA(subcell)  = {:10.5f} %'.format(Phi_subcell)+'\n')
        output.write('Phi_HEA(volume)   = {:10.5f} %'.format(Phi_volume)+'\n')
    if do_OSA:
        print('Phi_OSA          = {:10.5f} %'.format(Phi_OSA))
        output.write('Phi_OSA           = {:10.5f} %'.format(Phi_OSA)+'\n')
    if do_GPA:
        print('Phi_GPA(void)    = {:10.5f} %'.format(Phi_void))
        print('Phi_GPA(acc)     = {:10.5f} % with r_probe = {:10.5f}'.format(Phi_acc,probe_R))
        output.write('Phi_GPA(void)     = {:10.5f} %'.format(Phi_void)+'\n')
        output.write('Phi_GPA(acc)      = {:10.5f} % with r_probe = {:10.5f}'.format(Phi_acc,probe_R)+'\n')
    if do_PSD:
        form = 'Pores [A]        '+'{:10.5f}'*len(pores)
        print(form.format(*pores))
        output.write(form.format(*pores)+'\n')
        form = 'Distribution [%] '+'{:10.2f}'*len(pores)
        print(form.format(*distr))
        output.write(form.format(*distr)+'\n')
    output.close()



if __name__ == '__main__':
    # start from a cif file
    cif = 'nu1000.cif'
    # Choose which calculation to perform
    do_HEA = True
    do_OSA = True
    do_GPA = True
    do_PSD = True
    main(cif=cif,do_HEA=do_HEA,do_OSA=do_OSA,do_GPA=do_GPA,do_PSD=do_PSD)
