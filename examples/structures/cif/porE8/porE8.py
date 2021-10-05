from pore import porosity as p
from porE.hea.HEA import HEA
from pore import psd 
from porE.io.ase2pore import *
import glob

cifs = glob.glob('*.cif')
cifs.sort()
f_out = 'benchmark.out'
with open(f_out,'w',buffering=1) as o:
    o.write('file P(HEA) P(OSA) P(GPAvoid) P(GPAacc) PSD\n')
    for cif in cifs:
        # Perform HEA
        r1 = HEA(cif)
        # Transform cif into porE format
        ase2pore(cif)
        xyz = 'pypore.xyz'
        # Perform OSA
        r2 = p.osa(xyz)
        # perform PSD
        no_pores,tmp1,tmp2,tmp3,tmp4 = psd.get_psd(xyz,200,2000)
        pores           = tmp1[0:no_pores]
        distr           = tmp2[0:no_pores]
        pore_pos_cart   = tmp3[0:no_pores]
        pore_pos_frac   = tmp4[0:no_pores]
        # Perform GPA
        probe_R = 1.2
        grid_density = 5.0
        r3 = p.gpa_gridpera(xyz,probe_R,grid_density)
        # collect results
        r = [r1[0], r2[0], r3[0], r3[1]]
        #s = ' '.join([str(ri) for ri in r]+'\n')
        o.write('{:30.30s} {:15.3f} {:15.3f} {:15.3f} {:15.3f} {} {} {}\n'.format(cif,r[0],r[1],r[2],r[3],no_pores,pores,distr))
