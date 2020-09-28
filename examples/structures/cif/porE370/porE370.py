from pore import porosity as p
from porE.hea.HEA import HEA
from porE.io.ase2pore import *
import glob

cifs = glob.glob('*.cif')
cifs.sort()
f_out = 'benchmark.out'
with open(f_out,'w',buffering=10) as o:
    o.write('file P(HEA) P(OSA) P(GPAvoid) P(GPAacc)\n')
    for cif in cifs:
        # Perform HEA
        r1 = HEA(cif)
        # Transform cif into porE format
        ase2pore(cif)
        xyz = 'pypore.xyz'
        # Perform OSA
        r2 = p.osa(xyz)
        # Perform GPA
        probe_R = 1.2
        grid_density = 5.0
        # collect results
        r = [r1[0], r2[0], r3[0], r3[1]]
        s = ' '.join([str(ri) for ri in r])
        o.write(s)
