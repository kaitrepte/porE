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
        r = [r1[0]]
        s = ' '.join([str(ri) for ri in r])
        o.write(s)
