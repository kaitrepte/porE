import math
import os

ffile = open('testt.xyz','r')
lines = ffile.readlines()

os.system('touch new')
file2 = open('new','w')

a = [-10.4891999990,0.0000000000,10.4891999990]
b = [0.0000000000,10.4891999990,10.4891999990]
c = [-10.4891999990,10.4891999990,0.0000000000]

for u in range(0,len(lines),1):
	splitt = lines[u].split()
	for aa in range(0,2,1):
		for bb in range(0,2,1):
			for cc in range(0,2,1):
				file2.write(splitt[0]+' '+str(float(splitt[1]) + (aa)*a[0] + (bb)*b[0] + (cc)*c[0])+' '+ \
							  str(float(splitt[2]) + (aa)*a[1] + (bb)*b[1] + (cc)*c[1])+' '+ \
							  str(float(splitt[3]) + (aa)*a[2] + (bb)*b[2] + (cc)*c[2])+'\n')
	
