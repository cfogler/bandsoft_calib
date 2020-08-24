import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2 :
    print('Incorrect number of arguements. Instead use:')
    print('\tpython code.py [RunGlobalShift.txt]')
    exit(-1)


runNo = []
peak = []
std = []
with open(sys.argv[1],'r') as f:
    for line in f:
        parse = line.strip().split(" ")
        runNo.append(   int(parse[0])   )
        peak.append(    float(parse[3]) )
        std.append(     float(parse[4]) )

plt.errorbar(runNo,peak,yerr=std,linestyle='none',marker='o',color='blue')
plt.ylim([-6,3])
plt.xlabel('Run Number [a.u.]',fontsize=16)
plt.ylabel('Photon Peak Position [ns]',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.savefig('runshifts.pdf',bbox_inches='tight')
plt.show()
