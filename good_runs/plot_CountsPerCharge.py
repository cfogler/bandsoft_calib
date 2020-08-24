import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2 :
    print('Incorrect number of arguements. Instead use:')
    print('\tpython code.py [CountsPerCharge.txt]')
    exit(-1)


runNo = []
info = []
bkg = []
with open(sys.argv[1],'r') as f:
    for line in f:
        parse = line.strip().split(" ")
        if float(parse[2]) < 0 : continue
        runNo.append(   int(parse[0])   )
        info.append(    float(parse[3]) )
        bkg.append(    float(parse[4]) )

plt.figure(1)
plt.errorbar(runNo,info,linestyle='none',marker='o',color='green')
plt.xlabel('Run Number [a.u.]',fontsize=16)
plt.ylabel('Signal Counts per Charge [mC]',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.savefig('signalcharge.pdf',bbox_inches='tight')

plt.figure(2)
plt.errorbar(runNo,bkg,linestyle='none',marker='o',color='red')
plt.xlabel('Run Number [a.u.]',fontsize=16)
plt.ylabel('Background Counts per Charge [mC]',fontsize=16)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.savefig('bkgcharge.pdf',bbox_inches='tight')

plt.show()
