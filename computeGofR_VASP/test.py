#!/usr/bin/env python

import subprocess
import numpy as np

wc = subprocess.Popen(["wc", "XDATCAR"], stdout=subprocess.PIPE)
awk = subprocess.Popen(["awk", "{print $1}"], stdin=wc.stdout, stdout=subprocess.PIPE)
#print awk
output = float(awk.communicate()[0])
NsnapshotsTotal = int(np.ceil((output-8)/65))-1
print NsnapshotsTotal

for i in range(NsnapshotsTotal):
    print "i, startline",i, 9+65*i

print int(9+65*np.ceil(NsnapshotsTotal/2))
print int(65*np.ceil(NsnapshotsTotal/2))
print 65*int(np.ceil(NsnapshotsTotal/2))
print NsnapshotsTotal
print 65*int(NsnapshotsTotal/20)
print 65*int(np.ceil(NsnapshotsTotal/2/10))
print int(np.ceil(65*NsnapshotsTotal/2/10))


