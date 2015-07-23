import numpy as np
import sys
from numpy.linalg import *
import string

OO=64
tot=OO*3
st=int(sys.argv[1])
pos=open('P22.100.pos','r')
force=open('P22.100.for','r')
cel=open('P22.100.cel','r')

for cnt in range(st*(tot+1)+1): #skip prev steps and head(+1)
    pos.readline()
    force.readline()
for cnt in range(st*4+1):
    cel.readline()

op=open('10000.pos','w')
oc=open('10000.cel','w')
of=open('10000.for','w')
for cnt in range(OO):
    sp=[float(ele) for ele in string.split(pos.readline())]
    sf=[float(ele) for ele in string.split(force.readline())]
    print >>op, 'O', sp[0], sp[1], sp[2]
    print >>of, 'O', sf[0], sf[1], sf[2]

for cnt in range(OO*2):
    sp=[float(ele) for ele in string.split(pos.readline())]
    sf=[float(ele) for ele in string.split(force.readline())]
    print >>op, 'H', sp[0], sp[1], sp[2]
    print >>of, 'H', sf[0], sf[1], sf[2]
    
c=[]
for cnt in range(3):
    sc=[float(ele) for ele in string.split(cel.readline())]
    c.append(sc)

c=np.transpose(np.array(c))
for ele in c:
    print>>oc, ele[0],ele[1],ele[2]

inc=inv(c)
for ele in inc:
    print>>oc, ele[0],ele[1],ele[2]

