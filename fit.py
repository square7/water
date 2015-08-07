import numpy as np
from sklearn.linear_model import *
import string
import pylab
from sklearn.ensemble import *
from sklearn.svm import *
import sys

f=open(sys.argv[1],'r')
ene=[]
e0=[]
e1=[]
OO=[]
OH=[]
HH=[]
comb=[]

while 1:
    s=string.split(f.readline())
    if len(s)==0:
        break
    eDiff=float(s[0])-float(s[1])
    e0.append(s[0])
    e1.append(s[1])
    ene.append(eDiff)
    o=[0]
    h=[0]
    oh=[0]
    f.readline()
    for cnt in range(99):
        s=[float(e) for e in string.split(f.readline())]
        o.append(s[0])
        h.append(s[1])
        oh.append(s[2])
    OO.append(o)
    OH.append(oh)
    HH.append(h)
    
OO=np.array(OO)[:,20:90]
OH=np.array(OH)[:,20:90]
HH=np.array(HH)[:,20:90]
comb=np.concatenate((OO,HH,OH),axis=1)
alp=[1e-15,.01,1]
c=['r','b','k']
ax1=pylab.subplot(221)
ax2=pylab.subplot(222)
ax3=pylab.subplot(223)
ax4=pylab.subplot(224)
for a in range(3):
    f=Ridge(alp[a])
    #f=LinearRegression()
    f.fit(comb, -np.array(ene))
    ax1.plot(1+np.arange(70)*.05,f.coef_[:70],c[a])
    ax2.plot(1+np.arange(70)*.05,f.coef_[70:140],c[a])
    ax3.plot(1+np.arange(70)*.05,f.coef_[140:],c[a])
    #ax4.hist(f.predict(comb)+np.array(ene),bins=100)
    ax4.scatter(f.predict(comb)+np.array(e0).astype(float), np.array(e1).astype(float))


ax1.set_xlabel('OO(a.u.)')
ax2.set_xlabel('HH(a.u.)')
ax3.set_xlabel('OH(a.u.)')
ax4.set_xlabel('residue')
pylab.show()
