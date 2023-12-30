# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 17:11:15 2023

@author: lucky

Aggiunti Archi Laterali
"""

from math import *
from decimal import *
import time
import cmath
import numpy as np

from DxfUtility import *

from PlotDxfUtility import *

class Point:
  def __init__(self, x, y, z):
      self.x = x
      self.y = y
      self.z = z
      
tempo_iniziale = time.time()
Radicedi2=(2.0**0.50)


B=3.75

ExternalArches = True
b=0.50
#ExternalArches = False
if not ExternalArches:
    b=0
    
spessoreVOLTA=0.30
gammaVOLTA=16.00
spessoreCAPPA=0.0
gammaCAPPA=25.00
spessoreRINFIANCO=0.30
gammaRINFIANCO=14.00

qperm=0.04*25.00+0.06*18.00+0.30+0.30
qvar=3.00

t=b/2
r=B/2
R=r*Radicedi2

numeroSuddivisioni=8 #deve essere pari
infittimentoSpigoli=True
#infittimentoSpigoli=False


def GeneraZ(Ri, Re, xy):    
    #xysqr=np.square(x)+np.square(y)
    xysqr=np.sum(np.square(xy), axis=1)
    zi=np.sqrt(Ri**2-xysqr)
    ze=np.sqrt(Re**2-xysqr)   
    return zi, ze


xo=-B/2
xn=+B/2

X=np.linspace(xo, xn, numeroSuddivisioni+1)
Y=np.linspace(xo, xn, numeroSuddivisioni+1)

XY = np.array(np.meshgrid(X, Y)).T.reshape(-1,2, order='F')# Fortran-like index ordering

NumeroNodiLATO=len(X)
NumeroNodiXY=len(XY)

if infittimentoSpigoli:
    u=4
    xy11=[(X[1]+(u-1)*X[0])/u, (Y[1]+(u-1)*Y[0])/u]
    xy12=[((u-1)*X[-1]+X[-2])/u, ((u-1)*Y[-1]+Y[-2])/u]
    xy13=[((u-1)*X[-1]+X[-2])/u, (Y[1]+(u-1)*Y[0])/u]
    xy14=[(X[1]+(u-1)*X[0])/u, ((u-1)*Y[-1]+Y[-2])/u]
    
    u=2
    xy21=[(X[1]+(u-1)*X[0])/u, (Y[1]+(u-1)*Y[0])/u]
    xy22=[((u-1)*X[-1]+X[-2])/u, ((u-1)*Y[-1]+Y[-2])/u]
    xy23=[((u-1)*X[-1]+X[-2])/u, (Y[1]+(u-1)*Y[0])/u]
    xy24=[(X[1]+(u-1)*X[0])/u, ((u-1)*Y[-1]+Y[-2])/u]
    
    XY=np.vstack((XY, xy11, xy21, xy22, xy12, xy13, xy23, xy24, xy14))
    listDiagAdd = np.arange(NumeroNodiXY, NumeroNodiXY+8)
    
    u=2
    xy01x=[(X[1]+(u-1)*X[0])/u, Y[0]]
    xy02x=[((u-1)*X[-1]+X[-2])/u, Y[0]]
    
    xy03x=[(X[1]+(u-1)*X[0])/u, Y[-1]]
    xy04x=[((u-1)*X[-1]+X[-2])/u, Y[-1]]
    XY=np.vstack((XY, xy01x, xy02x, xy03x, xy04x))

    xy01y=[X[0], (Y[1]+(u-1)*Y[0])/u]
    xy02y=[X[0], ((u-1)*Y[-1]+Y[-2])/u]
    
    xy03y=[X[-1], (Y[1]+(u-1)*Y[0])/u]
    xy04y=[X[-1], ((u-1)*Y[-1]+Y[-2])/u]
    XY=np.vstack((XY, xy01y, xy02y, xy03y, xy04y))    
    listBrdAdd = np.arange(listDiagAdd[-1]+1, listDiagAdd[-1]+1+8)
 
NumeroNodiXYInf=len(XY)
#GENERA QUOTE NODI INTERNI
Re = R+spessoreVOLTA
Zi, Ze = GeneraZ(R, Re, XY)    

#GENERA  NODI ARCHI ESTERNI

XBorderlistDown = np.arange(0,NumeroNodiLATO)
if infittimentoSpigoli:
    XBorderlistDown = np.insert(XBorderlistDown, 1, listBrdAdd[0])
    XBorderlistDown = np.insert(XBorderlistDown, -1, listBrdAdd[1])
   
XBorderlistUp = np.arange(NumeroNodiLATO*(NumeroNodiLATO-1),NumeroNodiXY) 
if infittimentoSpigoli:
    XBorderlistUp = np.insert(XBorderlistUp, 1, listBrdAdd[2])
    XBorderlistUp = np.insert(XBorderlistUp,-1, listBrdAdd[3])

YBorderlistSx=np.linspace(0,(NumeroNodiLATO-1)*NumeroNodiLATO,NumeroNodiLATO,dtype=np.int32)
if infittimentoSpigoli:
    YBorderlistSx = np.insert(YBorderlistSx, 1, listBrdAdd[4])
    YBorderlistSx = np.insert(YBorderlistSx, -1, listBrdAdd[5])

YBorderlistDx=np.linspace(NumeroNodiLATO-1, NumeroNodiXY-1, NumeroNodiLATO,dtype=np.int32)
if infittimentoSpigoli:
    YBorderlistDx = np.insert(YBorderlistDx, 1, listBrdAdd[6])
    YBorderlistDx = np.insert(YBorderlistDx, -1, listBrdAdd[7])

if ExternalArches:
    XYLeftArch=XY[YBorderlistSx]+[-t, 0.00]
    XYRightArch=XY[YBorderlistDx]+[+t, 0.00]
    ZiArches=Zi[XBorderlistDown]
    ZeArches=Ze[XBorderlistDown]
    Zi = np.hstack((Zi.T, ZiArches, ZiArches))
    Ze = np.hstack((Ze.T, ZeArches, ZeArches))


nNodiLato=len(XYLeftArch)

listLeftArch = np.arange(NumeroNodiXYInf,NumeroNodiXYInf+nNodiLato)
listRightArch = np.arange(NumeroNodiXYInf+nNodiLato,NumeroNodiXYInf+2*nNodiLato)

XY = np.vstack((XY, XYLeftArch, XYRightArch))

NumeroNodiInterni=len(XY)

#GENERA  NODI ESTERNI
#delta=abs(np.sqrt(R**2-XY[1,0]**2)+XY[1,1])*0.99
delta=spessoreVOLTA*Radicedi2

XlistDown = np.arange(1,NumeroNodiLATO-1)
XlistDown = XBorderlistDown[1:-1]
XYExtrnl1=XY[XlistDown]+[0.0, -delta]
listExtrnlXDown = np.arange(NumeroNodiInterni, NumeroNodiInterni+len(XlistDown))                   
    
XlistUp = np.arange(NumeroNodiLATO*(NumeroNodiLATO-1)+1, NumeroNodiXY-1)
XlistUp = XBorderlistUp[1:-1]
XYExtrnl2=XY[XlistUp]+[0.0, +delta]
listExtrnlXUp = np.arange(listExtrnlXDown[-1]+1, listExtrnlXDown[-1]+1+len(XlistUp))
    
YlistSx=np.linspace(NumeroNodiLATO,(NumeroNodiLATO-2)*NumeroNodiLATO, NumeroNodiLATO-2,dtype=np.int32)
YlistSx = YBorderlistSx[1:-1]
if ExternalArches:
    YlistSx = listLeftArch[1:-1]
listExtrnlYSx = np.arange(listExtrnlXUp[-1]+1,listExtrnlXUp[-1]+1+len(YlistSx))
XYExtrnl3=XY[YlistSx]+[-delta, 0.0]

YlistDx=np.linspace(2*NumeroNodiLATO-1,(NumeroNodiLATO-1)*NumeroNodiLATO-1,NumeroNodiLATO-2,dtype=np.int32)
YlistDx = YBorderlistDx[1:-1] 
if ExternalArches:
    YlistDx = listRightArch[1:-1]
XYExtrnl4=XY[YlistDx]+[+delta, 0.0]
listExtrnlYDx = np.arange(listExtrnlYSx[-1]+1,listExtrnlYSx[-1]+1+len(YlistDx))

XY = np.vstack((XY, XYExtrnl1, XYExtrnl2, XYExtrnl3, XYExtrnl4))

Spigolilist=[0, NumeroNodiLATO-1, (NumeroNodiLATO-1)*NumeroNodiLATO, NumeroNodiXY-1,
             listLeftArch[0], listRightArch[0], listLeftArch[-1], listRightArch[-1]]

listBrd=np.hstack((XlistDown, XlistUp, YlistSx, YlistDx, Spigolilist))

deltaSpigoli=((np.ones((2,4))*delta)*np.array([-1,1,-1,1])).T
deltaSpigoli[1,1]=-deltaSpigoli[1,1]
deltaSpigoli[2,1]=-deltaSpigoli[2,1]
if ExternalArches:
    deltaSpigoli = np.vstack((deltaSpigoli, deltaSpigoli))

XYExtrnl5=XY[Spigolilist]+deltaSpigoli
XY = np.vstack((XY, XYExtrnl5))

listExtrnlSpigoli=np.arange(listExtrnlYDx[-1]+1,listExtrnlYDx[-1]+1+len(Spigolilist))
listExtrnl = np.hstack((listExtrnlXDown,listExtrnlXUp, listExtrnlYSx, listExtrnlYDx, listExtrnlSpigoli))  


NumeroNodiTotali=len(XY)
NumeroNodiEsterni=NumeroNodiTotali-NumeroNodiInterni

Xo=XY[XlistDown][:,0]
Zo1=np.sqrt(R**2-Xo*Xo)
Zo2=Ze[XlistDown]
ZoX=Zo2+(Zo2-Zo1)*delta/r
ZoY=ZoX
zo=Ze[0]+(Ze[0]-R)*delta*Radicedi2/R
Zo=[zo, zo, zo, zo]

ZiExtrnlX=Zi[XlistDown]-np.ones(len(XlistDown))*(+5.0)
ZiExtrnlS=Zi[Spigolilist]-np.ones(len(Spigolilist))*(+5.0)
ZiExtrnlY=Zi[YlistSx]-np.ones(len(YlistSx))*(+5.0)
if ExternalArches:
    #ZoYmed=(Ze[YlistSx]+Zi[YlistSx])/2
    ZoY=Ze[YlistSx]#ZoYmed*1.01  
    ZiExtrnlY=Zi[YlistSx]#ZoYmed*0.99  
    Zo=np.hstack((Zo, Zo))
    
ZeExtrnl=np.hstack((ZoX, ZoX, ZoY, ZoY, Zo)) 


#ZeExtrnl=Ze[listBrd]


Zi = np.hstack((Zi.T, ZiExtrnlX, ZiExtrnlX, ZiExtrnlY, ZiExtrnlY, ZiExtrnlS))
Ze = np.hstack((Ze.T, ZeExtrnl))

#Zi[listExtrnl]=np.ones(len(listExtrnl))*(-5.0)

q=np.zeros(NumeroNodiLATO, dtype=np.int32)
if infittimentoSpigoli:
    q[0]=1
    q[-1]=1
  
    
#GENERAZIONE BRANCHES X  
Branches=[]
nodoJ0=0
if infittimentoSpigoli:
    nodoI=0
    nodoJ=NumeroNodiXY+8
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=nodoJ0+1
    Branches.append((nodoI,nodoJ))
    nodoJ0 += 1
    
nodoI=nodoJ0
for ix in range(q[0], NumeroNodiLATO-1-q[-1]):
    nodoJ=nodoI+1
    Branches.append((nodoI,nodoJ))
    nodoI+=1
    
if infittimentoSpigoli:
    nodoJ=nodoJ=NumeroNodiXY+9
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=NumeroNodiLATO-1
    Branches.append((nodoI,nodoJ))    

for iy in range(1, NumeroNodiLATO-1):
    nodoI=iy*NumeroNodiLATO
    for ix in range(0, NumeroNodiLATO-1):
        nodoJ=nodoI+1
        Branches.append((nodoI,nodoJ))
        nodoI+=1

nodoJ0=(NumeroNodiLATO-1)*NumeroNodiLATO
if infittimentoSpigoli:
    nodoI=nodoJ0
    nodoJ=NumeroNodiXY+10
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=nodoJ0+1
    Branches.append((nodoI,nodoJ))
    nodoJ0 += 1
    
nodoI=nodoJ0 
for ix in range(q[0], NumeroNodiLATO-1-q[-1]):
    nodoJ=nodoI+1
    Branches.append((nodoI,nodoJ))
    nodoI+=1
    
if infittimentoSpigoli:
    nodoJ=nodoJ=NumeroNodiXY+11
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=NumeroNodiXY-1
    Branches.append((nodoI,nodoJ))    

BranchesXmaxIndex=len(Branches)-1

#GENERAZIONE BRANCHES Y  
nodoJ0=0
if infittimentoSpigoli:
    nodoI=0
    nodoJ=NumeroNodiXY+12
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=nodoJ0+NumeroNodiLATO
    Branches.append((nodoI,nodoJ))
    nodoJ0 += NumeroNodiLATO
    
nodoI=nodoJ0
for iy in range(q[0], NumeroNodiLATO-1-q[-1]):
    nodoJ=nodoI+NumeroNodiLATO
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    
if infittimentoSpigoli:
    nodoJ=NumeroNodiXY+13
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=(NumeroNodiLATO-1)*NumeroNodiLATO
    Branches.append((nodoI,nodoJ))

for ix in range(1, NumeroNodiLATO-1):
    nodoI=ix
    for iy in range(0, NumeroNodiLATO-1):
        nodoJ=nodoI+NumeroNodiLATO
        Branches.append((nodoI,nodoJ))
        nodoI=nodoJ


nodoJ=NumeroNodiLATO-1
if infittimentoSpigoli:
    nodoI=nodoJ
    nodoJ=NumeroNodiXY+14
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=2*NumeroNodiLATO-1
    Branches.append((nodoI,nodoJ))

nodoI=nodoJ
for iy in range(q[0], NumeroNodiLATO-1-q[-1]):
    nodoJ=nodoI+NumeroNodiLATO
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ

if infittimentoSpigoli:
    nodoJ=NumeroNodiXY+15
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoJ=NumeroNodiXY-1
    Branches.append((nodoI,nodoJ))

numeroBrachesXYo=len(Branches)

BranchesYmaxIndex=numeroBrachesXYo-1

listBranchesX=np.arange(0, BranchesXmaxIndex+1)
listBranchesY=np.arange(BranchesXmaxIndex+1, BranchesYmaxIndex+1)

#GENERAZIONE BRANCHES INFITTIMENTO SPIGOLI
listBranchesSpigoliAdd=[]
if infittimentoSpigoli:
   nodoI=listBrdAdd[0]
   nodoJ=listDiagAdd[0]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=listBrdAdd[4]
   Branches.append((nodoI,nodoJ))    
   nodoI=1
   nodoJ=listDiagAdd[1]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=NumeroNodiLATO
   Branches.append((nodoI,nodoJ)) 
   
   
   nodoI=listBrdAdd[1]
   nodoJ=listDiagAdd[4]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=listBrdAdd[6]
   Branches.append((nodoI,nodoJ))    
   nodoI=NumeroNodiLATO-2
   nodoJ=listDiagAdd[5]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=2*NumeroNodiLATO-1
   Branches.append((nodoI,nodoJ))    

   nodoI=listBrdAdd[2]
   nodoJ=listDiagAdd[7]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=listBrdAdd[5]
   Branches.append((nodoI,nodoJ))    
   nodoI=NumeroNodiLATO*(NumeroNodiLATO-1)+1
   nodoJ=listDiagAdd[6]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=NumeroNodiLATO*(NumeroNodiLATO-2)
   Branches.append((nodoI,nodoJ)) 

   nodoI=listBrdAdd[3]
   nodoJ=listDiagAdd[3]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=listBrdAdd[7]
   Branches.append((nodoI,nodoJ))    
   nodoI=NumeroNodiXY-2
   nodoJ=listDiagAdd[2]
   Branches.append((nodoI,nodoJ))    
   nodoI=nodoJ
   nodoJ=NumeroNodiLATO*(NumeroNodiLATO-1)-1
   Branches.append((nodoI,nodoJ)) 
   listBranchesSpigoliAdd=np.arange(numeroBrachesXYo, numeroBrachesXYo+16)


numeroBrachesXY=len(Branches)

#numeroBrachesDiagonali=numeroBrachesInterni-numeroBrachesXY




#GENERAZIONE BRANCHES DIAGONALI
nodoI=0
start=0
if infittimentoSpigoli:
    nodoJ=NumeroNodiXY
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ += 1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ=NumeroNodiLATO+1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    start=1
    
for i in range(start, NumeroNodiLATO-1-start):
    nodoJ=nodoI+1+NumeroNodiLATO
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ

nodoI=NumeroNodiLATO-1
if infittimentoSpigoli:
    nodoI=nodoJ
    nodoJ=NumeroNodiXY+2
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ += 1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ=NumeroNodiXY-1
    Branches.append((nodoI,nodoJ))
    
    nodoI=NumeroNodiLATO-1
    nodoJ=NumeroNodiXY+4
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ     
    nodoJ += 1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ   
    nodoJ=2*NumeroNodiLATO-2
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ
    nodoI=2*NumeroNodiLATO-2

for i in range(start, NumeroNodiLATO-1-start):
    nodoJ=nodoI-1+NumeroNodiLATO
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    
if infittimentoSpigoli:    
    nodoJlast=nodoI-1+NumeroNodiLATO    
    nodoJ=NumeroNodiXY+6
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ += 1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ=nodoJlast
    Branches.append((nodoI,nodoJ))

numeroBrachesVolta=len(Branches)
numeroBrachesDiagonali=numeroBrachesVolta-numeroBrachesXY
listBranchesDiagonali=np.arange(numeroBrachesXY, numeroBrachesVolta)

#GENERAZIONE ARCHES BRANCHES   
if ExternalArches:
    n=0
    for nodoJ in YBorderlistSx: 
        nodoI=listLeftArch[n]
        Branches.append((nodoI,nodoJ))
        n += 1
    n=0
    for nodoJ in YBorderlistDx: 
        nodoI=listRightArch[n]
        Branches.append((nodoI,nodoJ))
        n += 1
    
    listBranchesXArches=np.arange(numeroBrachesVolta, len(Branches))
    
    n=0
    for nodoJ in listLeftArch[1:]: 
        nodoI=listLeftArch[n]
        Branches.append((nodoI,nodoJ))
        n += 1
    n=0
    for nodoJ in listRightArch[1:]: 
        nodoI=listRightArch[n]
        Branches.append((nodoI,nodoJ))
        n += 1

numeroBrachesInterni=len(Branches)

#GENERAZIONE BRANCHES ESTERNI
n=0
for nodoJ in listExtrnl: 
    nodoI=listBrd[n]
    Branches.append((nodoI,nodoJ))
    n += 1
    
numeroBrachesTotali=len(Branches)

numeroBranchesExtrnl=len(listExtrnl)

numeroBranchesDiag=4
if ExternalArches:
    numeroBranchesDiag=8

numeroBranchesExtrnlXY=int((numeroBranchesExtrnl-numeroBranchesDiag)/2)

listBranchesEsterniY = np.arange(numeroBrachesInterni, numeroBrachesInterni+numeroBranchesExtrnlXY)
listBranchesEsterniX = np.arange(numeroBrachesInterni+numeroBranchesExtrnlXY, numeroBrachesInterni+2*numeroBranchesExtrnlXY)

listBranchesDiagonaliE=np.arange(numeroBrachesInterni+2*numeroBranchesExtrnlXY, numeroBrachesTotali)

listBranchesDiagonali=np.hstack((listBranchesDiagonali,listBranchesDiagonaliE))
NumeroBranchesEsterni=numeroBrachesTotali-numeroBrachesInterni

listBranchesEsterni = np.arange(numeroBrachesInterni, numeroBrachesTotali)

listBranchesX=np.hstack((listBranchesX, listBranchesXArches, listBranchesEsterniX))
listBranchesY=np.hstack((listBranchesY,  listBranchesEsterniY))

#GENERAZIONE FACES
Faces=[]
p=np.zeros(NumeroNodiLATO-1, dtype=np.int32)
if infittimentoSpigoli:
    p[0]=1
    p[-1]=1
    
    for w in range(0, 2):
        nodoVertice=w*(NumeroNodiLATO-1)
        nodoI=nodoVertice
        nodoJ=NumeroNodiXY+8+w
        nodoK=NumeroNodiXY+w*4
        Faces.append((nodoI,nodoJ,nodoK,nodoI))
        
        nodoI=nodoJ
        nodoJ=1+(NumeroNodiLATO-3)*w
        nodoL=nodoK
        nodoK=nodoL+1
        Faces.append((nodoI,nodoJ,nodoK,nodoL))    
        
        nodoI=nodoJ
        nodoJ=nodoI+NumeroNodiLATO
        Faces.append((nodoI,nodoJ,nodoK,nodoK))   
      
        nodoI=w*(NumeroNodiLATO-1)
        nodoJ=NumeroNodiXY+w*4
        nodoK=NumeroNodiXY+12+w*2
        Faces.append((nodoI,nodoJ,nodoK,nodoI))
           
        nodoI=nodoJ
        nodoJ += 1
        nodoL=nodoK
        nodoK=w*(NumeroNodiLATO-1)+NumeroNodiLATO 
        Faces.append((nodoI,nodoJ,nodoK,nodoL))    
           
        nodoI=nodoJ
        nodoJ=nodoVertice+NumeroNodiLATO+1-2*w
        Faces.append((nodoI,nodoJ,nodoK,nodoK))   
    
    for w in range(0, 2):
        nodoVertice=w*(NumeroNodiLATO-1)+NumeroNodiLATO*(NumeroNodiLATO-1)
        nodoI=nodoVertice
        nodoJ=NumeroNodiXY+10+w
        nodoK=NumeroNodiXY+7-w*4
        Faces.append((nodoI,nodoJ,nodoK,nodoI))
        
        nodoI=nodoJ
        nodoJ=nodoVertice+1-2*w
        nodoL=nodoK
        nodoK=nodoL-1
        Faces.append((nodoI,nodoJ,nodoK,nodoL))    
        
        nodoI=nodoJ
        nodoJ=nodoI-NumeroNodiLATO
        Faces.append((nodoI,nodoJ,nodoK,nodoK))   
      
        nodoI=nodoVertice
        nodoJ=NumeroNodiXY+7-w*4
        nodoK=NumeroNodiXY+13+w*2
        Faces.append((nodoI,nodoJ,nodoK,nodoI))
           
        nodoI=nodoJ
        nodoJ -= 1
        nodoL=nodoK
        nodoK=nodoVertice-NumeroNodiLATO 
        Faces.append((nodoI,nodoJ,nodoK,nodoL))    
           
        nodoI=nodoJ
        nodoJ=nodoVertice-NumeroNodiLATO+1-2*w
        Faces.append((nodoI,nodoJ,nodoK,nodoK))     
    

for iy in range(0, NumeroNodiLATO-1):
    nodoI=iy*NumeroNodiLATO+p[iy]
    nodoL=(iy+1)*NumeroNodiLATO+p[iy]
    for ix in range(p[iy], NumeroNodiLATO-1-p[iy]):
        nodoJ=nodoI+1
        nodoK=nodoL+1
        Faces.append((nodoI,nodoJ,nodoK,nodoL))
        nodoI=nodoJ
        nodoL=nodoK

#GENERAZIONE ARCHES FACES   
if ExternalArches:
    n=0
    for nodoI in listLeftArch[:-1]: 
        nodoJ = YBorderlistSx[n]
        nodoK=YBorderlistSx[n+1]
        nodoL=listLeftArch[n+1]
        Faces.append((nodoI,nodoJ,nodoK,nodoL))
        n += 1
    n=0
    for nodoI in listRightArch[:-1]: 
        nodoJ = YBorderlistDx[n]
        nodoK=YBorderlistDx[n+1]
        nodoL=listRightArch[n+1]
        Faces.append((nodoI,nodoJ,nodoK,nodoL))
        n += 1
    

#******************************************************************************   
#******************************************************************************
#PLOTTAGGIO *******************************************************************
#******************************************************************************
#******************************************************************************

fileName="VelaRect"
dxf_fileName=fileName+"-INPUT.dxf"
print ("\nInizio Salvataggio file grafico ", fileName+"-INPUT.dxf...")

dxf_file = open(dxf_fileName,"w")    
WriteIntestazioneDXF(dxf_file)

EtichettaNodo="NODI-Intradosso"
#layerNodiIntra="NODI-Intradosso"
layerNodiEstra="NODI-Estradosso"
layerNodiText="NODI-Index"
for n in range(0, NumeroNodiTotali):         
    P=Point(XY[n, 0],XY[n, 1], Zi[n])
    PointToDXF(dxf_file, P, EtichettaNodo)
    testo=str(n)
    TextToDXF(dxf_file, P,testo,0.05,layerNodiText)
    P=Point(XY[n, 0],XY[n, 1], Ze[n])
    PointToDXF(dxf_file, P, layerNodiEstra)
    
    
EtichettaBranch = "BRANCHE_INTERNO"

layerallBranchesText="BRANCHES-Index"
for n in range(0, numeroBrachesTotali): 
    branch=Branches[n]
    iNode=branch[0]    
    jNode=branch[1]
    Pi=Point(XY[iNode, 0],XY[iNode, 1], Zi[iNode])
    Pj=Point(XY[jNode, 0],XY[jNode, 1], Zi[jNode])
    LineToDXF(dxf_file, Pi, Pj, EtichettaBranch)
    testo=str(n)
    Pm=MediaP([Pi,Pj])
    TextToDXF(dxf_file, Pm,testo,0.02,layerallBranchesText)   


#faces
NumeroFaces=len(Faces)
layerFaces="Faces"
layerFacesText="Faces-Index"
for n in range(0, NumeroFaces): 
    face=Faces[n]
    iNode=face[0]    
    jNode=face[1]
    kNode=face[2]    
    lNode=face[3]    
    
    Pi=Point(XY[iNode, 0],XY[iNode, 1], Zi[iNode])
    Pj=Point(XY[jNode, 0],XY[jNode, 1], Zi[jNode])
    Pk=Point(XY[kNode, 0],XY[kNode, 1], Zi[kNode])
    Pl=Point(XY[lNode, 0],XY[lNode, 1], Zi[lNode])    
    FaceToDXF(dxf_file, Pi, Pj, Pk, Pl, layerFaces)
    
    testo=str(n)
    Pm=MediaP([Pi, Pj, Pk, Pl])
    TextToDXF(dxf_file, Pm,testo,0.002,layerFacesText) 
    
    Pi=Point(XY[iNode, 0],XY[iNode, 1], Ze[iNode])
    Pj=Point(XY[jNode, 0],XY[jNode, 1], Ze[jNode])
    Pk=Point(XY[kNode, 0],XY[kNode, 1], Ze[kNode])
    Pl=Point(XY[lNode, 0],XY[lNode, 1], Ze[lNode]) 
    FaceToDXF(dxf_file, Pi, Pj, Pk, Pl, "EstradossoFaces")   


CloseWrittenDXF(dxf_file)

print ("...Salvataggio file grafico concluso")
print ("Dati salvati nel file " + fileName + "-INPUT.DXF")
print ("Salvati ", str(NumeroNodiTotali), " Nodi, ",
       str(numeroBrachesTotali),  " branches e ", str(NumeroFaces),  " Faces")


def Volumi(indxs):    
    x=[0]*3
    y=[0]*3
    z1=[0]*3
    z2=[0]*3
    dzSumV=0.        
    dzSumC=0
    dzSumR=0.
    for i in range(0, 3):    
        x[i]=XY[indxs[i], 0]
        y[i]=XY[indxs[i], 1]
        z1[i]=Zi[indxs[i]]
        z2[i]=Ze[indxs[i]]
        #z3[i]=ZEstradossoCAPPA[indxs[i]]
        dzSumV += z2[i]-z1[i]
        #dzSumC += z3[i]-z2[i]
        dzSumR += R+spessoreVOLTA+spessoreRINFIANCO-z2[i]
    S=0.50*abs(y[0]*(x[1]-x[2]) +y[1]*(x[2]-x[0]) +y[2]*(x[0]-x[1]))
    VlmV = S*dzSumV/3
    VlmR = S*dzSumR/3
    return S, VlmV, VlmR 

fzV=np.zeros(NumeroNodiInterni)
fzR=np.zeros(NumeroNodiInterni)
fzP=np.zeros(NumeroNodiInterni)
i=0
conta=0
for face in Faces:    
    indexes=list(sorted(set(face), key=face.index))
    if len(indexes)==3:
        S, VlmV, VlmR = Volumi(indexes)
        for j in range (0,3):
            fzV[face[j]]+=-VlmV/3 
            fzR[face[j]]+=-VlmR/3   
            fzP[face[j]]+=-S/3 
        conta += 1
        print("sono passato ", conta)
    else:
        indexes=face[:-1]
        S, VlmV, VlmR = Volumi(indexes)
        indexes = list(face)
        indexes.pop(1)
        S1, VlmV1, VlmR1 = Volumi(indexes)
        VlmV += VlmV1
        VlmR += VlmR1
        S += S1
        for j in range (0,4):
            fzV[face[j]]+=-VlmV/4 
            fzR[face[j]]+=-VlmR/4    
            fzP[face[j]]+=-S/4 
i+=1    


fzV *= gammaVOLTA
fzR *= gammaRINFIANCO
fzP *= (qperm+qvar)
PesoTotaleV = 0.0
PesoTotaleR = 0.0
PesoTotaleP = 0.0
i=0
PesoTotaleV=-np.sum(fzV)
PesoTotaleR=-np.sum(fzR)
PesoTotaleP=-np.sum(fzP)
PesoTotale=PesoTotaleV+PesoTotaleR+PesoTotaleP

fx=np.zeros(NumeroNodiTotali)
fy=np.zeros(NumeroNodiTotali)
fz=fzV+fzR+fzP
fz=np.hstack((fz, np.zeros(NumeroNodiEsterni)))


print("\nPeso Proprio Volta =", PesoTotaleV)
print("Peso rifianco      =", PesoTotaleR)
print("Peso Altri Carichi permanenti e Variabili      =", PesoTotaleP)
print("Peso Totale      =", PesoTotale)






tipoNodo=np.ones(NumeroNodiInterni, dtype=int)
tipoNodo=np.hstack((tipoNodo, np.zeros(NumeroNodiEsterni, dtype=int)))

NodesLayer = ["EXTERNAL_NODES", "INTERNAL_NODES"]


# *****************************************************************************
# SCRITTURA FILE PER APPLICAZIONE RosatiTNA.PY
#******************************************************************************

dat_file = open((fileName+".dat"),"w")

print ("\nInizio Esportazione Dati nel file ", fileName+".dat...")

dat_file.write(str(NumeroNodiTotali)+"\n")

n=0

while (n<NumeroNodiTotali):
    x=format(XY[n, 0], '+.6f') 
    y=format(XY[n,1], '+.6f')
    zmin=format(Zi[n], '+.6f')
    zmax=format(Ze[n], '+.6f')  
    fxn=format(fx[n], '+.3f')
    fyn=format(fy[n], '+.3f')
    fzn=format(fz[n], '+.3f')
    row = str(n+1)+"  "+ str(x)+"  "+str(y)+"  "+str(zmin)+"  "+str(zmax)
    row += "  "+str(fxn)+"  "+str(fyn)+"  "+str(fzn)
    row += "  "+str(tipoNodo[n])+"  "+str(NodesLayer[tipoNodo[n]])
    dat_file.write(row +"\n")
    n+=1

dat_file.write(str(numeroBrachesTotali)+"  "+ str(NumeroBranchesEsterni)+"\n")

dminBrachesX=0.10#0.5
dminBrachesY=0.20#0.5
dminBrachesDiagonali=.30
dminBrachesInfittimentoSpigoli=0.05

branchesdmin=np.ones(numeroBrachesTotali)*dminBrachesX
branchesdmin[listBranchesY]=dminBrachesY
branchesdmin[listBranchesSpigoliAdd]=np.ones(len(listBranchesSpigoliAdd))*dminBrachesInfittimentoSpigoli

branchesdmin[listBranchesDiagonali]=np.ones(len(listBranchesDiagonali))*dminBrachesDiagonali

for n in range(0, numeroBrachesTotali): 
    branch=Branches[n]
    iNode=branch[0]    
    jNode=branch[1]
    row = str(iNode+1)+" "+str(jNode+1)+" "+format(branchesdmin[n], '+.3f')
    #str(branches[b,0]+1)+" "+str(branches[b,1]+1)+" "+str(branchesdmin[b,0])+"\n
    dat_file.write(row+"\n")

nNLt=len(listExtrnlXDown)

dat_file.write(str(listBranchesEsterni[:nNLt])+ "\n")
dat_file.write(str(listBranchesEsterni[nNLt:2*nNLt])+ "\n")
dat_file.write(str(listBranchesEsterni[2*nNLt:3*nNLt])+ "\n")
dat_file.write(str(listBranchesEsterni[3*nNLt:4*nNLt])+ "\n")
dat_file.write(str(listBranchesEsterni[4*nNLt:])+ "\n")
dat_file.write(str(np.hstack((listLeftArch[1:-1], listRightArch[1:-1])))+ "\n")
dat_file.write(str(np.hstack((listExtrnlYSx, listExtrnlYDx)))+ "\n")
dat_file.write(str(PesoTotale))
dat_file.close()

print ("...Esportazione conclusa")
print ("Dati salvati nel file " + fileName + ".DAT")
print ("Esportati ", str(NumeroNodiTotali), " Nodi e ",
       str(numeroBrachesTotali),  " branches", )
tempo_finale = time.time()
print ("Impiegati ", str(tempo_finale - tempo_iniziale), " secondi per la generazione dei dati e la creazione dei file.")

