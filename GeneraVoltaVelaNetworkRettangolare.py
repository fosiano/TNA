# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 17:11:15 2023

@author: lucky
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



# =============================================================================
# ExternalVoult = True
# b=0.50
# ExternalVoult = False
# 
# if not ExternalVoult:
#     b=0
#     
# =============================================================================
    
spessoreVOLTA=0.30
gammaVOLTA=16.00
spessoreCAPPA=0.0
gammaCAPPA=25.00
spessoreRINFIANCO=0.30
gammaRINFIANCO=14.00

qperm=0.04*25.00+0.06*18.00+0.30+0.30
qvar=3.00

#t=b/2
r=B/2
R=r*Radicedi2

numeroSuddivisioni=8 #deve essere pari
infittisiSpigoli=True
#infittisiSpigoli=False


def GeneraZ(Ri, Re, xy):    
    #xysqr=np.square(x)+np.square(y)
    xysqr=np.sum(np.square(xy), axis=1)
    zi=np.sqrt(Ri**2-xysqr)
    ze=np.sqrt(Re**2-xysqr)   
    return zi, ze

#delta=B/numeroSuddivisioni
xo=-B/2
xn=+B/2
X=np.linspace(xo, xn, numeroSuddivisioni+1)
Y=np.linspace(xo, xn, numeroSuddivisioni+1)
XY = np.array(np.meshgrid(X, Y)).T.reshape(-1,2, order='F')# Fortran-like index ordering
NumeroNodiLATO=len(X)
NumeroNodiXY=len(XY)

if infittisiSpigoli:
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

    u=2
    xy01=[(X[1]+(u-1)*X[0])/u, Y[0]]
    xy02=[((u-1)*X[-1]+X[-2])/u, Y[0]]
    
    xy03=[(X[1]+(u-1)*X[0])/u, Y[-1]]
    xy04=[((u-1)*X[-1]+X[-2])/u, Y[-1]]
    XY=np.vstack((XY, xy01, xy02, xy03, xy04))
    
    
NumeroNodiTOT=len(XY)

#GENERA QUOTE
Zi, Ze = GeneraZ(R, R+spessoreVOLTA, XY)

q=np.zeros(NumeroNodiLATO, dtype=np.int32)
if infittisiSpigoli:
    q[0]=1
    q[-1]=1
    
Branches=[]
nodoI=0
for iy in range(0, NumeroNodiLATO):
    nodoI=iy*NumeroNodiLATO+q[iy]
    for ix in range(q[iy], NumeroNodiLATO-1-q[iy]):
        nodoJ=nodoI+1
        Branches.append((nodoI,nodoJ))
        nodoI+=1


for ix in range(0, NumeroNodiLATO):
    nodoI=ix+q[ix]*NumeroNodiLATO
    for iy in range(q[ix], NumeroNodiLATO-1-q[ix]):
        nodoJ=nodoI+NumeroNodiLATO
        Branches.append((nodoI,nodoJ))
        nodoI=nodoJ

#DIAGONALI


nodoI=0
start=0
if infittisiSpigoli:
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
if infittisiSpigoli:
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
    
if infittisiSpigoli:    
    nodoJlast=nodoI-1+NumeroNodiLATO    
    nodoJ=NumeroNodiXY+6
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ += 1
    Branches.append((nodoI,nodoJ))
    nodoI=nodoJ 
    nodoJ=nodoJlast
    Branches.append((nodoI,nodoJ))

# nodoI=1
# nodoJ=NumeroNodiXY+1
# Branches.append((nodoI,nodoJ))   
# nodoI=nodoJ
# nodoJ=NumeroNodiLATO
# Branches.append((nodoI,nodoJ)) 
 

p=np.zeros(NumeroNodiLATO-1, dtype=np.int32)
if infittisiSpigoli:
    p[0]=1
    p[-1]=1
Faces=[]
for iy in range(0, NumeroNodiLATO-1):
    nodoI=iy*NumeroNodiLATO+p[iy]
    nodoL=(iy+1)*NumeroNodiLATO+p[iy]
    for ix in range(p[iy], NumeroNodiLATO-1-p[iy]):
        nodoJ=nodoI+1
        nodoK=nodoL+1
        Faces.append((nodoI,nodoJ,nodoK,nodoL))
        nodoI=nodoJ
        nodoL=nodoK

#******************************************************************************   
#******************************************************************************
#PLOTTAGGIO *******************************************************************
#******************************************************************************
#******************************************************************************

fileName="VelaRect"
dxf_fileName=fileName+"-INPUT.dxf"
dxf_file = open(dxf_fileName,"w")    
WriteIntestazioneDXF(dxf_file)

EtichettaNodo="NODI-Intradosso"
#layerNodiIntra="NODI-Intradosso"
layerNodiEstra="NODI-Estradosso"
layerNodiText="NODI-Index"
for n in range(0, NumeroNodiTOT):         
    P=Point(XY[n, 0],XY[n, 1], Zi[n])
    PointToDXF(dxf_file, P, EtichettaNodo)
    testo=str(n)
    TextToDXF(dxf_file, P,testo,0.05,layerNodiText)
    P=Point(XY[n, 0],XY[n, 1], Ze[n])
    PointToDXF(dxf_file, P, layerNodiEstra)
    
    
EtichettaBranch = "BRANCHE_INTERNO"
NumeroBranches=len(Branches)
layerallBranchesText="BRANCHES-Index"
for n in range(0, NumeroBranches): 
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
