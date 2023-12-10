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
      
def prova(a):
    return "ok", a

Radicedi2=(2.0**0.50)
B=3.75
b=0.5
t=b/2
r=B/2
R=r*Radicedi2
#ao=np.arctan(t*Radicedi2/(2*R+t*Radicedi2))

NsA=6
NsB=8

a45=np.pi/4.
daA=a45/NsA
daB=a45/NsB

dgamma=2*a45/(NsA-1)

def GeneraNodi(R,alfa, betaV):   
    z=R*np.sin(betaV)
    ro=z/np.tan(betaV)
    x=ro*np.cos(alfa)
    y=ro*np.sin(alfa)    
    return [x, y, z]

def betaArco(r, R, alfa):
    beta=np.arcsin(r/R*np.sqrt(1-(np.tan(np.abs(alfa)))**2))
    return beta

def GeneraNodiEstradosso(Re,x,y):
    z=np.sqrt(Re**2-(np.square(x)+np.square(y)))
    return z

def MediaP(listaP):
    n=len(listaP)
    xm=0
    ym=0
    zm=0
    for i in range(0, n):
        xm+=listaP[i].x
        ym+=listaP[i].y
        zm+=listaP[i].z
    xm=xm/n
    ym=ym/n
    zm=zm/n
    P=Point(xm,ym,zm)
    return P

#suddivisione degli alfa tale che gli archi sull'arco di peimetro siano di equilughezza
gammaV=np.linspace(-2*a45, 2*a45, 2*NsA-1)
alfaV=np.arctan(np.sin(gammaV))
# alfaV[1]=alfaV[0]+ao
# alfaV[-2]=alfaV[-1]-ao
#alfaV=np.linspace(0.0, a45, NsA)
betaV0=np.linspace(a45, 2*a45, NsB)#da 45° a 90° (90 compreso)


betaArc=betaArco(r, R, alfaV[-(NsA-1):])
# =============================================================================
# betaArc=betaArco(r, R, alfaV[:NsA-1])
# listbetaArc=[betaArc]
# for t in range(1,betaArc.size):
#     listbetaArc.append(betaArc[t:])
# listbetaArc.append([])
# betaArcFlipped=np.flipud(np.copy(betaArc))
# =============================================================================

betaV=np.copy(betaV0)
listaBetaV=[betaV]# max da 0° a 90° (90 compreso)
for i in range(0, NsA-1):
    betaV=np.insert(betaV,0, betaArc[i]) 
    listaBetaV.append(betaV)
    listaBetaV.insert(0,betaV)

i=0
for betaV in listaBetaV: 
     betaV=np.flipud(betaV)
     listaBetaV[i]=betaV
     i+=1

alfa=0.0
beta=2*a45
NodoColmo=GeneraNodi(R,alfa,beta)
NodiOrdinatiPerMeridiani=[]
i=0
for alfa in alfaV:
    betaV=listaBetaV[i]       
    NodiMeridiano=GeneraNodi(R,alfa,betaV[1:])   
    NodiOrdinatiPerMeridiani.append(np.array(NodiMeridiano))      
    i+=1


def NodiArcoPerimetro(NodiOrdinatiPerMeridiani):
    #Estendi nodi ARCO DI PERIMETRO
    NodiMeridiano = NodiOrdinatiPerMeridiani[0]
    xyzt=np.array(NodiMeridiano)[::, -1]
    k=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]:       
        xyz=np.array(NodiMeridiano)[::, -1]      
        xyzt=np.vstack((xyzt,xyz))              
        k+=1
    return xyzt#(np.array([x, y, z])).T


def EstendiArcoPerimetro(NodiArco, gap,indx, dz):
    #Estendi nodi ARCO DI PERIMETRO 
    NodiArcoCopy=np.copy(NodiArco)
    NodiArcoCopy.T[indx] += gap
    NodiArcoCopy.T[2] += dz
    #print(xyz)
    return NodiArcoCopy

NodiOrdinatiPerMeridiani180=[]
NodiOrdinatiPerMeridiani90=[]
NodiOrdinatiPerMeridiani270=[]
   
for NodiMeridiano in NodiOrdinatiPerMeridiani:
    NodiMeridiano90=np.copy(NodiMeridiano)
    NodiMeridiano90x=-np.copy(NodiMeridiano[1])
    NodiMeridiano90y=np.copy(NodiMeridiano[0])
    NodiMeridiano90[0]=NodiMeridiano90x
    NodiMeridiano90[1]=NodiMeridiano90y
    NodiOrdinatiPerMeridiani90.append(NodiMeridiano90)    
    
    NodiMeridiano180=np.copy(NodiMeridiano)
    NodiMeridiano180x=-np.copy(NodiMeridiano[0])
    NodiMeridiano180y=-np.copy(NodiMeridiano[1])
    NodiMeridiano180[0]=NodiMeridiano180x
    NodiMeridiano180[1]=NodiMeridiano180y        
    NodiOrdinatiPerMeridiani180.append(NodiMeridiano180)
        
    NodiMeridiano270=np.copy(NodiMeridiano)
    NodiMeridiano270x=np.copy(NodiMeridiano[1])
    NodiMeridiano270y=-np.copy(NodiMeridiano[0])
    NodiMeridiano270[0]=NodiMeridiano270x
    NodiMeridiano270[1]=NodiMeridiano270y    
    NodiOrdinatiPerMeridiani270.append(NodiMeridiano270)    


NodiOrdinatiPerMeridianiToT360=NodiOrdinatiPerMeridiani+NodiOrdinatiPerMeridiani90[1:]
NodiOrdinatiPerMeridianiToT360+=NodiOrdinatiPerMeridiani180[1:]+NodiOrdinatiPerMeridiani270[1:-1]

CoordinateNodidArcoPerimetro0=NodiArcoPerimetro(NodiOrdinatiPerMeridiani)
CoordinateNodidArcoPerimetro90=NodiArcoPerimetro(NodiOrdinatiPerMeridiani90)
CoordinateNodidArcoPerimetro180=NodiArcoPerimetro(NodiOrdinatiPerMeridiani180)
CoordinateNodidArcoPerimetro270=NodiArcoPerimetro(NodiOrdinatiPerMeridiani270)

TrueArcNodesY=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro0,t,0,  0.)
TrueArcNodesY180=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro180,-t,0,  0.)

FalseArcNodesY=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro0[1:-1],t+0.60, 0, +0.125)
FalseArcNodesY180=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro180[1:-1],-(t+0.60),0, +0.125)

FalseArcNodesX90=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro90,+0.60, 1, -3.)
FalseArcNodesX270=EstendiArcoPerimetro(CoordinateNodidArcoPerimetro270,-0.60, 1, -3.)

spessore=0.30
ZNodiEstradossoOrdinatiPerMeridiani=[]
i=0
for NodiMeridiano in NodiOrdinatiPerMeridianiToT360:       
    Z=GeneraNodiEstradosso(R+spessore,NodiMeridiano[0],NodiMeridiano[1])
    ZNodiEstradossoOrdinatiPerMeridiani.append(np.array(Z))   
    i+=1
ZnodoColmoEstradosso=NodoColmo[2]+spessore


    
def index(k, n, betaProgressiveLen):
    indx=betaProgressiveLen+n
    return indx

def betaProgressiveLen(NodiOrdinatiPerMeridiani):
    dimension=[0]
    dimProgressive=0
    dim=[]
    for NodiMeridiano in NodiOrdinatiPerMeridiani:
        d=len(NodiMeridiano[0])
        dim.append(d)
        dimProgressive+=d
        dimension.append(dimProgressive)
    return dimension, dim


def RiordinaNodiMeridiani(NodiOrdinatiPerMeridiani, betaProgressiveLen):
    #Riordina NODI E BRANCHES NodiOrdinatiPerMeridiani
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2]         
    xyzToT=np.array([[xo], [yo], [zo]])
    k=0
    MeridianBranches=[]
    ZEstradosso=[ZnodoColmoEstradosso]
    ZNodiEstradossoOrdinatiPerMeridiani
    for NodiMeridiano in NodiOrdinatiPerMeridiani:
        n=1       
        X=NodiMeridiano[0] 
        Y=NodiMeridiano[1] 
        Z=NodiMeridiano[2] 
        XYZ=np.array([X, Y, Z])
        xyzToT=np.hstack((xyzToT,XYZ))
        ZEstradosso=np.hstack((ZEstradosso,ZNodiEstradossoOrdinatiPerMeridiani[k]))
        nodoI=0
        
        for x in X:   
            nodoJ=index(k, n, btPrgrssvLn[k])
            MeridianBranches.append((nodoI,nodoJ))
            nodoI=nodoJ
            n+=1
        #MeridianBranches.append(branches)
        k+=1
        
    return xyzToT, ZEstradosso, MeridianBranches


def RiordinaParalleli(NodiOrdinatiPerMeridiani):
    #Riordina BRANCHES PARALLELI   
    ParallelBranches=[]     
    k=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani: 
        nextIndex=k+1
        try:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[nextIndex]
        except IndexError:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[0]
            nextIndex=0
        X=NodiMeridiano[0]
        Xsuccessivo=NodiMeridianoSuccessivo[0]                     
        nmax=min(len(X), len(Xsuccessivo))
        for n in range(1, nmax+1):   
            nodoI=index(k, n, btPrgrssvLn[k]) 
            nodoJ=index(nextIndex, n, btPrgrssvLn[nextIndex])
            ParallelBranches.append((nodoI,nodoJ))
            n+=1
        k+=1        
    return ParallelBranches 

def RiordinaArcoPerimetro(NodiOrdinatiPerMeridiani):
    #RIORDINA BRANCHES ARCO DI PERIMETRO
    ArcBranches=[]
    NodiMeridiano1=NodiOrdinatiPerMeridiani[0]
    n=len(NodiMeridiano1[0])
    k=0           
    nodoI0=index(0, n, btPrgrssvLn[0])
    nodoI=nodoI0
    for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]:        
        n=len(NodiMeridiano[0])
        nodoJ=index(k+1, n, btPrgrssvLn[k+1])      
        ArcBranches.append((nodoI,nodoJ))  
        nodoI=nodoJ        
        k+=1       
    ArcBranches.append((nodoJ,nodoI0)) 
    return ArcBranches

def RiordinaFace(NodiOrdinatiPerMeridiani):
    #RIORDINA FACE
    Faces=[]
    k=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani:#[:-1]:
        nodo1=0
        nodo2=0        
        nextIndex=k+1
        try:
            NodiOrdinatiPerMeridiani[nextIndex]
        except IndexError:
            nextIndex=0

        nodiVertici=[nodo1, nodo2]
        nFaces=min(btParzLn[k], btParzLn[nextIndex])
        for n in range(0, nFaces):  
            nodo3=index(nextIndex, n+1, btPrgrssvLn[nextIndex])
            nodo4=index(k, n+1, btPrgrssvLn[k])
            nodiVertici.append(nodo3) 
            nodiVertici.append(nodo4)
            Faces.append(nodiVertici)
            nodo1=nodo4
            nodo2=nodo3
            nodiVertici=[nodo1, nodo2] 
        n+=1    
        if n+1<btParzLn[nextIndex]+1:
            nodo3=index(nextIndex, n+1, btPrgrssvLn[nextIndex]) 
        else:
            nodo3=nodo2                              
        nodiVertici.append(nodo3)    
            
        if n+1<btParzLn[k]+1:
            nodo4=index(k, n+1, btPrgrssvLn[k])                 
        else:
            nodo4=nodo1       
        nodiVertici.append(nodo4)  

        Faces.append(nodiVertici)             
        k+=1
    return Faces



def RiordinoTrueArchi(NodiArchi,k0,n0):
    k=k0
    nodoJ=n0+1  
    nodoI=index(k, btParzLn[k], btPrgrssvLn[k])      
    Nnod=len(NodiArchi) 
    TrueArcBranches1=[(nodoI, nodoJ)]
    TrueArcBranches2=[]
    k+=1
    for n in range(1, Nnod):        
        nodoI=index(k, btParzLn[k], btPrgrssvLn[k])
        nodoJ += 1 
        TrueArcBranches1.append((nodoI,nodoJ)) 
        TrueArcBranches2.append((nodoJ-1, nodoJ))        
        k+=1
    TrueArcBranches=TrueArcBranches1+TrueArcBranches2
    newIndex=n0+len(TrueArcBranches1)
    return TrueArcBranches, newIndex

def RiordinoFalseArchiY(NodiArchi,n0, n1):         
    Nnod=len(NodiArchi) 
    nodoI = n1
    nodoJ = n0
    FalseArcBranches=[]
    for n in range(0, Nnod):        
        nodoI += 1
        nodoJ += 1 
        FalseArcBranches.append((nodoI,nodoJ))     
    newIndex=n0+len(FalseArcBranches)
    return FalseArcBranches, newIndex


def RiordinoFalseArchiX(NodiArchi,k0,n0):   
    k=k0
    nodoJ=n0       
    Nnod=len(NodiArchi) 
    FalseArcBranches=[]
    for n in range(0, Nnod):        
        nodoI=index(k, btParzLn[k], btPrgrssvLn[k])
        nodoJ += 1 
        FalseArcBranches.append((nodoI,nodoJ))     
        k+=1
    newIndex=n0+len(FalseArcBranches)
    return FalseArcBranches, newIndex


def AddArchiFace(NodiArchi,k0,n0):
    ArchiFaces=[]
    
    k=k0
    nodoJ=n0+1  
    nodoI=index(k, btParzLn[k], btPrgrssvLn[k])      
    Nnod=len(NodiArchi) 
    nodiVertici=[nodoI, nodoJ]
    nodoK = nodoJ
    k+=1
    for n in range(1, Nnod):        
        nodoL=index(k, btParzLn[k], btPrgrssvLn[k])
        nodoK += 1 
        nodiVertici.append(nodoK) 
        nodiVertici.append(nodoL) 
        ArchiFaces.append(nodiVertici)        
        nodoI=nodoL
        nodoJ=nodoK
        nodiVertici=[nodoI, nodoJ] 
        k+=1
    return ArchiFaces



fileName="Vela"

#DISEGNO DI TUTTI MERIDIALI RELATIVI A 360°
dxf_fileName=fileName+"-OUT.dxf"
dxf_file = open(dxf_fileName,"w")    
WriteIntestazioneDXF(dxf_file)

maxNumeroBeta=len(listaBetaV[0])
numeroParalleli = maxNumeroBeta-2 
plotNodoColmo(NodoColmo, dxf_file)
plotNodiMeridiani(NodiOrdinatiPerMeridianiToT360, NodoColmo, NsA, dxf_file)  
ParallelBranches=plotParalleli(NodiOrdinatiPerMeridianiToT360,numeroParalleli, dxf_file)
plotArcoPerimetro(NodiOrdinatiPerMeridianiToT360, dxf_file)
plotFace(NodiOrdinatiPerMeridianiToT360, NodoColmo, dxf_file)

G=NsA-1

plotTrueArchi(TrueArcNodesY, CoordinateNodidArcoPerimetro0, dxf_file)
plotTrueArchi(TrueArcNodesY180, CoordinateNodidArcoPerimetro180, dxf_file)


 
plotFalseArchi(FalseArcNodesY, TrueArcNodesY[1:-1], dxf_file)
plotFalseArchi(FalseArcNodesY180, TrueArcNodesY180[1:-1], dxf_file)
plotFalseArchi(FalseArcNodesX90[1:-1], CoordinateNodidArcoPerimetro90[1:-1], dxf_file)
plotFalseArchi(FalseArcNodesX270[1:-1], CoordinateNodidArcoPerimetro270[1:-1], dxf_file)

abc=CoordinateNodidArcoPerimetro90[-1]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x-0.5, y+0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=CoordinateNodidArcoPerimetro90[0]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x+0.5, y+0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=CoordinateNodidArcoPerimetro270[0]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x-0.5, y-0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=CoordinateNodidArcoPerimetro270[-1]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x+0.5, y-0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")




#**********************************************
abc=TrueArcNodesY[0]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x, y-0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=TrueArcNodesY[-1]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x, y+0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=TrueArcNodesY180[0]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x, y+0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

abc=TrueArcNodesY180[-1]
x=abc[0]
y=abc[1] 
z=abc[2] 
P0=Point(x, y, z)

P1=Point(x, y-0.5, z-10.00) 
LineToDXF(dxf_file, P0, P1, "ang") 
PointToDXF(dxf_file, P1, "ang")

#*************************************

CloseWrittenDXF(dxf_file)
# 
# =============================================================================

#RIORDINO NUMERAZIONI

btPrgrssvLn, btParzLn = betaProgressiveLen(NodiOrdinatiPerMeridianiToT360)
xyzNodi, ZEstradosso, Mrdnbranches = RiordinaNodiMeridiani(NodiOrdinatiPerMeridianiToT360, btPrgrssvLn)
Prlllbranches=RiordinaParalleli(NodiOrdinatiPerMeridianiToT360)
ArcBranches=RiordinaArcoPerimetro(NodiOrdinatiPerMeridianiToT360)

lastIndex=len(xyzNodi[0])-1
lastTrueVoltaIndex=lastIndex
xyzNodi=np.hstack((xyzNodi,TrueArcNodesY.T, TrueArcNodesY180.T))
xyzNodi=np.hstack((xyzNodi,FalseArcNodesY.T, FalseArcNodesY180.T))
xyzNodi=np.hstack((xyzNodi,FalseArcNodesX90[1:-1].T, FalseArcNodesX270[1:-1].T))

TrueArcNodesYZEstradosso=np.copy(TrueArcNodesY)

ZTrueArcNodesYEstradosso=np.sqrt((r+spessore)**2-(np.square(TrueArcNodesY.T[1])))
ZTrueArcNodesYEstradosso180=np.sqrt((r+spessore)**2-(np.square(TrueArcNodesY180.T[1])))
ZEstradosso=np.hstack((ZEstradosso, ZTrueArcNodesYEstradosso, ZTrueArcNodesYEstradosso180))

Faces=RiordinaFace(NodiOrdinatiPerMeridianiToT360)
ArchiFace = AddArchiFace(TrueArcNodesY,0,lastIndex)

allBranches=np.vstack((Mrdnbranches,Prlllbranches, ArcBranches))


TrueArcBranches, lastIndex = RiordinoTrueArchi(TrueArcNodesY,0,lastIndex)
allBranches=np.vstack((allBranches,TrueArcBranches))

k0=4*(NsA-1)
ArchiFace180 = AddArchiFace(TrueArcNodesY,k0,lastIndex)
allFaces=np.vstack((Faces, ArchiFace, ArchiFace180))

TrueArcBranches, lastIndex = RiordinoTrueArchi(TrueArcNodesY180,k0,lastIndex)
allBranches=np.vstack((allBranches,TrueArcBranches))

NumeroBranchesInterni = len(allBranches)


lastTrueNodesIndex=lastIndex
n1=lastIndex-(2*NsA-1)-(2*NsA-1)+1
FalseArcBranches, lastIndex = RiordinoFalseArchiY(FalseArcNodesY, lastIndex, n1 )
allBranches=np.vstack((allBranches, FalseArcBranches))
n1=lastIndex-(2*NsA-3)-(2*NsA-1)+1
FalseArcBranches, lastIndex = RiordinoFalseArchiY(FalseArcNodesY180, lastIndex, n1)
allBranches=np.vstack((allBranches, FalseArcBranches))

k0=2*NsA-1
FalseArcBranches, lastIndex = RiordinoFalseArchiX(FalseArcNodesX90[1:-1],k0,lastIndex)
allBranches=np.vstack((allBranches, FalseArcBranches))

k0=6*NsA-5
FalseArcBranches, lastIndex = RiordinoFalseArchiX(FalseArcNodesX270[1:-1],k0,lastIndex)
allBranches=np.vstack((allBranches, FalseArcBranches))


P=CoordinateNodidArcoPerimetro90[0]
NodesAdd=np.array([P[0]+0.5, P[1]+0.5, z-10.00]) 

k=2*(NsA-1)
nodoI = index(k, btParzLn[k], btPrgrssvLn[k])
nodoJ = lastIndex+1
branch=[nodoI, nodoJ]
NewBranches=[branch]


P=CoordinateNodidArcoPerimetro90[-1]
NodesAddi=np.array([P[0]-0.5, P[1]+0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))
k=4*(NsA-1)
nodoI = index(k, btParzLn[k], btPrgrssvLn[k])
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]


P=CoordinateNodidArcoPerimetro270[0]
NodesAddi=np.array([P[0]-0.5, P[1]-0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))

k=6*(NsA-1)
nodoI = index(k, btParzLn[k], btPrgrssvLn[k])
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]

P=CoordinateNodidArcoPerimetro270[-1]
NodesAddi=np.array([P[0]+0.5, P[1]-0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))

k=0
nodoI = index(k, btParzLn[k], btPrgrssvLn[k])
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]

P=TrueArcNodesY[0]
NodesAddi=np.array([P[0], P[1]-0.5, z-10.00])
NodesAdd=np.vstack((NodesAdd,NodesAddi))
nodoI =lastTrueVoltaIndex+1
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]

P=TrueArcNodesY[-1]
NodesAddi=np.array([P[0], P[1]+0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))
nodoI =lastTrueVoltaIndex+2*NsA-1
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]

P=TrueArcNodesY180[0]
NodesAddi=np.array([P[0],P[1]+0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))
nodoI += 1
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]


P=TrueArcNodesY180[-1]
NodesAddi=np.array([P[0], P[1]-0.5, z-10.00]) 
NodesAdd=np.vstack((NodesAdd,NodesAddi))
nodoI += 2*NsA-2
nodoJ += 1
branch=[nodoI, nodoJ]
NewBranches+=[branch]

#******************************************************************************
xyzNodi=np.hstack((xyzNodi,NodesAdd.T))
allBranches=np.vstack((allBranches, NewBranches))

NumeroNodiTotali = len(xyzNodi[0])
NumeroNodiInterni = lastTrueNodesIndex+1
NumeroNodiEsterni = NumeroNodiTotali- NumeroNodiInterni

NumeroBranchesTotali = len(allBranches)
NumeroBranchesEsterni = NumeroBranchesTotali - NumeroBranchesInterni

NumeroFacesTotali  = len(allFaces)

EtichettaNodo=[]
for n in range(0, NumeroNodiInterni):
    EtichettaNodo += ["NODO INTRADOSSO INTERNO"]
    
for n in range(0, NumeroNodiEsterni):
    EtichettaNodo += ["NODO ESTERNO ZMINIMO"]

EtichettaBranche=[]
for n in range(0, NumeroBranchesInterni):
    EtichettaBranche += ["BRANCHE INTERNO"]
    
for n in range(0, NumeroBranchesEsterni):
    EtichettaBranche += ["BRANCHE ESTERNO"]   


print("\nNumero Nodi Totali = ", NumeroNodiTotali)
print("Numero Nodi Interni = ", NumeroNodiInterni)
print("Numero Nodi Esterni = ", NumeroNodiEsterni)

print("\nNumero Branches Totali = ", NumeroBranchesTotali)
print("Numero Branches Interni = ", NumeroBranchesInterni)
print("Numero Branches Esterni = ", NumeroBranchesEsterni)

print("\nNumero faces Totali = ", NumeroFacesTotali, "\n")

dxf_fileName=fileName+"-OUT_RENUMBERED.dxf"
dxf_file = open(dxf_fileName,"w")    
WriteIntestazioneDXF(dxf_file)

NumeroNodi=NumeroNodiTotali
#layerNodiIntra="NODI-Intradosso"
layerNodiEstra="NODI-Estradosso"
layerNodiText="NODI-Index"
for n in range(0, NumeroNodiTotali):         
    P=Point(xyzNodi[0][n],xyzNodi[1][n],xyzNodi[2][n])
    PointToDXF(dxf_file, P, EtichettaNodo[n])
    testo=str(n)
    TextToDXF(dxf_file, P,testo,0.002,layerNodiText)
# =============================================================================
#     P=Point(xyzNodi[0][n],xyzNodi[1][n],ZEstradosso[n])
#     PointToDXF(dxf_file, P, layerNodiEstra)
# 
# =============================================================================
NumeroallBranches=len(allBranches)
layerallBranchesText="BRANCHES-Index"
for n in range(0, NumeroallBranches): 
    branch=allBranches[n]
    iNode=branch[0]    
    jNode=branch[1]
    Pi=Point(xyzNodi[0][iNode],xyzNodi[1][iNode],xyzNodi[2][iNode])
    Pj=Point(xyzNodi[0][jNode],xyzNodi[1][jNode],xyzNodi[2][jNode])
    LineToDXF(dxf_file, Pi, Pj, EtichettaBranche[n])
    testo=str(n)
    Pm=MediaP([Pi,Pj])
    TextToDXF(dxf_file, Pm,testo,0.002,layerallBranchesText)        

#faces
NumeroFaces=len(allFaces)
layerFaces="Faces"
layerFacesText="Faces-Index"
for n in range(0, NumeroFaces): 
    face=allFaces[n]
    iNode=face[0]    
    jNode=face[1]
    kNode=face[2]    
    lNode=face[3]    
    
    Pi=Point(xyzNodi[0][iNode],xyzNodi[1][iNode],xyzNodi[2][iNode])
    Pj=Point(xyzNodi[0][jNode],xyzNodi[1][jNode],xyzNodi[2][jNode])
    Pk=Point(xyzNodi[0][kNode],xyzNodi[1][kNode],xyzNodi[2][kNode])
    Pl=Point(xyzNodi[0][lNode],xyzNodi[1][lNode],xyzNodi[2][lNode])    
    
    FaceToDXF(dxf_file, Pi, Pj, Pk, Pl, layerFaces)
    testo=str(n)
    Pm=MediaP([Pi, Pj, Pk, Pl])
    TextToDXF(dxf_file, Pm,testo,0.002,layerFacesText) 
       
CloseWrittenDXF(dxf_file)