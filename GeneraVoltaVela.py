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


class Point:
  def __init__(self, x, y, z):
      self.x = x
      self.y = y
      self.z = z
      
def prova(a):
    return "ok", a

B=3.75
r=B/2
R=r*(2.0**0.50)


NsA=9
NsB=10

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

#suddivisione degli alfa tale che gli archi sull'arco di peimetro siano di equilughezza
gammaV=np.linspace(-2*a45, 2*a45, 2*NsA-1)
alfaV=np.arctan(np.sin(gammaV))
#alfaV=np.linspace(0.0, a45, NsA)
betaV0=np.linspace(a45, 2*a45, NsB)#da 45° a 90° (90 compreso)

#inseriti=[]


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
NodiOrdinatiPerMeridianiEstradosso=[]
i=0
for alfa in alfaV:
    betaV=listaBetaV[i]       
    NodiMeridiano=GeneraNodi(R,alfa,betaV[1:])   
    NodiOrdinatiPerMeridiani.append(np.array(NodiMeridiano))   
    
    NodiMeridianoEstradosso=np.copy(NodiMeridiano)
    Z=GeneraNodiEstradosso(R+0.5,NodiMeridiano[0],NodiMeridiano[1])
    NodiMeridianoEstradosso[2]=Z
    NodiOrdinatiPerMeridianiEstradosso.append(np.array(NodiMeridianoEstradosso))   
    i+=1
 
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

#e=np.hstack((c,d))


#******************************************************************************
#SALVA IN DXF INPUT - OUTPUT
#******************************************************************************
def PointToDXF(dxf_file, P, layer):  
    dxf_file.write("0\n")
    dxf_file.write("POINT\n")
    dxf_file.write("8\n")
    dxf_file.write(layer+"\n")
    dxf_file.write("10\n")  
    dxf_file.write(str(P.x)+"\n")
    dxf_file.write("20\n")
    dxf_file.write(str(P.y)+"\n")
    dxf_file.write("30\n")
    dxf_file.write(str(P.z)+"\n")

def LineToDXF(dxf_file, P1, P2, layer): 
    dxf_file.write("0\n")
    dxf_file.write("LINE\n")
    dxf_file.write("8\n")
    dxf_file.write(layer+"\n")
    dxf_file.write("10\n")  
    dxf_file.write(str(P1.x)+"\n")
    dxf_file.write("20\n")
    dxf_file.write(str(P1.y)+"\n")
    dxf_file.write("30\n")
    dxf_file.write(str(P1.z)+"\n")
    dxf_file.write("11\n")
    dxf_file.write(str(P2.x)+"\n")
    dxf_file.write("21\n")
    dxf_file.write(str(P2.y)+"\n")
    dxf_file.write("31\n")
    dxf_file.write(str(P2.z)+"\n")    

def FaceToDXF(dxf_file, P1, P2,P3,P4, layer): 
    dxf_file.write("0\n")
    dxf_file.write("3DFACE\n")
    dxf_file.write("8\n")
    dxf_file.write(layer+"\n")
    dxf_file.write("10\n")  
    dxf_file.write(str(P1.x)+"\n")
    dxf_file.write("20\n")
    dxf_file.write(str(P1.y)+"\n")
    dxf_file.write("30\n")
    dxf_file.write(str(P1.z)+"\n")
    dxf_file.write("11\n")
    dxf_file.write(str(P2.x)+"\n")
    dxf_file.write("21\n")
    dxf_file.write(str(P2.y)+"\n")
    dxf_file.write("31\n")
    dxf_file.write(str(P2.z)+"\n")    
    dxf_file.write("12\n")
    dxf_file.write(str(P3.x)+"\n")
    dxf_file.write("22\n")
    dxf_file.write(str(P3.y)+"\n")
    dxf_file.write("32\n")
    dxf_file.write(str(P3.z)+"\n") 
    dxf_file.write("13\n")
    dxf_file.write(str(P4.x)+"\n")
    dxf_file.write("23\n")
    dxf_file.write(str(P4.y)+"\n")
    dxf_file.write("33\n")
    dxf_file.write(str(P4.z)+"\n")    

def TextToDXF(dxf_file, P, text,textheight, layer): 
        dxf_file.write("0\n")
        dxf_file.write("TEXT\n")
        dxf_file.write("8\n")
        dxf_file.write(layer+"\n")     
        dxf_file.write("10\n")  
        dxf_file.write(str(P.x)+"\n")
        dxf_file.write("20\n")
        dxf_file.write(str(P.y)+"\n")
        dxf_file.write("30\n")
        dxf_file.write(str(P.z)+"\n")
        dxf_file.write("40\n")
        dxf_file.write(str(textheight)+"\n")
        dxf_file.write("1\n")
        dxf_file.write(text+"\n")
        
        dxf_file.write("72\n")
        dxf_file.write("1\n")        
        dxf_file.write("11\n")  
        dxf_file.write(str(P.x)+"\n")
        dxf_file.write("21\n")
        dxf_file.write(str(P.y)+"\n")
        dxf_file.write("31\n")
        dxf_file.write(str(P.z)+"\n")
        dxf_file.write("73\n")
        dxf_file.write("2\n")


def plotNodoColmo(NodoColmo):
    #PLOT NODO DI COLMO
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo)
    PointToDXF(dxf_file, Po,"ColmoNodes")

def index(k, n, betaProgressiveLen):
    indx=betaProgressiveLen[k]+n
    return indx

def betaProgressiveLen(NodiOrdinatiPerMeridiani):
    dimension=[0]
    dimProgressive=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani:     
        dimProgressive+=len(NodiMeridiano[0])
        dimension.append(dimProgressive)
    return dimension


def RiordinaNodiMeridiani(NodiOrdinatiPerMeridiani, betaProgressiveLen):
    #PLOT NODI E BRANCHES NodiOrdinatiPerMeridiani
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2]         
    xyzToT=np.array([[xo], [yo], [zo]])
    k=0
    MeridianBranches=[]
    for NodiMeridiano in NodiOrdinatiPerMeridiani:
        n=1       
        X=NodiMeridiano[0] 
        Y=NodiMeridiano[1] 
        Z=NodiMeridiano[2] 
        XYZ=np.array([X, Y, Z])
        xyzToT=np.hstack((xyzToT,XYZ))
        nodoI=0
        for x in X:   
            nodoJ=index(k, n, btPrgrssvLn)
            MeridianBranches.append((nodoI,nodoJ))
            nodoI=nodoJ
            n+=1
        #MeridianBranches.append(branches)
        k+=1
        
    return xyzToT, MeridianBranches


MeridianBranches=[]
def plotNodiMeridiani(NodiOrdinatiPerMeridiani):
    #PLOT NODI E BRANCHES NodiOrdinatiPerMeridiani
    layerL = "MeridianiBranches"
    layerN = "NodiInterni"
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo)    
    k=1
    
    for NodiMeridiano in NodiOrdinatiPerMeridiani: 
        if k==1 or k==2*NsA-1 or k==4*NsA-3 or k==6*NsA-5:
            layerL="DiagonaliBranches"
            layerN="DiagonaliNodes"
        if k==NsA or k==3*NsA-2 or k==5*NsA-4 or k==7*NsA-6:
            layerL="ColmoBranches"
            layerN="ColmoNodes"            
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]   
        n=0
        PP=Po
        nodoI=(k-1,-1)
        #PointToDXF(PP, layerN)
        branches=[]
        for X in x:              
            P=Point(x[n],y[n],z[n])
            PointToDXF(dxf_file, P, layerN)
            LineToDXF(dxf_file,  PP, P, layerL)
            nodoJ=(k-1,n)
            branches.append([(k-1,n),[nodoI,nodoJ]])
            nodoI=nodoJ
            n+=1
            PP=P                    
        #LineToDXF(PP, Po, layerL)
        MeridianBranches.append(branches)
        k+=1
        layerL = "MeridianiBranches"
        layerN = "NodiInterni"

   
def plotParalleli(NodiOrdinatiPerMeridiani,numeroParalleli):
    #PLOT BRANCHES PARALLELI   
    ParallelBranches=[]     
    k=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani: 
        nextIndex=k+1
        try:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[nextIndex]
        except IndexError:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[0]
            nextIndex=0
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]    
        xsuccessivo=NodiMeridianoSuccessivo[0]
        ysuccessivo=NodiMeridianoSuccessivo[1]   
        zsuccessivo=NodiMeridianoSuccessivo[2] 
        branches=[]
        for n in range(0, numeroParalleli):
            try:
                P2=Point(xsuccessivo[n],ysuccessivo[n],zsuccessivo[n])
            except IndexError:
                break
            try:    
                P1=Point(x[n],y[n],z[n])
            except IndexError:
                break      
            LineToDXF(dxf_file, P1, P2, "ParalleliBranches") 
            nodoI=(k,n)
            nodoJ=(nextIndex,n)
            branches.append([(k,n),[nodoI,nodoJ]])
        ParallelBranches.append(branches)       
        k+=1
    return ParallelBranches        
        
def plotArcoPerimetro(NodiOrdinatiPerMeridiani):
    #PLOT ARCO DI PERIMETRO
    ArcBranches=[]
    NodiMeridiano1=NodiOrdinatiPerMeridiani[0]
    x1=NodiMeridiano1[0][-1]
    y1=NodiMeridiano1[1][-1] 
    z1=NodiMeridiano1[2][-1] 
    n=len(NodiMeridiano1[0])-1
    P1=Point(x1, y1, z1)
    #Nm=len(NodiOrdinatiPerMeridiani)
    #branches=[]
    k=0
    nodoJ=(0,n)
    for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]:  
        nodoI=nodoJ
        x2=NodiMeridiano[0][-1]
        y2=NodiMeridiano[1][-1] 
        z2=NodiMeridiano[2][-1] 
        n=len(NodiMeridiano[0])-1
        nodoJ=(k+1,n)
        P2=Point(x2, y2, z2)
        LineToDXF(dxf_file, P1, P2, "ArcoBranches")          
        ArcBranches.append([(k,n),[nodoI,nodoJ]])  
        P1=P2
        k+=1
    P2=Point(x1, y1, z1)    
    LineToDXF(dxf_file, P1, P2, "ArcoBranches")
    return ArcBranches

#adeguare per PLOT solo nodi
def plotNodiEstradosso(NodiOrdinatiPerMeridianiEstradosso):
    #PLOT NODI NodiOrdinatiPerMeridiani
    layerN = "NodiEstradossoInterni"
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo+0.50)   
    PointToDXF(dxf_file, Po, "ColmoNodesEstradosso")
    k=1
    for NodiMeridiano in NodiOrdinatiPerMeridianiEstradosso: 
        if k==1 or k==2*NsA-1:
            layerN="DiagonaliNodesEstradosso"
        if k==NsA:
            layerN="ColmoNodesEstradosso"            
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]   
        n=0
        PP=Point(x[0],y[0],z[0])
        PointToDXF(dxf_file, PP, layerN)
        for X in x[:-1]:   
            n+=1
            P=Point(x[n],y[n],z[n])
            PointToDXF(dxf_file, P, layerN)
            PP=P                    
        k+=1
        layerN = "NodiInterniEstradosso"

# =============================================================================
# def MediaP(P1, P2,P3,P4):
#     x=(P1.x+P2.x+P3.x+P4.x)/4
#     y=(P1.y+P2.y+P3.y+P4.y)/4
#     z=(P1.z+P2.z+P3.z+P4.z)/4
#     P=Point(x,y,z)
#     return P
# =============================================================================


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

Faces=[]

def plotFace(NodiOrdinatiPerMeridiani):
    #PLOT FACE
    layerF = "Faces"
    layerT = "FacesName"
    textheight=0.02
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo)    
    k=0
    #n0=-1   
    for NodiMeridiano in NodiOrdinatiPerMeridiani:#[:-1]: 
        P1=Po
        P2=Po
        nextIndex=k+1
        try:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[nextIndex]
        except IndexError:
            NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[0]
            nextIndex=0
        #n1=n0
        #n2=n0
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]    
        xsuccessivo=NodiMeridianoSuccessivo[0]
        ysuccessivo=NodiMeridianoSuccessivo[1]   
        zsuccessivo=NodiMeridianoSuccessivo[2]   
        #n=0
        nodiVertici=[(k,-1)]
        nFaces=max(x.size, xsuccessivo.size)
        face=[]
        for n in range(0, nFaces):
            try:
                P3=Point(xsuccessivo[n],ysuccessivo[n],zsuccessivo[n])
                nodiVertici.append((nextIndex,n))
            except IndexError:
                P3=P2
            try:    
                P4=Point(x[n],y[n],z[n])
                nodiVertici.append((k,n))
            except IndexError:
                P4=P1                
            FaceToDXF(dxf_file, P1, P2,P3,P4, layerF)

            face.append([(k,n),nodiVertici])           
            #P=MediaP(P1, P2,P3,P4)
            P=MediaP([P1, P2, P3, P4])
            testo=str(k)+"."+str(n)
            TextToDXF(dxf_file, P,testo,textheight,layerT)
            #n+=1
            P1=P4
            P2=P3 
            nodiVertici=[(k,n),(nextIndex,n)]
        Faces.append(face)
        #FaceToDXF(PP, Po, layerL)
        k+=1


outfile="Vela"
dxf_file = open(outfile+"-OUT.dxf","w")
dxf_file.write("999\n")
dxf_file.write("DXF created from myself (Fortunato Siano)\n")
dxf_file.write("0\n")
dxf_file.write("SECTION\n")
dxf_file.write("2\n")
dxf_file.write("ENTITIES\n")

maxNumeroBeta=len(listaBetaV[0])
numeroParalleli = maxNumeroBeta-2 
#PLOT NODO DI COLMO
plotNodoColmo(NodoColmo)

#DISEGNO DI TUTTI MERIDIALI RELATIVI A 360°
plotNodiMeridiani(NodiOrdinatiPerMeridianiToT360)  
ParallelBranches=plotParalleli(NodiOrdinatiPerMeridianiToT360,numeroParalleli)
ArcBranches=plotArcoPerimetro(NodiOrdinatiPerMeridianiToT360)
plotFace(NodiOrdinatiPerMeridianiToT360)
# 
# =============================================================================

plotNodiEstradosso(NodiOrdinatiPerMeridianiEstradosso)

dxf_file.write("0\n")
dxf_file.write("ENDSEC\n")
dxf_file.write("0\n")
dxf_file.write("EOF\n")

dxf_file.close()

btPrgrssvLn=betaProgressiveLen(NodiOrdinatiPerMeridianiToT360)
print(index(0, 1, btPrgrssvLn))

xyzNodi, Mrdnbranches = RiordinaNodiMeridiani(NodiOrdinatiPerMeridianiToT360, btPrgrssvLn)


outfile="Velarenum"
dxf_file = open(outfile+"-OUT.dxf","w")
dxf_file.write("999\n")
dxf_file.write("DXF created from myself (Fortunato Siano)\n")
dxf_file.write("0\n")
dxf_file.write("SECTION\n")
dxf_file.write("2\n")
dxf_file.write("ENTITIES\n")

NumeroNodi=len(xyzNodi[0])
layerNodi="NODI-Intradosso"
layerNodiText="NODI-Index"
for n in range(0, NumeroNodi):         
    P=Point(xyzNodi[0][n],xyzNodi[1][n],xyzNodi[2][n])
    PointToDXF(dxf_file, P, layerNodi)
    testo=str(n)
    TextToDXF(dxf_file, P,testo,0.002,layerNodiText)
    
    
NumeroMrdnbranches=len(Mrdnbranches)
layerMrdnbranches="BRNCHES-Meridian"
layerMrdnbranchesText="BRNCHES-Index"
for n in range(0, NumeroMrdnbranches): 
    branch=Mrdnbranches[n]
    iNode=branch[0]    
    jNode=branch[1]
    Pi=Point(xyzNodi[0][iNode],xyzNodi[1][iNode],xyzNodi[2][iNode])
    Pj=Point(xyzNodi[0][jNode],xyzNodi[1][jNode],xyzNodi[2][jNode])
    LineToDXF(dxf_file, Pi, Pj, layerMrdnbranches)
    testo=str(n)
    
    #calcolare Pmedio
    Pm=MediaP([Pi,Pj])
    TextToDXF(dxf_file, Pm,testo,0.002,layerMrdnbranchesText)    
    

dxf_file.write("0\n")
dxf_file.write("ENDSEC\n")
dxf_file.write("0\n")
dxf_file.write("EOF\n")
dxf_file.close()
