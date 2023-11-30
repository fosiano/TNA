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


NsA=7
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



alfa=0.0
beta=2*a45
NodoColmo=GeneraNodi(R,alfa,beta)
NodiOrdinatiPerMeridiani=[]
NodiOrdinatiPerMeridianiEstradosso=[]
i=0
for alfa in alfaV:
# =============================================================================
#     r=NsA-1-abs(NsA-1-i)#rivedere in funzione di NsA
#     betaV=np.insert(betaV0,0, listbetaArc[r]) 
# =============================================================================   
    betaV=listaBetaV[i]  
    NodiMeridiano=GeneraNodi(R,alfa,betaV[:-1])      
    NodiOrdinatiPerMeridiani.append(np.array(NodiMeridiano))   
    
    NodiMeridianoEstradosso=np.copy(NodiMeridiano)
    Z=GeneraNodiEstradosso(R+0.5,NodiMeridiano[0],NodiMeridiano[1])
    NodiMeridianoEstradosso[2]=Z
    NodiOrdinatiPerMeridianiEstradosso.append(np.array(NodiMeridianoEstradosso))   
    i+=1

i=0
for NodiMeridiano in NodiOrdinatiPerMeridiani:  
    NodiMeridiano[0]=np.flipud(NodiMeridiano[0])
    NodiMeridiano[1]=np.flipud(NodiMeridiano[1])
    NodiMeridiano[2]=np.flipud(NodiMeridiano[2])  
    NodiMeridianoE=NodiOrdinatiPerMeridianiEstradosso[i]    
    NodiMeridianoE[0]=np.flipud(NodiMeridianoE[0])
    NodiMeridianoE[1]=np.flipud(NodiMeridianoE[1])
    NodiMeridianoE[2]=np.flipud(NodiMeridianoE[2])      
    
    
    
NodiOrdinatiPerMeridiani180=[]
NodiOrdinatiPerMeridiani90=[]
NodiOrdinatiPerMeridiani270=[]

for NodiMeridiano in NodiOrdinatiPerMeridiani:
    NodiMeridiano180=np.copy(NodiMeridiano)
    NodiMeridiano180[0]=-1*NodiMeridiano180[0]
    NodiOrdinatiPerMeridiani180.append(NodiMeridiano180)
    
    NodiMeridiano90=np.copy(NodiMeridiano)
    NodiMeridiano90x=np.copy(NodiMeridiano[0])
    NodiMeridiano90y=np.copy(NodiMeridiano[1])
    NodiMeridiano90[0]=NodiMeridiano90y
    NodiMeridiano90[1]=NodiMeridiano90x
    NodiOrdinatiPerMeridiani90.append(NodiMeridiano90)
    
    NodiMeridiano270=np.copy(NodiMeridiano90)
    NodiMeridiano270[1]=-1*NodiMeridiano270[1]
    NodiOrdinatiPerMeridiani270.append(NodiMeridiano270)
#******************************************************************************
#SALVA IN DXF INPUT - OUTPUT
#******************************************************************************
def PointToDXF(P, layer):  
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

def LineToDXF(P1, P2, layer): 
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

def FaceToDXF(P1, P2,P3,P4, layer): 
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

def TextToDXF(P, text,textheight, layer): 
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

outfile="Vela"
dxf_file = open(outfile+"-OUT.dxf","w")
dxf_file.write("999\n")
dxf_file.write("DXF created from myself (Fortunato Siano)\n")
dxf_file.write("0\n")
dxf_file.write("SECTION\n")
dxf_file.write("2\n")
dxf_file.write("ENTITIES\n")

def plotNodoColmo(NodoColmo):
    #PLOT NODO DI COLMO
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo)
    PointToDXF(Po,"ColmoNodes")

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
        if k==1 or k==2*NsA-1:
            layerL="DiagonaliBranches"
            layerN="DiagonaliNodes"
        if k==NsA:
            layerL="ColmoBranches"
            layerN="ColmoNodes"            
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]   
        n=0
        PP=Po
        #PointToDXF(PP, layerN)
        for X in x:              
            P=Point(x[n],y[n],z[n])
            PointToDXF(P, layerN)
            LineToDXF(PP, P, layerL)
            n+=1
            PP=P                    
        #LineToDXF(PP, Po, layerL)
        k+=1
        layerL = "NodiOrdinatiPerMeridianiBranches"
        layerN = "NodiInterni"

def plotParalleli(NodiOrdinatiPerMeridiani,numeroParalleli):
    #PLOT BRANCHES PARALLELI  
    NodiMeridiano1=NodiOrdinatiPerMeridiani[0]
    for j in range(0, numeroParalleli):
        x1=NodiMeridiano1[0][j]
        y1=NodiMeridiano1[1][j] 
        z1=NodiMeridiano1[2][j] 
        P1=Point(x1, y1, z1)
        vuoto=False
        for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]: 
            dim=NodiMeridiano[0].size
            if j<dim:
                x2=NodiMeridiano[0][j]
                y2=NodiMeridiano[1][j] 
                z2=NodiMeridiano[2][j] 
                P2=Point(x2, y2, z2)
                if not vuoto: 
                    LineToDXF(P1, P2, "ParalleliBranches")    
                P1=P2
                vuoto = False
            else:
                vuoto = True
            
def plotArcoPerimetro(NodiOrdinatiPerMeridiani):
    #PLOT ARCO DI PERIMETRO
    NodiMeridiano1=NodiOrdinatiPerMeridiani[0]
    x1=NodiMeridiano1[0][-1]
    y1=NodiMeridiano1[1][-1] 
    z1=NodiMeridiano1[2][-1] 
    P1=Point(x1, y1, z1)
    #Nm=len(NodiOrdinatiPerMeridiani)
    for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]:  
        x2=NodiMeridiano[0][-1]
        y2=NodiMeridiano[1][-1] 
        z2=NodiMeridiano[2][-1] 
        P2=Point(x2, y2, z2)
        LineToDXF(P1, P2, "ArcoBranches")    
        P1=P2


#adeguare per PLOT solo nodi
def plotNodiEstradosso(NodiOrdinatiPerMeridianiEstradosso):
    #PLOT NODI NodiOrdinatiPerMeridiani
    layerN = "NodiEstradossoInterni"
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo+0.50)   
    PointToDXF(Po, "ColmoNodesEstradosso")
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
        PointToDXF(PP, layerN)
        for X in x[:-1]:   
            n+=1
            P=Point(x[n],y[n],z[n])
            PointToDXF(P, layerN)
            PP=P                    
        k+=1
        layerN = "NodiInterniEstradosso"

def MediaP(P1, P2,P3,P4):
    x=(P1.x+P2.x+P3.x+P4.x)/4
    y=(P1.y+P2.y+P3.y+P4.y)/4
    z=(P1.z+P2.z+P3.z+P4.z)/4
    P=Point(x,y,z)
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
    for NodiMeridiano in NodiOrdinatiPerMeridiani[:-1]: 
        P1=Po
        P2=Po
        NodiMeridianoSuccessivo=NodiOrdinatiPerMeridiani[k+1]
        #n1=n0
        #n2=n0
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]    
        xsuccessivo=NodiMeridianoSuccessivo[0]
        ysuccessivo=NodiMeridianoSuccessivo[1]   
        zsuccessivo=NodiMeridianoSuccessivo[2]   
        n=0
        face=[]
        for X in x[:-1]:       
            P3=Point(xsuccessivo[n],ysuccessivo[n],zsuccessivo[n])
            P4=Point(x[n],y[n],z[n])
            FaceToDXF(P1, P2,P3,P4, layerF)
            nodi=[(k,n-1),(k+1,n-1),(k+1,n),(k,n)]
            testo=str(k)+"."+str(n)
            face.append([testo,nodi])
            
            P=MediaP(P1, P2,P3,P4)
            
            TextToDXF(P,testo,textheight,layerT)
            n+=1
            P1=P4
            P2=P3 
            
            
        if len(x)>len(xsuccessivo):
            P4=Point(x[n],y[n],z[n])
            FaceToDXF(P1, P2,P3,P4, layerF)          
            nodi=[(k,n-1),(k+1,n-1),(k+1,n-1),(k,n)]
            testo=str(k)+"."+str(n)
            face.append([testo,nodi])

            P=MediaP(P1, P2,P3,P4)
            TextToDXF(P,testo,textheight,layerT)
        else:
            P3=Point(xsuccessivo[n],ysuccessivo[n],zsuccessivo[n])
            P4=Point(x[n],y[n],z[n])  
            FaceToDXF(P1,P2,P3,P4, layerF)
            
            
            testo=str(k)+"."+str(n)
            nodi=[(k,n-1),(k+1,n-1),(k+1,n),(k,n)]
            face.append([testo,nodi])
            P=MediaP(P1, P2,P3,P4)
            TextToDXF(P,testo,textheight,layerT)
            P1=Point(xsuccessivo[n+1],ysuccessivo[n+1],zsuccessivo[n+1])            
            FaceToDXF(P4,P3,P1,P4, layerF)
            nodi=[(k,n),(k+1,n),(k+1,n+1),(k+1,n+1)]
            testo=str(k)+"."+str(n)
            face.append([testo,nodi])
            P=MediaP(P4,P3,P1,P4)
            TextToDXF(P,testo,textheight,layerT)
        Faces.append(face)
        #FaceToDXF(PP, Po, layerL)
        k+=1



maxNumeroBeta=len(listaBetaV[0])
numeroParalleli = maxNumeroBeta-2 
#PLOT NODO DI COLMO
plotNodoColmo(NodoColmo)



#***********************************************************
#ROTAZIONE 0
plotNodiMeridiani(NodiOrdinatiPerMeridiani)  
plotParalleli(NodiOrdinatiPerMeridiani,numeroParalleli)
plotArcoPerimetro(NodiOrdinatiPerMeridiani)

plotFace(NodiOrdinatiPerMeridiani)


#***********************************************************
#ROTAZIONE 180
plotNodiMeridiani(NodiOrdinatiPerMeridiani180)  
plotParalleli(NodiOrdinatiPerMeridiani180,numeroParalleli)
plotArcoPerimetro(NodiOrdinatiPerMeridiani180)
plotFace(NodiOrdinatiPerMeridiani180)


#***********************************************************
#ROTAZIONE 90
plotNodiMeridiani(NodiOrdinatiPerMeridiani90[1:-1])  
plotParalleli(NodiOrdinatiPerMeridiani90,numeroParalleli)
plotArcoPerimetro(NodiOrdinatiPerMeridiani90)
plotFace(NodiOrdinatiPerMeridiani90)

#***********************************************************
#ROTAZIONE 270
plotNodiMeridiani(NodiOrdinatiPerMeridiani270[1:-1])  
plotParalleli(NodiOrdinatiPerMeridiani270,numeroParalleli)
plotArcoPerimetro(NodiOrdinatiPerMeridiani270)
plotFace(NodiOrdinatiPerMeridiani270)

plotNodiEstradosso(NodiOrdinatiPerMeridianiEstradosso)

dxf_file.write("0\n")
dxf_file.write("ENDSEC\n")
dxf_file.write("0\n")
dxf_file.write("EOF\n")

dxf_file.close()