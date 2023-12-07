# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:50:29 2023

@author: lucky
"""
from DxfUtility import *

class Point:
  def __init__(self, x, y, z):
      self.x = x
      self.y = y
      self.z = z

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

def plotNodoColmo(NodoColmo, dxf_file):
    #PLOT NODO DI COLMO
    xo=NodoColmo[0]
    yo=NodoColmo[1]
    zo=NodoColmo[2] 
    Po=Point(xo,yo,zo)
    PointToDXF(dxf_file, Po,"ColmoNodes")      
 
def plotNodiMeridiani(NodiOrdinatiPerMeridiani, NodoColmo, NsA, dxf_file):
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
        for X in x:              
            P=Point(x[n],y[n],z[n])
            PointToDXF(dxf_file, P, layerN)
            LineToDXF(dxf_file,  PP, P, layerL)
            n+=1
            PP=P                    
        k+=1
        layerL = "MeridianiBranches"
        layerN = "NodiInterni"

   
def plotParalleli(NodiOrdinatiPerMeridiani,numeroParalleli, dxf_file):
    #PLOT BRANCHES PARALLELI      
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
        k+=1
      
        
def plotArcoPerimetro(NodiOrdinatiPerMeridiani, dxf_file):
    #PLOT ARCO DI PERIMETRO
    NodiMeridiano1=NodiOrdinatiPerMeridiani[0]
    x1=NodiMeridiano1[0][-1]
    y1=NodiMeridiano1[1][-1] 
    z1=NodiMeridiano1[2][-1] 

    P1=Point(x1, y1, z1)
    k=0
    for NodiMeridiano in NodiOrdinatiPerMeridiani[1:]:  

        x2=NodiMeridiano[0][-1]
        y2=NodiMeridiano[1][-1] 
        z2=NodiMeridiano[2][-1] 
        P2=Point(x2, y2, z2)
        LineToDXF(dxf_file, P1, P2, "ArcoBranches")          

        P1=P2
        k+=1
    P2=Point(x1, y1, z1)    
    LineToDXF(dxf_file, P1, P2, "ArcoBranches")

#adeguare per PLOT solo nodi
def plotNodiEstradosso(NodiOrdinatiPerMeridianiEstradosso, dxf_file):
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


def plotFace(NodiOrdinatiPerMeridiani, NodoColmo, dxf_file):
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
        x=NodiMeridiano[0]
        y=NodiMeridiano[1]   
        z=NodiMeridiano[2]    
        xsuccessivo=NodiMeridianoSuccessivo[0]
        ysuccessivo=NodiMeridianoSuccessivo[1]   
        zsuccessivo=NodiMeridianoSuccessivo[2]   
        nFaces=max(x.size, xsuccessivo.size)
        for n in range(0, nFaces):
            try:
                P3=Point(xsuccessivo[n],ysuccessivo[n],zsuccessivo[n])
            except IndexError:
                P3=P2
            try:    
                P4=Point(x[n],y[n],z[n])
            except IndexError:
                P4=P1                
            FaceToDXF(dxf_file, P1, P2,P3,P4, layerF)
     
            #P=MediaP(P1, P2,P3,P4)
            P=MediaP([P1, P2, P3, P4])
            testo=str(k)+"."+str(n)
            TextToDXF(dxf_file, P,testo,textheight,layerT)
            #n+=1
            P1=P4
            P2=P3 

        k+=1
        
        
def plotTrueArchi(NodiArchi, NodiAppoggio, dxf_file):
    #PLOT A
    layerL = "ArchiBranches"
    layerN = "NodiEsterni"
        
    k=0
    XYZ = NodiAppoggio[k]
    P0=Point(XYZ[0], XYZ[1], XYZ[2]) 
    x1=NodiArchi[0][0]
    y1=NodiArchi[0][1] 
    z1=NodiArchi[0][2] 
    P1=Point(x1, y1, z1)
    
    PointToDXF(dxf_file, P1, layerN)
    LineToDXF(dxf_file, P0, P1, layerL)

    Nnod=len(NodiArchi)    
    for k in range(1, Nnod): 
        x2=NodiArchi[k][0]
        y2=NodiArchi[k][1] 
        z2=NodiArchi[k][2] 
        P2=Point(x2, y2, z2)
        PointToDXF(dxf_file, P2, layerN)
        LineToDXF(dxf_file, P1, P2, layerL)           
        XYZ = NodiAppoggio[k]
        P0=Point(XYZ[0], XYZ[1], XYZ[2]) 
        LineToDXF(dxf_file, P0, P2, layerL)        
        P1=P2       
      
def plotFalseArchi(NodiArchi, NodiAppoggio, dxf_file):
    #PLOT A
    layerL = "ArchiBranchesEsterni"
    layerN = "NodiEsterniEsterni"
    k=0
    for Nodo in NodiArchi: 
        XYZ = NodiAppoggio[k]
        P1=Point(XYZ[0], XYZ[1], XYZ[2])       
        x2=NodiArchi[k][0]
        y2=NodiArchi[k][1] 
        z2=NodiArchi[k][2] 
        P2=Point(x2, y2, z2)
        PointToDXF(dxf_file, P2, layerN)    

        LineToDXF(dxf_file, P1, P2, layerL)

        k+=1
        
  
