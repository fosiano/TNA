#FILTRO .DXF TO .DAT BY FORTUNATO SIANO 16.09.2017
#CONVERSIONE FILE DI UN NETWORK CREATO IN DXF AL FILE DI INPUT PER LA PROCEDURA "ROSATI"
#IL DXF CONTIENE I BRANCHES CHE COLLEGANO I NODI ALLE QUTE ZMIN
#IL DXF CONTIENE I PUNTI CHE RIPORTANO LE ZMAX
from math import *
from decimal import *
import time
import cmath
#from numpy import *
from array import *
import numpy  as np

layerExternalBranches="EXT"#DA IMPLEMENTARE
file="RosatiGroinVault"

dminMeridiani=0.5
dminDiagonali=1.0
dminParalleli=0.75

tempo_iniziale = time.time()
print ("Elaborazione file " + file + ".dxf ....")

dxf_file = open((file+".dxf"),"r")

dat_file = open((file+".dat"),"w")

line=dxf_file.readline()
while (not (line == "ENTITIES\n")):
      line=dxf_file.readline()

branches=np.empty(shape=(0,2),dtype=np.uint16)
branchesdmin=np.empty(shape=(0,1))
branchesINnode=[]
LayerbranchesINnode=[] #IMPORTANTE DA FARE
Nn=-1
Nb=-1
xyz=[]
Nnn=-1
xyzmax=[]

NodesLayer=[]
flagLine=0
flagPoint=0
line=dxf_file.readline()
while (not (line == "EOF\n")):  
    if line == "LINE\n":
        flagLine=1
        Nb+=1  
        #print "LINE"
    if line == "POINT\n":
        flagPoint=1
        #print "POINT"
        
    Vc = line.split()
    if ((Vc[0] == "8")and(flagLine==1)): 
        layer=dxf_file.readline()
        if (layer=="MERIDIANI\n")or(layer=="EXT\n"):
            q=np.array([[dminMeridiani]])
        else:
            if layer=="DIAGONALI\n":
                q=np.array([[dminDiagonali]])
            else:
               q=np.array([[dminParalleli]]) 
        branchesdmin= np.append(branchesdmin,q,axis=0)
    if ((Vc[0] == "10")and(flagLine==1)):                        
        x1=np.around(float(dxf_file.readline()),9)
        dxf_file.readline()
        y1=np.around(float(dxf_file.readline()),9)
        dxf_file.readline()
        z1=np.around(float(dxf_file.readline()),9)
        
        xyz1=[x1,y1,z1]
        #NodesLayer1=[layer.split()]
        lyr=layer.split()[0]
        try:
            n = xyz.index(xyz1)
            branchesINnode[n]+=[Nb]
            LayerbranchesINnode[n]+=[lyr]
        except ValueError:
            Nn+=1
            n=Nn
            xyz.append(xyz1)
            NodesLayer.append(layer.split()[0])           
            branchesINnode.append([Nb])
            LayerbranchesINnode.append([lyr])           

        dxf_file.readline()
        x2=np.around(float(dxf_file.readline()),9)    
        dxf_file.readline()
        y2=np.around(float(dxf_file.readline()),9)
        dxf_file.readline()
        z2=np.around(float(dxf_file.readline()),9)
        xyz2=[x2,y2,z2]
        
        try:
            m = xyz.index(xyz2)
            branchesINnode[m]+=[Nb]
            LayerbranchesINnode[m]+=[lyr]
        except ValueError:
            Nn+=1
            m=Nn
            xyz.append(xyz2)
            NodesLayer.append(layer.split()[0])
            branchesINnode.append([Nb])
            LayerbranchesINnode.append([lyr]) 
        a=np.array([[n,m]])   
        branches= np.append(branches,a,axis=0)
        flagLine=0

#PER POINT ZMAX   
    if ((Vc[0] == "10")and(flagPoint==1)):                        
        xm=np.around(float(dxf_file.readline()),9)
        dxf_file.readline()
        ym=np.around(float(dxf_file.readline()),9)
        dxf_file.readline()
        zm=np.around(float(dxf_file.readline()),9)
        xyzm=[xm,ym,zm] 
        try:
            oo = xyzmax.index(xyzm)
            print (xyzm)
        except ValueError:
            Nnn+=1
            xyzmax.append(xyzm)
        flagPoint=0
#END PER POINT ZMAX        
    line=dxf_file.readline()
dxf_file.close()

Nn=Nn+1
Nb=Nb+1
Nnn=Nnn+1

def ordinaXY(elem):
    return [elem[0],elem[1]]

xyzord=[[xyz[i][0],xyz[i][1],xyz[i][2],i] for i in range(0,Nn)]
xyzord=sorted(xyzord, key=ordinaXY)
xyzmax=sorted(xyzmax, key=ordinaXY)
xyz=[[xyzord[i][0],xyzord[i][1],xyzord[i][2],xyzmax[i][2],xyzord[i][3]] for i in range(0,Nn)]

#########xyzMinMax=[[xyzord[i].append( xyzmax[i][2])] for i in range(0,Nn)]

def ordinaIndex(elem):
    return elem[4]
xyz=sorted(xyz, key=ordinaIndex)

##########[xyzord[i][0:2] for i in range (0,Nn)]

#PER IL CALCOLO DELL'AREA INCIDENTE SUL NODO
#def Anodo(nodo):

        
     
#    return A
'''        
        rop=ro[nb]
        
        
        alfap=alfa[nb]
        i=1
        while (i<nbn):
            nb=orderIndex[i]
            A+=ro[nb]*rop*np.sin(alfa[nb]-alfap)/2
            rop=ro[nb]
            alfap=alfa[nb]
            i+=1
        nb=orderIndex[0]
        A+=ro[nb]*rop*np.sin(alfa[nb]-alfap)/2     
'''       
        
   

Zmax=np.zeros(Nn)
fx=np.zeros(Nn)
fy=np.zeros(Nn)
fz=np.zeros(Nn)
tipoNodo=np.zeros(Nn,dtype=np.uint8)
Anodo=np.zeros(Nn)
#conta=np.zeros(Nn,dtype=np.uint16)
# *****************************************************************************
#DETERMINAZIONE DEI CARICHI*********************************
# *****************************************************************************

gamma=16.00
Zr=2.700 #quota massima del materiale di riempimento
gammar=15.00 #peso specifico del materiale di riempimento fino alla quota Zr in kN/mc
q=20.00 # carico uniforme in kN/mq
nodo=0
while (nodo<Nn):
    if len(branchesINnode[nodo])>1:
        tipoNodo[nodo]=1
    nodo+=1


nodo=0
Stot=0.0
while (nodo<Nn):
    print("nodo: "+str(nodo+1))
    #START 
    A=0.0
    nbn=len(branchesINnode[nodo])       
    if nbn>1:
        dx=np.zeros(nbn)
        dy=np.zeros(nbn)
        NodesIndex=np.zeros(nbn+1,dtype=np.uint32)
        brnchs=branchesINnode[nodo]   #branches nel nodo=nodo
        x1=xyz[nodo][0]
        y1=xyz[nodo][1]         
        pn=[x1,y1]
  
        #ORDINA I NODI COLLEGATI DA BRANCHES AL NODO "nodo", IN VERSO ORARIO NON CONSIDERANDO I NODI ESTERNI
        ib=0
        while (ib<nbn):
            nI=branches[brnchs[ib],0] #nodo iniziale del branche corrente
            nJ=branches[brnchs[ib],1] #nodo finale del branche corrente
            pib=[xyz[nI][0],xyz[nI][1]] #coordinate del nodo iniziale del branche corrente
            pjb=[xyz[nJ][0],xyz[nJ][1]] #coordinate del nodo finale del branche corrente
            if pib != pn: #se il nodo iniziale del branche corrente non è il nodo stesso
                pb=pib    #considera questo come secondo estremo
                ndsIndx=nI
            else:        #se il nodo iniziale del branche corrente è il nodo stesso
                pb=pjb   #considera questo come secondo estremo
                ndsIndx=nJ
            dx[ib]=pb[0]-pn[0]
            dy[ib]=pb[1]-pn[1]
            NodesIndex[ib]=ndsIndx #tiene traccia dell'indice del secondo estremo
            ib+=1              
        ro=np.sqrt(dx**2+dy**2)#??????????????????????????????????????????????????
        alfa=np.arctan2(dy,dx)#ordina in base alla tangente dell'angolo
        orderIndex=np.argsort(alfa)            
        NodesIndex[ib]=NodesIndex[0]
        orderIndex= np.append(orderIndex, orderIndex[0])
        NodesIndex=NodesIndex[orderIndex]
        #fine ordinamento  
        
        ib=0
        while (ib<nbn):
            n1=NodesIndex[ib]
            n1next=NodesIndex[ib+1]
            #print ("nodo: "+str(nodo+1)+"/ n1: "+str(n1+1)+"/ n1next: "+str(n1next+1))
            n1n=len(branchesINnode[n1])       
            if n1n>0:
                brnchs=branchesINnode[n1] 
                i1=0
                while (i1<n1n):
                    nI=branches[brnchs[i1],0] #nodo iniziale del branche corrente
                    nJ=branches[brnchs[i1],1] #nodo finale del branche corrente            
                    if nI != n1:
                        n2=nI
                    else:
                        n2=nJ
                    deltaX1=xyz[n1][0]-xyz[nodo][0]    
                    deltaY1=xyz[n1][1]-xyz[nodo][1]                       
                    deltaX2=xyz[n2][0]-xyz[nodo][0]    
                    deltaY2=xyz[n2][1]-xyz[nodo][1]
                    
                   # ba = np.array([deltaX1,deltaY1])
                    #bc = np.array([deltaX2,deltaY2])

                    #cosine_teta = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                    #teta = np.arccos(cosine_teta)
                    
                    a_phase = cmath.phase(complex(deltaX1,deltaY1))
                    b_phase = cmath.phase(complex(deltaX2,deltaY2))
                    teta= b_phase-a_phase
                    if n2==nodo:
                        teta=0.000
                    """
                    alfa1=np.arctan2(deltaY1,deltaX1)#rivedere
                    if (alfa1<0.00):
                        alfa1=2*np.pi+alfa1
                    alfa2=np.arctan2(deltaY2,deltaX2)
                    if (alfa2<0.00):
                        alfa2=2*np.pi+alfa2
                    
                    
                    teta=alfa2-alfa1
                    """
                    #print ("teta: "+str(nodo+1)+"/"+str(n1+1)+"/"+str(n2+1)+"/"+str(np.degrees(teta)))
                    n2OK=0
                    if ((teta>0.0)and(teta<np.pi))or((teta>-2*np.pi)and(teta<-np.pi)):
                        n2n=len(branchesINnode[n2])
                        if n2n>0:
                            brnchs2=branchesINnode[n2]
                            i2=0
                            while (i2<n2n):
                                nn=branches[brnchs2[i2],0] #nodo iniziale del branche corrente
                                nm=branches[brnchs2[i2],1] #nodo finale del branche corrente   
                                if nn != n2:
                                    n3=nn
                                else:
                                    n3=nm 
                                if (n3==nodo):
                                    print ("triangle: "+str(nodo+1)+"/"+str(n1+1)+"/"+str(n2+1)+"/"+str(n3+1)+"/"+str(np.degrees(teta)))
                                    n2OK=1
                                else:
                                    if (n3==n1next):
                                        print("quad: "+str(nodo+1)+"/"+str(n1+1)+"/"+str(n2+1)+"/"+str(n3+1)+"/"+str(nodo+1)+"/"+str(np.degrees(teta)))
                                        n2OK=1
                                i2+=1
                    if n2OK==1:
                        break
                    i1+=1
                ib+=1

    fz[nodo]=-5.00
    nodo+=1     
print("Ho finito: "+str(nodo))
      
"""
        ib=0
        while (ib<nbn):
            n1=NodesIndex[ib]  
            n2=NodesIndex[ib+1]
            if tipoNodo[n1]==0:
                ib+=1
                continue
            if tipoNodo[n2]==0:
                ib+=1
                continue        
            nbn1=len(branchesINnode[n1])
            nbn2=len(branchesINnode[n2])
            brnchsn1=branchesINnode[n1]
            brnchsn2=branchesINnode[n2]
            ib1=0            
            #connesso=0
            ne=[]#nodi esterni connessi sia a n1 che a n2              
            while (ib1<nbn1):
                if brnchsn1[ib1] in brnchsn2:
                   #connesso=1 #i due nodi n1 ed n2  sono connessi dal branch ib1
                   ne=[]
                   #connessocon=ib1                                    
                   break
                #trovare nodi comuni                 
                #trova i nodi esterni "ne" connessi da branch sia da n1 che da n2
                n1I=branches[brnchsn1[ib1],0]
                n1J=branches[brnchsn1[ib1],1]
                n1e=n1J
                if n1I != n1:
                    n1e=n1I   
                ib2=0
                while (ib2<nbn2):
                    #da completae con la memorizzazione di più nodi esterni trovati
                    n2I=branches[brnchsn2[ib2],0]
                    n2J=branches[brnchsn2[ib2],1] 
                    n2e=n2J
                    if n2I != n2:
                        n2e=n2I
                    if (n1e==n2e) and (not (n1e==nodo)):
                        ne+=[n1e]
                    ib2+=1                  
                ib1+=1#
            #n3=nodo  
            nNodilinked=len(ne)
            if nNodilinked>0:
                n=ne[0]
                n3=n
                if nNodilinked>1:#prendendo come esterno il nodo "ne" più vicino a "nodo"     
                    rnmin=np.sqrt((xyz[0][0]-x1)**2+(xyz[0][1]-y1)**2)            
                    jn=1
                    while (jn<nNodilinked):
                        n=ne[jn]
                        rn=np.sqrt((xyz[n][0]-x1)**2+(xyz[n][1]-y1)**2)       
                        if rn<rnmin:
                            n3=n
                            rnmin=rn                   
                        jn+=1 
                #print "nodo="+str(nodo+1)+", "+str(n1+1)+", "+ str(n2+1)+", "+str(np.array(ne)+1)+", scelto = "+str(n3+1) 
                i1=nodo
                i2=n1
                i3=n3
                A0=(xyz[i2][0]-xyz[i1][0])*(xyz[i3][1]-xyz[i1][1])
                A0+=-(xyz[i3][0]-xyz[i1][0])*(xyz[i2][1]-xyz[i1][1])
                A0=abs(A0)/2
                Anodo[nodo]+=A0/4
                i2=n3
                i3=n2
                A0=(xyz[i2][0]-xyz[i1][0])*(xyz[i3][1]-xyz[i1][1])
                A0+=-(xyz[i3][0]-xyz[i1][0])*(xyz[i2][1]-xyz[i1][1])
                A0=abs(A0)/2
                Anodo[nodo]+=A0/4
            else:
                i1=nodo
                i2=n1
                i3=n2
                A0=(xyz[i2][0]-xyz[i1][0])*(xyz[i3][1]-xyz[i1][1])
                A0+=-(xyz[i3][0]-xyz[i1][0])*(xyz[i2][1]-xyz[i1][1])
                A0=abs(A0)/2
                Anodo[nodo]+=A0/3   

            #print str(n1+1)+", "+ str(n2+1)+ ",An1="+str(An[0]) 
            #y2=(y1+xyz[n1][1])/2                      
            #x4=(x1+xyz[n2][0])/2
            #y4=(y1+xyz[n2][1])/2                        
            '''
            #ok="nook"   
            if connesso==0:      
                x3=(x1+xyz[n3][0])/2
                y3=(y1+xyz[n3][1])/2
            else:                             
                x3=(xyz[n1][0]+xyz[n2][0])/2
                y3=(xyz[n1][1]+xyz[n2][1])/2

            A1=(x3-x1)*(y2-y1)
            A1+=-(x2-x1)*(y3-y1)
            A1=abs(A1)/2
            A2=(x4-x1)*(y3-y1)
            A2+=-(x3-x1)*(y4-y1)
            A2=abs(A2)/2
            A+=A1+A2
            '''
            #A1=(x4-x1)*(y2-y1)
            #A1+=-(x2-x1)*(y4-y1)
            #A+=abs(A1)/2
                         
            #print 
            #ib+=1
        #da completare con calcolo baricentro e area della terna o quaderna di nodi collegati a "nodo"             

    #print str(nodo+1)+", " +str(NodesIndex+1)+", A="+str(A)
    #fz[nodo]=-Anodo[nodo]*((xyz[nodo][3]-xyz[nodo][2])*gamma+(Zr-xyz[nodo][3])*gammar+q)
    #print fz[n]
   # Stot+=Anodo[nodo]



print ("L'area Totale è pari a: "+str(Stot))
# *****************************************************************************
"""
# *****************************************************************************
# SCRITTURA FILE PER APPLICAZIONE ROSASATI.PY
#******************************************************************************
dat_file.write(str(Nn)+"\n")

n=0

while (n<Nn):
    x=float(xyz[n][0])
    y=float(xyz[n][1])
    zmin=float(xyz[n][2])
    zmax=float(xyz[n][3])
    x=format(x, '+.6f') 
    y=format(y, '+.6f')
    zmin=format(zmin, '+.6f')
    zmax=format(zmax, '+.6f')  
    fxn=format(float(fx[n]), '+.3f')
    fyn=format(float(fy[n]), '+.3f')
    fzn=format(float(fz[n]), '+.3f')
    dat_file.write(str(n+1)+"  "+ str(x)+"  "+str(y)+"  "+str(zmin)+"  "+str(zmax)+"  "+str(fxn)+"  "+str(fyn)+"  "+str(fzn)+"  "+str(tipoNodo[n])+"  "+str(NodesLayer[n])+"\n")
    n+=1

dat_file.write(str(Nb)+"\n")

b=0
while (b<Nb):
    dat_file.write(str(branches[b,0]+1)+" "+str(branches[b,1]+1)+" "+str(branchesdmin[b,0])+"\n")
    b+=1

dat_file.close()
print ("...Elaborazione conclusa")
print ("Dati salvati nel file " + file + ".DAT")
print ("Esportati " + str(Nn) + " Nodi e " + str(Nb) + " branches")
tempo_finale = time.time()
print ("Impiegati ", str(tempo_finale - tempo_iniziale), " secondi per l'esportazione.")
