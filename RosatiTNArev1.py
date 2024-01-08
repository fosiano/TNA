import time
import numpy as np
import xalglib as xal
#from scipy.optimize import show_options
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
from DxfUtility import *

class Point:
  def __init__(self, x, y, z):
      self.x = x
      self.y = y
      self.z = z

Radicedi2=(2.0**0.50)

def CalcolaStampa(th):
    RyLatoDown=np.sum(th[ExternalBranchesIndexesDown])    
    RyLatoUp=np.sum(th[ExternalBranchesIndexesUp])
    RxLatoLeft=np.sum(th[ExternalBranchesIndexesLeft])
    RxLatoRight=np.sum(th[ExternalBranchesIndexesRight])
    
       
    RDiagDownLeft=th[ExternalBranchesIndexesDiag[0]]/Radicedi2       
    RDiagDownRight=th[ExternalBranchesIndexesDiag[1]]/Radicedi2
    RDiagUpLeft=th[ExternalBranchesIndexesDiag[2]]/Radicedi2  
    RDiagUpRight=th[ExternalBranchesIndexesDiag[3]]/Radicedi2  
    if len(ExternalBranchesIndexesDiag) == 8:
        RDiagDownLeft+=th[ExternalBranchesIndexesDiag[4]]/Radicedi2       
        RDiagDownRight+=th[ExternalBranchesIndexesDiag[5]]/Radicedi2
        RDiagUpLeft+=th[ExternalBranchesIndexesDiag[6]]/Radicedi2  
        RDiagUpRight+=th[ExternalBranchesIndexesDiag[7]]/Radicedi2   
    
    
    RyDown = RyLatoDown + RDiagDownLeft + RDiagDownRight
    RyUp = RyLatoUp + RDiagUpLeft + RDiagUpRight
    RxLeft = RxLatoLeft + RDiagDownLeft + RDiagUpLeft
    RxRight = RxLatoRight + RDiagDownRight + RDiagUpRight
    
    Ry=RyDown-RyUp
    Rx=RxLeft-RxRight
    print("RyLatoDown = ", RyLatoDown, "RDiagDownLeft = ", RDiagDownLeft ,"RDiagDownRight = ",RDiagDownRight )
    print("RyLatoUp = ", RyLatoUp, "RDiagUpLeft = ", RDiagUpLeft,"RDiagUpRight = ", RDiagUpRight)
    
    
    print("RxLatoLeft = ", RxLatoLeft, "RDiagDownLeft = ", RDiagDownLeft,"RDiagUpLeft = ",  RDiagUpLeft)
    print("RxLatoright = ", RxLatoRight, "RDiagDownRight = ", RDiagDownRight,"RDiagUpRight = ", RDiagUpRight )
    print("RyDown = ", RyDown)
    print("RyUp = ", RyUp)
    print("RxLeft = ", RxLeft)
    print("RxRight = ", RxRight)
    print("Rx = ", Rx)
    print("Ry = ", Ry)


#show_options(solver= 'minimize', method='SLSQP', disp=True)
#d=densità di forza (ex ths = "th segnato")


file="VelaRect" #file di lavoro INPUT (.DAT)/OUTPUT (.DXF)

#file="RosatiGroinVault" #file di lavoro INPUT (.DAT)/OUTPUT (.DXF)

#file="ArcoRosatiFH" #file di lavoro INPUT (.DAT)/OUTPUT (.DXF)

tempo_iniziale = time.time()

use=0.0 # distanza minima della curva delle pressioni dall'estradosso come frazione dello spessore della volta 
usi=use # distanza minima della curva delle pressioni dall'intradosso come frazione dello spessore della volta 

#minimod=0.10
r1=0.01 #
tol=0.001 #tolleranza accettata per interruzione dell'iterazione nel caso di carichi orizzontali presenti

# These variables define stopping conditions for the optimizer xalglib.
# We use very simple condition - |g|<=epsg
epsg = 1.0e-8
epsf = 0.
epsx = 0.
maxits = 0


tolopt=1.0e-7 #tolleranza nelle procedure di ottimizzazione delle z (r) 
itermax=200
#******************************************************************************
#LETTURA FILE DATI INPUT - MEMORIZZAZIONE DATI PER ELABORAZIONE
#******************************************************************************
network = open(file+".dat","r")

print("file di LAVORO = " + file + ".dat")
time.sleep(0.1)

line=network.readline()

Nn=int(line) #numero dei nodi

#NB: gli indici dei vettori/matrici partono da 0!!!!
x=np.zeros(Nn)
y=np.zeros(Nn)
zBounds=np.zeros(shape=(Nn,2))
fx=np.zeros(Nn)
fy=np.zeros(Nn)
fz=np.zeros(Nn)
tipo=np.zeros(Nn, dtype=np.uint8)
NodesLayer=np.empty(Nn,dtype=np.dtype('U16'))

exttpl=()

#internal branches: represent a thrust force that is interior to the network
#external branches: represent the support reaction
#edge branches: represent forces that are on a free edge

Ni=0 #numero dei nodi interni (connected to internal branches)
Nr=0 #numero dei nodi esterni [(restrained)] (connected only to external branches)
n=0 #itera sui nodi
while (n<Nn):
      line=network.readline()
      dati = line.split()
      x[n]=float(dati[1]) #x
      y[n]=float(dati[2]) #y    
      zBounds[n,0]=float(dati[3]) #zmin     
      zBounds[n,1]=float(dati[4]) #zmax      
      fx[n]=float(dati[5]) #fx
      fy[n]=float(dati[6]) #fy   
      fz[n]=float(dati[7]) #fz
      flag=int(dati[8]) #flag external/internal/edge
      tipo[n]=flag  
      if (flag==1):      
          Ni=Ni+1
      if (flag==0):      
          Nr=Nr+1
          exttpl+=(n,)
      NodesLayer[n]=dati[9] #Layer del Nodo
      n+=1

Ne=Nn-Nr #numero dei nodi esterni liberi (connected to edge branches)

#branchesINnode = [[] for i in range(Nn)]

line=network.readline()
dati = line.split()
Nb=int(dati[0]) #numero dei branches
Nbe=int(dati[1]) #numero dei branches esterno
Nbi=Nb-Nbe #numero dei branches interni
#branches=np.zeros(shape=(Nb,2),dtype=np.uint16)
#branchesflag=np.zeros(Nb, dtype=np.uint8)
Mc=np.zeros(shape=(Nn,Nb))

dmin=np.zeros(Nb)

b=0 #itera sui branches
while (b<Nb):
      line=network.readline()
      dati = line.split()
      nI=int(dati[0])-1 #nodoI
      nJ=int(dati[1])-1#nodoI
      dmin[b]=float(dati[2])
      Mc[nI,b]=1.0
      Mc[nJ,b]=-1.0   
      #branches[b,0]= nI
      #branches[b,1]= nJ
      #branchesINnode[nI]+=[branches[b,0]]
      #branchesINnode[nJ]+=[branches[b,0]]
      b=b+1
      
dx=np.dot(Mc.T,x)
dy=np.dot(Mc.T,y)
lh=np.sqrt(dx**2+dy**2)

def leggiArray(file):
    line=file.readline()
    while not ("]" in line):
        lineAdd=file.readline()
        line += lineAdd
    ar=np.fromstring( line[1:-2], dtype=np.int32, sep=' ' )    
    return ar


ExternalBranchesIndexesDown = leggiArray(network) 
ExternalBranchesIndexesUp = leggiArray(network) 
ExternalBranchesIndexesLeft = leggiArray(network) 
ExternalBranchesIndexesRight = leggiArray(network) 
ExternalBranchesIndexesDiag = leggiArray(network) 

ArcNodeIndexes = leggiArray(network) 
ExternalArcNodeIndexes = leggiArray(network) 


line=network.readline()
PesoTotale=float(line)
#qqq=np.asarray(branchesINnode)
#print len(qqq[2])
#print len(branchesINnode[0])
#print b
network.close()
#******************************************************************************

ib = np.transpose(np.ones(Nb))
unit=(np.ones(Nb)).tolist()
def function1_grad(t, grad, param):
    global ib
    global unit
    fn=np.dot(ib, np.array(t))
    grad=unit
    return fn

def function1(t):
    global ib
    fo=np.dot(ib, np.array(t))
    return fo


def objectiverp(zr):
    a=zr[Nn]
    return a

def objectivern(zr):
    a=zr[Nn]
    return -a
      
def equilibrioZ(zr):
    #return np.squeeze(np.asarray(D*np.matrix(z).T))
    #return np.array(D*np.matrix(z).T).ravel()
    return np.ravel(D*np.matrix(zr).T)#COMMENTARE RAVEL????????????????????????????????????

def equalZ(zr):
   #a=np.zeros(len(zr))
   # a=zr.take(ArrayArcIndexes)
   # b=zr.take(ArrayFalseArcIndexes)
   
   a=zr.take(ArcNodeIndexes)
   b=zr.take(ExternalArcNodeIndexes  )
   c=a-b
   return c
        


jacobrp=np.zeros(Nn+1)

jacobrp[Nn]=1.00
def DERobjectiverp(zr):
    return jacobrp

jacobrn=np.zeros(Nn+1)
jacobrn[Nn]=-1.00
def DERobjectivern(zr):
    return jacobrn

def callbackFZrp(t):
  global Nfeval
  isZero=equilibrioZ(t)
  SommaQuadrati=np.sqrt(np.sum(isZero*isZero))
  print ("passo:"+str(Nfeval)+" / r = "+str(objectiverp(t))+" / sQFeq = "+str(SommaQuadrati))
  #print(d)
  Nfeval += 1

def callbackFZrn(t):
  global Nfeval
  print ("passo:"+str(Nfeval)+" / "+str(objectivern(t)))
  #print(d)
  Nfeval += 1

#******************************************************************************
ib2 = np.transpose(np.ones(Ni))
def objectivermed(zr):
    a=np.delete(zr, Nn, axis=0)
    a = np.delete(a, exttpl, axis=0)
    b=(a-zmed)**2
    cy=np.dot(ib2, b)/Ni
    return cy

'''
da studiare
jacobrmed=np.zeros(Nn+1)
jacobrmed[Nn]=-1.00
def DERobjectivermed(zr):
    a=np.delete(zr, Nn, axis=0)
    a = np.delete(a, exttpl, axis=0)
    b=2*(a-zmed)
    cy=np.dot(ib2, b)/Ni       
    return jacobrmed 
'''
 
def callbackFZmed(t):
  global Nfeval
  print ("passo:"+str(Nfeval)+" / "+str(objectivermed(t)))
  #print(d)
  Nfeval += 1

#______________________________________________________________________________
#
#START PROCEDURA ROSATI
#______________________________________________________________________________


print ("OTTIMIZZAZIONE d PER r=r0 ...")
time.sleep(0.1)
#******************************************************************************
#EVALUATION OF THRUSTS' HORIZONTAL COMPONENTS
#******************************************************************************
C=Mc*(dx/lh)
S=Mc*(dy/lh)

#print C
C = np.delete(C, exttpl, axis=0)
S = np.delete(S, exttpl, axis=0)
CS=np.append(C, S, axis=0)

fx = np.delete(fx, exttpl, axis=0)
fy = np.delete(fy, exttpl, axis=0)

fh=np.append(fx, fy, axis=0)
fh=np.matrix(fh).T

zero=np.zeros(2*Ni)
ct=zero.tolist() # 0 sta ad indicare che il vincolo consite in "equality" 


lowerbnd=dmin.tolist()
upperbnd=np.ones(Nb)*np.inf
upperbnd=upperbnd.tolist()

CS0=np.append(CS, np.matrix(zero).T, axis=1)
c=CS0.tolist()

#
# Now we are ready to actually optimize something:
# * first we create optimizer
# * we add linear constraints
# * we tune stopping conditions
# * and, finally, optimize and obtain results...
#
state = xal.minbleiccreate(lowerbnd)
xal.minbleicsetbc(state, lowerbnd, upperbnd)
xal.minbleicsetlc(state, c, ct, 2*Ni)
xal.minbleicsetcond(state, epsg, epsf, epsx, maxits)

xal.minbleicoptimize_g(state, function1_grad)
d0, rep = xal.minbleicresults(state)
# ...and evaluate these results
print ("... OTTIMIZZAZIONE d PER r=r0 COMPLETATA")
print ("f(x) = "+ str(function1(d0))+" - termination type code = "+str((rep.terminationtype)))
print ("number of iterations = "+str((rep.iterationscount))+" - number of gradient evaluations = "+str((rep.nfev)))

#print(rep.terminationtype) # expected 4

convergenzad=False
if rep.terminationtype==4:
    convergenzad=True
    VerticalLoads= np.array_equal(fh,np.matrix(zero).T)
    if not(VerticalLoads):
        print ("OTTIMIZZAZIONE d PER r=r1 ...")
        time.sleep(0.1)
        CS0=np.append(CS, -fh*r1, axis=1)
        c=CS0.tolist()
        state = xal.minbleiccreate(lowerbnd)
        xal.minbleicsetbc(state, lowerbnd, upperbnd)
        xal.minbleicsetlc(state, c, ct,2*Ni)
        xal.minbleicsetcond(state, epsg, epsf, epsx, maxits)
        xal.minbleicoptimize_g(state, function1_grad)
        d1, rep = xal.minbleicresults(state)
        print ("")        
        print ("... OTTIMIZZAZIONE d PER r=r1 COMPLETATA")
        print ("f(x) = "+ str(function1(d1))+" - termination type code = "+str((rep.terminationtype)))
        print ("number of iterations = "+str((rep.iterationscount))+" - number of gradient evaluations = "+str((rep.nfev)))
        d0=np.array(d0)
        convergenzad=convergenzad and rep.terminationtype==4
    else:
        d1=d0


#******************************************************************************
#SALVA IN DXF INPUT
#******************************************************************************
# =============================================================================
dxf_fileName=file+"_TNA-IN.dxf"
dxf_file = open(dxf_fileName,"w")    
WriteIntestazioneDXF(dxf_file)
#WriteNewLayerDXF(dxf_file, "prova", "6")
ScalaForze=0.30
b=0 #itera sui branches
while (b<Nb):
        n=np.argmax(Mc[:,[b]]) #nodoI
        m=np.argmin(Mc[:,[b]]) #nodoJ
        Pi=Point(x[n], y[n], zBounds[n,0])
        Pj=Point(x[m], y[m], zBounds[m,0])
        LineToDXF(dxf_file, Pi, Pj, "INTRADOSSO")        
        Pi=Point(x[n], y[n], zBounds[n,1])
        Pj=Point(x[m], y[m], zBounds[m,1])      
        LineToDXF(dxf_file, Pi, Pj, "ESTRADOSSO")       
        b+=1 
n=0 #itera sui nodes
while (n<Nn):     
    P=Point(x[n], y[n], zBounds[n,0])
    testo=str(n+1)
    PointToDXF(dxf_file, P, NodesLayer[n])
    TextToDXF(dxf_file, P,testo,0.05,"ID-NODES")    
    intensity=zBounds[n,1]-fz[n]*ScalaForze
    Pi=Point(x[n], y[n], zBounds[n,1])
    Pj=Point(x[n], y[n], intensity)      
    LineToDXF(dxf_file, Pi, Pj, "fz", (6, 1.00))
    n+=1    
CloseWrittenDXF(dxf_file)    
     
#******************************************************************************
#EVALUATION OF THE NODAL HEIGHTS
#******************************************************************************
if convergenzad:
    zmed=0.50*np.sum(zBounds,axis=1)
    zmedr=np.append(zmed, [1.0], axis=0)
    zmed= np.delete(zmed, exttpl, axis=0)
    '''
    dzBounds=np.diff(zBounds,axis=1)
    zBoundsRID=zBounds+np.append(usi*dzBounds,-use*dzBounds, axis=1)
    zBoundsRID[exttpl,0]=zBounds[exttpl,0]
    zBoundsRID[exttpl,1]=zBounds[exttpl,1]
    '''
    zBoundsRID=zBounds
    
    # bounds
    bnds=tuple(map(tuple, zBoundsRID))# spiegare map/tuple?????????
    bnds+= ((0.000000, None),)
        
    # equazione di equilibrio
    equZ = {'type': 'eq', 'fun': equilibrioZ}
    asseSimm = {'type': 'eq', 'fun': equalZ}
# =============================================================================
#     Bounds=np.vstack((zBounds, [0.000000, np.inf]))
#     zr_zero=np.zeros(len(bnds))
# =============================================================================
    
    
    fz=np.matrix(fz).T
    if not(VerticalLoads): 
        attenz="OK: Soluzione valida per carichi anche orizzontali poichè r1<rmin"
    else:
        attenz="OK: Soluzione valida per soli carichi verticali"       
    
    
    rprecedente=r1
    err=10
    #INIZIO PRCEDURA ITERATIVA
    print ("")    
    print ("OTTIMIZZAZIONE rmin")
    niterMIN=0
    d=np.array(d1)
    while (err>tol):
        niterMIN+=1    
        Do=np.dot(Mc*(d/lh), Mc.T)  
        D=np.append(Do, fz, axis=1)
        D=np.delete(D, exttpl, axis=0)
        
        
        
        #Shallowest configuration network (maximum thrust) ///////, jac=DERobjectiverp
        Nfeval=1
# =============================================================================
#         
#         
#         linear_constraint = LinearConstraint(D, zr_zero, zr_zero)
#         solutionZrmin = minimize(objectiverp, zmedr, method='trust-constr',bounds=Bounds, constraints=[linear_constraint], 
#                                  jac=DERobjectiverp,  
#                                  options={'disp': True, 'maxiter': itermax})
#         
#         
# =============================================================================
        solutionZrmin = minimize(objectiverp, zmedr, method='SLSQP',
                                 bounds=bnds, constraints=[equZ, asseSimm], 
                                 jac=DERobjectiverp, callback=callbackFZrp, 
                                 options={'disp': True, 'maxiter': itermax, 'ftol': tolopt})
  
        zrmin = solutionZrmin.x #Shallowest configuration network (maximum thrust)
        
        isZero=equilibrioZ(zrmin)
        SommaQuadrati=np.sqrt(np.sum(isZero*isZero))
        print("\nVERIFICA DELL'EQUILIBRIO LUNGO LA VERTICALE\n")
        
        #print("il vettore completo delle somma delle Fz nodali è pari a: \n\n", isZero, "\n")
        print("La Radice Quadrata della Somma dei Quadratri è pari a: ", SommaQuadrati, "\n")
        
        rmin=zrmin[Nn]
        OKrmin=solutionZrmin.success
        print("OKrmin = ",OKrmin)
        if not(OKrmin): break 
        if VerticalLoads: break 
        if rmin>r1:
            d=d0+rmin/r1*(d1-d0) 
            err=abs((rmin-rprecedente)/rprecedente)
            
            rprecedente=rmin
        else:
            attenz="r1 TR0PPO GRANDE"
            break
    
    zmin=np.delete(zrmin, Nn, axis=0)
    lzmin=np.dot(Mc.T,zmin)
    thmax=(1/rmin)*d
    tmax=thmax*np.sqrt(lh**2+lzmin**2)/lh
    tzmax=thmax*lzmin/lh
    Rz=np.sum(tzmax[Nbi:])
    print("Rz = ", Rz, "/",PesoTotale, " = ", Rz/PesoTotale)
    
    
    CalcolaStampa(thmax)
    
    
    
    rprecedente=r1
    err=10
    #INIZIO PROCEDURA ITERATIVA
    print ("")
    print ("OTTIMIZZAZIONE rmax")
    niterMAX=0
    d=np.array(d1)
    while (err>tol):  
        niterMAX+=1 
        Do=np.dot(Mc*(d/lh),Mc.T)  
        D=np.append(Do, fz, axis=1)
        D=np.delete(D, exttpl, axis=0)
        
        #deepest configuration network (minimum thrust) ///////, jac=DERobjectivern
        Nfeval=1
        solutionZrmax = minimize(objectivern,zmedr,method='SLSQP',bounds=bnds, constraints=[equZ, asseSimm], jac=DERobjectivern, callback=callbackFZrn, options={'disp': True, 'maxiter': itermax,  'ftol': tolopt})
        zrmax = solutionZrmax.x #deepest configuration network (minimum thrust)
        isZero=equilibrioZ(zrmax)
        SommaQuadrati=np.sqrt(np.sum(isZero*isZero))
        print("\nVERIFICA DELL'EQUILIBRIO LUNGO LA VERTICALE\n")
        
        #print("il vettore completo delle somma delle Fz nodali è pari a: \n\n", isZero, "\n")
        print("La Radice Quadrata della Somma dei Quadratri è pari a: ", SommaQuadrati, "\n")
        
        rmax=zrmax[Nn]
        OKrmax=solutionZrmax.success
        if not(OKrmax): break
        if VerticalLoads: break
        err=abs((rmax-rprecedente)/rprecedente)
        d=d0+rmax/r1*(d1-d0)
        rprecedente=rmax
     
    zmax=np.delete(zrmax, Nn, axis=0)
    lzmax=np.dot(Mc.T,zmax)
    #tmin=(1/rmax)*d*np.sqrt(lh**2+lzmax**2)/lh
    
    thmin=(1/rmax)*d
    tmin=thmin*np.sqrt(lh**2+lzmax**2)/lh
    tzmin=thmin*lzmax/lh
    Rz=np.sum(tzmin[Nbi:])
    print("Rz = ", Rz, "/",PesoTotale, " = ",Rz/PesoTotale)
    
    CalcolaStampa(thmin)
    
    
    rprecedente=r1
    err=10
    #INIZIO PRCEDURA ITERATIVA
    print ("")
    print ("OTTIMIZZAZIONE rmed")
    niterMED=0
    d=np.array(d1)
    while (err>tol):    
        niterMED+=1 
        Do=np.dot(Mc*(d/lh),Mc.T)  
        D=np.append(Do, fz, axis=1)
        D=np.delete(D, exttpl, axis=0)
           
        #CURVA CHE APPROSSIMA MEGLIO LA SUPERFICIE MEDIA DELLA VOLATA
        Nfeval=1
        solutionZrmed = minimize(objectivermed,zmedr,method='SLSQP',bounds=bnds, constraints=[equZ, asseSimm], callback=callbackFZmed, options={'disp': True, 'maxiter': itermax,'ftol': tolopt})
        zrmed = solutionZrmed.x #
        
        isZero=equilibrioZ(zrmed)
        SommaQuadrati=np.sqrt(np.sum(isZero*isZero))
        print("\nVERIFICA DELL'EQUILIBRIO LUNGO LA VERTICALE\n")
        
        #print("il vettore completo delle somma delle Fz nodali è pari a: \n\n", isZero, "\n")
        print("La Radice Quadrata della Somma dei Quadrati è pari a: ", SommaQuadrati, "\n")
        
        
        rmed=zrmed[Nn] 
        OKrmed=solutionZrmed.success
        if not(OKrmed): break
        if VerticalLoads: break
        err=abs((rmed-rprecedente)/rprecedente)
        d=d0+rmed/r1*(d1-d0)
        rprecedente=rmed  
    
    zmed=np.delete(zrmed, Nn, axis=0)
    lzmed=np.dot(Mc.T,zmed)
    tmed=(1/rmed)*d*np.sqrt(lh**2+lzmed**2)/lh
    
    thmed=(1/rmed)*d
    tmed=thmed*np.sqrt(lh**2+lzmed**2)/lh
    tzmed=thmed*lzmed/lh
    Rz=np.sum(tzmed[Nbi:])
    print("Rz = ", Rz, "/",PesoTotale, " = ",Rz/PesoTotale)
    
    CalcolaStampa(thmed)
    
    print ("")
    print ("file di INPUT = " + file + ".dat")
    print (attenz)
    print ("r+ = rmin = " + str(rmin) + " in " + str(niterMIN)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmin))
    print ("rmed = " + str(rmed) + " in " + str(niterMED)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmed))
    print ("r- = rmax = " + str(rmax) + " in "+ str(niterMAX)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmax))
    print ("file di OUTPUT = " + file + ".dxf")       
    tempo_finale = time.time()
    print ("Impiegati ", str(tempo_finale - tempo_iniziale), " secondi per la soluzione.")
    print("Convergenza generale: "+ str(OKrmed and  OKrmax and OKrmin))
    #************solutionZ = minimize(objectiver,zmedr,method='SLSQP',constraints=consz,options={'disp': True, 'eps': 1.0e-05, 'maxiter': 500, 'ftol': 1e-05})
    
    #______________________________________________________________________________
    #
    #END PROCEDURA ROSATI
    #______________________________________________________________________________
    
    
    tbmax=np.max([tmax, tmin,tmed])
    tbmin=np.min([tmax, tmin,tmed])
    deltat=(tbmax-tbmin)/9
           
    #******************************************************************************
    #SALVA IN DXF INPUT - OUTPUT
    #******************************************************************************
# =============================================================================
#     dxf_file = open(file+"_TNA-OUT.dxf","w")
#     
#     dxf_file.write("999\n")
#     dxf_file.write("DXF created from myself (Fortunato Siano)\n")
#     dxf_file.write("0\n")
#     dxf_file.write("SECTION\n")
#     dxf_file.write("2\n")
#     dxf_file.write("ENTITIES\n")
# =============================================================================
    force_file = open(file+"_TNA-OUT.frc","w")
    dxf_fileName=file+"_TNA-OUT.dxf"
    dxf_file = open(dxf_fileName,"w")    
    WriteIntestazioneDXF(dxf_file)
    
    
    b=0 #itera sui branches
    while (b<Nb):
            n=np.argmax(Mc[:,[b]]) #nodoI
            m=np.argmin(Mc[:,[b]]) #nodoJ
            #'''
            color=21+20*int(np.around((tmax[b]-tbmin)/deltat,decimals=0))
            Pi=Point(x[n], y[n], zrmin[n])
            Pj=Point(x[m], y[m], zrmin[m])
            LineToDXF(dxf_file, Pi, Pj, "MAXSPINTA", (color, tmax[b])) 
            
            color=21+20*int(np.around((tmin[b]-tbmin)/deltat,decimals=0))    
            Pi=Point(x[n], y[n], zrmax[n])
            Pj=Point(x[m], y[m], zrmax[m])
            LineToDXF(dxf_file, Pi, Pj, "MINSPINTA", (color, tmin[b])) 
          
            color=21+20*int(np.around((tmed[b]-tbmin)/deltat,decimals=0))  
            Pi=Point(x[n], y[n], zrmed[n])
            Pj=Point(x[m], y[m], zrmed[m])
            LineToDXF(dxf_file, Pi, Pj, "MEDIONETWORK", (color, tmed[b])) 
                          
            Pi=Point(x[n], y[n], zmedr[n])
            Pj=Point(x[m], y[m], zmedr[m])
            LineToDXF(dxf_file, Pi, Pj, "ASSE", (253, 1.00)) 
          
            Pi=Point(x[n], y[n], zBounds[n,0])
            Pj=Point(x[m], y[m], zBounds[m,0])
            LineToDXF(dxf_file, Pi, Pj, "INTRADOSSO")                
            
            Pi=Point(x[n], y[n], zBounds[n,1])
            Pj=Point(x[m], y[m], zBounds[m,1])      
            LineToDXF(dxf_file, Pi, Pj, "ESTRADOSSO") 
            
            
# =============================================================================
#EVENTUALE
#             Pi=Point(x[n], y[n], zBoundsRID[n,0])
#             Pj=Point(x[m], y[m], zBoundsRID[m,0])      
#             LineToDXF(dxf_file, Pi, Pj, "lim-inf") 
# 
#             Pi=Point(x[n], y[n], zBoundsRID[n,1])
#             Pj=Point(x[m], y[m], zBoundsRID[m,1])      
#             LineToDXF(dxf_file, Pi, Pj, "lim-sup") 
# =============================================================================
                      
            force_file.write(str(b)+", "+str(n)+", "+str(m)+", "+str(tmin[b])+", "+str(tmed[b])+", "+str(tmax[b])+"\n")              
                        
            b+=1
    
    n=0 #itera sui nodes
    while (n<Nn):  
            P=Point(x[n], y[n], zBounds[n,0])
            testo=str(n+1)
            PointToDXF(dxf_file, P, NodesLayer[n])
            TextToDXF(dxf_file, P,testo,0.05,"ID-NODES") 
          
            intensity=zBounds[n,1]-fz[n,0]*ScalaForze
            Pi=Point(x[n], y[n], zBounds[n,1])
            Pj=Point(x[n], y[n], intensity)      
            LineToDXF(dxf_file, Pi, Pj, "fz", (6, 1.00))

            n+=1
               
    
    CloseWrittenDXF(dxf_file)  
    
    force_file.close()
else:
    print ("L'ottimizzazione di d non converge")
#******************************************************************************
import matplotlib
import matplotlib.pyplot as plt
import csv

cmap = plt.cm.jet_r
norm = matplotlib.colors.Normalize(vmin=tbmin, vmax=tbmax)

color = np.empty((0, 5))
for n, thrust in enumerate(tmin):
    mbrColor = cmap(norm(thrust)) + (thrust,)
    color = np.append(color, [mbrColor], axis=0)
filename='Branches-Color.csv'
with open(filename, 'w', newline='') as csvfile:
    csvwriter=csv.writer(csvfile)
    csvwriter.writerows(color)    


'''
ini=np.transpose(np.ones(2*Ni))
def callbackF(t):
  global Nfeval
  uu=np.sqrt(np.dot(ini, (equilibrio(t))**2))
  print "passo:"+str(Nfeval)+" - "+str(funzione(t))+" - "+str(uu)
  #print(d)
  Nfeval += 1


def callbackSimplex(xk, **kwargs):
  global Niter
  uu=np.sqrt(np.dot(c.T, xk))
  ww=np.sqrt(np.dot(ini, (np.dot(A, xk))**2))
  print "passo:"+str(Niter)+" - "+str(uu)+" - "+str(ww)+" - "+str(np.min(xk))

  #print(d)
  Niter += 1
'''
