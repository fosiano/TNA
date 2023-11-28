import time
import numpy as np
import xalglib as xal

from scipy.optimize import minimize

#d=densità di forza (ex ths = "th segnato")

file="RosatiGroinVault" #file di lavoro INPUT (.DAT)/OUTPUT (.DXF)

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


tolopt=1.0e-10 #tolleranza nelle procedure di ottimizzazione delle z (r) 
itermax=250
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
Nb=int(line) #numero dei branches
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
  print ("passo:"+str(Nfeval)+" / "+str(objectiverp(t)))
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
ct=zero.tolist() # 0 sta ad indicare che il vinccolo consite in "equality" 


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
xal.minbleicsetlc(state, c, ct,2*Ni)
xal.minbleicsetcond(state, epsg, epsf, epsx, maxits)


#



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
        Do=np.dot(Mc*(d/lh),Mc.T)  
        D=np.append(Do, fz, axis=1)
        D=np.delete(D, exttpl, axis=0)
        
        #Shallowest configuration network (maximum thrust) ///////, jac=DERobjectiverp
        Nfeval=1
        solutionZrmin = minimize(objectiverp,zmedr,method='SLSQP',bounds=bnds, constraints=equZ, jac=DERobjectiverp, callback=callbackFZrp, options={'disp': True, 'maxiter': itermax,'ftol': tolopt})
        zrmin = solutionZrmin.x #Shallowest configuration network (maximum thrust)
        rmin=zrmin[Nn]
        OKrmin=solutionZrmin.success
        print("OKrmin = ",OKrmin)
        if not(OKrmin): break 
        if VerticalLoads: break 
        if rmin>r1:
            d=d0+rmin/r1*(d1-d0) 
            err=(rmin-rprecedente)/rprecedente    
            rprecedente=rmin
        else:
            attenz="r1 TR0PPO GRANDE"
            break
    
    
    zmin=np.delete(zrmin, Nn, axis=0)
    lzmin=np.dot(Mc.T,zmin)
    tmax=(1/rmin)*d*np.sqrt(lh**2+lzmin**2)/lh
    
    
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
        solutionZrmax = minimize(objectivern,zmedr,method='SLSQP',bounds=bnds, constraints=equZ, jac=DERobjectivern, callback=callbackFZrn, options={'disp': True, 'maxiter': itermax,  'ftol': tolopt})
        zrmax = solutionZrmax.x #deepest configuration network (minimum thrust)
        rmax=zrmax[Nn]
        OKrmax=solutionZrmax.success
        if not(OKrmax): break
        if VerticalLoads: break
        err=(rmax-rprecedente)/rprecedente
        d=d0+rmax/r1*(d1-d0)
        rprecedente=rmax
     
    zmax=np.delete(zrmax, Nn, axis=0)
    lzmax=np.dot(Mc.T,zmax)
    tmin=(1/rmax)*d*np.sqrt(lh**2+lzmax**2)/lh
    
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
        solutionZrmed = minimize(objectivermed,zmedr,method='SLSQP',bounds=bnds, constraints=equZ, callback=callbackFZmed, options={'disp': True, 'maxiter': itermax,'ftol': tolopt})
        zrmed = solutionZrmed.x #
        rmed=zrmed[Nn] 
        OKrmed=solutionZrmed.success
        if not(OKrmed): break
        if VerticalLoads: break
        err=(rmed-rprecedente)/rprecedente
        d=d0+rmed/r1*(d1-d0)
        rprecedente=rmed  
    
    zmed=np.delete(zrmed, Nn, axis=0)
    lzmed=np.dot(Mc.T,zmed)
    tmed=(1/rmed)*d*np.sqrt(lh**2+lzmed**2)/lh
    print ("")
    print ("file di INPUT = " + file + ".dat")
    print (attenz)
    print ("r+ = rmin = " + str(rmin) + " in " + str(niterMIN)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmin))
    print ("rmed = " + str(rmed) + " in " + str(niterMED)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmed))
    print ("r- = rmax = " + str(rmax) + " in "+ str(niterMAX)+" iterazioni, con convergenza nei singoli passi: "+ str(OKrmax))
    print ("file di OUTPUT = " + file + ".dxf")       
    tempo_finale = time.time()
    print ("Impiegati ", str(tempo_finale - tempo_iniziale), " secondi per la soluzione.")
    
    #************solutionZ = minimize(objectiver,zmedr,method='SLSQP',constraints=consz,options={'disp': True, 'eps': 1.0e-05, 'maxiter': 500, 'ftol': 1e-05})
    
    #______________________________________________________________________________
    #
    #END PROCEDURA ROSATI
    #______________________________________________________________________________
    
    
    tbmax=np.max(tmax)
    tbmin=np.min(tmin)
    deltat=(tbmax-tbmin)/9
    
    
    
    
        
    #******************************************************************************
    #SALVA IN DXF INPUT - OUTPUT
    #******************************************************************************
    dxf_file = open(file+"-OUT.dxf","w")
    force_file = open(file+"-OUT.frc","w")
    dxf_file.write("999\n")
    dxf_file.write("DXF created from myself (Fortunato Siano)\n")
    dxf_file.write("0\n")
    dxf_file.write("SECTION\n")
    dxf_file.write("2\n")
    dxf_file.write("ENTITIES\n")
    
    b=0 #itera sui branches
    while (b<Nb):
            n=np.argmax(Mc[:,[b]]) #nodoI
            m=np.argmin(Mc[:,[b]]) #nodoJ
            #'''
            color=21+20*int(np.around((tmax[b]-tbmin)/deltat,decimals=0))
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("MAXSPINTA\n")
            dxf_file.write("62\n")       
            dxf_file.write(str(color)+"\n") 
            dxf_file.write("48\n")       
            dxf_file.write(str(tmax[b])+"\n")                    
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zrmin[n])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zrmin[m])+"\n")

            color=21+20*int(np.around((tmin[b]-tbmin)/deltat,decimals=0))     
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("MINSPINTA\n")
            dxf_file.write("62\n")
            dxf_file.write(str(color)+"\n")
            dxf_file.write("48\n")       
            dxf_file.write(str(tmin[b])+"\n")  
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zrmax[n])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zrmax[m])+"\n")

            color=21+20*int(np.around((tmed[b]-tbmin)/deltat,decimals=0))    
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("MEDIONETWORK\n")
            dxf_file.write("62\n")
            dxf_file.write(str(color)+"\n")
            dxf_file.write("48\n")       
            dxf_file.write(str(tmed[b])+"\n")  
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zrmed[n])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zrmed[m])+"\n")
    
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("ASSE\n")
            dxf_file.write("62\n")
            dxf_file.write("253\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zmedr[n])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zmedr[m])+"\n")
            #'''
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("INTRADOSSO\n")
            dxf_file.write("62\n")
            dxf_file.write("256\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBounds[n,0])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zBounds[m,0])+"\n")    
            
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("ESTRADOSSO\n")
            dxf_file.write("62\n")
            dxf_file.write("256\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBounds[n,1])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zBounds[m,1])+"\n")     
            '''
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("lim-inf\n")
            dxf_file.write("62\n")
            dxf_file.write("13\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBoundsRID[n,0])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zBoundsRID[m,0])+"\n")    
            
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("lim-sup\n")
            dxf_file.write("62\n")
            dxf_file.write("14\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBoundsRID[n,1])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[m])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[m])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zBoundsRID[m,1])+"\n")           
            '''
            force_file.write(str(b)+", "+str(n)+", "+str(m)+", "+str(tmin[b])+", "+str(tmed[b])+", "+str(tmax[b])+"\n")              
            
            
            b+=1
    
    n=0 #itera sui nodes
    while (n<Nn):      
            dxf_file.write("0\n")
            dxf_file.write("POINT\n")
            dxf_file.write("8\n")
            #if tipo[n]==1:
            #    dxf_file.write("NODES-I\n")
            #else:
            #    dxf_file.write("NODES-E\n")        
            dxf_file.write(NodesLayer[n]+"\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBounds[n,0])+"\n")
      
            dxf_file.write("0\n")
            dxf_file.write("TEXT\n")
            dxf_file.write("8\n")
            dxf_file.write("ID-NODES\n")     
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBounds[n,0])+"\n")
            dxf_file.write("40\n")
            dxf_file.write("0.15\n")
            dxf_file.write("1\n")
            dxf_file.write(str(n+1)+"\n")
            
        
            dxf_file.write("0\n")
            dxf_file.write("LINE\n")
            dxf_file.write("8\n")
            dxf_file.write("fz\n")
            dxf_file.write("62\n")
            dxf_file.write("7\n")
            dxf_file.write("10\n")  
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("20\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("30\n")
            dxf_file.write(str(zBounds[n,1])+"\n")
            dxf_file.write("11\n")
            dxf_file.write(str(x[n])+"\n")
            dxf_file.write("21\n")
            dxf_file.write(str(y[n])+"\n")
            dxf_file.write("31\n")
            dxf_file.write(str(zBounds[n,1]-2.0*fz[n,0])+"\n")    
            #dxf_file.write(str(zBounds[n,1]-2.0*fz[n])+"\n")   
                  
            n+=1
            
    dxf_file.write("0\n")
    dxf_file.write("ENDSEC\n")
    dxf_file.write("0\n")
    dxf_file.write("EOF\n")
    
    dxf_file.close()
    force_file.close()
else:
    print ("L'ottimizzazione di d non converge")
#******************************************************************************


#PER IL CALCOLO DELL'AREA INCIDENTE SUL NODO
#x=np.array([5.0,6.0,7.0])
#y=np.array([8.0,9.0,10.0])
#alfa=np.arctan2(y,x)
#orderIndex=np.argsort(alfa)

#???
#q=np.array([3,1,2])
#q=np.matrix(q).T
#alfa1=np.append(alfa, q, axis=1)
#alfa1[np.argsort(alfa1.A[:, 0])]
# alfa1[np.lexsort(np.fliplr(alfa1).T)]
#non SERVE if alfa<0:alfa=(2*np.pi + alfa)


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
