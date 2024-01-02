# -*- coding: utf-8 -*-
"""
Created on Mon Jan  1 20:57:27 2024

@author: lucky
"""
import numpy as np

def ottimizza(r1,d1,tol,):
    rprecedente=r1
    err=10
    #INIZIO PRCEDURA ITERATIVA
    print ("")    
    print ("OTTIMIZZAZIONE rmin")
    niter=0
    d=np.array(d1)
    while (err>tol):
        niter+=1    
        Do=np.dot(Mc*(d/lh),Mc.T)  
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
 