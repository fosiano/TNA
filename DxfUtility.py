# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:41:26 2023

@author: lucky
"""

#******************************************************************************
#SEZIONE DXF INPUT - OUTPUT
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

def LineToDXF(dxf_file, P1, P2, layer, clrlnscl=(256, 1.00)): 
    dxf_file.write("0\n")
    dxf_file.write("LINE\n")
    dxf_file.write("8\n")
    dxf_file.write(layer+"\n")
    dxf_file.write("62\n")
    dxf_file.write(str(clrlnscl[0])+"\n")    
    dxf_file.write("48\n")
    dxf_file.write(str(clrlnscl[1])+"\n")        
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

def WriteIntestazioneDXF(dxf_file):     
    dxf_file.write("999\n")
    dxf_file.write("DXF created from myself (Fortunato Siano)\n")
    dxf_file.write("0\n")
    dxf_file.write("SECTION\n")
    dxf_file.write("2\n")
    dxf_file.write("ENTITIES\n")
    print(dxf_file.name+" aperto in Scrittura")
  
  
# =============================================================================
# def WriteNewLayerDXF(dxf_file, layerName, color):
#     dxf_file.write("0\n")
#     dxf_file.write("SECTION\n")
#     dxf_file.write("2\n")
#     dxf_file.write("TABLES\n")
#     dxf_file.write("0\n")
#     dxf_file.write("TABLE\n")
#     dxf_file.write("2\n")    
#     dxf_file.write("LAYER\n")
#     dxf_file.write("2\n")
#     dxf_file.write(layerName+"\n")
#     dxf_file.write("62\n")
#     dxf_file.write(color+"\n")
#     dxf_file.write("ENDTAB\n")
#     dxf_file.write("0\n")
#     dxf_file.write("ENDSEC\n")
#     dxf_file.write("0\n")
#     dxf_file.write("SECTION\n")
#     dxf_file.write("2\n")
#     dxf_file.write("ENTITIES\n")
#     print(dxf_file.name+" aperto in Scrittura")    
# =============================================================================
    
    
def CloseWrittenDXF(dxf_file): 
    dxf_file.write("0\n")
    dxf_file.write("ENDSEC\n")
    dxf_file.write("0\n")
    dxf_file.write("EOF\n")
    dxf_file.close()    
    print(dxf_file.name+" chiuso")
    
#******************************************************************************    
#******************************************************************************
