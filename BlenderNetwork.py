import bpy
import csv
import numpy as np
from numpy import genfromtxt
import mathutils
import math
from mathutils import Matrix 

# minuto 8:50

def importData():
    lowerVertices = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\VerticesLowerBounds.csv', delimiter = ',')
    upperVertices = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\VerticesUpperBounds.csv', delimiter = ',')
    edges = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\Edges.csv', delimiter = ',')
    faces = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\Faces.csv', delimiter = ',')

    #color = genfromtxt('Member-Color.csv', delimiter = ';')
    return lowerVertices, upperVertices, edges, faces #, color

def main(chunkSize, elementScale, redraw):
    lowerVertices, upperVertices, edges, faces = importData() #, color
    #print(edges)
    #print(vertices)
    print(faces)
    coll=bpy.data.collections.new("lowerBounds")
    bpy.context.scene.collection.children.link(coll)

    name = 'network' 
    AllVertices=[]
    for vertice in lowerVertices:
        P=vertice   
        AllVertices += [tuple(P)] 
    AllUpperVertices=[]
    for vertice in upperVertices:
        P=vertice   
        AllUpperVertices += [tuple(P)]           
    
    allEdges=[]
    for edge in edges:
        i=int(edge[0])
        j=int(edge[1])
        allEdges += [(i, j)] 
    print(allEdges)  
    
    allfaces=[]  
    for face in faces:
        i=int(face[0])
        j=int(face[1])
        k=int(face[2])
        l=int(face[3])
        allfaces += [(i, j, k, l)] 
    print(allfaces)         
    mesh = bpy.data.meshes.new(name) 
    mesh.from_pydata(AllVertices, allEdges, allfaces)
    
    mesh2 = bpy.data.meshes.new('iforiofo') 
    mesh2.from_pydata(AllUpperVertices, allEdges, allfaces)    
    obj_network = bpy.data.objects.new(name, mesh)
    coll.objects.link(obj_network)    
    obj_network1 = bpy.data.objects.new(name, mesh2)
    coll.objects.link(obj_network1)             
    verts=[]
    edge=[]    
    for n, e in enumerate(edges):
        i=int(e[0])
        j=int(e[1])
        #Pi=vertices[i,:]
        #Pj=vertices[j,:]       

      
        
        
      
    #mesh.from_pydata(verts, edge, [])
    # Create Object and link to scene

    #bpy.context.scene.collection.objects.link(obj_network)
    #print(Pi)
    
main(10, 0.03, False)