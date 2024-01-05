import bpy
import csv
import numpy as np
from numpy import genfromtxt
import mathutils
import math
from mathutils import Matrix 



def makeMaterial(name, diffuse, metallic, roughness):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = diffuse
    mat.metallic = metallic
    mat.roughness = roughness    
    return mat

def setMaterial(ob, mat):
    me = ob.data
    me.materials.append(mat)


def importData():
    lowerVertices = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\VerticesLowerBounds.csv', delimiter = ',')
    upperVertices = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\VerticesUpperBounds.csv', delimiter = ',')
    edges = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\Edges.csv', delimiter = ',')
    faces = genfromtxt('E:\Archivio Lavori\Sviluppo\TNA\Faces.csv', delimiter = ',')

    #color = genfromtxt('Member-Color.csv', delimiter = ';')
    return lowerVertices, upperVertices, edges, faces #, color

def calculateRotationMatrix(Pi, Pj):
    R = np.identity(4)
    ix = Pi[0]
    iy = Pi[1]
    iz = Pi[2]
    
    jx = Pj[0]
    jy = Pj[1]
    jz = Pj[2]   
    
    dx = jx-ix
    dy = jy-iy
    dz = jz-iz
    
    length = math.sqrt(dx**2+dy**2+dz**2) #Magnitude of vector (lehgth of member)
    
    
    if (abs(dx)<0.001 and abs(dy)<0.001):
        #Element is vertical - offset in positive global x to define local z-x plane
        i_offset = np.array([ix+1, iy, iz]) #Offset node i by 1m in positive global x-direction
        j_offset = np.array([jx+1, jy, jz]) #Offset node j by 1m in positive global x-direction
    else:
        #Element is not vertical - offset in positive global z to define local x-z plane
        i_offset = np.array([ix, iy, iz+1]) #Offset node i by 1m in positive global z-direction
        j_offset = np.array([jx, jy, jz+1]) #Offset node j by 1m in positive global z-direction
    
    node_k = i_offset + 0.50*(j_offset - i_offset) # Point in the local x-y plane
    
        
    #Local z-vector in global RF running along the member
    local_z_vector = np.array(Pj) - np.array(Pi) # Vector along local z-axis
    local_z_unit = local_z_vector/length #Local unit vector defininf local z-axis
    
    #Local x-vector in global RF using Gram-Schmidt process
    vector_in_plane = node_k - Pi #Vector in the x-y plane
    local_x_vector = vector_in_plane - np.dot(vector_in_plane, local_z_unit)*local_z_unit #local x-vector in RF (Gram-Schmidt)
    magX = math.sqrt(local_x_vector[0]**2+local_x_vector[1]**2+local_x_vector[2]**2) #length of local x-vector
    local_x_unit = local_x_vector/magX #Local unit vector defining local x-axis  
    
    #Local y-vector in global RF using matrix cross product
    local_y_unit = np.cross(local_z_unit, local_x_unit) ##Local unit vector defining local y-axis  
    
    #combine reference frame into standard rotation matrix for the element vector x, y, z: column, 1, 2, 3
    rotationMatrix = np.array([local_x_unit, local_y_unit, local_z_unit]).T
    
    R[0:3, 0:3] = rotationMatrix
    
    
    return Matrix(R), length


#Function to return translation component of transformation matrix
def calculateTranslationMatrix(Pi, Pj):
    offset = Pi + 0.50*(Pj - Pi) #translation to centre of eleemnt
    transMatrix = Matrix.Translation(offset) #Construct a 4x4 transformation matrix that encodes this transformation
    return transMatrix

def generateElement(scale, length, num, forceMag = 1000.): #, color[n, 4])
    #Add a mesh primitive
    bpy.ops.mesh.primitive_cylinder_add(
        vertices = 4,
        radius = scale,
        depth = length,
        align = 'WORLD',
        location = (0., 0., 0.),
        scale = (1, 1, 1))
    bpy.ops.object.shade_smooth()
    cylinder = bpy.context.object #Get handle to object
    cylinder.name = 'Element ' + str(num)+'_'+str(round(forceMag/1000.)) + 'kN' #Rename object in outliner
    
    return cylinder

def generateNode(v, scale, num):
    #Add a mesh primitive
    bpy.ops.mesh.primitive_ico_sphere_add(
        radius = 2.*scale,
        subdivisions = 1,
        enter_editmode=False,
        align = 'WORLD',
        location = (v[0], v[1], v[2]),
        scale = (1, 1, 1))
    bpy.ops.object.shade_smooth()
    ico = bpy.context.object #Get handle to object
    ico.name = 'Node ' + str(num) #Rename object in outliner
    return ico


def combineElements(prefix, collection):
    #https://blender.stackexchange.com/questions/13986/how-to-join-objects-with-python
    obs=[]
    #collection = bpy.data.collections["COLORED-NETWORK"]
    for ob in collection.objects:# bpy.context.scene.objects:
        if ob.name.startswith(prefix): obs.append(ob)
    
    print('\n',obs) 
    #active_object = obs[0]
     
    #selected_objects = obs
    
    with bpy.context.temp_override(active_object=obs[0], selected_editable_objects=obs):
        bpy.ops.object.join() 

    return obs[0]
    
    
def main(chunkSize, elementScale, redraw):
    lowerVertices, upperVertices, edges, faces = importData() #, color
    #print(edges)
    #print(vertices)
    #print(faces)

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
    #print(allEdges)  
    
    allfaces=[]  
    for face in faces:        
        i=int(face[0])
        j=int(face[1])
        k=int(face[2])
        l=int(face[3])
        allfaces += [(i, j, k, l)] 
    #print(allfaces)   
    lowColl=bpy.data.collections.new("lowerBounds")
    bpy.context.scene.collection.children.link(lowColl)
    
    
    name = 'lowerBoundsNetwork'    
    LOWnetmesh = bpy.data.meshes.new(name) 
    LOWnetmesh.from_pydata(AllUpperVertices, allEdges, [])  
    obj_LOWnetwork = bpy.data.objects.new(name, LOWnetmesh)
    lowColl.objects.link(obj_LOWnetwork)      
   
   
    name = 'loweBoundsSurface'      
    LOWsurfmesh = bpy.data.meshes.new(name) 
    #LOWmesh.from_pydata(AllVertices, allEdges, allfaces)
    LOWsurfmesh.from_pydata(AllVertices, [], allfaces)
    obj_LOWsurface = bpy.data.objects.new(name, LOWsurfmesh)
    lowColl.objects.link(obj_LOWsurface)        
    
    
    
    upperColl=bpy.data.collections.new("upperBounds")
    bpy.context.scene.collection.children.link(upperColl)
    
    name = 'upperBoundsNetwork' 
    UPnetmesh = bpy.data.meshes.new(name) 
    UPnetmesh.from_pydata(AllUpperVertices, allEdges, [])  
    obj_UPnetwork = bpy.data.objects.new(name, UPnetmesh)
    upperColl.objects.link(obj_UPnetwork)         
    
    name = 'upperBoundsSurface' 
    UPsurfmesh = bpy.data.meshes.new(name)
    UPsurfmesh.from_pydata(AllUpperVertices, [], allfaces)  
    obj_UPsurfice = bpy.data.objects.new(name, UPsurfmesh)
    upperColl.objects.link(obj_UPsurfice)
    
                 

    #DegreeTutors
    
    NetworkCollectionName='COLORED-NETWORK'
    Active_collection=bpy.context.view_layer.active_layer_collection.collection
    coll=bpy.data.collections.new(NetworkCollectionName) #Create a new collection
    bpy.context.scene.collection.children.link(coll) #link collection to master scene collection #???
    
    collBranches=bpy.data.collections.new("Branches") 
    bpy.data.collections[NetworkCollectionName].children.link(collBranches)

    collNodes=bpy.data.collections.new("Nodes") 
    bpy.data.collections[NetworkCollectionName].children.link(collNodes)
    
    #Generate elements
    
    #Create material
    matName = "elementMat"
    c = [0.95, 0.05, 0.05]
    eleMat = makeMaterial(matName, (c[0], c[1], c[2], 1), 0.5, 0.3) #single material for all element
    
    cnt=0
    for n, e in enumerate(edges):
        i=int(e[0])
        j=int(e[1])
        Pi=lowerVertices[i,:]
        Pj=lowerVertices[j,:]       

        #Calculate position and orientation
        R, length = calculateRotationMatrix(Pi, Pj)
        T = calculateTranslationMatrix(Pi, Pj)
        transformationMatrix  = T @ R
        
        #Generate eleemnt and apply tranformation
        element = generateElement(elementScale, length, n) #, color[n, 4])
        element.matrix_world = transformationMatrix
        
        #Link object with 'COLORED-NETWORK' collection
        collBranches.objects.link(element) 
        Active_collection.objects.unlink(element)  
        #if n==0: ob_to_unlink = element
        
        #Apply material to current element
        setMaterial(element, eleMat)        
        cnt+=1
        #print(cnt)        
        if cnt == chunkSize:# and n<=len(edges)- chunkSize:
            cnt=0            
            group = combineElements('Element', collBranches)
            group.name ='Group_' + group.name  
            print('\n', f'Completed element: {n+1} of {len(edges)}')
            if redraw:
                bpy.ops.wm.redraw_timer(type = 'DRAW_WIN_SWAP', iterations = 1)#Update 3D viewport 
            #bpy.context.scene.collection.objects.unlink(group)
   
    #Generate nodes same pattern as above elements)
    matName = "nodeMat"
    nodemat= makeMaterial(matName, (0.2, 0.8, 0.8, 1), 0.8, 0.5) #single material for all nodes
    cnt=0
    for n, v in enumerate(lowerVertices):
        node = generateNode(v, elementScale, n)        
        collNodes.objects.link(node)   
        Active_collection.objects.unlink(node) 
        setMaterial(node, nodemat)
        cnt+=1
        if cnt == chunkSize:# and n<=len(lowerVertices)- chunkSize:
            cnt=0            
            group = combineElements('Node', collNodes)
            group.name ='Group_' + group.name
            print('\n', f'Completed node: {n+1} of {len(lowerVertices)}')
            if redraw:
                bpy.ops.wm.redraw_timer(type = 'DRAW_WIN_SWAP', iterations = 1)#Update 3D viewport 
            #bpy.context.scene.collection.objects.unlink(group)
    
main(20, 0.005, False)