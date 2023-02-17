from pickle import NEWTRUE
from tkinter import INSERT
import pyvista as pv
import numpy as np
import math

'''
Please go to the bottom section, "Script Execution," to change how the code is executed.
'''

#####################################
# Helper Functions/Global Variables #
#####################################

lawsonStack = []
skinnyQueue = []
segmentQueue = []

facesDict = {}
edgesDict = {}

edgesPSLG = set()
pointsList = []

#Margin of error to account for floating point errors
#This value can be adjusted based on the dimensions of mesh
ERRORMARGIN = 0.00001

#Need a way to load PSLG info from a file?
#But we can worry about that later


#May need to check this function again in the future, for now this is acceptable
def EdgesPSLGFiller(faces, edges):
    facenum = 0
    prevpoint = -1
    firstpoint = -1
    for n in faces:
        #print(n)
        #print(facenum)
        if facenum == 0:
            if firstpoint != -1:
               e = frozenset((prevpoint, firstpoint))
               edgesPSLG.add(e)
            facenum = n
            prevpoint = -1
            firstpoint = -1

        else:

            if prevpoint >= 0:
                e = frozenset((prevpoint,n))
                edgesPSLG.add(e)
            if firstpoint == -1:
                firstpoint = n

            facenum -= 1
            prevpoint = n
    
    if firstpoint != -1:
        e = frozenset((prevpoint, firstpoint))
        edgesPSLG.add(e)

    
    facenum = 0
    prevpoint = -1
    for n in edges:
        if facenum == 0:
            facenum = n
            prevpoint = -1
        else:
            if prevpoint >= 0:
                e = frozenset((prevpoint, n))
                edgesPSLG.add(e)
            prevpoint = n
            facenum -= 1

def InitialFacesEdgesPointsSetup(c_delaunay):
    cd_faces = c_delaunay.faces.reshape(-1,4)
    for p in c_delaunay.points:
        pointsList.append(p.tolist())

    for f in cd_faces:
        triangle = tuple(f)
        facesDict[triangle] = 1
        InsertFacesIntoEdges(triangle)

def InsertFacesIntoEdges(f):
    global edgesDict
    if frozenset((f[1],f[2])) in edgesDict.keys():
        edgesDict[frozenset((f[1],f[2]))].append(f)
    else:
        edgesDict[frozenset((f[1],f[2]))] = [f]

    if frozenset((f[2],f[3])) in edgesDict.keys():
        edgesDict[frozenset((f[2],f[3]))].append(f)
    else:
        edgesDict[frozenset((f[2],f[3]))] = [f]

    if frozenset((f[1],f[3])) in edgesDict.keys():
        edgesDict[frozenset((f[1],f[3]))].append(f)
    else:
        edgesDict[frozenset((f[1],f[3]))] = [f]

def EdgesEncroached(edge, point):
    global pointsList

    p1,p2 = ExtractEdgePoints(edge)
    mp = Midpoint(p1,p2)
    radsq = DistSquared(p1,mp)

    #If a point is given, it's a point trying to be inserted
    #Otherwise, it's the initial check
    if point != None:
        if DistSquared(point,mp) <= radsq:
            return True
        else:
            return False
    else:
        counter = 0
        for p in pointsList:
            if counter in edge:
                pass
            elif DistSquared(p,mp) <= radsq:
                #print("point has distance smaller than the radius of the circle")
                #For the initial check
                #The encroaching point must be visible from the edge
                if PointVisibleToEdge(edge,p,mp):
                    #print("point is visible to edge")
                    return True
            counter += 1

        return False

def PointVisibleToEdge(edge, point, midpoint):

    global facesDict, edgesDict, edgesPSLG, pointsList

    ptest1, ptest2 = ExtractEdgePoints(edge) 
    #print("Checking original edge given to PointVisibleFromEdge")
    #print((ptest1,ptest2))

    for e in edgesPSLG:
        if e != edge:
            p1, p2 = ExtractEdgePoints(e)
            l1inf = False
            l2inf = False
            if p1[0] == p2[0]:
                l1inf = True
            if point[0] == midpoint[0]:
                l2inf = True
            if l1inf and l2inf:
                continue

            #print("Checking PointVisibleFromEdge Function")
            #print((p1, p2, point, midpoint))

            smallx = min(p1[0],p2[0])
            largex = max(p1[0],p2[0])
            smally = min(p1[1],p2[1]) 
            largey = max(p1[1],p2[1])

            slope1 = None
            slope2 = None
            intercept1 = None
            intercept2 = None

            if not l1inf:
                slope1 = (p1[1] - p2[1]) / (p1[0] - p2[0])
                intercept1 = p1[1] - (slope1 * p1[0])
            if not l2inf:
                slope2 = (point[1]-midpoint[1]) / (point[0]-midpoint[0])            
                intercept2 = point[1] - (slope2 * point[0])

            if slope1 != slope2:
                if slope1 == None:
                    if min(point[0],midpoint[0]) < p1[0] and max(point[0],midpoint[0]) > p1[0]:
                        yinter = slope2 * smallx + intercept2
                        if yinter > smally and yinter < largey:
                            return False
                elif slope2 == None:
                    if smallx < point[0] and largex > point[0]:
                        yinter = slope1 * point[0] + intercept1
                        if yinter > min(point[1],midpoint[1]) and yinter < max(point[1],midpoint[1]):
                            return False
                else:
                    xinter = (intercept2 - intercept1) / (slope1 - slope2)
                    if xinter > smallx and xinter < largex:
                        if xinter > min(point[0],midpoint[0]) and xinter < max(point[0],midpoint[0]):
                            return False
            
    return True

def ExtractEdgePoints(e):
    global pointsList
    pid1 = -1
    pid2 = -1
    for c in e:
        if pid1 >= 0:
            pid2 = c
        else:
            pid1 = c
    p1 = pointsList[pid1]
    p2 = pointsList[pid2]
    return p1,p2

def ExtractEdgeIDs(e):
    pid1 = -1
    pid2 = -1
    for c in e:
        if pid1 >= 0:
            pid2 = c
        else:
            pid1 = c

    return pid1,pid2

def Midpoint(p1,p2):
    return [(p1[0]+p2[0])/2, (p1[1]+p2[1])/2, 0]

def DistSquared(p1,p2):
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2

def checkCCW(pid1,pid2,pid3):
    global pointsList

    p1 = pointsList[pid1] 
    p2 = pointsList[pid2] 
    p3 = pointsList[pid3] 
    ax = p1[0]
    ay = p1[1]
    bx = p2[0]
    by = p2[1]
    cx = p3[0]
    cy = p3[1]

    if (bx - ax) * (cy - ay) - (cx - ax) * (by - ay) > 0:
        return True
    else:
        return False

def CircumcenterCoords(t):
    ax = pointsList[t[1]][0]
    ay = pointsList[t[1]][1]
    bx = pointsList[t[2]][0]
    by = pointsList[t[2]][1]
    cx = pointsList[t[3]][0]
    cy = pointsList[t[3]][1]

    '''
    print(t) 
    print(ax)
    print(ay) 
    print(bx) 
    print(by) 
    print(cx)  
    print(cy)
    '''
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    #print(d)
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d

    return ux, uy

def AdjTirangles(t):
    #t is a triangle tuple/array with 4 elements
    global facesDict, edgesDict

    e1 = frozenset((t[1],t[2]))
    e2 = frozenset((t[2],t[3]))
    e3 = frozenset((t[1],t[3])) 

    adjarr = []
    edgearr = [e1,e2,e3]

    for e in edgearr:
        edgetri = edgesDict[e] 
        for et in edgetri:
            if et != t and facesDict[et] == 1:
                adjarr.append((e,et))

    #return an array with elements as (edgeshared, adjtriangle)
    return adjarr

def DisplayMeshOutput():
    global facesDict
    facesarr = []
    for t in facesDict:
        #print(t) 
        #print(facesDict[t])
        if facesDict[t] == 1:
            for i in t:
                facesarr.append(i)
            

    #print(facesarr)

    displaymesh = pv.PolyData(pointsList, facesarr) 
    displaymesh.plot(cpos='xy', show_edges=True)

def CheckDelaunay():
    for t in facesDict:
        if facesDict[t] == 1:
            ccx,ccy = CircumcenterCoords(t)
            radsq = DistSquared((ccx,ccy),pointsList[t[1]])
            for p in pointsList:
                if DistSquared((p[0],p[1]),(ccx,ccy)) < radsq - ERRORMARGIN:
                    return False 
    return True
                
def QueueSkinnyTriangles():
    global facesDict, skinnyQueue, pointsList
    bbound = math.sqrt(2)
    for f in facesDict:
        if facesDict[f] == 1:
            ccx,ccy = CircumcenterCoords(f)
            radsq = DistSquared((ccx,ccy), pointsList[f[1]]) 
            e1 = (pointsList[f[1]],pointsList[f[2]])
            e2 = (pointsList[f[2]],pointsList[f[3]])
            e3 = (pointsList[f[1]],pointsList[f[3]])
            minseg = min(DistSquared(e1[0],e1[1]),DistSquared(e2[0],e2[1]),DistSquared(e3[0],e3[1]))
            bratio = math.sqrt(radsq/minseg)
            if bratio > bbound:
                skinnyQueue.append(f)

def EdgeInBadTriangle(triori,edge,badtriangles):
    global edgesDict, facesDict

    for f in edgesDict[edge]:
        if facesDict[f] == 1 and f != triori:
            if f in badtriangles:
                return True

    return False

def SetupPSLGCD(dataname):
    str1 = "PolyData_Files\\" + dataname + "_pslg.vtk"
    str2 = "PolyData_Files\\" + dataname + "_cd.vtk"

    PSLG = pv.read(str1)
    PSLG.plot(cpos='xy', show_edges=True)
    tess = pv.read(str2)
    tess.plot(cpos='xy', show_edges=True)

    return PSLG, tess

############################################
# Lawson/Bowyer-Watson Algorithm Functions #
############################################

def LawsonFlip(triangle, ccenter, radsq, adjtriarr, checkadj = False):
    #We will be using the ERRORMARGIN here as it is possible that
    #due to floating point errors, an edge is flipped that results in
    #obviously incorrect triangles.
    #Thus, to avoid this we will give some leeway for points that 
    #potentially lie close enough to be on the circumcircle and therefore
    #is a valid configuration.

    #Return True if a flip happened, Return False if it did not

    global lawsonStack, segmentQueue, facesDict, edgesDict, edgesPSLG, pointsList

    #print("Checking adjtriarr input")
    #print(adjtriarr)
    for sharededge, adjtri in adjtriarr:

        #We never want to flip an edge that's a constraint
        #So if this somehow is the case, skip
        if sharededge in edgesPSLG:
            continue

        checkpointid = -1
        #print("Checking adjtri outside accessing adjtri index")
        #print(adjtri)
        if adjtri[1] not in sharededge:
            checkpointid = adjtri[1] 
        elif adjtri[2] not in sharededge:
            checkpointid = adjtri[2] 
        else:
            checkpointid = adjtri[3]

        connecterid = -1
        if triangle[1] not in sharededge:
            connecterid = triangle[1] 
        elif triangle[2] not in sharededge:
            connecterid = triangle[2] 
        else:
            connecterid = triangle[3]
        
        distsq = DistSquared(pointsList[checkpointid],ccenter) 

        #If circumcircle is encroached by point
        if distsq < radsq - ERRORMARGIN:
            #Remove the offending triangles
            facesDict[triangle] = 0
            facesDict[adjtri] = 0

            #Perform Lawson flip
            for p in sharededge:

                #Add new triangles
                newtri = (3,checkpointid,connecterid,p)
                if not checkCCW(checkpointid,connecterid,p):
                    newtri = (3,p,connecterid,checkpointid)
                facesDict[newtri] = 1
                InsertFacesIntoEdges(newtri)

                #Add new triangles to LawsonStack
                lawsonStack.append(newtri)

            #report that a flip took place
            return True
    
    #If currently checking original triangle, also have to check adjacents
    #If one return True, this entire function returns True
    #If checking adjacents, check ends here.
    if not checkadj:
        for sharededge, adjtri in adjtriarr:
            #print("Checking adjacent triangle lawsons")
            newadjtriarr = [(sharededge, triangle)]
            #print(newadjtriarr)
            newccx, newccy = CircumcenterCoords(adjtri) 
            newradsq = DistSquared((newccx,newccy),pointsList[adjtri[1]])
            checkflip = LawsonFlip(adjtri,(newccx,newccy),newradsq,newadjtriarr,True)
            if checkflip:
                return True

    return False

def LawsonAlg(point = None):
    
    global lawsonStack, segmentQueue, facesDict, edgesDict, edgesPSLG, pointsList

    #point == None is for the initial first run

    notempty = True
    firstpass = True
    
    counter = 0
    while notempty:
        for e in edgesPSLG:

            #Check if the edge is encroached
            #If it is, add to the segment split queue
            if EdgesEncroached(e, point):
                segmentQueue.append(e)

        #If after checking no segments need to be split,
        #set flag to False, which ends function

        if len(segmentQueue) == 0:
            notempty = False

            #Given a point, it must not encroach any edges
            #This only works if that's the case
            #If it encroached, firstpass will be False
            #once the segment queue is emptied
            if firstpass and point != None:
                return True

        for e in segmentQueue:
            #Perform midpoint insertion and LawsonAlg here

            #Find the midpoint
            ep1,ep2 = ExtractEdgePoints(e) 
            mp = Midpoint(ep1,ep2) 

            #Add the midpoint to the pointsList
            pointsList.append(mp)

            #Get the midpoint id
            mpid = len(pointsList) - 1

            #Add subsegments to edgesPSLG
            #Remove the original segment from edgesPSLG

            pid1,pid2 = ExtractEdgeIDs(e)

            edgesPSLG.remove(e)
            edgesPSLG.add(frozenset((mpid,pid1)))
            edgesPSLG.add(frozenset((mpid,pid2)))

            #Remove old triangles and draw new ones
            #Track the boundaries here
            boundaries = []

            #Find the old triangles, they should be the only ones active
            for t in edgesDict[e]:
                if facesDict[t] == 1:
                    te1 = frozenset((t[1],t[2]))
                    te2 = frozenset((t[2],t[3]))
                    te3 = frozenset((t[1],t[3])) 
                    if te1 != e:
                        boundaries.append(te1) 
                    if te2 != e:
                        boundaries.append(te2)
                    if te3 != e:
                        boundaries.append(te3)
                    
                    #Remove the old triangles
                    facesDict[t] = 0

            #Iterate through the list and create new triangles
            for b in boundaries:
                pid1, pid2 = ExtractEdgeIDs(b)

                #Ensure new triangle is ordered CCW
                newtri = None
                if checkCCW(pid1,pid2,mpid):
                    newtri = (3, pid1, pid2, mpid)
                else:
                    newtri = (3, mpid, pid2, pid1)

                #Add new triangle to the Lawson algorithm stack
                lawsonStack.append(newtri)
                
                #Perform maintenance
                facesDict[newtri] = 1
                InsertFacesIntoEdges(newtri)

            #Check each triangle and perform Lawson flip if necessary
            while len(lawsonStack) > 0:
                #Extract a triangle
                t = lawsonStack.pop()

                #Check if triangle still exists in the first place
                #It may have been removed in a flip
                if facesDict[t] == 1:
                    #Find the circumcenter
                    ccx, ccy = CircumcenterCoords(t) 
                    #Find the radius by connecting circumcenter to one of the triangle points
                    radsq = DistSquared((ccx,ccy),pointsList[t[1]])

                    #Find all triangles surrounding the given triangle
                    adjtri = AdjTirangles(t) 

                    #Attempt to perform Lawson Flip
                    LawsonFlip(t,(ccx,ccy),radsq,adjtri)

                    #Need to write functions for:

        #clear out segmentQueue after reaching the end
        #if need to split more, loop will add more in
        segmentQueue = []

        #DisplayMeshOutput()

        firstpass = False
        
    return False

def BowyerWatsonAlg(point):

    global facesDict,pointsList

    #Store triangles to remove here
    badtriangles = []

    #Add triangles to badtriangles
    #If their circumcircle is encroached by the point
    for f in facesDict:
        if facesDict[f] == 1:
            ccx,ccy = CircumcenterCoords(f)
            radsq = DistSquared((ccx,ccy),pointsList[f[1]])
            if radsq > DistSquared((ccx,ccy),point): 
                badtriangles.append(f)

    #Find the bounding edges of the bad triangles
    boundaries = []

    for t in badtriangles:
        e1 = frozenset((t[1],t[2]))
        e2 = frozenset((t[2],t[3]))
        e3 = frozenset((t[1],t[3]))

        edges = [e1,e2,e3]

        for e in edges:
            if not EdgeInBadTriangle(t,e,badtriangles):
                boundaries.append(e)
        
    #Remove all bad triangles from facesDict
    for t in badtriangles:
        facesDict[t] = 0

    #Create and add new triangles and new point
    pointsList.append(point)
    pidp = len(pointsList) - 1
    #print("checking boundaries extracted in bowyer watson")
    #print(boundaries)
    for b in boundaries:
        pid1,pid2 = ExtractEdgeIDs(b)
        newtri = (3, pid1, pid2, pidp)
        if not checkCCW(pid1,pid2,pidp): 
            newtri = (3, pidp, pid2, pid1)
        #print(newtri)
        facesDict[newtri] = 1
        InsertFacesIntoEdges(newtri)

#######################
# Ruppert's Algorithm #
#######################

def RuppertAlg(pslg_data,c_delaunay_mesh):
    global skinnyQueue, segmentQueue, facesDict

    #Set up the necessary data structures
    EdgesPSLGFiller(pslg_data.faces, pslg_data.lines)
    InitialFacesEdgesPointsSetup(c_delaunay_mesh)

    #Perform initial segment encroachment check
    LawsonAlg() 

    #print("Ended initial lawson alg in ruppert alg")

    #Remove skinny triangles
    while True:
        QueueSkinnyTriangles()
        if len(skinnyQueue) == 0:
            #If after queueing triangles the queue is still empty
            #Function is done, exit
            break
        else:
            #Iterate through the skinnyQueue
            counter = 0
            while counter < len(skinnyQueue):
                if facesDict[skinnyQueue[counter]] == 1:
                    #Create the center point and check if it encroaches any segments
                    ccx, ccy = CircumcenterCoords(skinnyQueue[counter])
                    newpoint = [ccx,ccy,0.0]
                    #print("Checking newpoint in bowyer-watson alg")
                    #print(newpoint)
                    insertpoint = LawsonAlg(newpoint)

                    #If it encroaches, LawsonAlg handles situation
                    #Do not advance as this triangle may still be skinny
                    if insertpoint:
                        #Insert the point into the mesh using Bowyer-Watson algorithm
                        #This will for sure eliminate the skinny triangle, so advance position
                        BowyerWatsonAlg(newpoint)
                        #print("Added new point with Bowyer-Watson")
                        #DisplayMeshOutput()
                        counter += 1
                else:
                    #If triangle no longer exists, ignore and move on
                    counter += 1

            #At the end of the skinnyQueue, remove all items
            skinnyQueue = []

####################
# Script Execution #
####################

PSLG, tess = SetupPSLGCD("segmented_square")

RuppertAlg(PSLG,tess)
DisplayMeshOutput()
print(CheckDelaunay())

