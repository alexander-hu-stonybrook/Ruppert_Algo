
import pyvista as pv

convexpoints = [
        [1,1,0],
        [4,-1,0],
        [10,1,0],
        [9,6,0],
        [4,8,0],
        [0,7,0],
        [1,6,0],
        [5,6,0],
        [4,2,0]
    ]
convexface = [9, 0, 1, 2, 3, 4, 5, 6, 7, 8]
convex = pv.PolyData(convexpoints,convexface)

trianglepoints = [[7,3,0],[7,5,0],[8,5,0],[9,3,0]]
triangleface = [4, 0, 1, 2, 3]
triangle = pv.PolyData(trianglepoints,triangleface)

line1points = [[3,7,0],[5,7,0],[6,6,0]]
line1edge = [3, 0, 1, 2]
line1 = pv.PolyData(line1points, lines=line1edge)

line2points = [[4,1,0],[6,3,0]]
line2edge = [2, 0, 1]
line2 = pv.PolyData(line2points, lines=line2edge)

PSLG = convex + triangle + line1 + line2
PSLG.plot(cpos="xy", show_edges=True)
print(PSLG.points)
print(PSLG.faces)
print(PSLG.lines)

tess = PSLG.delaunay_2d(edge_source=PSLG)
tess.plot(cpos="xy", show_edges=True)

PSLG.save("PolyData_Files\concave_poly_pslg.vtk")
tess.save("PolyData_Files\concave_poly_cd.vtk")

'''
readpslg = pv.read("PolyData_Files\convex_poly_pslg.vtk")
readpslg.plot(cpos="xy", show_edges=True)
print(readpslg.points)
print(readpslg.faces)
print(readpslg.lines)

readcd = pv.read("PolyData_Files\convex_poly_cd.vtk")
readcd.plot(cpos="xy", show_edges=True)
'''