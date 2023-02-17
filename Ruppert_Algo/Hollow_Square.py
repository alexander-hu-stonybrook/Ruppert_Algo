import pyvista as pv

points1 = [[1,0,0],[1,1,0],[0,1,0],[0,0,0]]
lines1 = [5, 0, 1, 2, 3, 0]

points4 = [[0.25,0.5,0],[0.25,0.25,0],[0.5,0.25,0],[0.5,0.5,0]]
faces4 = [4, 0, 3, 2, 1]
points5 = [[0.7,0.6,0],[0.7,0.7,0],[0.6,0.7,0],[0.6,0.6,0]]
faces5 = [4, 0, 3, 2, 1]

rect1 = pv.PolyData(points1, lines=lines1)
rect4 = pv.PolyData(points4, faces4)
rect5 = pv.PolyData(points5, faces5)

PSLG = rect4 + rect1 + rect5

tess = PSLG.delaunay_2d(edge_source=PSLG)

PSLG.save("PolyData_Files\hollow_square_pslg.vtk")
tess.save("PolyData_Files\hollow_square_cd.vtk")