import pyvista as pv

points1 = [[1,0,0],[1,1,0],[0,1,0],[0,0,0]]
points2 = [[0.25,0.5,0],[0.25,0.25,0]]
points3 = [[0.9, 0.05, 0.0],[0.9, 0.6, 0.0],[0.7, 0.7, 0.0]]

lines1 = [5, 0, 1, 2, 3, 0]
lines2 = [2, 0, 1]
lines3 = [3, 0, 1, 2]

rect1 = pv.PolyData(points1, lines=lines1)

segment1 = pv.PolyData(points2, lines=lines2)
segment2 = pv.PolyData(points3, lines=lines3)

PSLG = rect1 + segment1 + segment2

tess = PSLG.delaunay_2d(edge_source=PSLG)

PSLG.save("PolyData_Files\segmented_square_pslg.vtk")
tess.save("PolyData_Files\segmented_square_cd.vtk")