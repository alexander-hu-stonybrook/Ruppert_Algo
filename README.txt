Alexander Hu
ID: 110544388
CSE528 Final Project Code Execution README

Requirements to run code:

1. Python3 (I used Python3.9)
2. PyVista installed as a package in Python3 (pip install pyvista should install the package)

Code was written and tested using Visual Studio Community Version 2022 via the Python IDE option.
The solution file is left in the zip folder in case you wish to access it that way.
However, it should be runnable with any python method.

To run code, simply run Ruppert_Algo.py.
All the other files were used to generate the test data.

When running code, simply close the popup windows that appear to show the 
next one. There should be 3 pop up windows:

1. The PSLG input.
2. The constrained delaunay mesh using the PSLG.
3. The resulting refined delaunay mesh.

If you wish to change which test data to use in the PolyData_Files folder,
open Ruppert_Algo.py and scroll down to the bottom, then change the string
argument "file_name" in SetupPSLGCD("file_name") to the relevant file in
the PolyData_Files folder.