Introduction.
=============

 This is a fairly simple and somewhat crooked script designed for visual analysis of the distributionof the electrostatic potential over various electron density isosurfaces. 
 An external Multiwfn (http://sobereva.com/multiwfn/) program is used to generate ED and ESP. (J. Comput. Chem., 33, 580-592 (2012))


Requirements and Installation.
==============================

 In order to run the script, you need:

 1) Python.
    Im using Python3.8 x64 and strongly advise you to download it from the official website (https://www.python.org/).

 2) (The Visualization Toolkit) VTK’s python interface.
    The current version (9.1.0) of VTK, available via pip, unfortunately completely breaks the Mayavi installation process 
    (therefore, the method of installation via pip proposed on the official Mayavi website is likely to lead to an error, but you might still try).
    You need to use version 9.0.1. 
    Just in case, I attach a whl file for win64 and linux, which will install the necessary working version. Just write on the command line:
    "pip install vtk-9.0.1-cp_38-cp_38-win_amd64.whl"

 3) PyQt5; PyOpenGl; Mayavi.
    It is installed easily using pip, which comes bundled with python, you just need to write on the command line:
    "pip install PyQt5 PyOpenGl Mayavi"
    and the necessary components themselves will be downloaded and installed

 4) Multiwfn.
    Any version will do, you can download it for free from the official website



Necessary code modification.
============================

 "Multiwfnpath" variable should be changed if the path to the multiwfn differs from 'Multiwfn' (if you did not specify the path or for some other reason)


  
Working with the script. Part 1, generating cubes.
==================================================
 
 Supported wavefunction files: .fchk; .wfn; .wfx

 Main launch method: python VisMap4.0.py NameOfWFNFile -nproc=N -mode=old/new -vis=y/n

 Decryption and default values (if the string is not in the input):
  
  -nproc=N
   (The number of threads that multivan will use to calculate cubes/crit.points on the surface. 
   Default value is 4)
  
  -mode=old/new
   (If the cubes have already been generated, they are saved as NameOfWFNFile (without extension) + _ESP.cub (for ESP)
                                                               NameOfWFNFile (without extension) + _Dens.cub (for density)
   old mode first tries to find these files and, if they are not found, calls the Multiwfn to generate
   new mode does not check the folder and immediately calls the Multiwfn to generate
   Default value is old)

  -vis=y/n
   (To call the visualizer or not. It is necessary if you want to calculate cubes on some external server and do not want to put a visualizer there. 
   If mode n is selected, then the function with the visualizer simply will not be called.
   Default value is y)



Working with the script. Part 2, visualizer and CP of ESP on the surface.
=========================================================================
 
 Commands of the "interactive" mode:
 
 isosurf value (i.e. float number)
  (By default, a map is built for the EP 0.001 isosurface. You can specify some other number - then the isosurface is regenerated, this is fast.)

 gen
  (Looking for a file NameOfWFNFile (without extension)_sa_*isoval*.txt
   If it is not there, it calls Multiwfn to generate it.
   Reads all the points and adds them to the visualizer.)
 
 scan
  (Looking for a file NameOfWFNFile (without extension)_sa_*isoval*.txt
   If it is not there, it calls Multiwfn to generate it.
   Reads all the points, checks if they are within a radius of 5.0 a.u. from
   the atoms (H, O, F, S, Cl, Se, Br, Te, I, Po, At) and adds them to the visualizer.)

 kill N M
  (Since there are too many dots even with a scan, I added the ability to erase dots.
   When calling this command, all points with an ESP value equal to N ± M are destroyed.)

 exit
  (Exit from the "interactive" mode – you can work with Mayavi window.)



Examples - see the attached pdf file.
=====================================








