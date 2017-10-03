# IVA_python
Inversion Velocity Analysis Code in 1D/2D

This is a pythonified version of the Code I wrote during my PhD thesis [3] (full Fortran 90). The most computer-expensive part will remain in Fortran and I need to make all f2py, python, mpi work together...
Work in progress...

This work is in the very general framework of seismic imaging. The particular imaging technique studied here is *Migration Velocity Analysis* using model extension with the subsurface-offset as extension parameter \cite{symes}.

The main idea tested in this code is to replace the migration step by an inversion. Two strategies are (or should be in a near future) implemented :
* direct inversion: an asymptotic inverse of the extended Born modelling operator is used to derive an extended reflectivity image. Here, the inverse proposied in [1] is studied.
* iterative inversion: the extended reflectivity model is obtained as solution of an inner inverse problem (namely *iterative migration*, or *least-squares migration*.


[1] Chauris and Cocher, 2017.

[2] Cocher et al, 2017.

[3] Cocher, 2017.
