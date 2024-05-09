# NH-TOPOPT-FORTRAN

This code allows the generation of mechanical metamaterials based on a process of topological optimisation and numerical homogenisation.

This is essentially based on the integration of the work of authors such as Liu and Tovar(2018) and Andreassen et al.(2011) who efficiently developed topological optimisation processes and the work of Andreassen et al.(2013), Bendøe et al.(1988) and Yvonnet(2019) focused on the generation of cellular structures and characterisation by means of computational homogenisation.

The code was developed in FORTRAN and requires the Lapack, BLAS, OpenMP and Metis libraries for its operation. Additionally, we used an adapted version of the MMA algorithm developed by Svanberg (1987) and the MA87 of the HSL, which is a library that uses a direct method to solve large symmetric linear systems of AX=B equations.

This code was used in the article "How does the initial cell configuration influence the final topology in a metamaterial generation process?" published in the Latin American Journal of Solids and Structures.


Note: it is clarified that the code is free, its use is allowed for academic purposes and no warranty or license is given for commercial purposes.


- Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B., & Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization, 43(1), 1-16.
- Andreassen, E., & Andreasen, C. (2013). How to determine composite material properties using numerical homogenization. Computational Materials Science, 83, 488-495.
- Bendsøe, M.P., & Kikuchi, N. (1988). Generating optimal topologies in structural design using a homogenization method. Computer Methods in Applied Mechanics and Engineering, 71(2), 197-224.
- Liu, K., & Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. Structural and Multidisciplinary Optimization, 50(6), 1174-1196.
- Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. International journal for numerical methods in engineering, 24(2), 359-373.
- Yvonnet, J. (2019). Computational Homogenization of Heterogeneous Materials with Finite Elements. Springer.
