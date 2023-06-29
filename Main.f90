program Main
    ! Modules
    use TopologyOptimizationModule
    use ParaviewModule
    implicit none
    ! Derived variables
    type(Optimization)                                        :: OptimizationModel
    type(PostprocessingInfo)                                  :: Postprocessing
    !---------------------------------------------------------------------------------!
    !                   UNIVERSIDAD INDUSTRIAL DE SANTANDER - UIS                     !
    !                     PHD PROGRAM IN MECHANICAL ENGINEERING                       !
    !                              SCHOOL OF MECHANICAL                               !
    !---------------------------------------------------------------------------------!
    ! These routines were developed for academic, pedagogical and research purposes   !
    ! only. No guarantee is given for their results, nor are they authorised for      !
    ! commercial use.                                                                 !
    !---------------------------------------------------------------------------------!
    ! 1.1 Enter structure information
    call SetDimAnalysis(OptimizationModel,'2D')                 ! Set 2D or 3D
    call SetAnalysisType(OptimizationModel,'PlaneStress')       ! Set PlaneStress(2D), PlainStrain(2D), HookeLaw(3D)
    call SetElementType(OptimizationModel,'quad4')              ! tria3 tria6, quad4, quad8, tetra4, tetra10, hexa8, hexa20
    call SetThickness(OptimizationModel,50.0d0)                 ! Only 2D
    call SetYoungModulus(OptimizationModel,210000.0d0)          ! Young modulus
    call SetPoissonModulus(OptimizationModel,0.3d0)             ! Poisson modulus
    call SetGaussAprox(OptimizationModel,3)                     ! can use up to 5 gauss points
    call SetStrainNumericalH(OptimizationModel,0.1d0)           ! Strain used for numerical homogenization
    call SetCoordinates(OptimizationModel,'DataStructure/Coordinates.txt')
    call SetConnectivity(OptimizationModel,'DataStructure/Connectivity.txt')
    call SetCornerNodes(OptimizationModel,'DataStructure/CornerNodes.txt')
    call SetEdgeNodes(OptimizationModel,'DataStructure/EdgeNodes.txt')
    call SetFaceNodes(OptimizationModel,'DataStructure/FaceNodes.txt')
    call GetMasterSlaveNodes(OptimizationModel)
    call SetBondaryConditions(OptimizationModel,'DataStructure/BoundaryConditions.txt')
    ! 1.2 Enter Optimization parameters
    call SetVolFraction(OptimizationModel,0.5d0)                ! limiting volume fraction (lower limit of optimization)
    call SetMutationRate(OptimizationModel,0.1d0)               ! rate of change/mutation of optimization
    call SetFilterRadius(OptimizationModel,2.0d0*1.5)           ! smoothing filter radius (1.5 times de FE size)
    call SetMaxIterations(OptimizationModel,100)                ! maximum number of iterations
    call SetPenalFactor(OptimizationModel,3.0d0)                ! SIMP method penalty factor
    ! 1.3 Plotting selection                                       (3.0 for 2D and 5.0 for 3D)
    !    Strain   -   StrainEqv
    !    Stress   - StressVonMises
    ! DensEnergy  -    TopOpt
    call SetFilterValue(PostProcessing,0.2d0)                   ! doenst plot values below the filter
    call PlottingSelection(Postprocessing,'TopOpt')             ! Select only ONE
    call SetPostprocessingPath(Postprocessing,'ParaviewPostprocessing')
    ! 1.4 Star optimization process
    call TopologyOptimization(OptimizationModel)
    ! 1.5 Paraview Postprocessing (Select)
    call GenerateParaviewFiles(Postprocessing,OptimizationModel)
end program Main