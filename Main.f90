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
    write(unit=*, fmt=*) '1. Structure Properties'
    call SetAnalysisType(OptimizationModel,'PlaneStress')       ! Set PlaneStress(2D), PlaneStrain(2D), HookeLaw(3D)
    call SetElementType(OptimizationModel,'quad4')              ! tria3 tria6, quad4, quad8, tetra4, tetra10, hexa8, hexa20
    call SetThickness(OptimizationModel,50.0d0)                 ! Only 2D
    call SetYoungModulus(OptimizationModel,210000.0d0)          ! Young modulus
    call SetPoissonModulus(OptimizationModel,0.3d0)             ! Poisson modulus
    call SetGaussAprox(OptimizationModel,2)                     ! can use up to 5 gauss points
    call SetStrainNumericalH(OptimizationModel,0.1d0)           ! Strain used for numerical homogenization
    call ReadFiles(OptimizationModel)
    call PreAssemblyRoutine(OptimizationModel)
    ! 1.2 Enter Optimization parameters
    call SetVolFraction(OptimizationModel,0.35d0)              ! limiting volume fraction (lower limit of optimization)
    call SetMutationRate(OptimizationModel,0.1d0)             ! rate of change/mutation of optimization
    call SetFilterRadius(OptimizationModel,2.0d0*1.5)         ! smoothing filter radius (1.5 times de FE size)
    call SetMaxIterations(OptimizationModel,150)               ! maximum number of iterations
    call SetPenalFactor(OptimizationModel,3.0d0)              ! SIMP method penalty factor
    ! 1.3 Plotting selection
    !    Strain   -   StrainEqv
    !    Stress   - StressVonMises
    ! DensEnergy  -    TopOpt
    call SetFilterValue(PostProcessing,0.2d0)                   ! doenst plot values below the filter
    call PlottingSelection(Postprocessing,'TopOpt')             ! Select only ONE
    call SetPostprocessingPath(Postprocessing,'Paraview')
    ! 1.4 Star optimization process
    call TopologyOptimizationProcess(OptimizationModel)
    call ProcessingResults(OptimizationModel)
    ! 1.5 Paraview Postprocessing (Select)
    call ParaviewPostprocessing(Postprocessing,OptimizationModel)
end program Main
