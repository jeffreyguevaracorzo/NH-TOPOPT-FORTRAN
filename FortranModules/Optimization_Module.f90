module Optimization_Module
    use MMA_Module
    use HN_Module
    implicit none
    ! ----------------------------------------------------------------- !
    !     this variable groups all the information of the structure     !
    ! ----------------------------------------------------------------- !
    type, extends(Structure)                                    :: Optimization
        ! control variables
        integer                                                 :: Iteration
        integer                                                 :: MaxIterations
        double precision                                        :: PostProcesingFilter
        double precision                                        :: FilterRadius
        double precision                                        :: PenalFactor
        double precision                                        :: VolFraction
        double precision                                        :: MutationRate
        double precision                                        :: Change
        double precision, dimension(:), allocatable             :: DensityVector    
        double precision, dimension(:), allocatable             :: FinalDensityVector
        double precision, dimension(:), allocatable             :: Compilance
        double precision, dimension(:), allocatable             :: DiffVolume
        double precision, dimension(:), allocatable             :: DiffCompilance
        ! Additional for OC
        double precision                                        :: L1
        double precision                                        :: L2
        ! Additional for MMA
        integer                                                 :: mm
        integer                                                 :: nn
        double precision, dimension(:,:), allocatable           :: xmin
        double precision, dimension(:,:), allocatable           :: xmax
        double precision, dimension(:,:), allocatable           :: XMMA
        double precision, dimension(:,:), allocatable           :: xold1
        double precision, dimension(:,:), allocatable           :: xold2
        double precision, dimension(:,:), allocatable           :: low
        double precision, dimension(:,:), allocatable           :: upp
        double precision                                        :: f0val
        double precision, dimension(:,:), allocatable           :: df0dx
        double precision, dimension(:,:), allocatable           :: fval
        double precision, dimension(:,:), allocatable           :: dfdx
        double precision                                        :: a0 
        double precision, dimension(:,:), allocatable           :: a 
        double precision, dimension(:,:), allocatable           :: c 
        double precision, dimension(:,:), allocatable           :: d 
    contains
        procedure                                               :: MMAParameters
        procedure                                               :: SetMaxIterations
        procedure                                               :: SetFilterRadius
        procedure                                               :: SetPenalFactor
        procedure                                               :: SetVolFraction
        procedure                                               :: SetMutationRate
        procedure                                               :: UploadOptimizationParameters
        procedure                                               :: TopologyOptimizationProcess
    end type 

contains
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for TOP      !
    ! ----------------------------------------------------------------- !
    ! Additional routines for MMA Algorithm
    ! 1. setting MMA parameters
    subroutine MMAParameters(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision                                             :: xminval
        double precision                                             :: xmaxval
        Self%mm = 1
        Self%nn = Self%Ne
        allocate(Self%xmin(Self%nn,1))
        allocate(Self%xmax(Self%nn,1))
        allocate(Self%XMMA(Self%nn,1))
        allocate(Self%xold1(Self%nn,1))
        allocate(Self%xold2(Self%nn,1))
        allocate(Self%low(Self%nn,1))
        allocate(Self%upp(Self%nn,1))
        allocate(Self%fval(Self%mm,1))
        allocate(Self%dfdx(Self%mm,Self%nn))
        allocate(Self%df0dx(Self%nn,1))
        allocate(Self%a(Self%mm,1),Self%c(Self%mm,1),Self%d(Self%mm,1))
        xmaxval = 1.0d0
        xminval = 0.001d0
        Self%xmax = xmaxval
        Self%xmin = xminval
        Self%xold1(:,1) = Self%DensityVector
        Self%xold2(:,1) = Self%DensityVector
        Self%upp = xmaxval
        Self%low = xmaxval
        Self%a0 = 1.0d0
        Self%a = 0.0
        Self%c = 1000.0
        Self%d = 0.0
    end subroutine MMAParameters
    ! topology optimization subroutines
    ! 1. Input max iterations
    subroutine SetMaxIterations(Self,MaxIterations)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        integer, intent(in)                                          :: MaxIterations
        Self%MaxIterations = MaxIterations
    end subroutine SetMaxIterations
    ! 2. Input filter radius
    subroutine SetFilterRadius(Self,FilterRadius)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: FilterRadius
        Self%FilterRadius = FilterRadius
    end subroutine SetFilterRadius
    ! 3. Input penal factor
    subroutine SetPenalFactor(Self,PenalFactor)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: PenalFactor
        Self%PenalFactor = PenalFactor
    end subroutine SetPenalFactor
    ! 4. Input max volumne fraction
    subroutine SetVolFraction(Self,VolFraction)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: VolFraction
        Self%VolFraction = VolFraction
    end subroutine SetVolFraction
    ! 5. Input mutation/changing rate
    subroutine SetMutationRate(Self,Mutation)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                 :: Mutation
        Self%MutationRate = Mutation
    end subroutine SetMutationRate
    ! 6. Updating Optimization parameters
    subroutine UploadOptimizationParameters(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! Optimality criteria parameters
        if (self%DimAnalysis.eq.2) then
            self%L1 = 0.0d0
            self%L2 = 100000.0d0
        elseif(self%DimAnalysis.eq.3) then
            self%L1 = 0.0d0
            self%L2 = 1000000000.0d0
        end if
        self%Change = 1.0d0
        Self%Iteration = 0
        ! Density vector, compilance and diffcompilance
        allocate(self%DensityVector(self%Ne));
        self%DensityVector = 1.0d0
        allocate(self%Compilance(self%Ne));
        self%Compilance = 0.0d0
        allocate(self%DiffCompilance(self%Ne));
        self%DiffCompilance = 0.0d0
        ! Volume of all elements
        call VolumeAllElement(Self)
        ! Assign the MMA Params.
        call MMAParameters(Self)
    end subroutine UploadOptimizationParameters
    ! ----------------------------------------------------------------- !
    !    subroutines that define the topology optimization procedure    !
    ! ----------------------------------------------------------------- !
    ! 1. Density filter
    subroutine DensityFilter(self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        ! internal variables
        integer                                                      :: i,j,ne
        double precision                                             :: suma
        double precision                                             :: fac
        double precision, dimension(:), allocatable                  :: DiffCompilanceNew
        double precision, dimension(:), allocatable                  :: DiffVolumeNew
        double precision, dimension(:), allocatable                  :: Radius
        double precision, dimension(:,:), allocatable                :: RadPromE,PosPromE
        allocate(DiffCompilanceNew(self%Ne));      DiffCompilanceNew = 0.0d0
        allocate(DiffVolumeNew(Self%Ne));              DiffVolumeNew = 0.0d0
        allocate(RadPromE(self%Ne,Self%DimAnalysis));       RadPromE = 0.0d0
        allocate(PosPromE(self%Ne,Self%DimAnalysis));       PosPromE = 0.0d0
        ! Get elemento position
        do i = 1, Self%Ne, 1
            PosPromE(i,:) = Sum(Self%Coordinates(Self%ConnectivityN(i,:),:),1)/Self%Npe
        end do
        ! Star Filtering
        do i = 1, Self%Ne, 1
            suma = 0.0d0
            do j = 1, Self%DimAnalysis, 1
                RadPromE(:,j) = PosPromE(:,j) - PosPromE(i,j)
            end do
            Radius = sqrt(sum(RadPromE**2,2))
            do j = 1, Self%Ne, 1
                fac = Self%FilterRadius - Radius(j)
                suma = suma + max(0.0d0,fac)
                DiffCompilanceNew(i) = DiffCompilanceNew(i) + (max(0.0d0,fac))*Self%DensityVector(j)*Self%DiffCompilance(j)
                DiffVolumeNew(i) = DiffVolumeNew(i) + (max(0.0d0,fac))*Self%DensityVector(j)*Self%VolumePerElement(j)
            end do
            DiffCompilanceNew(i) = DiffCompilanceNew(i)/(Self%DensityVector(i)*suma)
            DiffVolumeNew(i) = DiffVolumeNew(i)/(Self%DensityVector(i)*suma)
        end do
        ! Diff compilance update
        Self%DiffCompilance = DiffCompilanceNew
        Self%DiffVolume = DiffVolumeNew/(Self%DensityVector*Self%VolumePerElement)
        deallocate(DiffCompilanceNew,DiffVolumeNew,Radius,RadPromE)
    end subroutine DensityFilter
    ! 2. Optimization algorithm
    ! 2.1. Optimality criteria (OC)
    subroutine OptimalityCriteria(self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        ! internal variables
        double precision                                             :: move
        double precision                                             :: xminval
        double precision                                             :: v1
        double precision                                             :: v2
        double precision                                             :: L1
        double precision                                             :: L2
        double precision                                             :: Lmid
        double precision, dimension(:), allocatable                  :: one
        double precision, dimension(:), allocatable                  :: DensityVectorNew
        double precision, dimension(:), allocatable                  :: DensityVectorOld
        write(unit=*, fmt=*) '** Using the optimality criteria method (OC)'
        allocate(one(Self%ne));                               one = 1.0d0
        allocate(DensityVectorNew(Self%ne));     DensityVectorNew = 0.0d0
        DensityVectorOld = Self%DensityVector
        xminval = 0.01d0
        move = Self%MutationRate
        L1 = Self%L1
        L2 = Self%L2
        do while((L2-L1)/(L2+L1).gt.1e-3)
            Lmid = 0.5d0*(L2 + L1)
            DensityVectorNew = min(DensityVectorOld+move,DensityVectorOld*sqrt((-1.0d0*Self%DiffCompilance)/Lmid))
            DensityVectorNew = min(one,DensityVectorNew)
            DensityVectorNew = max(DensityVectorOld-move,DensityVectorNew)
            DensityVectorNew = max(xminval,DensityVectorNew)
            v1 = sum(DensityVectorNew*Self%VolumePerElement)
            v2 = sum(Self%VolumePerElement*Self%VolFraction)
            if ((v1-v2).gt.0.0d0) then
                L1 = Lmid
            else
                L2 = Lmid
            end if
        end do
        ! updating
        Self%DensityVector = DensityVectorNew
        Self%Change = maxval(abs(DensityVectorOld-DensityVectorNew))
        deallocate(DensityVectorNew,DensityVectorOld,one)
    end subroutine OptimalityCriteria
    ! 2.2. Method of Moving Asymptotes (MMA)
    subroutine MethodMovingAsympotes(Self)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision                                             :: v1
        double precision                                             :: v2
        write(unit=*, fmt=*) '** using then Method of Moving Asympotes (MMA)'
        Self%xMMA(:,1) = Self%DensityVector
        v1 = sum(Self%DensityVector*Self%VolumePerElement)
        v2 = sum(Self%VolumePerElement*Self%VolFraction)
        Self%fval(1,1) = (v1 - v2)
        Self%f0val = sum(Self%Compilance)

        Self%dfdx(1,:) = Self%DiffVolume
        Self%df0dx(:,1) = Self%DiffCompilance
        ! applying the MMA solver interface of Svanberg
        call MMA_Solver_Interface(self%mm,Self%nn,Self%Iteration,Self%XMMA,Self%xmin,Self%xmax,Self%xold1,Self%xold2, &
                                  Self%f0val,Self%df0dx,Self%fval,Self%dfdx,Self%low,Self%upp,Self%a0,Self%a,Self%c,Self%d)
        ! Replace solutions
        Self%xold2 = Self%xold1
        Self%xold1(:,1) = Self%DensityVector
        Self%DensityVector = Self%XMMA(:,1)
        Self%Change = maxval(abs(Self%DensityVector-Self%xold1(:,1)))
    end subroutine MethodMovingAsympotes
    ! 3. Compilance and Diff Compilance (here goes the objetive function)
    subroutine GetCompilance(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! internal variables
        integer                                                      :: i
        double precision, dimension(:), allocatable                  :: ObjFun
        double precision, dimension(:), allocatable                  :: Emod,Vmod,Kmod,Gmod
        double precision, dimension(:,:,:), allocatable              :: Diso
        double precision, dimension(:,:,:), allocatable              :: penalization
        ! ------------------------------------------------------------------ !
        ! note: It is only necessary to change the formula that defines the  !
        !       objective function.                                          !
        ! Tensor variable convention (el, Comp-X, Comp-Y)                    !
        ! ------------------------------------------------------------------ !
        ! ** Isotropy forced penalization
        if (Self%DimAnalysis.eq.2) then
            allocate(Diso(Self%Ne,3,3)); Diso = 0.0d0
            allocate(penalization(Self%Ne,3,3)); penalization = 0.0d0
            Diso(:,1,1) = (self%DTensor(:,1,1)+self%DTensor(:,2,2))/2.0d0; Diso(:,2,2) = Diso(:,1,1)
            Diso(:,1,2) = (self%DTensor(:,1,2)+self%DTensor(:,2,1))/2.0d0; Diso(:,2,1) = Diso(:,1,2)
        elseif (Self%DimAnalysis.eq.3) then
            allocate(Diso(Self%Ne,6,6)); Diso = 0.0d0
            allocate(penalization(Self%Ne,6,6)); penalization = 0.0d0
            Diso(:,1,1) = (self%DTensor(:,1,1)+self%DTensor(:,2,2)+self%DTensor(:,3,3))/3.0d0
            Diso(:,2,2) = Diso(:,1,1); Diso(:,3,3) = Diso(:,1,1)
            Diso(:,4,4) = (self%DTensor(:,4,4)+self%DTensor(:,5,5)+self%DTensor(:,6,6))/3.0d0
            Diso(:,5,5) = Diso(:,4,4); Diso(:,6,6) = Diso(:,4,4)
            Diso(:,1,2) = (self%DTensor(:,1,2)+self%DTensor(:,1,3)+self%DTensor(:,2,3))/3.0d0
            Diso(:,1,3) = Diso(:,1,2); Diso(:,2,3) = Diso(:,1,2)
            Diso(:,2,1) = Diso(:,1,2); Diso(:,3,1) = Diso(:,1,2); Diso(:,3,2) = Diso(:,1,2)
        end if
        ! --- comment the penalization to avoid this ---
        Penalization = abs(Self%DTensor - Diso)
        Self%DTensor = Self%DTensor + Penalization
        ! ** Mechanical properties
        ! Emod - Young modulus
        ! Vmod - Poisson Modulus
        ! Kmod - Bulk Modulus
        ! Gmod - Shear Modulus
        if (Self%DimAnalysis.eq.2) then
            Emod = (self%DTensor(:,1,1)+self%DTensor(:,2,2))/2.0d0
            Vmod = (self%DTensor(:,1,2)+self%DTensor(:,2,1))/(2*Emod)
            Kmod = Emod/(3.0d0*(1.0d0-Vmod))
            Gmod = Emod/(2.0d0*(1.0d0+Vmod))
        elseif (Self%DimAnalysis.eq.3) then
            Emod = (self%DTensor(:,1,1)+self%DTensor(:,2,2)+self%DTensor(:,3,3))/3.0d0
            Vmod = (self%DTensor(:,1,2)+self%DTensor(:,2,3)+self%DTensor(:,1,3) &
                    +self%DTensor(:,2,1)+self%DTensor(:,3,2)+self%DTensor(:,3,1))/(2*Emod)
            Kmod = Emod/(3.0d0*(1.0d0-Vmod))
            Gmod = Emod/(2.0d0*(1.0d0+Vmod))
        end if
        ! -- some objetive functions --
        !ObjFun = Kmod                      ! Max compresion hid
        !ObjFun = Kmod/Gmod                 ! Max comp and min shear
        !ObjFun = Kmod/Vmod                 ! Max comp and min poisson --unstable--
        !ObjFun = vmod                      ! min poisson ratio
        ObjFun = Kmod
        self%Compilance = (self%DensityVector**self%PenalFactor)*ObjFun
        self%DiffCompilance = -self%PenalFactor*(self%DensityVector**(self%PenalFactor-1.0d0))*ObjFun
        deallocate(ObjFun,Diso,penalization)
    end subroutine GetCompilance
    ! 4. Printing information
    subroutine PrintStatus(Self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        write(unit=*, fmt=*) '--- IteNo.', Self%Iteration,'ObjFun.', Sum(self%Compilance),'Change', self%Change, &
                                'FracVol.', Sum(self%VolumePerElement*self%DensityVector)/Sum(self%VolumePerElement)
    end subroutine PrintStatus
    ! 5. Topology Optimization process
    subroutine TopologyOptimizationProcess(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        call UploadOptimizationParameters(self)
        write(unit=*, fmt=*) '---- topology optimization process ----'
        do while ((Self%Change.gt.0.001).and.(Self%Iteration.lt.Self%MaxIterations))
            Self%Iteration = Self%Iteration + 1
            write(unit=*, fmt=*) '-- Solve system with FEM'
            call UploadStructure(self,Self%DensityVector,Self%PenalFactor)
            call HSL87Solver(self)
            call NumericalHomogenization(self,Self%DensityVector,Self%PenalFactor)
            write(unit=*, fmt=*) '-- Getting compilance'
            call GetCompilance(self)
            write(unit=*, fmt=*) '-- Filtering'
            call DensityFilter(self)
            write(unit=*, fmt=*) '-- Getting new solution'
            ! choose the Opt. Algorhtm
            !call OptimalityCriteria(Self)
            call MethodMovingAsympotes(Self)
            ! Print Optimization status
            call PrintStatus(Self)
        end do
    end subroutine TopologyOptimizationProcess
end module Optimization_Module