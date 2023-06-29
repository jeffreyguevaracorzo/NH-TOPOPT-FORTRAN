module TopologyOptimizationModule
    ! this module needs the FEAModule
    use HomogenizationModule
    implicit none
    ! ----------------------------------------------------------------- !
    !     this variable groups all the information of the structure     !
    ! ----------------------------------------------------------------- !
    type, extends(Structure)                                    :: Optimization
        integer                                                 :: MaxIterations
        double precision                                        :: FilterRadius
        double precision                                        :: PenalFactor
        double precision                                        :: VolFraction
        double precision                                        :: MutationRate
        double precision                                        :: L1
        double precision                                        :: L2
        double precision                                        :: Change
        double precision, dimension(:), allocatable             :: DensityVector
        double precision, dimension(:), allocatable             :: Compilance
        double precision, dimension(:), allocatable             :: DiffCompilance
        double precision, dimension(:,:), allocatable           :: PosPromE
    contains
        procedure                                               :: SetMaxIterations
        procedure                                               :: SetFilterRadius
        procedure                                               :: SetPenalFactor
        procedure                                               :: SetVolFraction
        procedure                                               :: SetMutationRate
        procedure                                               :: UploadOptParameters
        procedure                                               :: TopologyOptimization
    end type 

contains
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for TOP      !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine SetMaxIterations(Self,MaxIterations)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        integer, intent(in)                                          :: MaxIterations
        Self%MaxIterations = MaxIterations
    end subroutine SetMaxIterations
    ! listo
    subroutine SetFilterRadius(Self,FilterRadius)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: FilterRadius
        Self%FilterRadius = FilterRadius
    end subroutine SetFilterRadius
    ! listo
    subroutine SetPenalFactor(Self,PenalFactor)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: PenalFactor
        Self%PenalFactor = PenalFactor
    end subroutine SetPenalFactor
    ! listo
    subroutine SetVolFraction(Self,VolFraction)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: VolFraction
        Self%VolFraction = VolFraction
    end subroutine SetVolFraction
    ! listo
    subroutine SetMutationRate(Self,Mutation)
        implicit none
        class(Optimization), intent(inout)                           :: Self
        double precision, intent(in)                                             :: Mutation
        Self%MutationRate = Mutation
    end subroutine SetMutationRate
    ! listo
    subroutine UploadOptParameters(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! Optimality criteria parameters
        if (self%DimAnalysis.eq.'2D') then
            self%L1 = 0.0d0
            self%L2 = 100000.0d0
        elseif(self%DimAnalysis.eq.'3D') then
            self%L1 = 0.0d0
            self%L2 = 1000000000.0d0
        end if
        self%Change = 1.0d0
        ! Density vector, compilance and diffcompilance
        allocate(self%DensityVector(self%Ne));
        self%DensityVector = Self%VolFraction
        allocate(self%Compilance(self%Ne));
        self%Compilance = 0.0d0
        allocate(self%DiffCompilance(self%Ne));
        self%DiffCompilance = 0.0d0
        ! Volume of all elements
        call VolumeAllElement(self)
    end subroutine UploadOptParameters
    ! ----------------------------------------------------------------- !
    !    subroutines that define the topology optimization procedure    !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine DensityFilter(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! internal variables
        integer                                                      :: i,j
        double precision                                                         :: suma
        double precision                                                         :: fac
        double precision, dimension(:), allocatable                              :: DiffCompilanceNew
        double precision, dimension(:), allocatable                              :: Radius
        double precision, dimension(:,:), allocatable                            :: RadPromE
        ! process
        allocate(DiffCompilanceNew(self%Ne));      DiffCompilanceNew = 0.0d0
        allocate(RadPromE(self%Ne,3));                      RadPromE = 0.0d0
        ! Position of each element
        if (allocated(self%PosPromE)) then;            self%PosPromE = 0.0d0
        else; allocate(self%PosPromE(self%Ne,3));      self%PosPromE = 0.0d0
        end if
        ! Get elemento position
        do i = 1, self%Ne, 1
            self%PosPromE(i,:) = Sum(self%Coordinates(self%ConnectivityN(i,:),:),1)/self%Npe
        end do
        ! Star Filtering
        do i = 1, self%Ne, 1
            suma = 0.0d0
            RadPromE(:,1) = self%PosPromE(:,1) - self%PosPromE(i,1) 
            RadPromE(:,2) = self%PosPromE(:,2) - self%PosPromE(i,2)
            RadPromE(:,3) = self%PosPromE(:,3) - self%PosPromE(i,3)
            Radius = sqrt(sum(RadPromE**2,2))
            do j = 1, self%Ne, 1
                fac = self%FilterRadius - Radius(j)
                suma = suma + max(0.0d0,fac)
                DiffCompilanceNew(i) = DiffCompilanceNew(i) + (max(0.0d0,fac))*self%DensityVector(j)*self%DiffCompilance(j)
            end do
            DiffCompilanceNew(i) = DiffCompilanceNew(i)/(self%DensityVector(i)*suma)
        end do
        ! Diff compilance update
        self%DiffCompilance = DiffCompilanceNew
        deallocate(DiffCompilanceNew,Radius,RadPromE)
    end subroutine DensityFilter
    ! listo
    subroutine OptimalityCriteria(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! internal variables
        double precision                                                         :: move
        double precision                                                         :: L1
        double precision                                                         :: L2
        double precision                                                         :: Lmid
        double precision, dimension(:), allocatable                              :: one
        double precision, dimension(:), allocatable                              :: DensityVectorNew
        double precision, dimension(:), allocatable                              :: DensityVectorOld
        ! process
        allocate(one(self%Ne));                           one = 1.0d0
        allocate(DensityVectorNew(self%Ne)); DensityVectorNew = 0.0d0
        allocate(DensityVectorOld(self%Ne)); DensityVectorOld = self%DensityVector
        move = self%MutationRate
        L1 = self%L1
        L2 = self%L2
        ! star optimality process
        do while((L2-L1)/(L2+L1).gt.1e-3)
            Lmid = 0.5d0*(L2+L1)
            DensityVectorNew = max(0.01d0*one,max(DensityVectorOld-move,min(one,min(DensityVectorOld+move, &
                            DensityVectorOld*sqrt((-1.0)*self%DiffCompilance/Lmid)))))
            ! Volume constrain
            if ((sum(DensityVectorNew*self%VolumePerElement)-sum(self%VolumePerElement*self%VolFraction)).gt.0.0d0) then
                L1 = Lmid
            else
                L2 = Lmid
            end if
        end do
        ! updating
        self%DensityVector = DensityVectorNew
        call WriteDoublePrecisionVector('DensityVector.txt',self%DensityVector)
        self%Change = maxval(abs(DensityVectorOld-DensityVectorNew))
        deallocate(DensityVectorNew,DensityVectorOld,one)
    end subroutine OptimalityCriteria
    ! listo
    subroutine GetCompilance(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! internal variables
        integer                                                      :: i
        double precision, dimension(:), allocatable                              :: ObjFun
        ! process
        ! select or define the objective function
        ! Tensor (el, dimx, dimy)
        ObjFun = self%DTensor(:,1,1)+self%DTensor(:,2,2) 
        ! getting the compilance and the diffcompilance
        self%Compilance = (self%DensityVector**self%PenalFactor)*ObjFun
        self%DiffCompilance = -self%PenalFactor*(self%DensityVector**(self%PenalFactor-1.0d0))*ObjFun
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/ObjFunction.txt',ObjFun)
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/Compilance.txt',Self%Compilance)
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/DiffCompilance.txt',Self%DiffCompilance)
        deallocate(ObjFun)
    end subroutine GetCompilance
    ! listo
    subroutine TopologyOptimization(self)
        implicit none
        class(Optimization), intent(inout)                           :: self
        ! internal variables
        integer                                                      :: ite = 0
        ! inicializar parametros
        call UploadOptParameters(self)
        ! inicia proceso de optimizacion
        write(unit=*, fmt=*) '---- topology optimization process ----'
        do while (self%Change.gt.0.001)
            ite = ite + 1
            ! Actualizar estructura
            write(unit=*, fmt=*) 'Solver system with FEM'
            call UploadStructure(self,Self%DensityVector,Self%PenalFactor)
            ! Aplicar HSL87
            call HSL87Solver(self)
            call NumericalHomogenization(self,Self%DensityVector,Self%PenalFactor)
            call ProcessingResults(self)
            ! Calcular compilance y diff compilance
            write(unit=*, fmt=*) 'Getting compilance'
            call GetCompilance(self)
            ! Filtro de densidad
            write(unit=*, fmt=*) 'Filtering'
            call DensityFilter(self)
            ! Criterio de optimalidad
            write(unit=*, fmt=*) 'Getting new solution'
            call OptimalityCriteria(self)
            write(unit=*, fmt=*) 'IteNo.', ite, 'ObjFun.',Sum(self%Compilance), &
                                 'FracVol.',Sum(self%VolumePerElement*self%DensityVector)/Sum(self%VolumePerElement), &
                                 'Change',self%Change
            if (ite.ge.self%MaxIterations) then; exit; end if
        end do
    end subroutine TopologyOptimization
end module TopologyOptimizationModule