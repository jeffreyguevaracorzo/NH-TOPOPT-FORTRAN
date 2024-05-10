module HN_Module
    use Base_Module
    use Solver_Module
    implicit none
    ! ----------------------------------------------------------------- !
    !     this variable groups all the information of the structure     !
    ! ----------------------------------------------------------------- !
    type                                                        :: Structure        ! Name of the type of derived structure
        character(len=20)                                       :: AnalysisType     ! PlaneStress(2D), PlainStrain(2D), HookeLaw(3D)
        character(len=20)                                       :: ElementType      ! tri3 tri6, cuad4, cuad8, tetra4, tetra10, hexa8, hexa20
        integer                                                 :: DimAnalysis      ! 2D or 3D
        integer                                                 :: QuadGauss        !
        integer                                                 :: N                ! Number of nodes
        integer                                                 :: Ne               ! Number of elements
        integer                                                 :: Npe              ! Number of nodes per element
        integer                                                 :: NBc              ! Number of nodes with restrictions
        integer                                                 :: MN_Corner        ! Number of master corner nodes
        integer                                                 :: MN_Edge          ! Number of master edge nodes
        integer                                                 :: MN_Face          ! Number of master face nodes
        integer, Dimension(:), allocatable                      :: FreeD            ! Free degrees of freedom
        integer, Dimension(:), allocatable                      :: FixedD           ! Fixed degrees of freedom
        integer, dimension(:,:), allocatable                    :: CornerNodes      ! Master and slave corner nodes
        integer, dimension(:,:), allocatable                    :: EdgeNodes        ! Master and slave edge nodes
        integer, dimension(:,:), allocatable                    :: FaceNodes        ! Master and slave face nodes
        integer, dimension(:,:), allocatable                    :: ConnectivityN    ! Node connectivity per element
        integer, dimension(:,:), allocatable                    :: ConnectivityD    ! Dof Connectivity per element
        integer, dimension(:,:), allocatable                    :: ConnectivityNPBC ! Node connectivity per element (with PBC)
        integer, dimension(:,:), allocatable                    :: ConnectivityDPBC ! Dof Connectivity per element (with PBC)
        integer, dimension(:,:), allocatable                    :: BoundaryC        ! Boundary Conditions
        double precision                                        :: YoungModulus     ! 
        double precision                                        :: PoissonModulus   !
        double precision                                        :: Thickness        ! Only for 2D analysis (in 3D thickness = 0.0)
        double precision                                        :: StrainH          ! Strain used for the numerical homogenization
        double precision, dimension(:), allocatable             :: GaussPoint       ! Gauss Approximation points
        double precision, dimension(:), allocatable             :: GaussWeights     ! Gauss Approximation weights
        double precision, dimension(:), allocatable             :: VolumePerElement ! Volume
        double precision, dimension(:,:), allocatable           :: Coordinates      ! Coordinates of nodes
        double precision, dimension(:,:,:), allocatable         :: BLocal           ! Matrix of derivatives of shape functions
        double precision, dimension(:,:,:), allocatable         :: KLocal           ! Local stiffness matrix per element
        double precision, dimension(:,:,:), allocatable         :: FLocal           ! All the 3 or 6 force states
        double precision, dimension(:,:,:), allocatable         :: DTensor          ! Elastic Tensor per element
        ! System of equations [K]{u}={f} in sparse format               
        integer, dimension(:), allocatable                      :: index_i_KGlobal  ! Indices for rows of the global stiffness matrix
        integer, dimension(:), allocatable                      :: index_j_KGlobal  ! Indices for columns of the global stiffness matrix
        integer, dimension(:), allocatable                      :: Rows_KGlobal     !
        integer, dimension(:), allocatable                      :: Cols_KGlobal     !
        integer, dimension(:,:,:), allocatable                  :: Location_KGlobal !
        double precision, dimension(:), allocatable             :: value_KGlobal    ! Value of the global stiffness matrix at position i,j
        double precision, dimension(:,:), allocatable           :: value_FGlobal    ! Global load vector (for Sparse system)
        double precision, dimension(:,:), allocatable           :: value_UGlobal    ! Global displacement vector (for Sparse system)
        ! Results
        double precision, dimension(:,:), allocatable           :: UGlobal          ! Global displacement vector
        double precision, dimension(:,:), allocatable           :: UGlobal_0        ! Ideal Global displacement 
        double precision, dimension(:,:), allocatable           :: StrainEnergyE    ! Strain energy per element
        double precision, dimension(:,:,:), allocatable         :: StrainE          ! Strain per element
        double precision, dimension(:,:,:), allocatable         :: StressE          ! Stress per element
    contains
        procedure                                               :: SetAnalysisType
        procedure                                               :: SetElementType
        procedure                                               :: SetYoungModulus
        procedure                                               :: SetPoissonModulus
        procedure                                               :: SetThickness
        procedure                                               :: SetGaussAprox
        procedure                                               :: SetStrainNumericalH
        procedure                                               :: SetCoordinates
        procedure                                               :: SetConnectivity
        procedure                                               :: SetCornerNodes
        procedure                                               :: SetEdgeNodes
        procedure                                               :: SetFaceNodes
        procedure                                               :: GetMasterSlaveNodes
        procedure                                               :: SetBondaryConditions
        procedure                                               :: UploadStructure
        procedure                                               :: HSL87Solver
        procedure                                               :: ProcessingResults
    end type
    
contains
    ! ----------------------------------------------------------------- !
    !           base subroutines for reading and operations             !
    ! ----------------------------------------------------------------- !
    ! listo
    Subroutine ReadingfileInteger(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                         :: Path
        integer                                              :: ios,iounit,i,j
        integer, intent(inout)                               :: Nrow
        integer, intent(in)                                  :: Ncol
        integer, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileInteger
    Subroutine ReadingfileDP(Path,Nrow,Ncol,Matrix)
        implicit none
        character(len=*), intent(in)                                  :: Path
        integer                                                       :: ios,iounit,i,j
        integer, intent(inout)                                        :: Nrow
        integer, intent(in)                                           :: Ncol
        double precision, dimension(:,:), allocatable, intent(inout)  :: Matrix
        open(unit=iounit, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=iounit, fmt=*) !title
            read(unit=iounit, fmt=*) Nrow ; allocate(Matrix(Nrow,Ncol))
            read(unit=iounit, fmt=*) !references
            do i = 1, Nrow, 1
                read(unit=iounit, fmt=*) (Matrix(i,j), j = 1, Ncol, 1)
            end do
        close(iounit)
    end subroutine ReadingfileDP
    ! 2. Area
    function Area(Coordinates,Type) result(AreaAprox)
        character(len=*), intent(in)                                :: Type
        double precision, dimension(:,:), allocatable, intent(in)   :: Coordinates
        ! internal variables
        double precision, dimension(:), allocatable                 :: vector1, vector2, vector3, vector4
        double precision, dimension(:), allocatable                 :: Av1, Av2
        double precision                                            :: AreaAprox
        allocate(vector1(3),vector2(3),vector3(3),vector4(3),Av1(3), Av2(3))
        if ((Type.eq.'tetra4').or.(Type.eq.'tetra10')) then
            vector1 = Coordinates(2,:)-Coordinates(1,:)
            vector2 = Coordinates(3,:)-Coordinates(1,:)
            Av1 = CrossProduct(vector1,vector2)
            AreaAprox = abs(Norm(Av1))/2
        elseif ((Type.eq.'hexa8').or.(Type.eq.'hexa20')) then
            vector1 = Coordinates(1,:)-Coordinates(2,:)
            vector2 = Coordinates(3,:)-Coordinates(2,:)
            vector3 = Coordinates(1,:)-Coordinates(4,:)
            vector4 = Coordinates(3,:)-Coordinates(4,:)
            Av1 = CrossProduct(vector1,vector2)
            Av2 = CrossProduct(vector3,vector4)
            AreaAprox = abs(Norm(Av1))/2 + abs(Norm(Av2))/2
        end if
    end function Area
    ! 3. Volume
    subroutine VolumeAllElement(self)
        implicit none
        class(Structure), intent(inout)                              :: self
        ! internal variables
        integer                                                      :: i
        double precision                                             :: Area
        double precision, dimension(:), allocatable                  :: v1,v2,v3,v4,Centroid
        double precision, dimension(:,:), allocatable                :: Coordinates
        ! process
        allocate(self%VolumePerElement(self%Ne))
        self%VolumePerElement = 0.0d0
        ! Process
        do i = 1, self%Ne, 1
            Area = 0.0d0
            if ((self%ElementType.eq.'tria3').or.(self%ElementType.eq.'tria6')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(3,:) - Coordinates(2,:),0.0d0]
                v2 = [Coordinates(1,:) - Coordinates(2,:),0.0d0]
                Area = norm2(CrossProduct(v1,v2))/2.0d0
                self%VolumePerElement(i) = Area*self%Thickness
            elseif ((self%ElementType.eq.'quad4').or.(self%ElementType.eq.'quad8')) then
                Coordinates = Self%Coordinates(Self%ConnectivityN(i,:),:)
                v1 = [Coordinates(2,:) - Coordinates(1,:),0.0d0]
                v2 = [Coordinates(4,:) - Coordinates(1,:),0.0d0]
                v3 = [Coordinates(4,:) - Coordinates(3,:),0.0d0]
                v4 = [Coordinates(2,:) - Coordinates(3,:),0.0d0]
                Area = norm2(CrossProduct(v1,v2))/2.0d0 + norm2(CrossProduct(v3,v4))/2.0d0
                self%VolumePerElement(i) = Area*self%Thickness
            elseif ((self%ElementType.eq.'tetra4').or.(self%ElementType.eq.'tetra10')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                v1 = Coordinates(3,:) - Coordinates(1,:)
                v2 = Coordinates(2,:) - Coordinates(1,:)
                v3 = Coordinates(4,:) - Coordinates(1,:)
                self%VolumePerElement(i) = abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            elseif ((self%ElementType.eq.'hexa8').or.(self%ElementType.eq.'hexa20')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                Centroid = sum(Coordinates,1)/size(Coordinates(:,1))
                ! 1
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(2,:) - Centroid
                v3 = Coordinates(5,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 2
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(5,:) - Centroid
                v3 = Coordinates(6,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 3
                v1 = Coordinates(3,:) - Centroid
                v2 = Coordinates(7,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 4
                v1 = Coordinates(3,:) - Centroid
                v2 = Coordinates(4,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 5
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(6,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 6
                v1 = Coordinates(2,:) - Centroid
                v2 = Coordinates(3,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 7
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(5,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 8
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(4,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 9
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(2,:) - Centroid
                v3 = Coordinates(3,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 10
                v1 = Coordinates(1,:) - Centroid
                v2 = Coordinates(3,:) - Centroid
                v3 = Coordinates(4,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 11
                v1 = Coordinates(5,:) - Centroid
                v2 = Coordinates(6,:) - Centroid
                v3 = Coordinates(7,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
                ! 12
                v1 = Coordinates(5,:) - Centroid
                v2 = Coordinates(7,:) - Centroid
                v3 = Coordinates(8,:) - Centroid
                self%VolumePerElement(i) = self%VolumePerElement(i) + abs(dot_product(CrossProduct(v1,v2),v3))/6.0d0
            end if
            deallocate(Coordinates)
        end do
    end subroutine VolumeAllElement
    ! 4. Sorting
    subroutine Sort(vector,n)
        implicit none
        integer, intent(in)                 :: n
        integer, allocatable, intent(inout) :: vector(:)
        ! internal variables
        integer                             :: i, j, temp
        integer, allocatable                :: dumb(:)
        ! Ordenado de valores
        do i = 1, n-1
            do j = 1, n-i
                if (vector(j) > vector(j+1)) then
                    temp = vector(j)
                    vector(j) = vector(j+1)
                    vector(j+1) = temp
                endif
            end do
        end do
        ! eliminar 0's y repetidos
        do i = 2, n, 1
            if (Vector(i-1).eq.vector(i)) then 
            vector(i-1) = 0
            end if
        end do
        vector = pack(vector,vector.gt.0)
    end subroutine Sort
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for FEA      !
    ! ----------------------------------------------------------------- !
    ! 1. Input Analysis type
    subroutine SetAnalysisType(Self,AnalysisType)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: AnalysisType
        Self%AnalysisType = AnalysisType
        if (Self%AnalysisType.eq.'PlaneStress'.or.Self%AnalysisType.eq.'PlaneStrain') then
            Self%DimAnalysis = 2
        elseif (Self%AnalysisType.eq.'HookeLaw') then 
            Self%DimAnalysis = 3
        else
            stop "ERROR, Setting Analysis type"
        end if
    end subroutine SetAnalysisType
    ! 2. Input Element type
    subroutine SetElementType(Self,ElementType)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: ElementType
        Self%ElementType = ElementType
        if ( Self%ElementType.eq.'tria3' ) then; Self%Npe = 3 ; end if
        if ( Self%ElementType.eq.'tria6' ) then; Self%Npe = 6 ; end if
        if ( Self%ElementType.eq.'quad4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'quad8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'tetra4' ) then; Self%Npe = 4 ; end if
        if ( Self%ElementType.eq.'tetra10' ) then; Self%Npe = 10 ; end if
        if ( Self%ElementType.eq.'hexa8' ) then; Self%Npe = 8 ; end if
        if ( Self%ElementType.eq.'hexa20' ) then; Self%Npe = 20 ; end if
        !write(unit=*, fmt=*) 'ElementType   ',Self%ElementType
    end subroutine SetElementType
    ! 3. Input Thickness (only for 3D cases)
    subroutine SetThickness(Self,Thickness)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: Thickness
        Self%Thickness = Thickness
    end subroutine SetThickness
    ! 4. Input Young modulus
    subroutine SetYoungModulus(Self,YoungModulus)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: YoungModulus
        Self%YoungModulus = YoungModulus
    end subroutine SetYoungModulus
    ! 5. Input Poisson modulus
    subroutine SetPoissonModulus(Self,PoissonModulus)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: PoissonModulus
        Self%PoissonModulus = PoissonModulus
    end subroutine SetPoissonModulus
    ! 6. Input Gauss approximation
    subroutine SetGaussAprox(Self,Gauss)
        implicit none
        class(Structure), intent(inout)                             :: Self 
        integer, intent(in)                                         :: Gauss
        double precision                                            :: LB, UB, Coef1, Coef2
        if (Gauss.eq.0.or.Gauss.gt.5) stop "ERROR, GuassQuadrature (greater than 5 or <= to zero)" 
        Self%QuadGauss = Gauss
        select case (Gauss)
            case (1)
                allocate(Self%GaussPoint(1));       Self%GaussPoint = [0.00000000]
                allocate(Self%GaussWeights(1));   Self%GaussWeights = [2.00000000]
            case (2)
                allocate(Self%GaussPoint(2));       Self%GaussPoint = [-0.57735026,0.57735026]
                allocate(Self%GaussWeights(2));   Self%GaussWeights = [1.00000000,1.00000000]
            case (3)
                allocate(Self%GaussPoint(3));       Self%GaussPoint = [-0.77459666,0.00000000,0.77459666]
                allocate(Self%GaussWeights(3));   Self%GaussWeights = [0.55555555,0.88888888,0.55555555]
            case (4)
                allocate(Self%GaussPoint(4));       Self%GaussPoint = [-0.86113631,-0.33998104,0.33998104,0.86113631]
                allocate(Self%GaussWeights(4));   Self%GaussWeights = [0.34785484,0.65214515,0.65214515,0.34785484]
            case (5)
                allocate(Self%GaussPoint(5));       Self%GaussPoint = [-0.90617984,-0.53846931,0.00000000,0.53846931,0.90617984]
                allocate(Self%GaussWeights(5));   Self%GaussWeights = [0.23692688,0.47862867,0.56888888,0.47862867,0.23692688]
        end select
        LB = 0.0d0;  UB = 1.0d0
        Coef1 = (UB - LB)/2.0d0
        Coef2 = (UB + LB)/2.0d0
        ! changing coordinates and points due to boundary changes
        if (Self%ElementType.eq.'tria3'.or.Self%ElementType.eq.'tria6'.or. &
            Self%ElementType.eq.'tetra4'.or.Self%ElementType.eq.'tetra10') then
            Self%GaussPoint = Self%GaussPoint*Coef1 + Coef2
            Self%GaussWeights = Self%GaussWeights*Coef1
        end if
    end subroutine SetGaussAprox
    ! 7. Input Strain for numerical Homogenization
    subroutine SetStrainNumericalH(Self,Strain)
        implicit none
        class(Structure), intent(inout)                                 :: Self 
        double precision, intent(in)                                    :: Strain
        Self%StrainH = Strain
        !write(unit=*, fmt=*) 'Strain for Homogenization   ',Self%StrainH
    end subroutine SetStrainNumericalH
    ! 8. Input Coordinates
    subroutine SetCoordinates(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        call ReadingfileDP(Path,Self%N,Self%DimAnalysis,Self%Coordinates)
        write(unit=*, fmt=*) '- Coordinates'
        !call FilePrinting(Self%Coordinates,'DataResults/.InternalData/CoordinatesLecture.txt')
    end subroutine SetCoordinates
    ! 9. Input Connectivity
    subroutine SetConnectivity(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        integer                                                     :: i,j,k
        call ReadingfileInteger(Path,Self%Ne,Self%Npe,Self%ConnectivityN)
        allocate(Self%ConnectivityD(Self%Ne,Self%Npe*Self%DimAnalysis))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = Self%DimAnalysis-1, 0, -1
                    Self%ConnectivityD(i,j*Self%DimAnalysis - k) = Self%ConnectivityN(i,j)*Self%DimAnalysis - k
                end do
            end do
        end do
        write(unit=*, fmt=*) '- Connectivity'
        !call FilePrinting(Self%ConnectivityN,'DataResults/.InternalData/ConnectivityNLecture.txt')
        !call FilePrinting(Self%ConnectivityD,'DataResults/.InternalData/ConnectivityDLecture.txt')
    end subroutine SetConnectivity
    ! 10. Input Cornernodes
    subroutine SetCornerNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: i,j,ios
        if (Self%DimAnalysis.eq.2) then; j=3; end if
        if (Self%DimAnalysis.eq.3) then; j=7; end if
        write(unit=*, fmt=*) '- Corner Nodes'
        open(unit=1, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error corner nodes"
            allocate(Self%CornerNodes(j,2)); Self%CornerNodes = 0
            read(unit=1, fmt=*) ! master nodes
            read(unit=1, fmt=*) Self%CornerNodes(1,1)
            read(unit=1, fmt=*) ! slave nodes
            read(unit=1, fmt=*) (Self%CornerNodes(i,2),i=1,j,1)
        close(1)
        !call FilePrinting(Self%CornerNodes,'DataResults/.InternalData/CornerNodes.txt')
    end subroutine SetCornerNodes
    ! 11. Input Edgenodes
    subroutine SetEdgeNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: i,j,ios
        write(unit=*, fmt=*) '- Edge Nodes'
        select case (Self%DimAnalysis)
            case (2)
                open(unit=1, file=Path, iostat=ios, status="old", action="read")
                    read(unit=1, fmt=*) !header
                    read(unit=1, fmt=*) j
                    allocate(Self%EdgeNodes(j*2,2)) ; Self%EdgeNodes = 0
                    read(unit=1, fmt=*) !master nodes
                    read(unit=1, fmt=*) !x
                    read(unit=1, fmt=*) (Self%EdgeNodes(i,1),i=1,j,1)
                    read(unit=1, fmt=*) !y
                    read(unit=1, fmt=*) (Self%EdgeNodes(j+i,1),i=1,j,1)
                    read(unit=1, fmt=*) !z(ignore)
                    read(unit=1, fmt=*) ! -
                    read(unit=1, fmt=*) !slave nodes
                    read(unit=1, fmt=*) !x
                    read(unit=1, fmt=*) (Self%EdgeNodes(i,2),i=1,j,1)
                    read(unit=1, fmt=*) !y
                    read(unit=1, fmt=*) (Self%EdgeNodes(j+i,2),i=1,j,1)
                    read(unit=1, fmt=*) !z(ignore)
                    read(unit=1, fmt=*) ! -
                close(1)
            case (3)
                open(unit=1, file=Path, iostat=ios, status="old", action="read")
                    read(unit=1, fmt=*) !header
                    read(unit=1, fmt=*) j
                    allocate(Self%EdgeNodes(j*9,2)) ; Self%EdgeNodes = 0
                    read(unit=1, fmt=*) !master nodes
                    read(unit=1, fmt=*) !x
                    read(unit=1, fmt=*) (Self%EdgeNodes(i,1),i=1,j,1)
                    read(unit=1, fmt=*) !y
                    read(unit=1, fmt=*) (Self%EdgeNodes(3*j+i,1),i=1,j,1)
                    read(unit=1, fmt=*) !z(ignore)
                    read(unit=1, fmt=*) (Self%EdgeNodes(6*j+i,1),i=1,j,1)
                    read(unit=1, fmt=*) !slave nodes
                    read(unit=1, fmt=*) !x
                    read(unit=1, fmt=*) (Self%EdgeNodes(i,2),i=1,3*j,1)
                    read(unit=1, fmt=*) !y
                    read(unit=1, fmt=*) (Self%EdgeNodes(3*j+i,2),i=1,3*j,1)
                    read(unit=1, fmt=*) !z(ignore)
                    read(unit=1, fmt=*) (Self%EdgeNodes(6*j+i,2),i=1,3*j,1)
                close(1)
        end select
        !call FilePrinting(Self%EdgeNodes,'DataResults/.InternalData/EdgeNodes.txt')
    end subroutine SetEdgeNodes
    ! 12. Input Facenodes
    subroutine SetFaceNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: dim,i,j,ios
        write(unit=*, fmt=*) '- Face Nodes'
        select case (Self%DimAnalysis)
            case (2)
                allocate(Self%FaceNodes(1,2))
                Self%FaceNodes = 0
            case (3)
                open(unit=1, file=Path, iostat=ios, status="old", action="read")
                    if ( ios /= 0 ) stop "Error opening file name"
                    read(unit=1, fmt=*) !header
                    read(unit=1, fmt=*) j
                    allocate(Self%FaceNodes(j*3,2)) ; Self%FaceNodes = 0
                    read(unit=1, fmt=*) !master node
                    read(unit=1, fmt=*) !xy
                    read(unit=1, fmt=*) (Self%FaceNodes(i,1),i=1,j,1)
                    read(unit=1, fmt=*) !xz
                    read(unit=1, fmt=*) (Self%FaceNodes(j+i,1),i=1,j,1)
                    read(unit=1, fmt=*) !yz
                    read(unit=1, fmt=*) (Self%FaceNodes(2*j+i,1),i=1,j,1)
                    read(unit=1, fmt=*) !slave nodes
                    read(unit=1, fmt=*) !xy
                    read(unit=1, fmt=*) (Self%FaceNodes(i,2),i=1,j,1)
                    read(unit=1, fmt=*) !xz
                    read(unit=1, fmt=*) (Self%FaceNodes(j+i,2),i=1,j,1)
                    read(unit=1, fmt=*) !yz
                    read(unit=1, fmt=*) (Self%FaceNodes(2*j+i,2),i=1,j,1)
                close(1)
        end select
        !call FilePrinting(Self%FaceNodes,'DataResults/.InternalData/FaceNodes.txt')
    end subroutine SetFaceNodes
    ! 13. Get Master-Slave nodes
    subroutine GetMasterSlaveNodes(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal values
        integer                                                         :: i,j,k,l,dim
        integer, dimension(:), allocatable                              :: Array
        integer, dimension(:,:), allocatable                            :: CondArray
        write(unit=*, fmt=*) '- Assembly Master-Slave Nodes'
        ! corner nodes
        Self%CornerNodes(:,1) = Self%CornerNodes(1,1)
        ! edge nodes
        select case (self%DimAnalysis)
            case (2)
                k = size(Self%EdgeNodes(:,1))/2
                allocate(Array(k))
                Array = Self%EdgeNodes(1:k,1)
                ! x axis
                do i = 1, k, 1 !master
                    do j = 1, k, 1 !slave
                        if (Self%Coordinates(Array(i),1).eq.Self%Coordinates(Self%EdgeNodes(j,2),1)) then
                            Self%EdgeNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                ! y axis
                Array = Self%EdgeNodes((k+1):2*k,1)
                do i = 1, k, 1 !master
                    do j = k+1, 2*k, 1 !slave
                        if (Self%Coordinates(Array(i),2).eq.Self%Coordinates(Self%EdgeNodes(j,2),2)) then
                            Self%EdgeNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                deallocate(Array)
            case (3)
                k = size(Self%EdgeNodes(:,1))/9
                allocate(Array(k))
                ! x axis
                Array = Self%EdgeNodes(1:k,1)
                do i = 1, k, 1 !master
                    do j = 1, 3*k, 1 !slave
                        if (Self%Coordinates(Array(i),1).eq.Self%Coordinates(Self%EdgeNodes(j,2),1)) then
                            Self%EdgeNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                ! y axis
                Array = Self%EdgeNodes((3*k+1):4*k,1)
                do i = 1, k, 1 !master
                    do j = 3*k+1, 6*k, 1 !slave
                        if (Self%Coordinates(Array(i),2).eq.Self%Coordinates(Self%EdgeNodes(j,2),2)) then
                            Self%EdgeNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                ! z axis
                Array = Self%EdgeNodes((6*k+1):7*k,1)
                do i = 1, k, 1 !master
                    do j = 6*k+1, 9*k, 1 !slave
                        if (Self%Coordinates(Array(i),3).eq.Self%Coordinates(Self%EdgeNodes(j,2),3)) then
                            Self%EdgeNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                deallocate(Array)
        end select
        ! face nodes
        select case (self%DimAnalysis)
            case (2)
                !ignore
                Self%FaceNodes = 0
            case (3)
                k = size(Self%FaceNodes(:,1))/3
                allocate(Array(k))
                ! xy plane
                Array = Self%FaceNodes(1:k,1)
                do i = 1, k, 1 !master
                    do j = 1, k, 1 !slave
                        if (all(Self%Coordinates(Array(i),[1,2]).eq.Self%Coordinates(Self%FaceNodes(j,2),[1,2]))) then
                            Self%FaceNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                ! xz plane
                Array = Self%FaceNodes(k+1:2*k,1)
                do i = 1, k, 1 !master
                    do j = k+1, 2*k, 1 !slave
                        if (all(Self%Coordinates(Array(i),[1,3]).eq.Self%Coordinates(Self%FaceNodes(j,2),[1,3]))) then
                            Self%FaceNodes(j,1) = Array(i)
                        end if
                    end do
                end do
                ! yz plane
                Array = Self%FaceNodes(2*k+1:3*k,1)
                do i = 1, k, 1 !master
                    do j = 2*k+1, 3*k, 1 !slave
                        if (all(Self%Coordinates(Array(i),[2,3]).eq.Self%Coordinates(Self%FaceNodes(j,2),[2,3]))) then
                            Self%FaceNodes(j,1) = Array(i)
                        end if
                    end do
                end do
        end select
        ! node reduction
        ! reduce edge-corner
        do i = 1, size(Self%CornerNodes(:,1)), 1
            do j = 1, size(Self%EdgeNodes(:,1)), 1
                if (Self%CornerNodes(i,2).eq.Self%EdgeNodes(j,1)) then
                    Self%EdgeNodes(j,1) = Self%CornerNodes(i,1)
                end if
            end do
        end do
        ! reduce face-edge
        do i = 1, size(Self%EdgeNodes(:,1)), 1
            do j = 1, size(Self%FaceNodes(:,1)), 1
                if (Self%EdgeNodes(i,2).eq.Self%FaceNodes(j,1)) then
                    Self%FaceNodes(j,1) = Self%EdgeNodes(i,1)
                end if
            end do
        end do
        ! assembly of connectivities with periodic boundary conditions
        i = size(Self%CornerNodes(:,1))
        j = size(Self%EdgeNodes(:,1))
        k = Size(Self%FaceNodes(:,1))
        allocate(CondArray(i+j+k,2))
        CondArray(1:i,:) = Self%CornerNodes
        CondArray((i+1):(i+j),:) = Self%EdgeNodes
        CondArray((i+j+1):(i+j+k),:) = Self%FaceNodes
        allocate(Self%ConnectivityNPBC(Self%Ne,Self%Npe))
        allocate(Self%ConnectivityDPBC(Self%Ne,self%DimAnalysis*Self%Npe))
        Self%ConnectivityNPBC = Self%ConnectivityN
        Self%ConnectivityDPBC = 0
        l = size(CondArray(:,1))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = 1, l, 1
                    if (Self%ConnectivityNPBC(i,j).eq.CondArray(k,2)) then
                        Self%ConnectivityNPBC(i,j) = CondArray(k,1)
                    end if
                end do
            end do
        end do
        ! get the dof
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = self%DimAnalysis-1, 0, -1
                    Self%ConnectivityDPBC(i,j*self%DimAnalysis-k) = Self%ConnectivityNPBC(i,j)*self%DimAnalysis - k
                end do
            end do
        end do
        !call FilePrinting(Self%ConnectivityNPBC,'DataResults/.InternalData/ConnectivityNPBC.txt')
        !call FilePrinting(Self%ConnectivityDPBC,'DataResults/.InternalData/ConnectivityDPBC.txt')
    end subroutine GetMasterSlaveNodes
    ! 14. Input Boundary Conditions
    subroutine SetBondaryConditions(Self,Path)
        implicit none
        class(Structure), intent(inout)                             :: Self
        character(len=*), intent(in)                                :: Path
        ! internal variables
        integer                                                     :: i,j,k,l,Freedof,Fixeddof
        logical                                                     :: r1,r2,r3,r4
        integer, dimension(:), allocatable                          :: InPosition
        logical, dimension(:), allocatable                          :: InLogical
        ! ---- 1st part ----
        ! Get the degress of freedom
        call ReadingfileInteger(Path,Self%NBc,(Self%DimAnalysis+1),Self%BoundaryC)
        l = 0
        do i = 1, Self%N, 1
            r1 = any(Self%CornerNodes(:,2).eq.i)    ! if it is a slave node (corner)
            r2 = any(Self%EdgeNodes(:,2).eq.i)      ! if it is a slave node (edges)
            r3 = any(Self%FaceNodes(:,2).eq.i)      ! if it is a slave node (faces)
            r4 = any(Self%BoundaryC(:,1).eq.i)      ! if it is a constrain node
            if (any([r1,r2,r3,r4])) then            ! GL fixed
                l = l + 1
            end if
        end do
        Fixeddof = l*Self%DimAnalysis
        Freedof = (Self%N-l)*Self%DimAnalysis
        allocate(Self%FreeD(Freedof));     Self%FreeD = 0
        allocate(Self%FixedD(Fixeddof));  Self%FixedD = 0
        k = 1
        l = 1
        do i = 1, Self%N, 1
            r1 = any(Self%CornerNodes(:,2).eq.i)    ! if it is a slave node (corner)
            r2 = any(Self%EdgeNodes(:,2).eq.i)      ! if it is a slave node (edges)
            r3 = any(Self%FaceNodes(:,2).eq.i)      ! if it is a slave node (faces)
            r4 = any(Self%BoundaryC(:,1).eq.i)      ! if it is a constrain node
            if (any([r1,r2,r3,r4])) then  ! fixed
                do j = Self%DimAnalysis-1, 0, -1
                    Self%FixedD(k*Self%DimAnalysis-j) = i*Self%DimAnalysis - j
                end do
                k = k + 1
            else    ! Free
                do j = Self%DimAnalysis-1, 0, -1
                    Self%FreeD(l*Self%DimAnalysis-j) = i*Self%DimAnalysis - j
                end do
                l = l + 1
            end if
        end do
        ! printing
        !call FilePrinting(Self%FreeD,'V','DataResults/.InternalData/FreeDof.txt')
        !call FilePrinting(Self%FixedD,'V','DataResults/.InternalData/FixedDof.txt')
    end subroutine SetBondaryConditions
    ! 15 Read Files
    subroutine ReadFiles(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        call SetCoordinates(Self,'DataStructure/Coordinates.txt')
        call SetConnectivity(Self,'DataStructure/Connectivity.txt')
        call SetCornerNodes(Self,'DataStructure/CornerNodes.txt')
        call SetEdgeNodes(Self,'DataStructure/EdgeNodes.txt')
        call SetFaceNodes(Self,'DataStructure/FaceNodes.txt')
        call GetMasterSlaveNodes(Self)
        call SetBondaryConditions(Self,'DataStructure/BoundaryConditions.txt')
    end subroutine ReadFiles
    ! ----------------------------------------------------------------- !
    !              subroutines that define the FEM procedure            !
    ! ----------------------------------------------------------------- !
    ! 1. Strain
    subroutine StrainDisplacement(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal values
        integer                                                         :: dim,nl,i,j,k
        double precision, dimension(:), allocatable                     :: Vector
        double precision, dimension(:,:,:), allocatable                 :: Array
        if (Self%DimAnalysis.eq.2) then
            dim = 2
            nl = 3
            allocate(Array(nl,dim,dim))
            Array = 0.0d0
            Array(1,1,1) = Self%StrainH
            Array(2,2,2) = Self%StrainH
            Array(3,2,1) = 0.5d0*Self%StrainH
            Array(3,1,2) = 0.5d0*Self%StrainH
        elseif (Self%DimAnalysis.eq.3) then 
            dim = 3
            nl = 6
            allocate(Array(nl,dim,dim))
            Array = 0.0d0
            Array(1,1,1) = Self%StrainH
            Array(2,2,2) = Self%StrainH
            Array(3,3,3) = Self%StrainH
            Array(4,2,1) = 0.5d0*Self%StrainH
            Array(4,1,2) = 0.5d0*Self%StrainH
            Array(5,3,2) = 0.5d0*Self%StrainH
            Array(5,2,3) = 0.5d0*Self%StrainH
            Array(6,3,1) = 0.5d0*Self%StrainH
            Array(6,1,3) = 0.5d0*Self%StrainH
        end if
        if(allocated(Self%UGlobal_0)) then
            Self%UGlobal_0 = 0
        else
            allocate(Self%UGlobal_0(Self%N*dim,nl))
            Self%UGlobal_0 = 0
        end if 
        do i = 1, Self%N, 1
            do j = 1, nl, 1
                Vector = matmul(Self%Coordinates(i,1:dim),Array(j,:,:))
                do k = dim-1,0,-1
                    Self%UGlobal_0(dim*i-k,j) = Vector(dim-k)
                end do
            end do
        end do
        deallocate(Array,Vector)
    end subroutine StrainDisplacement
    ! 2. Derivative of the shape functions
    subroutine DiffFormFunction(Self,DiffFunction,e,n,z)
        ! dN1/de     dN2/de     dN3/de      ....     dNn/de
        ! dN1/dn     dN2/dn     dN3/dn      ....     dNn/dn
        ! dN1/dz     dN2/dz     dN3/dz      ....     dNn/dz (caso 3D)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: e, n
        double precision, intent(inout), optional                       :: z
        double precision, dimension(:,:), allocatable, intent(out)      :: DiffFunction
        if (Self%ElementType.eq.'tria3') then
            allocate(DiffFunction(2,3))
            DiffFunction(1,:) = [-1.0d0,1.0d0,0.0d0]
            DiffFunction(2,:) = [-1.0d0,0.0d0,1.0d0]
        elseif (Self%ElementType.eq.'tria6') then
            allocate(DiffFunction(2,6))
            DiffFunction(1,:) = [4.0d0*e+4.0d0*n-3.0d0,4.0d0*e-1.0d0,0.0d0,4.0d0-4.0d0*n-8.0d0*e,4.0d0*n,-4.0d0*n]
            DiffFunction(2,:) = [4.0d0*e+4.0d0*n-3.0d0,0.0d0,4.0d0*n-1.0d0,-4.0d0*e,4.0d0*e,4.0d0-8.0d0*n-4.0d0*e]
        elseif (Self%ElementType.eq.'quad4') then
            allocate(DiffFunction(2,4))
            DiffFunction(1,:) = [n/4.0d0-1.0d0/4.0d0,1.0d0/4.0d0-n/4.0d0,n/4.0d0+1.0d0/4.0d0,-n/4.0d0-1.0d0/4.0d0]
            DiffFunction(2,:) = [e/4.0d0-1.0d0/4.0d0,-e/4.0d0-1.0d0/4.0d0,e/4.0d0+1.0d0/4.0d0,1.0d0/4.0d0-e/4.0d0]
        elseif (Self%ElementType.eq.'quad8') then
            allocate(DiffFunction(2,8))
            DiffFunction(1,:) = [-(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+1.0d0))/4.0d0, &
                            ((n-1.0d0)*(n-e+1.0d0))/4.0d0-(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n-1.0d0))/4.0d0, &
                            (e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)+((n+1.0d0)*(e-n+1.0d0))/4.0d0,  &
                            e*(n-1.0d0),1.0d0/2.0d0-n**2.0d0/2.0d0,-e*(n+1.0d0),n**2.0d0/2.0d0-1.0d0/2.0d0]
            DiffFunction(2,:) = [-(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)-(e/4.0d0-1.0d0/4.0d0)*(e+n+1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n-e+1.0d0)+(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)+(e/4.0d0+1.0d0/4.0d0)*(e+n-1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(e-n+1.0d0)-(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            e**2.0d0/2.0d0-1.0d0/2.0d0,-n*(e+1.0d0),1.0d0/2.0d0-e**2.0d0/2.0d0,n*(e-1.0d0)]
        elseif (Self%ElementType.eq.'tetra4') then
            allocate(DiffFunction(3,4))
            DiffFunction(1,:) = [1.0d0,0.0d0,0.0d0,-1.0d0]
            DiffFunction(2,:) = [0.0d0,1.0d0,0.0d0,-1.0d0]
            DiffFunction(3,:) = [0.0d0,0.0d0,1.0d0,-1.0d0]
        elseif (Self%ElementType.eq.'tetra10') then
            allocate(DiffFunction(3,10))
            DiffFunction(1,:) = [4.0d0*e-1.0d0,0.0d0,0.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,4.0d0*n,0.0d0,4.0d0*z, &
                                4.0d0-4.0d0*n-4.0d0*z-8.0d0*e,-4.0d0*n,0.0d0]
            DiffFunction(2,:) = [0.0d0,4.0d0*n-1.0d0,0.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,4.0d0*e,4.0d0*z,0.0d0,-4.0d0*e, &
                                4.0d0-8.0d0*n-4.0d0*z-4.0d0*e,0.0d0]
            DiffFunction(3,:) = [0.0d0,0.0d0,4.0d0*z-1.0d0,4.0d0*e+4.0d0*n+4.0d0*z-3.0d0,0.0d0,4.0d0*n,4.0d0*e,-4.0d0*e, &
                                -4.0d0*n,0.0d0]
        elseif (Self%ElementType.eq.'hexa8') then
            allocate(DiffFunction(3,8))
            DiffFunction(1,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)+((n-1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+z))/8.0d0, &
                             (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)-((n+1.0d0)*(e+n+z))/8.0d0, &
                            -(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)-((n-1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                             (e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)+((n-1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                            -(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)-((n+1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                             (e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)+((n+1.0d0)*(e+n+z-2.0d0))/8.0d0]
            DiffFunction(2,:) = [(e/8.0d0-1.0d0/8.0d0)*(e+n+z)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(e+n+z)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(e+n+z)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(e+n+z)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(e+n+z-2.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(e+n+z-2.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(e+n+z-2.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0), &
                             (e/8.0d0-1.0d0/8.0d0)*(e+n+z-2.0d0)+(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)]
            DiffFunction(3,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0),-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                             (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0),-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0),(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0),(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)]
        elseif (Self%ElementType.eq.'hexa20') then
            allocate(DiffFunction(3,20))
            !  line 1
            DiffFunction(1,1) = - (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,2) = (e*n*z*(n - 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,3) = - (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,4) = (e*n*z*(n + 1.0d0)*(z - 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(1,5) = (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,6) = - (e*n*z*(n - 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,7) = (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 + (n*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,8) = - (e*n*z*(n + 1.0d0)*(z + 1.0d0))/8.0d0 - (n*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(1,9) = -(e*(n - 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,10) = (e*(n + 1.0d0)*(z - 1.0d0))/2.0d0 
            DiffFunction(1,11) = -(e*(n + 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,12) = (e*(n - 1.0d0)*(z + 1.0d0))/2.0d0 
            DiffFunction(1,13) = -((n**2 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,14) = ((n**2 - 1.0d0)*(z - 1.0d0))/4.0d0 
            DiffFunction(1,15) = -((n**2 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,16) = ((n**2 - 1.0d0)*(z + 1.0d0))/4.0d0 
            DiffFunction(1,17) = -((z**2 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,18) = ((z**2 - 1.0d0)*(n - 1.0d0))/4.0d0 
            DiffFunction(1,19) = -((z**2 - 1.0d0)*(n + 1.0d0))/4.0d0 
            DiffFunction(1,20) = 0.0d0 
            !  line 2
            DiffFunction(2,1) = - (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,2) = (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,3) = - (e*n*z*(e + 1.0d0)*(z - 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,4) = (e*n*z*(e - 1.0d0)*(z - 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(2,5) = (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,6) = - (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,7) = (e*n*z*(e + 1.0d0)*(z + 1.0d0))/8.0d0 + (e*z*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,8) = - (e*n*z*(e - 1.0d0)*(z + 1.0d0))/8.0d0 - (e*z*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(2,9) = -(e**2/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,10) = (e**2/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,11) = -(e**2/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,12) = (e**2/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,13) = -2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,14) = 2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z - 1.0d0) 
            DiffFunction(2,15) = -2.0d0*n*(e/4.0d0 + 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,16) = 2.0d0*n*(e/4.0d0 - 1.0d0/4.0d0)*(z + 1.0d0) 
            DiffFunction(2,17) = -(e/4.0d0 - 1.0d0/4.0d0)*(z**2 - 1.0d0) 
            DiffFunction(2,18) = (e/4.0d0 + 1.0d0/4.0d0)*(z**2 - 1.0d0) 
            DiffFunction(2,19) = -(e/4.0d0 + 1.0d0/4.0d0)*(z**2 - 1.0d0) 
            DiffFunction(2,20) = 0.0d0 
            !  line 3
            DiffFunction(3,1) = - (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,2) = (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n - 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,3) = - (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,4) = (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n + 1.0d0)*(z - 1.0d0))/8.0d0 
            DiffFunction(3,5) = (e*n*z*(e - 1.0d0)*(n - 1.0d0))/8.0d0 + (e*n*(e - 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,6) = - (e*n*z*(e + 1.0d0)*(n - 1.0d0))/8.0d0 - (e*n*(e + 1.0d0)*(n - 1.0d0)*(z + 1.0d0))/8.0d0
            DiffFunction(3,7) = (e*n*z*(e + 1.0d0)*(n + 1.0d0))/8.0d0 + (e*n*(e + 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,8) = - (e*n*z*(e - 1.0d0)*(n + 1.0d0))/8.0d0 - (e*n*(e - 1.0d0)*(n + 1.0d0)*(z + 1.0d0))/8.0d0 
            DiffFunction(3,9) = -(e**2/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,10) = (e**2/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,11) = -(e**2/4.0d0 - 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,12) = (e**2/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,13) = -(e/4.0d0 - 1.0d0/4.0d0)*(n**2 - 1.0d0) 
            DiffFunction(3,14) = (e/4.0d0 + 1.0d0/4.0d0)*(n**2 - 1.0d0) 
            DiffFunction(3,15) = -(e/4.0d0 + 1.0d0/4.0d0)*(n**2 - 1.0d0) 
            DiffFunction(3,16) = (e/4.0d0 - 1.0d0/4.0d0)*(n**2 - 1.0d0) 
            DiffFunction(3,17) = -2.0d0*z*(e/4.0d0 - 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,18) = 2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n - 1.0d0) 
            DiffFunction(3,19) = -2.0d0*z*(e/4.0d0 + 1.0d0/4.0d0)*(n + 1.0d0) 
            DiffFunction(3,20) = 0.0d0 
        end if
    end subroutine DiffFormFunction
    ! 3. Elasticity tensor
    subroutine ElasticityTensor(Self,ETensor)
        implicit none
        class(Structure), intent(inout)                                 :: Self    
        double precision, dimension(:,:), allocatable, intent(inout)    :: ETensor
        ! internal variables
        double precision                                                :: E,V,Constant
        E = Self%YoungModulus
        V = Self%PoissonModulus    
        if (Self%AnalysisType.eq.'PlaneStress') then
            allocate(ETensor(3,3))
            Constant = E/(1.0d0-V**2d0)
            ETensor(1,:) = [1.0d0,V,0.0d0]
            ETensor(2,:) = [V,1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'PlaneStrain') then
            allocate(ETensor(3,3))
            Constant = E*(1.0d0-v)/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0,V/(1.0d0-v),0.0d0]
            ETensor(2,:) = [V/(1.0d0-v),1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v))]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'HookeLaw') then
            allocate(ETensor(6,6))
            Constant = E*(1.0d0-v)/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0,V/(1.0d0-v),V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(2,:) = [V/(1.0d0-v),1.0d0,V/(1.0d0-v),0.0d0,0.0d0,0.0d0]
            ETensor(3,:) = [V/(1.0d0-v),V/(1.0d0-v),1.0d0,0.0d0,0.0d0,0.0d0]
            ETensor(4,:) = [0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0,0.0d0]
            ETensor(5,:) = [0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v)),0.0d0]
            ETensor(6,:) = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/(2.0d0*(1.0d0-v))]
            ETensor = Constant*ETensor
        end if
    end subroutine ElasticityTensor
    ! 4. pressembly routine
    subroutine PreAssemblyRoutine(self)
        use omp_lib
        implicit none
        integer, dimension(:,:), allocatable                        :: Matrix
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: i,j,k,i1,IndexRow,IndexCol
        integer, dimension(:), allocatable                          :: InPosition
        logical, dimension(:), allocatable                          :: InLogical
        integer, dimension(:,:), allocatable                        :: Node_Interaction,Elem_Interaction
        ! -------------------------------------------------------------------------------------- !
        ! note: In this part, a preliminary assembly of the stiffness matrix (in sparse form)    !
        !       is made. The idea is to locate the position of each element of the matrix in     !
        !       the sparse vector so that the assembly of the global matrix in each iteration    !
        !       can be faster.                                                                   !
        ! -------------------------------------------------------------------------------------- !
        write(unit=*, fmt=*) '- Assembly routine'
        ! get de element and DoF incidence for each DoF (only the free degres)
        j = size(Self%FreeD)
        ! maximum of 20 elements per node
        allocate(Elem_Interaction(j,50))
        Elem_Interaction = 0
        ! considering the quad20 30*20 
        allocate(Node_Interaction(j,600))
        Node_Interaction = 0

        !$omp parallel do private(InLogical,InPosition)
        do i = 1, size(Self%FreeD), 1
            ! element interaction
            InLogical = any(Self%ConnectivityDPBC.eq.Self%FreeD(i),2)
            InPosition = pack([(k,k=1,Self%Ne)],InLogical)
            j = count(InLogical)
            Elem_Interaction(i,1:j) = InPosition 
            ! DoF interaction
            InPosition = reshape(Self%ConnectivityDPBC(Elem_Interaction(i,1:j),:),[j*(Self%Npe)*Self%DimAnalysis])
            ! letting only free degres
            do k = 1, size(InPosition), 1       
                if(any(InPosition(k).eq.Self%FixedD)) InPosition(k) = 0
            end do
            ! Sorting ang eliminating 0s
            call sort(InPosition,j*Self%Npe*Self%DimAnalysis)
            InLogical = InPosition.ge.Self%FreeD(i)
            InPosition = pack(InPosition,InLogical)
            j = size(InPosition)
            Node_Interaction(i,1:j) = InPosition
        end do
        !$omp end parallel do        

        !call FilePrinting(Elem_Interaction,'DataResults/.InternalData/Elem_Interaction.txt')
        !call FilePrinting(Node_Interaction,'DataResults/.InternalData/Node_Interaction.txt')

        ! pre-assembly
        i = count(Node_Interaction.ne.0)
        allocate(Self%index_i_KGlobal(i)) ; Self%index_i_KGlobal = 0
        allocate(Self%index_j_KGlobal(i)) ; Self%index_j_KGlobal = 0
        allocate(Self%Location_KGlobal(Self%Ne,Self%Npe*Self%DimAnalysis,Self%Npe*Self%DimAnalysis))
        Self%Location_KGlobal = 0
        k = 1
        do i = 1, size(Self%FreeD), 1                           ! Col
            do j = 1, count(Node_Interaction(i,:).ne.0), 1      ! row
                Self%index_i_KGlobal(k) = findloc(Self%FreeD,Node_Interaction(i,j),1)   ! row-indx
                Self%index_j_KGlobal(k) = i                                             ! col-indx
                do i1 = 1, count(Elem_Interaction(i,:).ne.0), 1
                    IndexRow = findloc(Self%ConnectivityDPBC(Elem_Interaction(i,i1),:),Node_Interaction(i,j),1)
                    IndexCol = findloc(Self%ConnectivityDPBC(Elem_Interaction(i,i1),:),Node_Interaction(i,1),1)
                    if (IndexCol.eq.0.or.IndexRow.eq.0) cycle
                    Self%Location_KGlobal(Elem_Interaction(i,i1),IndexRow,IndexCol) = k
                end do
                k = k+1
            end do
        end do
    end subroutine PreAssemblyRoutine
    ! 5. Local force vectors
    subroutine SetFlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                         :: Self
        double precision, intent(inout)                                         :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)              :: DensityVector 
        ! internal variables
        integer                                                                 :: el,i,j,k,l,m
        double precision                                                        :: e,n,z,w1,w2,w3,DetJacobian,Fac
        double precision, dimension(:,:), allocatable                           :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable                           :: DiffN,DiffNXY,ElementCoordinates
        double precision, dimension(:,:), allocatable                           :: Epsilon_0
        if (Self%DimAnalysis.eq.2) then
            if (allocated(Self%FLocal)) then
                Self%FLocal = 0.0d0
            else
                allocate(Self%FLocal(Self%Ne,Self%Npe*2,3));      Self%FLocal = 0.0d0;
            end if
            allocate(Be(3,2*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(2,2));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(2,2));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,2));      ElementCoordinates = 0.0d0;
            allocate(Epsilon_0(3,3))
            ! assumed strain for computation of the state of loading
            Epsilon_0 = 0.0d0
            Epsilon_0(1,1)=Self%StrainH
            Epsilon_0(2,2)=Self%StrainH
            Epsilon_0(3,3)=Self%StrainH
            do el = 1, Self%Ne, 1
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        ! Jacobian
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = abs(Determinant(Jacobian))
                        DiffNXY = matmul(InvJacobian,DiffN)
                        Fac = DetJacobian*w1*w2*(Self%Thickness)
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            Be(1,2*k-1) = DiffNxy(1,k)
                            Be(2,2*k) = DiffNxy(2,k)
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! F Local
                        Self%FLocal(el,:,:) = Self%FLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Epsilon_0)))
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            if (allocated(Self%FLocal)) then
                Self%FLocal = 0.0d0
            else
                allocate(Self%FLocal(Self%Ne,Self%Npe*3,6));      Self%FLocal = 0.0d0;
            end if
            allocate(Be(6,3*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(3,3));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(3,3));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,3));      ElementCoordinates = 0.0d0;
            allocate(Epsilon_0(6,6))
            ! assumed strain for computation of the state of loading
            Epsilon_0 = 0.0d0
            Epsilon_0(1,1)=Self%StrainH
            Epsilon_0(2,2)=Self%StrainH
            Epsilon_0(3,3)=Self%StrainH
            Epsilon_0(4,4)=Self%StrainH
            Epsilon_0(5,5)=Self%StrainH
            Epsilon_0(6,6)=Self%StrainH
            do el = 1, Self%Ne, 1
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            ! Jacobian
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = abs(Determinant(Jacobian))
                            DiffNXY = matmul(InvJacobian,DiffN)
                            Fac = DetJacobian*w1*w2*w3
                            ! Be
                            Be = 0.0d0
                            do l = 1, size(DiffN,2), 1
                                Be(1,3*l-2) = DiffNxy(1,l)
                                Be(2,3*l-1) = DiffNxy(2,l)
                                Be(3,3*l) = DiffNxy(3,l)
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                Be(6,3*l-2) = DiffNxy(3,l)
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! K Local
                            Self%FLocal(el,:,:) = Self%FLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Epsilon_0)))
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates,Epsilon_0)
    end subroutine SetFlocal
    ! 6. Local Stiffness matrix
    subroutine SetKlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)      :: DensityVector 
        ! internal variables
        integer                                                         :: el,i,j,k,l,m
        double precision                                                :: e,n,z,w1,w2,w3,DetJacobian,Fac
        double precision, dimension(:,:), allocatable                   :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable                   :: DiffN,DiffNXY,ElementCoordinates
        ! ---- 1st part ----
        ! in this section the local stiffness matrix (dense) of each of the elements is calculated using gauss 
        ! quadrature and each of them is stored in the array KLocal[element,[K_x,K_y]].
        !open(unit=1,file='DataStructure/AdditionalData/StiffnessMatrixArray.txt',iostat=ios,status="replace",action="write")
        ! this is just to check
        if (Self%DimAnalysis.eq.2) then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*2,Self%Npe*2)); Self%KLocal = 0.0d0;
                allocate(Self%BLocal(Self%Ne,3,2*Self%Npe));          Self%BLocal = 0.0d0;
            end if
            allocate(Be(3,2*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(2,2));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(2,2));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,2));      ElementCoordinates = 0.0d0;
            do el = 1, Self%Ne, 1
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:2)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        call DiffFormFunction(Self,DiffN,e,n)
                        ! Jacobian
                        Jacobian = matmul(DiffN,ElementCoordinates)
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = abs(Determinant(Jacobian))
                        DiffNXY = matmul(InvJacobian,DiffN)
                        Fac = DetJacobian*w1*w2*(Self%Thickness)
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            Be(1,2*k-1) = DiffNxy(1,k)
                            Be(2,2*k) = DiffNxy(2,k)
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! Be Local
                        Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Fac*Be
                        ! K Local
                        Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Be)))
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.3) then
            if (allocated(Self%KLocal).and.allocated(Self%BLocal)) then
                Self%KLocal = 0.0d0
                Self%BLocal = 0.0d0
            else
                allocate(Self%KLocal(Self%Ne,Self%Npe*3,Self%Npe*3)); Self%KLocal = 0.0d0;
                allocate(Self%BLocal(Self%Ne,6,3*Self%Npe));          Self%BLocal = 0.0d0;
            end if
            allocate(Be(6,3*Self%Npe));                                    Be = 0.0d0;
            allocate(Jacobian(3,3));                                 Jacobian = 0.0d0;
            allocate(InvJacobian(3,3));                           InvJacobian = 0.0d0;
            allocate(ElementCoordinates(Self%Npe,3));      ElementCoordinates = 0.0d0;
            do el = 1, Self%Ne, 1
                call ElasticityTensor(Self,D)
                D = (DensityVector(el)**PenalFactor)*D
                ElementCoordinates = Self%Coordinates(Self%ConnectivityN(el,:),1:3)
                do i = 1, Self%QuadGauss, 1
                    e = Self%GaussPoint(i)
                    w1 = Self%GaussWeights(i)
                    do j = 1, Self%QuadGauss, 1
                        n = Self%GaussPoint(j)
                        w2 = Self%GaussWeights(j)
                        do k = 1, Self%QuadGauss, 1
                            z = Self%GaussPoint(k)
                            w3 = Self%GaussWeights(k)
                            call DiffFormFunction(Self,DiffN,e,n,z)
                            ! Jacobian
                            Jacobian = matmul(DiffN,ElementCoordinates)
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = abs(Determinant(Jacobian))
                            DiffNXY = matmul(InvJacobian,DiffN)
                            Fac = DetJacobian*w1*w2*w3
                            ! Be
                            Be = 0.0d0
                            do l = 1, size(DiffN,2), 1
                                Be(1,3*l-2) = DiffNxy(1,l)
                                Be(2,3*l-1) = DiffNxy(2,l)
                                Be(3,3*l) = DiffNxy(3,l)
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                Be(6,3*l-2) = DiffNxy(3,l)
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! Be Local
                            Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + Fac*Be
                            ! K Local
                            Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + Fac*(matmul(transpose(Be),matmul(D,Be)))
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if     
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates)
    end subroutine SetKlocal
    ! 7. Assembly Load vector
    subroutine GetFGlobalSparse(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: i
        if(allocated(Self%value_FGlobal)) then
                Self%value_FGlobal = 0.0d0
        else
            if (Self%DimAnalysis.eq.2) then
                allocate(Self%value_FGlobal(Self%DimAnalysis*Self%N,3)); 
                Self%value_FGlobal = 0.0d0
            elseif (Self%DimAnalysis.eq.3) then
                allocate(Self%value_FGlobal(Self%DimAnalysis*Self%N,6)); 
                Self%value_FGlobal = 0.0d0
            end if
        end if 
        do i = 1, Self%Ne, 1
            Self%value_FGlobal(Self%ConnectivityDPBC(i,:),:) = Self%value_FGlobal(Self%ConnectivityDPBC(i,:),:) &
                                                             + Self%FLocal(i,:,:)
        end do
        !call FilePrinting(Self%value_FGlobal,'DataResults/.InternalData/value_FGlobal.txt')
    end subroutine GetFGlobalSparse
    ! 8. Assembly Stiffness matrix
    subroutine GetKGlobalSparse(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: n,i,j,k
        logical, dimension(:), allocatable                          :: InLogical
        if (allocated(Self%value_KGlobal)) then
            Self%value_KGlobal = 0.0
        else
            n = size(Self%index_i_KGlobal)
            allocate(Self%value_KGlobal(n))
            Self%value_KGlobal = 0.0
        end if
        do i = 1, Self%Ne, 1
            do j = 1, Self%DimAnalysis*Self%Npe, 1
                do k = 1, Self%DimAnalysis*Self%Npe, 1
                    n = Self%Location_KGlobal(i,j,k)
                    if(n.eq.0) cycle
                    Self%value_KGlobal(n) = Self%value_KGlobal(n) + Self%KLocal(i,j,k)
                end do
            end do
        end do
        ! eliminating null or zero elements (StiffnessMatrix=0.0)
        InLogical = Self%value_KGlobal.ne.0.0
        Self%Rows_KGlobal = pack(Self%index_i_KGlobal,InLogical)
        Self%Cols_KGlobal = pack(Self%index_j_KGlobal,InLogical)
        Self%value_KGlobal = pack(Self%value_KGlobal,InLogical)
    end subroutine GetKGlobalSparse
    ! 9. Upload Data
    subroutine UploadStructure(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                             :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)  :: DensityVector 
        ! 1. Getting Local Stiffnes matrix
        !write(unit=*, fmt=*) '- Local Stiffness matrix'
        call SetFlocal(Self,DensityVector,PenalFactor)
        ! 2. Getting Local Stiffnes matrix
        !write(unit=*, fmt=*) '- Global Stiffness matrix'
        call SetKlocal(Self,DensityVector,PenalFactor)
        ! 3. Getting Local Stiffnes matrix
        !write(unit=*, fmt=*) '- Local Force Vector'
        call GetFGlobalSparse(Self)
        ! 4. Getting Local Stiffnes matrix
        !write(unit=*, fmt=*) '- Global Force Vector'
        call GetKGlobalSparse(Self)
        Self%value_FGlobal = Self%value_FGlobal(Self%FreeD,:)
    end subroutine UploadStructure
    ! 10. MA87 Sparse Solver
    subroutine HSL87Solver(Self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        ! internal variables
        if(allocated(Self%UGlobal).and.allocated(Self%value_UGlobal)) then
            Self%UGlobal = 0.0d0
            Self%value_UGlobal = 0.0d0
        else
            if (Self%DimAnalysis.eq.2) then
                allocate(Self%UGlobal(Self%N*Self%DimAnalysis,3))
                allocate(Self%value_UGlobal(size(Self%FreeD),3))
                Self%UGlobal = 0.0d0; Self%value_UGlobal = 0.0d0;
            elseif (Self%DimAnalysis.eq.3) then
                allocate(Self%UGlobal(Self%N*Self%DimAnalysis,6))
                allocate(Self%value_UGlobal(size(Self%FreeD),6));
                Self%UGlobal = 0.0d0; Self%value_UGlobal = 0.0d0;
            end if
        end if
        !write(unit=*, fmt=*) '- Applying MA87 Sparse Solver'
        Self%value_UGlobal = SparseSystemSolver(Self%Rows_KGlobal,Self%Cols_KGlobal,Self%value_KGlobal,Self%value_FGlobal)
        deallocate(Self%value_KGlobal,Self%value_FGlobal)
        Self%UGlobal(Self%FreeD,:) = Self%value_UGlobal
        !write(unit=*, fmt=*) '- Extending solution to slave nodes'
        call ExtendSolutionPBC(Self)
        !call FilePrinting(Self%UGlobal,'DataResults/.InternalData/UGlobalBase.txt')
        Call StrainDisplacement(Self)
        Self%UGlobal = Self%UGlobal_0 - Self%UGlobal
        !call FilePrinting(Self%UGlobal,'DataResults/.InternalData/UGlobal.txt')
    end subroutine HSL87Solver
    ! 11. Extending solution to Slave nodes
    subroutine ExtendSolutionPBC(self)
        implicit none
        class(Structure), intent(inout)                             :: Self
        integer                                                     :: i,j,ind1,ind2
        ! corner
        do i = 1, size(Self%CornerNodes(:,1)), 1
            do j = Self%DimAnalysis-1,0,-1
                ind1 = Self%CornerNodes(i,1)*Self%DimAnalysis-j    !master
                ind2 = Self%CornerNodes(i,2)*Self%DimAnalysis-j    !slave
                if ((Self%CornerNodes(i,1).eq.0).or.(Self%CornerNodes(i,2).eq.0)) cycle
                Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
            end do
        end do
        ! edge
        do i = 1, size(Self%EdgeNodes(:,1)), 1
            do j = Self%DimAnalysis-1,0,-1
                ind1 = Self%EdgeNodes(i,1)*Self%DimAnalysis-j      !master
                ind2 = Self%EdgeNodes(i,2)*Self%DimAnalysis-j      !slave
                if ((Self%EdgeNodes(i,1).eq.0).or.(Self%EdgeNodes(i,2).eq.0)) cycle
                Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
            end do
        end do
        ! face
        if (Self%DimAnalysis.eq.3) then
            do i = 1, size(Self%FaceNodes(:,1)), 1
                do j = Self%DimAnalysis-1,0,-1
                    ind1 = Self%FaceNodes(i,1)*Self%DimAnalysis-j  !master
                    ind2 = Self%FaceNodes(i,2)*Self%DimAnalysis-j  !slave
                    if ((Self%FaceNodes(i,1).eq.0).or.(Self%FaceNodes(i,2).eq.0)) cycle
                    Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
                end do
            end do
        end if
    end subroutine ExtendSolutionPBC
    ! 12. Numerical homogenization 
    subroutine NumericalHomogenization(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)      :: DensityVector 
        ! internal values
        integer                                                         :: i,j,Vtotal
        double precision, dimension(:,:), allocatable                   :: D,Dtotal,Displacement
        if (allocated(Self%DTensor)) then
            Self%DTensor = 0.0d0
            if (Self%DimAnalysis.eq.2) allocate(Dtotal(3,3))  
            if (Self%DimAnalysis.eq.3) allocate(Dtotal(6,6)) 
            Dtotal = 0.0d0
        else
            if (Self%DimAnalysis.eq.2) then
                allocate(Self%DTensor(Self%Ne,3,3));            Self%DTensor = 0.0d0
                allocate(Dtotal(3,3));                                Dtotal = 0.0d0
            elseif (Self%DimAnalysis.eq.3) then
                allocate(Self%DTensor(Self%Ne,6,6));            Self%DTensor = 0.0d0
                allocate(Dtotal(6,6));                                Dtotal = 0.0d0
            else
                stop "ERROR Numerical homogenization"
            end if
        end if
        Vtotal = sum(Self%VolumePerElement)
        do i = 1, Self%Ne, 1
            Displacement = Self%UGlobal(Self%ConnectivityD(i,:),:)
            Self%DTensor(i,:,:) = (matmul(transpose(Displacement),matmul(Self%KLocal(i,:,:),Displacement)))
            Self%DTensor(i,:,:) = (1.0d0/Vtotal)*Self%DTensor(i,:,:)
            if (Self%DimAnalysis.eq.2) then
                Self%DTensor(i,1:2,3)=0.0d0
                Self%DTensor(i,3,1:2)=0.0d0
            elseif (Self%DimAnalysis.eq.3) then
                Self%DTensor(i,1:3,4:6)=0.0d0
                Self%DTensor(i,4:6,1:3)=0.0d0
                Self%DTensor(i,4,5:6)=0.0d0
                Self%DTensor(i,5:6,4)=0.0d0
                Self%DTensor(i,5,6)=0.0d0
                Self%DTensor(i,6,5)=0.0d0
            end if
            Dtotal = Dtotal + Self%DTensor(i,:,:)
        end do
        !call FilePrinting(Dtotal,'DataResults/.InternalData/HomogenizedTensor.txt')
        !deallocate(D,Dtotal)
    end subroutine NumericalHomogenization
    ! 13. Processing numerical results
    subroutine ProcessingResults(Self)
        ! this subroutine calculates stresses per element, per node, deformations and strain energy.
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal variables
        integer                                                         :: i,j,k,nl
        integer, dimension(:), allocatable                              :: index
        double precision                                                :: energy
        double precision, dimension(:), allocatable                     :: Strain, Stress
        double precision, dimension(:,:), allocatable                   :: D
        call ElasticityTensor(Self,D)
        !write(unit=*, fmt=*) '- Processing Results'
        ! allocating
        if (Self%DimAnalysis.eq.2) then
            nl = 3
            allocate(Strain(nl),Stress(nl))
            Strain = 0.0d0; Stress = 0.0d0;
            if (allocated(Self%StrainE)) then
                Self%StrainE = 0.0d0
            else
                allocate(Self%StrainE(nl,Self%Ne,nl))        
                Self%StrainE = 0.0d0
            end if
            if (allocated(Self%StressE)) then
                Self%StressE = 0.0d0
            else
                allocate(Self%StressE(nl,Self%Ne,nl))
                Self%StressE = 0.0d0
            end if
        elseif (Self%DimAnalysis.eq.3) then
            nl = 6
            allocate(Strain(nl),Stress(nl))
            Strain = 0.0d0; Stress = 0.0d0;
            if (allocated(Self%StrainE)) then
                Self%StrainE = 0.0d0
            else
                allocate(Self%StrainE(nl,Self%Ne,nl))
                Self%StrainE = 0.0d0
            end if
            if (allocated(Self%StressE)) then
                Self%StressE = 0.0d0
            else
                allocate(Self%StressE(nl,Self%Ne,nl))
                Self%StressE = 0.0d0
            end if
        end if
        if (allocated(Self%StrainEnergyE)) then
            Self%StrainEnergyE = 0.0d0
        else
            allocate(Self%StrainEnergyE(nl,Self%Ne)); 
            Self%StrainEnergyE = 0.0d0
        end if
        ! star processing
        ! [LoadState,Element/Node,Stress/Srain-Direction]
        ! (calculation per element)
        do i = 1, nl, 1
            do j = 1, Self%Ne, 1
                ! 1. Strain     Strain_e = Be*Ue        
                ! exx - eyy - ezz(3D) - exy - eyz(3D) - exz(3D)
                Strain = matmul(Self%BLocal(j,:,:),Self%UGlobal(Self%ConnectivityD(j,:),i))
                Self%StrainE(i,j,:) = Strain
                ! 2. Stress      Sigma_e = D*Strain_e   
                ! Sxx - Syy - Szz(3D) - Sxy - Syz(3D) - Sxz(3D)
                Stress = matmul(D,Strain)
                Self%StressE(i,j,:) = Stress
                ! 3. Energy    Ee = 0.5*Sigma_e*Strain_e
                energy = 0.5d0*dot_product(Strain,Stress) 
                Self%StrainEnergyE(i,j) = energy
            end do
        end do
    end subroutine ProcessingResults
end module HN_Module
