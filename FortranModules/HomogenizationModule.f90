module HomogenizationModule
    implicit none
    ! ----------------------------------------------------------------- !
    !     this variable groups all the information of the structure     !
    ! ----------------------------------------------------------------- !
    type                                                        :: Structure        ! Name of the type of derived structure
        character(len=20)                                       :: DimAnalysis      ! 2D or 3D
        character(len=20)                                       :: AnalysisType     ! PlaneStress(2D), PlainStrain(2D), HookeLaw(3D)
        character(len=20)                                       :: ElementType      ! tri3 tri6, cuad4, cuad8, tetra4, tetra10, hexa8, hexa20
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
        integer, dimension(:,:), allocatable                    :: Node_Interaction
        integer, dimension(:,:), allocatable                    :: Element_Interaction 
        double precision, dimension(:), allocatable             :: value_KGlobal    ! Value of the global stiffness matrix at position i,j
        double precision, dimension(:,:), allocatable           :: value_FGlobal    ! Global load vector (for Sparse system)
        double precision, dimension(:,:), allocatable           :: value_UGlobal    ! Global displacement vector (for Sparse system)
        ! Results
        double precision, dimension(:,:), allocatable           :: UGlobal          ! Global displacement vector
        double precision, dimension(:,:), allocatable           :: UGlobal_0        ! Ideal Global displacement 
        double precision, dimension(:,:), allocatable           :: StrainEnergyE    ! Strain energy per element
        double precision, dimension(:,:), allocatable           :: StrainEnergyN    ! Strain energy per node
        double precision, dimension(:,:,:), allocatable         :: StrainE          ! Strain per element
        double precision, dimension(:,:,:), allocatable         :: StrainN          ! Strain per node
        double precision, dimension(:,:,:), allocatable         :: StressE          ! Stress per element
        double precision, dimension(:,:,:), allocatable         :: StressN          ! Stress per node
    contains
        procedure                                               :: SetDimAnalysis
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
    subroutine ReadDoublePrecisionFile(Path,nr,nc,Array)
        implicit none
        character(len=*), intent(in)                                    :: Path
        integer                                                         :: i, j, ios
        integer, intent(out)                                            :: nr
        integer, intent(in)                                             :: nc
        double precision, dimension(:,:), allocatable, intent(out)      :: Array
        !read file
        open(unit=1, file=Path, iostat=ios, status='old', Action='read')
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=1, fmt=*)
            read(unit=1, fmt=*) nr
            allocate(Array(nr,nc))
            read(unit=1, fmt=*)
            do i = 1, nr, 1
                read(unit=1, fmt=*) (Array(i,j), j=1, nc, 1)
            end do 
        close(1)
    end subroutine ReadDoublePrecisionFile
    ! listo
    subroutine ReadIntegerFile(Path,nr,nc,Array)
        implicit none
        character(len=*), intent(in)                                    :: Path
        integer                                                         :: i, j, ios
        integer, intent(out)                                            :: nr
        integer, intent(in)                                             :: nc
        integer, dimension(:,:), allocatable, intent(out)               :: Array
        !read file
        open(unit=1, file=Path, iostat=ios, status='old', Action='read')
            if ( ios /= 0 ) stop "Error opening file name"
            read(unit=1, fmt=*)
            read(unit=1, fmt=*) nr
            allocate(Array(nr,nc))
            read(unit=1, fmt=*)
            do i = 1, nr, 1
                read(unit=1, fmt=*) (Array(i,j), j=1, nc, 1)
            end do 
        close(1)
    end subroutine ReadIntegerFile
    ! listo
    function Area(Coordinates,Type) result(AreaAprox)
        character(len=*), intent(in)                                    :: Type
        double precision, dimension(:,:), allocatable, intent(in)       :: Coordinates
        double precision                                                :: AreaAprox
        ! internal variables
        double precision, dimension(:), allocatable                     :: vector1, vector2, vector3, vector4
        double precision, dimension(:), allocatable                     :: Av1, Av2
        allocate(vector1(3),vector2(3),vector3(3),vector4(3))
        allocate(Av1(3),Av2(3))
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
    ! listo
    subroutine VolumeAllElement(self)
        implicit none
        class(Structure), intent(inout)                                 :: self
        ! internal variables
        integer                                                         :: i
        double precision                                                :: Area
        double precision, dimension(:,:), allocatable                   :: Coordinates
        ! process
        allocate(self%VolumePerElement(self%Ne))
        self%VolumePerElement = 0.0d0
        ! Process
        do i = 1, self%Ne, 1
            if ((self%ElementType.eq.'tria3').or.(self%ElementType.eq.'tria6')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                Area = norm2(CrossProduct(Coordinates(3,:)-Coordinates(2,:), &
                                          Coordinates(1,:)-Coordinates(2,:)))/2
                self%VolumePerElement(i) = Area*self%Thickness
            elseif ((self%ElementType.eq.'quad4').or.(self%ElementType.eq.'quad8')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                Area = norm2(CrossProduct(Coordinates(2,:)-Coordinates(1,:), &
                                          Coordinates(4,:)-Coordinates(1,:)))/2 + &
                       norm2(CrossProduct(Coordinates(4,:)-Coordinates(3,:), &
                                          Coordinates(2,:)-Coordinates(3,:)))/2
                self%VolumePerElement(i) = Area*self%Thickness
            elseif ((self%ElementType.eq.'tetra4').or.(self%ElementType.eq.'tetra10')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                self%VolumePerElement(i) = abs(dot_product(CrossProduct(Coordinates(3,:)-Coordinates(1,:), &
                                                                         Coordinates(2,:)-Coordinates(1,:))/2, &
                                                                         Coordinates(4,:)-Coordinates(1,:))/2)
            elseif ((self%ElementType.eq.'hexa8').or.(self%ElementType.eq.'hexa20')) then
                Coordinates = self%Coordinates(self%ConnectivityN(i,:),:)
                self%VolumePerElement(i) = abs(dot_product(CrossProduct(Coordinates(1,:)-Coordinates(3,:), &
                                                                         Coordinates(4,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2) + &
                                            abs(dot_product(CrossProduct(Coordinates(2,:)-Coordinates(3,:), &
                                                                         Coordinates(1,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2) + &
                                            abs(dot_product(CrossProduct(Coordinates(2,:)-Coordinates(3,:), &
                                                                         Coordinates(6,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2) + &
                                            abs(dot_product(CrossProduct(Coordinates(6,:)-Coordinates(3,:), &
                                                                         Coordinates(7,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2) + &
                                            abs(dot_product(CrossProduct(Coordinates(7,:)-Coordinates(3,:), &
                                                                         Coordinates(8,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2) + &
                                            abs(dot_product(CrossProduct(Coordinates(8,:)-Coordinates(3,:), &
                                                                         Coordinates(4,:)-Coordinates(3,:))/2, &
                                                                         Coordinates(5,:)-Coordinates(3,:))/2)
            end if
        end do
        call WriteDoublePrecisionVector('DataStructure/AdditionalData/VolumePerElement.txt',self%VolumePerElement)
    end subroutine VolumeAllElement
    ! listo
    function CrossProduct(v1, v2) result(v3)
        double precision, dimension(3), intent(in)                      :: v1, v2
        double precision, dimension(3)                                  :: v3
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end function CrossProduct
    ! listo
    function Determinant(A) result(Det)
        implicit none
        double precision                                                :: Det
        double precision, intent(in)                                    :: A(:,:)
        if (size(A(1,:)).eq.2) then
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        elseif (size(A(1,:)).eq.3) then
            det = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)&
                 -A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)&
                 +A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)
        else
            write(*,*) 'ERROR, Determinant'
            call exit
        end if
    end function Determinant
    ! listo
    function Inverse(A) result(Ainv)
        ! this function use lapack library
        implicit none
        double precision,intent(in)                                     :: A(:,:)
        double precision                                                :: Ainv(size(A,1),size(A,2))
        ! internal variables
        double precision                                                :: work(size(A,1))
        integer                                                         :: n,info,ipiv(size(A,1))
        Ainv = A
        n = size(A,1)
        call DGETRF(n,n,Ainv,n,ipiv,info)
        !call SGETRF(n,n,Ainv,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        call DGETRI(n,Ainv,n,ipiv,work,n,info)
        !call SGETRI(n,Ainv,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function Inverse
    ! listo
    function Norm(vec) result(Ans)
        implicit none
        double precision, dimension(:), allocatable, intent(in)         :: vec
        double precision                                                :: sum_of_squares, Ans
        integer                                                         :: i
        sum_of_squares = 0.0d0
        do i = 1, size(vec)
            sum_of_squares = sum_of_squares + vec(i)**2
        end do
        Ans = sqrt(sum_of_squares)
    end function Norm
    ! listo
    subroutine Sort(vector,n)
        implicit none
        integer, intent(in)                                             :: n
        integer, allocatable, intent(inout)                             :: vector(:)
        ! internal variables
        integer                                                         :: i, j, temp
        integer, allocatable                                            :: dumb(:)
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
    ! escritura
    subroutine WriteIntegerArray(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                    :: path
        integer, dimension(:,:), allocatable, intent(in)                :: Array
        ! internal value
        integer                                                         :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array,1), 1
                write(unit=1, fmt=*) (Array(i,j),j=1,size(Array,2),1)
            end do    
        close(1)
    end subroutine WriteIntegerArray
    ! escritura
    subroutine WriteDoublePrecisionArray(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                    :: path
        double precision, dimension(:,:), allocatable, intent(in)       :: Array
        ! internal value
        integer                                                         :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array,1), 1
                write(unit=1, fmt=*) (Array(i,j),j=1,size(Array,2),1)
            end do    
        close(1)
    end subroutine WriteDoublePrecisionArray
    ! escritura
    subroutine WriteIntegerVector(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                    :: path
        integer, dimension(:), allocatable, intent(in)                  :: Array
        ! internal value
        integer                                                         :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array), 1
                write(unit=1, fmt=*) Array(i)
            end do    
        close(1)
    end subroutine WriteIntegerVector
    ! escritura
    subroutine WriteDoublePrecisionVector(path,Array)
        implicit none
        ! input
        character(len=*), intent(in)                                    :: path
        double precision, dimension(:), allocatable, intent(in)         :: Array
        ! internal value
        integer                                                         :: i,j,ios
        ! process
        open(unit=1, file=path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(Array), 1
                write(unit=1, fmt=*) Array(i)
            end do    
        close(1)
    end subroutine WriteDoublePrecisionVector
    ! ----------------------------------------------------------------- !
    !       subroutines to define the information required for FEA      !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine SetDimAnalysis(Self,DimAnalysis)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: DimAnalysis
        Self%DimAnalysis = DimAnalysis
        !write(unit=*, fmt=*) 'DimAnalysis   ',Self%DimAnalysis
    end subroutine
    ! listo
    subroutine SetAnalysisType(Self,AnalysisType)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: AnalysisType
        Self%AnalysisType = AnalysisType
        !write(unit=*, fmt=*) 'AnalysisType   ',Self%AnalysisType
    end subroutine SetAnalysisType
    ! listo
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
    ! listo
    subroutine SetThickness(Self,Thickness)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: Thickness
        Self%Thickness = Thickness
        !write(unit=*, fmt=*) 'Thickness   ', Self%Thickness
    end subroutine SetThickness
    ! listo
    subroutine SetYoungModulus(Self,YoungModulus)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: YoungModulus
        Self%YoungModulus = YoungModulus
        !write(unit=*, fmt=*) 'YoungModulus   ',Self%YoungModulus
    end subroutine SetYoungModulus
    ! listo
    subroutine SetPoissonModulus(Self,PoissonModulus)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(in)                                    :: PoissonModulus
        Self%PoissonModulus = PoissonModulus
        !write(unit=*, fmt=*) 'PoissonModulus   ',Self%PoissonModulus
    end subroutine SetPoissonModulus
    ! listo
    subroutine SetGaussAprox(Self,Gauss)
        implicit none
        class(Structure), intent(inout)                                 :: Self 
        integer, intent(in)                                             :: Gauss
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
        !write(unit=*, fmt=*) 'GaussAproximation   ',Self%QuadGauss
    end subroutine SetGaussAprox
    ! listo
    subroutine SetStrainNumericalH(Self,Strain)
        implicit none
        class(Structure), intent(inout)                                 :: Self 
        double precision, intent(in)                                    :: Strain
        Self%StrainH = Strain
        !write(unit=*, fmt=*) 'Strain for Homogenization   ',Self%StrainH
    end subroutine SetStrainNumericalH
    ! listo
    subroutine SetCoordinates(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        call ReadDoublePrecisionFile(Path,Self%N,3,Self%Coordinates)
    end subroutine SetCoordinates
    ! listo
    subroutine SetConnectivity(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: i, j, k, dim
        ! calculate ConnectivityD
        if ( Self%DimAnalysis.eq.'2D' ) then; dim = 2 ; end if
        if ( Self%DimAnalysis.eq.'3D' ) then; dim = 3 ; end if
        call ReadIntegerFile(Path,Self%Ne,Self%Npe,Self%ConnectivityN)
        allocate(Self%ConnectivityD(Self%Ne,Self%Npe*dim))
        do i = 1, Self%Ne, 1
            do j = 1, Self%Npe, 1
                do k = dim-1, 0, -1
                    Self%ConnectivityD(i,j*dim - k) = Self%ConnectivityN(i,j)*dim - k
                end do
            end do
        end do
        !call WriteIntegerArray('DataStructure/AdditionalData/ConnectivityN.txt',Self%ConnectivityN)
        !call WriteIntegerArray('DataStructure/AdditionalData/ConnectivityD.txt',Self%ConnectivityD)
    end subroutine SetConnectivity
    ! listo
    subroutine SetCornerNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: i,j,ios
        if (Self%DimAnalysis.eq.'2D') then; j=3; end if
        if (Self%DimAnalysis.eq.'3D') then; j=7; end if 
        open(unit=1, file=Path, iostat=ios, status="old", action="read")
            if ( ios /= 0 ) stop "Error corner nodes"
            allocate(Self%CornerNodes(j,2)); Self%CornerNodes = 0
            read(unit=1, fmt=*) ! master nodes
            read(unit=1, fmt=*) Self%CornerNodes(1,1)
            read(unit=1, fmt=*) ! slave nodes
            read(unit=1, fmt=*) (Self%CornerNodes(i,2),i=1,j,1)
        close(1)
        !call WriteIntegerArray('VerCornerNodes.txt',Self%CornerNodes)
    end subroutine SetCornerNodes
    ! listo
    subroutine SetEdgeNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: dim,i,j,ios
        if (Self%DimAnalysis.eq.'2D') then; dim = 2; end if
        if (Self%DimAnalysis.eq.'3D') then; dim = 3; end if
        select case (dim)
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
        !call WriteIntegerArray('VerEdgeNodes.txt',Self%EdgeNodes)
    end subroutine SetEdgeNodes
    ! listo
    subroutine SetFaceNodes(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: dim,i,j,ios
        if (Self%DimAnalysis.eq.'2D') then; dim = 2; end if
        if (Self%DimAnalysis.eq.'3D') then; dim = 3; end if
        select case (dim)
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
        !call WriteIntegerArray('VerFaceNodes.txt',Self%FaceNodes)
    end subroutine SetFaceNodes
    ! listo
    subroutine GetMasterSlaveNodes(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal values
        integer                                                         :: i,j,k,l,dim
        integer, dimension(:), allocatable                              :: Array
        integer, dimension(:,:), allocatable                            :: CondArray
        if (Self%DimAnalysis.eq.'2D') then; dim = 2; end if
        if (Self%DimAnalysis.eq.'3D') then; dim = 3; end if
        ! corner nodes
        Self%CornerNodes(:,1) = Self%CornerNodes(1,1)
        ! edge nodes
        select case (dim)
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
        select case (dim)
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
        allocate(Self%ConnectivityDPBC(Self%Ne,dim*Self%Npe))
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
                do k = dim-1, 0, -1
                    Self%ConnectivityDPBC(i,j*dim-k) = Self%ConnectivityNPBC(i,j)*dim - k
                end do
            end do
        end do
        !call WriteIntegerArray('DataStructure/AdditionalData/CornerNode.txt',Self%CornerNodes)
        !call WriteIntegerArray('DataStructure/AdditionalData/EdgeNode.txt',Self%EdgeNodes)
        !call WriteIntegerArray('DataStructure/AdditionalData/FaceNode.txt',Self%FaceNodes)
        !call WriteIntegerArray('DataStructure/AdditionalData/ConnectivityNPBC.txt',Self%ConnectivityNPBC)
        !call WriteIntegerArray('DataStructure/AdditionalData/ConnectivityDPBC.txt',Self%ConnectivityDPBC)
    end subroutine
    ! listo
    subroutine SetBondaryConditions(Self,Path)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        character(len=*), intent(in)                                    :: Path
        ! internal variables
        integer                                                         :: i,j,k,l,dim,Freedof,Fixeddof
        logical                                                         :: r1,r2,r3,r4
        integer, dimension(:), allocatable                              :: InPosition
        logical, dimension(:), allocatable                              :: InLogical
        if (Self%DimAnalysis.eq.'2D') then; dim = 2 ; end if
        if (Self%DimAnalysis.eq.'3D') then; dim = 3 ; end if
        ! ---- 1st part ----
        ! Get the degress of freedom
        call ReadIntegerFile(Path,Self%NBc,(dim+1),Self%BoundaryC)
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
        Fixeddof = l*dim
        Freedof = (Self%N-l)*dim
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
                do j = dim-1, 0, -1
                    Self%FixedD(k*dim-j) = i*dim - j
                end do
                k = k + 1
            else    ! Free
                do j = dim-1, 0, -1
                    Self%FreeD(l*dim-j) = i*dim - j
                end do
                l = l + 1
            end if
        end do
        ! printing
        !call WriteIntegerVector('DataStructure/AdditionalData/FreeDof.txt',Self%FreeD)
        !call WriteIntegerVector('DataStructure/AdditionalData/FixedDof.txt',Self%FixedD)
        ! ---- 2nd part ----
        ! get de element and DoF incidence for each DoF (only the free degres)
        ! (only considering connectivityPBC)
        j = size(Self%FreeD)
        ! consider a maximum of 20 elements per node
        allocate(Self%Element_Interaction(j,20))
        Self%Element_Interaction = 0
        ! considering the quad20 20*20 
        allocate(Self%Node_Interaction(j,400))
        Self%Node_Interaction = 0
        do i = 1, size(Self%FreeD), 1
            ! element (wtih PBC)
            InLogical = any(Self%ConnectivityDPBC.eq.Self%FreeD(i),2)
            InPosition = pack([(k,k=1,Self%Ne)],InLogical)
            j = count(InLogical)
            Self%Element_Interaction(i,1:j) = InPosition 
            ! other DoF (with PBC)
            InPosition = reshape(Self%ConnectivityDPBC(Self%Element_Interaction(i,1:j),:),[j*(Self%Npe)*dim])
            call sort(InPosition,j*Self%Npe*dim)
            j = size(InPosition)
            Self%Node_Interaction(i,1:j) = InPosition
        end do
        ! printing
        !call WriteIntegerArray('DataStructure/AdditionalData/FreeDofInteraction_Dof.txt',Self%Node_Interaction)
        !call WriteIntegerArray('DataStructure/AdditionalData/FreeDofInteraction_Element.txt',Self%Element_Interaction)
    end subroutine SetBondaryConditions
    ! ----------------------------------------------------------------- !
    !              subroutines that define the FEM procedure            !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine StrainDisplacement(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal values
        integer                                                         :: dim,nl,i,j,k
        double precision, dimension(:), allocatable                     :: Vector
        double precision, dimension(:,:,:), allocatable                 :: Array
        if (Self%DimAnalysis.eq.'2D') then
            dim = 2
            nl = 3
            allocate(Array(nl,dim,dim))
            Array = 0.0d0
            !xx
            Array(1,1,1) = Self%StrainH
            !yy
            Array(2,2,2) = Self%StrainH
            !xy
            Array(3,2,1) = 0.5d0*Self%StrainH
            Array(3,1,2) = 0.5d0*Self%StrainH
        elseif (Self%DimAnalysis.eq.'3D') then 
            dim = 3
            nl = 6
            allocate(Array(nl,dim,dim))
            Array = 0.0d0
            !xx
            Array(1,1,1) = Self%StrainH
            !yy
            Array(2,2,2) = Self%StrainH
            !zz
            Array(3,3,3) = Self%StrainH
            !xy
            Array(4,2,1) = 0.5d0*Self%StrainH
            Array(4,1,2) = 0.5d0*Self%StrainH
            !yz
            Array(5,3,2) = 0.5d0*Self%StrainH
            Array(5,2,3) = 0.5d0*Self%StrainH
            !xz
            Array(6,3,1) = 0.5d0*Self%StrainH
            Array(6,1,3) = 0.5d0*Self%StrainH
        end if
        allocate(Self%UGlobal_0(Self%N*dim,nl))
        Self%UGlobal_0 = 0
        do i = 1, Self%N, 1
            do j = 1, nl, 1
                Vector = matmul(Self%Coordinates(i,1:dim),Array(j,:,:))
                do k = dim-1,0,-1
                    Self%UGlobal_0(dim*i-k,j) = Vector(dim-k)
                end do
            end do
        end do
        deallocate(Array,Vector)
        !call WriteDoublePrecisionArray('DataStructure/AdditionalData/XoDisplacement.txt',Self%UGlobal_0)
    end subroutine StrainDisplacement
    ! listo
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
            DiffFunction(1,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0)+((n-1.0d0)*(z-3.0d0)*(e+n+z))/8.0d0, &
                            -((n-1.0d0)*(2.0d0*e+z-3.0d0)*(e+n+z))/8.0d0-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0) &
                            -2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            ((n+1.0d0)*(e+n+z)*(2.0d0*e+2.0d0*n+z-3.0d0))/8.0d0 &
                            +(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -((n+1.0d0)*(2.0d0*n+z-3.0d0)*(e+n+z))/8.0d0-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0), &
                            ((n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0)*(e+n+z-2.0d0))/8.0d0 &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -((n-1.0d0)*(2.0d0*n+z+1.0d0)*(e+n+z-2.0d0))/8.0d0-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0)+((n+1.0d0)*(z+1.0d0)*(e+n+z-2.0d0))/8.0d0, &
                            -((n+1.0d0)*(2.0d0*e+z+1.0d0)*(e+n+z-2.0d0))/8.0d0-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0) &
                            -2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0)+(e*(n-1.0d0)*(e+n+z))/2.0d0, &
                            -((n**2.0d0-1.0d0)*(e+n+z))/4.0d0-(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0)-(e*(n+1.0d0)*(e+n+z))/2.0d0, &
                            ((n**2.0d0-1.0d0)*(e+n+z))/4.0d0+(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            0.0d0, &
                            ((n**2.0d0-1.0d0)*(e+n+z-2.0d0))/4.0d0+(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0)+(e*(n+1.0d0)*(e+n+z-2.0d0))/2.0d0, &
                            -((n**2.0d0-1.0d0)*(e+n+z-2.0d0))/4.0d0-(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(((e+n+z-1.0d0)**2.0d0-1.0d0)*(n-1.0d0))/4.0d0 &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (((e+n+z-1.0d0)**2.0d0-1.0d0)*(n-1.0d0))/4.0d0 &
                            +(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(((e+n+z-1.0d0)**2.0d0-1.0d0)*(n+1.0d0))/4.0d0 &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (((e+n+z-1.0d0)**2.0d0-1.0d0)*(n+1.0d0))/4.0d0 &
                            +(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
            DiffFunction(2,:) = [(e/8.0d0-1.0d0/8.0d0)*(z-3.0d0)*(e+n+z)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(2.0d0*e+z-3.0d0)*(e+n+z) &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(e+n+z)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0) &
                            +2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -(e/8.0d0-1.0d0/8.0d0)*(2.0d0*n+z-3.0d0)*(e+n+z)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0) &
                            -2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            (e/8.0d0-1.0d0/8.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0)*(e+n+z-2.0d0) &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +2.0d0*(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(2.0d0*n+z+1.0d0)*(e+n+z-2.0d0) &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0) &
                            -2.0d0*(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(z+1.0d0)*(e+n+z-2.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(2.0d0*e+z+1.0d0)*(e+n+z-2.0d0) &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z)+(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0)-2.0d0*n*(e/4.0d0+1.0d0/4.0d0)*(e+n+z), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z)-(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0)+2.0d0*n*(e/4.0d0-1.0d0/4.0d0)*(e+n+z), &
                            0.0d0, &
                            (e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0)+2.0d0*n*(e/4.0d0+1.0d0/4.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(e+n+z-2.0d0)+(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0)-2.0d0*n*(e/4.0d0-1.0d0/4.0d0)*(e+n+z-2.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            +(e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*((e+n+z-1.0d0)**2.0d0-1.0d0) &
                            +(e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
            DiffFunction(3,:) = [(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(z-3.0d0)+(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+z-3.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+z-3.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*n+z-3.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z), &
                            (e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+z+1.0d0) &
                            +(e/8.0d0-1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(2.0d0*n+z+1.0d0)-(e/8.0d0+1.0d0/8.0d0)*(n-1.0d0)*(e+n+z-2.0d0), &
                            (e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(z+1.0d0)+(e/8.0d0+1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            -(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(2.0d0*e+z+1.0d0)-(e/8.0d0-1.0d0/8.0d0)*(n+1.0d0)*(e+n+z-2.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n-1.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            0.0d0, &
                            (e/4.0d0+1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            (e**2.0d0/4.0d0-1.0d0/4.0d0)*(n+1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n**2.0d0-1.0d0), &
                            -(e/4.0d0-1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0+1.0d0/4.0d0)*(n-1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            -(e/4.0d0+1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0), &
                            (e/4.0d0-1.0d0/4.0d0)*(n+1.0d0)*(2.0d0*e+2.0d0*n+2.0d0*z-2.0d0)]
        end if
    end subroutine DiffFormFunction
    ! listo
    subroutine ElasticityTensor(Self,ETensor)
        implicit none
        class(Structure), intent(inout)                                         :: Self    
        double precision, dimension(:,:), allocatable, intent(inout)            :: ETensor
        ! internal variables
        double precision                                                        :: E,V,Constant
        E = Self%YoungModulus
        V = Self%PoissonModulus    
        if (Self%AnalysisType.eq.'PlaneStress') then
            allocate(ETensor(3,3))
            Constant = E/(1.0d0-V**2d0)
            ETensor(1,:) = [1.0d0,V,0.0d0]
            ETensor(2,:) = [V,1.0d0,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0-V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'PlaneStrain') then
            allocate(ETensor(3,3))
            Constant = E/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,0.0d0]
            ETensor(3,:) = [0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        elseif (Self%AnalysisType.eq.'HookeLaw') then
            allocate(ETensor(6,6))
            Constant = E/((1.0d0+V)*(1.0d0-2.0d0*V))
            ETensor(1,:) = [1.0d0-V,V,V,0.0d0,0.0d0,0.0d0]
            ETensor(2,:) = [V,1.0d0-V,V,0.0d0,0.0d0,0.0d0]
            ETensor(3,:) = [V,V,1.0d0-V,0.0d0,0.0d0,0.0d0]
            ETensor(4,:) = [0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0,0.0d0,0.0d0]
            ETensor(5,:) = [0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0,0.0d0]
            ETensor(6,:) = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,(1.0d0-2.0d0*V)/2.0d0]
            ETensor = Constant*ETensor
        end if
    end subroutine ElasticityTensor
    ! listo
    subroutine SetFlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                         :: Self
        double precision, intent(inout)                                         :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)              :: DensityVector 
        ! internal variables
        integer                                                                 :: el,i,j,k,l,m,ios,dim
        double precision                                                        :: e,n,z,w1,w2,w3,DetJacobian
        double precision, dimension(:,:), allocatable                           :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable                           :: DiffN,DiffNXY,ElementCoordinates
        double precision, dimension(:,:), allocatable                           :: Epsilon_0
        write(unit=*, fmt=*) '(Local Load vector)'
        if (Self%DimAnalysis.eq.'2D') then
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
                        do k = 1, 2, 1
                            do l = 1, 2, 1
                                Jacobian(k,l) = dot_product(DiffN(k,:),ElementCoordinates(:,l))
                            end do
                        end do
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        ! Be
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            !xx
                            Be(1,2*k-1) = DiffNxy(1,k)
                            !Be(1,2*k) = 0.0d0
                            !yy
                            !Be(2,2*k-1) = 0.0d0
                            Be(2,2*k) = DiffNxy(2,k)
                            !xy
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! F Local
                        Self%FLocal(el,:,:) = Self%FLocal(el,:,:) + DetJacobian*(matmul(transpose(Be),matmul(D,Epsilon_0)))*&
                                                                    w1*w2*(Self%Thickness)
                        deallocate(DiffN)
                    end do
                end do
                deallocate(D)
            end do
        elseif(Self%DimAnalysis.eq.'3D') then
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
                            do l = 1, 3, 1
                                do m = 1, 3, 1
                                    Jacobian(l,m) = dot_product(DiffN(l,:),ElementCoordinates(:,m))
                                end do
                            end do
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            ! Be
                            Be = 0.0d0
                            do l = 1, size(DiffN,2), 1
                                !xx
                                Be(1,3*l-2) = DiffNxy(1,l)
                                !Be(1,3*l-1) = 0.0d0
                                !Be(1,3*l) = 0.0d0
                                !yy
                                !Be(2,3*l-2) = 0.0d0
                                Be(2,3*l-1) = DiffNxy(2,l)
                                !Be(2,3*l) = 0.0d0
                                !zz
                                !Be(3,3*l-2) = 0.0d0
                                !Be(3,3*l-1) = 0.0d0
                                Be(3,3*l) = DiffNxy(3,l)
                                !xy
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                !Be(4,3*l) = 0.0d0
                                !yz
                                !Be(5,3*l-2) = 0.0d0
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                !xz
                                Be(6,3*l-2) = DiffNxy(3,l)
                                !Be(6,3*l-1) = 0.0d0
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! K Local
                            Self%FLocal(el,:,:) = Self%FLocal(el,:,:) + DetJacobian*(matmul(transpose(Be),&
                                                                              matmul(D,Epsilon_0)))*w1*w2*w3
                            deallocate(DiffN)
                        end do
                    end do
                end do
                deallocate(D)
            end do
        end if
        close(1)        
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates,Epsilon_0)
        !assembling
        if (Self%DimAnalysis.eq.'2D') then
            dim = 2
            allocate(Self%value_FGlobal(dim*Self%N,3))
            Self%value_FGlobal = 0.0d0
        elseif (Self%DimAnalysis.eq.'3D') then
            dim = 3
            allocate(Self%value_FGlobal(dim*Self%N,6)) 
            Self%value_FGlobal = 0.0d0
        end if
        do i = 1, Self%Ne, 1
            Self%value_FGlobal(Self%ConnectivityDPBC(i,:),:) = Self%value_FGlobal(Self%ConnectivityDPBC(i,:),:) &
                                                               + Self%FLocal(i,:,:)
        end do
        call WriteDoublePrecisionArray('DataStructure/AdditionalData/StrainLoads.txt',Self%value_FGlobal)
    end subroutine SetFlocal
    ! listo
    subroutine SetKlocal(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)      :: DensityVector 
        ! internal variables
        integer                                                         :: el,i,j,k,l,m,ios
        integer                                                         :: row,col,loc1,loc2
        double precision                                                :: e,n,z,w1,w2,w3,DetJacobian,Kl,Kite
        double precision, dimension(:,:), allocatable                   :: Be,Jacobian,InvJacobian,D
        double precision, dimension(:,:), allocatable                   :: DiffN,DiffNXY,ElementCoordinates
        ! ---- 1st part ----
        ! in this section the local stiffness matrix (dense) of each of the elements is calculated using gauss 
        ! quadrature and each of them is stored in the array KLocal[element,[K_x,K_y]].
        !open(unit=1,file='DataStructure/AdditionalData/StiffnessMatrixArray.txt',iostat=ios,status="replace",action="write")
        ! this is just to check
        write(unit=*, fmt=*) '(Local Stiffness matrix)'
        if (Self%DimAnalysis.eq.'2D') then
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
                        do k = 1, 2, 1
                            do l = 1, 2, 1
                                Jacobian(k,l) = dot_product(DiffN(k,:),ElementCoordinates(:,l))
                            end do
                        end do
                        InvJacobian = Inverse(Jacobian)
                        DetJacobian = Determinant(Jacobian)
                        DiffNXY = matmul(InvJacobian,DiffN)
                        ! Be
                        Be = 0.0d0
                        do k = 1, size(DiffN,2), 1
                            Be(1,2*k-1) = DiffNxy(1,k)
                            !Be(1,2*k) = 0.0d0
                            !Be(2,2*k-1) = 0.0d0
                            Be(2,2*k) = DiffNxy(2,k)
                            Be(3,2*k-1) = DiffNxy(2,k)
                            Be(3,2*k) = DiffNxy(1,k)
                        end do
                        ! Be Local
                        Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + DetJacobian*Be*w1*w2
                        ! K Local
                        Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + DetJacobian*(matmul(transpose(Be),matmul(D,Be)))*&
                                                                    w1*w2*(Self%Thickness)
                        deallocate(DiffN)
                    end do
                end do
                !write(unit=1, fmt=*) reshape(Self%KLocal(el,:,:),[(Self%Npe*2)*(Self%Npe*2)])
                deallocate(D)
            end do
            !write(unit=*, fmt=*) 'salida de tensor de elasticidad', D(:,:)
        elseif(Self%DimAnalysis.eq.'3D') then
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
                            do l = 1, 3, 1
                                do m = 1, 3, 1
                                    Jacobian(l,m) = dot_product(DiffN(l,:),ElementCoordinates(:,m))
                                end do
                            end do
                            InvJacobian = Inverse(Jacobian)
                            DetJacobian = Determinant(Jacobian)
                            DiffNXY = matmul(InvJacobian,DiffN)
                            ! Be
                            Be = 0.0d0
                            do l = 1, size(DiffN,2), 1
                                Be(1,3*l-2) = DiffNxy(1,l)
                                !Be(1,3*l-1) = 0.0d0
                                !Be(1,3*l) = 0.0d0
                                !Be(2,3*l-2) = 0.0d0
                                Be(2,3*l-1) = DiffNxy(2,l)
                                !Be(2,3*l) = 0.0d0
                                !Be(3,3*l-2) = 0.0d0
                                !Be(3,3*l-1) = 0.0d0
                                Be(3,3*l) = DiffNxy(3,l)
                                Be(4,3*l-2) = DiffNxy(2,l)
                                Be(4,3*l-1) = DiffNxy(1,l)
                                !Be(4,3*l) = 0.0d0
                                !Be(5,3*l-2) = 0.0d0
                                Be(5,3*l-1) = DiffNxy(3,l)
                                Be(5,3*l) = DiffNxy(2,l)
                                Be(6,3*l-2) = DiffNxy(3,l)
                                !Be(6,3*l-1) = 0.0d0
                                Be(6,3*l) = DiffNxy(1,l)
                            end do
                            ! Be Local
                            Self%BLocal(el,:,:) = Self%BLocal(el,:,:) + DetJacobian*Be*w1*w2*w3
                            ! K Local
                            Self%KLocal(el,:,:) = Self%KLocal(el,:,:) + DetJacobian*(matmul(transpose(Be),matmul(D,Be)))*&
                                                                            w1*w2*w3
                            deallocate(DiffN)
                        end do
                    end do
                end do
                !write(unit=1, fmt=*) reshape(Self%KLocal(el,:,:),[Self%Npe*3**2])
                deallocate(D)
            end do
        end if
        close(1)        
        deallocate(Be,Jacobian,InvJacobian,DiffNXY,ElementCoordinates)
        !write(unit=*, fmt=*) 'max value stiffness matrix',maxval(Self%KLocal)
        !write(unit=*, fmt=*) 'min value stiffness matrix',minval(self%KLocal)
        
        ! ---- 2nd part ----
        ! in this section the sparse stiffness matrix (free degrees only) is assembled for further solution
        write(unit=*, fmt=*) '(Sparse Stiffness matrix)'
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseIndex.txt',iostat=ios,status="replace",action="write")
        open(unit=2,file='DataStructure/AdditionalData/SparseSystem/KSparseValue.txt',iostat=ios,status="replace",action="write")
        Kl = 0.0d0
        m = 0
        do j = 1, size(Self%FreeD), 1                                   ! column
            col = Self%FreeD(j)
            do i = 1, size(Self%Node_Interaction,2), 1                  ! row
                row = Self%Node_Interaction(j,i)
                if (row.eq.0) exit                                      ! exit (doenst interact)
                if (any(Self%FixedD.eq.row)) cycle                      ! discard (fixed DoF)
                if (col.gt.row) cycle                                   ! discard (bottom diagonal only)
                m = m + 1
                do k = 1,size(Self%Element_Interaction,2), 1            ! element
                    el = Self%Element_Interaction(j,k)
                    if (el.eq.0) exit                                   ! exit (no more elements)
                    loc1 = findloc(Self%ConnectivityDPBC(el,:),row,1)
                    if (loc1.eq.0) cycle                                ! discard (doesnt exist in element)
                    loc2 = findloc(Self%ConnectivityDPBC(el,:),col,1)
                    Kite = Self%KLocal(el,loc1,loc2)
                    ! if (Kite.eq.0) then; cycle; end if
                    Kl = Kl + Kite
                end do
                ! adding the elements to the vectors
                write(unit=1, fmt=*) findloc(Self%FreeD,row,dim=1), findloc(Self%FreeD,col,dim=1)
                write(unit=2, fmt=*) Kl
                Kl = 0.0d0                                             ! reset
            end do
        end do
        close(2)
        close(1)
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/Resume.txt',iostat=ios,status="replace",action="write")
            write(unit=1,fmt=*) size(Self%FreeD), m 
        close(1)
    end subroutine SetKlocal
    ! listo
    subroutine UploadStructure(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                             :: Self
        double precision, intent(inout)                                         :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)              :: DensityVector 
        write(unit=*, fmt=*) 'Uploading Structure'
        ! 1. Stiffnes matrix
        write(unit=*, fmt=*) '1. Stiffness matrix'
        call SetKlocal(Self,DensityVector,PenalFactor)
        ! 1. load vector
        write(unit=*, fmt=*) '2. Load vector'
        call SetFlocal(Self,DensityVector,PenalFactor)
        Self%value_FGlobal = Self%value_FGlobal(Self%FreeD,:)
        !call WriteDoublePrecisionArray('DataStructure/AdditionalData/SparseSystem/LoadVectorSparse.txt',Self%value_FGlobal)
    end subroutine UploadStructure
    ! listo
    subroutine SparseSolver(Self)
        use hsl_ma87_double
        use hsl_mc68_double
        use hsl_mc69_double
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! HSL-MA87 variables
        !integer, parameter                                              :: wp = kind(0d0)
        type(mc68_control)                                              :: control68
        type(mc68_info)                                                 :: info68
        type(ma87_keep)                                                 :: keep
        type(ma87_control)                                              :: control
        type(ma87_info)                                                 :: info
        integer                                                         :: dim
        !---------------------------
        !integer                                                     :: i
        !integer, dimension(:), allocatable                          :: index
        !------------------------
        integer                                                         :: n, ne, nrhs, lmap, flag
        integer, dimension(:), allocatable                              :: crow, ccol, ptr, row, order, map
        double precision, dimension(:), allocatable                     :: cval, val
        double precision, dimension(:,:), allocatable                   :: x
        ! Organizing variables
        if (Self%DimAnalysis.eq.'2D') then; nrhs = 3; end if
        if (Self%DimAnalysis.eq.'3D') then; nrhs = 6; end if
        n = size(self%value_FGlobal(:,1));    ne = size(Self%index_j_KGlobal);
        allocate(crow(ne));                 crow = self%index_i_KGlobal;
        allocate(ccol(ne));                 ccol = self%index_j_KGlobal;
        allocate(cval(ne));                 cval = self%value_KGlobal;
        allocate(x(n,nrhs));                   x = self%value_FGlobal;
        ! -----------------------------------------------------------
        ! cheque para saber si es positivo definido
        !index = pack([(i,i=1,ne)],crow.eq.ccol)
        !write(unit=*, fmt=*) 'el sistema es positivo definido?', all(cval(index).gt.0.0)
        ! -----------------------------------------------------------
        ! Convert to HSL standard format
        allocate(ptr(n+1))
        call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, n, n, ne, crow, ccol, &
            ptr, row, flag, val_in=cval, val_out=val, lmap=lmap, map=map)
        call stop_on_bad_flag("mc69_coord_convert", flag)
        ! Call mc68 to find a fill reducing ordering (1=AMD)
        allocate(order(n))
        call mc68_order(1, n, ptr, row, order, control68, info68)
        call stop_on_bad_flag("mc68_order", info68%flag)
        ! Analyse
        call ma87_analyse(n, ptr, row, order, keep, control, info)
        call stop_on_bad_flag("analyse", info%flag)
        ! Factor
        call ma87_factor_solve(n, ptr, row, val, order, keep, control, info, &
        nrhs, n, x)
        call stop_on_bad_flag("factor_solve", info%flag)
        Self%value_UGlobal = x
        !call WriteDoublePrecisionArray('DataStructure/AdditionalData/SparseSystemDisplacement.txt',Self%value_UGlobal)
        ! Finalize
        call ma87_finalise(keep, control)
    end subroutine SparseSolver
    ! listo
    subroutine stop_on_bad_flag(context, flag)
        character(len=*), intent(in)                                    :: context
        integer, intent(in)                                             :: flag
        if(flag.eq.0) return
        write(*,*) "Failure during ", context, " with flag = ", flag
        stop
    end subroutine stop_on_bad_flag
    ! listo
    subroutine HSL87Solver(Self)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal variables
        double precision, dimension(:,:), allocatable                   :: Deformation
        integer                                                         :: i,j,ind1,ind2,dim,n,ne,ios,ns
        ! procedure
        if (Self%DimAnalysis.eq.'2D') then;     dim=2; ns=3
        elseif (Self%DimAnalysis.eq.'3D') then; dim=3; ns=6
        end if
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/Resume.txt',iostat=ios, status="old", action="read")
            read(1,*) n, ne
        close(1)
        allocate(Self%index_i_KGlobal(ne));              Self%index_i_KGlobal = 0;
        allocate(Self%index_j_KGlobal(ne));              Self%index_j_KGlobal = 0;
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseIndex.txt',iostat=ios, status="old", action="read")
            do i = 1, ne, 1
                read(1,*) Self%index_i_KGlobal(i), Self%index_j_KGlobal(i)
            end do
        close(1)
        allocate(Self%value_KGlobal(ne));                  Self%value_KGlobal = 0.0d0;
        open(unit=1,file='DataStructure/AdditionalData/SparseSystem/KSparseValue.txt',iostat=ios, status="old", action="read")
            do i = 1, ne, 1
                read(1,*) Self%value_KGlobal(i)
            end do
        close(1)
        if (allocated(Self%UGlobal)) then
            Self%UGlobal = 0.0d0;
            Self%value_UGlobal = 0.0d0;
        else
             allocate(Self%UGlobal(Self%N*dim,ns));              Self%UGlobal = 0.0d0;
            allocate(Self%value_UGlobal(n,ns));            Self%value_UGlobal = 0.0d0;
        end if
        ! Sparse solver HSL-Ma87
        write(unit=*, fmt=*) '3. Applying Solver'
        call SparseSolver(Self)
        Self%UGlobal(Self%FreeD,:) = Self%value_UGlobal
        ! ---- extension to the slaves PBC ---- !
        ! corner
        do i = 1, size(Self%CornerNodes(:,1)), 1
            do j = dim-1,0,-1
                ind1 = Self%CornerNodes(i,1)*dim-j    !master
                ind2 = Self%CornerNodes(i,2)*dim-j    !slave
                if ((Self%CornerNodes(i,1).eq.0).or.(Self%CornerNodes(i,2).eq.0)) cycle
                Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
            end do
        end do
        ! edge
        do i = 1, size(Self%EdgeNodes(:,1)), 1
            do j = dim-1,0,-1
                ind1 = Self%EdgeNodes(i,1)*dim-j      !master
                ind2 = Self%EdgeNodes(i,2)*dim-j      !slave
                if ((Self%EdgeNodes(i,1).eq.0).or.(Self%EdgeNodes(i,2).eq.0)) cycle
                Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
            end do
        end do
        ! face
        if (dim.eq.3) then
            do i = 1, size(Self%FaceNodes(:,1)), 1
                do j = dim-1,0,-1
                    ind1 = Self%FaceNodes(i,1)*dim-j  !master
                    ind2 = Self%FaceNodes(i,2)*dim-j  !slave
                    if ((Self%FaceNodes(i,1).eq.0).or.(Self%FaceNodes(i,2).eq.0)) cycle
                    Self%UGlobal(ind2,:) = Self%UGlobal(ind1,:)
                end do
            end do
        end if
        call WriteDoublePrecisionArray('NumericalResults/Deformation.txt',Self%UGlobal)
        ! deallocating the sparse system, it is no longer necessary
        deallocate(Self%value_FGlobal,Self%value_KGlobal,Self%index_i_KGlobal,Self%index_j_KGlobal)
    end subroutine HSL87Solver
    ! listo
    subroutine NumericalHomogenization(Self,DensityVector,PenalFactor)
        implicit none
        class(Structure), intent(inout)                                 :: Self
        double precision, intent(inout)                                 :: PenalFactor
        double precision, dimension(:), allocatable, intent(inout)      :: DensityVector 
        ! internal values
        integer                                                         :: i,j,nl
        double precision                                                :: Vtotal
        double precision, dimension(:,:), allocatable                   :: D,Epsilon,Dtotal
        ! this section calculates the volume
        if (allocated(Self%VolumePerElement).eqv..false.) then
            call VolumeAllElement(Self)
        end if
        Vtotal = sum(Self%VolumePerElement)
        !Vtotal = 1.0
        ! this section calculates the Ideal displacement
        if (allocated(Self%UGlobal_0).eqv..false.) then
            Call StrainDisplacement(Self)
        end if
        ! the displacement is obtained
        Self%UGlobal = Self%UGlobal_0 - Self%UGlobal
        ! the homogenised tensor is calculated (per element)
        if (Self%DimAnalysis.eq.'2D') nl = 3
        if (Self%DimAnalysis.eq.'3D') nl = 6
        if (allocated(Self%DTensor)) then
            Self%DTensor = 0.0d0
        else            
            allocate(Self%DTensor(Self%Ne,nl,nl))
            Self%DTensor = 0.0d0
        end if
        allocate(Dtotal(nl,nl)) ;                     Dtotal = 0.0d0
        allocate(Epsilon(nl,nl)) ;                   Epsilon = 0.0d0
        do i = 1, Self%Ne, 1
            call ElasticityTensor(Self,D)
            D = (DensityVector(i)**PenalFactor)*D
            do j = 1, size(D,1), 1
                Epsilon(:,j) = matmul(Self%BLocal(i,:,:),Self%UGlobal(Self%ConnectivityD(i,:),j))
            end do
            Self%DTensor(i,:,:) = 1.0d0/Vtotal*(matmul(transpose(Epsilon),matmul(D,Epsilon)))
            deallocate(D)
        end do
        ! total
        do i = 1, Self%Ne, 1
            Dtotal = Dtotal + Self%DTensor(i,:,:)
        end do
        if (Self%DimAnalysis.eq.'2D') then
            Dtotal(1:2,3)=0.0d0
            Dtotal(3,1:2)=0.0d0
        elseif (Self%DimAnalysis.eq.'3D') then
            Dtotal(1:3,4:6)=0.0d0
            Dtotal(4,5:6)=0.0d0
            Dtotal(5,6)=0.0d0
            Dtotal(4:6,1:3)=0.0d0
            Dtotal(6:4,5)=0.0d0
            Dtotal(6,5)=0.0d0
        end if
        call WriteDoublePrecisionArray('NumericalResults/HomogenizedTensor.txt',Dtotal)
        deallocate(Epsilon,Dtotal)
    end subroutine NumericalHomogenization
    ! listo
    subroutine ProcessingResults(Self)
        ! this subroutine calculates stresses per element, per node, deformations and strain energy.
        implicit none
        class(Structure), intent(inout)                                 :: Self
        ! internal variables
        integer                                                         :: i,j,nl
        integer, dimension(:), allocatable                              :: index
        double precision                                                :: energy
        double precision, dimension(:), allocatable                     :: Strain, Stress   ! (load case, elm/node, direction)
        double precision, dimension(:,:), allocatable                   :: D
        call ElasticityTensor(Self,D)
        write(unit=*, fmt=*) '4. Processing Results'
        ! allocating
        if (Self%DimAnalysis.eq.'2D') then
            nl = 3
            allocate(Strain(nl),Stress(nl))
            Strain = 0.0d0; Stress = 0.0d0;
            if (allocated(Self%StrainE).and.allocated(Self%StrainN)) then
                Self%StrainE = 0.0d0; Self%StrainN = 0.0d0
            else
                allocate(Self%StrainE(nl,Self%Ne,nl),Self%StrainN(nl,Self%N,nl))        
                Self%StrainE = 0.0d0; Self%StrainN = 0.0d0
            end if
            if (allocated(Self%StressE).and.allocated(Self%StressN)) then
                Self%StressE = 0.0d0; Self%StressN = 0.0d0
            else
                allocate(Self%StressE(nl,Self%Ne,nl),Self%StressN(nl,Self%N,nl))
                Self%StressE = 0.0d0; Self%StressN = 0.0d0
            end if
        elseif (Self%DimAnalysis.eq.'3D') then
            nl = 6
            allocate(Strain(nl),Stress(nl))
            Strain = 0.0d0; Stress = 0.0d0;
            if (allocated(Self%StrainE).and.allocated(Self%StrainN)) then
                Self%StrainE = 0.0d0; Self%StrainN = 0.0d0
            else
                allocate(Self%StrainE(nl,Self%Ne,nl),Self%StrainN(nl,Self%N,nl))        
                Self%StrainE = 0.0d0; Self%StrainN = 0.0d0
            end if
            if (allocated(Self%StressE).and.allocated(Self%StressN)) then
                Self%StressE = 0.0d0; Self%StressN = 0.0d0
            else
                allocate(Self%StressE(nl,Self%Ne,nl),Self%StressN(nl,Self%N,nl))
                Self%StressE = 0.0d0; Self%StressN = 0.0d0
            end if
        end if
        if (allocated(Self%StrainEnergyE).and.Allocated(Self%StrainEnergyN)) then
            Self%StrainEnergyE = 0.0d0; Self%StrainEnergyN = 0.0d0
        else
            allocate(Self%StrainEnergyE(nl,Self%Ne),Self%StrainEnergyN(nl,Self%N)); 
            Self%StrainEnergyE = 0.0d0; Self%StrainEnergyN = 0.0d0
        end if
        ! star processing
        ! (calculation per element)
        do i = 1, Self%Ne, 1
            do j = 1, nl, 1
                ! 1. Strain     Strain_e = Be*Ue        
                ! exx - eyy - ezz(3D) - exy - eyz(3D) - exz(3D)
                Strain = matmul(Self%BLocal(i,:,:),Self%UGlobal(Self%ConnectivityD(i,:),j))
                Self%StrainE(j,i,:) = Strain
                ! 2. Stress      Sigma_e = D*Strain_e   
                ! Sxx - Syy - Szz(3D) - Sxy - Syz(3D) - Sxz(3D)
                Stress = matmul(D,Strain)
                Self%StressE(j,i,:) = Stress
                ! 3. Energy    Ee = 0.5*Sigma_e*Strain_e
                energy = 0.5d0*dot_product(Strain,Stress) 
                Self%StrainEnergyE(j,i) = energy
            end do
        end do
        ! (node-weighted calculation)
        do i = 1, Self%N, 1
            do j = 1, nl, 1
                ! locating the nodes 
                index = pack([(j,j=1,Self%Ne)],any(Self%ConnectivityN.eq.i,2))
                ! 1. Strain
                Strain = sum(Self%StrainE(j,index,:),1)
                Self%StrainN(j,i,:) = Strain/size(index)
                ! 2. Stress
                Stress = sum(Self%StressE(j,index,:),1)
                Self%StressN(j,i,:) = Stress/size(index)
                ! 3. Energy
                Energy = sum(Self%StrainEnergyE(j,index))
                Self%StrainEnergyN(j,i) = Energy/size(index)
            end do
        end do
    end subroutine ProcessingResults
end module HomogenizationModule
