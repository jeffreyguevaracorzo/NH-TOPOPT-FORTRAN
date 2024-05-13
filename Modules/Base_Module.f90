module Base_Module
    implicit none
    ! Lapack and BLAS are needed in some routines/functions
    ! -------------------------------------------------------------------------------------------------!
    ! Matematical constants
    double precision, parameter                             :: Pi = 3.14159265358979323846d0
    double precision, parameter                             :: Euler = 2.71828182845904523536d0
    double precision, parameter                             :: Aureo = 1.618033988749895d0
    double precision, parameter                             :: EulerMascheroni = 0.57721566490153286060d0
    double precision, parameter                             :: Feigenbaum = 4.66920160910299067185d0
    double precision, parameter                             :: Catalan = 0.91596559417721901505d0
    double precision, parameter                             :: Apery = 1.20205690315959428539d0
    ! Physical/Fundamental constants
    double precision, parameter                             :: LightSpeed = 299792458.0d0
    double precision, parameter                             :: GravConst = 6.674E-11
    double precision, parameter                             :: PlankConst = 6.626E-34
    double precision, parameter                             :: BoltzmannConst = 1.381E-23
    double precision, parameter                             :: AvogadroConst = 6.02214076E-23
    double precision, parameter                             :: FaradayConst = 96485.33212d0
    ! ----------------- ALL INTERFACES -----------------
    ! Algebra
    interface Inverse
        module procedure InverseReal
        module procedure InverseDP
    end interface Inverse
    interface TransposeMatrix
        module procedure TransposeMatrixInteger
        module procedure TransposeMatrixReal
        module procedure TransposeMatrixDP
    end interface TransposeMatrix
    interface Determinant
        module procedure DeterminantReal
        module procedure DeterminantDP
    end interface Determinant
    interface CrossProduct
        module procedure CrossProductInteger
        module procedure CrossProductReal
        module procedure CrossProductDP
    end interface CrossProduct
    interface DotProduct
        module procedure DotProductInteger
        module procedure DotProductReal
        module procedure DotProductDP
    end interface DotProduct
    interface Norm
        module procedure NormRealVector
        module procedure NormDPVector
        module procedure NormRealMatrix
        module procedure NormDPMatrix
    end interface Norm 
    interface UnitVector
        module procedure UnitVectorReal
        module procedure UnitVectorDP
    end interface UnitVector
    interface Add2Vector
        module procedure AddInteger2Vector
        module procedure AddReal2Vector
        module procedure AddDP2Vector
    end interface Add2Vector
    Interface AddVector2Matrix
        module procedure AddIntegerVector2Matrix
        module procedure AddRealVector2Matrix
        module procedure AddDPVector2Matrix
    end interface AddVector2Matrix
    Interface AddMatrix2Matrix
        module procedure AddIntegerMatrix2Matrix
        module procedure AddRealMatrix2Matrix
        module procedure AddDPMatrix2Matrix
    end interface AddMatrix2Matrix
    ! Additional base functions
    interface ConsolePrinting
        module procedure PrintIntegerVectorInConsole
        module procedure PrintRealVectorInConsole
        module procedure PrintDPVectorInConsole
        module procedure PrintIntegerMatrixInConsole
        module procedure PrintRealMatrixInConsole
        module procedure PrintDPMatrixInConsole
    end interface ConsolePrinting
    interface FilePrinting
        module procedure PrintIntegerVectorInFile
        module procedure PrintRealVectorInFile
        module procedure PrintDPVectorInFile
        module procedure PrintIntegerMatrixInFile
        module procedure PrintRealMatrixInFile
        module procedure PrintDPMatrixInFile
    end interface FilePrinting
    interface ScatterPlot
        module procedure ScatterPlotReal
        module procedure ScatterPlotDP
    end interface ScatterPlot
    interface LinePlot
        module procedure LinePlotReal
        module procedure LinePlotDP
    end interface LinePlot

    contains
    ! ----------------- ALL FUNCTIONS AND SUBROUTINES -----------------
    ! 1. Algebra/Vector operations
    ! 1.1. Inverse
    function InverseReal(Matrix) result(InverseMatrix)
        ! this function use lapack library
        implicit none
        real, dimension(:,:), allocatable, intent(in)       :: Matrix
        real, dimension(:,:), allocatable                   :: InverseMatrix
        real, dimension(:), allocatable                     :: work
        integer                                             :: n,info
        integer, dimension(:), allocatable                  :: ipiv
        allocate(InverseMatrix(size(Matrix,1),size(Matrix,2)))
        allocate(work(size(Matrix,1)))
        allocate(ipiv(size(Matrix,1)))
        InverseMatrix = Matrix
        n = size(Matrix,1)
        call SGETRF(n,n,InverseMatrix,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        call SGETRI(n,InverseMatrix,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function InverseReal
    function InverseDP(Matrix) result(InverseMatrix)
        ! this function use lapack library
        implicit none
        double precision, dimension(:,:), allocatable, intent(in)       :: Matrix
        double precision, dimension(:,:), allocatable                   :: InverseMatrix
        double precision, dimension(:), allocatable                     :: work
        integer                                                         :: n,info
        integer, dimension(:), allocatable                              :: ipiv
        allocate(InverseMatrix(size(Matrix,1),size(Matrix,2)))
        allocate(work(size(Matrix,1)))
        allocate(ipiv(size(Matrix,1)))
        InverseMatrix = Matrix
        n = size(Matrix,1)
        call DGETRF(n,n,InverseMatrix,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        call DGETRI(n,InverseMatrix,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
    end function InverseDP
    ! 1.2. Transpose
    function TransposeMatrixInteger(Matrix) result(TransMatrix)
        implicit none
        integer                                             :: i,j,k
        integer, dimension(:,:), allocatable, intent(in)    :: Matrix
        integer, dimension(:,:), allocatable                :: TransMatrix
        j = size(Matrix,1)   ! columnas
        k = Size(Matrix,2)   ! filas
        allocate(TransMatrix(j,k))
        do i = 1, j, 1
            TransMatrix(j,:) = TransMatrix(:,j)
        end do
    end function TransposeMatrixInteger
    function TransposeMatrixReal(Matrix) result(TransMatrix)
        implicit none
        integer                                          :: i,j,k
        Real, dimension(:,:), allocatable, intent(in)    :: Matrix
        Real, dimension(:,:), allocatable                :: TransMatrix
        j = size(Matrix,1)   ! columnas
        k = Size(Matrix,2)   ! filas
        allocate(TransMatrix(j,k))
        do i = 1, j, 1
            TransMatrix(j,:) = TransMatrix(:,j)
        end do
    end function TransposeMatrixReal
    function TransposeMatrixDP(Matrix) result(TransMatrix)
        implicit none
        integer                                                      :: i,j,k
        double precision, dimension(:,:), allocatable, intent(in)    :: Matrix
        double precision, dimension(:,:), allocatable                :: TransMatrix
        j = size(Matrix,1)   ! columnas
        k = Size(Matrix,2)   ! filas
        allocate(TransMatrix(j,k))
        do i = 1, j, 1
            TransMatrix(j,:) = TransMatrix(:,j)
        end do
    end function TransposeMatrixDP
    ! 1.3. Determinant
    function DeterminantReal(Matrix) result(DetValue)
        ! this function use lapack library
        implicit none
        real, dimension(:,:), allocatable, intent(in)       :: Matrix
        integer                                             :: i
        integer                                             :: n,info
        real                                                :: DetValue
        real, dimension(:), allocatable                     :: work
        integer, dimension(:), allocatable                  :: ipiv
        allocate(work(size(Matrix,1)))
        allocate(ipiv(size(Matrix,1)))
        n = size(Matrix,1)
        call SGETRF(n,n,Matrix,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        DetValue = 1.0
        do i = 1, n, 1
            DetValue = DetValue*Matrix(i,i)
            if (ipiv(i) .ne. i) DetValue = (-1.0)*DetValue
        end do
    end function DeterminantReal
    function DeterminantDP(Matrix) result(DetValue)
        ! this function use lapack library
        implicit none
        double precision, dimension(:,:), allocatable, intent(in)   :: Matrix
        integer                                                     :: i
        integer                                                     :: n,info
        double precision                                            :: DetValue
        double precision, dimension(:), allocatable                 :: work
        integer, dimension(:), allocatable                          :: ipiv
        allocate(work(size(Matrix,1)))
        allocate(ipiv(size(Matrix,1)))
        n = size(Matrix,1)
        call DGETRF(n,n,Matrix,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        DetValue = 1.0
        do i = 1, n, 1
            DetValue = DetValue*Matrix(i,i)
            if (ipiv(i) .ne. i) DetValue = (-1.0d0)*DetValue
        end do
    end function DeterminantDP
    ! 1.4. Cross Product (vectors)
    function CrossProductInteger(Vector1,Vector2) result(CPVector)
        implicit none
        integer, dimension(:), allocatable, intent(in)         :: Vector1
        integer, dimension(:), allocatable, intent(in)         :: Vector2
        integer, dimension(:), allocatable                     :: CPVector
        allocate(CPVector(3))
        CPVector(1) = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
        CPVector(2) = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
        CPVector(3) = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)
    end function CrossProductInteger
    function CrossProductReal(Vector1,Vector2) result(CPVector)
        implicit none
        real, dimension(:), allocatable, intent(in)         :: Vector1
        real, dimension(:), allocatable, intent(in)         :: Vector2
        real, dimension(:), allocatable                     :: CPVector
        allocate(CPVector(3))
        CPVector(1) = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
        CPVector(2) = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
        CPVector(3) = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)
    end function CrossProductReal
    function CrossProductDP(Vector1,Vector2) result(CPVector)
        implicit none
        double precision, dimension(:), allocatable, intent(in)         :: Vector1
        double precision, dimension(:), allocatable, intent(in)         :: Vector2
        double precision, dimension(:), allocatable                     :: CPVector
        allocate(CPVector(3))
        CPVector(1) = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
        CPVector(2) = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
        CPVector(3) = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)
    end function CrossProductDP
    ! 1.5. Dot Product (vectors)
    function DotProductInteger(Vector1,Vector2) result(DPVector)
        implicit none
        integer, dimension(:), allocatable, intent(in)  :: Vector1
        integer, dimension(:), allocatable, intent(in)  :: Vector2
        integer, dimension(:), allocatable              :: DPVector
        DPVector = sum(Vector1*Vector2)
    end function DotProductInteger
    function DotProductReal(Vector1,Vector2) result(DPVector)
        implicit none
        real, dimension(:), allocatable, intent(in)  :: Vector1
        real, dimension(:), allocatable, intent(in)  :: Vector2
        real, dimension(:), allocatable              :: DPVector
        DPVector = sum(Vector1*Vector2)
    end function DotProductReal
    function DotProductDP(Vector1,Vector2) result(DPVector)
        implicit none
        double precision, dimension(:), allocatable, intent(in)  :: Vector1
        double precision, dimension(:), allocatable, intent(in)  :: Vector2
        double precision, dimension(:), allocatable              :: DPVector
        DPVector = sum(Vector1*Vector2)
    end function DotProductDP
    ! 1.6. Norm
    function NormRealVector(Vector) result(NormValue)
        implicit none
        real, dimension(:), allocatable, intent(in)         :: Vector
        real                                                :: NormValue
        NormValue = sqrt(sum(Vector**2))
    end function NormRealVector
    function NormDPVector(Vector) result(NormValue)
        implicit none
        double precision, dimension(:), allocatable, intent(in)         :: Vector
        double precision                                                :: NormValue
        NormValue = sqrt(sum(Vector**2))
    end function NormDPVector
    function NormRealMatrix(Matrix) result(NormValue)
        implicit none
        integer                                             :: i,j
        real, dimension(:,:), allocatable, intent(in)       :: Matrix
        real                                                :: NormValue
        real                                                :: Result = 0.0
        do i = 1, size(Matrix(:,1)), 1
            do j = 1, size(Matrix(1,:)), 1
                Result = Result + abs(Matrix(i,j))
            end do
        end do
        NormValue = Result
    end function NormRealMatrix
    function NormDPMatrix(Matrix) result(NormValue)
        implicit none
        integer                                                         :: i,j
        double precision, dimension(:,:), allocatable, intent(in)       :: Matrix
        double precision                                                :: NormValue
        double precision                                                :: Result = 0.0d0
        do i = 1, size(Matrix(:,1)), 1
            do j = 1, size(Matrix(1,:)), 1
                Result = Result + abs(Matrix(i,j))
            end do
        end do
        NormValue = Result
    end function NormDPMatrix
    ! 1.7. Unitary Vector
    function UnitVectorReal(Vector) result(UnitVector)
        implicit none
        real, dimension(:), allocatable, intent(in)         :: Vector
        real, dimension(:), allocatable                     :: UnitVector
        UnitVector = Vector/(sqrt(sum(Vector**2)))
    end function UnitVectorReal
    function UnitVectorDP(Vector) result(UnitVector)
        implicit none
        double precision, dimension(:), allocatable, intent(in)         :: Vector
        double precision, dimension(:), allocatable                     :: UnitVector
        UnitVector = Vector/(sqrt(sum(Vector**2)))
    end function UnitVectorDP
    ! 1.8. Add element to vector
    function AddInteger2Vector(Vector,IntegerValue,Location) Result(ExtVector)
        implicit none
        character, intent(in)                                   :: Location
        integer, dimension(:), allocatable, intent(in)          :: Vector
        integer, intent(in)                                     :: IntegerValue
        integer, dimension(:), allocatable                      :: ExtVector
        allocate(ExtVector(Size(Vector)+1))
        if (Location.eq.'B') then    ! al comienzo
            ExtVector(1) = IntegerValue
            ExtVector(2:(size(Vector)+1)) = Vector
        elseif (Location.eq.'E') then ! al final
            ExtVector(1:size(vector)) = Vector
            ExtVector(size(Vector)+1) = IntegerValue
        else
            stop "ERROR Location in Add2Vector"
        end if
    end function AddInteger2Vector
    function AddReal2Vector(Vector,RealValue,Location) Result(ExtVector)
        implicit none
        character, intent(in)                                :: Location
        real, dimension(:), allocatable, intent(in)          :: Vector
        real, intent(in)                                     :: RealValue
        real, dimension(:), allocatable                      :: ExtVector
        allocate(ExtVector(Size(Vector)+1))
        if (Location.eq.'B') then    ! al comienzo
            ExtVector(1) = RealValue
            ExtVector(2:(size(Vector)+1)) = Vector
        elseif (Location.eq.'E') then ! al final
            ExtVector(1:size(vector)) = Vector
            ExtVector(size(Vector)+1) = RealValue
        else
            stop "ERROR Location in Add2Vector"
        end if
    end function AddReal2Vector
    function AddDP2Vector(Vector,DPValue,Location) Result(ExtVector)
        implicit none
        character, intent(in)                                            :: Location
        double precision, dimension(:), allocatable, intent(in)          :: Vector
        double precision, intent(in)                                     :: DPValue
        double precision, dimension(:), allocatable                      :: ExtVector
        allocate(ExtVector(Size(Vector)+1))
        if (Location.eq.'B') then    ! al comienzo
            ExtVector(1) = DPValue
            ExtVector(2:(size(Vector)+1)) = Vector
        elseif(Location.eq.'E') then ! al final
            ExtVector(1:size(vector)) = Vector
            ExtVector(size(Vector)+1) = DPValue
        else
            stop "ERROR Location in Add2Vector"
        end if
    end function AddDP2Vector
    ! 1.9. Add Vector to Matrix
    function AddIntegerVector2Matrix(Matrix,Vector,Location) Result(ExtMatrix)
        implicit none
        integer                                                 :: i
        character, intent(in)                                   :: Location
        integer, dimension(:,:), allocatable, intent(in)        :: Matrix
        integer, dimension(:), allocatable, intent(in)          :: Vector
        integer, dimension(:,:), allocatable                    :: ExtMatrix
        if(size(Vector).eq.size(Matrix,2)) then     ! agregar fila
            allocate(ExtMatrix((size(Matrix,1)+1),size(Matrix,2)))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(1+i,:) = Matrix(i,:)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(size(Matrix,1)+1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(i,:) = Matrix(i,:)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        elseif(size(Vector).eq.size(Matrix,1)) then ! agregar columna
            allocate(ExtMatrix(size(Matrix,1),size(Matrix,2)+1))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(:,1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,1+i) = Matrix(:,i)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(:,size(Matrix,2)+1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,i) = Matrix(:,i)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        else
            stop "ERROR Dimension in AddVector2Matrix"
        end if
    end function AddIntegerVector2Matrix
    function AddRealVector2Matrix(Matrix,Vector,Location) Result(ExtMatrix)
        implicit none
        integer                                              :: i
        character, intent(in)                                :: Location
        real, dimension(:,:), allocatable, intent(in)        :: Matrix
        real, dimension(:), allocatable, intent(in)          :: Vector
        real, dimension(:,:), allocatable                    :: ExtMatrix
        if(size(Vector).eq.size(Matrix,2)) then     ! agregar fila
            allocate(ExtMatrix(size(Matrix,1)+1,size(Matrix,2)))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(1+i,:) = Matrix(i,:)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(size(Matrix,1)+1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(i,:) = Matrix(i,:)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        elseif(size(Vector).eq.size(Matrix,1)) then ! agregar columna
            allocate(ExtMatrix(size(Matrix,1),size(Matrix,2)+1))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(:,1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,1+i) = Matrix(:,i)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(:,size(Matrix,2)+1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,i) = Matrix(:,i)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        else
            stop "ERROR Dimension in AddVector2Matrix"
        end if
    end function AddRealVector2Matrix
    function AddDPVector2Matrix(Matrix,Vector,Location) Result(ExtMatrix)
        implicit none
        integer                                                     :: i
        character, intent(in)                                       :: Location
        double precision, dimension(:,:), allocatable, intent(in)   :: Matrix
        double precision, dimension(:), allocatable, intent(in)     :: Vector
        double precision, dimension(:,:), allocatable               :: ExtMatrix
        if(size(Vector).eq.size(Matrix,2)) then     ! agregar fila
            allocate(ExtMatrix(size(Matrix,1)+1,size(Matrix,2)))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(1+i,:) = Matrix(i,:)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(size(Matrix,1)+1,:) = Vector
                do i = 1, size(Matrix,1), 1
                    ExtMatrix(i,:) = Matrix(i,:)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        elseif(size(Vector).eq.size(Matrix,1)) then ! agregar columna
            allocate(ExtMatrix(size(Matrix,1),size(Matrix,2)+1))
            ExtMatrix = 0.0
            if (Location.eq.'B') then
                ExtMatrix(:,1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,1+i) = Matrix(:,i)
                end do
            elseif(Location.eq.'E') then  
                ExtMatrix(:,size(Matrix,2)+1) = Vector
                do i = 1, size(Matrix,2), 1
                    ExtMatrix(:,i) = Matrix(:,i)
                end do
            else
                stop "ERROR Location in AddVector2Matrix"
            end if
        else
            stop "ERROR Dimension in AddVector2Matrix"
        end if
    end function AddDPVector2Matrix
    ! 1.10 Add Matrix to MAtrix
    function AddIntegerMatrix2Matrix(BaseMatrix,AddMatrix,Location) Result(ExtMatrix)
        implicit none
        integer                                                 :: i,ColBaseMatrix,RowBaseMatrix,ColAddMatrix,RowAddMatrix
        character, intent(in)                                   :: Location
        integer, dimension(:,:), allocatable, intent(in)        :: BaseMatrix
        integer, dimension(:,:), allocatable, intent(in)        :: AddMatrix
        integer, dimension(:,:), allocatable                    :: ExtMatrix
        RowBaseMatrix = size(BaseMatrix,1)
        ColBaseMatrix = size(BaseMatrix,2) 
        RowAddMatrix = size(AddMatrix,1)
        ColAddMatrix = size(AddMatrix,2)
        if (Location.eq.'U') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowAddMatrix,:) = AddMatrix
            ExtMatrix((RowAddMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = BaseMatrix
        elseif (Location.eq.'D') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowBaseMatrix,:) = BaseMatrix
            ExtMatrix((RowBaseMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = AddMatrix
        elseif (Location.eq.'L') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension" 
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColAddMatrix) = AddMatrix
            ExtMatrix(:,(ColAddMatrix+1):(ColBaseMatrix+ColAddMatrix)) = BaseMatrix
        elseif (Location.eq.'R') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColBaseMatrix) = BaseMatrix
            ExtMatrix(:,(ColBaseMatrix+1):(ColBaseMatrix+ColAddMatrix)) = AddMatrix
        end if
    end function AddIntegerMatrix2Matrix
    function AddRealMatrix2Matrix(BaseMatrix,AddMatrix,Location) Result(ExtMatrix)
        implicit none
        integer                                                 :: i,ColBaseMatrix,RowBaseMatrix,ColAddMatrix,RowAddMatrix
        character, intent(in)                                   :: Location
        real, dimension(:,:), allocatable, intent(in)           :: BaseMatrix
        real, dimension(:,:), allocatable, intent(in)           :: AddMatrix
        real, dimension(:,:), allocatable                       :: ExtMatrix
        RowBaseMatrix = size(BaseMatrix,1)
        ColBaseMatrix = size(BaseMatrix,2) 
        RowAddMatrix = size(AddMatrix,1)
        ColAddMatrix = size(AddMatrix,2)
        if (Location.eq.'U') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowAddMatrix,:) = AddMatrix
            ExtMatrix((RowAddMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = BaseMatrix
        elseif (Location.eq.'D') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowBaseMatrix,:) = BaseMatrix
            ExtMatrix((RowBaseMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = AddMatrix
        elseif (Location.eq.'L') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension" 
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColAddMatrix) = AddMatrix
            ExtMatrix(:,(ColAddMatrix+1):(ColBaseMatrix+ColAddMatrix)) = BaseMatrix
        elseif (Location.eq.'R') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColBaseMatrix) = BaseMatrix
            ExtMatrix(:,(ColBaseMatrix+1):(ColBaseMatrix+ColAddMatrix)) = AddMatrix
        end if
    end function AddRealMatrix2Matrix
    function AddDPMatrix2Matrix(BaseMatrix,AddMatrix,Location) Result(ExtMatrix)
        implicit none
        integer                                                     :: i,ColBaseMatrix,RowBaseMatrix,ColAddMatrix,RowAddMatrix
        character, intent(in)                                       :: Location
        double precision, dimension(:,:), allocatable, intent(in)   :: BaseMatrix
        double precision, dimension(:,:), allocatable, intent(in)   :: AddMatrix
        double precision, dimension(:,:), allocatable               :: ExtMatrix
        RowBaseMatrix = size(BaseMatrix,1)
        ColBaseMatrix = size(BaseMatrix,2) 
        RowAddMatrix = size(AddMatrix,1)
        ColAddMatrix = size(AddMatrix,2)
        if (Location.eq.'U') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowAddMatrix,:) = AddMatrix
            ExtMatrix((RowAddMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = BaseMatrix
        elseif (Location.eq.'D') then
            if (ColBaseMatrix.ne.ColAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix((RowAddMatrix+RowBaseMatrix),ColBaseMatrix))
            ExtMatrix(1:RowBaseMatrix,:) = BaseMatrix
            ExtMatrix((RowBaseMatrix+1):(RowAddMatrix+RowBaseMatrix),:) = AddMatrix
        elseif (Location.eq.'L') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension" 
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColAddMatrix) = AddMatrix
            ExtMatrix(:,(ColAddMatrix+1):(ColBaseMatrix+ColAddMatrix)) = BaseMatrix
        elseif (Location.eq.'R') then
            if (RowBaseMatrix.ne.RowAddMatrix) stop "ERROR Matrix dimension"
            allocate(ExtMatrix(RowBaseMatrix,(ColBaseMatrix+ColAddMatrix)))
            ExtMatrix(:,1:ColBaseMatrix) = BaseMatrix
            ExtMatrix(:,(ColBaseMatrix+1):(ColBaseMatrix+ColAddMatrix)) = AddMatrix
        end if
    end function AddDPMatrix2Matrix
    ! Some Additional base functions
    ! 1. printing
    ! 1.1. ...in console
    subroutine PrintIntegerVectorInConsole(IntegerVector,Direction)
        implicit none
        integer                                             :: i
        character, intent(in)                               :: Direction
        integer, dimension(:), allocatable, intent(in)      :: IntegerVector
        if (Direction.eq.'H') then
            write(unit=*, fmt=*) IntegerVector
        elseif (Direction.eq.'V') then
            do i = 1, size(IntegerVector), 1
                write(unit=*, fmt=*) IntegerVector(i)
            end do
        else
            stop "ERROR Direction Console printing"
        end if
    end subroutine PrintIntegerVectorInConsole
    subroutine PrintRealVectorInConsole(RealVector,Direction)
        implicit none
        integer                                             :: i
        character, intent(in)                               :: Direction
        real, dimension(:), allocatable, intent(in)         :: RealVector
        if (Direction.eq.'H') then
            write(unit=*, fmt=*) RealVector
        elseif (Direction.eq.'V') then
            do i = 1, size(RealVector), 1
                write(unit=*, fmt=*) RealVector(i)
            end do
        else
            stop "ERROR Direction Console printing"
        end if
    end subroutine PrintRealVectorInConsole
    subroutine PrintDPVectorInConsole(DPVector,Direction)
        implicit none
        integer                                                   :: i
        character, intent(in)                                     :: Direction
        double precision, dimension(:), allocatable, intent(in)   :: DPVector
        if (Direction.eq.'H') then
            write(unit=*, fmt=*) DPVector
        elseif (Direction.eq.'V') then
            do i = 1, size(DPVector), 1
                write(unit=*, fmt=*) DPVector(i)
            end do
        else
            stop "ERROR Direction Console printing"
        end if
    end subroutine PrintDPVectorInConsole
    subroutine PrintIntegerMatrixInConsole(IntegerMatrix)
        implicit none
        integer                                             :: i
        integer, dimension(:,:), allocatable, intent(in)    :: IntegerMatrix
        do i = 1, size(IntegerMatrix,1), 1
            write(unit=*, fmt=*) IntegerMatrix(i,:)
        end do
    end subroutine PrintIntegerMatrixInConsole
    subroutine PrintRealMatrixInConsole(RealMatrix)
        implicit none
        integer                                          :: i
        real, dimension(:,:), allocatable, intent(in)    :: RealMatrix
        do i = 1, size(RealMatrix,1), 1
            write(unit=*, fmt=*) RealMatrix(i,:)
        end do
    end subroutine PrintRealMatrixInConsole
    subroutine PrintDPMatrixInConsole(DPMatrix)
        implicit none
        integer                                                      :: i
        double precision, dimension(:,:), allocatable, intent(in)    :: DPMatrix
        do i = 1, size(DPMatrix,1), 1
            write(unit=*, fmt=*) DPMatrix(i,:)
        end do
    end subroutine PrintDPMatrixInConsole
    ! 1.2. ...in file
    subroutine PrintIntegerVectorInFile(IntegerVector,Direction,Path)
        implicit none
        integer                                             :: i,ios,iounit
        character, intent(in)                               :: Direction
        character(len=*), intent(in)                        :: Path
        integer, dimension(:), allocatable, intent(in)      :: IntegerVector
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Direction.eq.'H') then
                write(unit=iounit, fmt=*) IntegerVector
            elseif (Direction.eq.'V') then
                do i = 1, size(IntegerVector), 1
                    write(unit=iounit, fmt=*) IntegerVector(i)
                end do
            else
                stop "ERROR Direction Console printing"
            end if
        close(iounit)
    end subroutine PrintIntegerVectorInFile
    subroutine PrintRealVectorInFile(RealVector,Direction,Path)
        implicit none
        integer                                             :: i,ios,iounit
        character, intent(in)                               :: Direction
        character(len=*), intent(in)                        :: Path
        real, dimension(:), allocatable, intent(in)         :: RealVector
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Direction.eq.'H') then
                write(unit=iounit, fmt=*) RealVector
            elseif (Direction.eq.'V') then
                do i = 1, size(RealVector), 1
                    write(unit=iounit, fmt=*) RealVector(i)
                end do
            else
                stop "ERROR Direction Console printing"
            end if
        close(iounit)
    end subroutine PrintRealVectorInFile
    subroutine PrintDPVectorInFile(DPVector,Direction,Path)
        implicit none
        integer                                                   :: i,ios,iounit
        character, intent(in)                                     :: Direction
        character(len=*), intent(in)                              :: Path
        double precision, dimension(:), allocatable, intent(in)   :: DPVector
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Direction.eq.'H') then
                write(unit=iounit, fmt=*) DPVector
            elseif (Direction.eq.'V') then
                do i = 1, size(DPVector), 1
                    write(unit=iounit, fmt=*) DPVector(i)
                end do
            else
                stop "ERROR Direction Console printing"
            end if
        close(iounit)
    end subroutine PrintDPVectorInFile
    subroutine PrintIntegerMatrixInFile(IntegerMatrix,Path)
        implicit none
        integer                                                   :: i,ios,iounit
        character(len=*), intent(in)                              :: Path
        integer, dimension(:,:), allocatable, intent(in)          :: IntegerMatrix
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(IntegerMatrix,1), 1
                write(unit=iounit, fmt=*) IntegerMatrix(i,:)
            end do
        close(iounit)
    end subroutine PrintIntegerMatrixInFile
    subroutine PrintRealMatrixInFile(RealMatrix,Path)
        implicit none
        integer                                          :: i,ios,iounit
        character(len=*), intent(in)                     :: Path
        real, dimension(:,:), allocatable, intent(in)    :: RealMatrix
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(RealMatrix,1), 1
                write(unit=iounit, fmt=*) RealMatrix(i,:)
            end do
        close(iounit)
    end subroutine PrintRealMatrixInFile
    subroutine PrintDPMatrixInFile(DPMatrix,Path)
        implicit none
        integer                                                     :: i,ios,iounit
        character(len=*), intent(in)                                :: Path
        double precision, dimension(:,:), allocatable, intent(in)   :: DPMatrix
        open(unit=iounit, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            do i = 1, size(DPMatrix,1), 1
                write(unit=iounit, fmt=*) DPMatrix(i,:)
            end do
        close(iounit)
    end subroutine PrintDPMatrixInFile
    ! 2. Plotting with gnuplot (functions still under construction :D)
    ! 2.1 Scatter plot
    subroutine ScatterPlotReal(Xvalues,Yvalues)
        implicit none
        real, dimension(:), allocatable, intent(in)         :: Xvalues
        real, dimension(:), allocatable, intent(in)         :: Yvalues
        CHARACTER(200)                                      :: CommandLine
        integer                                             :: i,iounit,ios
        character(len=*), parameter                         :: name = 'DumbFile.txt'
        open(unit=iounit, file=name, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Size(Xvalues).ne.size(Yvalues)) stop "ScatterPlot ERROR Dimensions" 
            do i = 1, size(Xvalues), 1
                write(unit=iounit, fmt=*) Xvalues(i),Yvalues(i)
            end do
            CommandLine = 'gnuplot -persist << EOF' // NEW_LINE('A') // &
            'set terminal x11' // NEW_LINE('A') // &
            'plot "DumbFile.txt" using 1:2 with points' // NEW_LINE('A') // &
            'EOF'
            call SYSTEM(CommandLine)
        close(iounit,status="delete")
    end subroutine ScatterPlotReal
    subroutine ScatterPlotDP(Xvalues,Yvalues)
        implicit none
        double precision, dimension(:), allocatable, intent(in)     :: Xvalues
        double precision, dimension(:), allocatable, intent(in)     :: Yvalues
        CHARACTER(200)                                              :: CommandLine
        integer                                                     :: i,iounit,ios
        character(len=*), parameter                                 :: name = 'DumbFile.txt'
        open(unit=iounit, file=name, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Size(Xvalues).ne.size(Yvalues)) stop "ScatterPlot ERROR Dimensions" 
            do i = 1, size(Xvalues), 1
                write(unit=iounit, fmt=*) Xvalues(i),Yvalues(i)
            end do
            CommandLine = 'gnuplot -persist << EOF' // NEW_LINE('A') // &
            'set terminal x11' // NEW_LINE('A') // &
            'plot "DumbFile.txt" using 1:2 with points' // NEW_LINE('A') // &
            'EOF'
            call SYSTEM(CommandLine)
        close(iounit,status="delete")
    end subroutine ScatterPlotDP
    ! 2.2 Line plot
    subroutine LinePlotReal(Xvalues,Yvalues)
        implicit none
        real, dimension(:), allocatable, intent(in)         :: Xvalues
        real, dimension(:), allocatable, intent(in)         :: Yvalues
        CHARACTER(200)                                      :: CommandLine
        integer                                             :: i,iounit,ios
        character(len=*), parameter                         :: name = 'DumbFile.txt'
        open(unit=iounit, file=name, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Size(Xvalues).ne.size(Yvalues)) stop "ScatterPlot ERROR Dimensions" 
            do i = 1, size(Xvalues), 1
                write(unit=iounit, fmt=*) Xvalues(i),Yvalues(i)
            end do
            CommandLine = 'gnuplot -persist << EOF' // NEW_LINE('A') // &
            'set terminal x11' // NEW_LINE('A') // &
            'plot "DumbFile.txt" using 1:2 with lines' // NEW_LINE('A') // &
            'EOF'
            call SYSTEM(CommandLine)
        close(iounit,status="delete")
    end subroutine LinePlotReal
    subroutine LinePlotDP(Xvalues,Yvalues)
        implicit none
        double precision, dimension(:), allocatable, intent(in)     :: Xvalues
        double precision, dimension(:), allocatable, intent(in)     :: Yvalues
        CHARACTER(200)                                              :: CommandLine
        integer                                                     :: i,iounit,ios
        character(len=*), parameter                                 :: name = 'DumbFile.txt'
        open(unit=iounit, file=name, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name"
            if (Size(Xvalues).ne.size(Yvalues)) stop "ScatterPlot ERROR Dimensions" 
            do i = 1, size(Xvalues), 1
                write(unit=iounit, fmt=*) Xvalues(i),Yvalues(i)
            end do
            CommandLine = 'gnuplot -persist << EOF' // NEW_LINE('A') // &
            'set terminal x11' // NEW_LINE('A') // &
            'plot "DumbFile.txt" using 1:2 with lines' // NEW_LINE('A') // &
            'EOF'
            call SYSTEM(CommandLine)
        close(iounit,status="delete")
    end subroutine LinePlotDP
    ! ...
end module Base_Module
