module Paraview_Module
    implicit none
    type                                            :: PostprocessingInfo
        character(len=20)                           :: PostprocessingPath
        character(len=20)                           :: Selection
        double precision                            :: FilterValue
    contains
        procedure                                   :: SetPostprocessingPath
        procedure                                   :: SetFilterValue
        procedure                                   :: PlottingSelection
        procedure                                   :: ParaviewPostprocessing
    end type
    
contains
    ! ----------------------------------------------------------------- !
    !                           base subroutines                        !
    ! ----------------------------------------------------------------- !
    ! Generation of Geometry file 
    subroutine ParaviewGeometryFile(Path,Coordinates,Connectivity,Element)
        implicit none
        ! input
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        integer, dimension(:,:), allocatable, intent(in)           :: Connectivity
        double precision, dimension(:,:), allocatable, intent(in)  :: Coordinates
        ! internal 
        integer                                                    :: i, ios
        ! Generation process
        ! Geometry File
        open(unit=1, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_geometryfile"
            write(unit=1,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'This is the 1st description line of the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'node id given'
            write(unit=1,fmt='(a)') 'element id given'
            write(unit=1,fmt='(a)') 'extents'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') '-100.00000+00 100.00000+00'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') '2D uns-elements (description line for part 1)'
            write(unit=1,fmt='(a)') 'coordinates'
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt=*) size(Coordinates,1)
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,1)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                write(unit=1,fmt=*) coordinates(i,2)
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Coordinates,1), 1
                if (size(Coordinates,2).eq.2) then
                    write(unit=1,fmt=*) 0
                elseif (Size(Coordinates,2).eq.3) then
                    write(unit=1,fmt=*) coordinates(i,3)
                else
                    stop "ERROR Dimension-Paraview"
                end if
            end do
            write(unit=1,fmt='(a)') ''
            write(unit=1,fmt='(a)') Element
            write(unit=1,fmt=*) size(Connectivity,1)
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) i
            end do
            write(unit=1,fmt='(a)') ''
            do i = 1, size(Connectivity,1), 1
                write(unit=1,fmt=*) Connectivity(i,:)
            end do
        close(1)
    end subroutine ParaviewGeometryFile

    Subroutine ParaviewScalarFilePerElement(Path,Result,Element)
        implicit none
        ! input
        character(len=*), intent(in)                               :: Path
        character(len=*), intent(in)                               :: Element
        double precision, dimension(:), allocatable, intent(in)    :: Result
        ! Internal
        integer                                                    :: i, ios
        ! Scalar File
        open(unit=1, file=Path, iostat=ios, status="replace", action="write")
            if ( ios /= 0 ) stop "Error opening file name_coordinatesfile"
            write(unit=1,fmt='(a)') 'Per_element scalar values for the EnSight Gold geometry example'
            write(unit=1,fmt='(a)') 'part'
            write(unit=1,fmt='(a)') '1'
            write(unit=1,fmt='(a)') Element
            do i = 1, size(Result), 1
                write(unit=1,fmt=*) Result(i)
            end do
        close(1)
    end subroutine ParaviewScalarFilePerElement
    ! ----------------------------------------------------------------- !
    !          subroutines to define the information required           !
    ! ----------------------------------------------------------------- !
    ! listo
    subroutine SetPostprocessingPath(self,path)
        implicit none
        character(len=*), intent(in)                :: path
        class(PostprocessingInfo), intent(inout)    :: self
        self%PostprocessingPath = path
    end subroutine
    ! listo
    subroutine PlottingSelection(self,info)
        implicit none
        character(len=*), intent(in)                :: info
        class(PostprocessingInfo), intent(inout)    :: self 
        self%Selection = info
    end subroutine PlottingSelection
    ! listo
    subroutine SetFilterValue(self,FilterValue)
        implicit none
        double precision, intent(in)                            :: FilterValue
        class(PostprocessingInfo), intent(inout)    :: self 
        Self%FilterValue = FilterValue
    end subroutine SetFilterValue
    ! ----------------------------------------------------------------- !
    !             subroutines to define the Postprocessing              !
    ! ----------------------------------------------------------------- !
    ! Rutina para generar archivos de paraview
    subroutine ParaviewPostprocessing(self1,self2)
        use Optimization_Module
        implicit none
        class(PostprocessingInfo), intent(inout)        :: self1
        class(Optimization), intent(inout)              :: self2
        ! internal variables
        integer                                         :: i,j,ios
        integer, dimension(:), allocatable              :: Index
        integer, dimension(:,:), allocatable            :: Connectivity
        double precision, dimension(:), allocatable     :: ResultVector
        double precision, dimension(:,:), allocatable   :: Coordinates
        logical, dimension(:), allocatable              :: InLogical
        character(len=14), dimension(:), allocatable    :: EstadosCarga
        character(len=80)                               :: path
        character(len=1)                                :: ind
        ! --------------------------------------------------------------!
        ! Applying filter (deleting elements with low density)
        Coordinates = Self2%Coordinates
        j = size(Self2%DensityVector)
        InLogical = Self2%DensityVector.gt.self1%FilterValue
        Index = Pack([(i,i=1,j)],InLogical)
        Connectivity = Self2%ConnectivityN(Index,:)
        ! ---------------------------------------------------------------!
        ! star postprocessing
        if (self2%DimAnalysis.eq.2) then
            allocate(EstadosCarga(3))
            EstadosCarga(1) = '/LoadCase1-Exx'
            EstadosCarga(2) = '/LoadCase2-Eyy'
            EstadosCarga(3) = '/LoadCase3-Exy'
        elseif (self2%DimAnalysis.eq.3) then
            allocate(EstadosCarga(6))
            EstadosCarga(1) = '/LoadCase1-Exx'
            EstadosCarga(2) = '/LoadCase2-Eyy'
            EstadosCarga(3) = '/LoadCase3-Ezz'
            EstadosCarga(4) = '/LoadCase4-Exy'
            EstadosCarga(5) = '/LoadCase5-Eyz'
            EstadosCarga(6) = '/LoadCase6-Exz'
        end if
        write(unit=*, fmt=*) 'PostProcessing FEM results - Paraview'
        ! Generation of .geom file
        path = trim(self1%PostprocessingPath) // '/Geometry.geom'
        write(unit=*, fmt=*) '** Generating Geometry file'
        call ParaviewGeometryFile(path,Coordinates,Connectivity,self2%ElementType)
        ! note: 
        ! Deformation - DeformationEqv 
        !    Stress   - StressVonMises 
        !    Strain   -   StrainEqv
        !  DensEnergy
        if (self1%Selection.eq.'Stress') then !listo
            select case (self2%DimAnalysis)
                case (2)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    do j = 1, 3, 1
                        ! Generation of .case file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//EstadosCarga(j)//'-Stress.case'
                        open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                            if ( ios /= 0 ) stop "Error opening file name_casefile"
                            write(unit=1,fmt='(a)') 'FORMAT'
                            write(unit=1,fmt='(a)') 'type: ensight gold'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'GEOMETRY'
                            write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'VARIABLE'
                            path = 'Elem-StressXX'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressXX' // ' ' // path
                            path = 'Elem-StressYY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressYY' // ' ' // path
                            path = 'Elem-StressXY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressXY' // ' ' // path
                        close(1)
                        ! Generation of .esca file
                        ! ** per element
                        path = trim(self1%PostprocessingPath) //'/Elem-StressXX'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,1)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressYY'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressXY'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,3)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
                case (3)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    do j = 1, 6, 1
                        ! Generation of .case file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//EstadosCarga(j)//'-Stress.case'
                        open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                            if ( ios /= 0 ) stop "Error opening file name_casefile"
                            write(unit=1,fmt='(a)') 'FORMAT'
                            write(unit=1,fmt='(a)') 'type: ensight gold'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'GEOMETRY'
                            write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'VARIABLE'
                            path = 'Elem-StressXX'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressXX' // ' ' // path
                            path = 'Elem-StressYY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressYY' // ' ' // path
                            path = 'Elem-StressZZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressZZ' // ' ' // path
                            path = 'Elem-StressXY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressXY' // ' ' // path
                            path = 'Elem-StressYZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressYZ' // ' ' // path
                            path = 'Elem-StressXZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StressXZ' // ' ' // path
                        close(1)
                        ! Generation of .esca file
                        ! ** per element
                        path = trim(self1%PostprocessingPath) //'/Elem-StressXX'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,1)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressYY'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressZZ'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,3)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressXY'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,4)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressYZ'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,5)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath) //'/Elem-StressXZ'// ind //'.esca'
                        ResultVector = self2%StressE(j,Index,6)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
            end select
        elseif (self1%Selection.eq.'StressVonMises') then !listo
            select case (self2%DimAnalysis)
                case (2)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    ! Generation of .case file
                    path = trim(self1%PostprocessingPath)//'/StressVonMises.case'
                    open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                        if ( ios /= 0 ) stop "Error opening file name_casefile"
                        write(unit=1,fmt='(a)') 'FORMAT'
                        write(unit=1,fmt='(a)') 'type: ensight gold'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'GEOMETRY'
                        write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'VARIABLE'
                        path = 'Elem-StressVonMises1.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXX' // ' ' // path
                        path = 'Elem-StressVonMises2.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesYY' // ' ' // path
                        path = 'Elem-StressVonMises3.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXY' // ' ' // path
                    close(1)
                    do j = 1, 3, 1
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//'/Elem-StressVonMises'// ind //'.esca'
                        ResultVector = sqrt(self2%StressE(j,Index,1)**2-self2%StressE(j,Index,1)*self2%StressE(j,Index,2) &
                                            +self2%StressE(j,Index,2)**2+3*self2%StressE(j,Index,3)**2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
                case (3)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    do j = 1, 6, 1
                        ! Generation of .case file
                        path = trim(self1%PostprocessingPath)//'/StressVonMises.case'
                        open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                            if ( ios /= 0 ) stop "Error opening file name_casefile"
                            write(unit=1,fmt='(a)') 'FORMAT'
                            write(unit=1,fmt='(a)') 'type: ensight gold'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'GEOMETRY'
                            write(unit=1,fmt='(a)') 'model:     ' //'Geometry.geom'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'VARIABLE'
                            path = 'Elem-StressVonMises1.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXX' // ' ' // path
                            path = 'Elem-StressVonMises2.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesYY' // ' ' // path
                            path = 'Elem-StressVonMises3.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesZZ' // ' ' // path
                            path = 'Elem-StressVonMises4.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXY' // ' ' // path
                            path = 'Elem-StressVonMises5.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesYZ' // ' ' // path
                            path = 'Elem-StressVonMises6.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXZ' // ' ' // path
                        close(1)
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//'/Elem-StressVonMises'// ind //'.esca'
                        ResultVector = sqrt(self2%StressE(j,Index,1)**2 + self2%StressE(j,Index,2)**2 &
                                        + self2%StressE(j,Index,3)**2 - (self2%StressE(j,Index,1)*self2%StressE(j,Index,2) &
                                        + self2%StressE(j,Index,2)*self2%StressE(j,Index,3) &
                                        + self2%StressE(j,Index,1)*self2%StressE(j,Index,3)) + 3*(self2%StressE(j,Index,4)**2 &
                                        + self2%StressE(j,Index,5)**2 + self2%StressE(j,Index,6)**2))
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
            end select
        elseif (self1%Selection.eq.'Strain') then !listo
            select case (self2%DimAnalysis)
                case (2)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    do j = 1, 3, 1
                        ! Generation of .case file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//EstadosCarga(j)//'-Strain.case'
                        open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                            if ( ios /= 0 ) stop "Error opening file name_casefile"
                            write(unit=1,fmt='(a)') 'FORMAT'
                            write(unit=1,fmt='(a)') 'type: ensight gold'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'GEOMETRY'
                            write(unit=1,fmt='(a)') 'model:     ' //'Geometry.geom'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'VARIABLE'
                            path = 'Elem-StrainXX'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainXX' // ' ' // path
                            path = 'Elem-StrainYY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainYY' // ' ' // path
                            path = 'Elem-StrainXY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainXY' // ' ' // path
                        close(1)
                        ! Generation of .esca file
                        ! ** per element
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainXX'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,1)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainYY'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainXY'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,3)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
                case (3)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    do j = 1, 6, 1
                        ! Generation of .case file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//EstadosCarga(j)//'-Strain.case'
                        open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                            if ( ios /= 0 ) stop "Error opening file name_casefile"
                            write(unit=1,fmt='(a)') 'FORMAT'
                            write(unit=1,fmt='(a)') 'type: ensight gold'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'GEOMETRY'
                            write(unit=1,fmt='(a)') 'model:     ' //'Geometry.geom'
                            write(unit=1,fmt='(a)') ''
                            write(unit=1,fmt='(a)') 'VARIABLE'
                            path = 'Elem-StrainXX'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainXX' // ' ' // path
                            path = 'Elem-StrainYY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainYY' // ' ' // path
                            path = 'Elem-StrainZZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainZZ' // ' ' // path
                            path = 'Elem-StrainXY'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainXY' // ' ' // path
                            path = 'Elem-StrainYZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainYZ' // ' ' // path
                            path = 'Elem-StrainXZ'// ind //'.esca'
                            write(unit=1,fmt='(a)') 'scalar per element:     ' // 'Elem-StrainXZ' // ' ' // path
                        close(1)
                        ! ** per element
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainXX'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,1)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainYY'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainZZ'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,3)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainXY'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,4)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainYZ'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,5)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                        path = trim(self1%PostprocessingPath)//'/Elem-StrainXZ'// ind //'.esca'
                        ResultVector = self2%StrainE(j,Index,6)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
            end select
        elseif (self1%selection.eq.'StrainEqv') then !list
            select case (self2%DimAnalysis)
                case (2)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    ! Generation of .case file
                    path = trim(self1%PostprocessingPath)//'/StrainEqv.case'
                    open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                        if ( ios /= 0 ) stop "Error opening file name_casefile"
                        write(unit=1,fmt='(a)') 'FORMAT'
                        write(unit=1,fmt='(a)') 'type: ensight gold'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'GEOMETRY'
                        write(unit=1,fmt='(a)') 'model:     ' //'Geometry.geom'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'VARIABLE'
                        path = 'StrainEqv1.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StrainEqv-XX' // ' ' // path
                        path = 'StrainEqv2.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StrainEqv-YY' // ' ' // path
                        path = 'StrainEqv3.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StrainEqv-XY' // ' ' // path
                    close(1)
                    do j = 1, 3, 1
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//'/StrainEqv'// ind //'.esca'
                        ResultVector = sqrt(self2%StrainE(j,Index,1)**2 + self2%StrainE(j,Index,2)**2 &
                                            + 2*self2%StrainE(j,Index,3)**2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
                case (3)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    ! Generation of .case file
                    path = trim(self1%PostprocessingPath)//'/StrainEqv.case'
                    open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                        if ( ios /= 0 ) stop "Error opening file name_casefile"
                        write(unit=1,fmt='(a)') 'FORMAT'
                        write(unit=1,fmt='(a)') 'type: ensight gold'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'GEOMETRY'
                        write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'VARIABLE'
                        path = 'StrainEqv1.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXX' // ' ' // path
                        path = 'StrainEqv2.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesYY' // ' ' // path
                        path = 'StrainEqv3.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesZZ' // ' ' // path
                        path = 'StrainEqv4.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXY' // ' ' // path
                        path = 'StrainEqv5.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesYZ' // ' ' // path
                        path = 'StrainEqv6.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'StressVonMisesXZ' // ' ' // path
                    close(1)
                    do j = 1, 6, 1
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath)//'/StrainEqv'// ind //'.esca'
                        ResultVector = sqrt(self2%StrainE(j,Index,1)**2 + self2%StrainE(j,Index,2)**2 &
                                            + self2%StrainE(j,Index,3)**2 + 2*self2%StrainE(j,:,4)**2 &
                                            + 2*self2%StrainE(j,Index,5)**2 + 2*self2%StrainE(j,Index,6)**2)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
            end select
        elseif (self1%Selection.eq.'DensEnergy') then !listo
            select case (Self2%DimAnalysis)
                case (2)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    ! Generation of .case file
                    path = trim(self1%PostprocessingPath)//'/DensEnergy.case'
                    open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                        if ( ios /= 0 ) stop "Error opening file name_casefile"
                        write(unit=1,fmt='(a)') 'FORMAT'
                        write(unit=1,fmt='(a)') 'type: ensight gold'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'GEOMETRY'
                        write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'VARIABLE'
                        path = 'DensEnergy1.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-XX' // ' ' // path
                        path = 'DensEnergy2.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-YY' // ' ' // path
                        path = 'DensEnergy3.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-XY' // ' ' // path
                    close(1)
                    do j = 1, 3, 1
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath) //'/DensEnergy'// ind //'.esca'
                        ResultVector = self2%StrainEnergyE(j,Index)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
                case (3)
                    write(unit=*, fmt=*) '** Generating .case and .esca files'
                    ! Generation of .case file
                    path = trim(self1%PostprocessingPath)//'/DensEnergy.case'
                    open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                        if ( ios /= 0 ) stop "Error opening file name_casefile"
                        write(unit=1,fmt='(a)') 'FORMAT'
                        write(unit=1,fmt='(a)') 'type: ensight gold'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'GEOMETRY'
                        write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                        write(unit=1,fmt='(a)') ''
                        write(unit=1,fmt='(a)') 'VARIABLE'
                        path = 'DensEnergy1.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-XX' // ' ' // path
                        path = 'DensEnergy2.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-YY' // ' ' // path
                        path = 'DensEnergy3.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-ZZ' // ' ' // path
                        path = 'DensEnergy4.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-XY' // ' ' // path
                        path = 'DensEnergy5.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-YZ' // ' ' // path
                        path = 'DensEnergy6.esca'
                        write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensEnergy-XZ' // ' ' // path
                    close(1)
                    do j = 1, 6, 1
                        ! Generation of .esca file
                        write(unit=ind, fmt='(I1)') j
                        path = trim(self1%PostprocessingPath) //'/DensEnergy'// ind //'.esca'
                        ResultVector = self2%StrainEnergyE(j,Index)
                        call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
                    end do
            end select
        elseif (Self1%selection.eq.'TopOpt') then ! Listo
            ! Generation of .case file
            write(unit=*, fmt=*) '** Generating .case and .esca files'
            path = trim(self1%PostprocessingPath)//'/TopOpt.case'
            open(unit=1, file = Path, iostat=ios, status="replace", action="write")
                if ( ios /= 0 ) stop "Error opening file name_casefile"
                write(unit=1,fmt='(a)') 'FORMAT'
                write(unit=1,fmt='(a)') 'type: ensight gold'
                write(unit=1,fmt='(a)') ''
                write(unit=1,fmt='(a)') 'GEOMETRY'
                write(unit=1,fmt='(a)') 'model:     ' // 'Geometry.geom'
                write(unit=1,fmt='(a)') ''
                write(unit=1,fmt='(a)') 'VARIABLE'
                path = '/TopOpt.esca'
                write(unit=1,fmt='(a)') 'scalar per element:     ' // 'DensityDistribution' // ' ' // path
            close(1)
            ! Generation of .esca file
            path = trim(self1%PostprocessingPath) //'/TopOpt.esca'
            ResultVector = self2%DensityVector(Index)
            call ParaviewScalarFilePerElement(Path,ResultVector,self2%ElementType)
        end if
    end subroutine ParaviewPostprocessing
end module Paraview_Module