! COPYRIGHT (c) 2009 Council for the Central Laboratory
!               of the Research Councils
! Original date 20 October 2009. Version 1.0.0.

! Fortran 95 version of the mc34 package.

! 18 May 2010 Version 1.1.0 -jhogg
!             Create hsl_mc34_integer
!             Change from logical Hermitian to integer sym_type to cope with
!             skew symmetric matrices as well as Hermitian and symmetric ones.

! to change precision:
!    change _double, kind(0.0d0)
! For complex version:
!    change real to complex
! For Hermitian case, take conjugate for upper triangle entries
! For integer version:
!    change real to integer, kind(0)

   module hsl_mc34_single
   implicit none
   private
   public mc34_expand

   integer, parameter :: wp = kind(0.0)

   interface mc34_expand
      module procedure mc34_expand_single
   end interface

   contains

      subroutine mc34_expand_single(n,row,ptr,iw,a,sym_type)

!  this subroutine generates the expanded structure for a
!   matrix a with a symmetric sparsity pattern given the structure 
!   for the lower triangular part.  diagonal entries need not be present.

      integer, intent(in) :: n  ! holds the order of a.

      integer, intent(inout) :: row(*) ! must be set by the user to
!       hold the row indices of the lower triangular part of a.
!       the entries of a single column must be
!       contiguous. the entries of column j must precede those of column
!       j+1, and there must be no wasted space between
!       columns. row indices within a column may be in any order.  on
!       exit, it will have the same meaning but will be changed to hold
!       the row indices of the entries in the expanded structure.  diagonal
!       entries need not be present. the new row indices added in the
!       upper triangular part will be in order for each column and will
!       precede the row indices for the lower triangular part which will
!       remain in the input order.

     integer, intent(inout) ::ptr(n+1)  !  must be set
!       by the user so that ptr(j) is the position in row
!       of the first entry in column j and
!       ptr(n+1) must be set to one more than the total number of
!       entries.  on exit, ptr(j) will have the same meaning but
!       will be changed to point to the position of the first entry of
!       column j in the expanded structure. the new value of
!       ptr(n+1) will be one greater than the number of entries in
!       the expanded structure.

     integer :: iw(n) ! workspace

     real(wp), optional, intent(inout) :: a(*) 
!       if present, a(1:ptr(n+1)-1) must be set by the user so that
!       a(k) holds the value of the entry in row(k). 
!       on exit, a will hold the values of the entries in the expanded 
!       structure corresponding to the output values of row.

     integer, optional, intent(in) :: sym_type 
!      if present with value 1, matrix is skew symmetric.
!      if present with value 2, matrix is hermitian.
!      otherwise matrix is symmetric.

      integer :: ckp1 ! used as running pointer
      integer :: i,i1,i2,ii,ipkp1,ipos
      integer :: j,jstart 
      integer :: lenk ! number of entries in col. j of original structure
      integer :: ndiag ! number diagonal entries present
      integer :: newtau ! number of entries in expanded storage
      integer :: oldtau ! number of entries in symmetric storage
      integer :: r_sym_type ! real sym_type value (used as argument is optional)

      oldtau = ptr(n+1) - 1
      iw(1:n) = 0

! iw(j) set to total number entries in col. j of expanded mx.
      ndiag = 0
      do j = 1,n
        i1 = ptr(j)
        i2 = ptr(j+1) - 1
        iw(j) = iw(j) + i2 - i1 + 1
        do ii = i1,i2
          i = row(ii)
          if (i /= j) then
            iw(i) = iw(i) + 1
          else
            ndiag = ndiag + 1
          end if
        end do
      end do

      newtau = 2*oldtau - ndiag
! ipkp1 points to position  after end of column being currently processed
      ipkp1 = oldtau + 1
! ckp1 points to position  after end of same column in expanded structure
      ckp1 = newtau + 1
! go through the array in the reverse order placing lower triangular
!     elements in  appropriate slots.
      do j = n,1,-1
        i1 = ptr(j)
        i2 = ipkp1
        lenk = i2 - i1
! jstart is running pointer to position in new structure
        jstart = ckp1
! set ikp1 for next column
        ipkp1 = i1
        i2 = i2 - 1
! run through columns in reverse order
! lower triangular part of col. moved to end of same column in expanded form
        if (present(a)) then
          do ii = i2,i1,-1
            jstart = jstart - 1
            a(jstart) = a(ii)
            row(jstart) = row(ii)
          end do
        else
          do ii = i2,i1,-1
            jstart = jstart - 1
            row(jstart) = row(ii)
          end do
        end if
! ptr is set to position of first entry in lower triangular part of
!     column j in expanded form
        ptr(j) = jstart
! set ckp1 for next column
        ckp1 = ckp1 - iw(j)
! reset iw(j) to number of entries in lower triangle of column.
        iw(j) = lenk
      end do
 
! again sweep through the columns in the reverse order, this
!     time when one is handling column j the upper triangular
!     elements a(j,i) are put in position.
        do j = n,1,-1
          i1 = ptr(j)
          i2 = ptr(j) + iw(j) - 1
! run down column in order
! note that i is always greater than or equal to j
          if (present(a)) then
            r_sym_type = 0 ! symmetric
            if(present(sym_type)) r_sym_type = sym_type
            select case(r_sym_type)
            case(1) ! skew symmetric
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = -a(ii)
                row(ipos) = j
              end do
            case default ! symmetric or hermitian
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = a(ii)
                row(ipos) = j
              end do
            end select
          else
            do ii = i1,i2
              i = row(ii)
              if (i == j) cycle
              ptr(i) = ptr(i) - 1
              ipos = ptr(i)
              row(ipos) = j
            end do
          end if
        end do
      ptr(n+1) = newtau + 1

      end subroutine mc34_expand_single
   end module hsl_mc34_single
! COPYRIGHT (c) 2010 Science and Technology Facilities Council
! Original date 18 January 2011. Version 1.0.0
!
! Written by: Jonathan Hogg, John Reid, and Sue Thorne
!
! Version 1.4.2
! For version history, see ChangeLog
!
module hsl_mc69_single
   implicit none

   private
   public :: HSL_MATRIX_UNDEFINED,                             &
      HSL_MATRIX_REAL_RECT, HSL_MATRIX_CPLX_RECT,              &
      HSL_MATRIX_REAL_UNSYM, HSL_MATRIX_CPLX_UNSYM,            &
      HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_CPLX_HERM_PSDEF,   &
      HSL_MATRIX_REAL_SYM_INDEF, HSL_MATRIX_CPLX_HERM_INDEF,   &
      HSL_MATRIX_CPLX_SYM,                                     &
      HSL_MATRIX_REAL_SKEW, HSL_MATRIX_CPLX_SKEW
   public :: mc69_cscl_clean, mc69_verify, mc69_print, mc69_csclu_convert, &
      mc69_coord_convert, mc69_set_values, mc69_csrlu_convert, &
      mc69_cscl_convert, mc69_cscu_convert, mc69_csru_convert, &
      mc69_csrl_convert

   integer, parameter :: wp = kind(0.0)
   real(wp), parameter :: zero = 0.0_wp

   ! matrix types : real
   integer, parameter :: HSL_MATRIX_UNDEFINED      =  0 ! undefined/unknown
   integer, parameter :: HSL_MATRIX_REAL_RECT      =  1 ! real rectangular
   integer, parameter :: HSL_MATRIX_REAL_UNSYM     =  2 ! real unsymmetric
   integer, parameter :: HSL_MATRIX_REAL_SYM_PSDEF =  3 ! real symmetric pos def
   integer, parameter :: HSL_MATRIX_REAL_SYM_INDEF =  4 ! real symmetric indef
   integer, parameter :: HSL_MATRIX_REAL_SKEW      =  6 ! real skew symmetric

   ! matrix types : complex
   integer, parameter :: HSL_MATRIX_CPLX_RECT      = -1 ! complex rectangular
   integer, parameter :: HSL_MATRIX_CPLX_UNSYM     = -2 ! complex unsymmetric
   integer, parameter :: HSL_MATRIX_CPLX_HERM_PSDEF= -3 ! hermitian pos def
   integer, parameter :: HSL_MATRIX_CPLX_HERM_INDEF= -4 ! hermitian indef
   integer, parameter :: HSL_MATRIX_CPLX_SYM       = -5 ! complex symmetric
   integer, parameter :: HSL_MATRIX_CPLX_SKEW      = -6 ! complex skew symmetric

   ! Error flags
   integer, parameter :: MC69_SUCCESS                =  0
   integer, parameter :: MC69_ERROR_ALLOCATION       = -1
   integer, parameter :: MC69_ERROR_MATRIX_TYPE      = -2
   integer, parameter :: MC69_ERROR_N_OOR            = -3
   integer, parameter :: MC69_ERROR_M_NE_N           = -4
   integer, parameter :: MC69_ERROR_PTR_1            = -5
   integer, parameter :: MC69_ERROR_PTR_MONO         = -6
   integer, parameter :: MC69_ERROR_ROW_BAD_ORDER    = -7
   integer, parameter :: MC69_ERROR_ROW_OOR          = -8
   integer, parameter :: MC69_ERROR_ROW_DUP          = -9
   integer, parameter :: MC69_ERROR_ALL_OOR          = -10
   integer, parameter :: MC69_ERROR_MISSING_DIAGONAL = -11
   integer, parameter :: MC69_ERROR_IMAG_DIAGONAL    = -12
   integer, parameter :: MC69_ERROR_MISMATCH_LWRUPR  = -13
   integer, parameter :: MC69_ERROR_UPR_ENTRY        = -14
   integer, parameter :: MC69_ERROR_VAL_MISS         = -15
   integer, parameter :: MC69_ERROR_LMAP_MISS        = -16


   ! warning flags
   integer, parameter :: MC69_WARNING_IDX_OOR          = 1
   integer, parameter :: MC69_WARNING_DUP_IDX          = 2
   integer, parameter :: MC69_WARNING_DUP_AND_OOR      = 3
   integer, parameter :: MC69_WARNING_MISSING_DIAGONAL = 4
   integer, parameter :: MC69_WARNING_MISS_DIAG_OORDUP = 5

!            Possible error returns:

! MC69_ERROR_ALLOCATION         Allocation error
! MC69_ERROR_MATRIX_TYPE        Problem with matrix_type
! MC69_ERROR_N_OOR              n < 0 or m < 0
! MC69_ERROR_PTR_1              ptr(1) < 1
! MC69_ERROR_PTR_MONO           Error in ptr (not monotonic)
! MC69_ERROR_VAL_MISS           Only one of val_in and val_out is present
! MC69_ERROR_ALL_OOR            All the variables in a column are out-of-range
! MC69_ERROR_IMAG_DIAGONAL      Hermitian case and diagonal not real
! MC69_ERROR_ROW_BAD_ORDER      Entries within a column are not sorted by
!                               increasing row index 
! MC69_ERROR_MISMATCH_LWRUPR    Symmetric, skew symmetric or Hermitian: 
!                               entries in upper and lower
!                               triangles do not match
! MC69_ERROR_MISSING_DIAGONAL   Pos def and diagonal entries missing 
! MC69_ERROR_ROW_OOR            Row contains out-of-range entries      
! MC69_ERROR_ROW_DUP            Row contains duplicate entries         
! MC69_ERROR_M_NE_N             Square matrix and m .ne. n        

!           Possible warnings:

! MC69_WARNING_IDX_OOR          Out-of-range variable indices
! MC69_WARNING_DUP_IDX          Duplicated variable indices
! MC69_WARNING_DUP_AND_OOR      out of range and duplicated indices
! MC69_WARNING_MISSING_DIAGONAL Indefinite case and diagonal entries missing
! MC69_WARNING_MISS_DIAG_OORDUP As MC69_WARNING_MISSING_DIAGONAL, and 
!                               out-of-range and/or duplicates


!!!!!!!!!!!!!!!!!!!!!!!!
! Internal types
!!!!!!!!!!!!!!!!!!!!!!!!
type dup_list
   integer :: src
   integer :: dest
   type(dup_list), pointer :: next => null()
end type dup_list

interface mc69_verify
   module procedure mc69_verify_single
end interface mc69_verify
interface mc69_print
   module procedure mc69_print_single
end interface
interface mc69_cscl_clean
   module procedure mc69_cscl_clean_single
end interface mc69_cscl_clean
interface mc69_set_values
   module procedure mc69_set_values_single
end interface mc69_set_values
interface mc69_cscl_convert
   module procedure mc69_cscl_convert_single
end interface mc69_cscl_convert
interface mc69_cscu_convert
   module procedure mc69_cscu_convert_single
end interface mc69_cscu_convert
interface mc69_csclu_convert
   module procedure mc69_csclu_convert_single
end interface mc69_csclu_convert
interface mc69_csrl_convert
   module procedure mc69_csrl_convert_single
end interface mc69_csrl_convert
interface mc69_csru_convert
   module procedure mc69_csru_convert_single
end interface mc69_csru_convert
interface mc69_csrlu_convert
   module procedure mc69_csrlu_convert_single
end interface mc69_csrlu_convert
interface mc69_coord_convert
   module procedure mc69_coord_convert_single
end interface mc69_coord_convert

contains

!
! To verify that a matrix is in HSL format, or identify why it is not
!
subroutine mc69_verify_single(lp, matrix_type, m, n, ptr, row, flag, more, val)
   integer, intent(in) :: lp ! output unit
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr ! column starts
   integer, dimension(*), intent(in) :: row ! row indices.
     ! Entries within each column must be sorted in order of 
     ! increasing row index. no duplicates and/or out of range entries allowed.
   integer, intent(out) :: flag ! return code
   integer, intent(out) :: more ! futher error information (or set to 0)
   real(wp), dimension(:), optional, intent(in) :: val ! matrix values,if any

   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   logical :: diag ! flag for detection of diagonal
   integer :: i
   integer :: j
   integer :: k
   integer :: last ! last row index
   logical :: lwronly
   integer, dimension(:), allocatable :: ptr2
   integer :: st

   context = 'mc69_verify'
   flag = MC69_SUCCESS

   more = 0 ! ensure more is not undefined.

   ! Check matrix_type
   select case(matrix_type)
   case(0)
      ! Undefined matrix. OK, do nothing
   case(1:4,6)
      ! Real matrix. OK, do nothing
   case(-6:-1)
      ! Complex matrix. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   case default
      ! Out of range value. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   end select

   ! Check m and n are valid; skip remaining tests if n=0
   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,lp,flag)
      return
   end if
   if(abs(matrix_type).ne.HSL_MATRIX_REAL_RECT .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,lp,flag)
      return
   endif
   if(n == 0 .or. m == 0) return

   ! Check ptr(1) is valid
   if(ptr(1) < 1) then
      more = ptr(1)
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,lp,flag)
      return
   end if

   ! Check ptr is monotonically increasing
   do i = 1, n
      if(ptr(i+1).ge.ptr(i)) cycle ! ptr is monotonically increasing, good.
      flag = MC69_ERROR_PTR_MONO
      more = i + 1
      call mc69_print_flag(context,lp,flag)
      return
   end do

   ! Count number of entries in each row. Also:
   ! * Check ordering of entries
   ! * Check for out-of-range entries
   ! * Check for duplicate entries
   ! * Lack of diagonal entries in real pos. def. case.

   ! ptr2(k+2) is set to hold the number of entries in row k
   allocate(ptr2(m+2), stat=st)
   if(st.ne.0) goto 100
   ptr2(:) = 0

   lwronly = (abs(matrix_type).ne.HSL_MATRIX_UNDEFINED) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_RECT) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_UNSYM)
   do col = 1, n
      last = -1
      diag = .false.
      do j = ptr(col), ptr(col+1)-1
         k = row(j)
         ! Check out-of-range
         if(k.lt.1 .or. k.gt.m) then
            flag = MC69_ERROR_ROW_OOR
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(lwronly .and. k.lt.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. k.eq.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check duplicate
         if(k.eq.last) then
            flag = MC69_ERROR_ROW_DUP
            more = j-1
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check order
         if(k.lt.last) then
            flag = MC69_ERROR_ROW_BAD_ORDER
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check for diagonal
         diag = diag .or. (k.eq.col)
         ! Increase count for row k
         ptr2(k+2) = ptr2(k+2) + 1
         ! Store value for next loop
         last = k
      end do
      ! If marked as positive definite, check if diagonal was present
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF .and. .not.diag) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         more = col
         call mc69_print_flag(context,lp,flag)
         return
      endif
   end do

   if(present(val)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr(j)
            ! Note: column cannot be empty as previously checked pattern ok
            if(real(val(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               more = j
               call mc69_print_flag(context,lp,flag)
               return
            end if 
         end do
      end select
   endif

   return

   100 continue ! Allocation error
   flag = MC69_ERROR_ALLOCATION
   more = st
   call mc69_print_flag(context,lp,flag)
   return

end subroutine mc69_verify_single

!****************************************

!
! Pretty prints a matrix as best it can
!
subroutine mc69_print_single(lp, lines, matrix_type, m, n, ptr, row, val, cbase)
   integer, intent(in) :: lp ! unit to print on
   integer, intent(in) :: lines ! max number of lines to use (ignored if -ive)
   integer, intent(in) :: matrix_type ! type of matrix
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of cols in matrix
   integer, dimension(n+1), intent(in) :: ptr ! column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices
   real(wp), dimension(ptr(n+1)-1), optional, intent(in) :: val ! matrix vals
   logical, optional, intent(in) :: cbase ! is true, rebase for C

   integer :: col, j, k
   integer :: llines
   integer, dimension(:,:), allocatable :: dmat
   character(len=5) :: mfrmt, nfrmt, nefrmt
   character(len=12) :: negfrmt, valfrmt, emptyfrmt
   integer ::  rebase

   if(lp.lt.0) return ! invalid unit

   ! Check if we need to rebase for C consistent output
   rebase = 0
   if(present(cbase)) then
      if(cbase) rebase = 1
   endif

   ! Calculate number of lines to play with
   llines = huge(llines)
   if(lines.gt.0) llines = lines

   ! Print a summary statement about the matrix
   mfrmt = digit_format(m)
   nfrmt = digit_format(n)
   nefrmt = digit_format(ptr(n+1)-1)

   select case(matrix_type)
   case(HSL_MATRIX_UNDEFINED)
      write(lp, "(a)", advance="no") &
         "Matrix of undefined type, dimension "
   case(HSL_MATRIX_REAL_RECT)
      write(lp, "(a)", advance="no") &
         "Real rectangular matrix, dimension "
   case(HSL_MATRIX_REAL_UNSYM)
      write(lp, "(a)", advance="no") &
         "Real unsymmetric matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric positive definite matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_INDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric indefinite matrix, dimension "
   case(HSL_MATRIX_REAL_SKEW)
      write(lp, "(a)", advance="no") &
         "Real skew symmetric matrix, dimension "
   case default
      write(lp, "(a,i5)") &
         "Unrecognised matrix_type = ", matrix_type
      return
   end select
   write(lp, mfrmt, advance="no") m
   write(lp, "(a)", advance="no") "x"
   write(lp, nfrmt, advance="no") n
   write(lp, "(a)", advance="no") " with "
   write(lp, nefrmt, advance="no") ptr(n+1)-1
   write(lp, "(a)") " entries."

   ! immediate return if m = 0 or n = 0
   if (m == 0 .or. n == 0) return

   if(((present(val) .and. n.lt.10) .or. (.not.present(val) .and. n.lt.24)) &
         .and. m+1.le.llines) then
      ! Print the matrix as if dense
      allocate(dmat(m, n))
      dmat(:,:) = 0
      do col = 1, n
         do j = ptr(col), ptr(col+1) - 1
            k = row(j)
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) then
               dmat(col, k) = -j
            endif
            dmat(k, col) = j
         end do
      end do

      select case(n)
      case(:6)
         valfrmt = "(1x,es12.4)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,a12)"
      case(7)
         valfrmt = "(1x,es10.2)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,a10)"
      case(8:)
         valfrmt = "(1x,es8.2)"
         negfrmt = "(1x,es8.1)"
         emptyfrmt = "(1x,a8)"
      end select

      do k = 1, m
         write(lp,mfrmt,advance="no") k-rebase
         write(lp,"(':')",advance="no")
         if(present(val)) then
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,emptyfrmt,advance="no") '                '
               elseif(dmat(k,j).gt.0) then
                  if(val(dmat(k,j)).gt.zero) then
                     write(lp,valfrmt,advance="no") val(dmat(k,j))
                  else
                     write(lp,negfrmt,advance="no") val(dmat(k,j))
                  endif
               else
                  ! in upper triangle
                  select case(matrix_type)
                  case(HSL_MATRIX_REAL_SYM_INDEF,HSL_MATRIX_REAL_SYM_PSDEF)
                     if(val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") val(-dmat(k,j))
                     endif
                  case(HSL_MATRIX_REAL_SKEW)
                     if(-val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") -val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") -val(-dmat(k,j))
                     endif
                  end select
               endif
            end do
         else ! pattern only
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,"('  ')",advance="no")
               else
                  write(lp,"(1x,'x')",advance="no")
               endif
            end do
         end if
         write(lp,"()")
      end do
   else
      ! Output first x entries from each column
      llines = llines - 1 ! Discount first info line
      if(llines.le.2) return
      write(lp, "(a)") "First 4 entries in columns:"
      llines = llines - 1
      do col = 1, llines
         write(lp, "(a)", advance="no") "Col "
         write(lp, nfrmt, advance="no") col-rebase
         write(lp, "(':')", advance="no")
         do j = ptr(col), min(ptr(col+1)-1, ptr(col)+3)
            write(lp, "('  ')", advance="no")
            write(lp, mfrmt, advance="no") row(j)-rebase
            if(present(val)) then
               write(lp, "(' (',es12.4,')')", advance="no") val(j)
            endif
         end do
         write(lp, "()")
      end do
   endif
end subroutine mc69_print_single

!****************************************

character(len=5) function digit_format(x)
   integer, intent(in) :: x

   integer :: ndigit

   ndigit = int(log10(real(x))) + 1
   if(ndigit<10) then
      write(digit_format,"('(i',i1,')')") ndigit
   else
      write(digit_format,"('(i',i2,')')") ndigit
   endif
end function digit_format

!****************************************

! Convert a matrix held in lower compressed sparse column format to standard
! HSL format in-place. Duplicates are summed and out-of-range entries are
! remove.  For symmetric, skew-symmetric and Hermitian matrices, entries in 
! the upper triangle are removed and discarded.
!
! Note: this routine is in place! cscl_convert is the out of place version.
subroutine mc69_cscl_clean_single(matrix_type, m, n, ptr, row, flag, &
      val, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m, n ! matrix dimensions
   integer, dimension(n+1), intent(inout) :: ptr ! column pointers
   integer, dimension(*), intent(inout) :: row ! row indices
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(inout) :: val ! values
   integer, optional, intent(out) :: lmap ! size of map
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means val_out(i)=val(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag ! number of columns with a diagonal entry
   integer :: idup ! number of duplicates summed
   integer :: ioor ! number of out-of-range entries
   integer :: j
   integer :: k
   integer :: k1
   integer :: k2
   integer :: l
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer, allocatable :: temp(:) ! work array

   context = 'mc69_cscl_clean'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type) > 1 .and. m /= n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if
 
   if(ptr(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate work array

   allocate(temp(m),stat=st)
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,lp,flag)
      return
   endif
   temp(:) = 0

   ! count duplicates and out-of-range indices
   idup = 0; ioor = 0; idiag = 0
   k = 0 ! last diagonal seen
   l = 1
   do col = 1, n
      if(ptr(col+1).lt.ptr(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      if(abs(matrix_type) > 2) l = col
      do i = ptr(col), ptr(col+1)-1
         j = row(i)
         if(j<l .or. j>m) then
            ioor = ioor + 1
            row(i) = m+1
         else
            if(temp(j)==col) idup = idup + 1
            temp(j) = col
         end if
         if(j.eq.col .and. k<col) then
            idiag = idiag + 1
            k = col
         endif
      end do
   end do
   if(present(ndup)) ndup = idup
   if(present(noor)) noor = ioor
   deallocate(temp,stat=st)
   

   k = 0
   l = ptr(n+1) 
   if(present(map)) then
      deallocate(map,stat=st)
      lmap = l - 1 - ioor + idup
      allocate(map(l-1+idup*2),stat=st)
      if(st /= 0) then
         flag = MC69_ERROR_ALLOCATION
         call mc69_print_flag(context,lp,flag)
         return
      endif
      do i = 1,ptr(n+1)-1
         map(i) = i
      end do
      if(present(val)) then
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:), val=val(k1))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
               val(k) = val(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
                  val(k) = val(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  val(k) = val(k) + val(i)
                  l = l + 2
               end if
            end do
         end do
      else 
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  l = l + 2
               end if
            end do
         end do
      end if
      l = ptr(n+1) - 1
      ptr(n+1) = k + 1
      ! Move duplicate pair in map forward
      do i = 1, idup*2
         map(k+i) = map(l+i)
      end do
   else if(present(val)) then
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1, val=val(k1))
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
            val(k) = val(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
               val(k) = val(i)
            else
               val(k) = val(k) + val(i)
               l = l + 2
            end if
         end do
      end do
      ptr(n+1) = k + 1
   else 
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1)
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
            end if
         end do
      end do
      ptr(n+1) = k + 1
   end if

   select case(matrix_type)
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      ! Check for positive diagonal entries
      do j = 1,n
         k = ptr(j)
         ! Note: we're OK if the column is empty - row(k) is still a diagonal
         ! entry, however we must be careful that we've not gone past the
         ! end of the matrix
         if(k.ge.ptr(n+1)) exit
         if(row(k)/=j)then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
         if(present(val))then
            if(val(k)<=zero)then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end if 
      end do
   end select

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

end subroutine mc69_cscl_clean_single

!****************************************

!
! Converts CSC with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_cscl_convert_single(matrix_type, m, n, ptr_in, row_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscl_convert_single

!****************************************

!
! Converts CSC (with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices) to HSL standard format.
! Also used for symmetric, skew-symmetric and Hermitian matrices in upper
! CSR format
!
subroutine mc69_cscl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 or 1, differs for csc/csr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer :: minidx

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   nullify(dup, duphead)

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! ensure output arrays are not allocated

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)

   idup = 0; ioor = 0; idiag = 0

   allocate(row_out(ptr_in(n+1)-1),stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Allocate map for worst case where all bar one are duplicates
      allocate(map(2*ptr_in(n+1)-2),stat=st)
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            map(k) = multiplier*i
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, map=map(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               call cleanup_dup(duphead)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  allocate(dup,stat=st)
                  if(st.ne.0) goto 100
                  dup%next => duphead
                  duphead => dup
                  dup%src = map(i)
                  dup%dest = k-1
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               map(k) = map(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      allocate(val_out(ptr_in(n+1)-1),stat=st)
      if(st.ne.0) goto 100
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         select case(matrix_type)
         case(HSL_MATRIX_REAL_SKEW)
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = multiplier*val_in(i)
               k = k + 1
            end do
         case default
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = val_in(i)
               k = k + 1
            end do
         end select
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, &
               val=val_out(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, sum then drop from pattern
                  idup = idup + 1
                  val_out(i-1) = val_out(i-1) + val_out(i)
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               val_out(k) = val_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         endif
      end do
      ptr_out(n+1) = k
   else ! pattern only
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i)
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  cycle
               endif
               if(row_out(i).eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_cscl_convert_main

!****************************************

!
! Converts CSC with only upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_cscu_convert_single(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: row_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscu_convert_single

!****************************************

!
! Converts CSC with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csclu_convert_single(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csclu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, -1, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csclu_convert_single

!****************************************

subroutine mc69_csclu_convert_main(context, multiplier, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 for upr source, +1 for lwr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nlwr ! number of strictly lower triangular entries
   integer :: npre
   integer :: nupr ! number of strictly upper triangular entries
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0
   nlwr = 0; nupr = 0

   !
   ! First pass, count number of entries in each row of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   do col = 1, n
      if(ptr_in(col+1).lt.ptr_in(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = nlwr + nupr + idiag
      if(ptr_in(col+1).eq.ptr_in(col)) cycle ! no entries in column
      do i = ptr_in(col), ptr_in(col+1)-1
         j = row_in(i)
         if(j.lt.1 .or. j.gt.n) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) then
            ioor = ioor + 1
            cycle
         endif
         if(j.gt.col) then
            nlwr = nlwr + 1
            cycle
         endif
         if(j.lt.col) then
            nupr = nupr + 1
         elseif(j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = col
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(nlwr+nupr+idiag.eq.npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Check number in lower triangle = number in upper triangle
   if(nlwr .ne. nupr) then
      flag = MC69_ERROR_MISMATCH_LWRUPR
      call mc69_print_flag(context,nout,flag)
      return
   endif

   ! Determine column starts for transposed matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle ! ignore oor and lwr triangle entries
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.ge.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.gt.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csclu_convert_main

!****************************************

!
! Converts CSR with only lower entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrl_convert_single(matrix_type, m, n, ptr_in, col_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup)
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 1 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrl_convert_single

!****************************************

subroutine mc69_csrl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! 1 for lwr triangle source, -1 for upr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: row ! current row
   integer :: col ! current col
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: npre ! number of out-of-range entries prior to this column
   integer :: st ! stat parameter
   integer :: maxv ! maximum value before out of range

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   maxv = n
   do row = 1, m
      if(ptr_in(row+1).lt.ptr_in(row)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = ioor
      if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
      do i = ptr_in(row), ptr_in(row+1)-1
         j = col_in(i)
         if(j.lt.1 .or. j.gt.maxv) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.row .and. j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = row
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(ioor-npre.ne.0 .and. ptr_in(row+1)-ptr_in(row).eq.ioor-npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Determine column starts for matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle ! ignore oor and upr triangle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do row = 1, m
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.ge.row) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case default
         maxv= n
         do row = 1, m
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
            if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.gt.maxv) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csrl_convert_main

!****************************************

!
! Converts CSR with upper entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_csru_convert_single(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: col_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if lp not present)

   context = 'mc69_csru_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csru_convert_single

!****************************************

!
! Converts CSR with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrlu_convert_single(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrlu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, 1, matrix_type, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrlu_convert_single

!****************************************

!
! Converts COOR format to CSC with only lower entries present for
! (skew-)symmetric problems. Entries within each column ordered by increasing
! row index (HSL standard format)
!
subroutine mc69_coord_convert_single(matrix_type, m, n, ne, row, col, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of columns in matrix
   integer, intent(in) :: ne ! number of input nonzero entries
   integer, dimension(:), intent(in) :: row(ne) ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(in) :: col(ne) ! column indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(out) :: ptr_out(n+1) ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means
      ! val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k, l1, l2
   integer :: l
   integer :: ne_new
   integer :: nout ! output unit (set to -1 if lp not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   context = 'mc69_coord_convert'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type).ge.HSL_MATRIX_REAL_UNSYM .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix
   ! matrix. Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   do l = 1, ne
      i = row(l)
      j = col(l)
      if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) then
         ioor = ioor + 1
         cycle
      endif

      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. i.eq.j) then
         ioor = ioor + 1
         cycle
      endif
   
      select case (abs(matrix_type))
      case (HSL_MATRIX_REAL_SYM_PSDEF:)
         if(i.ge.j) then
            ptr_out(j+1) = ptr_out(j+1) + 1
         else
            ptr_out(i+1) = ptr_out(i+1) + 1
         end if
      case default
          ptr_out(j+1) = ptr_out(j+1) + 1
      end select
   end do


   ! Determine column starts for transposed expanded matrix such 
   ! that column i starts at ptr_out(i)
   ne_new = 0
   ptr_out(1) = 1
   do i = 2, n+1
      ne_new = ne_new + ptr_out(i)
      ptr_out(i) = ptr_out(i) + ptr_out(i-1)
   end do

   ! Check whether all entries out of range
   if(ne.gt.0 .and. ne_new.eq.0) then
      flag = MC69_ERROR_ALL_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! expanded matrix
   !
   allocate(row_out(ne_new), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      if(allocated(map)) deallocate(map,stat=st)
      allocate(map(2*ne), stat=st)
      if(st.ne.0) goto 100
      map(:) = 0
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = -l
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = l
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            map(k) = l
         end do
      end select
   elseif(present(val_out)) then
      allocate(val_out(ne_new), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = -val_in(l)
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = val_in(l)
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            val_out(k) = val_in(l)
         end do
      end select


   else
      ! neither val_out or map present
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
          do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
         end do
      end select
   endif

   do j=n,2,-1
      ptr_out(j) = ptr_out(j-1)
   end do
   ptr_out(1) = 1

   !
   ! Third pass, in place sort and removal of duplicates
   ! Also checks for diagonal entries in pos. def. case.
   ! Uses a modified insertion sort for best speed on already ordered data
   !
   idup=0
   if(present(map)) then
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.gt.1) call sort ( row_out(l1:l2),l, map=map(l1:l2) )
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         map(k) = map(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   else if(present(val_out)) then
      ! ADD
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1) call sort( row_out(l1:l2),l,val=val_out(l1:l2) )
         ! ADD
      end do 

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         val_out(k) = val_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              val_out(k-1) = val_out(k-1)+val_out(i)
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1)call sort ( row_out(l1:l2),l)
         ! ADD
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   
   endif


   ! Check for missing diagonals in pos def and indef cases
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
      do j = 1,n
         if(ptr_out(j).lt.ptr_out(n+1)) then
            if(row_out(ptr_out(j)) .ne. j) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if
         end if
      end do
   end if

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
      if(ioor > 0) flag = MC69_WARNING_IDX_OOR
      if(idup > 0) flag = MC69_WARNING_DUP_IDX
      if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. &
            idiag<n .and. ioor>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n .and.&
            idup>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n) then
         flag = MC69_WARNING_MISSING_DIAGONAL
      end if
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup
   return

   100 continue
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
   end if
   return

end subroutine mc69_coord_convert_single

!*************************************************

!
! This subroutine will use map to translate the values of val to val_out
!
subroutine mc69_set_values_single(matrix_type, lmap, map, val, ne, val_out)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: lmap
   integer, dimension(lmap), intent(in) :: map
   real(wp), dimension(*), intent(in) :: val
   integer, intent(in) :: ne
   real(wp), dimension(ne), intent(out) :: val_out

   integer :: i, j, k

   select case(matrix_type)
   case default
      !
      ! Rectangular, Unsymmetric or Symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + val(k)
      end do
   case(HSL_MATRIX_REAL_SKEW)
      !
      ! Skew symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = sign(1.0,real(map(i)))*val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + sign(1.0,real(map(i+1)))*val(k)
      end do
   end select
end subroutine mc69_set_values_single

!*************************************************

subroutine mc69_print_flag(context,nout,flag)
   integer, intent (in) :: flag, nout
   character (len=*), optional, intent(in) :: context

   if(nout < 0) return
   if(flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context), &
         '. Error flag = ', flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context), &
         '. Warning flag = ', flag
   end if

   select case(flag)
   !
   ! Errors
   !
   case(MC69_ERROR_ALLOCATION)
      write (nout,'(a)') ' Allocation error'
   case(MC69_ERROR_MATRIX_TYPE)
       write (nout,'(a)') ' matrix_type has invalid value'
   case(MC69_ERROR_N_OOR)
      write (nout,'(a)') ' m or n is out-of-range'
   case(MC69_ERROR_ALL_OOR)
      write (nout,'(a)') ' All entries in a column out-of-range'
   case(MC69_ERROR_PTR_MONO)
      write (nout,'(a)') ' ptr not monotonic'
   case(MC69_ERROR_PTR_1)
      write (nout,'(a)') ' ptr(1) < 1'
   case(MC69_ERROR_IMAG_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is not real'
   case(MC69_ERROR_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries are not positive'
   case(MC69_ERROR_VAL_MISS)
      write (nout,'(a)') ' Only one of val and val_out is present'
   case(MC69_ERROR_LMAP_MISS)
      write (nout,'(a)') ' Only one of lmap and map is present'
   case(MC69_ERROR_UPR_ENTRY)
      write (nout,'(a)') ' Entry in upper triangle'
   case(MC69_ERROR_M_NE_N)
      write (nout,'(a)') ' m is not equal to n'
   !
   ! Warnings
   !
   case(MC69_WARNING_IDX_OOR)
      write (nout,'(a)') ' out-of-range indices detected'
   case(MC69_WARNING_DUP_IDX)
      write (nout,'(a)') ' duplicate entries detected'
   case(MC69_WARNING_DUP_AND_OOR)
      write (nout,'(a)') &
         ' out-of-range indices detected and duplicate entries detected'
   case(MC69_WARNING_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is missing'
   case(MC69_WARNING_MISS_DIAG_OORDUP)
      write (nout,'(a)') ' one or more diagonal entries is missing and'
      write (nout,'(a)') ' out-of-range and/or duplicate entries detected'
   end select

end subroutine mc69_print_flag

!************************************************************************

!
!   Sort an integer array by heapsort into ascending order.
!
subroutine sort( array, n, map, val )
   integer, intent(in) :: n       ! Size of array to be sorted
   integer, dimension(n), intent(inout) :: array ! Array to be sorted
   integer, dimension(n), optional, intent(inout) :: map
   real(wp), dimension(n), optional, intent(inout) :: val ! Apply same
      ! permutation to val

   integer :: i
   integer :: temp
   real(wp) :: vtemp
   integer :: root

   if(n.le.1) return ! nothing to do

   !
   ! Turn array into a heap with largest element on top (as this will be pushed
   ! on the end of the array in the next phase)
   !
   ! We start at the bottom of the heap (i.e. 3 above) and work our way
   ! upwards ensuring the new "root" of each subtree is in the correct
   ! location
   root = n / 2
   do root = root, 1, -1
      call pushdown(root, n, array, val=val, map=map)
   end do

   !
   ! Repeatedly take the largest value and swap it to the back of the array
   ! then repeat above code to sort the array
   !
   do i = n, 2, -1
      ! swap a(i) and head of heap a(1)
      temp = array(1)
      array(1) = array(i)
      array(i) = temp
      if(present(val)) then
         vtemp = val(1)
         val(1) = val(i)
         val(i) = vtemp
      endif
      if(present(map)) then
         temp = map(1)
         map(1) = map(i)
         map(i) = temp
      endif
      call pushdown(1,i-1, array, val=val, map=map)
   end do
end subroutine sort

!****************************************

! This subroutine will assume everything below head is a heap and will
! push head downwards into the correct location for it
subroutine pushdown(root, last, array, val, map)
   integer, intent(in) :: root
   integer, intent(in) :: last
   integer, dimension(last), intent(inout) :: array
   real(wp), dimension(last), optional, intent(inout) :: val
   integer, dimension(last), optional, intent(inout) :: map

   integer :: insert ! current insert position
   integer :: test ! current position to test
   integer :: root_idx ! value of array(root) at start of iteration
   real(wp) :: root_val ! value of val(root) at start of iteration
   integer :: root_map ! value of map(root) at start of iteration

   ! NB a heap is a (partial) binary tree with the property that given a
   ! parent and a child, array(child)>=array(parent).
   ! If we label as
   !                      1
   !                    /   \
   !                   2     3
   !                  / \   / \
   !                 4   5 6   7
   ! Then node i has nodes 2*i and 2*i+1 as its children

   if(present(val) .and. present(map)) then ! both val and map
      root_idx = array(root)
      root_val = val(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test);
         val(insert) = val(test);
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
      map(insert) = root_map
   elseif(present(val)) then ! val only, not map
      root_idx = array(root)
      root_val = val(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         val(insert) = val(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
   elseif(present(map)) then ! map only, not val
      root_idx = array(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating mapue up
         array(insert) = array(test)
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root mapue into location found
      array(insert) = root_idx
      map(insert) = root_map
   else ! neither map nor val
      root_idx = array(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
   endif

end subroutine pushdown

!****************************************

subroutine cleanup_dup(duphead)
   type(dup_list), pointer :: duphead ! NB: can't have both intent() and pointer

   type(dup_list), pointer :: dup

   do while(associated(duphead))
      dup => duphead%next
      deallocate(duphead)
      duphead => dup
   end do
end subroutine cleanup_dup

end module hsl_mc69_single
