module MMA_Module
    use MMA_Routines
    implicit none

    interface MMA_Solver_Interface
        module procedure MMA_Solver_Real
        module procedure MMA_Solver_DP
    end interface MMA_Solver_Interface

contains

    ! interface of MMA 
    subroutine MMA_Solver_Real(m,n,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d)
        implicit none
        ! input or output
        integer, intent(in)                                         :: m            ! 
        integer, intent(in)                                         :: n            !
        integer, intent(in)                                         :: loop         !
        real, dimension(:,:), allocatable, intent(inout)            :: xval         !
        real, dimension(:,:), allocatable, intent(inout)            :: xmin         !
        real, dimension(:,:), allocatable, intent(inout)            :: xmax         !
        real, dimension(:,:), allocatable, intent(inout)            :: xold1        !
        real, dimension(:,:), allocatable, intent(inout)            :: xold2        !
        real, intent(inout)                                         :: f0val        !
        real, dimension(:,:), allocatable, intent(inout)            :: df0dx        !
        real, dimension(:,:), allocatable, intent(inout)            :: fval         !
        real, dimension(:,:), allocatable, intent(inout)            :: dfdx         !
        real, dimension(:,:), allocatable, intent(inout)            :: low          !
        real, dimension(:,:), allocatable, intent(inout)            :: upp          !
        real                                                        :: a0           !
        real, dimension(:,:), allocatable, intent(inout)            :: a            !
        real, dimension(:,:), allocatable, intent(inout)            :: c            !
        real, dimension(:,:), allocatable, intent(inout)            :: d            !
        ! internal variables
        real, dimension(:,:), allocatable                           :: xmma         ! 
        real, dimension(:,:), allocatable                           :: ymma         ! 
        real                                                        :: zmma         ! 
        real, dimension(:,:), allocatable                           :: lam          ! 
        real, dimension(:,:), allocatable                           :: xsi          ! 
        real, dimension(:,:), allocatable                           :: eta          ! 
        real, dimension(:,:), allocatable                           :: mu           ! 
        real                                                        :: zet          ! 
        real, dimension(:,:), allocatable                           :: s            !
        !---------------------------------------------------------------------------!

        ! MMA subproblem is solved
        call MMA_Sub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp, &
                    m,n,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,a0,a,c,d)
        ! Updating of some vectors
        xval  = xmma
    end subroutine MMA_Solver_Real

    subroutine MMA_Solver_DP(m,n,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d)
        implicit none
        ! input or output
        integer, intent(in)                                          :: m            ! 
        integer, intent(in)                                          :: n            !
        integer, intent(in)                                          :: loop         !
        double precision, dimension(:,:), allocatable, intent(inout) :: xval         !
        double precision, dimension(:,:), allocatable, intent(inout) :: xmin         !
        double precision, dimension(:,:), allocatable, intent(inout) :: xmax         !
        double precision, dimension(:,:), allocatable, intent(inout) :: xold1        !
        double precision, dimension(:,:), allocatable, intent(inout) :: xold2        !
        double precision, intent(inout)                              :: f0val        !
        double precision, dimension(:,:), allocatable, intent(inout) :: df0dx        !
        double precision, dimension(:,:), allocatable, intent(inout) :: fval         !
        double precision, dimension(:,:), allocatable, intent(inout) :: dfdx         !
        double precision, dimension(:,:), allocatable, intent(inout) :: low          !
        double precision, dimension(:,:), allocatable, intent(inout) :: upp          !
        double precision                                             :: a0           !
        double precision, dimension(:,:), allocatable, intent(inout) :: a            !
        double precision, dimension(:,:), allocatable, intent(inout) :: c            !
        double precision, dimension(:,:), allocatable, intent(inout) :: d            !
        ! internal variables
        double precision, dimension(:,:), allocatable                :: xmma         ! 
        double precision, dimension(:,:), allocatable                :: ymma         ! 
        double precision                                             :: zmma         ! 
        double precision, dimension(:,:), allocatable                :: lam          ! 
        double precision, dimension(:,:), allocatable                :: xsi          ! 
        double precision, dimension(:,:), allocatable                :: eta          ! 
        double precision, dimension(:,:), allocatable                :: mu           ! 
        double precision                                             :: zet          ! 
        double precision, dimension(:,:), allocatable                :: s            !
        !----------------------------------------------------------------------------!

        ! MMA subproblem is solved
        call MMA_Sub(xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp, &
                     m,n,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,a0,a,c,d)
        xval  = xmma
    end subroutine MMA_Solver_DP  
end module MMA_Module