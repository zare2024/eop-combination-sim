! lsq.f90
! Least-squares adjustment (Gauss-Markov model) through normal equations.
!
!   observation equations   l + v = A x        weights P = diag(w)
!   normal equations        N x = b            N = A^T P A,   b = A^T P l
!   solution                x = N^-1 b         cofactor matrix Q = N^-1
!   residual sum            v^T P v = l^T P l - x^T b
!
! A system is kept as N, b and l^T P l rather than as A. That is the form in
! which systems of different techniques can later be stacked (combination on
! the normal equation level).

module lsq
  use constants, only: dp
  implicit none
  private
  public :: normal_eq_t, neq_init, neq_add, neq_solve, neq_vtpv

  type :: normal_eq_t
    integer  :: n_par = 0, n_obs = 0
    real(dp), allocatable :: n(:, :)     ! normal matrix A^T P A
    real(dp), allocatable :: b(:)        ! right-hand side A^T P l
    real(dp) :: ltpl = 0.0_dp            ! l^T P l
  end type normal_eq_t

  ! LAPACK routines for symmetric positive definite matrices (Cholesky)
  interface
    subroutine dpotrf(uplo, n, a, lda, info)
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n, lda
      double precision, intent(inout) :: a(lda, *)
      integer,          intent(out)   :: info
    end subroutine dpotrf
    subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n, nrhs, lda, ldb
      double precision, intent(in)    :: a(lda, *)
      double precision, intent(inout) :: b(ldb, *)
      integer,          intent(out)   :: info
    end subroutine dpotrs
    subroutine dpotri(uplo, n, a, lda, info)
      character,        intent(in)    :: uplo
      integer,          intent(in)    :: n, lda
      double precision, intent(inout) :: a(lda, *)
      integer,          intent(out)   :: info
    end subroutine dpotri
  end interface

contains

  subroutine neq_init(neq, n_par)
    type(normal_eq_t), intent(out) :: neq
    integer,           intent(in)  :: n_par
    neq%n_par = n_par
    neq%n_obs = 0
    neq%ltpl  = 0.0_dp
    allocate(neq%n(n_par, n_par), neq%b(n_par))
    neq%n = 0.0_dp
    neq%b = 0.0_dp
  end subroutine neq_init

  ! Add one observation: partials a, reduced observation l (observed minus computed), weight w
  subroutine neq_add(neq, a, l, w)
    type(normal_eq_t), intent(inout) :: neq
    real(dp),          intent(in)    :: a(:), l, w
    integer :: i, j
    do j = 1, neq%n_par
      do i = 1, neq%n_par
        neq%n(i, j) = neq%n(i, j) + w * a(i) * a(j)
      end do
    end do
    neq%b     = neq%b + w * a * l
    neq%ltpl  = neq%ltpl + w * l * l
    neq%n_obs = neq%n_obs + 1
  end subroutine neq_add

  ! Solve N x = b by Cholesky decomposition N = U^T U; q receives N^-1.
  ! info > 0 means N is not positive definite: some parameter is not estimable.
  subroutine neq_solve(neq, x, q, info)
    type(normal_eq_t), intent(in)  :: neq
    real(dp),          intent(out) :: x(:), q(:, :)
    integer,           intent(out) :: info
    integer :: i, j, n

    n = neq%n_par
    q = neq%n
    call dpotrf('U', n, q, n, info)          ! q now holds U
    if (info /= 0) return
    x = neq%b
    call dpotrs('U', n, 1, q, n, x, n, info) ! forward and back substitution
    if (info /= 0) return
    call dpotri('U', n, q, n, info)          ! N^-1 from U, upper triangle only
    do j = 1, n
      do i = j + 1, n
        q(i, j) = q(j, i)
      end do
    end do
  end subroutine neq_solve

  ! Weighted sum of squared residuals for the solution x, without forming v
  pure function neq_vtpv(neq, x) result(vtpv)
    type(normal_eq_t), intent(in) :: neq
    real(dp),          intent(in) :: x(:)
    real(dp) :: vtpv
    vtpv = neq%ltpl - dot_product(x, neq%b)
  end function neq_vtpv

end module lsq
