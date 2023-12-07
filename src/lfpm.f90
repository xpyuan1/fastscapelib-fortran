! -------------------------- functions ------------------------------

module lfpm_functions
contains

  function m_inv (matrix) ! (m11,m12,m21,m22)
  implicit none
  real(kind=8), dimension(1:4) :: matrix
  real(kind=8), dimension(1:4) :: m_inv
  real(kind=8) :: det
  det=matrix(1)*matrix(4)-matrix(2)*matrix(3)
  m_inv=(/ matrix(4)/det,-matrix(2)/det,-matrix(3)/det,matrix(1)/det /)
  end function

  function m_mul_m (matrix1, matrix2) ! (m11,m12,m21,m22)
  implicit none
  real(kind=8), dimension(1:4) :: matrix1,matrix2
  real(kind=8), dimension(1:4) :: m_mul_m
  m_mul_m=(/ matrix1(1)*matrix2(1)+matrix1(2)*matrix2(3), &
             matrix1(1)*matrix2(2)+matrix1(2)*matrix2(4), &
             matrix1(3)*matrix2(1)+matrix1(4)*matrix2(3), &
             matrix1(3)*matrix2(2)+matrix1(4)*matrix2(4) /)
  end function m_mul_m

  function m_mul_v (matrix, vector)
  implicit none
  real(kind=8), dimension(1:4) :: matrix
  real(kind=8), dimension(1:2) :: vector
  real(kind=8), dimension(1:2) :: m_mul_v
  m_mul_v=(/ matrix(1)*vector(1)+matrix(2)*vector(2), &
             matrix(3)*vector(1)+matrix(4)*vector(2) /)
  end function m_mul_v

end module lfpm_functions


! -------------------------- subroutine ------------------------------
subroutine lfpm (h,nx,ny,p)
! This f90 program is modified from Hergarten's work and his c++ codes.
! Reference: Hergarten and Robl, GMD, 2022
! Paper DOI: https://doi.org/10.5194/gmd-15-2063-2022
! Codes availability: http://hergarten.at/openlem/lfpm.php

use lfpm_functions

implicit none

!!! Note: Take care of the axis directions
integer, intent(in) :: nx,ny
integer :: i,j
real(kind=8), intent(in) :: h(nx,ny)     ! elevation 
real(kind=8) :: qv(nx,ny)    ! vapor flux
real(kind=8) :: qc(nx,ny)    ! cloud water flux
real(kind=8) :: ptot(nx,ny)  ! total precipitation
real(kind=8), intent(out) :: p(nx,ny)     ! effective precipitation (total minus evaporation)
real(kind=8) :: psi(nx)      ! effective precipitation rate in each loop

real(kind=8) :: lc        ! length scale of condensation
real(kind=8) :: lf        ! length scale of fallout
real(kind=8) :: ll        ! length scale of long-range transport
real(kind=8) :: ld        ! length scale of dispersion
real(kind=8) :: refheight ! reference elevation
real(kind=8) :: evap      ! evaporation fraction at sea level
real(kind=8) :: qin       ! influx from the boundary (such as sea)

real(kind=8) :: beta      ! beta0 at sea level

!!! Note: Take care of the axis directions
real(kind=8) :: diag(ny,4)     ! main diagonal of the matrix
real(kind=8) :: upper(ny-1,4)  ! directly right to the main diagonal
real(kind=8) :: lower(ny-1,4)  ! directly below the main diagonal
real(kind=8) :: rhs(ny,2)      ! right-hand side of the linear equation system

real(kind=8) :: f
real(kind=8) :: fm(4)

!do j=ny,1,-1
!  print*, j, ( h(i,j), i=1,nx )
!enddo
!print*, "----------------- Check Topography ------------------"

! Set parameters (unit: grid, 1 grid = 625 m)
! lc=lf=25 km
lc=40.0d0
lf=40.0d0
! play with ll and qin
! ll=350 km
ll=600.0d0
! ld=25 km
ld=40.0d0
! reference height 1 km
refheight=1.6d0
! evaporation
evap=0.5d0
! flux input (qin=0.8 km2/yr)
qin=2.0d0 ! unit: grid^2/year

beta=(1-lc/ll)*(ll/lf-1)

! Set the boundary condition for influx (from sea)
do j=1,ny !!! Note: Parallel to "coastline", normal to advection direction (i=1)
  f=lf/(ll-lf)*exp(h(nx,j)/refheight)
  rhs(j,1)=qin/(1+f)
  rhs(j,2)=qin*f/(1+f)
  qv(nx,j)=rhs(j,1)
  qc(nx,j)=rhs(j,2)
  ptot(nx,j)=rhs(j,2)/lf
  psi(j)=1-evap*exp(-h(nx,j)/refheight)
  p(nx,j)=psi(j)*rhs(j,2)/lf
enddo

! We have stored the real rhs(1,1:2) and rhs(nx,1:2) for j=1.
! But beacause of homogeneous Neumann boundary, here we set
! them as (0,0) only used in the next loop of j=2.
rhs(1,1:2)=(/0.0d0,0.0d0/)
rhs(ny,1:2)=(/0.0d0,0.0d0/)

! Loop in the direction of advection (i), each uses the fluxes of the previous.
do i=nx-1,1,-1

  ! Set up a sparse matrix consisting of 2x2 blocks.
  do j=1,ny
    f=exp(-h(i,j)/refheight)
    psi(j)=1-f*evap
    f=f*beta/lc
    if (j.gt.1.and.j.lt.ny) then
      diag(j,1:4)=(/2*ld+1+1/lc,-f-(1-psi(j))/lf,-1/lc,2*ld+1+f+1/lf/) ! (m11,m12,m21,m22)
    elseif (j.eq.1) then
      diag(j,1:4)=(/1.0d0,0.0d0,0.0d0,1.0d0/) ! E(1,0,0,1) from homogeneous Neumann boundary
    elseif (j.eq.ny) then
      diag(j,1:4)=(/1.0d0,0.0d0,0.0d0,1.0d0/) ! E(1,0,0,1) from homogeneous Neumann boundary
    endif
!    print*, j,i,diag(i,:)
  enddo
  do j=1,ny-1
    if (j.gt.1) then
      upper(j,1:4)=(/-ld,0.0d0,0.0d0,-ld/) ! (m11,m12,m21,m22)
    else
      upper(j,1:4)=(/-1.0d0,0.0d0,0.0d0,-1.0d0/) ! -E(1,0,0,1) from homogeneous Neumann boundary
    endif
!    print*, j,i,upper(i,:)
  enddo
  do j=1,ny-1
    if (j.lt.ny-1) then
      lower(j,1:4)=(/-ld,0.0d0,0.0d0,-ld/) ! (m11,m12,m21,m22)
    else
      lower(j,1:4)=(/-1.0d0,0.0d0,0.0d0,-1.0d0/) ! -E(1,0,0,1) from homogeneous Neumann boundary
    endif
!    print*, j,i,lower(i,:)
  enddo

  ! Gauss algorithm for converting the matrix to an upper tridiagonal matrix
  do j=1,ny-1
    fm=m_mul_m(m_inv(diag(j,1:4)), lower(j,1:4))
    diag(j+1,1:4)=diag(j+1,1:4)-m_mul_m(fm, upper(j,1:4))
    rhs(j+1,1:2)=rhs(j+1,1:2)-m_mul_v(fm, rhs(j,1:2))
  enddo

  ! Solve the tridiagonal system from the bottom
  rhs(ny,1:2)=m_mul_v(m_inv(diag(ny,1:4)), rhs(ny,1:2))
  rhs(ny-1,1:2)=rhs(ny-1,1:2)-m_mul_v(upper(ny-1,1:4), rhs(ny,1:2))
  rhs(ny-1,1:2)=m_mul_v(m_inv(diag(ny-1,1:4)), rhs(ny-1,1:2))
  do j=ny-2,1,-1
    rhs(j,1:2)=rhs(j,1:2)-m_mul_v(upper(j,1:4), rhs(j+1,1:2))
    rhs(j,1:2)=m_mul_v(m_inv(diag(j,1:4)), rhs(j,1:2))
  enddo
  do j=1,ny
    qv(i,j)=rhs(j,1)
    qc(i,j)=rhs(j,2)
    ptot(i,j)=rhs(j,2)/lf
    p(i,j)=psi(j)*rhs(j,2)/lf
  enddo

  ! We have stored the real rhs(1,1:2) and rhs(nx,1:2) for j=j.
  ! But beacause of homogeneous Neumann boundary, here we set
  ! them as (0,0) only used in the next loop of j=j+1.
  rhs(1,1:2)=(/0.0d0,0.0d0/)
  rhs(ny,1:2)=(/0.0d0,0.0d0/)

enddo

!! Output
!print*, "------ qv:"
!do j=ny,1,-1
!  print*, j, ( qv(i,j), i=1,nx )
!enddo
!print*, "------ qc:"
!do j=ny,1,-1
!  print*, j, ( qc(i,j), i=1,nx )
!enddo
!print*, "------ ptot:"
!do j=ny,1,-1
!  print*, j, ( ptot(i,j), i=1,nx )
!enddo
!print*, "------ p:"
!do j=ny,1,-1
!  print*, j, ( p(i,j), i=1,nx )
!enddo

end subroutine lfpm
