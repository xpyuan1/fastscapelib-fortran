!--------------------------------------------------------------------------------------------

subroutine StreamPowerLaw ()

  ! subroutine to solve the stream power law equation following the FastScape method described
  ! in Braun and Willett, Geomorphology, 2015 and Yuan et al., JGR, 2019

  use FastScapeContext

  implicit none

  integer :: ij,ijk,ijr,k,ijr1
  double precision :: dx,dy,fact,tol,err
  double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
  double precision, dimension(:), allocatable :: ht,kfint,dh,hp
  double precision, dimension(:), allocatable :: elev
  double precision, dimension(:), allocatable :: water,lake_water_volume,lake_sediment
  integer, dimension(:), allocatable :: lake_sill

! define the parameters for calculating correlation coefficients (Xiaoping Yuan, Sep, 2023)
  double precision, dimension(:), allocatable :: crx,cry1,cry2,topomean
  double precision, dimension(:,:), allocatable :: crtopo,crprecip,cruplift,crerate,cretot
  double precision mean_x,mean_y1,mean_y2,cov,var_x,var_y1,var_y2,rpe,rue,width
  integer size_x,size_y1,size_y2,np,i,j
  allocate (topomean(nx))
  allocate (crtopo(nx,ny),crprecip(nx,ny),cruplift(nx,ny),crerate(nx,ny),cretot(nx,ny))

  allocate (ht(nn),kfint(nn),dh(nn),hp(nn))
  allocate (elev(nn))
  allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),lake_sill(nn))

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! set g, dimensionless parameter for sediment transport and deposition
  ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
  if (g1.ge.0.d0) then
    g=g1
    if (g2.gt.0.d0) where ((h-b).gt.1.d0) g=g2
  endif

  ! set kf / kfsed
  kfint=kf
  if (kfsed.gt.0.d0) where ((h-b).gt.1.d0) kfint=kfsed

  if (count(mstack==0).ne.0) print*,'incomplete stack',count(mstack==0),nn

  ! modified by Jean Braun (20/11/2022) to allow for relative versus
  ! absolute tolerance
  tol = tol_rel*maxval(abs(h)) + tol_abs
  err=2.d0*tol

  ! store the elevation at t
  ht=h

  ! Gauss-Seidel iteration
  nGSStreamPowerLaw=0

  lake_sediment=0.d0
  lake_sill=0.d0
  dh=0.d0
  hp=h
  
  do while (err.gt.tol.and.nGSStreamPowerLaw.lt.nGSStreamPowerLawMax-1)
    nGSStreamPowerLaw=nGSStreamPowerLaw+1

    where (bounds_bc)
      elev=ht
    elsewhere
!      elev=ht+(dh-(ht-hp))*g*dx*dy/a
      elev=ht+(dh-(ht-hp))*g*dx*dy/agrid
    endwhere

      ! apply modified stream power law using lake surface (hwater)

      if (abs(n-1.d0).lt.tiny(n)) then

        do ij=nn,1,-1
          ijk=mstack(ij)
          ijr1=rec(ijk)
          if (ijr1.eq.ijk) then
            water(ijk)=ht(ijk)
            lake_sill(ijk)=ijk
            lake_water_volume(ijk)=0.d0
          else
            w_rcv=water(ijr1)
            if (elev(ijk).gt.w_rcv) then
              if (mnrec(ijk).gt.0) then
                if (h(ijk).ge.sealevel.or..not.runMarine) then
                  f = elev(ijk)
                  df = 1.d0
                  do k=1,mnrec(ijk)
                    if (ht(ijk).ge.ht(mrec(k,ijk))) then
                      fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)
                      f = f + fact*h(mrec(k,ijk))
                      df = df + fact
                    endif
                  enddo
                  h(ijk)=f/df
                endif
              endif
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
              if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
            else
              h(ijk)=elev(ijk)
              lake_sill(ijk)=lake_sill(ijr1)
              if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
              lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
            endif
            water(ijk)=max(w_rcv,h(ijk))
          endif
        enddo

      else

        do ij=nn,1,-1
          ijk=mstack(ij)
          ijr1=rec(ijk)
          if (ijr1.eq.ijk) then
            water(ijk)=ht(ijk)
            lake_sill(ijk)=ijk
            lake_water_volume(ijk)=0.d0
          else
            w_rcv=water(ijr1)
            if (elev(ijk).gt.w_rcv) then
              if (mnrec(ijk).gt.0) then
                if (ht(ijk).ge.sealevel.or..not.runMarine) then
                  omega=0.875d0/n
                  tolp=1.d-3
                  errp=2.d0*tolp
                  h0=elev(ijk)
                  do while (errp.gt.tolp)
                    f=h(ijk)-h0
                    df=1.d0
                    do k=1,mnrec(ijk)
                      if (ht(ijk).gt.ht(mrec(k,ijk))) then
                        fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)**n
                        f=f+fact*max(0.d0,h(ijk)-h(mrec(k,ijk)))**n
                        df=df+fact*n*max(0.d0,h(ijk)-h(mrec(k,ijk)))**(n-1.d0)
                      endif
                    enddo
                    hn=h(ijk)-f/df
                    errp=abs(hn-h(ijk))
                    h(ijk)=h(ijk)*(1.d0-omega)+hn*omega
                  enddo
                endif
              endif
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
              if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
            else
              h(ijk)=elev(ijk)
              lake_sill(ijk)=lake_sill(ijr1)
              if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
              lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
            endif
            water(ijk)=max(w_rcv,h(ijk))
          endif
        enddo

      endif

      err=sqrt(sum((h-hp)**2)/nn)
    ! Jean Braun modification 18/11/2022: moved the computation of redistribution of sediment in lakes
    ! following Sebastian Wolf's suggestion; this ensures mass conservation in multi-minima cases
    ! guess/update the elevation at t+Dt (k)
    hp=h

    ! calculate erosion/deposition at each node
    dh=ht-hp

    ! sum the erosion in stack order
    do ij=1,nn
      ijk=mstack(ij)
      ijr1=rec(ijk)
      if (ijr1.ne.ijk) then
        dh(ijk)=dh(ijk)-(ht(ijk)-hp(ijk))
        if (lake_sill(ijk).eq.ijk) then
          if (dh(ijk).le.0.d0) then
            lake_sediment(ijk)=0.d0
          else
            lake_sediment(ijk)=dh(ijk)
          endif
        endif
        dh(ijk)=dh(ijk)+(ht(ijk)-hp(ijk))
        do k=1,mnrec(ijk)
          ijr=mrec(k,ijk)
          dh(ijr)=dh(ijr)+dh(ijk)*mwrec(k,ijk)
        enddo
      else
        lake_sediment(ijk)=dh(ijk)
      endif
    enddo

      if (maxval(g).lt.tiny(g)) err=0.d0

    enddo

    b=min(h,b)

    do ij=1,nn
      if (lake_sill(ij).ne.0) then
        if (lake_water_volume(lake_sill(ij)).gt.0.d0) h(ij)=h(ij) &
        +max(0.d0,min(lake_sediment(lake_sill(ij)),lake_water_volume(lake_sill(ij))))/ &
        lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))
      endif
    enddo

    ! stores total erosion, erosion rate and flux for output
    etot=etot+ht-h
    erate=(ht-h)/dt
    Sedflux=ht-h
    !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0

    open (30,file='erate_M.txt',status='unknown',position='append')
    write (30,*) step,minval(erate),sum(erate)/real(nn),maxval(erate)
    close (30)
    open (40,file='etot_M.txt',status='unknown',position='append')
    write (40,*) step,minval(etot),sum(etot)/real(nn),maxval(etot)
    close (40)
    open (50,file='Sedflux_M.txt',status='unknown',position='append')
    write (50,*) step,sum(Sedflux)
    close (50)

deallocate (topomean)
deallocate (crtopo,crprecip,cruplift,crerate,cretot)

    deallocate (ht,kfint,dh,hp,elev,water,lake_water_volume,lake_sediment,lake_sill)

    return

  end subroutine StreamPowerLaw

  !--------------------------------------------------------------------------------------------

  subroutine StreamPowerLawSingleFlowDirection ()

    ! subroutine to solve the stream power law equation following the FastScape method described
    ! in Braun and Willett, Geomorphology, 2015

    use FastScapeContext

    implicit none

    integer :: ij,ijk,ijr
    double precision :: dx,dy,fact,tol,err
    double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
    double precision, dimension(:), allocatable :: ht,kfint,dh,hp
    double precision, dimension(:), allocatable :: elev
    double precision, dimension(:), allocatable :: water,lake_water_volume,lake_sediment
    integer, dimension(:), allocatable :: lake_sill

! define the parameters for calculating correlation coefficients (Xiaoping Yuan, Sep, 2023)
  double precision, dimension(:), allocatable :: crx,cry1,cry2,topomean,topomin,topomax
  double precision, dimension(:,:), allocatable :: crtopo,crprecip,cruplift,crerate,cretot
  double precision mean_x,mean_y1,mean_y2,cov,var_x,var_y1,var_y2,rpe,rue,width,topominmax
  integer size_x,size_y1,size_y2,np,i,j,nxmid
  allocate (topomean(nx),topomin(nx),topomax(nx))
  allocate (crtopo(nx,ny),crprecip(nx,ny),cruplift(nx,ny),crerate(nx,ny),cretot(nx,ny))

    allocate (ht(nn),kfint(nn),dh(nn),hp(nn))
    allocate (elev(nn))
    allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),lake_sill(nn))

    dx=xl/(nx-1)
    dy=yl/(ny-1)

    ! set g, dimensionless parameter for sediment transport and deposition
    ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
    if (g1.ge.0.d0) then
      g=g1
      if (g2.gt.0.d0) where ((h-b).gt.1.d0) g=g2
    endif

    ! set kf / kfsed
    kfint=kf
    if (kfsed.gt.0.d0) where ((h-b).gt.1.d0) kfint=kfsed

    ! modified by Jean Braun (20/11/2022) to allow for relative versus
    ! absolute tolerance
    tol = tol_rel*maxval(abs(h)) + tol_abs 
    err=2.d0*tol

    ! store the elevation at t
    ht=h

    ! Gauss-Seidel iteration
    nGSStreamPowerLaw=0

    lake_sediment=0.d0
    dh=0.d0
    hp=h

    do while (err.gt.tol.and.nGSStreamPowerLaw.lt.nGSStreamPowerLawMax-1)
      nGSStreamPowerLaw=nGSStreamPowerLaw+1
      where (bounds_bc)
        elev=ht
      elsewhere
!        elev=ht+(dh-(ht-hp))*g*dx*dy/a
        elev=ht+(dh-(ht-hp))*g*dx*dy/agrid
      endwhere

        ! apply modified stream power law using lake surface (hwater)

        if (abs(n-1.d0).lt.tiny(n)) then

          do ij=1,nn
            ijk=stack(ij)
            ijr=rec(ijk)
            if (ijr.eq.ijk) then
              water(ijk)=ht(ijk)
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
            else
              w_rcv=water(ijr)
              if (elev(ijk).gt.w_rcv) then
                if (h(ijk).ge.sealevel.or..not.runMarine) then
                f = elev(ijk)
                df = 1.d0
                ! todo: check if we don't need those checks for single flow
!                if (ht(ijk).ge.ht(ijr)) then
                  fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)
                  f = f + fact*h(ijr)
                  df = df + fact
!                endif
                h(ijk)=f/df
!                h(ijk)=min(f/df,minval(h(don(1:ndon(ijk),ijk))))
                  endif
                lake_sill(ijk)=ijk
                lake_water_volume(ijk)=0.d0
                if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
              else
                h(ijk)=elev(ijk)
                lake_sill(ijk)=lake_sill(ijr)
                if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
              endif
              water(ijk)=max(w_rcv,h(ijk))
            endif
          enddo

        else

          do ij=1,nn
            ijk=stack(ij)
            ijr=rec(ijk)
            if (ijr.eq.ijk) then
              water(ijk)=ht(ijk)
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
            else
              w_rcv=water(ijr)
              if (elev(ijk).gt.w_rcv) then
                if (ht(ijk).ge.sealevel.or..not.runMarine) then
                omega=0.875d0/n
                tolp=1.d-3
                errp=2.d0*tolp
                h0=elev(ijk)
                do while (errp.gt.tolp)
                  f=h(ijk)-h0
                  df=1.d0
                  if (ht(ijk).gt.ht(ijr)) then
                    fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)**n
                    f=f+fact*max(0.d0,h(ijk)-h(ijr))**n
                    df=df+fact*n*max(0.d0,h(ijk)-h(ijr))**(n-1.d0)
                  endif
                  hn=h(ijk)-f/df
                  errp=abs(hn-h(ijk))
                  h(ijk)=h(ijk)*(1.d0-omega)+hn*omega
                enddo
                endif
                lake_sill(ijk)=ijk
                lake_water_volume(ijk)=0.d0
                if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
              else
                h(ijk)=elev(ijk)
                lake_sill(ijk)=lake_sill(ijr)
                if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
              endif
              water(ijk)=max(w_rcv,h(ijk))
            endif
          enddo

        endif

        err=sqrt(sum((h-hp)**2)/nn)

      ! Jean Braun modification 18/11/2022: moved the computation of redistribution of sediment in lakes
      ! following Sebastian Wolf's suggestion; this ensures mass conservation in multi-minima cases
      ! guess/update the elevation at t+Dt (k)
      hp=h

       ! calculate erosion/deposition at each node
      dh=ht-hp

      ! sum the erosion in stack order
      do ij=nn,1,-1
        ijk=stack(ij)
        ijr=rec(ijk)
        if (ijr.ne.ijk) then
          dh(ijk)=dh(ijk)-(ht(ijk)-hp(ijk))
          if (lake_sill(ijk).eq.ijk) then
            if (dh(ijk).le.0.d0) then
              lake_sediment(ijk)=0.d0
            else
              lake_sediment(ijk)=min(dh(ijk),lake_water_volume(ijk))
              dh(ijk) = dh(ijk)-lake_water_volume(ijk) !remove the sediment that is going to be deposited in the lake from the dh stack
              if (dh(ijk)<0.d0) dh(ijk)=0.d0
            endif
          endif
          dh(ijk)=dh(ijk)+(ht(ijk)-hp(ijk))
          dh(ijr)=dh(ijr)+dh(ijk)
        else
          lake_sediment(ijk)=dh(ijk)
        endif
      enddo

        if (maxval(g).lt.tiny(g)) err=0.d0
      enddo

      b=min(h,b)

      do ij=1,nn
        if (lake_sill(ij).ne.0) then
          if (lake_water_volume(lake_sill(ij)).gt.0.d0) h(ij)=h(ij) &
          +max(0.d0,min(lake_sediment(lake_sill(ij)),lake_water_volume(lake_sill(ij))))/ &
          lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))
        endif
      enddo

      ! stores total erosion, erosion rate and flux for output
      etot=etot+ht-h
      where(etot.lt.0.d0) etot=0.d0
      erate=(ht-h)/dt
      Sedflux=ht-h
      !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0

      open (30,file='erate.txt',status='unknown',position='append')
      write (30,*) step,minval(erate),sum(erate)/real(nn),maxval(erate)
      close (30)
      open (40,file='etot.txt',status='unknown',position='append')
      write (40,*) step,minval(etot),sum(etot)/real(nn),maxval(etot)
      close (40)
      open (50,file='Sedflux.txt',status='unknown',position='append')
      write (50,*) step,sum(Sedflux)
      close (50)

! To calculate the correlation coefficient that measures
! linear correlation between two sets of data
! https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
! calculate correlation coefficients (Xiaoping Yuan, Sep, 2023)
do j=1,ny
do i=1,nx
ij=i+(j-1)*nx
crtopo(i,j)=h(ij)
crprecip(i,j)=precip(ij)
cruplift(i,j)=u(ij)
crerate(i,j)=erate(ij)
cretot(i,j)=etot(ij)
enddo
enddo

! obtain min, mean, and max of topo along x-axial
do i=1,nx
topomin(i)=minval(crtopo(i,:))
topomax(i)=maxval(crtopo(i,:))
topomean(i)=sum(crtopo(i,:))/real(ny)
enddo

! find the location of divide for the pro-side mountain
topominmax=maxval(topomin)
do i=1,nx
  if (abs(topomin(i)-topominmax).lt.1.d-3) nxmid=i
enddo

! --------------- elevation > 200 m -------------------
! count only erosion, uplift, and precipitation for elevation > 200 m
np=0
!i=320 refers 200 km from left boundary
do i=320,nx-320
if (topomean(i).gt.200.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nx-320
if (topomean(i).gt.200.d0) then
  np=np+1
  crx(np)=sum(crerate(i,:))/real(ny)
  cry1(np)=sum(crprecip(i,:))/real(ny)
  cry2(np)=sum(cruplift(i,:))/real(ny)
endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (60,file='1_correlation_200.txt',status='unknown',position='append')
write (60,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (60)

deallocate (crx,cry1,cry2)


! --------------- elevation > 500 m -------------------
! count only erosion, uplift, and precipitation for elevation > 500 m
np=0
do i=320,nx-320
if (topomean(i).gt.500.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nx-320
if (topomean(i).gt.500.d0) then
  np=np+1
  crx(np)=sum(crerate(i,:))/real(ny)
  cry1(np)=sum(crprecip(i,:))/real(ny)
  cry2(np)=sum(cruplift(i,:))/real(ny)
endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (70,file='1_correlation_500.txt',status='unknown',position='append')
write (70,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (70)

deallocate (crx,cry1,cry2)


! --------------- elevation > 1000 m -------------------
! count only erosion, uplift, and precipitation for elevation > 1000 m
np=0
do i=320,nx-320
if (topomean(i).gt.1000.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nx-320
if (topomean(i).gt.1000.d0) then
  np=np+1
  crx(np)=sum(crerate(i,:))/real(ny)
  cry1(np)=sum(crprecip(i,:))/real(ny)
  cry2(np)=sum(cruplift(i,:))/real(ny)
endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (80,file='1_correlation_1000.txt',status='unknown',position='append')
write (80,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (80)

deallocate (crx,cry1,cry2)


! --------------- elevation > 200 m on the pro-side mountain ------------------
! count only erosion, uplift, and precipitation for elevation > 200 m
np=0
do i=320,nxmid
  if (topomax(i).gt.200.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nxmid
  if (topomax(i).gt.200.d0) then
    np=np+1
    crx(np)=sum(crerate(i,:))/real(ny)
    cry1(np)=sum(crprecip(i,:))/real(ny)
    cry2(np)=sum(cruplift(i,:))/real(ny)
  endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (90,file='2_correlation_200.txt',status='unknown',position='append')
write (90,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (90)

deallocate (crx,cry1,cry2)


! --------------- elevation > 500 m on the pro-side mountain ------------------
! count only erosion, uplift, and precipitation for elevation > 500 m
np=0
do i=320,nxmid
  if (topomax(i).gt.500.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nxmid
  if (topomax(i).gt.500.d0) then
    np=np+1
    crx(np)=sum(crerate(i,:))/real(ny)
    cry1(np)=sum(crprecip(i,:))/real(ny)
    cry2(np)=sum(cruplift(i,:))/real(ny)
  endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (100,file='2_correlation_500.txt',status='unknown',position='append')
write (100,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (100)

deallocate (crx,cry1,cry2)


! --------------- elevation > 1000 m on the pro-side mountain ------------------
! count only erosion, uplift, and precipitation for elevation > 1000 m
np=0
do i=320,nxmid
  if (topomax(i).gt.1000.d0) np=np+1
enddo

allocate (crx(np),cry1(np),cry2(np))

np=0
do i=320,nxmid
  if (topomax(i).gt.1000.d0) then
    np=np+1
    crx(np)=sum(crerate(i,:))/real(ny)
    cry1(np)=sum(crprecip(i,:))/real(ny)
    cry2(np)=sum(cruplift(i,:))/real(ny)
  endif
enddo

size_x = size(crx)
size_y1 = size(cry1)
size_y2 = size(cry2)

mean_x = sum(crx) / size_x
var_x = sum((crx(1:size_x) - mean_x)*(crx(1:size_x) - mean_x))

mean_y1 = sum(cry1) / size_y1
cov = sum((crx(1:size_x) - mean_x)*(cry1(1:size_x) - mean_y1))
var_y1 = sum((cry1(1:size_x) - mean_y1)*(cry1(1:size_x) - mean_y1))
rpe = (cov / sqrt(var_x)) / sqrt(var_y1)

mean_y2 = sum(cry2) / size_y2
cov = sum((crx(1:size_x) - mean_x)*(cry2(1:size_x) - mean_y2))
var_y2 = sum((cry2(1:size_x) - mean_y2)*(cry2(1:size_x) - mean_y2))
rue = (cov / sqrt(var_x)) / sqrt(var_y2)

open (110,file='2_correlation_1000.txt',status='unknown',position='append')
write (110,*) step,dt,rpe,rue,np,real(np)*dx/1.d3,maxval(h),maxval(topomean)
close (110)

deallocate (crx,cry1,cry2)


deallocate (topomean,topomin,topomax)
deallocate (crtopo,crprecip,cruplift,crerate,cretot)

      deallocate (ht,kfint,dh,hp,elev,water,lake_water_volume,lake_sediment,lake_sill)

      return

    end subroutine StreamPowerLawSingleFlowDirection

    !--------------------------------------------------------------------------------------------
