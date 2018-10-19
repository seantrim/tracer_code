subroutine viscosity_contrast(level,diff,T0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,diff,a,T0(-span:nx+span,-span:nz+span)
integer*4 :: i,k,level,nxx,nzz
if (iterate.eq.0) then
 nxx=nx
 nzz=nz
elseif (iterate.eq.1) then
 nxx=nx_grid(level)
 nzz=nz_grid(level)
end if

diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a) REDUCTION(max:diff)
do i=0,nxx-1
 do k=1,nzz-1
  a1=max(T0(i,k-1),T0(i+1,k-1),T0(i+1,k),T0(i,k+1),T0(i+1,k+1))
  a2=min(T0(i,k-1),T0(i+1,k-1),T0(i+1,k),T0(i,k+1),T0(i+1,k+1))
  a3=max(a1,T0(i,k))/min(a1,T0(i,k))
  a4=max(a2,T0(i,k))/min(a2,T0(i,k))
  a=max(a3,a4)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
end


subroutine compute_smoothing_parameters
 use basics
 use arrays
implicit none
integer*4 :: i
real*8 :: A(1:ntype)
!!!!!!!!!!!!!!!!!!!!!!!!!!Temperature Smoothing
if ((ivis.eq.0).or.(RaT.eq.0.d0)) then
 delta_T=2.d0 !!no smoothing required
 delta_T_init=2.d0
elseif (ivis.ge.1) then
 if (visT.ne.1.d0) then
  delta_T=dlog(delta_visT)/dlog(visT)
  delta_T_init=dlog(delta_visT_init)/dlog(visT)
 else
  delta_T=2.d0
  delta_T_init=2.d0
 end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!End Temperature Smoothing

!!!!!!!!!!!!!!!!!!!!!!!!!!Composition Smoothing
if (comp.eq.1) then
 allocate(delta_C(1:ntype))
 delta_C=2.d0 !!default -- no smoothing
 do i=1,ntype
  if (conduct_factor(i).ne.1.d0) then
   delta_C(i)=dabs(dlog(delta_visC(i))/dlog(1.d0/conduct_factor(i))) 
  end if
 end do
 write(*,*) "conductivity delta_C=",delta_C
 if (ivis.eq.2) then
  do i=1,ntype
   if (visC(i).ne.1.d0) then
    delta_C(i)=dabs(dlog(delta_visC(i))/dlog(visC(i))) !!!!!vectorize
   end if
  end do
 end if
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!End Composition Smoothing
write(*,*) "T Smoothing Parameter:", delta_T
if (comp.eq.1) write(*,*) "C Smoothing Parameters:",delta_C
end

subroutine smoothness(diff,T0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,a5,diff,a,T0(-span:nx+span,-span:nz+span)
integer*4 :: i,k
diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a5,a) REDUCTION(max:diff)
do i=0,nx-1
 do k=1,nz-1
  a1=dabs(T0(i,k-1)-T0(i,k))
  a2=dabs(T0(i+1,k-1)-T0(i,k))
  a3=dabs(T0(i+1,k)-T0(i,k))
  a4=dabs(T0(i,k+1)-T0(i,k))
  a5=dabs(T0(i+1,k+1)-T0(i,k))
  a=max(a1,a2,a3,a4,a5)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
end

subroutine smoothness_residual(level,diff,r0)
!$ use OMP_LIB
 use basics
implicit none
real*8 :: a1,a2,a3,a4,a5,diff,a,r0(1:nx-1,1:nz-1)
integer*4 :: i,k,level
diff=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,a1,a2,a3,a4,a5,a) REDUCTION(max:diff)
do i=1,nx_grid(level)-2
 do k=2,nz_grid(level)-2
  a1=dabs(r0(i,k-1)-r0(i,k))
  a2=dabs(r0(i+1,k-1)-r0(i,k))
  a3=dabs(r0(i+1,k)-r0(i,k))
  a4=dabs(r0(i,k+1)-r0(i,k))
  a5=dabs(r0(i+1,k+1)-r0(i,k))
  a=max(a1,a2,a3,a4,a5)
  if (a.gt.diff) diff=a
 end do
end do
!$OMP END PARALLEL DO
diff=diff/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8) !!scale since r0 array has not been normalized
end


real*8 function smoother_space(i,k,T0) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz
real*8 :: T0(-span:nx+span,-span:nz+span)
Tr=dot_product(D1(-span1:span1),T0(i-span1:i+span1,k))/dr
Ts=dot_product(D1(-span1:span1),T0(i,k-span1:k+span1))/ds
Trr=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2
Tss=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2

Tx=Tr/xr(i)
Tz=Ts/zs(k)
Txx=(Trr-Tx*xrr(i))/xr(i)**2.d0
Tzz=(Tss-Tz*zss(k))/zs(k)**2.d0
smoother_space=Txx+Tzz
end

real*8 function vis_space(i,k,T0,level) !!T0 represents the composition field
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz
real*8 :: T0(-span:nx+span,-span:nz+span)

i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
k1=k*2**(level-1)

Tr=dot_product(D1(-span1:span1),T0(i-span1:i+span1,k))/dr_grid(level)
Ts=dot_product(D1(-span1:span1),T0(i,k-span1:k+span1))/ds_grid(level)
Trr=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2_grid(level)
Tss=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2_grid(level)

Tx=Tr/xr(i1)
Tz=Ts/zs(k1)
Txx=(Trr-Tx*xrr(i1))/xr(i1)**2.d0
Tzz=(Tss-Tz*zss(k1))/zs(k1)**2.d0
vis_space=Txx+Tzz
end

subroutine smoother_vis(level,interpolate) !!smoother for the viscosity field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,n,level,interpolate,kmin,kmax
integer*4 :: nn,ntemp,Niter,trigger,option
real*8 :: vis_space,temp,diff,c,ratio_min
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,la,lb

if ((level.eq.1).and.(interpolate.eq.0)) then !!filter out unresolved frequencies of the physical viscosity
 option=1
 Niter=6
 lb=-5.d0
 la=-2.d0
else          !!!for interpolation of viscosity to coarser grids during multigrid: only retain the lower half of the frequency spectrum
 option=0
 Niter=vis_iter
 lb=-2.d0
 la=-1.d0
end if

if ((interpolate.eq.0).or.(level.eq.1)) then
 Nvis_store=0
 nn_store=1
end if

if (RaT.eq.0.d0) then !!isothermal convection: viscosity not fixed at horizontal boundaries
 kmin=0
 kmax=nz_grid(level)
else                  !!thermal/thermochemical convection: viscosity is fixed by T values (assumed fixed -- insulating bottom not included yet)
 kmin=1
 kmax=nz_grid(level)-1
end if

n=0
trigger=0
do
do nn=nn_store,Niter !!for Richardson's method

  call viscosity_gradients(level)
! call viscosity_contrast(level,diff,vis_grid(level,-span:nx+span,-span:nz+span))
  call test_eigenvalues(level)
  ratio_min=minval(ratio(level,1:nx_grid(level)-1,1:nz_grid(level)-1))

 if (option.eq.1) then
  if (((Nvis_store.ge.1).and.(ratio_min.gt.eig_ratio_min)).or.(ivis.eq.0).or.(n.gt.20)) then
   eig_fine=ratio_min
   n_eta=n
!   write(*,*) "smoother_vis:",time,n,nn,ratio_min
   trigger=1
   exit
  end if
 elseif (option.eq.0) then
  if (interpolate.eq.0) then !!!prep for smoother iterations
   if ((ratio_min.gt.eig_ratio).or.(ivis.eq.0)) then
!    write(*,*) "smoother_vis:",level,n,nn,ratio_min,interpolate
    trigger=1
    nn_store=nn
    exit
   end if
  elseif (interpolate.eq.1) then !!!prep for restriction of viscosity to coarser grid
    if ((Nvis_store.ge.1).or.(ivis.eq.0)) then
!     write(*,*) "smoother_vis:",level,n,nn,ratio_min,interpolate
     trigger=1
     exit
    end if
  end if
 end if

 ntemp=Niter-nn+1 !!!start with smaller time steps to kill off high frequencies first
 c=2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8)))

 T0=vis_grid(level,-span:nx+span,-span:nz+span) !!feed in sharp viscosity field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
 do i=0,nx_grid(level)
  do k=kmin,kmax
   temp=vis_space(i,k,T0,level)
   dts=dt_vis(level,i,k)*0.99d0
   vis_grid(level,i,k)=T0(i,k)+c*dts*temp
  end do
 end do
!$OMP END PARALLEL DO
 call viscosityBCs(level)
end do

 if (trigger.eq.1) exit
 Nvis_store=Nvis_store+1
 nn_store=1
 if (option.eq.1) then !!!Once high unresolved frequencies are damped, start damping resolved frequencies to increase the eigenvalue ratio
  nn_store=1
 Niter=6
 lb=-2.d0
 la=-1.d0
 end if
n=n+1
end do
end

subroutine smoother_res(level,res)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: res,vis_space
real*8 :: temp(0:nx,0:nz)
real*8 :: T0(-span:nx+span,-span:nz+span)
 temp=0.d0
 T0=vis_grid(level,-span:nx+span,-span:nz+span) !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k)
 do i=0,nx_grid(level)
  do k=1,nz_grid(level)-1
   temp(i,k)=vis_space(i,k,T0,level)
  end do
 end do
!$OMP END PARALLEL DO
res=sum(dabs(temp(0:nx_grid(level),1:nz_grid(level)-1)))/real(nx_grid(level)+1,8)/real(nz_grid(level)-1,8)
end

subroutine smoother_time !!smoother for the C field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nn,tt,ntemp,Niter,trigger
real*8 :: smoother_space,temp,diff,c
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,la,lb
 Niter=10  !!!first pass: filter out unresolved frequencies and the upper half of the resolved spectrum to smooth C while supressing spurious oscillations
 lb=-5.d0
 la=-1.d0

 Cvis=Cnew
 Cbuoy=Cnew
do tt=1,ntype
nn=0
trigger=0
do
do n=1,Niter
 call smoothness(diff,Cvis(tt,:,:))
 if (diff.lt.delta_C(tt)) then
  trigger=1
  exit
 end if

 if (nn.eq.0) then !!!filter out unresolved frequencies on first pass
  ntemp=Niter-n+1 !!!start with smaller time steps to kill off high frequencies first
  c=2.d0/(-lb-la+(lb-la)*dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8)))
 else
  c=1.0d0  !!continue with more efficient smoothing once highest frequencies are damped
 end if

T0=Cvis(tt,-span:nx+span,-span:nz+span) !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do i=0,nx
 do k=0,nz
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantC
  Cvis(tt,i,k)=T0(i,k)+c*dts*temp
  if (Cvis(tt,i,k).lt.0.d0) Cvis(tt,i,k)=0.d0
  if (Cvis(tt,i,k).gt.1.d0) Cvis(tt,i,k)=1.d0
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_comp(Cvis(tt,:,:))
end do
 if (trigger.eq.1) exit
 Niter=1
 nn=nn+1
end do
 n_comp(tt)=nn
! write(*,*) "C smoother:",nn,n,diff,tt
end do
end

subroutine smoother_time_T(Tin) !!smoother for the temperature field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nsmall,Nsmoother,ntemp
real*8 :: smoother_space,temp,diff
real*8 :: T0(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,Tin(-span:nx+span,-span:nz+span)
 Tvis=Tin
 Tbuoy=Tin
if (RaT.eq.0) return
n=0
do
 call smoothness(diff,Tvis)
 if (diff.lt.delta_T) then
  n_temp=n
!  write(*,*) "T smoother:",n,diff,tstep
  exit
 end if

T0=Tvis !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do i=0,nx
 do k=0,nz
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantT
  Tvis(i,k)=T0(i,k)+dts*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs(Tvis)
n=n+1
end do
end

subroutine smoother_initial_T !!smoother for the initial temperature field
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,nsmall,Nsmoother,ntemp
real*8 :: smoother_space,temp,diff
real*8 :: T0(-span:nx+span,-span:nz+span),Tdumb(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: time1,time2,Tin(-span:nx+span,-span:nz+span)
n=0
do
 call smoothness(diff,T)
 if (diff.lt.delta_T_init) then
  write(*,*) "initial T smoother:",n,diff,tstep
  exit
 end if

T0=T !!feed in sharp C field to be smoothed
!$OMP PARALLEL DO PRIVATE(i,k,temp,dts)
do i=0,nx
 do k=0,nz
  temp=smoother_space(i,k,T0)
  dts=dts_array(i,k)*courantT
  T(i,k)=T0(i,k)+dts*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs(T)
n=n+1
end do

if (tracer.eq.1) then
 call extendT(T) !!for interpolation of T to Ttr during equalizing
 call compute_derivatives(Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp),DERr_Textend)
 Tdumb=Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
!$OMP PARALLEL DO PRIVATE(ID)
 do ID=1,ntr !!!interpolate temperature values to repositioned tracers
  call interpolate(rtr(ID),str(ID),DERr_Textend,Tdumb,Ttr(ID))
 end do
!$OMP END PARALLEL DO
 call tracers_to_corners(rtr,str,Ttr,T)
end if
end


subroutine stability_RK4_smoother
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: eigA_smoother,eigD_smoother,dtA,dtD,eR,eI
Ts_max=Tr_base/ds     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2 !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr
Trr_max=-Trr_base/dr2

!$OMP PARALLEL DO PRIVATE(i,k,eR,eI)
do i=0,nx !!using local time steps for smoother to account for grid spacing
 do k=0,nz
  eR=eigD_smoother(i,k)
  eI=eigA_smoother(i,k)
  dts_array(i,k)=2.d0*eR/(eR**2.d0+eI**2.d0)
 end do
end do
!$OMP END PARALLEL DO
end

real*8 function eigA_smoother(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
!!!!!!mathematical advection from grid spacing
  eigA_smoother=dabs(xrr(i)/(xr(i)**3.d0))*Tr_max+dabs(zss(k)/(zs(k)**3.d0))*Ts_max
end


real*8 function eigD_smoother(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
  eigD_smoother=dabs(Trr_max/(xr(i)**2.d0))+dabs(Tss_max/(zs(k)**2.d0))
end


subroutine stability_vis_smoother(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eigA_grid,eigD_grid,eR,eI
Ts_max=Tr_base/ds_grid(level)     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr_grid(level)
Trr_max=-Trr_base/dr2_grid(level)

!$OMP PARALLEL DO PRIVATE(i,k,eR,eI)
do i=0,nx_grid(level) !!using local time steps for smoother to account for grid spacing
 do k=1,nz_grid(level)-1
  eR=eigD_grid(level,i,k)
  eI=eigA_grid(level,i,k)
  dt_vis(level,i,k)=2.d0*eR/(eR**2.d0+eI**2.d0)
 end do
end do
!$OMP END PARALLEL DO
end

real*8 function eigA_grid(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
!!!!!!mathematical advection from grid spacing
  eigA_grid=dabs(xrr(i1)/(xr(i1)**3.d0))*Tr_max+dabs(zss(k1)/(zs(k1)**3.d0))*Ts_max
end


real*8 function eigD_grid(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
  i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
  k1=k*2**(level-1)
  eigD_grid=dabs(Trr_max/(xr(i1)**2.d0))+dabs(Tss_max/(zs(k1)**2.d0))
end


