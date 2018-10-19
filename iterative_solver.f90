subroutine optimize_grid(nx,nz,grids,n_points_min)
implicit none
integer*4 :: n_in,grids,nx,nz,nxbig
integer*4 :: n,n_points,n_coarse,n_temp,n_coarse_temp,n_points_temp
integer*4 :: n_min,n_max,n_points_max,n_coarse_max,n_points_min,n_coarse_min
n_points_max=min(nx,nz)   !!!max # of grid levels on coarsest grid
n_coarse_min=1
n_coarse_max=grids-1   !!!max # of coarse grids
n_min=n_points_min*2**n_coarse_min
n_max=n_points_max*2**n_coarse_max 

if (nz.le.nx) then !!# of coarse grids constrained by smaller of nx and nz
 nxbig=1
 n_in=nz
else
 nxbig=0
 n_in=nx
end if
if (n_min.gt.n_in) then
 write(*,*) "Grid Cannot Be Optimized With Current Settings"
 write(*,*) "Current settings allow a minimum grid size of",n_min,"."
end if
if (n_max.lt.n_in) then
 write(*,*) "Grid Cannot Be Optimized With Current Settings"
 write(*,*) "Current settings allow a maximum grid size of",n_max,"."
end if
n_temp=n_max !!start with upper limit
do n_points=n_points_min,n_points_max!!!range for # of grid levels on coarsest grid
 do n_coarse=n_coarse_min,n_coarse_max!!!range for total # of coarse grids 
  n=n_points*2**n_coarse !!# of grid levels on fine grid
  if ((n.ge.n_in).and.(n.lt.n_temp)) then
   n_temp=n
   n_coarse_temp=n_coarse
   n_points_temp=n_points
  end if
 end do
end do
!write(*,*) n_temp,n_points_temp,n_coarse_temp
grids=n_coarse_temp+1
if (nxbig.eq.0) then
 nx=n_temp
else
 nz=n_temp
end if

!!!given the # of coarse grids, figure out the resolution for the larger of nx and nz
if (nz.le.nx) then !!# of coarse grids constrained by smaller of nx and nz
 n_in=nx
else
 n_in=nz
end if
n_points_temp=nint(real(n_in,8)/real(2**n_coarse_temp,8),4)
n_temp=n_points_temp*2**n_coarse_temp
!write(*,*) n_temp,n_points_temp,n_coarse_temp
if (nxbig.eq.0) then
 nz=n_temp
else
 nx=n_temp
end if
write(*,*) "# of multigrid levels=",grids
end

subroutine print_field(level,field,fname)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: field(1:nx-1,1:nz-1)
character*6 :: fname
open(unit=999,file=fname)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   write(999,*) xg(i),zg(k),field(i,k)
 end do
end do
 close(999)
end

subroutine multigrid
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: level,dir,multiply,INFO,q,p,ii
real*8 :: ROWCND,COLCND,AMAX,diff,res_store,time1,time2
!$ time1=omp_get_wtime()
error(1,:,:)=SF(:,:)
 call compute_RHS(1)
 call compute_reference_residual(1)

if (non_Newtonian.eq.0) then !!if Newtonian rheology, these quantities do not need to be computed as often
 do level=1,grids
  call compute_viscosity(level)
 end do
 if ((time.eq.0.d0).or.(ivis.gt.0)) then
  call Compute_Matrix_coarse
  call DGBTRF(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse,2*kl+ku+1,IPIV,INFO) !!compute LU factors
 end if
 do level=1,grids
  call iterate_local_time_steps(level)
 end do
end if

q=1
do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!begin multigrid cycle
!!!!!!!!!!!!!!!!!!!! For non-Newtonian rheology
if (non_Newtonian.eq.1) then
  SF(:,:)=error(1,:,:)
  call compute_velocities
  call compute_strain_rate_invariant
 do level=1,grids !!compute viscosity and viscosity dependent local time steps for each multigrid level
  call compute_viscosity(level)
  call iterate_local_time_steps(level)
 end do
  call Compute_Matrix_coarse
  call DGBTRF(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse,2*kl+ku+1,IPIV,INFO) !!compute LU factors
end if
!!!!!!!!!!!!!!!!!!!! End for non-Newtoninan rheology
 do p=1,nsteps
  level=multi(p)
  if (multi(p+1).gt.multi(p)) then
   dir=0
  else
   dir=1 
  end if
  if (multi(p-1).lt.multi(p)) then
   call compute_RHS(level)
   if ((p.le.grids).and.(q.eq.1)) then
    call compute_reference_residual(level)
   end if
  end if
  if (level.lt.grids) then
   call iterate_stream(level,Nmulti,dir,p)
   if (level.eq.1) write(*,'(i5,2(g20.8))') q,residual_mag(p),res_ref(level)
!  write(*,*) q,level,p,residual_mag(p),res_ref(level),residual_mag(p)/res_ref(level)
  else
   !!!!!!!!!!!!!!!set up linear system
   call Compute_coarse_vector
   call DGBTRS('N',ngrid_coarse,kl,ku,1,Matrix_coarse,2*kl+ku+1,IPIV,B_coarse,ngrid_coarse,INFO)
   call Stream_grid_coarse
   call enforceBCs_stream(level,error(level,:,:))
!   call compute_residual(level,error(level,:,:),p)
   !!!!!!!!!!!!!!!end set up linear system
!   write(*,*) q,level,p,residual_mag(p),res_ref(level),residual_mag(p)/res_ref(level)
  end if
  if (level.eq.1) then
   if (non_Newtonian.eq.1) then !!make sure that viscosity is consistent with current stream function before judging convergence
    SF(:,:)=error(1,:,:)
    call compute_velocities
    call compute_strain_rate_invariant
    call compute_viscosity(level)
    call compute_residual(level,error(level,:,:),p)
   end if
   if ((residual_mag(1)/res_ref(1)).le.tolerance) exit
  end if
  if (multi(p+1).lt.multi(p)) then
   call prolong(level)
  end if
 end do
 if ((residual_mag(1)/res_ref(1)).le.tolerance) exit
!!!!!!!!!!!!!!!!!!!!!!!end multigrid cycle
q=q+1
 if (q.gt.max_cycles) then
  write(*,*) "Multigrid Convergence Problem: Too Many Multigrid Cycles Needed"
  stop
 end if
end do
 SF(:,:)=error(1,:,:)
 write(*,'(2(i5),2(g20.8))') tstep,q,residual_mag(1),res_ref(1)
!$ time2=omp_get_wtime()
!$ tmultigrid=tmultigrid+(time2-time1)
end

subroutine prolong(level_coarse)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: c1,c2,c3,c4
integer*4 :: level_coarse,i,k,lc,lf
lc=level_coarse
lf=level_coarse-1
!$OMP PARALLEL DO PRIVATE(i,k)
do i=1,nx_grid(lc)-1 !!prolong
 do k=1,nz_grid(lc)-1
!   error(lf,2*i,2*k)=error(lf,2*i,2*k)+error(lc,i,k)*4.d0/9.d0   !!arithmetic avg
   error(lf,2*i,2*k)=error(lf,2*i,2*k)+error(lc,i,k)              !!bilinear
 end do
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,k,c1,c2)
do i=2,nx_grid(lf)-2,2
 do k=1,nz_grid(lf)-1,2
!   error(lf,i,k)=error(lf,i,k)+(error(lc,i/2,(k-1)/2)+error(lc,i/2,(k+1)/2))*4.d0/9.d0 !!arithmetic avg
   error(lf,i,k)=error(lf,i,k)+(error(lc,i/2,(k-1)/2)+error(lc,i/2,(k+1)/2))/2.d0     !!bilinear
!   c1=min(vis_grid(lf,i,k-1),vis_grid(lf,i,k))/max(vis_grid(lf,i,k-1),vis_grid(lf,i,k))
!   c2=min(vis_grid(lf,i,k+1),vis_grid(lf,i,k))/max(vis_grid(lf,i,k+1),vis_grid(lf,i,k))
!   error(lf,i,k)=error(lf,i,k)+(c1*error(lc,i/2,(k-1)/2)+c2*error(lc,i/2,(k+1)/2))/(c1+c2)     !!viscosity dependent
 end do
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,k,c1,c2)
do k=2,nz_grid(lf)-2,2
 do i=1,nx_grid(lf)-1,2
!   error(lf,i,k)=error(lf,i,k)+(error(lc,(i-1)/2,k/2)+error(lc,(i+1)/2,k/2))*4.d0/9.d0 !!arithmetic avg
   error(lf,i,k)=error(lf,i,k)+(error(lc,(i-1)/2,k/2)+error(lc,(i+1)/2,k/2))/2.d0     !!bilinear
!   c1=min(vis_grid(lf,(i-1),k),vis_grid(lf,i,k))/max(vis_grid(lf,(i-1),k),vis_grid(lf,i,k))
!   c2=min(vis_grid(lf,(i+1),k),vis_grid(lf,i,k))/max(vis_grid(lf,(i+1),k),vis_grid(lf,i,k))
!   error(lf,i,k)=error(lf,i,k)+(c1*error(lc,(i-1)/2,k/2)+c2*error(lc,(i+1)/2,k/2))/(c1+c2)     !!viscosity dependent
 end do
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,k,c1,c2,c3,c4)
do i=1,nx_grid(lf)-1,2
 do k=1,nz_grid(lf)-1,2
!   error(lf,i,k)=error(lf,i,k)+&
!&(error(lc,(i-1)/2,(k-1)/2)+error(lc,(i+1)/2,(k-1)/2)+error(lc,(i+1)/2,(k+1)/2)+error(lc,(i-1)/2,(k+1)/2))*4.d0/9.d0 !!arithmetic avg
   error(lf,i,k)=error(lf,i,k)+&
&(error(lc,(i-1)/2,(k-1)/2)+error(lc,(i+1)/2,(k-1)/2)+error(lc,(i+1)/2,(k+1)/2)+error(lc,(i-1)/2,(k+1)/2))/4.d0 !!bilinear

!   c1=min(vis_grid(lf,(i-1),(k-1)),vis_grid(lf,i,k))/max(vis_grid(lf,(i-1),(k-1)),vis_grid(lf,i,k))
!   c2=min(vis_grid(lf,(i+1),(k-1)),vis_grid(lf,i,k))/max(vis_grid(lf,(i+1),(k-1)),vis_grid(lf,i,k))
!   c3=min(vis_grid(lf,(i+1),(k+1)),vis_grid(lf,i,k))/max(vis_grid(lf,(i+1),(k+1)),vis_grid(lf,i,k))
!   c4=min(vis_grid(lf,(i-1),(k+1)),vis_grid(lf,i,k))/max(vis_grid(lf,(i-1),(k+1)),vis_grid(lf,i,k))
!   error(lf,i,k)=error(lf,i,k)+&
!&(c1*error(lc,(i-1)/2,(k-1)/2)+c2*error(lc,(i+1)/2,(k-1)/2)+c3*error(lc,(i+1)/2,(k+1)/2)+c4*error(lc,(i-1)/2,(k+1)/2))/&
!&(c1+c2+c3+c4) !!viscosity dependent
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(lf,error(lf,:,:))
end

subroutine scale_factor(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,iq,kq,level
real*8 :: stream_space,factor,temp,Sval
real*8 :: dot(1:nx_grid(level)-1,1:nz_grid(level)-1),mag(1:nx_grid(level)-1,1:nz_grid(level)-1)
scaling=1

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,error(level,:,:),level,iq,kq,Sval)
  mag(i,k)=temp
  dot(i,k)=RHS(level,i,k)*temp
 end do
end do
!$OMP END PARALLEL DO
factor=sum(dot)/sum(mag**2.d0)
write(*,*) level,factor
scaling=0
end

subroutine compute_residual(level,sol,p)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,p,iq,kq
real*8 ::sol(-span:nx+span,-span:nz+span),stream_space,Sval
!$OMP PARALLEL DO PRIVATE(i,k)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  residual(i,k)=stream_space(i,k,sol,level,iq,kq,Sval)
 end do
end do
!$OMP END PARALLEL DO
residual_mag(p)=sum(dabs(residual(1:nx_grid(level)-1,1:nz_grid(level)-1)))/real(nx_grid(level)-1,8)/real(nz_grid(level)-1,8)

!write(*,*) level,residual_mag(p)
if ((level.eq.1).and.(residual_mag(p).gt.(1.d10*res_ref(level)))) then
 write(*,*) "Multigrid Convergence Problem: Abnormally High Residual Detected -- stopping"
 stop
end if
end

subroutine compute_reference_residual(level)
 use basics
 use arrays
implicit none
integer*4 :: level
res_ref(level)=sum(dabs(RHS(level,1:nx_grid(level)-1,1:nz_grid(level)-1)))/real(nx_grid(level)-2,8)/real(nz_grid(level)-2,8)
!write(*,*) res_ref
end

subroutine viscosity_gradients(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,ii,kk
real*8 :: visf,vis_r,vis_rr,vis_s,vis_ss,vis_rs

if (ivis.eq.0) then
   vis_r=0.d0
   vis_rr=0.d0
   vis_s=0.d0
   vis_ss=0.d0
   vis_rs=0.d0
   vis_x(level,:,:)=0.d0
   vis_z(level,:,:)=0.d0
   vis_xx(level,:,:)=0.d0
   vis_zz(level,:,:)=0.d0
   vis_xz(level,:,:)=0.d0
else
!$OMP PARALLEL DO PRIVATE(i,k,vis_r,vis_rr,vis_s,vis_ss,vis_rs,ii,kk)
do i=1,nx_grid(level)-1
 ii=i*2**(level-1)
 do k=1,nz_grid(level)-1
   kk=k*2**(level-1)
   vis_r=dot_product(D1(-span1:span1),vis_grid(level,i-span1:i+span1,k))/dr_grid(level)
   vis_rr=dot_product(D2(-span1:span1),vis_grid(level,i-span1:i+span1,k))/dr2_grid(level)
   vis_s=dot_product(D1(-span1:span1),vis_grid(level,i,k-span1:k+span1))/ds_grid(level)
   vis_ss=dot_product(D2(-span1:span1),vis_grid(level,i,k-span1:k+span1))/ds2_grid(level)
   vis_rs=sum(Drs(-span1:span1,-span1:span1)*vis_grid(level,i-span1:i+span1,k-span1:k+span1))/drds_grid(level)
   vis_x(level,i,k)=vis_r/xr(ii)
   vis_z(level,i,k)=vis_s/zs(kk)
   vis_xx(level,i,k)=(vis_rr-vis_x(level,i,k)*xrr(ii))/xr(ii)**2.d0
   vis_zz(level,i,k)=(vis_ss-vis_z(level,i,k)*zss(kk))/zs(kk)**2.d0
   vis_xz(level,i,k)=vis_rs/xr(ii)/zs(kk)
 end do
end do
!$OMP END PARALLEL DO
end if
end

subroutine compute_viscosity(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,ii,kk,q,i1,k1,interpolate
real*8 :: visf
if (level.eq.1) then

!$OMP PARALLEL DO PRIVATE(i,k)
do i=-span,nx_grid(level)+span
 do k=-span,nz_grid(level)+span
  vis_grid(level,i,k)=visf(i,k)
 end do
end do
!$OMP END PARALLEL DO
!  call viscosity_gradients(level)

  interpolate=0
  call smoother_vis(level,interpolate)  !!smooth computational viscosity to unresolved frequencies
else
  vis_gridf(level-1,:,:)=vis_grid(level-1,:,:) !!!store original viscosity
  vis_xf(level-1,:,:)=vis_x(level-1,:,:)
  vis_zf(level-1,:,:)=vis_z(level-1,:,:)
  vis_xxf(level-1,:,:)=vis_xx(level-1,:,:)
  vis_zzf(level-1,:,:)=vis_zz(level-1,:,:)
  vis_xzf(level-1,:,:)=vis_xz(level-1,:,:)
  interpolate=1
  call smoother_vis(level-1,interpolate) !!!!!smooth viscosity field before interpolating to coarse grid

 call restrict_coefficients(level,vis_grid(level-1,1:nx-1,1:nz-1),vis_grid(level,1:nx-1,1:nz-1))
 call viscosityBCs_initialize(level)
 if (level.lt.grids) then
  interpolate=0
  call smoother_vis(level,interpolate)  !!smooth computational viscosity to remove high frequencies
 else
  call viscosity_gradients(level)
 end if

  vis_grid(level-1,:,:)=vis_gridf(level-1,:,:) !!restore original viscosity
  vis_x(level-1,:,:)=vis_xf(level-1,:,:)
  vis_z(level-1,:,:)=vis_zf(level-1,:,:)
  vis_xx(level-1,:,:)=vis_xxf(level-1,:,:)
  vis_zz(level-1,:,:)=vis_zzf(level-1,:,:)
  vis_xz(level-1,:,:)=vis_xzf(level-1,:,:)
end if
end

subroutine viscosityBCs_initialize(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,ii,kk,level

if (level.gt.1) then
do k=1,nz_grid(level)-1 !!vertical boundaries
 kk=2*k !!fine grid
 vis_grid(level,0,k)=vis_grid(level-1,0,kk)
 vis_grid(level,nx_grid(level),k)=vis_grid(level-1,nx_grid(level-1),kk)
end do

do i=0,nx_grid(level) !!horizontal boundaries
 ii=2*i !!fine grid
 vis_grid(level,i,0)=vis_grid(level-1,ii,0)
 vis_grid(level,i,nz_grid(level))=vis_grid(level-1,ii,nz_grid(level-1))
end do
end if

do i=1,span !!!side ghost points
 vis_grid(level,-i,0:nz_grid(level))=vis_grid(level,i,0:nz_grid(level))
 vis_grid(level,nx_grid(level)+i,0:nz_grid(level))=vis_grid(level,nx_grid(level)-i,0:nz_grid(level))
end do

if (RaT.ne.0.d0) then
 do k=1,span !!!top/bottom ghost points !!based on T/C BCs
  vis_grid(level,-span:nx_grid(level)+span,-k)=1.d0/vis_grid(level,-span:nx_grid(level)+span,k)/visT**2.d0
  vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)+k)=1.d0/vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)-k)
 end do
else
 do k=1,span !!!top/bottom ghost points !!based on T/C BCs
  vis_grid(level,-span:nx_grid(level)+span,-k)=vis_grid(level,-span:nx_grid(level)+span,k)
  vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)+k)=vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)-k)
 end do
end if
end

subroutine viscosityBCs(level)
 use basics
 use arrays
implicit none
integer*4 :: i,k,ii,kk,level

if (level.gt.1) then
!do k=1,nz_grid(level)-1 !!vertical boundaries
! kk=2*k !!fine grid
! vis_grid(level,0,k)=vis_grid(level-1,0,kk)
! vis_grid(level,nx_grid(level),k)=vis_grid(level-1,nx_grid(level-1),kk)
!end do

do i=0,nx_grid(level) !!horizontal boundaries
 ii=2*i !!fine grid
 vis_grid(level,i,0)=vis_grid(level-1,ii,0)
 vis_grid(level,i,nz_grid(level))=vis_grid(level-1,ii,nz_grid(level-1))
end do
end if

do i=1,span !!!side ghost points: insulating
 vis_grid(level,-i,0:nz_grid(level))=vis_grid(level,i,0:nz_grid(level))
 vis_grid(level,nx_grid(level)+i,0:nz_grid(level))=vis_grid(level,nx_grid(level)-i,0:nz_grid(level))
end do

if (RaT.ne.0.d0) then
 do k=1,span !!!top/bottom ghost points !!based on T/C BCs
  vis_grid(level,-span:nx_grid(level)+span,-k)=1.d0/vis_grid(level,-span:nx_grid(level)+span,k)/visT**2.d0
  vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)+k)=1.d0/vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)-k)
 end do
else
 do k=1,span !!!top/bottom ghost points !!based on T/C BCs
  vis_grid(level,-span:nx_grid(level)+span,-k)=vis_grid(level,-span:nx_grid(level)+span,k)
  vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)+k)=vis_grid(level,-span:nx_grid(level)+span,nz_grid(level)-k)
 end do
end if
end


real*8 function stream_space(i,k,S0,level,iq,kq,Sval)
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,level,i1,k1,di,dk,iq,kq
real*8 :: Tr,Cr
real*8 :: Tx,Cx
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz
real*8 :: S0(-span:nx+span,-span:nz+span)
real*8 :: term1,term2,term3,term4,term5,term6,Sval

i1=i*2**(level-1)  !!!fine grid coordinates
k1=k*2**(level-1)

if (build.eq.0) then !!for smoother iterations
 Sr=dot_product(D1(-span1:span1),S0(i-span1:i+span1,k))/dr_grid(level)
 Ss=dot_product(D1(-span1:span1),S0(i,k-span1:k+span1))/ds_grid(level)
 Srr=dot_product(D2(-span1:span1),S0(i-span1:i+span1,k))/dr2_grid(level)
 Sss=dot_product(D2(-span1:span1),S0(i,k-span1:k+span1))/ds2_grid(level)
 Srrr=dot_product(D3(-span:span),S0(i-span:i+span,k))/dr3_grid(level)
 Ssss=dot_product(D3(-span:span),S0(i,k-span:k+span))/ds3_grid(level)
 Srrrr=dot_product(D4(-span:span),S0(i-span:i+span,k))/dr4_grid(level)
 Sssss=dot_product(D4(-span:span),S0(i,k-span:k+span))/ds4_grid(level)

 Srs=sum(Drs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds_grid(level)
 Srss=sum(Drss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/drds2_grid(level)
 Srrs=sum(Drrs(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds_grid(level)
 Srrss=sum(Drrss(-span1:span1,-span1:span1)*S0(i-span1:i+span1,k-span1:k+span1))/dr2ds2_grid(level)
elseif (build.eq.1) then !!for building the coarse matrix: extract matrix operator coefficients
 di=iq-i
 dk=kq-k

 if (k.eq.kq) then
  Sr=D1(di)*Sval/dr_grid(level)
  Srr=D2(di)*Sval/dr2_grid(level)
  Srrr=D3(di)*Sval/dr3_grid(level)
  Srrrr=D4(di)*Sval/dr4_grid(level)
 else
  Sr=0.d0
  Srr=0.d0
  Srrr=0.d0
  Srrrr=0.d0
 end if

 if (i.eq.iq) then
  Ss=D1(dk)*Sval/ds_grid(level)
  Sss=D2(dk)*Sval/ds2_grid(level)
  Ssss=D3(dk)*Sval/ds3_grid(level)
  Sssss=D4(dk)*Sval/ds4_grid(level)
 else
  Ss=0.d0
  Sss=0.d0
  Ssss=0.d0
  Sssss=0.d0
 end if

 Srs=Drs(di,dk)*Sval/drds_grid(level)
 Srss=Drss(di,dk)*Sval/drds2_grid(level)
 Srrs=Drrs(di,dk)*Sval/dr2ds_grid(level)
 Srrss=Drrss(di,dk)*Sval/dr2ds2_grid(level)
end if

 Sx=Sr/xr(i1)
 Sz=Ss/zs(k1)
 Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
 Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
 Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
 Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
 Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
 Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*Szz-Sz*zssss(k1))/zs(k1)**4.d0

 Sxz=Srs/xr(i1)/zs(k1)
 Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
 Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
 Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

if (ivis.eq.0) then
 term1=0.d0
 term2=vis_grid(level,i,k)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(level,i,k)*Sxxzz
 term4=0.d0
 term5=0.d0
 term6=0.d0
else
 term1=(vis_xx(level,i,k)-vis_zz(level,i,k))*(Sxx-Szz)
 term2=vis_grid(level,i,k)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(level,i,k)*Sxxzz
 term4=4.d0*vis_xz(level,i,k)*Sxz
 term5=2.d0*vis_x(level,i,k)*(Sxxx+Sxzz)
 term6=2.d0*vis_z(level,i,k)*(Szzz+Sxxz)
end if

if (build.eq.0) then !!for computing residuals during smoother iterations
 if (scaling.eq.0) then !!for smoother iterations
  stream_space=term1+term2+term3+term4+term5+term6-RHS(level,i,k) !!add terms
  stream_space=-stream_space
 elseif (scaling.eq.1) then !!for computing scaling factor for multigrid convergence
  stream_space=term1+term2+term3+term4+term5+term6
 end if
elseif (build.eq.1) then !!for constructing the matrix on the coarsest grid
 stream_space=term1+term2+term3+term4+term5+term6
end if
end

subroutine compute_RHS(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,j,k,level,ii,kk
real*8 :: Tr,Tx,Cr,Cx
real*8 :: a1,a2,a3
if (level.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,k,Tr,Tx,Cr,Cx,ii,kk,a1,a2,a3)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   RHS(level,i,k)=0.d0
   if (RaT.ne.0) then
    Tr=dot_product(D1(-span1:span1),Tbuoy(i-span1:i+span1,k))/dr
    Tx=Tr/xr(i)
    RHS(level,i,k)=RaT*Tx
   end if
   if (comp.eq.1) then !!add in compositional buoyancy
    do j=1,ntype
     Cr=dot_product(D1(-span1:span1),Cbuoy(j,i-span1:i+span1,k))/dr
     Cx=Cr/xr(i)
     RHS(level,i,k)=RHS(level,i,k)-RaC(j)*Cx
    end do
   end if
 end do
end do
!$OMP END PARALLEL DO
else
 call restrict(level,residual,RHS(level,:,:))
end if
end

subroutine restrict(level,array_in,array_out)
use basics
use arrays
implicit none
integer*4 :: i,k,ii,kk,level,di,dk
real*8 :: corners,sides,centre,array_in(1:nx-1,1:nz-1),array_out(1:nx-1,1:nz-1),c(-1:1,-1:1)
!$OMP PARALLEL DO PRIVATE(i,k,ii,kk,corners,sides,centre,c,di,dk)
do i=1,nx_grid(level)-1
 ii=2*i !!fine grid coordinates from one level up
 do k=1,nz_grid(level)-1
   kk=2*k
    corners=(array_in(ii-1,kk-1)+array_in(ii-1,kk+1)+array_in(ii+1,kk+1)+array_in(ii-1,kk+1))/16.d0
    sides=(array_in(ii-1,kk)+array_in(ii,kk-1)+array_in(ii+1,kk)+array_in(ii,kk+1))/8.d0
    centre=array_in(ii,kk)/4.d0
    array_out(i,k)=(corners+sides+centre) !!bilinear

!    do di=-1,1
!     do dk=-1,1
!      c(di,dk)=min(vis_grid(level-1,ii+di,kk+dk),vis_grid(level-1,ii,kk))/&
!      &max(vis_grid(level-1,ii+di,kk+dk),vis_grid(level-1,ii,kk))
!     end do
!    end do
!    array_out(i,k)=sum(c*array_in(ii-1:ii+1,kk-1:kk+1))/sum(c) !!viscosity dependent

!    array_out(i,k)=sum(array_in(ii-1:ii+1,kk-1:kk+1)/9.d0)   !!arithmetic avg    

!    array_out(i,k)=array_in(ii,kk) !!simple injection    
 end do
end do
!$OMP END PARALLEL DO
end

subroutine restrict_coefficients(level,array_in,array_out)
use basics
use arrays
implicit none
integer*4 :: i,k,ii,kk,level,p,nsample,di,dk
real*8 :: corners,sides,centre,array_in(1:nx-1,1:nz-1),array_out(1:nx-1,1:nz-1),c(-1:1,-1:1)
!$OMP PARALLEL DO PRIVATE(i,k,ii,kk,corners,sides,centre,di,dk,c)
do i=1,nx_grid(level)-1
 ii=2*i !!fine grid coordinates from one level up
 do k=1,nz_grid(level)-1
   kk=2*k
!    corners=(array_in(ii-1,kk-1)*array_in(ii-1,kk+1)*array_in(ii+1,kk+1)*array_in(ii-1,kk+1))**(1.d0/16.d0)
!    sides=(array_in(ii-1,kk)*array_in(ii,kk-1)*array_in(ii+1,kk)*array_in(ii,kk+1))**(1.d0/8.d0)
!    centre=array_in(ii,kk)**(1.d0/4.d0)
!    array_out(i,k)=corners*sides*centre !!exponential

    corners=(array_in(ii-1,kk-1)+array_in(ii-1,kk+1)+array_in(ii+1,kk+1)+array_in(ii-1,kk+1))/16.d0
    sides=(array_in(ii-1,kk)+array_in(ii,kk-1)+array_in(ii+1,kk)+array_in(ii,kk+1))/8.d0
    centre=array_in(ii,kk)/4.d0
    array_out(i,k)=(corners+sides+centre)

!    do di=-1,1
!     do dk=-1,1
!      c(di,dk)=min(vis_grid(level-1,ii+di,kk+dk),vis_grid(level-1,ii,kk))/&
!      &max(vis_grid(level-1,ii+di,kk+dk),vis_grid(level-1,ii,kk))
!     end do
!    end do
!    array_out(i,k)=sum(c*array_in(ii-1:ii+1,kk-1:kk+1))/sum(c) !!viscosity dependent


!     array_out(i,k)=sum(array_in(ii-1:ii+1,kk-1:kk+1)/9.d0)   !!arithmetic avg    

!     array_out(i,k)=array_in(ii,kk) !!simple injection
 end do
end do
!$OMP END PARALLEL DO
end

subroutine test_eigenvalues(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: eI,eR
real*8 :: eigI,eigR
Ts_max=Tr_base/ds_grid(level)     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary
Tsss_max=Trrr_base/ds3_grid(level) !!for real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tssss_max=Trrrr_base/ds4_grid(level)
Tr_max=Tr_base/dr_grid(level)
Trr_max=-Trr_base/dr2_grid(level)
Trrr_max=Trrr_base/dr3_grid(level)
Trrrr_max=Trrrr_base/dr4_grid(level)
Trs_max=-Trs_base/drds_grid(level)
Trrss_max=Trrss_base/dr2ds2_grid(level)
Trss_max=Trrs_base/drds2_grid(level)
Trrs_max=Trrs_base/dr2ds_grid(level)

ratio=1.d99
if (ivis.eq.0) return
!$OMP PARALLEL DO PRIVATE(i,k,eI,eR)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   eI=eigI(level,i,k)
   eR=eigR(level,i,k)
   ratio(level,i,k)=dabs(eR/eI)
 end do
end do
!$OMP END PARALLEL DO
end

subroutine iterate_local_time_steps(level)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,level
real*8 :: dtA,dtD,eI,eR
real*8 :: eigI,eigR
!real*8 :: eR_array(1:grids,1:nx,1:nz)
Ts_max=Tr_base/ds_grid(level)     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2_grid(level) !!even # of derivative operators=real, odd #=imaginary
Tsss_max=Trrr_base/ds3_grid(level) !!for real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tssss_max=Trrrr_base/ds4_grid(level)
Tr_max=Tr_base/dr_grid(level)
Trr_max=-Trr_base/dr2_grid(level)
Trrr_max=Trrr_base/dr3_grid(level)
Trrrr_max=Trrrr_base/dr4_grid(level)
Trs_max=-Trs_base/drds_grid(level)
Trrss_max=Trrss_base/dr2ds2_grid(level)
Trss_max=Trrs_base/drds2_grid(level)
Trrs_max=Trrs_base/dr2ds_grid(level)

!eR_array=0.d0
!$OMP PARALLEL DO PRIVATE(i,k,dtA,dtD,eI,eR)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
   eI=eigI(level,i,k)
   if (iterate_order.eq.1) then
     eR=eigR(level,i,k)
!     eR_array(level,i,k)=eR
     dt_stream(level,i,k)=-courant_stream*2.d0*eR/(eR**2.d0+eI**2.d0)    
   elseif (iterate_order.eq.4) then
    dtD=-2.8d0/eigR(level,i,k)
    if ((dtD.gt.0.d0).and.(ei.ne.0.d0)) then
     dtA=2.8d0/eI
    else
     dtA=-2.8d0/eI
    end if
    if (ei.eq.0.d0) then
     dt_stream(level,i,k)=courant_stream*dtD
    else
     dt_stream(level,i,k)=courant_stream/(1.d0/dtA+1.d0/dtD)
    end if
   end if
 end do
end do
!$OMP END PARALLEL DO
!write(*,*) level,maxval(eR_array(level,:,:)),minval(eR_array(level,:,:))
end

real*8 function eigI(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: term1,term2,term3,term4,term5,term6
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz

i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
k1=k*2**(level-1)

!!!!!!!!!!!!!!!!!!!!!feed in imaginary parts of the eigenvalues (odd # of derivative operators)
Sr=Tr_max
Ss=Ts_max
Srr=0.d0
Sss=0.d0
Srrr=Trrr_max
Ssss=Tsss_max
Srrrr=0.d0
Sssss=0.d0
Srs=0.d0
Srrs=Trrs_max
Srss=Trss_max
Srrss=0.d0
!!!!!!!!!!!!!!!!!!!!

Sx=Sr/xr(i1)  !!!!!!!!!!!!!adding absolute values here to be conservative
Sz=Ss/zs(k1)
Sxx=(Srr+dabs(Sx*xrr(i1)))/xr(i1)**2.d0
Szz=(Sss+dabs(Sz*zss(k1)))/zs(k1)**2.d0
Sxxx=(Srrr+dabs(3.d0*Sxx*xr(i1)*xrr(i1))+dabs(Sx*xrrr(i1)))/dabs(xr(i1)**3.d0)
Szzz=(Ssss+dabs(3.d0*Szz*zs(k1)*zss(k1))+dabs(Sz*zsss(k1)))/dabs(zs(k1)**3.d0)
Sxxxx=(Srrrr+dabs(6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0)+(dabs(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*dabs(Sxx))+&
&dabs(Sx*xrrrr(i1)))/xr(i1)**4.d0
Szzzz=(Sssss+dabs(6.d0*Szzz*zss(k1)*zs(k1)**2.d0)+(dabs(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*dabs(Szz))+&
&dabs(Sz*zssss(k1)))/zs(k1)**4.d0

Sxz=Srs/xr(i1)/zs(k1)
Sxxz=(Srrs+dabs(Sxz*zs(k1)*xrr(i1)))/dabs(zs(k1))/xr(i1)**2.d0
Sxzz=(Srss+dabs(Sxz*xr(i1)*zss(k1)))/dabs(xr(i1))/zs(k1)**2.d0
Sxxzz=(Srrss+dabs(Sxzz*xrr(i1))*zs(k1)**2.d0+dabs(Sxxz*zss(k1))*xr(i1)**2.d0+dabs(Sxz*xrr(i1)*zss(k1)))/(xr(i1)*zs(k1))**2.d0

 term1=(vis_xx(level,i,k)+vis_zz(level,i,k))*(Sxx+Szz)
 term2=vis_grid(level,i,k)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(level,i,k)*Sxxzz
 term4=4.d0*vis_xz(level,i,k)*Sxz
 term5=2.d0*vis_x(level,i,k)*(Sxxx+Sxzz)
 term6=2.d0*vis_z(level,i,k)*(Szzz+Sxxz)

 eigI=dabs(term1)+dabs(term2)+dabs(term3)+dabs(term4)+dabs(term5)+dabs(term6)
end


real*8 function eigR(level,i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k,level,i1,k1
real*8 :: term1,term2,term3,term4,term5,term6
real*8 :: Sr,Ss,Srr,Sss,Srrr,Ssss,Srrrr,Sssss
real*8 :: Sx,Sz,Sxx,Szz,Sxxx,Szzz,Sxxxx,Szzzz
real*8 :: Srs,Srrs,Srss,Srrss
real*8 :: Sxz,Sxxzz,Sxxz,Sxzz

i1=i*2**(level-1) !!fine grid coordinates for grid metrics (may need to restrict the metrics eventually)
k1=k*2**(level-1)

!!!!!!!!!!!!!!!!!!!!!feed in real parts of the eigenvalues (even # of derivative operators)
Sr=0.d0
Ss=0.d0
Srr=Trr_max
Sss=Tss_max
Srrr=0.d0
Ssss=0.d0
Srrrr=Trrrr_max
Sssss=Tssss_max
Srs=Trs_max
Srrs=0.d0
Srss=0.d0
Srrss=Trrss_max
!!!!!!!!!!!!!!!!!!!!

Sx=Sr/xr(i1)
Sz=Ss/zs(k1)
Sxx=(Srr-Sx*xrr(i1))/xr(i1)**2.d0
Szz=(Sss-Sz*zss(k1))/zs(k1)**2.d0
Sxxx=(Srrr-3.d0*Sxx*xr(i1)*xrr(i1)-Sx*xrrr(i1))/xr(i1)**3.d0
Szzz=(Ssss-3.d0*Szz*zs(k1)*zss(k1)-Sz*zsss(k1))/zs(k1)**3.d0
Sxxxx=(Srrrr-6.d0*Sxxx*xrr(i1)*xr(i1)**2.d0-(3.d0*xrr(i1)**2.d0+4.d0*xr(i1)*xrrr(i1))*Sxx-Sx*xrrrr(i1))/xr(i1)**4.d0
Szzzz=(Sssss-6.d0*Szzz*zss(k1)*zs(k1)**2.d0-(3.d0*zss(k1)**2.d0+4.d0*zs(k1)*zsss(k1))*Szz-Sz*zssss(k1))/zs(k1)**4.d0

Sxz=Srs/xr(i1)/zs(k1)
Sxxz=(Srrs-Sxz*zs(k1)*xrr(i1))/zs(k1)/xr(i1)**2.d0
Sxzz=(Srss-Sxz*xr(i1)*zss(k1))/xr(i1)/zs(k1)**2.d0
Sxxzz=(Srrss-Sxzz*xrr(i1)*zs(k1)**2.d0-Sxxz*zss(k1)*xr(i1)**2.d0-Sxz*xrr(i1)*zss(k1))/(xr(i1)*zs(k1))**2.d0

 term1=(vis_xx(level,i,k)-vis_zz(level,i,k))*(Sxx-Szz)
 term2=vis_grid(level,i,k)*(Sxxxx+Szzzz)
 term3=2.d0*vis_grid(level,i,k)*Sxxzz
 term4=4.d0*vis_xz(level,i,k)*Sxz
 term5=2.d0*vis_x(level,i,k)*(Sxxx+Sxzz)
 term6=2.d0*vis_z(level,i,k)*(Szzz+Sxxz)

  eigR=-(term1+term2+term3+term4+term5+term6)
end


subroutine iterate_stream(level,Niter,dir,p)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID,n,Niter,tt,level,dir,nn,iq,kq,p
real*8 :: temp,la,lb
real*8 :: S0(-span:nx+span,-span:nz+span),S1(-span:nx+span,-span:nz+span),S2(-span:nx+span,-span:nz+span)
real*8 :: time1,time2,dto2,dto6
real*8 :: stream_space,c(1:nx,1:nz),cmin,ntemp,Sval

do n=1,Niter
 if (iterate_order.eq.1) then
  ntemp=Niter-n+1 !!!start with smaller time steps to kill off high frequencies first
  la=-1.d0

temp=dcos(real(2*ntemp-1,8)*pii/real(2*Niter,8))
!$OMP PARALLEL DO PRIVATE(i,k,lb)
  do i=1,nx_grid(level)-1
   do k=1,nz_grid(level)-1
    if (ratio(level,i,k).gt.1.d99) then
     lb=-2.d0
    else
     lb=-2.d0*ratio(level,i,k)**2.d0/(1.d0+ratio(level,i,k)**2.d0)
    end if
    c(i,k)=2.d0/(-lb-la+(lb-la)*temp)
   end do
  end do
!$OMP END PARALLEL DO


 elseif (iterate_order.eq.4) then
  cmin=1.6d0/2.8d0
  if (Niter.gt.1) then
   c=(cmin-1.d0)/real(1-Niter,8)*real(n-1,8)+cmin
  else
   c=cmin
  end if
 end if
if (dir.eq.1)  c=1.0d0
S0=error(level,:,:) !!initial stream function

if (iterate_order.eq.1) then !!Euler's method

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,S0,level,iq,kq,Sval)
  error(level,i,k)=S0(i,k)+c(i,k)*dt_stream(level,i,k)*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,error(level,:,:))

elseif (iterate_order.eq.4) then !!RK4

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,S0,level,iq,kq,Sval)
  error(level,i,k)=temp
  S1(i,k)=S0(i,k)+0.5d0*c(i,k)*dt_stream(level,i,k)*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,S1)

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,S1,level,iq,kq,Sval)
  error(level,i,k)=error(level,i,k)+2.d0*temp
  S2(i,k)=S0(i,k)+0.5d0*c(i,k)*dt_stream(level,i,k)*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,S2)

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,S2,level,iq,kq,Sval)
  error(level,i,k)=error(level,i,k)+2.d0*temp
  S1(i,k)=S0(i,k)+c(i,k)*dt_stream(level,i,k)*temp
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,S1)

!$OMP PARALLEL DO PRIVATE(i,k,temp)
do i=1,nx_grid(level)-1
 do k=1,nz_grid(level)-1
  temp=stream_space(i,k,S1,level,iq,kq,Sval)
  error(level,i,k)=S0(i,k)+(c(i,k)*dt_stream(level,i,k)/6.d0)*(error(level,i,k)+temp)
 end do
end do
!$OMP END PARALLEL DO
 call enforceBCs_stream(level,error(level,:,:))

end if

end do
 call compute_residual(level,error(level,:,:),p)
end

subroutine enforceBCs_stream(level,S0) !!!!update this
 use basics
implicit none
integer*4 :: i,level
real*8 :: S0(-span:nx+span,-span:nz+span)
S0(0,-span:nz_grid(level)+span)=0.d0   !!impermeable boundaries
S0(nx_grid(level),-span:nz_grid(level)+span)=0.d0
S0(-span:nx_grid(level)+span,0)=0.d0
S0(-span:nx_grid(level)+span,nz_grid(level))=0.d0
 do i=1,span  !!!antisymmetry for free-slip sidewalls
  S0(-i,0:nz_grid(level))=-S0(i,0:nz_grid(level))      !!antsymmetry: free-slip
  S0(nx_grid(level)+i,0:nz_grid(level))=-S0(nx_grid(level)-i,0:nz_grid(level))
 end do
if (Vbc.eq.1) then !!symmetry for rigid top/bottom
 do i=1,span
  S0(-span:nx_grid(level)+span,-i)=S0(-span:nx_grid(level)+span,i)
  S0(-span:nx_grid(level)+span,nz_grid(level)+i)=S0(-span:nx_grid(level)+span,nz_grid(level)-i)
 end do
elseif (Vbc.eq.0) then !!antisymmetry for free-slip top/bottom
 do i=1,span
  S0(-span:nx_grid(level)+span,-i)=-S0(-span:nx_grid(level)+span,i)
  S0(-span:nx_grid(level)+span,nz_grid(level)+i)=-S0(-span:nx_grid(level)+span,nz_grid(level)-i)
 end do
end if
end

integer*4 function pf_coarse(ip,kp) !!input is ip,kp -- output is p value
 use basics
 use arrays
implicit none
integer*4 :: ip,kp,N
if (Database.eq.1) then
 N=nx_grid(grids)+1
 pf_coarse=N*kp+ip+1
elseif (Database.eq.2) then
 N=nz_grid(grids)+1
 pf_coarse=N*ip+kp+1
elseif (Database.eq.3) then !!update or take out
 pf_coarse=parray(ip,kp)
end if
end

subroutine indices_coarse(p,ip,kp)  !!input is p, output is ip,kp
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp,N
if (Database.eq.1) then
 N=nx_grid(grids)+1
 ip=p-N*((p-1)/N)-1
 kp=(p-1)/N
elseif (Database.eq.2) then
 N=nz_grid(grids)+1
 kp=p-N*((p-1)/N)-1
 ip=(p-1)/N
elseif (Database.eq.3) then !!!!!update or remove
 ip=ipvec(p)
 kp=kpvec(p)
end if
end

subroutine Compute_Matrix_coarse
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: ip,kp,p,iq,kq,q,pf_coarse,iiq,kkq,INFO,limit(-span:span),j,temp
integer*4 :: nx_coarse,nz_coarse,kl_ku_1
real*8 :: stream_space,xval,zval,S0(-span:nx+span,-span:nz+span)
real*8 :: ROWCND,COLCND,AMAX,ta1,ta2,Sval
!$ ta1=omp_get_wtime()

nx_coarse=nx_grid(grids)
nz_coarse=nz_grid(grids)

!!!!!!!!!!!!!!!!!define the stencil shape for shorter loops
limit(-span)=0
limit(-span1:-1)=span1
limit(0)=span
limit(1:span1)=span1
limit(span)=0
!!!!!!!!!!!!!!!! end stencil

if (analysis.eq.1) then !!!figure of matrix bandwidth on startup
ku=0
kl=0
!$OMP PARALLEL DO PRIVATE(p,q,temp,j,iq,kq,ip,kp) REDUCTION(max:ku,kl)
 do p=1,ngrid_coarse
  call indices_coarse(p,ip,kp)
   do j=-span,span
    iq=j+ip
    if ((iq.gt.0).and.(iq.lt.nx_coarse)) then
     do kq=-limit(j)+kp,kp+limit(j)
       if ((kq.gt.0).and.(kq.lt.nz_coarse)) then
        q=pf_coarse(iq,kq)
        temp=q-p
        if (temp.gt.ku) ku=temp
        if (-temp.gt.kl) kl=-temp
       end if
     end do
    end if
   end do
 end do
!$OMP END PARALLEL DO
write(*,*) "Coarse Matrix:",ku,kl,ngrid_coarse
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct Matrix Once Bandwidth is known
if (analysis.eq.0) then
kl_ku_1=kl+ku+1
build=1 !!tells stream_space function not to include RHS terms in the matrix elements
Matrix_coarse=0.d0
!$OMP PARALLEL DO PRIVATE(p,q,j,ip,kp,iq,kq,iiq,kkq,xval,zval,Sval)
do p=1,ngrid_coarse
 call indices_coarse(p,ip,kp)
  if ((ip.gt.0).and.(ip.lt.nx_coarse).and.(kp.gt.0).and.(kp.lt.nz_coarse)) then  !!!for interior values of stream function
   do j=-span,span
    iq=j+ip
    if ((iq.ge.0).and.(iq.le.nx_coarse)) then
      iiq=iq
      xval=1.d0
     do kq=-limit(j)+kp,kp+limit(j)
      if ((kq.ge.0).and.(kq.le.nz_coarse)) then
       kkq=kq
       zval=1.d0
       Sval=1.d0
       q=pf_coarse(iiq,kkq)
       Matrix_coarse(kl_ku_1+p-q,q)=stream_space(ip,kp,S0,grids,iq,kq,Sval)
      end if
     end do
    end if
   end do
  end if
end do
!$OMP END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!handle boundary conditions
do p=1,ngrid_coarse
 call indices_coarse(p,ip,kp)
 if ((ip.gt.0).and.(ip.lt.nx_coarse).and.(kp.gt.0).and.(kp.lt.nz_coarse)) then  !!!for ghost values of the streamfunction
  do j=-span,span
   iq=j+ip
   do kq=-limit(j)+kp,kp+limit(j)
    if ((iq.lt.0).or.(iq.gt.nx_coarse).or.(kq.lt.0).or.(kq.gt.nz_coarse)) then
     if (iq.lt.0) then  !!free-slip left sidewall: antisymmetry
      iiq=-iq
      xval=-1.d0
     elseif (iq.gt.nx_coarse) then!!free-slip right sidewall: antisymmetry
      iiq=2*nx_coarse-iq
      xval=-1.d0
     else
      iiq=iq
      xval=1.d0
     end if

     if (kq.lt.0) then !!!bottom
      kkq=-kq
      if (Vbc.eq.0) then !!free-slip: antisymmetry
       zval=-1.d0
      elseif (Vbc.eq.1) then !!rigid: symmetry
       zval=1.d0
      end if
     elseif (kq.gt.nz_coarse) then  !!!!top
      kkq=2*nz_coarse-kq
      if (Vbc.eq.0) then !!free-slip: antisymmetry
       zval=-1.d0
      elseif (Vbc.eq.1) then !!rigid: symmetry
       zval=1.d0
      end if
     else
      kkq=kq
      zval=1.d0
     end if
     Sval=xval*zval
     q=pf_coarse(iiq,kkq)
     Matrix_coarse(kl_ku_1+p-q,q)=Matrix_coarse(kl_ku_1+p-q,q)+stream_space(ip,kp,S0,grids,iq,kq,Sval)
    end if
   end do
  end do
 end if
end do

!$OMP PARALLEL DO PRIVATE(kp,p)
do kp=0,nz_coarse                !!impermeable boundaries
   p=pf_coarse(0,kp)
   Matrix_coarse(kl_ku_1,p)=1.d0
   p=pf_coarse(nx_coarse,kp)
   Matrix_coarse(kl_ku_1,p)=1.d0
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(ip,p)
do ip=0,nx_coarse                !!impermeable boundaries
   p=pf_coarse(ip,0)
   Matrix_coarse(kl_ku_1,p)=1.d0
   p=pf_coarse(ip,nz_coarse)
   Matrix_coarse(kl_ku_1,p)=1.d0
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!end handle boundary conditions

 call DGBEQU(ngrid_coarse,ngrid_coarse,kl,ku,Matrix_coarse(kl+1:2*kl+ku+1,1:ngrid_coarse),kl+ku+1,RowB,ColB,ROWCND,COLCND,AMAX,INFO)
! write(*,*) INFO,AMAX,ROWCND,COLCND
!$OMP PARALLEL DO PRIVATE(q,p)
 do q=1,ngrid_coarse  !!scale matrix
  do p=max(1,q-kl),min(ngrid_coarse,q+ku)
   Matrix_coarse(kl_ku_1+p-q,q)=RowB(p)*Matrix_coarse(kl_ku_1+p-q,q)
  end do
 end do
!$OMP END PARALLEL DO

build=0 !!subsequent calls to stream_space will take RHS terms into account (as required for smoothing iterations)
end if
!$ ta2=omp_get_wtime()
!$ tcoarse_matrix=tcoarse_matrix+ta2-ta1
end

subroutine Compute_coarse_vector !!build RHS for linear system on coarsest grid
 use basics
 use arrays
implicit none
integer*4 :: p,i,k,pf_coarse
B_coarse=0.d0 !!ensure boundary values are zero to ensure stream function will be zero on model boundaries
do i=1,nx_grid(grids)-1
 do k=1,nz_grid(grids)-1
  p=pf_coarse(i,k)
  B_coarse(p)=RHS(grids,i,k)*RowB(p)
 end do
end do
end

subroutine Stream_grid_coarse
 use basics
 use arrays
implicit none
integer*4 :: p,ip,kp
do p=1,ngrid_coarse
 call indices_coarse(p,ip,kp)
 error(grids,ip,kp)=B_coarse(p)
end do
end
