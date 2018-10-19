module basics
integer*4 :: nx,nz,ngrid,order,span,span1,spanT,local_time,subD,supD,restart,Database,initial_temp
integer*4 :: Nsnap,Nout,fcount,Tbot,Vbc,tstep,RK_update,torder
integer*4 :: comp,tpc,ntr,initial_comp,Dnumber,span_interp
integer*4 :: tracer,equalize,ntype,kent,nequalize,add_layer,gridX,gridZ,nbound
integer*4 :: iterate,grids,iterate_order,ngrid_coarse,kl,ku,analysis,build,iq_save,kq_save
integer*4 :: mingrid,nsteps,initial_temp_smooth,smooth_vis,scaling,vis_iter
integer*4, allocatable :: nx_grid(:),nz_grid(:),multi(:)
integer*4 :: non_Newtonian
real*8, allocatable :: dr_grid(:),ds_grid(:),dr2_grid(:),ds2_grid(:),dr3_grid(:),ds3_grid(:)
real*8, allocatable :: dr4_grid(:),ds4_grid(:),drds_grid(:),drds2_grid(:),dr2ds_grid(:),dr2ds2_grid(:)
real*8 :: aspect,spaceX,spaceZ,dr,ds,dt,time,Ax,Az,tend,tout,tsnap
real*8 :: RaT,V,ivis,visT,visP,H,yield,yield_gradient,yield_depth
real*8 :: Tr_max,Trr_max,Trrr_max,Trrrr_max,Ts_max,Tss_max,Tsss_max,Tssss_max
real*8 :: Trs_max,Trrs_max,Trss_max,Trrss_max
real*8 :: Tr_base,Trr_base,Trrr_base,Trrrr_base,Ts_base,Tss_base,Tsss_base,Tssss_base
real*8 :: Trs_base,Trrs_base,Trss_base,Trrss_base
real*8 :: visr,visrr,viss,visss,visx,visxx,visz,viszz,visrs,visxz
real*8 :: dlayer,filter,dts,courant,courantT,courantC,zent,dent,courant_stream
real*8 :: tenergy,tcd,tbSF,tsSF,tequalize,tadvect_tracer,tconvert,tmultigrid
real*8 :: tcoarse_matrix
real*8, parameter :: pii=3.1415926535897932d0

real*8 :: dr2,ds2,dr3,ds3,dr4,ds4,drds,drds2,dr2ds,dr2ds2
real*8 :: dro2,dso2,dro10,dso10,dist_r,dist_s,dist_factor
real*8 :: thickness,dt_manual,tolerance,delta_T,delta_T_init,delta_visT,delta_visT_init,eig_ratio
integer*4 :: max_cycles,Nvis_store,nn_store,Nmulti
real*8 :: diff_max,eig_ratio_min
real*8, allocatable :: res_ref(:),residual_mag(:),delta_visC(:),delta_C(:)
real*8 :: eig_fine
integer*4 :: n_eta,n_temp
integer*4, allocatable :: n_comp(:)
integer*4, allocatable :: seed_i(:),seed_j(:),seed_k(:),seed_n(:)
real*8, allocatable :: Cinit(:)
end module basics

subroutine storage
 use basics
 use arrays
implicit none
!allocate(IPIV(1:ngrid),Row(1:ngrid),Col(1:ngrid),Matrix(1:ngrid,1:ngrid))
allocate(Bvec(1:ngrid),vis(-span:nx+span,-span:nz+span))
allocate(xg(-(span+span_interp):nx+span+span_interp),zg(-(span+span_interp):nz+span+span_interp))
allocate(T(-span:nx+span,-span:nz+span))
allocate(SF(-span:nx+span,-span:nz+span),u(-span_interp:nx+span_interp,-span_interp:nz+span_interp))
allocate(w(-span_interp:nx+span_interp,-span_interp:nz+span_interp))
allocate(D1(-span:span),D2(-span:span),D3(-span:span),D4(-span:span))
allocate(Drs(-span:span,-span:span),Drss(-span:span,-span:span),Drrs(-span:span,-span:span),Drrss(-span:span,-span:span))
allocate(xr(-span_interp:nx+span_interp),xrr(-span_interp:nx+span_interp),xrrr(0:nx),xrrrr(0:nx))
allocate(zs(-span_interp:nz+span_interp),zss(-span_interp:nz+span_interp),zsss(0:nz),zssss(0:nz))
allocate(dt_array(0:nx,0:nz),dts_array(0:nx,0:nz))
if (Database.eq.3) allocate(parray(0:nx,0:nz),ipvec(1:ngrid),kpvec(1:ngrid))
allocate(D(1:Dnumber,-span_interp:span_interp))
allocate(Tvis(-span:nx+span,-span:nz+span),Tbuoy(-span:nx+span,-span:nz+span))
if (comp.eq.1) then
 allocate(Cvis(1:ntype,-span:nx+span,-span:nz+span),Cbuoy(1:ntype,-span:nx+span,-span:nz+span))
 allocate(Cnew(1:ntype,-span:nx+span,-span:nz+span)) !!assume type 0 is ambient -- do not compute explicitly
 Cnew=0.d0
 allocate(mass(1:ntype))
 allocate(n_comp(1:ntype))
end if
if (tracer.eq.1) then
 allocate(Ttr(1:ntr),Ttr0(1:ntr),Ttr1(1:ntr))
 if (torder.gt.2) allocate(Ttr2(1:ntr))
 allocate(Tratio(-span:nx+span,-span:nz+span))
 allocate(DERr_tracer(1:order-1,0:nx,-span_interp:nz+span_interp))
 allocate(tracer_space_array(-span_interp:nx+span_interp,-span_interp:nz+span_interp))
 allocate(Textend(-spanT:nx+spanT,-spanT:nz+spanT))
 allocate(DERr_Textend(1:order-1,0:nx,-span_interp:nz+span_interp))
end if
if ((comp.eq.1).or.(tracer.eq.1)) then
 allocate(Ttype(1:ntr),count_type(0:ntype,-span:nx+span,-span:nz+span))
 allocate(nempty(-span:nx+span,-span:nz+span))
 allocate(rtr(1:ntr),str(1:ntr))
 allocate(DERr_u(1:order-1,0:nx,-span_interp:nz+span_interp))
 allocate(DERr_w(1:order-1,0:nx,-span_interp:nz+span_interp))
 allocate(DERxr(1:order-1,0:nx),DERzs(1:order-1,0:nz))
 allocate(rtr0(1:ntr),rtr1(1:ntr))
 allocate(str0(1:ntr),str1(1:ntr))
 if (torder.gt.2) allocate(rtr2(1:ntr),str2(1:ntr))
 if (local_time.eq.1) allocate(dt_tr(1:ntr))
end if
 allocate(conduct(-spanT:nx+spanT,-spanT:nz+spanT))
 allocate(conduct_r(-span_interp:nx+span_interp,-span_interp:nz+span_interp))
 allocate(conduct_s(-span_interp:nx+span_interp,-span_interp:nz+span_interp))
 allocate(strain(-span:nx+span,-span:nz+span))
if (iterate.eq.1) then
 allocate(dt_vis(1:grids,0:nx,0:nz))
 allocate(dt_stream(1:grids,1:nx-1,1:nz-1),residual(1:nx-1,1:nz-1))
 allocate(vis_x(1:grids,1:nx-1,1:nz-1),vis_z(1:grids,1:nx-1,1:nz-1),vis_xx(1:grids,1:nx-1,1:nz-1))
 allocate(vis_zz(1:grids,1:nx-1,1:nz-1),vis_xz(1:grids,1:nx-1,1:nz-1))
 allocate(vis_xf(1:grids,1:nx-1,1:nz-1),vis_zf(1:grids,1:nx-1,1:nz-1),vis_xxf(1:grids,1:nx-1,1:nz-1))
 allocate(vis_zzf(1:grids,1:nx-1,1:nz-1),vis_xzf(1:grids,1:nx-1,1:nz-1))
 allocate(RHS(1:grids,1:nx-1,1:nz-1),error(1:grids,-span:nx+span,-span:nz+span))
 allocate(vis_grid(1:grids,-span:nx+span,-span:nz+span))
 allocate(vis_gridf(1:grids,-span:nx+span,-span:nz+span))
 allocate(ratio(1:grids-1,1:nx-1,1:nz-1))
end if
if (iterate.eq.0) then
 grids=1
 allocate(vis_grid(1:grids,-span:nx+span,-span:nz+span))
 allocate(dt_vis(1:grids,0:nx,0:nz))
 allocate(vis_x(1:grids,1:nx-1,1:nz-1),vis_z(1:grids,1:nx-1,1:nz-1),vis_xx(1:grids,1:nx-1,1:nz-1))
 allocate(vis_zz(1:grids,1:nx-1,1:nz-1),vis_xz(1:grids,1:nx-1,1:nz-1))
 allocate(ratio(1:grids,1:nx-1,1:nz-1))
end if
if (equalize.eq.1) then
 allocate(Cinit(1:ntype))
end if
end

subroutine compute_basics
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: j,Nthreads,ran_int
Dnumber=order  !!numder of derivatives needed in tracer velocity interpolator (plus one)
span1=order/2   !!stencil span for 1st and 2nd derivatives
span=order/2+1  !!full stencil span (due to 4th derivatives)
span_interp=(order+Dnumber-2)/2 !!for interpolation to tracer positions
spanT=span_interp+span1
ngrid=(nx+1)*(nz+1)
dr=aspect/real(nx,8)
ds=1.d0/real(nz,8)
dr2=dr**2.d0
ds2=ds**2.d0
dr3=dr**3.d0
ds3=ds**3.d0
dr4=dr**4.d0
ds4=ds**4.d0
drds=dr*ds
drds2=dr*ds2
dr2ds=dr2*ds
dr2ds2=dr2*ds2
dro2=0.5d0*dr
dso2=0.5d0*ds
dro10=0.1d0*dr
dso10=0.1d0*ds
Ax=(datan(spaceX*(0.5d0))-datan(spaceX*(-0.5d0)))**(-1.d0)
Az=(datan(spaceZ*(0.5d0))-datan(spaceZ*(-0.5d0)))**(-1.d0)
 allocate(factorial(0:Dnumber),dr_power(0:Dnumber),ds_power(0:Dnumber))
 call compute_factorials
write(*,*) "nx,nz=",nx,nz
if ((comp.eq.1).or.(tracer.eq.1)) then
 ntr=(nx+1)*(nz+1)*tpc !!# of total tracers
 write(*,*) "ntr=",ntr
end if
tmultigrid=0.d0
tenergy=0.d0
tcd=0.d0
tbSF=0.d0
tsSF=0.d0
tequalize=0.d0
tadvect_tracer=0.d0
nequalize=0
tconvert=0.d0
tcoarse_matrix=0.0d0
dist_r=dist_factor*dr
dist_s=dist_factor*ds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (iterate.eq.0) grids=1
 allocate(nx_grid(1:grids),nz_grid(1:grids))
 allocate(dr_grid(1:grids),ds_grid(1:grids),dr2_grid(1:grids),ds2_grid(1:grids),dr3_grid(1:grids),ds3_grid(1:grids))
 allocate(dr4_grid(1:grids),ds4_grid(1:grids),drds_grid(1:grids),drds2_grid(1:grids),dr2ds_grid(1:grids))
 allocate(dr2ds2_grid(1:grids))
 do j=1,grids
  nx_grid(j)=nx/2**(j-1)
  nz_grid(j)=nz/2**(j-1)
  dr_grid(j)=aspect/real(nx_grid(j),8)
  ds_grid(j)=1.d0/real(nz_grid(j),8)
  dr2_grid(j)=dr_grid(j)**2.d0
  ds2_grid(j)=ds_grid(j)**2.d0
  dr3_grid(j)=dr_grid(j)**3.d0
  ds3_grid(j)=ds_grid(j)**3.d0
  dr4_grid(j)=dr_grid(j)**4.d0
  ds4_grid(j)=ds_grid(j)**4.d0
  drds_grid(j)=dr_grid(j)*ds_grid(j)
  drds2_grid(j)=dr_grid(j)*ds2_grid(j)
  dr2ds_grid(j)=dr2_grid(j)*ds_grid(j)
  dr2ds2_grid(j)=dr2_grid(j)*ds2_grid(j)
 end do
 ngrid_coarse=(nx_grid(grids)+1)*(nz_grid(grids)+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (iterate.eq.1) then
 ku=0  !!for banded matrix storage and routines
 kl=0
 analysis=1
 call Compute_Matrix_coarse !!figure out bandwidth before allocating
 analysis=0
 allocate(Matrix_coarse(1:2*kl+ku+1,1:ngrid_coarse),B_coarse(1:ngrid_coarse),IPIV(1:ngrid_coarse))
 allocate(RowB(1:ngrid_coarse),ColB(1:ngrid_coarse))
 allocate(res_ref(1:grids))
 build=0
 scaling=0
end if
smooth_vis=1
!!!!!!!!!!!!!!!!!! Initial random seeds for mzran for each OpenMP thread
Nthreads=1 !!serial case
!$ Nthreads=omp_get_max_threads()
allocate(seed_i(0:Nthreads-1),seed_j(0:Nthreads-1),seed_k(0:Nthreads-1),seed_n(0:Nthreads-1))
 seed_i(0)=521288629
 seed_j(0)=362436069
 seed_k(0)=16163801
 seed_n(0)=1131199299
do j=1,Nthreads-1
 ran_int=seed_i(j-1)-seed_k(j-1)
 if (ran_int.lt.0) ran_int=ran_int+2147483579
 seed_i(j)=seed_j(j-1)
 seed_j(j)=seed_k(j-1)
 seed_k(j)=ran_int
 seed_n(j)=69069*seed_n(j-1)+1013904243
end do
end

subroutine initial_temperature
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: Tu,Tl,Tr,Ts,Q,u0,v0,x,z,temperature
 do i=0,nx
 x=xg(i)
  do k=0,nz
   z=zg(k)
   T(i,k)=temperature(x,z)
  end do
 end do

T(:,0)=1.d0  !!BCs
T(:,nz)=0.d0
 call enforceBCs(T)
 if (tracer.eq.0) call smoother_time_T(T)
end

subroutine stats
 use basics
 use arrays
implicit none
integer*4 :: i,k,empty_cells,tracer_min,tracer_max,tt
real*8 :: Tz_surf(0:nx),Tz_bot(0:nx),ent(1:ntype),zgrid
real*8 :: integrate1Dx,integrate2D,Tavg,vrms,Nu,Cavg(1:ntype),entrainment,Nu_bot,RaT_avg
real*8 :: q1,q2,q3,q4,strain_max,mobility
 if (iterate.eq.1) vis(:,:)=vis_grid(1,:,:)
 call compute_strain_rate_invariant
 strain_max=maxval(strain(0:nx,nz))
if (comp.eq.1) then
 do tt=1,ntype
  Cavg(tt)=integrate2D(Cnew(tt,0:nx,0:nz))/aspect
 end do
 if (tstep.eq.0) then
  call entrainment_setup
  zent=zgrid(real(kent,8)*ds)
  write(*,*) "Entrainment Calculation:",kent,zent
  do tt=1,ntype
   mass(tt)=Cavg(tt)
  end do
 end if
 do tt=1,ntype
  ent(tt)=entrainment(Cnew(tt,0:nx,0:nz),tt)
 end do
end if
if ((comp.eq.1).or.(tracer.eq.1)) then
 empty_cells=count(nempty(0:nx,0:nz).eq.0)
 tracer_min=minval(nempty(0:nx,0:nz))
 tracer_max=maxval(nempty(0:nx,0:nz))
end if
vrms=(integrate2D(u(0:nx,0:nz)**2.d0+w(0:nx,0:nz)**2.d0)/(aspect))**0.5d0
mobility=((integrate1Dx(u(0:nx,nz)**2.d0+w(0:nx,nz)**2.d0)/aspect)**0.5d0)/vrms
if (RaT.ne.0.d0) then
 Tavg=integrate2D(T(0:nx,0:nz))/aspect  !!don't include ghost points
 RaT_avg=integrate2D(RaT/conduct(0:nx,0:nz)/vis(0:nx,0:nz))/aspect !!!!!!!!!!!
 do i=0,nx
  Tz_surf(i)=dot_product(D1(-span1:span1),T(i,nz-span1:nz+span1))/ds/zs(nz)
  Tz_bot(i)=dot_product(D1(-span1:span1),T(i,-span1:span1))/ds/zs(0)
 end do
 Nu=integrate1Dx(-conduct(0:nx,nz)*Tz_surf(0:nx))/integrate1Dx(T(0:nx,0))
 Nu_bot=integrate1Dx(-conduct(0:nx,0)*Tz_bot(0:nx))/aspect
 q1=-conduct(0,nz)*dot_product(D1(-span1:span1),T(0,nz-span1:nz+span1))/ds/zs(nz)
 q2=-conduct(nx,nz)*dot_product(D1(-span1:span1),T(nx,nz-span1:nz+span1))/ds/zs(nz)
 q3=-conduct(nx,0)*dot_product(D1(-span1:span1),T(nx,-span1:span1))/ds/zs(0)
 q4=-conduct(0,0)*dot_product(D1(-span1:span1),T(0,-span1:span1))/ds/zs(0)
 if (local_time.eq.0) then
  if (comp.eq.0) then
   if (tracer.eq.0) then
    write(303,'(12(g20.8))') time,Tavg,vrms,Nu,Nu_bot,RaT_avg,mobility,strain_max,q1,q2,q3,q4
   elseif (tracer.eq.1) then
    write(303,'(12(g20.8),3(i6))') time,Tavg,vrms,Nu,Nu_bot,RaT_avg,mobility,strain_max,q1,q2,q3,q4,&
&empty_cells,tracer_min,tracer_max
   end if
  else
   write(303,'(12(g20.8),3(i6),10(g20.8))') time,Tavg,vrms,Nu,Nu_bot,RaT_avg,mobility,strain_max,q1,q2,q3,q4,&
&empty_cells,tracer_min,tracer_max,Cavg(1:ntype),ent(1:ntype)
  end if
 else
  if (tracer.eq.0) then
   write(303,'(i10,11(g20.8))') tstep,Tavg,vrms,Nu,Nu_bot,RaT_avg,mobility,strain_max,q1,q2,q3,q4
  else
   write(303,'(i10,11(g20.8),3(i6))') tstep,Tavg,vrms,Nu,Nu_bot,RaT_avg,mobility,strain_max,q1,q2,q3,q4,&
&empty_cells,tracer_min,tracer_max
  end if
 end if
else
 if (local_time.eq.0) then
  if (comp.eq.0) then
   write(303,'(4(g20.8))') time,vrms,mobility,strain_max
  else
   write(303,'(4(g20.8),3(i6),10(g20.8))') time,vrms,mobility,strain_max,empty_cells,tracer_min,tracer_max,&
&Cavg(:),ent(:)
  end if
 else
  write(303,'(i10,3(g20.8))') tstep,vrms,mobility,strain_max
 end if
end if
!!!!!!!!!!!!!!!!!!viscosity smoother info
   if ((RaT.ne.0.d0).and.(comp.eq.0)) then
    write(304,'(g20.8,i5,g20.8,i5)') time,n_eta,eig_fine,n_temp
   elseif ((RaT.ne.0.d0).and.(comp.eq.1)) then
    write(304,'(g20.8,i5,g20.8,11(i5))') time,n_eta,eig_fine,n_temp,n_comp(1:ntype)
   elseif ((RaT.eq.0.d0).and.(comp.eq.1)) then
    write(304,'(g20.8,i5,g20.8,10(i5))') time,n_eta,eig_fine,n_comp(1:ntype)
   end if
if ((time.eq.0.d0).and.(equalize.eq.1)) Cinit=Cavg 
end

real*8 function entrainment(f,tt)
 use basics
 use arrays
implicit none
integer*4 :: i,k,kmin,kmax,tt
real*8 :: integrate1Dx,integrate1Dz_ent
real*8 :: fintx(0:nz),fint
real*8 :: f(0:nx,0:nz),mass0
do k=0,nz
 fintx(k)=integrate1Dx(f(0:nx,k))
end do
fint=integrate1Dz_ent(fintx(0:nz))
entrainment=fint/aspect/mass(tt)
end

real*8 function integrate2D(f)
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: integrate1Dx,integrate1Dz
real*8 :: fintx(0:nz),fint
real*8 :: f(0:nx,0:nz)
do k=0,nz
 fintx(k)=integrate1Dx(f(0:nx,k))
end do
fint=integrate1Dz(fintx(0:nz))
integrate2D=fint
end

real*8 function integrate1Dx(f)
 use basics
 use arrays
implicit none
integer*4 :: i,j
real*8 :: f(0:nx),fint,hx,derivative
!!!!original
fint=0.d0
do i=0,nx-1
 fint=fint+0.5d0*(f(i)*xr(i)+f(i+1)*xr(i+1))*dr
end do
!!!!end original

!hx=dx/2.d0
!fint=f(i)*h
!do j=2,stats_order,2
! derivative=dot_product(D(j,-span_stats,span_stats),f(i-span_stats,i+span_stats))/dr**real(j,8)
! fint=fint+2.d0*derivative*h**real(j+1)/factorial(j+1)
!end do

integrate1Dx=fint
end

real*8 function integrate1Dz(f)
 use basics
 use arrays
implicit none
integer*4 :: k
real*8 :: f(0:nz),fint
fint=0.d0
do k=0,nz-1
 fint=fint+0.5d0*(f(k)*zs(k)+f(k+1)*zs(k+1))*ds
end do
integrate1Dz=fint
end

real*8 function integrate1Dz_ent(f)
 use basics
 use arrays
implicit none
integer*4 :: k
real*8 :: f(0:nz),fint
fint=0.d0
do k=kent,nz-1
 fint=fint+0.5d0*(f(k)*zs(k)+f(k+1)*zs(k+1))*ds
end do
integrate1Dz_ent=fint
end

subroutine compute_factorials
 use basics
 use arrays
implicit none
integer*4 :: j
factorial(0)=1.d0
dr_power(0)=1.d0
ds_power(0)=1.d0
do j=1,Dnumber
 factorial(j)=real(j,8)*factorial(j-1)
 dr_power(j)=dr**real(j,8)
 ds_power(j)=ds**real(j,8)
end do
end
