program hindenburg
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: t1,t2,tmin
integer*4 :: i,j,k,N,level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BEGIN INPUTS
!!!!!!!!!!!!! grid options
aspect=1.d0       !!aspect ratio of Cartesian domain
nx=200            !!# of cells in the horizontal
nz=200            !!# of cells in the vertical
Database=2        !!ordering of streamfunction and buoyancy vectors -- 1 (large nz), 2 (large nx), or 3 (experiment)
order=2           !!space order: even numbers
gridX=0           !!0=uniform,1=arctan
gridZ=0           !!0=uniform,1=arctan,2=boundary layer
thickness=0.1d0   !!gridZ=2: thickness of grid boundary layer
nbound=8          !!# of cell corners in boundary layer grid
spaceX=1.d0       !!stretching factor for arctan function
spaceZ=1.d0       !!stretching factor for arctan and boundary layer grids
!!!!!!!!!!!!! end grid options

!!!!!!!!!!!!!velocity solver options
iterate=1         !!how to solve for the velocities?: 0=direct method, 1=iterative method
non_Newtonian=1   !!is the fluid non-Newtonian? (only works if iterate=1)
iterate_order=1   !!order of iterative smoother: 1=RK1 (done), 4=RK4 (work in progress)
tolerance=1.d-3   !!tolerance for multigrid -- smaller values indicate more accurate solutions
Nmulti=70         !!# of iterations on each multigrid level
mingrid=4         !!min grid size (<nx and <nz)
grids=2           !!max # of multigrid levels (currently between 2 and 6)
eig_ratio=4.0d0   !!min ratio of diffusive to advective eigenvalues desired
vis_iter=12       !!max # of viscosity smoother iterations for multigrid
max_cycles=1000   !!max # of multigrid cycles allowed for convergence before code stops
!!!!!!!!!!!!!end velocity solver options

!!!!!!!!!!!!! time stepping options
torder=2              !!order of RK method: RK2 only if temperature tracers are used - RK4 for field method and/or tracers (2 and 4 currently available)
tend=1.d-4          !!if local time steps are NOT used, run controlled by time cutoff 
tsnap=tend/5.d0      !!if local time steps are NOT used, snapshots are taken at time intervals
tmin=tend/2.d0        !!if local_time=0, begin time averaging at this time for stats
N=1000                !!if local time steps ARE used, run controlled by iteration cutoff
Nsnap=100             !!if local time steps ARE used, snapshots controlled by iteration count
local_time=0          !!only use with RK4 without updates (for field methods primarily)
dt_manual=0.5d0       !!manual time step size for isothermal convection
courant=0.99d0        !!factor that multiplies the time step size computed from stability analysis
courantT=0.99d0       !!for temperature smoother
courantC=0.99d0       !!for composition smoother
courant_stream=0.99d0 !!for iterative solver
RK_update=1           !!update the velocities between RK4 stages? (use 0 for steady state problems, 1 for time accurate runs)
!!!!!!!!!!!!! end time stepping options

!!!!!tracer options
comp=0                  !!switch for composition: 0=off,1=on
ntype=0                 !!# of tracer types:=0 for uniform composition with T tracers, >=1 for heterogenious composition
equalize=1              !!redistribute tracers?
tracer=1                !!temperature tracers?
tpc=30                  !!# of tracers per cell corner
dist_factor=1.d-4       !! <<1 controls min. dist between tracers and side boundaries during equalize
!!!!!end tracer options

!!!!!!!!!!!!!!!!!!!!Initial and Boundary Conditions
initial_temp=-1         !!initial T: -1=user defined, 0 for conductive, 1 for convective cell, 2 alternative conductive
initial_comp=0          !!initial C: -1=user defined, 0=dense layer at base, 1=sinusoidal deflection, 2=step function
initial_temp_smooth=1   !!smooth the initial temperature field? (controlled by delta_vis in the viscosity options)
delta_visT_init=2.0d0    !!controls smoothing of the initial T field based upon viscosity dependence on T
restart=1               !!if 0 use analytic initial T and C, if 1 use T (and possibly C -- see "add_layer" tracer option) from previous run
add_layer=0             !!add heterogenious composition to a previous isocompositional run when restarting?
dlayer=0.1d0            !!for initial_comp=0 or 1: layer thickness
dent=0.1d0              !!reference height for entrainment calculation

!! reflecting sidewalls assumed
Tbot=0 !!Bottom thermal BC: 0=isothermal, 1=insulating
Vbc=0  !!Top and Bottom velocity BCs: 0=free-slip,1=rigid
!!!!!!!!!!!!!!!!!!!!End Initial and Boundary Conditions

!!!!!!!!!!!!!Rayleigh Numbers and Buoyancy
RaT=(1.d3)/3.d0   !!Thermal Rayleigh number
if (comp.eq.1) then
 allocate(RaC(1:ntype))
 RaC(1)=0.d0           !!Compositional Rayleigh Number(s): >0 for increased density, <0 for decreased density
! RaC(2)=-1.d4         !!create as many terms as there are enriched compositions
! RaC(3)=-3.d5
end if
!!!!!!!!!!!!!End Rayleigh Numbers and Buoyancy

!!!!!!!!!!Internal heating and Thermal Conductivity Options
H=0.d0                  !!internal heating rate of ambient material
if (comp.eq.1) then
 allocate(conduct_factor(1:ntype),Htype(1:ntype))
 conduct_factor(1)=1.d0  !!multiplier for nondimesnional conductivity within each compositional component (C is smoothed by delta_visC)
! conduct_factor(2)=1.d0 !!create as many terms as there are enriched compositions

 Htype(1)=0.d0           !!set H value within components
! Htype(2)=0.d0        !!create as many terms as there are enriched compositions
end if
!!!!!!!!!!!Internal heating and Thermal Conductivity Options

!!!!!!!!!viscosity options
eig_ratio_min=4.d0             !!min ratio of diffusive to advective eigenvalues in streamfunction eqn
ivis=3                         !!viscosity: 0=isoviscous, 1=T and z dependent, 2=dependence on z, T and C, 3=T dependence and plastic yielding
visP=1.d0                      !!for depth dependence (>1 for an increase with depth)
visT=3.d4                      !!for temperature dependence (>1 for decrease)
delta_visT=100.d0               !!max allowed viscosity contrast due to T between adjacent grid points (T used to compute viscosity field is smoothed until this criterion is satisfied)
yield=1.5d1                     !!yield stress at the surface
yield_gradient=0.d0            !!yield stress gradient (linear increase with depth)
yield_depth=0.1d0                !!max depth that yielding can occur
if (comp.eq.1) then            !!composition dependence (>1 for decrease, <1 for increase)
 allocate(visC(1:ntype),delta_visC(1:ntype))
 visC(1)=1.d0 
! visC(2)=1.d-2  !!create as many terms as there are enriched compositions 
! visC(3)=1.d1
 delta_visC(1)=2.0d0            !!max allowed viscosity/conductivity contrast due to C between adjacent grid points (multiplicative factor)
! delta_visC(2)=2.d0
! delta_visC(3)=1.d0
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END INPUTS

 if (equalize.eq.1) open(unit=1024,file="equalize.dat")
 call compute_smoothing_parameters
 if (iterate.eq.1) call multigrid_setup
 call compute_basics
 call storage
 if (Database.eq.3) call Database3
 call FD_coefficients
 call build_grid
 call derivative_eigenvalues
 call stability_RK4_smoother
 if (iterate.eq.1) then
  do level=1,grids
   call stability_vis_smoother(level)
  end do
 else
   grids=1
   level=1
   call stability_vis_smoother(level)
 end if
 if (restart.eq.1) then
  call smoother_time_T(T)
  if (comp.ne.0) then
   call smoother_time
  end if
 call load_conductivity_array
 end if
 if (restart.eq.0) call initial_temperature
 if ((restart.eq.0).and.((comp.eq.1).or.(tracer.eq.1))) then
  call initialize_tracers
 end if
!!!!!Initial temperature smoothing here
 if ((restart.eq.0).and.(initial_temp_smooth.eq.1)) call smoother_initial_T
 if (((restart.eq.1).and.(add_layer.eq.1)).and.(comp.eq.1)) then
  if (tracer.eq.0) then
   call initialize_tracers
  else
   call add_tracers
  end if
 end if
 tstep=-1
 if (iterate.eq.0) then
  call build_SF_matrix !!find bandwith and banded matrix size
 elseif (iterate.eq.1) then
  if (restart.eq.0) SF=0.d0 !!initial guess for stream function
  error(2:grids,:,:)=0.d0
 end if
! call tracer_snapshots !!!!time sink for large jobs
 call cpu_time(t1)
!$ t1=omp_get_wtime()
 open(unit=303,file="stats.dat")
 open(unit=304,file="viscosity_smoother.dat")
fcount=0
Nout=Nsnap
tout=0.d0
time=0.d0
if (local_time.eq.0) N=1000000 !!!max iterations for time accurate run
do tstep=0,N
 if (iterate.eq.0) then
  if ((ivis.ge.1).or.(tstep.eq.0)) then
   call build_SF_matrix  
  end if
  call solve_SF_equation
 elseif (iterate.eq.1) then
  call multigrid
 end if
 call compute_velocities
 call stats
 if (local_time.eq.0) then  !!snapshots controlled by time
  if (time.ge.tout) then
   call snapshots
   tout=tout+tsnap
   fcount=fcount+1
   if (time.ge.tend) exit
  end if
 else                       !!snapshots controlled by iteration
  if (Nout.ge.Nsnap) then
   call snapshots
   Nout=1
   fcount=fcount+1
  else
   Nout=Nout+1
  end if
 end if
 call energy_time
 if (local_time.eq.0) then
  time=time+dt
 end if
end do
 close(303)
 close(304)
 call cpu_time(t2)
!$ t2=omp_get_wtime()
 write(*,*) "Run Time=",t2-t1,"s"
 write(*,*) "energy_time=",tenergy,"s"
 if ((tracer.eq.1).or.(comp.eq.1)) write(*,*) "advect_tracer=",tadvect_tracer,"s"
 if (iterate.eq.0) then
  write(*,*) "build_SF_matrix=",tbSF,"s"
  write(*,*) "solve_SF_equation=",tsSF,"s"
 elseif (iterate.eq.1) then
  write(*,*) "Multigrid=",tmultigrid,"s"
  write(*,*) "Compute_Matrix_coarse=",tcoarse_matrix,"s"
 end if
 if ((tracer.eq.1).or.(comp.eq.1)) then
  write(*,*) "compute_derivatives=",tcd,"s"
  write(*,*) "tracers_to_corners=",tconvert,"s"
 end if
 if (equalize.eq.1) then
  write(*,*) "equalize=",tequalize,"s"
 end if
  write(*,*) "time steps=",tstep
 call print_restart
! call tracer_snapshots
 if ((local_time.eq.0).and.(RaT.ne.0.d0).and.(tstep.gt.10)) call time_average(tmin,tend)
end program hindenburg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine snapshots
 use basics
 use arrays
implicit none
real*8 :: xgrid,zgrid
integer*4 :: i,k,id
 character*6 :: fname
write(fname,'(a1,i5)') 'T',10000+fcount
open(unit=101,file=fname)
if (RaT.ne.0.d0) then
 do i=0,nx
  do k=0,nz
   if ((comp.eq.0).and.(tracer.eq.0)) then
    write(101,'(8(g14.5E3))') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),vis(i,k),Tvis(i,k)
   elseif ((comp.eq.0).and.(tracer.eq.1)) then
    write(101,'(9(g14.5E3),i6)') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),strain(i,k),vis(i,k),Tvis(i,k),&
&nempty(i,k)
   else
 write(101,'(9(g14.5E3),i6,20(g14.5E3))') xg(i),zg(k),T(i,k),SF(i,k),u(i,k),w(i,k),strain(i,k),vis(i,k),Tvis(i,k),nempty(i,k),&
&Cnew(1:ntype,i,k),Cvis(1:ntype,i,k)
   end if
  end do
  write(101,*) " "
 end do
else                !!!!isothermal convection (pure compositional convection)
 do i=0,nx
  do k=0,nz
 write(101,'(7(g14.5E3),i6,20(g14.5E3))') xg(i),zg(k),SF(i,k),u(i,k),w(i,k),strain(i,k),vis(i,k),nempty(i,k),&
&Cnew(1:ntype,i,k),Cvis(1:ntype,i,k)
  end do
  write(101,*) " "
 end do
end if
 close(101)
end

subroutine tracer_snapshots
 use basics
 use arrays
implicit none
real*8 :: xgrid,zgrid
integer*4 :: id
if ((comp.eq.1).or.(tracer.eq.1)) then
 if (tstep.eq.-1) then
  open(unit=666,file="tracers.dat")
 else
  open(unit=666,file="tracers.dat",position='append')
 end if
 if ((comp.eq.1).and.(tracer.eq.1)) then
    do id=1,ntr
     write(666,'(3(f9.5),i4)') xgrid(rtr(id)),zgrid(str(id)),Ttr(id),Ttype(id)
    end do
 elseif (comp.eq.1) then
    do id=1,ntr
     write(666,'(2(f9.5),i4)') xgrid(rtr(id)),zgrid(str(id)),Ttype(id)
    end do
 elseif (tracer.eq.1) then
    do id=1,ntr
     write(666,'(3(f9.5))') xgrid(rtr(id)),zgrid(str(id)),Ttr(id)
    end do
 end if
    write(666,*) " "
    write(666,*) " "
    close(666)
end if
end

subroutine print_restart
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID
open(unit=101,file="Trestart",form='unformatted')
if (RaT.ne.0.d0) then
 if (comp.eq.0) then
  write(101) xg,zg,T,SF
 else
  write(101) xg,zg,T,Cnew,SF
 end if
else
 if (comp.eq.0) then
  write(101) xg,zg,SF
 else
  write(101) xg,zg,Cnew,SF
 end if
end if
 close(101)
if ((comp.eq.1).or.(tracer.eq.1)) then
 open(unit=101,file="restart_tracers",form='unformatted')
  if (tracer.eq.0) then
   write(101) rtr,str,Ttype
  elseif (comp.eq.0) then
   write(101) rtr,str,Ttr
  else
   write(101) rtr,str,Ttype,Ttr
  end if
 close(101)
end if
end

subroutine test_velocities
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: x,z
do i=0,nx
 x=xg(i)
 do k=0,nz
 z=zg(k)
  u(i,k)=-V*dsin(pii*x/aspect)*dcos(pii*z)
  w(i,k)=(V/aspect)*dcos(pii*x/aspect)*dsin(pii*z)
 end do
end do
end

subroutine FD_coefficients
 use basics
 use arrays
implicit none
integer*4 :: i,m,n,bigN,bigM,k
real*8 :: x0,c1,c2,c3
real*8, allocatable :: a(:),delta(:,:,:)
 bigM=max(4,Dnumber)  !!need up to 4th derivatives in SF equation - may need more for tracer velocity interpolation
 bigN=order+bigM-2 !!assuming bigM is even
! bigN=order+2
 allocate(a(0:bigN),delta(-1:bigM,0:bigN,0:bigN))
 delta=0.d0
 x0=0.d0
 a(0)=0.d0
 k=1
 do i=1,bigN,2
  a(i)=real(k,8)
  a(i+1)=-real(k,8)
  k=k+1
 end do

 delta(0,0,0)=1.d0
 c1=1.d0
 do n=1,bigN
  c2=1.d0
  do i=0,n-1
   c3=a(n)-a(i)
   c2=c2*c3
   if (n.le.bigM) delta(n,n-1,i)=0.d0
   do m=0,min(n,bigM)
    delta(m,n,i)=((a(n)-x0)*delta(m,n-1,i)-real(m,8)*delta(m-1,n-1,i))/c3
   end do
  end do
  do m=0,min(n,bigM)
   delta(m,n,n)=(c1/c2)*(real(m,8)*delta(m-1,n-1,n-1)-(a(n-1)-x0)*delta(m,n-1,n-1))
  end do
 c1=c2
 end do

 D1=0.d0
 D2=0.d0
 D1(0)=delta(1,order,0)
 D2(0)=delta(2,order,0)
 k=1
 do i=1,order,2
  D1(k)=delta(1,order,i)
  D1(-k)=delta(1,order,i+1)
  D2(k)=delta(2,order,i)
  D2(-k)=delta(2,order,i+1)
 k=k+1
 end do

!!!!!!!!for stats integrals
 D(1:Dnumber,-span_interp:span_interp)=0.d0
 do m=2,Dnumber,2
  D(m-1,0)=delta(m-1,order+m-2,0)
  D(m,0)=delta(m,order+m-2,0)
  k=1
  do i=1,order+m-2,2
   D(m-1,k)=delta(m-1,order+m-2,i)
   D(m-1,-k)=delta(m-1,order+m-2,i+1)
   D(m,k)=delta(m,order+m-2,i)
   D(m,-k)=delta(m,order+m-2,i+1)
  k=k+1
  end do
 end do
 do m=1,Dnumber,2
  D(m,0)=0.d0     !!!enforce symemetry and antisymmetry properties
  do i=1,span_interp
   D(m,i)=-D(m,-i)
   D(m+1,i)=D(m+1,-i)
  end do
 end do
!!!!!!!!

 D3=0.d0
 D4=0.d0
 D3(0)=delta(3,order+2,0)
 D4(0)=delta(4,order+2,0)
 k=1
 do i=1,order+2,2
  D3(k)=delta(3,order+2,i)
  D3(-k)=delta(3,order+2,i+1)
  D4(k)=delta(4,order+2,i)
  D4(-k)=delta(4,order+2,i+1)
 k=k+1
 end do
 

 D1(0)=0.d0     !!!enforce symemetry and antisymmetry properties
 D3(0)=0.d0
 do i=1,span
  D1(i)=-D1(-i)
  D2(i)=D2(-i)
  D3(i)=-D3(-i)
  D4(i)=D4(-i)
 end do

 Drs=0.d0
 Drrss=0.d0
 Drss=0.d0
 Drrs=0.d0
 do k=-span1,span1
  do i=-span1,span1
   Drs(i,k)=D1(k)*D1(i)
   Drss(i,k)=D2(k)*D1(i)
   Drrs(i,k)=D1(k)*D2(i)
   Drrss(i,k)=D2(k)*D2(i)
  end do
 end do

! write(*,'(100(g12.3,x))') D1
! write(*,'(100(g12.3,x))') D2
! write(*,'(100(g12.3,x))') D3
! write(*,'(100(g12.3,x))') D4
! write(*,*) ""
! do m=1,Dnumber
!  write(*,'(100(g12.3,x))') D(m,-span_interp:span_interp)
! end do
end
