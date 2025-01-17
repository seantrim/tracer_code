real*8 function heating(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,j,k
if (comp.eq.0) then
 heating=H
else
  heating=H*(1.d0-sum(Cvis(1:ntype,i,k)))  !!contribution from ambient material
  do j=1,ntype
   heating=Htype(j)*Cvis(j,i,k)+heating
  end do
end if
end

real*8 function conductivity(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,j,k

if (comp.eq.0) then
 conductivity=1.d0
else
 do j=1,ntype
  conductivity=dexp(-Cvis(j,i,k)*dlog(1.d0/conduct_factor(j)))
 end do
end if
end

subroutine load_conductivity_array
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: conductivity
 do i=0,nx
  do k=0,nz
   conduct(i,k)=conductivity(i,k)
  end do
 end do
 do i=1,spanT
  conduct(-i,0:nz)=conduct(i,0:nz)
  conduct(nx+i,0:nz)=conduct(nx-i,0:nz)
 end do
 do k=1,spanT
  conduct(-spanT:nx+spanT,-k)=conduct(-spanT:nx+spanT,k)
  conduct(-spanT:nx+spanT,nz+k)=conduct(-spanT:nx+spanT,nz-k)
 end do

 do i=-span_interp,nx+span_interp
  do k=-span_interp,nz+span_interp
   conduct_r(i,k)=dot_product(D1(-span1:span1),conduct(i-span1:i+span1,k))/dr
   conduct_s(i,k)=dot_product(D1(-span1:span1),conduct(i,k-span1:k+span1))/ds
  end do
 end do

end

real*8 function energy_space(i,k,T0)
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz,kx,kz
real*8 :: T0(-span:nx+span,-span:nz+span),heating
Tr=dot_product(D1(-span1:span1),T0(i-span1:i+span1,k))/dr
Ts=dot_product(D1(-span1:span1),T0(i,k-span1:k+span1))/ds
Trr=dot_product(D2(-span1:span1),T0(i-span1:i+span1,k))/dr2
Tss=dot_product(D2(-span1:span1),T0(i,k-span1:k+span1))/ds2

kx=conduct_r(i,k)/xr(i)
kz=conduct_s(i,k)/zs(k)
Tx=Tr/xr(i)
Tz=Ts/zs(k)
Txx=(Trr-Tx*xrr(i))/xr(i)**2.d0
Tzz=(Tss-Tz*zss(k))/zs(k)**2.d0
energy_space=conduct(i,k)*(Txx+Tzz)+(kx*Tx+kz*Tz)-u(i,k)*Tx-w(i,k)*Tz+heating(i,k)
end

real*8 function tracer_space(i,k) !!Lagrangian frame
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: Tr,Ts,Trr,Tss
real*8 :: Tx,Tz,Txx,Tzz,kx,kz,heating
Tr=dot_product(D1(-span1:span1),Textend(i-span1:i+span1,k))/dr
Ts=dot_product(D1(-span1:span1),Textend(i,k-span1:k+span1))/ds
Trr=dot_product(D2(-span1:span1),Textend(i-span1:i+span1,k))/dr2
Tss=dot_product(D2(-span1:span1),Textend(i,k-span1:k+span1))/ds2

kx=conduct_r(i,k)/xr(i)
kz=conduct_s(i,k)/zs(k)
Tx=Tr/xr(i)
Tz=Ts/zs(k)
Txx=(Trr-Tx*xrr(i))/xr(i)**2.d0
Tzz=(Tss-Tz*zss(k))/zs(k)**2.d0
tracer_space=conduct(i,k)*(Txx+Tzz)+(kx*Tx+kz*Tz)+heating(i,k)
end

subroutine load_tracer_space_array
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: tracer_space
  call extendT(Tratio)
 do i=-span_interp,nx+span_interp
  do k=-span_interp,nz+span_interp
   tracer_space_array(i,k)=tracer_space(i,k)
  end do
 end do
end

subroutine energy_time
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: i,k,ID
real*8 :: energy_space,temp,tempr,temps,tempT,xr_interp,zs_interp
real*8 :: T0(-span:nx+span,-span:nz+span),T1(-span:nx+span,-span:nz+span),T2(-span:nx+span,-span:nz+span)
real*8, allocatable :: C0(:,:,:)
real*8 :: time1,time2,dto6,ta1,ta2,dtt,dtstage
!$ time1=omp_get_wtime()
if (local_time.eq.0) then
 call stability_RK4
elseif(local_time.eq.1) then
 call local_time_steps
end if
if (equalize.eq.1) then
 call equalize_tracers
end if
if (comp.eq.1) then
 allocate(C0(1:ntype,-span:nx+span,-span:nz+span))
 C0=Cnew
end if
if (tracer.eq.0) then
 T0=T
elseif (tracer.eq.1) then
 T0=Tratio
 Ttr0=Ttr
end if
if ((comp.eq.1).or.(tracer.eq.1)) then
 rtr0=rtr
 str0=str
end if

if (torder.ge.1) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 1st stage
 if (RaT.ne.0.d0) then
  if (tracer.eq.1) then
   call load_tracer_space_array
   call compute_derivatives(tracer_space_array,DERr_tracer)
  end if
  if (tracer.eq.0) then
  do i=0,nx
   do k=0,nz
    temp=energy_space(i,k,T0)
    T(i,k)=temp
    if (local_time.eq.1) dt=dt_array(i,k)
    T1(i,k)=T0(i,k)+0.5d0*dt*temp
   end do
  end do
  call enforceBCs(T1)
  end if
 end if
 if ((comp.eq.1).or.(tracer.eq.1)) then
  call compute_derivatives(u,DERr_u)
  call compute_derivatives(w,DERr_w)
 end if
 if ((comp.eq.1).or.(tracer.eq.1)) then
!$ ta1=omp_get_wtime()
   if (local_time.eq.0) then
    if (torder.eq.2) then
     dtstage=dt
    elseif (torder.eq.4) then
     dtstage=0.5d0*dt
    end if
   end if
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp,i,k) FIRSTPRIVATE(dtstage)
  do ID=1,ntr
   if (local_time.eq.1) then
    i=int(rtr0(ID)/dr,4)
    k=int(str0(ID)/ds,4)
    dt_tr(ID)=dt_array(i,k)
    if (torder.eq.2) then
     dtstage=dt_tr(ID)
    elseif (torder.eq.4) then
     dtstage=0.5d0*dt_tr(ID)
    end if
   end if
   call interpolate(rtr0(ID),str0(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr0(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   call interpolate(rtr0(ID),str0(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str0(ID),zs_interp)
    temps=temps/zs_interp
   end if
   if (torder.eq.2) then
    rtr(ID)=tempr
    str(ID)=temps
    if (tracer.eq.1) then
     call interpolate(rtr0(ID),str0(ID),DERr_tracer,tracer_space_array,tempT)
     Ttr(ID)=tempT
    end if
   elseif (torder.eq.4) then
    rtr(ID)=tempr
    str(ID)=temps
    if (tracer.eq.1) then
     call interpolate(rtr0(ID),str0(ID),DERr_tracer,tracer_space_array,tempT)
     Ttr(ID)=tempT
    end if
   end if
   rtr1(ID)=rtr0(ID)+dtstage*tempr    !!!stuff for all RK methods
   str1(ID)=str0(ID)+dtstage*temps
   if (tracer.eq.1) then
    Ttr1(ID)=Ttr0(ID)+dtstage*tempT
   end if
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
 end if
end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 1st Stage

if (torder.ge.2) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 2nd Stage
   if (tracer.eq.1) then
    call tracers_to_corners(rtr1,str1,Ttr1,T0)
    call load_tracer_space_array
    call compute_derivatives(tracer_space_array,DERr_tracer)
   end if
 if (RK_update.eq.1) then
  if ((comp.eq.1).and.(tracer.eq.0)) call tracers_to_corners(rtr1,str1,Ttr1,T0)
  if (tracer.eq.0) then
   call smoother_time_T(T1)
  end if
  call update_velocities
  if ((comp.eq.1).or.(tracer.eq.1)) then
   call compute_derivatives(u,DERr_u)
   call compute_derivatives(w,DERr_w)
  end if
 end if
 if (RaT.ne.0.d0) then
  if (tracer.eq.0) then
  do i=0,nx
   do k=0,nz
    temp=energy_space(i,k,T1)
    T(i,k)=T(i,k)+2.d0*temp
    if (local_time.eq.1) dt=dt_array(i,k)
    T2(i,k)=T0(i,k)+0.5d0*dt*temp
   end do
  end do
  call enforceBCs(T2)
  end if
 end if
 if ((comp.eq.1).or.(tracer.eq.1)) then
!$ ta1=omp_get_wtime()
   if (local_time.eq.0) then
    if (torder.eq.2) then
     dtstage=0.5d0*dt
    elseif (torder.eq.4) then
     dtstage=0.5d0*dt
    end if
   end if
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp) FIRSTPRIVATE(dtstage)
  do ID=1,ntr
   if (local_time.eq.1) dtstage=0.5d0*dt_tr(ID)
   call interpolate(rtr1(ID),str1(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr1(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   call interpolate(rtr1(ID),str1(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str1(ID),zs_interp)
    temps=temps/zs_interp
   end if
   if (torder.eq.2) then
    rtr(ID)=rtr(ID)+tempr
    str(ID)=str(ID)+temps
    if (tracer.eq.1) then
     call interpolate(rtr1(ID),str1(ID),DERr_tracer,tracer_space_array,tempT)
     Ttr(ID)=Ttr(ID)+tempT
    end if
   elseif (torder.eq.4) then
    rtr(ID)=rtr(ID)+2.d0*tempr
    str(ID)=str(ID)+2.d0*temps
    if (tracer.eq.1) then
     call interpolate(rtr1(ID),str1(ID),DERr_tracer,tracer_space_array,tempT)
     Ttr(ID)=Ttr(ID)+2.d0*tempT
    end if
   end if
   if (torder.eq.2) then !!final stage for RK2
    rtr(ID)=rtr0(ID)+dtstage*rtr(ID)
    str(ID)=str0(ID)+dtstage*str(ID)
    if (tracer.eq.1) then
     Ttr(ID)=Ttr0(ID)+dtstage*Ttr(ID)
    end if
   else
    rtr2(ID)=rtr0(ID)+dtstage*tempr
    str2(ID)=str0(ID)+dtstage*temps
    if (tracer.eq.1) then
     Ttr2(ID)=Ttr0(ID)+dtstage*tempT
    end if
   end if
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
 end if
end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 2nd Stage

if (torder.ge.3) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 3rd stage
  if (tracer.eq.1) then
   if (comp.eq.1) Cnew=C0 !!reset initial composition values for second RK4 stage in case of empty cells
   call tracers_to_corners(rtr2,str2,Ttr2,T0)
   call load_tracer_space_array
   call compute_derivatives(tracer_space_array,DERr_tracer)
  end if
 if (RK_update.eq.1) then
  if ((comp.eq.1).and.(tracer.eq.0)) then
   Cnew=C0 !!reset initial composition values for second RK4 stage in case of empty cells
   call tracers_to_corners(rtr2,str2,Ttr2,T0)
  end if
  if (tracer.eq.0) then
   call smoother_time_T(T2)
  end if
  call update_velocities
  if ((comp.eq.1).or.(tracer.eq.1)) then
   call compute_derivatives(u,DERr_u)
   call compute_derivatives(w,DERr_w)
  end if
 end if
 if (RaT.ne.0.d0) then
  if (tracer.eq.0) then
  do i=0,nx
   do k=0,nz
    temp=energy_space(i,k,T2)
    T(i,k)=T(i,k)+2.d0*temp
    if (local_time.eq.1) dt=dt_array(i,k)
    T1(i,k)=T0(i,k)+dt*temp  
   end do
  end do
  call enforceBCs(T1)
  end if
 end if
 if ((comp.eq.1).or.(tracer.eq.1)) then
  if (local_time.eq.0) dtt=dt
!$ ta1=omp_get_wtime()
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp) FIRSTPRIVATE(dtt)
  do ID=1,ntr
   if (local_time.eq.1) dtt=dt_tr(ID)
   call interpolate(rtr2(ID),str2(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr2(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   rtr(ID)=rtr(ID)+2.d0*tempr
   rtr1(ID)=rtr0(ID)+dtt*tempr
   call interpolate(rtr2(ID),str2(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str2(ID),zs_interp)
    temps=temps/zs_interp
   end if
   str(ID)=str(ID)+2.d0*temps
   str1(ID)=str0(ID)+dtt*temps
   if (tracer.eq.1) then
    call interpolate(rtr2(ID),str2(ID),DERr_tracer,tracer_space_array,tempT)
    Ttr(ID)=Ttr(ID)+2.d0*tempT
    Ttr1(ID)=Ttr0(ID)+dtt*tempT
   end if
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
 end if
end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 3rd Stage

if (torder.ge.4) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Begin 4th stage
  if (tracer.eq.1) then
   if (comp.eq.1) Cnew=C0 !!reset initial composition values for third RK4 stage in case of empty cells
   call tracers_to_corners(rtr1,str1,Ttr1,T0)
   call load_tracer_space_array
   call compute_derivatives(tracer_space_array,DERr_tracer)
  end if
 if (RK_update.eq.1) then
  if ((comp.eq.1).and.(tracer.eq.0)) then
   Cnew=C0 !!reset initial composition values for third RK4 stage in case of empty cells
   call tracers_to_corners(rtr1,str1,Ttr1,T0)
  end if
  if (tracer.eq.0) then
   call smoother_time_T(T1)
  end if
  call update_velocities
  if ((comp.eq.1).or.(tracer.eq.1)) then
   call compute_derivatives(u,DERr_u)
   call compute_derivatives(w,DERr_w)
  end if
 end if
 if (RaT.ne.0.d0) then
  if (tracer.eq.0) then
  do i=0,nx
   do k=0,nz
    temp=energy_space(i,k,T1)
    if (local_time.eq.1) dt=dt_array(i,k)
    T(i,k)=T0(i,k)+(dt/6.d0)*(T(i,k)+temp)
   end do
  end do
  call enforceBCs(T)
  end if
 end if
 if ((comp.eq.1).or.(tracer.eq.1)) then
!$ ta1=omp_get_wtime()
  if (local_time.eq.0) dto6=dt/6.d0
!$OMP PARALLEL DO PRIVATE(tempr,temps,tempT,xr_interp,zs_interp) FIRSTPRIVATE(dto6)
  do ID=1,ntr
   if (local_time.eq.1) dto6=dt_tr(ID)/6.d0
   call interpolate(rtr1(ID),str1(ID),DERr_u,u,tempr)
   if (gridX.ne.0) then
    call interpolate_xr(rtr1(ID),xr_interp)
    tempr=tempr/xr_interp
   end if
   rtr(ID)=rtr0(ID)+dto6*(rtr(ID)+tempr)
   call interpolate(rtr1(ID),str1(ID),DERr_w,w,temps)
   if (gridZ.ne.0) then
    call interpolate_zs(str1(ID),zs_interp)
    temps=temps/zs_interp
   end if
   str(ID)=str0(ID)+dto6*(str(ID)+temps)
   if (tracer.eq.1) then
    call interpolate(rtr1(ID),str1(ID),DERr_tracer,tracer_space_array,tempT)
    Ttr(ID)=Ttr0(ID)+dto6*(Ttr(ID)+tempT)
   end if
  end do
!$OMP END PARALLEL DO
!$ ta2=omp_get_wtime()
!$ tadvect_tracer=tadvect_tracer+ta2-ta1
 end if
end if!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End 4th stage

!!!!!!wrap up: convert tracers to corner values
 if (comp.eq.1) Cnew=C0 !!reset initial composition values for final RK4 stage in case of empty cells
 if ((comp.eq.1).or.(tracer.eq.1)) call tracers_to_corners(rtr,str,Ttr,T0)
 if (tracer.eq.1) T=Tratio
 if (tracer.eq.0) call smoother_time_T(T)
!$ time2=omp_get_wtime()
!$ tenergy=tenergy+time2-time1
end

subroutine enforceBCs(T0)
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: T0(-span:nx+span,-span:nz+span)
do i=1,span
 T0(-i,0:nz)=T0(i,0:nz) !!symmetry
 T0(nx+i,0:nz)=T0(nx-i,0:nz)
end do
do k=1,span
 if (Tbot.eq.0) then !!isothermal
  T0(:,-k)=-(T0(:,k)-1.d0)+1.d0 !!antisymmetry
 elseif (Tbot.eq.1) then
  T0(:,-k)=T0(:,k) !!symmetry
 end if
 T0(:,nz+k)=-T0(:,nz-k)
end do
 if (Tbot.eq.0) T0(:,0)=1.d0
T0(:,nz)=0.d0
end


subroutine extendT(T0) !!define Textend for interpolation
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: T0(-span:nx+span,-span:nz+span)
 Textend(0:nx,0:nz)=T0(0:nx,0:nz)
 do i=1,spanT
  Textend(-i,0:nz)=Textend(i,0:nz) !!symmetry
  Textend(nx+i,0:nz)=Textend(nx-i,0:nz)
 end do
 do k=1,spanT
  if (Tbot.eq.0) then !!isothermal
   Textend(:,-k)=-(Textend(:,k)-1.d0)+1.d0 !!antisymmetry
  elseif (Tbot.eq.1) then
   Textend(:,-k)=Textend(:,k) !!symmetry
  end if
  Textend(:,nz+k)=-Textend(:,nz-k)
 end do
  if (Tbot.eq.0) Textend(:,0)=1.d0
 Textend(:,nz)=0.d0
end

subroutine derivative_eigenvalues
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer :: i,j,k,n_angle,angle,angle2,jj
real*8 :: kr,ks,d_angle
complex*16 :: T_test,Tr,Trr,Trrr,Trrrr,Ts,Tss,Tsss,Tssss,Trs,Trrs,Trrss,temp

!n_angle=2**(20)+1   
n_angle=2**(11)+1
d_angle=2.d0*pii/real(n_angle-1,8)
Tr_max=0.d0
Trr_max=0.d0
Trrr_max=0.d0
Trrrr_max=0.d0
Ts_max=0.d0
Tss_max=0.d0
Tsss_max=0.d0
Tssss_max=0.d0
Trs_max=0.d0
Trrs_max=0.d0
Trss_max=0.d0
Trrss_max=0.d0

i=0 !!central position of waveform
k=0


!$OMP PARALLEL DO PRIVATE(angle,j,jj,kr,ks,Tr,Trr,Trrr,Trrrr,Trs,Trrs,Trrss,temp) &
!$OMP& REDUCTION(max:Tr_max,Trr_max,Trrr_max,Trrrr_max,Trs_max,Trrs_max,Trrss_max)
do angle=1,n_angle
 ks=0.d0
 kr=real(angle-1,8)*d_angle
 Tr=(0.d0,0.d0)
 Trr=(0.d0,0.d0)
 Trrr=(0.d0,0.d0)
 Trrrr=(0.d0,0.d0)
 do j=-span1,span1
  temp=T_test(j+i,k,kr,ks)
  Tr=Tr+D1(j)*temp
  Trr=Trr+D2(j)*temp
 end do
 do j=-span,span
  temp=T_test(j+i,k,kr,ks)
  Trrr=Trrr+D3(j)*temp
  Trrrr=Trrrr+D4(j)*temp
 end do
 Tr=Tr/T_test(i,k,kr,ks)     !!divide out central value of test temperature to get eigenvalues
 Trr=Trr/T_test(i,k,kr,ks)
 Trrr=Trrr/T_test(i,k,kr,ks)
 Trrrr=Trrrr/T_test(i,k,kr,ks)
 if (dabs(dimag(Tr)).gt.Tr_max) Tr_max=dabs(dimag(Tr))
 if (dabs(real(Trr,8)).gt.Trr_max) Trr_max=dabs(real(Trr,8))
 if (dabs(dimag(Trrr)).gt.Trrr_max) Trrr_max=dabs(dimag(Trrr))
 if (dabs(real(Trrrr,8)).gt.Trrrr_max) Trrrr_max=dabs(real(Trrrr,8))
 do angle2=1,n_angle
  ks=real(angle2-1,8)*d_angle
  Trs=(0.d0,0.d0)
  Trrs=(0.d0,0.d0)
  Trrss=(0.d0,0.d0)
  do j=-span1,span1
   do jj=-span1,span1
    temp=T_test(j+i,jj+k,kr,ks)
    Trs=Trs+Drs(j,jj)*temp
    Trrs=Trrs+Drrs(j,jj)*temp
    Trrss=Trrss+Drrss(j,jj)*temp
   end do
  end do
  if (dabs(real(Trs,8)).gt.Trs_max) Trs_max=dabs(real(Trs,8))
  if (dabs(dimag(Trrs)).gt.Trrs_max) Trrs_max=dabs(dimag(Trrs))
  if (dabs(real(Trrss,8)).gt.Trrss_max) Trrss_max=dabs(real(Trrss,8))
 end do
end do
!$OMP END PARALLEL DO
Tr_base=Tr_max  !!!base values that can be scaled by the denominators of FD formulas (scaling depends on grid size)
Trr_base=Trr_max
Trrr_base=Trrr_max
Trrrr_base=Trrrr_max
Ts_base=Ts_max
Tss_base=Tss_max
Tsss_base=Tsss_max
Tssss_base=Tssss_max
Trs_base=Trs_max
Trrss_base=Trrss_max
Trss_base=Trss_max
Trrs_base=Trrs_max
end

complex*16 function T_test(i,k,kr,ks)  !!for Fourier analysis
 use basics
implicit none
integer*4 :: i,k
real*8 :: r,kr,s,ks
r=real(i,8) !!position of test waveform
s=real(k,8)
T_test=cmplx(dcos(kr*r),dsin(kr*r))*cmplx(dcos(ks*s),dsin(ks*s))
end

subroutine stability_RK4
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: eigA,eigA_max,eigD,eigD_max,dtA,dtD

!!!!!!!scale derivative eigenvalues
Ts_max=Tr_base/ds     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2 !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr
Trr_max=-Trr_base/dr2


!!!!!!!!!!!!!!! find worst case scenario for stability eigenvalues (based on transformed energy equation)
if (RaT.ne.0.d0) then
 eigA_max=0.d0
 eigD_max=0.d0
 do i=0,nx
  do k=0,nz
   if (eigA(i,k).gt.eigA_max) eigA_max=eigA(i,k)
   if (eigD(i,k).gt.eigD_max) eigD_max=eigD(i,k)
  end do
 end do

 if (torder.eq.2) then
  if (eigA_max.ge.eigD_max) then
   write(*,*) "Advective Eigenvalues Exceed Diffusive Eigenvalues -- this does not work for RK2"
   stop
  end if
  dtA=1.0d0/eigA_max          !!diffusive stability bound of 2.0 for RK2 (unstable for pure advection so I'm being conservative with that bound)
  dtD=2.0d0/eigD_max
 elseif (torder.eq.4) then
  dtA=2.8d0/eigA_max          !!approx. 2.8 stability bound for RK4 (conditionally stable for both diffusion and advection)
  dtD=2.8d0/eigD_max
 end if

 if (eigA_max.eq.0.d0) then !!if no mathematical advection
  dt=dtD
 else
  dt=1.d0/(1.d0/dtA+1.d0/dtD)
 end if
else !!for pure compositional convection -- no stability requirement for C -- time steps must be small enough for accuracy (and so tracers do not exit the domain)
 dt=dt_manual
end if
dt=courant*dt
end

subroutine local_time_steps
 use basics
 use arrays
implicit none
integer*4 :: i,k
real*8 :: eigA,eigA_max,eigD,eigD_max,dtA,dtD,eA,eD
!!!!!!!scale derivative eigenvalues
Ts_max=Tr_base/ds     !!scale according to the denominators of finite difference formulas
Tss_max=-Trr_base/ds2 !!even # of derivative operators=real, odd #=imaginary. For real values keep signs (sign=(imaginary unit)^(# of derivatives)), for imaginary eigenvalues take absolute values
Tr_max=Tr_base/dr
Trr_max=-Trr_base/dr2

do i=0,nx
 do k=0,nz
  if (torder.eq.4) then
   dtA=2.8d0/eigA(i,k)          !!approx. 2.8 stability bound
   dtD=2.8d0/eigD(i,k)
  elseif (torder.eq.2) then
   eA=eigA(i,k)
   eD=eigD(i,k)
   if (eA.ge.eD) then
    write(*,*) "Advective Eigenvalues Exceed Diffusive Eigenvalues -- this does not work for RK2"
    stop
   end if
   dtA=1.0d0/eA          !!unstable for pure advection - taking a conservative limit
   dtD=2.0d0/eD          !!2.0 diffusive stability bound
  end if
   if (eigA(i,k).eq.0) then
    dt_array(i,k)=dtD
   elseif (eigD(i,k).eq.0) then
    dt_array(i,k)=dtA
   else
    dt_array(i,k)=1.d0/(1.d0/dtA+1.d0/dtD)
   end if
 end do
end do
end

real*8 function eigA(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
 if (tracer.eq.0) then !!field method for solving the energy eqn -- physical advection terms are present
  eigA=dabs(-conduct(i,k)*xrr(i)/(xr(i)**3.d0)-u(i,k)/xr(i)+conduct_r(i,k)/(xr(i)**2.d0))*Tr_max+&
&dabs(-conduct(i,k)*zss(k)/(zs(k)**3.d0)-w(i,k)/zs(k)+conduct_s(i,k)/(zs(k)**2.d0))*Ts_max
 elseif (tracer.eq.1) then !!Lagrangian method -- no physical advection terms (some mathematical advection from nonuniform grids and variable conductivity though)
  eigA=dabs(-conduct(i,k)*xrr(i)/(xr(i)**3.d0)+conduct_r(i,k)/(xr(i)**2.d0))*Tr_max+&
&dabs(-conduct(i,k)*zss(k)/(zs(k)**3.d0)+conduct_s(i,k)/(zs(k)**2.d0))*Ts_max
 end if
end

real*8 function eigD(i,k)
 use basics
 use arrays
implicit none
integer*4 :: i,k
  eigD=conduct(i,k)*(dabs(Trr_max/(xr(i)**2.d0))+dabs(Tss_max/(zs(k)**2.d0)))
end
