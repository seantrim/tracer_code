module tracers
 type trace_cell
  integer*4, allocatable :: id(:)
 end type trace_cell
end module tracers

subroutine mzran(ran_real,rank)
use basics
implicit none
integer*4 :: ran_int,rank
real*8 :: ran_real
ran_int=seed_i(rank)-seed_k(rank)
if (ran_int.lt.0) ran_int=ran_int+2147483579
seed_i(rank)=seed_j(rank)
seed_j(rank)=seed_k(rank)
seed_k(rank)=ran_int
seed_n(rank)=69069*seed_n(rank)+1013904243
ran_int=ran_int+seed_n(rank)
ran_real=0.5d0+(0.2328306d-9)*real(ran_int,8)
end

subroutine shuffle_serial(a,istart,N,rank) !!for small arrays (can be called by more than one thread simultaneously)
implicit none
integer*4 :: i,N,a(istart:N),rank,r,temp,istart
real*8 :: ran_real
do i=istart,N-1
 call mzran(ran_real,rank)
 r=floor(ran_real*real(N-i,8),4)+i+1   !!random integer from i+1 to N
 temp=a(i)
 a(i)=a(r)
 a(r)=temp
end do
end

subroutine shuffle(a,N) !!for a large array "a" of size "N" ------------------------------------ Maybe use random reals instead of integers to reduce repaeted random values????
!$ use OMP_LIB
implicit none
integer*4 :: i,N,rank,a(1:N),b(1:N),r(1:N),mloc_array(1:1),mloc,Nthreads,Nshare
integer*4, allocatable :: Nmin(:),Nmax(:),mloc_rank(:)
real*8 :: ran_real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!create array of random integers (repeats possible)
!$OMP PARALLEL PRIVATE(rank,i)
rank=omp_get_thread_num()
!$OMP DO
do i=1,N
 call mzran(ran_real,rank)
! r(i)=nint(ran_real*real(N-1,8),4)+1   !!random integer from 1 to N !!!!!!!!!!!!!!!!!!!adapt to use floor instead of nint (so boundary values are equally likely compared to interior values)
 r(i)=floor(ran_real*real(N,8),4)+1 !!!!!!!!!!!!!use this one (untested)
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end create array of random integers

!!!!!!!!!!!!!!!!!!!!!!!!!figure out how to divide work between processes
Nthreads=omp_get_max_threads()
Nshare=nint(real(N,8)/real(Nthreads,8),4)
allocate(Nmin(0:Nthreads-1),Nmax(0:Nthreads-1),mloc_rank(0:Nthreads-1))
!$OMP PARALLEL DO
do rank=0,Nthreads-1
 Nmin(rank)=rank*Nshare
 Nmax(rank)=(rank+1)*Nshare-1
end do
!$OMP END PARALLEL DO
Nmin(0)=1
Nmax(Nthreads-1)=N
!!!!!!!!!!!!!!!!!!!!!!!!!end figure out how to divide work between processes

do i=1,N
!$OMP PARALLEL PRIVATE(rank,mloc_array)
 rank=omp_get_thread_num()
 mloc_array=minloc(r(Nmin(rank):Nmax(rank)))
 mloc_rank(rank)=mloc_array(1)+Nmin(rank)-1
!$OMP END PARALLEL
 mloc_array=minloc(r(mloc_rank))
 rank=mloc_array(1)-1 !!rank with smallest minvalue
 mloc=mloc_rank(rank)
 b(i)=a(mloc)
 r(mloc)=N+1 !!once random number is used exclude from the search
end do
a=b
end

subroutine equalize_tracers
!$ use OMP_LIB
 use basics
 use arrays
 use tracers
implicit none
integer*4 :: ID,p,i,k,pmax,pp,tt,ii,kk
real*8 :: perturb_r,perturb_s,r,s,tempT,tempC,ran_i,ran_k,time1,time2
real*8 :: ran_q
real*8 :: Tdumb(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
integer*4 :: nexcess(0:nx,0:nz),nfree,nadd,ncells,q,qmax
real*8 :: ran_real
integer*4 :: rank,counter(0:nx,0:nz),type_save,escape,minexcess(0:ntype),tpc_type(0:ntype),count_save
type (trace_cell) tracer_in_cell(0:nx,0:nz)
integer*4 :: Nthreads,Nshare,nx_loop(0:nx),nz_loop(0:nz),count_display(0:ntype),Nideal(0:ntype),tracer_count(0:ntype)
integer*4, allocatable :: nx_min(:),nx_max(:),donor(:,:,:),index_save(:,:),index_save2(:,:),work(:),nfree_rank(:,:)
logical :: logic_test(0:nx,0:nz)
!$ time1=omp_get_wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Testing
if (comp.eq.1) then
 do tt=1,ntype
  Nideal(tt)=nint(Cinit(tt)*real(ntr,8)) 
 end do
  Nideal(0)=ntr-sum(Nideal(1:ntype))
else
  Nideal(0)=ntr
end if
if (time.eq.0.d0) write(*,*) ntr,Nideal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End testing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Which cells have excess tracers?? we want nexcess=0 ideally
if (comp.eq.0) then
 tpc_type(0)=tpc
!$OMP PARALLEL DO PRIVATE(i,k)
 do i=0,nx
  do k=0,nz
   nexcess(i,k)=0
   if (nempty(i,k).gt.0) then
    nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
   end if
  end do
 end do
!$OMP END PARALLEL DO
elseif (comp.eq.1) then
 logic_test=sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0
 count_save=count(logic_test)
 if (count_save.gt.0) then
  tpc_type(0)=nint(real(sum(count_type(0,0:nx,0:nz),MASK=logic_test),8)/&
  &real(count_save,8))
 else
  tpc_type(0)=0
 end if
 do tt=1,ntype
  logic_test=Cnew(tt,0:nx,0:nz).eq.1.d0
  count_save=count(logic_test)
  if (count_save.gt.0) then
   tpc_type(tt)=nint(real(sum(count_type(tt,0:nx,0:nz),MASK=Cnew(tt,0:nx,0:nz).eq.1.d0),8)/real(count_save,8))
  else
   tpc_type(tt)=0
  end if
 end do
!$OMP PARALLEL DO PRIVATE(i,k,tt)
 do i=0,nx
  do k=0,nz
!!are the tracers all of one type? Don't modify interface or near interface cells or cells near empty cells
   nexcess(i,k)=0
   do tt=0,ntype
    if ((all(count_type(tt,i-1:i+1,k-1:k+1).eq.nempty(i-1:i+1,k-1:k+1))).and.&
&(nempty(i,k).gt.0)) then
!     nexcess(i,k)=nempty(i,k)-tpc  !!# of excess tracers
     nexcess(i,k)=nempty(i,k)-tpc_type(tt)  !!# of excess tracers
     exit
    end if
   end do
  end do
 end do
!$OMP END PARALLEL DO
end if
nfree=sum(nexcess,MASK=nexcess.gt.0) !!# of freed tracers to play with
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Which cells have excess tracers?? we want nexcess=0 ideally

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Is equalize worth doing? If not exit
if (all(nexcess(0:nx,0:nz).eq.0)) then
 return
else
 nequalize=nequalize+1
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Is equalize worth doing?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Initial Prep
Nthreads=omp_get_max_threads()
Nshare=nx/Nthreads
allocate(nx_min(0:Nthreads-1),nx_max(0:Nthreads-1),index_save(0:Nthreads-1,0:ntype),nfree_rank(0:Nthreads-1,0:ntype),&
&index_save2(0:Nthreads-1,0:ntype))
do rank=0,Nthreads-1
 nx_min(rank)=rank*Nshare
 nx_max(rank)=(rank+1)*Nshare-1
end do
nx_max(Nthreads-1)=nx

!$OMP PARALLEL DO PRIVATE(i,k)
do i=0,nx
 do k=0,nz
  allocate(tracer_in_cell(i,k)%id(0:nempty(i,k)))
  tracer_in_cell(i,k)%id(0)=0
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Initial Prep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Organize tracer IDs by CV location (i,k)
counter=0                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Can Parallelization be more efficient??
!$OMP PARALLEL PRIVATE(rank,ID,i,k)
rank=omp_get_thread_num()
do ID=1,ntr
 i=nint(rtr(ID)/dr,4)
 if ((nx_min(rank).le.i).and.(i.le.nx_max(rank))) then
  k=nint(str(ID)/ds,4)
  counter(i,k)=counter(i,k)+1             !!!!!!!!!!!!!!!!!!!!!maybe reuse nempty(i,k) since final version of counter(i,k) is identical -- would save memory
  tracer_in_cell(i,k)%id(counter(i,k))=ID
 end if
end do
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Organize tracer IDs by CV location (i,k)

!!!!!!!!!!!!!!!!!!!!!ADD SHUFFLE HERE TO TRACER_IN_CELL ARRAY
!$OMP PARALLEL PRIVATE(i,k,rank)
rank=omp_get_thread_num()
!$OMP DO
do i=0,nx
 do k=0,nz
  call shuffle_serial(tracer_in_cell(i,k)%id(1:counter(i,k)),1,counter(i,k),rank)
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!END SHUFFLE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Prep for T interpolation
if (tracer.eq.1) then
 call extendT(Tratio) !!for interpolation of T to Ttr during equalizing
 call compute_derivatives(Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp),DERr_Textend)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Prep for T interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Look for donor tracers: NEWEST
allocate(donor(0:Nthreads-1,0:nfree,0:ntype))   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Save on memory here???? 0:nfree may be too big a range
donor=0  !!default values
index_save=0 !!staring value
!!!!!!!!!!!!!!!!Look near edges first
do i=0,nx,nx  !!!left and right sides
!$OMP PARALLEL PRIVATE(rank,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do k=0,nz
 if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
  do q=1,counter(i,k)
    ID=tracer_in_cell(i,k)%id(q)
    if ((str(ID).lt.dist_s).or.(str(ID).gt.(1.d0-dist_s)).or.&
&(rtr(ID).lt.dist_r).or.(rtr(ID).gt.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
    end if
  end do
 end if
end do
!$OMP END DO
!$OMP END PARALLEL
end do

do k=0,nz,nz  !!!bottom and top boundaries
!$OMP PARALLEL PRIVATE(rank,i,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do i=1,nx-1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! corners already done
 if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
  do q=1,counter(i,k)
    ID=tracer_in_cell(i,k)%id(q)
    if ((str(ID).lt.dist_s).or.(str(ID).gt.(1.d0-dist_s)).or.&
&(rtr(ID).lt.dist_r).or.(rtr(ID).gt.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
    end if
  end do
 end if
end do
!$OMP END DO
!$OMP END PARALLEL
end do
!!!!!!!!!!!!!!!!End Look near edges first

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!shuffle donor tracers near the boundaries
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 call shuffle_serial(donor(rank,1:index_save(rank,tt),tt),1,index_save(rank,tt),rank)
end do
!$OMP END PARALLEL
index_save2=index_save !!index_save should be saved here for later shuffle of interior donor tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end shuffle tracers near the boundaries

!!!!!!!!!!!!!!!!!!!!!!!!!!wrap up edge CVs (now look for excess tracers in edge CV interiors, if any)
do i=0,nx,nx !!left and right side interiors
!$OMP PARALLEL PRIVATE(rank,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
 do k=0,nz
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     if (((str(ID).ge.dist_s).and.(str(ID).le.(1.d0-dist_s))).and.&
&((rtr(ID).ge.dist_r).and.(rtr(ID).le.(aspect-dist_r)))) then !!too many tracers?? These tracers will be redistributed
      index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
      donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
      nexcess(i,k)=nexcess(i,k)-1
      if (nexcess(i,k).eq.0) exit
     end if
   end do
  end if
 end do
!$OMP END DO
!$OMP END PARALLEL
end do

do k=0,nz,nz !!bottom and top boundary CV interiors
!$OMP PARALLEL PRIVATE(rank,i,q,ID)
rank=omp_get_thread_num()
!$OMP DO
 do i=1,nx-1 !!already did corners
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     if (((str(ID).ge.dist_s).and.(str(ID).le.(1.d0-dist_s))).and.&
&(rtr(ID).ge.dist_r).and.(rtr(ID).le.(aspect-dist_r))) then !!too many tracers?? These tracers will be redistributed:-------------------------- Cut down if statement (no corners)
      index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
      donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
      nexcess(i,k)=nexcess(i,k)-1
      if (nexcess(i,k).eq.0) exit
     end if
   end do
  end if
 end do
!$OMP END DO
!$OMP END PARALLEL
end do
!!!!!!!!!!!!!!!!!!!!!!!!!End wrap up edge cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Search remaining (interior) CVs for donor tracers
!$OMP PARALLEL PRIVATE(rank,i,k,q,ID)
rank=omp_get_thread_num()
!$OMP DO
do i=1,nx-1 !!left and right side interiors
 do k=1,nz-1
  if (nexcess(i,k).gt.0) then !!only take donors from rich unrestricted CVs
   do q=1,counter(i,k)
     ID=tracer_in_cell(i,k)%id(q)
     index_save(rank,Ttype(ID))=index_save(rank,Ttype(ID))+1
     donor(rank,index_save(rank,Ttype(ID)),Ttype(ID))=ID
     nexcess(i,k)=nexcess(i,k)-1
     if (nexcess(i,k).eq.0) exit
   end do
  end if
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Search remaining (interior) CVs for donor tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Test to see if we got all the donor tracers
nfree_rank=0
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 nfree_rank(rank,tt)=count(donor(rank,:,tt).gt.0)
end do
!$OMP END PARALLEL
if (sum(nfree_rank).ne.nfree) then
 write(*,*) "Error: # of donor tracers found inconsistent with nfree value"
 stop
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End test to see if we got all the donor tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Look for donor tracers: NEWEST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************Shuffle donor tracers that are not near the boundaries
!$OMP PARALLEL PRIVATE(rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 call shuffle_serial(donor(rank,index_save2(rank,tt)+1:index_save(rank,tt),tt),index_save2(rank,tt)+1,index_save(rank,tt),rank)
end do
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***********************End Shuffle donor tracers that are not near the boundaries


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Add needed tracers to lacking cells via lottery: NEW
!$OMP PARALLEL DO !!!!!prep for randomized loops
do i=0,nx
 nx_loop(i)=i
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
do k=0,nz
 nz_loop(k)=k
end do
!$OMP END PARALLEL DO

index_save=0 !!to keep track of which tracer IDs need to be used for T interpolation
escape=0
do tt=0,ntype
 tracer_count(tt)=sum(count_type(tt,0:nx,0:nz))
end do
if (all(tracer_count.eq.Nideal)) escape=1

do
if (escape.eq.0) then
 do tt=0,ntype
  tracer_count(tt)=sum(count_type(tt,0:nx,0:nz))
 end do
 if (all(tracer_count.eq.Nideal)) escape=1
end if
if (tpc_type(0).gt.0) then
 if (comp.eq.0) then
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
 else
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 end if
else
 minexcess(0)=0
end if
do p=1,ntype
 if (tpc_type(p).gt.0) then
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 else
  minexcess(p)=0
 end if
end do
if (all(minexcess.eq.0)) exit
if (sum(index_save(:,:)).ge.nfree) exit
index_save2=index_save
rank=0
 call shuffle_serial(nx_loop,0,nx,rank) !!randomize loop order
 call shuffle_serial(nz_loop,0,nz,rank) !!randomize loop order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Find poorest CV locations
!$OMP PARALLEL PRIVATE(i,k,ii,kk,rank,r,s,perturb_r,perturb_s,ID,tt,p,type_save)
rank=omp_get_thread_num()
!$OMP DO
do ii=0,nx   !!find which cells lack tracers
 i=nx_loop(ii) !!!use randomized i values
 do kk=0,nz
  k=nz_loop(kk) !!!use randomized k values
  type_save=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
  if (escape.eq.0) then
   do p=0,ntype !!!!!!!!!!!!!!!!decide which type of tracer to donate in order to get optimum balance
    if ((tracer_count(p)).gt.Nideal(p)) then
     tt=p
     exit
    end if
   end do
  else
   tt=type_save !!if ideal values are met then do not swap tracer types
  end if
  if ((nexcess(i,k).eq.minexcess(type_save)).and.(index_save(rank,tt).lt.nfree_rank(rank,tt)).and.(nexcess(i,k).lt.0)) then
   r=dr*(real(i,8)) !!cell corners
   s=ds*(real(k,8))
   call mzran(perturb_r,rank)
   call mzran(perturb_s,rank)
   perturb_r=(perturb_r-0.5d0)*dr
   perturb_s=(perturb_s-0.5d0)*ds
   index_save(rank,tt)=index_save(rank,tt)+1
   ID=donor(rank,index_save(rank,tt),tt)
   rtr(ID)=r+perturb_r
   str(ID)=s+perturb_s
   if (rtr(ID).le.0.d0) rtr(ID)=dist_r
   if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
   if (str(ID).le.0.d0) str(ID)=dist_s
   if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
   if (comp.eq.1) then
    Ttype(ID)=type_save !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
   end if
   count_type(tt,i,k)=count_type(tt,i,k)-1
   count_type(Ttype(ID),i,k)=count_type(Ttype(ID),i,k)+1
   nexcess(i,k)=nexcess(i,k)+1
  end if
 end do
end do
!$OMP END DO
!$OMP END PARALLEL
 if (all(index_save.eq.index_save2)) exit !!if no repositioning happened on the last iteration
end do

!!!!!!!!!!!!!!!!!!Use any remaining donor tracers in serial (alternating between ranks to get even usage of donors)
rank=0
do
if (tpc_type(0).gt.0) then
 if (comp.eq.0) then
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
 else
  minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 end if
else
 minexcess(0)=0
end if
do p=1,ntype
 if (tpc_type(0).gt.0) then
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 else
  minexcess(p)=0
 end if
end do
if (all(minexcess.eq.0)) exit
if (sum(index_save(:,:)).ge.nfree) exit
index_save2=index_save
rank=0
 call shuffle_serial(nx_loop,0,nx,rank) !!randomize loop order
 call shuffle_serial(nz_loop,0,nz,rank) !!randomize loop order
do ii=0,nx   !!find which cells lack tracers
 i=nx_loop(ii) !!!use randomized i values
 do kk=0,nz
  k=nz_loop(kk) !!!use randomized k values
  tt=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
  if ((nexcess(i,k).eq.minexcess(tt)).and.(sum(index_save(:,tt)).lt.sum(nfree_rank(:,tt))).and.(nexcess(i,k).lt.0)) then
   do
    rank=rank+1
    if (rank.eq.Nthreads) rank=0
    if (index_save(rank,tt).lt.nfree_rank(rank,tt)) exit
   end do
   r=dr*(real(i,8)) !!cell corners
   s=ds*(real(k,8))
   call mzran(perturb_r,rank)
   call mzran(perturb_s,rank)
   perturb_r=(perturb_r-0.5d0)*dr
   perturb_s=(perturb_s-0.5d0)*ds
   index_save(rank,tt)=index_save(rank,tt)+1
   ID=donor(rank,index_save(rank,tt),tt)
   rtr(ID)=r+perturb_r
   str(ID)=s+perturb_s
   if (rtr(ID).le.0.d0) rtr(ID)=dist_r
   if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
   if (str(ID).le.0.d0) str(ID)=dist_s
   if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
   if (comp.eq.1) then
    Ttype(ID)=maxloc(count_type(0:ntype,i,k),DIM=1)-1 !!which tracer type is in that cell? Needed to subtract 1 since index starts at 0 in count_type array
   end if
   nexcess(i,k)=nexcess(i,k)+1
  end if
 end do
end do
 if (all(index_save.eq.index_save2)) exit !!if no repositioning happened on the last iteration
end do
!!!!!!!!!!!!!!!!!!End use any remaining donor tracers in serial 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Find poorest CV locations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interpolate T values of repositioned tracers
if (tracer.eq.1) then
Tdumb=Textend(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
!$OMP PARALLEL PRIVATE(ID,p,rank,tt)
rank=omp_get_thread_num()
do tt=0,ntype
 do p=1,index_save(rank,tt) !!!interpolate temperature values to repositioned tracers
  ID=donor(rank,p,tt)
  call interpolate(rtr(ID),str(ID),DERr_Textend,Tdumb,Ttr(ID))
 end do
end do
!$OMP END PARALLEL
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Interpolate T values of repositioned tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!End Add needed tracers to lacking cells via lottery: NEW

if (comp.eq.0) then
 minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=(nexcess(0:nx,0:nz).le.0)) !!ambient composition
else
 minexcess(0)=minval(nexcess(0:nx,0:nz),MASK=((sum(Cnew(1:ntype,0:nx,0:nz),DIM=1).eq.0.d0).and.(nexcess(0:nx,0:nz).le.0))) !!ambient composition
 do p=1,ntype
  minexcess(p)=minval(nexcess(0:nx,0:nz),MASK=((Cnew(p,0:nx,0:nz).eq.1.d0).and.(nexcess(0:nx,0:nz).le.0))) !!enriched compositions
 end do
end if

do tt=0,ntype
 count_display(tt)=sum(count_type(tt,0:nx,0:nz))
end do
if (comp.eq.0) then
 write(1024,'(g20.8,15(i10))') time,tpc,minexcess,sum(nexcess(0:nx,0:nz),MASK=nexcess(0:nx,0:nz).le.0),&
 &sum(index_save),nfree,count_display
else
 write(1024,'(g20.8,15(i10))') time,tpc_type,minexcess,sum(nexcess(0:nx,0:nz),MASK=nexcess(0:nx,0:nz).le.0),&
 &sum(index_save),nfree,count_display
end if
!$ time2=omp_get_wtime()
smooth_vis=0 !!viscosity smoothing not required yet
call tracers_to_corners(rtr,str,Ttr,T)
smooth_vis=1
!$ tequalize=tequalize+time2-time1
end

subroutine entrainment_setup
 use basics
implicit none
integer*4 :: k
real*8 :: scell,zgrid
  do k=0,nz
   scell=ds*(real(k,8))
    if (zgrid(scell).le.dent) kent=k
  end do
end

subroutine initialize_tracers
 use basics
 use arrays
implicit none
integer*4 :: i,k,p,id,nseed,composition
integer*4, allocatable :: seed(:)
real*8 :: rcell,scell,perturb_r,perturb_s,xtr,ztr
real*8 :: xgrid,zgrid,temperature
 call random_seed(size=nseed)
 allocate(seed(1:nseed))
 seed=1
 call random_seed(PUT=seed)
 id=1 !!tracer counter: should range from 1 to ntr
 do i=0,nx !!go through each cell corner
 rcell=dr*(real(i,8))
  do k=0,nz
   scell=ds*(real(k,8))
   if (comp.eq.1) then
    if (zgrid(scell).le.dent) kent=k
   end if
   do p=1,tpc
    call random_number(perturb_r)
    call random_number(perturb_s)
    perturb_r=(perturb_r-0.5d0)*dr
    perturb_s=(perturb_s-0.5d0)*ds
    rtr(id)=rcell+perturb_r   
    str(id)=scell+perturb_s
    if (rtr(ID).le.0.d0) rtr(ID)=dist_r
    if (rtr(ID).ge.aspect) rtr(ID)=aspect-dist_r
    if (str(ID).le.0.d0) str(ID)=dist_s
    if (str(ID).ge.1.d0) str(ID)=1.d0-dist_s
    xtr=xgrid(rtr(id))       !!!!!maybe interpolate these?? minor detail... 
    ztr=zgrid(str(id))
    if (comp.eq.1) then
     Ttype(id)=composition(xtr,ztr)
    else
     Ttype(id)=0  !!!uniform composition
    end if
    if (tracer.eq.1) then
     Ttr(id)=temperature(xtr,ztr)
    end if
    id=id+1
   end do
  end do
 end do
 call tracers_to_corners(rtr,str,Ttr,T)!!!!!!!!!!convert tracers to initial C field here
end

subroutine add_tracers
 use basics
 use arrays
implicit none
integer*4 :: i,k,id,nseed,tt,composition
integer*4, allocatable :: seed(:)
real*8 :: scell,xtr,ztr
real*8 :: xgrid,zgrid,temperature
do id=1,ntr
 xtr=xgrid(rtr(id))       !!!!!maybe interpolate these?? minor detail... 
 ztr=zgrid(str(id))
 Ttype(id)=composition(xtr,ztr)
end do
 if (comp.eq.1) then
  Cnew=0.d0
  do k=0,nz
   scell=ds*(real(k,8))
   if (zgrid(scell).le.dent) kent=k
   do i=0,nx
    do tt=1,ntype
     if (composition(xg(i),zg(k)).eq.tt) Cnew(tt,i,k)=1.d0
    end do
   end do
  end do
 end if
 call tracers_to_corners(rtr,str,Ttr,T)!!!!!!!!!!convert tracers to initial C field here
end

subroutine tracers_to_corners(rtr_,str_,Ttr_,T0) !!T0 input values are for empty cells
!$ use OMP_LIB
 use basics
 use arrays
implicit none
real*8 :: r,s,T0(-span:nx+span,-span:nz+span),T1(-span:nx+span,-span:nz+span)
real*8 :: temp,rtr_(1:ntr),str_(1:ntr),Ttr_(1:ntr),time1,time2
integer*4 :: id,i,k,tt
integer*4 :: rank,Nshare,Nthreads,id_min,id_max
real*8, allocatable :: Twork(:,:,:)
integer*4, allocatable :: count_work(:,:,:,:)
!$ time1=omp_get_wtime()
if (tracer.eq.1) Tratio=T0 !!for empty cells
 T1=0.d0
 count_type=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!parallel
Nthreads=1
!$ Nthreads=omp_get_max_threads()
Nshare=ntr/Nthreads
if (tracer.eq.1) then
 allocate(Twork(0:Nthreads-1,-span:nx+span,-span:nz+span))
 Twork=0.d0
end if
allocate(count_work(0:Nthreads-1,0:ntype,-span:nx+span,-span:nz+span))
count_work=0
!$OMP PARALLEL PRIVATE(id,id_min,id_max,r,s,i,k,rank)
rank=0
!$ rank=omp_get_thread_num()
id_min=(1+rank*Nshare)
if (rank.eq.(Nthreads-1)) then
 id_max=max((rank+1)*Nshare,ntr)
else
 id_max=(rank+1)*Nshare
end if
do id=id_min,id_max
 r=rtr_(id)
 s=str_(id)
 i=nint(r/dr,4)
 k=nint(s/ds,4)
 if ((r.le.0.d0).or.(r.ge.aspect).or.(s.le.0.d0).or.(s.ge.1.d0)) then
  write(*,*) "Tracer advected outside of domain - stopping"
  stop
 end if
 if (tracer.eq.1) Twork(rank,i,k)=Twork(rank,i,k)+Ttr_(id)
 count_work(rank,Ttype(id),i,k)=count_work(rank,Ttype(id),i,k)+1
end do
!$OMP END PARALLEL
!$OMP PARALLEL DO PRIVATE(i,k,tt)
do i=0,nx
 do k=0,nz
  if (tracer.eq.1) T1(i,k)=sum(Twork(:,i,k))
  do tt=0,ntype
   count_type(tt,i,k)=sum(count_work(:,tt,i,k))
  end do
 end do
end do
!$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!original
!do id=1,ntr
! r=rtr_(id)
! s=str_(id)
! i=nint(r/dr,4)
! k=nint(s/ds,4)
! if ((r.le.0.d0).or.(r.ge.aspect).or.(s.le.0.d0).or.(s.ge.1.d0)) then
!  write(*,*) "Tracer advected outside of domain - stopping"
!  stop
! end if
! if (tracer.eq.1) T1(i,k)=T1(i,k)+Ttr_(id)
! count_type(Ttype(id),i,k)=count_type(Ttype(id),i,k)+1
!end do
!write(*,*) sum(T1) !!!!test value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end original

  do tt=0,ntype
   do i=1,span
    count_type(tt,-i,0:nz)=count_type(tt,i,0:nz)
    count_type(tt,nx+i,0:nz)=count_type(tt,nx-i,0:nz)
   end do
   do k=1,span
    count_type(tt,0:nx,-k)=count_type(tt,0:nx,k)
    count_type(tt,0:nx,nz+k)=count_type(tt,0:nx,nz-k)
   end do
  end do

nempty=sum(count_type,DIM=1)
!$OMP PARALLEL DO PRIVATE(i,k,tt)
do i=0,nx
 do k=0,nz
  if (nempty(i,k).gt.0) then !!if there are no tracers close to the cell corner use previous value for C
   if (comp.eq.1) then
    do tt=1,ntype !!assume tracer type 0 has zero mass
     if (count_type(tt,i,k).eq.nempty(i,k)) then
      Cnew(tt,i,k)=1.d0
     elseif (count_type(tt,i,k).eq.0) then
      Cnew(tt,i,k)=0.d0
     else
      Cnew(tt,i,k)=real(count_type(tt,i,k),8)/real(nempty(i,k),8)
     end if
    end do
   end if
   if (tracer.eq.1) Tratio(i,k)=T1(i,k)/real(nempty(i,k),8)
  end if
 end do
end do
!$OMP END PARALLEL DO
 if (comp.eq.1) then
  do tt=1,ntype
   call enforceBCs_comp(Cnew(tt,:,:))
  end do
  if (smooth_vis.eq.1) call smoother_time
 end if
 if (tracer.eq.1) then
  call enforceBCs(Tratio)
  if (smooth_vis.eq.1) call smoother_time_T(Tratio)
 end if
 call load_conductivity_array
!$ time2=omp_get_wtime()
!$ tconvert=tconvert+time2-time1
end

subroutine enforceBCs_comp(C0)
 use basics
implicit none
integer*4 :: i
real*8 :: C0(-span:nx+span,-span:nz+span)
do i=1,span
 C0(-i,0:nz)=C0(i,0:nz)      !!symmetry: conservation of mass
 C0(nx+i,0:nz)=C0(nx-i,0:nz)
 C0(:,-i)=C0(:,i)      
 C0(:,nz+i)=C0(:,nz-i)
end do
end 


subroutine velocityBCs
 use basics
 use arrays
implicit none
integer*4 :: i,k
do i=1,span_interp         !!symmetry:free slip
 w(-i,0:nz)=w(i,0:nz)      
 w(nx+i,0:nz)=w(nx-i,0:nz)
 if (Vbc.eq.0) then        !!symmetry:free slip
  u(0:nx,-i)=u(0:nx,i)      
  u(0:nx,nz+i)=u(0:nx,nz-i)
 elseif (Vbc.eq.1) then    !!antisymmetry for rigid boundaries
  u(0:nx,-i)=-u(0:nx,i)      
  u(0:nx,nz+i)=-u(0:nx,nz-i)
 end if
end do
do k=1,span_interp         !!antisymmetry for impermeable boundaries
 w(:,-k)=-w(:,k)           
 w(:,nz+k)=-w(:,nz-k)
 u(-k,:)=-u(k,:)       
 u(nx+k,:)=-u(nx-k,:)
end do
u(0,:)=0.d0                !!impermeable boundaries
u(nx,:)=0.d0
w(:,0)=0.d0
w(:,nz)=0.d0
if (Vbc.eq.1) then !!rigid BCs
 u(:,0)=0.d0
 u(:,nz)=0.d0
end if
end

subroutine interpolate(r,s,DERr,f,finterp)
 use basics
 use arrays
implicit none
integer*4 :: i,j,jj,k,ID
real*8 :: hr,hs,fRinterp(-span_interp:span_interp+1),finterp
real*8 :: a1,r,s
real*8 :: f(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: DERs(1:order-1)       !!holds s derivatives of fRinterp
real*8 :: DERr(1:order-1,0:nx,-span_interp:nz+span_interp)

i=int(r/dr,4)
k=int(s/ds,4) 
!write(*,*) r,s,i,k     
hr=r-real(i,8)*dr

if (hr.le.dro2) then !!forward Taylor expansion
 do jj=-span_interp,span_interp+1
  fRinterp(jj)=f(i,k+jj)
 end do
 do j=1,order-1 !!j=derivative counter
  a1=hr**real(j,8)/factorial(j)
  do jj=-span_interp,span_interp+1
   fRinterp(jj)=fRinterp(jj)+DERr(j,i,k+jj)*a1
  end do
 end do
else                       !!backward Taylor expansion
 hr=dr-hr
 do jj=-span_interp,span_interp+1
  fRinterp(jj)=f(i+1,k+jj)
 end do
 do j=1,order-1 !!j=derivative counter
  a1=(-hr)**real(j,8)/factorial(j)
  do jj=-span_interp,span_interp+1
   fRinterp(jj)=fRinterp(jj)+DERr(j,i+1,k+jj)*a1
  end do
 end do
end if

hs=s-real(k,8)*ds
if (hs.le.dso2) then !!forward Taylor expansion
 do jj=1,order-1
  DERs(jj)=dot_product(D(jj,-span_interp:span_interp),fRinterp(0-span_interp:0+span_interp))/ds_power(jj)
 end do
 finterp=fRinterp(0)
 do j=1,order-1 !!j=derivative counter
  a1=hs**real(j,8)/factorial(j)
  finterp=finterp+DERs(j)*a1
 end do
else                       !!backward Taylor expansion
 hs=ds-hs
 do jj=1,order-1
  DERs(jj)=dot_product(D(jj,-span_interp:span_interp),fRinterp(1-span_interp:1+span_interp))/ds_power(jj)
 end do
 finterp=fRinterp(1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hs)**real(j,8)/factorial(j)
  finterp=finterp+DERs(j)*a1
 end do
end if
end

subroutine compute_derivatives(f,DERr)
!$ use OMP_LIB
 use basics
 use arrays
implicit none
integer*4 :: jj,i,k
real*8 :: f(-span_interp:nx+span_interp,-span_interp:nz+span_interp)
real*8 :: DERr(1:order-1,0:nx,-span_interp:nz+span_interp)
real*8 :: time1,time2
!$ time1=omp_get_wtime()
do jj=1,order-1
 do k=-span_interp,nz+span_interp
  do i=0,nx
   DERr(jj,i,k)=dot_product(D(jj,-span_interp:span_interp),f(i-span_interp:i+span_interp,k))/dr_power(jj)
  end do
 end do
end do
!$ time2=omp_get_wtime()
!$ tcd=tcd+time2-time1
end

subroutine compute_metric_interpolator_derivatives
 use basics
 use arrays
implicit none
integer*4 :: i,k,jj
do jj=1,order-1
 do i=0,nx
  DERxr(jj,i)=dot_product(D(jj,-span_interp:span_interp),xr(i-span_interp:i+span_interp))/dr**real(jj,8)
 end do
end do
do jj=1,order-1
 do k=0,nz
  DERzs(jj,k)=dot_product(D(jj,-span_interp:span_interp),zs(k-span_interp:k+span_interp))/ds**real(jj,8)
 end do
end do
end

subroutine interpolate_xr(r,xr_interp)
 use basics
 use arrays
implicit none
integer*4 :: i,j,jj,ID
real*8 :: hr,xr_interp
real*8 :: a1,r
if (gridX.ne.0) then
i=int(r/dr,4)
hr=r-real(i,8)*dr

if (hr.le.dro2) then !!forward Taylor expansion
 xr_interp=xr(i)
 do j=1,order-1 !!j=derivative counter
  a1=hr**real(j,8)/factorial(j)
  xr_interp=xr_interp+DERxr(j,i)*a1
 end do
else                       !!backward Taylor expansion
 hr=dr-hr
 xr_interp=xr(i+1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hr)**real(j,8)/factorial(j)
  xr_interp=xr_interp+DERxr(j,i+1)*a1
 end do
end if

else
 xr_interp=1.d0
end if
end

subroutine interpolate_zs(s,zs_interp)
 use basics
 use arrays
implicit none
integer*4 :: j,jj,k,ID
real*8 :: hs,zs_interp
real*8 :: a1,s
if (gridZ.ne.0) then
k=int(s/ds,4)
hs=s-real(k,8)*ds

if (hs.le.dso2) then !!forward Taylor expansion
 zs_interp=zs(k)
 do j=1,order-1 !!j=derivative counter
  a1=hs**real(j,8)/factorial(j)
  zs_interp=zs_interp+DERzs(j,k)*a1
 end do
else                       !!backward Taylor expansion
 hs=ds-hs
 zs_interp=zs(k+1)
 do j=1,order-1 !!j=derivative counter
  a1=(-hs)**real(j,8)/factorial(j)
  zs_interp=zs_interp+DERzs(j,k+1)*a1
 end do
end if

else
 zs_interp=1.d0
end if
end


