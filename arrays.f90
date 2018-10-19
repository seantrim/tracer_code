module arrays
real*8, allocatable :: Matrix(:,:),Bvec(:),Row(:),Col(:),MatrixB(:,:),RowB(:),ColB(:)  !!p,q indexing
integer*4, allocatable :: IPIV(:),IPIVB(:),ipvec(:),kpvec(:),parray(:,:),nempty(:,:),Ttype(:),count_type(:,:,:)
real*8, allocatable :: xg(:),zg(:),T(:,:),u(:,:),w(:,:),SF(:,:),vis(:,:) !!i,k indexing
real*8, allocatable :: xr(:),xrr(:),xrrr(:),xrrrr(:),zs(:),zss(:),zsss(:),zssss(:)
real*8, allocatable :: D1(:),D2(:),D3(:),D4(:),Drs(:,:),Drss(:,:),Drrs(:,:),Drrss(:,:),dt_array(:,:),dts_array(:,:)
real*8, allocatable :: factorial(:),D(:,:),dr_power(:),ds_power(:)
real*8, allocatable :: rtr(:),str(:),Ttr(:),DERr_u(:,:,:),DERr_w(:,:,:),DERxr(:,:),DERzs(:,:)
real*8, allocatable :: rtr0(:),rtr1(:),rtr2(:),str0(:),str1(:),str2(:),Cvis(:,:,:),Cbuoy(:,:,:),Tvis(:,:),Tbuoy(:,:)
real*8, allocatable :: Tratio(:,:),Ttr0(:),Ttr1(:),Ttr2(:),DERr_tracer(:,:,:),tracer_space_array(:,:),Textend(:,:)
real*8, allocatable :: DERr_Textend(:,:,:),Cnew(:,:,:),RaC(:),visC(:),mass(:)
real*8, allocatable :: dt_tr(:),dt_stream(:,:,:),dt_vis(:,:,:),dt_fine(:,:)
real*8, allocatable :: conduct(:,:),conduct_r(:,:),conduct_s(:,:),conduct_factor(:),Htype(:)
real*8, allocatable :: strain(:,:),RHS(:,:,:),residual(:,:),error(:,:,:)
real*8, allocatable :: vis_x(:,:,:),vis_z(:,:,:),vis_xx(:,:,:),vis_zz(:,:,:),vis_xz(:,:,:),Matrix_coarse(:,:),B_coarse(:)
real*8, allocatable :: vis_grid(:,:,:),vis_gridf(:,:,:)
real*8, allocatable :: vis_xf(:,:,:),vis_zf(:,:,:),vis_xxf(:,:,:),vis_zzf(:,:,:),vis_xzf(:,:,:)
real*8, allocatable :: ratio(:,:,:)
end module arrays
