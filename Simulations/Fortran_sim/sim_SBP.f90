! cd "/mnt/c/Data/copula2/fortran_sim"
! cd "/scratch/burda/copula2/fortran_sim"

! nvfortran random.f90 sim_SBP.f90 -o sim_SBP.exe -Mbounds
! nvfortran random.f90 sim_SBP.f90 -o sim_SBP.exe -O3 
! ./sim_SBP.exe

!====================================================================
!                    MODULE global_data
!====================================================================
MODULE global_data
 implicit none

 integer, parameter :: dp = selected_real_kind(p=12,r=60) 
 integer, parameter :: sp = selected_real_kind(p=6,r=37)  
 integer, parameter :: wp = dp

 !--WSL2
 character(len=29), parameter :: dirr = '/mnt/c/Data/copula2/data_sim/' 
 character(len=31), parameter :: dirw = '/mnt/c/Data/copula2/output_sim/'
 !--Trillium
 !character(len=31), parameter :: dirr = '/scratch/burda/copula2/data_sim/'
 !character(len=33), parameter :: dirw = '/scratch/burda/copula2/output_sim/'

 !--Constants
 integer,  parameter :: TT = 10000 ! data size @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 integer             :: Tn = TT 
 integer,  parameter :: dim = 3 ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 integer,  parameter :: smax = 6 !6 ! max tree depth

 integer,  parameter :: MC = 11000 !!!!!!!!!
 integer,  parameter :: burnin = 1000

 real(wp), parameter :: prior_a = real(smax/2) ! 1._wp !!!!!!!
 real(wp), parameter :: prior_b = 1._wp

 integer, parameter  :: NSim = TT
 real(wp), parameter :: NSim_wp = real(NSim)

 !--Variables
 integer  :: ttps(0:smax)  ! 2^s (two to the power of s) as an integer 
 real(wp) :: ttpsr(0:smax) 
 real(wp) :: Csim(Nsim,dim)

 real(wp), allocatable :: um(:,:), ln_um(:,:), ln_1_um(:,:),            &
 & Bp(:,:), hir(:,:,:), htreer(:,:,:), nsh(:,:), vsh(:,:), v_nsh(:,:),  & 
 & Ssh(:,:), Hsh(:,:), pish(:,:), prsu(:,:), prs(:,:), un(:),           &
 & lnGam1(:,:,:), lnGam2(:,:,:), lnGam3(:,:,:), Bdn(:,:,:),             &
 & spsh(:), PrMn(:), sPrMn(:)

 integer,  allocatable :: hi(:,:,:), htree(:,:,:), hp(:,:), dt(:),      &
 & pn(:,:), pc(:), ones(:), sst(:), NCd(:)
 integer :: Anodes, m, kf, nf
 character(len=8)  :: ks, ns
 character(len=20) :: rname 
 character(len=70) :: file_name

 !--Density over 2D grid
 integer, parameter    :: gsz = 100 ! grid size @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
 real(wp), parameter   :: lb(dim) = 0.01_wp, ub(dim) = 0.99_wp
 real(wp), allocatable :: gr(:,:), ln_gr(:,:), ln_1_gr(:,:),   &
 & lnGam_tree1(:,:,:), lnGam_tree2(:,:,:), lnGam_tree3(:,:,:), &
 & lnBe(:), lnBep(:), lnBep2(:), cdens_SBP(:)
 
END MODULE global_data

!====================================================================
!                    MODULE arrays
!====================================================================
MODULE arrays
 use global_data
 implicit none
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE arrays_allocate
 !-------------------------------------------------------------------- 
 SUBROUTINE arrays_allocate

  allocate(lnGam1(0:smax,TT,dim), lnGam2(0:smax,TT,dim),              &  
  & lnGam3(0:smax,TT,dim), Bdn(0:smax,TT,dim))                        

  allocate(um(TT,dim), ln_um(TT,dim), ln_1_um(TT,dim),                & 
  & Bp(0:smax,TT), hir(0:smax,TT,dim), hi(0:smax,TT,dim),             & 
  & htree(0:smax,TT,dim), htreer(0:smax,TT,dim), hp(0:smax,TT),       &
  & dt(0:smax), pn(smax,TT), pc(smax), ones(TT), sst(TT),             & 
  & nsh(0:smax,TT), vsh(0:smax,TT), v_nsh(0:smax,TT), Ssh(0:smax,TT), &
  & Hsh(smax,TT), pish(0:smax,TT), prsu(0:Smax,TT), prs(0:smax,TT),   &
  & spsh(0:smax), un(TT))

  allocate(gr(gsz**dim,dim), ln_gr(gsz**dim,dim), ln_1_gr(gsz**dim,dim))

  allocate(lnGam_tree1(0:smax,TT,dim), lnGam_tree2(0:smax,TT,dim),    &
  & lnGam_tree3(0:smax,TT,dim), lnBe(gsz**dim), lnBep(gsz**dim),      &
  & cdens_SBP(gsz**dim)) 

 END SUBROUTINE arrays_allocate

 !--------------------------------------------------------------------
 !                    SUBROUTINE arrays_deallocate
 !-------------------------------------------------------------------- 
 SUBROUTINE arrays_deallocate

  deallocate(lnGam1, lnGam2, lnGam3, Bdn)

  deallocate(um, ln_um, ln_1_um, Bp, hir, hi, htree, htreer, hp, dt, pn, & 
  & pc, ones, sst, nsh, vsh, v_nsh, Ssh, Hsh, pish, prsu, prs, spsh, un)

  deallocate(gr, ln_gr, ln_1_gr, lnGam_tree1, lnGam_tree2, lnGam_tree3, &
  & lnBe, lnBep, cdens_SBP)

 END SUBROUTINE arrays_deallocate

END MODULE arrays

!====================================================================
!                    MODULE functions
!====================================================================
MODULE functions
 use global_data; use random
 implicit none
 contains

 !--------------------------------------------------------------------
 !                    SUBROUTINE set_up_tree
 !--------------------------------------------------------------------   
 SUBROUTINE set_up_tree
  integer :: j, d, s, k
  logical :: newnode

   !--Create hi, htree, hp, dt, pn
   do concurrent(s=0:smax, j=1:Tn, d=1:dim)
      hir(s,j,d) = um(j,d)*ttpsr(s)
   end do
   hi(:,1:Tn,:) = ceiling(hir(:,1:Tn,:))
   hir(:,1:Tn,:) = real(hi(:,1:Tn,:))
   hir(:,Tn+1:TT,:) = 0._wp   

   dt = 0; hp = 0 ! initialize
   s = 0 ! level 0
   dt(s) = 1
   htree(s,dt(s),:) = 1 ! only one base node at level 0
   hp(0,:) = 1 ! all units point to the first (and only) tree node

   !--hp is the full path of h(s) for each individual j
   do s=1,smax
      do j=1,Tn
         newnode = .true.
         treeloop: do k=1,dt(s)
            if (all(hi(s,j,:)==htree(s,k,:))) then ! match with an existing node
               hp(s,j) = k
               newnode = .false.
               EXIT treeloop
            end if
         end do treeloop
         if (newnode) then
            dt(s) = dt(s) + 1
            htree(s,dt(s),:) = hi(s,j,:)             
            hp(s,j) = dt(s)            
         end if
      end do
   end do

   do s=0,smax
      htreer(s,1:dt(s),:) = real(htree(s,1:dt(s),:))
   end do 

   pn = 0
   do s=1,smax
      do j=1,Tn
         pn(s,hp(s,j)) = hp(s-1,j) ! parent node h of j at level s
      end do
   end do

   do s=1,smax
      pc(s) = maxval(pn(s,:))
   end do

   sst = smax ! initialize

   ! write(*,*) 'dt:'
   ! write(*,'(11I8)') dt ! dimensions of tree (# of unique nodes at each level s)
   ! write(*,*)

   ! write(*,*) 'pn:'
   ! do s=1,5 !smax-1
   !    write(*,'(40I4)') pn(s,1:dt(s)) 
   ! end do
   ! write(*,*)

   ! write(*,*) 'pc:'
   ! write(*,'(100I4)') pc
   ! write(*,*)

   ! write(*,*) 'hp for j=1:20'
   ! do s=0,smax
   !    write(*,'(40I4)') hp(s,1:20) ! for the first 10 observations
   ! end do
   ! write(*,*)

  END SUBROUTINE set_up_tree

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_Bernstein
 !--------------------------------------------------------------------   
 SUBROUTINE eval_Bernstein
  integer :: j, s, d, k

   do s=0,smax
      lnGam1(s,:,:) = log_gamma(ttpsr(s) + 1._wp)
   end do

   do d=1,dim
      do s=0,smax       
         do j=1,Tn
            lnGam2(s,j,d) = log_gamma(hir(s,j,d))
         end do
      end do
   end do

   do d=1,dim
      do s=0,smax
         do j=1,Tn
            lnGam3(s,j,d) = log_gamma(ttpsr(s) - hir(s,j,d) + 1._wp)
         end do
      end do
   end do

   ln_um = log(um)
   ln_1_um = log(1.-um)

   do d=1,dim
      do s=0,Smax
         do j=1,Tn
            Bdn(s,j,d) = lnGam1(s,j,d) - lnGam2(s,j,d) - lnGam3(s,j,d) + &
                       & (hir(s,j,d) - 1.)*ln_um(j,d) + &
                       & (ttpsr(s) - hir(s,j,d))*ln_1_um(j,d)
         end do
      end do
   end do

   do s=0,smax
      do j=1,Tn
         Bp(s,j) = sum(Bdn(s,j,:)) !Bdn is in logs so we sum 
      end do                       !to get the log of the product
   end do

   ! write(*,*) "Bp:"
   ! do j=1,10
   !    write(*,'(11F8.4)') (Bp(s,j), s=0,smax)
   ! end do
   ! STOP

  END SUBROUTINE eval_Bernstein

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_sst
 !--------------------------------------------------------------------   
 SUBROUTINE update_sst
  integer :: s, j
 
   do concurrent(s=1:smax, j=1:Tn) ! Bp is in logs
      prsu(s,j) = Bp(s,j) + log(pish(s,hp(s,j))) 
   end do

   prsu(:,1:Tn) = exp(prsu(:,1:Tn)) ! Pr(s) unnormalized
   prsu(0,:) = 0. ! level 0, vacuous parameter

   do j=1,Tn
      prs(:,j) = prsu(:,j)/sum(prsu(:,j)) ! Pr(s) normalized
   end do

! !if (m>100) then
! do j=1,Tn
!    write(*,*) j
!    write(*,'(11F8.3)') Bp(:,j) 
!    !write(*,'(11F8.3)') (log(pish(s,hp(s,j))), s=0,smax)
!    !write(*,'(11F8.3)') - log(spsh(:)) 
!    !write(*,'(11F8.3)') (log(pish(s,hp(s,j)))-log(spsh(s)), s=0,smax)   
!    write(*,'(11F8.3)') (pish(s,hp(s,j)), s=0,smax)     
!    write(*,'(11F8.3)') (spsh(s), s=0,smax)     
!    write(*,'(11F8.3)') (pish(s,hp(s,j))/spsh(s), s=0,smax)   
!    write(*,'(11F8.3)') prsu(:,j)
!    write(*,'(11F8.3)') prs(:,j)   
! read(*,*)    
! end do
! !end if

   do j=1,Tn
      do s=smax,1,-1 ! go backwards, otherwise the sum cumulates the already cumulated
         prs(s,j) = sum(prs(0:s,j)) ! cumulative Pr(s)
      end do
   end do

   call random_number(un(1:Tn))

   do j=1,Tn
      prs(:,j) = prs(:,j) - un(j)
   end do
   where (prs<0._wp) prs = 1._wp

   do j=1,Tn
      sst(j:j) = minloc(prs(:,j)) - 1 ! include level 0, minloc picks position in 0:Smax 
   end do

   !write(*,*) 'sst:'
   !write(*,'(40I4)') sst
   !read(*,*)
 
  END SUBROUTINE update_sst
 !--------------------------------------------------------------------
 !                    SUBROUTINE update_nv
 !--------------------------------------------------------------------   
 SUBROUTINE update_nv
  integer :: j, d, s, k

   !--Update nsh
   nsh = 0._wp
   do j=1,Tn
      s = sst(j) ! only the currently allocated s, i.e. sstar
      nsh(s, hp(s,j)) = nsh(s, hp(s,j)) + 1._wp
   end do

   ! write(*,*) 'nsh for k=1:20'
   ! do s=0,smax
   !    write(*,'(20F8.4)') nsh(s,1:20) ! for the first 20 columns
   ! end do
   ! write(*,*)

   ! write(*,*) 'proportion of sample allocated to level s:'
   ! write(*,'(20F8.4)') (sum(nsh(s,:))/dble(Tn), s=1,smax)

   !--Update vsh
   vsh = 0._wp
   do j=1,Tn
      do s=0,sst(j)
         vsh(s, hp(s,j)) = vsh(s, hp(s,j)) + 1._wp
      end do
   end do

   !  write(*,*) 'vsh for k=1:20'
   !  do s=0,smax
   !     write(*,'(20F8.4)') vsh(s,1:20) ! for the first 10 columns
   !  end do
   !  write(*,*)

   v_nsh = vsh - nsh

  END SUBROUTINE update_nv
  
 !--------------------------------------------------------------------
 !                    SUBROUTINE update_S
 !--------------------------------------------------------------------   
 SUBROUTINE update_S
  integer :: s, k

   Ssh = 0._wp
   do s=0,smax-1
      do k=1,dt(s)
         Ssh(s,k) = random_beta(1.0 + nsh(s,k), &
                              & prior_a + v_nsh(s,k), .true.) 
      end do
   end do
   Ssh(0,1) = 0._wp ! a-priori no stopping at level s=0 ! ##########################################################   
   Ssh(smax,1:dt(smax)) = 1._wp

   ! write(*,*) "Ssh:"
   ! do s=1,smax-1
   !    write(*,'(20F8.2)') Ssh(s,1:dt(s))
   ! end do   
   ! read(*,*)

   ! write(*,*) "Average Ssh for each level:"
   ! write(*,'(20F8.4)') (sum(Ssh(s,1:dt(s)))/dble(dt(s)), s=1,smax)
   ! read(*,*)

  END SUBROUTINE update_S

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_H
 !--------------------------------------------------------------------   
 SUBROUTINE update_H
  real(wp) :: Ddraws(Tn), params(Tn)
  integer  :: s, k, nps
  logical  :: active(Tn)
      
   Hsh = 0._wp
   do s=1,smax
      do k=1,pc(s)
         active = .false.
         where (pn(s,1:dt(s))==k) active(1:dt(s)) = .true.
         nps = sum(ones(1:Tn), mask=active) !number of active nodes to sample
         if (nps==1) then ! only one child node, no need to sample
            where(active(1:dt(s))) Hsh(s,1:dt(s)) = 1.
            !params = 0.            
         else
            params(1:dt(s)) = pack(vsh(s,1:dt(s)), mask=active(1:dt(s)))
            params(1:nps) = params(1:nps) + prior_b
            !params(nps+1:Tn) = 0.
            Ddraws(1:nps) = random_dirichlet(nps, params(1:nps))
            Hsh(s,1:dt(s)) = unpack(Ddraws(1:nps), mask=active(1:dt(s)), field=Hsh(s,1:dt(s)))
         end if

   ! write(*,*) s, pc(s), k
   ! write(*,*) 'active:'
   ! write(*,'(30L7)') active(1:dt(s))
   ! write(*,'(A5,I4)') 'nps =', nps
   ! write(*,*) 'nsh:'
   ! write(*,'(30F7.1)') nsh(s,1:dt(s))
   ! write(*,*) 'vsh:'
   ! write(*,'(30F7.1)') vsh(s,1:dt(s))
   ! write(*,*) 'pn:'
   ! write(*,'(30I7)') pn(s,1:dt(s))
   ! write(*,*) 'params:'
   ! write(*,'(30F7.1)') params(1:dt(s))
   ! write(*,*) 'Ddraws:'
   ! write(*,'(30F7.4)') Ddraws(1:nps), sum(Ddraws(1:nps))
   ! write(*,*) 'Hsh:'
   ! write(*,'(30F7.4)') Hsh(s,1:dt(s))
   ! write(*,*)

      end do
   end do
        
  END SUBROUTINE update_H
  
 !--------------------------------------------------------------------
 !                    SUBROUTINE update_pi
 !--------------------------------------------------------------------   
 SUBROUTINE update_pi
  integer  :: s, k

   pish = 0._wp ! initialize

   !--Evaluate pish
   s = 1 ! level 1
   do k=1,dt(s)
      pish(s,k) = (1._wp - Ssh(s-1,pn(s,k)))*Hsh(s,k)
   end do

   do s=2,smax ! levels > 1
      do k=1,dt(s)
         pish(s,k) = pish(s-1,pn(s,k))*(1._wp - Ssh(s-1,pn(s,k)))*Hsh(s,k)
         !pish(s,k) = (1._wp - Ssh(s-1,pn(s,k)))*Hsh(s,k) !@@@@@@@@@@@@@
      end do
   end do
   
   do s=1,smax
      do k=1,dt(s)
         pish(s,k) = pish(s,k)*Ssh(s,k)
      end do
   end do

   !--sum of pish for each s, i.e. total probability for level s
   spsh(0) = 1. ! vacuous value
   do s=1,smax
      spsh(s) = sum(pish(s,1:dt(s)))
   end do

   !--Normalize
   pish(1:smax,:) = pish(1:smax,:)/sum(spsh(1:smax))

   do s=1,smax
      spsh(s) = sum(pish(s,1:dt(s)))
   end do
   if (mod(m,MC)==0) write(*,'(20F8.4)') spsh(1:smax)

   !write(*,*) 'spsh:'
   !write(*,'(20F8.4)') spsh(1:smax) 
   !write(*,'(F8.2)') sum(spsh(1:smax)) ! should sum up to 1

   ! write(*,*) 'dt:'
   ! write(*,'(11I4)') dt
   ! write(*,*)

   ! write(*,*) 'pish:'
   ! do s=0,smax
   !    write(*,'(A3,I3)') 's =', s
   !    write(*,'(30F7.4)') pish(s,1:dt(s))
   !    write(*,'(A5,I3,A12,F10.5)') 'level', s, ' sum(pish) =', sum(pish(s,1:dt(s)))
   ! end do
   !  write(*,'(A11,F10.5)') 'sum(pish) =', sum(sum(pish,dim=1))
   !  read(*,*)

 END SUBROUTINE update_pi

 !--------------------------------------------------------------------
 !                    SUBROUTINE cop_sim_SBP
 !--------------------------------------------------------------------   
 SUBROUTINE cop_sim_SBP
  real(wp) :: uni, a, b
  integer  :: s, k, i, j, q, d

   !--Set up vector of multinomial probabilities
   i = 0
   PrMn = 0._wp
   do s=1,smax
      do k=1,dt(s)
         i=i+1
         PrMn(i) = pish(s,k)
      end do
   end do
   PrMn = PrMn/sum(PrMn)

   !--Draw from the multinomial w/ probs PrMn
   NCd = 0
   do i=1,Anodes
      sPrMn(i) = sum(PrMn(1:i)) ! cumulative sum
   end do   

   do j=1,NSim
      call random_number(uni)
      iloop: do i=1,Anodes
         if (uni<sPrMn(i)) then
            NCd(i) = NCd(i) + 1 ! # of copula draws for each node
            EXIT iloop
         end if
      end do iloop
   end do

   !--Simulate
   i = 0
   q = 0
   do s=1,smax
      do k=1,dt(s)
         i = i + 1
         do j=1,NCd(i)
            q = q + 1
 !k is just an ordinal count of the node.
 !at level s, dimension d, h is given by htree(s,k,d)
            do d=1,dim
               a = real(htree(s,k,d))
               b = ttpsr(s) - a + 1._wp
               Csim(q,d) = random_beta(a, b, .true.)  

!if (Csim(q,d)==0._wp) Csim(q,d) = tiny(1._wp) 
! if (Csim(q,d)==0._wp) then
!   write(*,*) 'Csim = 0:'
!   write(*,*) d, s, k, htree(s,k,d), ttpsr(s)
!   write(*,*) a, b
! end if
            end do
         end do
      end do
   end do   

  END SUBROUTINE cop_sim_SBP

 !--------------------------------------------------------------------
 !                    SUBROUTINE set_up_grid
 !--------------------------------------------------------------------   
 SUBROUTINE set_up_grid
  real(wp) :: stp(dim)
  integer :: i, j, k, idx(dim), ct, d  

   do d = 1, dim
       stp(d) = (ub(d) - lb(d)) / real(gsz - 1)
   end do

   idx = 0
   ct = 0

   do i = 0, gsz**dim - 1
      ct = ct + 1      
      k = i
      do d = 1, dim
         idx(d) = mod(k, gsz) ! Decode linear index i into dim indices
         k = k / gsz
      end do
      do d = 1, dim
         gr(ct,d) = lb(d) + real(idx(d)) * stp(d) ! Fill in grid points
      end do
   end do

   ln_gr = log(gr)
   ln_1_gr = log(1.-gr)

   open(11,file=dirw//'gr.out',action='write',status='replace') 
   do i=1,gsz**dim
      write(11,'(10F10.5)') gr(i,:)
   end do
   close(11)

  END SUBROUTINE set_up_grid

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_lnGam_tree
 !--------------------------------------------------------------------   
 SUBROUTINE eval_lnGam_tree
  integer :: s, d, k

   do s=0,smax
      lnGam_tree1(s,1:dt(s),:) = lnGam1(s,1:dt(s),:)
   end do

   do d=1,dim
      do s=0,smax       
         do k=1,dt(s)
            lnGam_tree2(s,k,d) = log_gamma(htreer(s,k,d))
         end do
      end do
   end do

   do d=1,dim
      do s=0,smax
         do k=1,dt(s)
            lnGam_tree3(s,k,d) = log_gamma(ttpsr(s) - htreer(s,k,d) + 1._wp)
         end do
      end do
   end do

  END SUBROUTINE eval_lnGam_tree  
  
 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_cop_dens_SBP
 !--------------------------------------------------------------------   
 SUBROUTINE eval_cop_dens_SBP
  integer :: s, k, d

   cdens_SBP = 0._wp

   do s=0,smax
      do k=1,dt(s)
  
         lnBep = 0._wp
         do d=1,dim            
            lnBe = lnGam_tree1(s,k,d) - & 
                 & lnGam_tree2(s,k,d) - lnGam_tree3(s,k,d) + &
                 & (htree(s,k,d) - 1._wp)*ln_gr(:,d) + &
                 & (ttpsr(s) - htree(s,k,d))*ln_1_gr(:,d)

            lnBep = lnBep + lnBe !lnBe is in logs so we sum across all dim
         end do
         cdens_SBP  = cdens_SBP + exp(lnBep + log(pish(s,k)))

      end do           
   end do

   !write(*,'(20F10.3)') cdens_SBP
   !read(*,*)

  END SUBROUTINE eval_cop_dens_SBP

!  !--------------------------------------------------------------------
!  !                    SUBROUTINE eval_KL_SBP_SN
!  !--------------------------------------------------------------------   
!  SUBROUTINE eval_KL_SBP_SN
!   real(wp) :: eps

!    eps = 1.d-6 ! add small epsilon to avoid division by zero
!    KLrat  = (cdens_SN + eps)/(cdens_SBP + eps)
!    KLvals = log(KLrat)*cdens_SN
!    KL_SN  = sum(KLvals)/dble(gsz**2)

!  END SUBROUTINE eval_KL_SBP_SN   

END MODULE functions

!====================================================================
!                    PROGRAM sim_SBP
!====================================================================
PROGRAM sim_SBP
 use global_data; use arrays; use functions
 implicit none
 integer, allocatable :: seed(:) 
 integer :: j, s, nsd, cop
 logical :: simulate, evaluate

   !--Set seed
   call random_seed(size = nsd)
   allocate(seed(nsd))
   seed = 123456
   call random_seed(put = seed)

   write(*,'(A13,I7)') 'sample size:', TT

   cop = 18
!do cop = 1,10
   do kf = 1,10 ! simulated data batch number
      write(*,'(A14,I3)') 'batch number:', kf

      !--Allocate arrays and runtime variables
      call arrays_allocate
      do s=0,Smax
         ttps(s) = 2**s ! two to the power of s
      end do
      ttpsr = real(ttps)
      Hsh   = 1._wp
      ones  = 1

      !--Read data
      call set_file_name
      call read_cop_dat

      !--Set up
      call set_up_tree
      call eval_Bernstein
      simulate = .false. ! if .true. then set kf loop to 1,1
      evaluate = .true.  

      !--Main loop
      do m=1,MC 

         !if (mod(m,1000)==0) write(*,*) m
         if (m>1) call update_sst
         call update_nv    
         call update_S 
         call update_pi   

      end do

      if (simulate) then ! simulate from copula

         call update_H
         call update_pi   

         !--Obtain # of active nodes in the tree for s>0
         Anodes = sum(dt) - 1 ! exclude level s=0
         !--multinomial and cumulative probabilities for sampling
         allocate(PrMn(Anodes), sPrMn(Anodes), NCd(Anodes)) 
         call cop_sim_SBP
         deallocate(PrMn, sPrMn, NCd)   
         call write_cop_sim_SBP 

      end if

      if (evaluate) then ! evaluate SBP copula on a grid

         call set_up_grid ! for SBP copula      
         call eval_lnGam_tree
         if (.not.(simulate)) call update_H
         call update_pi
         call eval_cop_dens_SBP
         call write_cop_dens_SBP

      end if

      call arrays_deallocate

   end do
!end do

contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE set_file_name
 !--------------------------------------------------------------------   
 SUBROUTINE set_file_name

   !--Choose true copula type
   select case(cop)
   case(1)
      rname = 'Gaussian1'
   case(2)
      rname = 'Gaussian2'
   case(3)
      rname = 'Studentt1'
   case(4)
      rname = 'Studentt2'
   case(5)
      rname = 'Studentt3'
   case(6)
      rname = 'Studentt4'
   case(7)
      rname = 'Gumbel1'     ! dim = 2, case 1 (a=2)
   case(8)
      rname = 'Gumbel2'     ! dim = 2, case 1 (a=2)
   case(9)
      rname = 'mix1'
   case(10)
      rname = 'mix2'
   case(11)
      rname = 'Gumbel-hd31' ! dim = 3, case 1 (a=2)
   case(12)
      rname = 'Gumbel-hd32' ! dim = 3, case 2 (a=4)     
   case(13)
      rname = 'Gumbel-hd41' ! dim = 4, case 1 (a=2)
   case(14)
      rname = 'Gumbel-hd42' ! dim = 4, case 2 (a=4)
   case(15)
      rname = 'mixture-hd31' ! dim = 3, case 1
   case(16)
      rname = 'mixture-hd32' ! dim = 3, case 2
   case(17)
      rname = 'mixture-hd41' ! dim = 4, case 1
   case(18)
      rname = 'mixture-hd42' ! dim = 4, case 2         
   end select

   if (kf==1) write(*,*) rname

   !--Set up file name labels
   nf = TT 
   write(ks,'(I0)') kf
   write(ns,'(I0)') nf

 END SUBROUTINE set_file_name

 !--------------------------------------------------------------------
 !                    SUBROUTINE read_cop_dat
 !--------------------------------------------------------------------   
 SUBROUTINE read_cop_dat
  integer :: i

   file_name = trim(rname)//'_dat_'//trim(ns)//'_'//trim(ks)
   open(10, file=dirr//trim(file_name)//'.csv', action='read')
   do i=1,TT
      read(10,*) um(i,:)
   end do
   close(10)  

 END SUBROUTINE read_cop_dat

 !--------------------------------------------------------------------
 !                    SUBROUTINE write_cop_sim_SBP
 !--------------------------------------------------------------------   
 SUBROUTINE write_cop_sim_SBP
  integer :: i, d

   file_name = trim(rname)//'_csim_'//trim(ns)//'_'//trim(ks)    
   open(11,file=dirw//trim(file_name)//'.out', action='write', &
       & status='replace')
   do i=1,Nsim
      write(11,'(10F15.7)') (Csim(i,d), d=1,dim)
   end do
   close(11)   

 END SUBROUTINE write_cop_sim_SBP 

 !--------------------------------------------------------------------
 !                    SUBROUTINE write_cop_dens_SBP
 !--------------------------------------------------------------------   
 SUBROUTINE write_cop_dens_SBP
  integer :: i, d

   file_name = trim(rname)//'_cdens_'//trim(ns)//'_'//trim(ks)      
   open(12,file=dirw//trim(file_name)//'.out', action='write', &
       & status='replace')
   do i=1,gsz**dim
      write(12,'(F15.5)') cdens_SBP(i)
   end do
   close(12)  

 END SUBROUTINE write_cop_dens_SBP

END PROGRAM sim_SBP
