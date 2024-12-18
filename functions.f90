MODULE functions

 use global_data
 use random
 use sorts
 use normals

 implicit none
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE set_up_grid
 !--------------------------------------------------------------------   
 SUBROUTINE set_up_grid
  real(wp) :: lb1, ub1, lb2, ub2, stp1, stp2
  integer  :: i, j, ct, d

   lb1 = 0.0001_wp; ub1 = 0.9999_wp
   lb2 = 0.0001_wp; ub2 = 0.9999_wp

   !--Set up gr
   stp1 = (ub1-lb1)/real(gsz-1) 
   stp2 = (ub2-lb2)/real(gsz-1) 

   ct = 0
   do i=1,gsz
      do j=1,gsz
         ct = ct + 1
         gr(ct,1) = lb1 + real((j-1))*stp1
         gr(ct,2) = lb2 + real((i-1))*stp2
      end do
   end do

   !--Set up grdat
   do d=1,dim
      grdat(:,d) = sum(um(:,d))/real(TT) ! average value
   end do
   grdat(:,d1) = gr(:,1) ! grid variables
   grdat(:,d2) = gr(:,2)
   ln_grdat = log(grdat)
   ln_1_grdat = log(1.-grdat)

  END SUBROUTINE set_up_grid

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

   dt = 0; hp = 0 
   s = 0 
   dt(s) = 1
   htree(s,dt(s),:) = 1 
   hp(0,:) = 1 

   do s=1,smax
      do j=1,Tn
         newnode = .true.
         treeloop: do k=1,dt(s)
            if (all(hi(s,j,:)==htree(s,k,:))) then 
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
         pn(s,hp(s,j)) = hp(s-1,j)
      end do
   end do

   do s=1,smax
      pc(s) = maxval(pn(s,:))
   end do

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
         Bp(s,j) = sum(Bdn(s,j,:)) 
      end do                       
   end do

  END SUBROUTINE eval_Bernstein

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
 !                    SUBROUTINE update_sst
 !--------------------------------------------------------------------   
 SUBROUTINE update_sst
  integer :: s, j
 
   sst = smax/2 

   do concurrent(s=0:smax, j=1:Tn)
      prsu(s,j) = Bp(s,j) + log(pish(s,hp(s,j))) 
   end do

   prsu(:,1:Tn) = exp(prsu(:,1:Tn)) 

   do j=1,Tn
      prs(:,j) = prsu(:,j)/sum(prsu(:,j)) 
   end do

   do j=1,Tn
      do s=smax,1,-1 
         prs(s,j) = sum(prs(0:s,j)) 
      end do
   end do

   call random_number(un(1:Tn))

   do j=1,Tn
      prs(:,j) = prs(:,j) - un(j)
   end do
   where (prs<0._wp) prs = 1._wp

   do j=1,Tn
      sst(j:j) = minloc(prs(:,j)) - 1 
   end do
 
  END SUBROUTINE update_sst
 !--------------------------------------------------------------------
 !                    SUBROUTINE update_nv
 !--------------------------------------------------------------------   
 SUBROUTINE update_nv
  integer :: j, d, s, k

   !--Update nsh
   nsh = 0._wp
   do j=1,Tn
      s = sst(j)
      nsh(s, hp(s,j)) = nsh(s, hp(s,j)) + 1._wp
   end do

   !--Update vsh
   vsh = 0._wp
   do j=1,Tn
      do s=0,sst(j)
         vsh(s, hp(s,j)) = vsh(s, hp(s,j)) + 1._wp
      end do
   end do

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
         Ssh(s,k) = random_beta(1._wp + nsh(s,k), &
                               & prior_a + v_nsh(s,k), .true.) 
      end do
   end do
   Ssh(0,1) = 0._wp 
   Ssh(smax,1:dt(smax)) = 1._wp

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
         nps = sum(ones(1:Tn), mask=active) 

         if (nps==1) then 
            where(active(1:dt(s))) Hsh(s,1:dt(s)) = 1.       
         else
            params(1:dt(s)) = pack(vsh(s,1:dt(s)), mask=active(1:dt(s)))
            params(1:nps) = params(1:nps) + prior_b
            Ddraws(1:nps) = random_dirichlet(nps, params(1:nps))
            Hsh(s,1:dt(s)) = unpack(Ddraws(1:nps), mask=active(1:dt(s)), field=Hsh(s,1:dt(s)))
         end if

      end do
   end do
        
  END SUBROUTINE update_H

 !--------------------------------------------------------------------
 !                    SUBROUTINE update_pi
 !--------------------------------------------------------------------   
 SUBROUTINE update_pi
  integer  :: s, k
  real(wp) :: spsh

   pish = 0._wp 

   pish(0,1) = Ssh(0,1) 

   !--Only the product of the path
   s = 1 
   do k=1,dt(s)
      pish(s,k) = (1._wp - Ssh(s-1,pn(s,k)))*Hsh(s,k)
   end do

   do s=2,smax 
      do k=1,dt(s)
         pish(s,k) = pish(s-1,pn(s,k))*(1._wp - Ssh(s-1,pn(s,k)))*Hsh(s,k)
      end do
   end do

   !--Multiply by Ssh
   do s=1,smax
      do k=1,dt(s)
         pish(s,k) = pish(s,k)*Ssh(s,k)
      end do
   end do

 END SUBROUTINE update_pi

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_cop_dens
 !--------------------------------------------------------------------   
 SUBROUTINE eval_cop_dens
  integer :: s, k, d, i

   cdens = 0._wp
   cdens2 = 0._wp   

   do s=0,smax
      do k=1,dt(s)

         do d=1,dim  
            do i=1,gsz**2
               lnBe(i,d) = lnGam_tree1(s,k,d) - & 
                         & lnGam_tree2(s,k,d) - lnGam_tree3(s,k,d) + &
                         & (htree(s,k,d) - 1._wp)*ln_grdat(i,d) + &
                         & (ttpsr(s) - htree(s,k,d))*ln_1_grdat(i,d)
            end do
         end do
         do i=1,gsz**2
            lnBep(i)  = sum(lnBe(i,:)) 
            lnBep2(i) = lnBe(i,d1) + lnBe(i,d2)
         end do
         cdens  = cdens + exp(lnBep + log(pish(s,k)))
         cdens2 = cdens2 + exp(lnBep2 + log(pish(s,k)))
      end do           
   end do

  END SUBROUTINE eval_cop_dens

 !--------------------------------------------------------------------
 !                    SUBROUTINE copula_simulate
 !--------------------------------------------------------------------   
 SUBROUTINE copula_simulate
  real(wp) :: uni, a, b
  integer  :: s, k, i, j, q, d

   !--Set up vector of multinomial probabilities
   i = 0
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
            NCd(i) = NCd(i) + 1 
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

            do d=1,dim
               a = real(htree(s,k,d))
               b = ttpsr(s) - a + 1._wp
               Csim(q,d) = random_beta(a, b, .true.)  
            end do

         end do      
      end do
   end do   

  END SUBROUTINE copula_simulate

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_lprsim
 !--------------------------------------------------------------------   
 SUBROUTINE eval_lprsim
  integer :: i
   
   Msim = exp(Msim) 
   do i=1,Nsim ! AOR portfolio weights
      lprsim(i) = 0.345_wp*Msim(i,1) + 0.331_wp*Msim(i,2) + &
                & 0.17_wp*Msim(i,3) + 0.06_wp*Msim(i,4) +   &
                & 0.058_wp*Msim(i,5) + 0.02_wp*Msim(i,6) +  & 
                & 0.009_wp*Msim(i,7)
   end do
   lprsim = log(lprsim) 

  END SUBROUTINE eval_lprsim   

END MODULE functions
