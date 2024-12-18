MODULE marginals

 use global_data
 use functions
 use random
 
 implicit none
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE sample_eta_ui
 !--------------------------------------------------------------------   
 SUBROUTINE sample_eta_ui(cd)  
  integer, intent(in) :: cd 
  real(wp) :: a, b, c
  integer  :: k, j

   do concurrent (k=1:Kmax) 
      ns(k) = sum(ones_Kmax, mask=(zt(:,cd)==k))
   end do

   c = 0._wp
   do k=1,Kmax
      a = ns(k) + 1._wp
      c = c + ns(k)
      b = real(Tn) - c + nu(cd)
      eta(k,cd) = random_beta(a, b, .true.)
   end do

   rho(1,cd) = eta(1,cd)
   do k=2,Kmax
      evec(k-1,cd) = 1._wp - eta(k-1,cd)
      rho(k,cd) = product(evec(1:k-1,cd))*eta(k,cd)
   end do

   !--draw ut(t,cd) from uniform[0,rho(zt(t,cd))]
   call random_number(ut(1:Tn,cd))
   do concurrent(j=1:Tn)
      ut(j,cd) = ut(j,cd)*rho(zt(j,cd),cd)
   end do 

   !--# of mixture components can't exceed # of observations
   rho(Tn+1:TT,cd) = 0._wp 

 END SUBROUTINE sample_eta_ui

 !--------------------------------------------------------------------
 !                    SUBROUTINE sample_lambda
 !--------------------------------------------------------------------   
 SUBROUTINE sample_lambda(cd)  
  integer, intent(in) :: cd
  real(wp) :: a, b, sxi
  logical  :: pt
  integer  :: t, k, i 

   !--Evaluate ht and xi
   ht(1:init,cd) = sum((lret(1:Tn,cd)-mu(cd))**2)/Tn_wp

   do t=init+1,Tn
      ht(t,cd) = omega(cd) + & 
               & dot_product(alpha(:,cd),(lret(t-Pdim:t-1,cd)-mu(cd))**2) + &
               & dot_product(beta(:,cd),ht(t-Pdim:t-1,cd))
      xi(t,cd) = (lret(t,cd)-mu(cd))/sqrt(ht(t,cd)) 
   end do

   do k=1,Tn
      a = (prior_lama + ns(k))/2._wp
      sxi = sum(xi(:,cd)**2, mask=(zt(:,cd)==k))
      b = (prior_lamb + sxi)/2._wp 
      b = 1._wp/b 
      lambda(k,cd) = random_gamma(a,b) 
   end do

   lambda_inv(:,cd) = 1._wp/lambda(:,cd) 

 END SUBROUTINE sample_lambda

 !--------------------------------------------------------------------
 !                    FUNCTION Norm_pdf
 !--------------------------------------------------------------------   
 PURE FUNCTION Norm_pdf(xin, var)
  !--The mean of xi is zero
  real(wp) :: Norm_pdf
  real(wp), intent(in) :: xin, var

   Norm_pdf = srtp_inv*exp(-0.5_wp*xin**2/var)/sqrt(var)

 END FUNCTION Norm_pdf

 !--------------------------------------------------------------------
 !                    SUBROUTINE sample_zt
 !--------------------------------------------------------------------   
 SUBROUTINE sample_zt(cd)  
  integer, intent(in) :: cd 
  integer :: t, k

   do t=1,Tn
      do k=1,Kmax
         mPrMn(t,k) = Norm_pdf(xi(t,cd), lambda_inv(k,cd))
      end do      
   end do

   do t=1,Tn
      do k=1,Kmax
         if (ut(t,cd)>(rho(k,cd)) ) mPrMn(t,k) = 0._wp
      end do
   end do

   call random_number(munif(1:Tn))

   do t=1,Tn

      smPrMn(t,1) = mPrMn(t,1)
      do k=2,Kmax      
         smPrMn(t,k) = smPrMn(t,k-1) + mPrMn(t,k) 
      end do
      !--reweight them to sum up to 1
      smPrMn(t,:) = smPrMn(t,:)/smPrMn(t,Kmax)

      !--multinomial sampling
      smPrMn(t,:) = smPrMn(t,:) - munif(t)
      where (smPrMn(t,:)<0._wp) smPrMn(t,:) = 1._wp
      zt(t:t,cd) = minloc(smPrMn(t,:))

   end do

   do t=1,Tn
      if (zt(t,cd)==0) then
         write(*,*) 'zt = 0:', t, cd
         read(*,*)
      end if
   end do

 END SUBROUTINE sample_zt   

 !--------------------------------------------------------------------
 !                    FUNCTION ln_mkernel
 !--------------------------------------------------------------------   
 REAL(wp) FUNCTION ln_mkernel(cd)
  integer, intent(in) :: cd
  real(wp):: lnmk(Tn)
  integer :: t, k

   !--Evaluate ht
   ht(1:init,cd) = sum((lret(1:Tn,cd)-mu(cd))**2)/Tn_wp 

   do t=init+1,Tn
      ht(t,cd) = omega(cd) + & 
               & dot_product(alpha(:,cd),(lret(t-Pdim:t-1,cd)-mu(cd))**2) + &
               & dot_product(beta(:,cd),ht(t-Pdim:t-1,cd))
   end do

   do concurrent(t=1:Tn)
      lnmk(t) = -0.5_wp*log(ht(t,cd)/lambda(zt(t,cd),cd)) &
              & -lambda(zt(t,cd),cd)*(lret(t,cd)-mu(cd))**2/(2._wp*ht(t,cd)) 
   end do        

   ln_mkernel = sum(lnmk)

 END FUNCTION ln_mkernel

 !--------------------------------------------------------------------
 !                    SUBROUTINE sample_mu
 !--------------------------------------------------------------------   
 SUBROUTINE sample_mu(cd)
  integer, intent(in) :: cd
  real(wp), parameter :: RWi_mu = 0.001_wp
  real(wp) :: ln_num, ln_den, ln_alp, alp, un
  real(wp) :: rn, mu_c
  integer  :: j, k, p

   ln_den = ln_mkernel(cd) 
   mu_c = mu(cd)
 
   rn = random_normal()   
   mu(cd) = mu(cd) + RWi_mu*rn
   ln_num = ln_mkernel(cd)
   ln_alp = ln_num - ln_den

   if (ln_alp>0._wp) then
      alp = 1._wp
   else
      alp = exp(ln_alp)
   end if

   call random_number(un) 
   if (un < alp) then 
      ac_mu(cd) = ac_mu(cd) + 1._wp
   else 
      mu(cd) = mu_c
   end if

   tm_mu(cd) = tm_mu(cd) + 1._wp
   acr_mu(cd) = ac_mu(cd)/tm_mu(cd)

 END SUBROUTINE sample_mu   

 !--------------------------------------------------------------------
 !                    SUBROUTINE sample_garch
 !--------------------------------------------------------------------   
 SUBROUTINE sample_garch(cd)
  integer, intent(in) :: cd
  real(wp), parameter :: RWi_ab = 0.01_wp
  real(wp) :: ln_num, ln_den, ln_alp, alp, un
  real(wp) :: omega_c, alpha_c(Pdim), beta_c(Qdim)
  real(wp) :: rno, rna(Pdim), rnb(Qdim)
  integer  :: j, k, p

   ln_den  = ln_mkernel(cd)
   omega_c = omega(cd)
   alpha_c = alpha(:,cd)
   beta_c  = beta(:,cd)
 
   rno = random_normal()
   do j=1,Pdim
      rna(j) = random_normal()   
   end do
   do j=1,Qdim
      rnb(j) = random_normal()
   end do

   omega(cd)   = omega(cd) + omega(cd)*RWi_ab*rno   
   alpha(:,cd) = alpha(:,cd) + alpha(:,cd)*RWi_ab*rna 
   beta(:,cd)  = beta(:,cd) + beta(:,cd)*RWi_ab*rnb  
   beta(:,cd)  = abs(beta(:,cd))
   ln_num = ln_mkernel(cd)
   ln_alp = ln_num - ln_den

   if (ln_alp>0._wp) then
      alp = 1._wp
   else
      alp = exp(ln_alp)
   end if

   call random_number(un) 
   if (un < alp) then 
      ac_ab(cd) = ac_ab(cd) + 1._wp
   else 
      omega(cd) = omega_c
      alpha(:,cd) = alpha_c 
      beta(:,cd) = beta_c
   end if

   tm_ab(cd) = tm_ab(cd) + 1._wp
   acr_ab(cd) = ac_ab(cd)/tm_ab(cd)

 END SUBROUTINE sample_garch

 !--------------------------------------------------------------------
 !                    SUBROUTINE eval_mcdf
 !--------------------------------------------------------------------   
 SUBROUTINE eval_mcdf(cd)
  integer, intent(in) :: cd
  integer :: t, k

   lret_mu(:,cd) = lret(:,cd)-mu(cd)

   do concurrent(t=1:Tn,k=1:Tn)
      Ncdf(t,k,cd) = lret_mu(t,cd)/sqrt(ht(t,cd)*lambda_inv(zt(k,cd),cd))
   end do      

   do concurrent(t=1:Tn,k=1:Tn)
      Ncdf(t,k,cd) = normal_01_cdf(real(Ncdf(t,k,cd)))
   end do

   do concurrent(t=1:Tn)
      um(t,cd) = sum(Ncdf(t,:,cd))/Tn_wp
   end do

 END SUBROUTINE eval_mcdf

 !--------------------------------------------------------------------
 !                    FUNCTION mcdf1
 !--------------------------------------------------------------------   
 FUNCTION mcdf1(x, y, cd)
  implicit none
  real(wp), intent(in) :: x, y
  integer, intent(in)  :: cd   
  real(wp)             :: mcdf1
  real(wp) :: xmu, cdfk(Tn)
  integer  :: i, k
  
   xmu = x - mu(cd)

   do concurrent(k=1:Tn)
      cdfk(k) = x/sqrt(ht(k,cd)*lambda_inv(zt(k,cd),cd))
   end do      

   do k=1,Tn
      cdfk(k) = normal_01_cdf(real(cdfk(k)))
   end do

   mcdf1 = sum(cdfk)/Tn_wp - y 

   END FUNCTION mcdf1

 !--------------------------------------------------------------------
 !                    FUNCTION mcdf_pred
 !--------------------------------------------------------------------   
 FUNCTION mcdf_pred(x, y, cd)
  implicit none
  real(wp), intent(in) :: x, y
  integer, intent(in)  :: cd   
  real(wp)             :: mcdf_pred
  real(wp) :: xmu, ht_pred, cdfk(Tn)
  integer  :: i, k
   
   xmu = x - mu(cd)

   !--Predict ht
   ht_pred = omega(cd) + & 
           & dot_product(alpha(:,cd),(lret(Tn-Pdim+1:Tn,cd)-mu(cd))**2) + &
           & dot_product(beta(:,cd),ht(Tn-Pdim+1:Tn,cd))

   do concurrent(k=1:Tn)
      cdfk(k) = x/sqrt(ht_pred*lambda_inv(zt(k,cd),cd))
   end do      

   do k=1,Tn
      cdfk(k) = normal_01_cdf(real(cdfk(k)))
   end do

   mcdf_pred = sum(cdfk)/Tn_wp - y 

 END FUNCTION mcdf_pred

 !--------------------------------------------------------------------
 !                    SUBROUTINE plot_marginals
 !--------------------------------------------------------------------  
   SUBROUTINE plot_marginals
      real(wp) :: pt
      integer :: i
    
       open(19,file=dirw//'sc1_mcdf1.out',action='write',status='replace')
       do cd=1,dim
          pt = -10._wp
          do i=1,200
             write(19,*) cd, pt, mcdf1(pt, 0._wp, cd)
             pt = pt + 0.1_wp
          end do
       end do
       close(19)
    
     END SUBROUTINE plot_marginals     

END MODULE marginals   
 