program SBP

 use global_data
 use arrays
 use timers
 use functions
 use marginals
 use rootfindmod
 use MPI_communications
 use newuoa_mod

 implicit none
 integer, allocatable :: seed(:)  
 real(wp) :: dev ! counter of density evaluations
 real(wp) :: lnk
 integer  :: istat, nsd, i, j, k, s, m, ct, ctTot
 real(wp) :: lret_full(TT,7), fout

   !--Initialize MPI
#if defined (_MPI)
   call MPI_init(ierror) 
   call MPI_comm_size(MPI_comm_world, Nranks, ierror)
   call MPI_comm_rank(MPI_comm_world, rank, ierror)

   if (rank==0) then
      root = .true.
      write(*,'(A16,I3,A6)') 'MPI enabled with', Nranks, ' ranks'
   else
      root = .false.
   end if

   if (wp==sp) then 
      MPI_type = MPI_real
   elseif (wp==dp) then
      MPI_type = MPI_double_precision
   else
      write(*,*) 'Precision neither single nor double, check wp, sp, dp'
   end if

#else
   write(*,'(A15)') 'No MPI detected'
   rank = 0
   root = .true.
   Nranks = 1
#endif  
   Nranks_wp = real(Nranks)

   !--Initialize different seed for each rank
   call random_seed(size = nsd)
   allocate(seed(nsd))
   call random_seed(get = seed)
   seed = (rank+1)*seed
   call random_seed(put = seed)

   !--Allocate arrays 
   call arrays_allocate
   if (eval_grd) call grid_arrays_allocate

   !--Initialize runtime variables
   do s=0,Smax
      ttps(s) = 2**s ! two to the power of s
   end do
   ttpsr = real(ttps)
   ones = 1
   dev = 0._wp
   ac_mu = 0._wp; tm_mu = 0._wp; ac_ab = 0._wp; tm_ab = 0._wp !RW counters
   m = 0

   !--Initialize runtime variables for marginals   
   nu  = 1._wp 
   eta = 1._wp/real(Kmax) 
   rho = 1._wp 
   ut  = 1._wp  
   xi  = 1._wp  

   do k=1,T1 
      zt(k,:) = 1 
   end do
   zt(T1+1:TT,:) = 0 

   Niter = ceiling(dim_wp/Nranks_wp)

   !--Read data
   open(10,file=dirr//'log_returns.txt',action='read') 
   do j=1,TT
      read(10,*) lret_full(j,:) 
   end do
   close(10)

   lret(:,1) = lret_full(:,1)   
   lret(:,2) = lret_full(:,2)
   lret(:,3) = lret_full(:,3)
   lret(:,4) = lret_full(:,4)
   lret(:,5) = lret_full(:,5)   
   lret(:,6) = lret_full(:,6)  
   lret(:,7) = lret_full(:,7)        

   !--Initialize
   if (root) then
      call timing_start       
      if (write_sim) open(16,file=dirw//'sc1_CMsim.out',action='write', &
                         & status='replace')  
      open(17,file=dirw//'sc1_VE.out',action='write',status='replace')    
      close(17)
   end if

   !-----------------------------------------------------------------
   do Tn=T1,TT 
      if (root) write(*,*) 'Tn =', Tn
      Tn_wp = real(Tn)
      if (Tn>T1) zt(Tn,:) = 1 

      !--Each MPI rank samples marginals 
      do itr=1,Niter ! distributed among ranks with MPI, otherwise root
         cd = (itr-1)*Nranks + rank + 1 ! each rank gets its current dimension
         if (cd <= dim) then ! in case Niter*Nranks > dim (with MPI)

            if (Tn==T1) then
               !--Initialize the marginals
               lambda_inv(:,cd) = 1._wp ! variance of scaled mixture components
               lambda(:,cd) = 1._wp ! inverse variance of scaled mixture components
      
               mu(cd) = sum(lret(:,cd))/real(TT) 
      
               omega(cd) = 0.01_wp
               alpha(:,cd) = 0.1_wp
               beta(:,cd) = 0.85_wp        
      
               dtheta(1) = omega(cd)
               dtheta(1+1:1+Pdim) = alpha(:,cd)
               dtheta(1+Pdim+1:1+Pdim+Qdim) = beta(:,cd)
      
               call newuoa(calfun, dtheta, fout, rhobeg=5.0d-2, maxfun=1000000)
      
               omega(cd) = dtheta(1) 
               alpha(:,cd) = dtheta(1+1:1+Pdim)
               beta(:,cd) = dtheta(1+Pdim+1:1+Pdim+Qdim)

            end if

            if (Tn==T1) then 
               ctTot = 1000 ! # of iterations through marginal params
            else
               ctTot = 100
            end if

            do j=1,ctTot

               call sample_eta_ui(cd) 
               call sample_lambda(cd) 
               call sample_zt(cd)
               call sample_mu(cd)

               call sample_garch(cd)       

            end do 
            write(*,'(I4,2F10.5)') cd, acr_mu(cd), acr_ab(cd)
            
         end if

         call eval_mcdf(cd)

      end do ! end of loop over marginals, distributed among ranks if _MPI == true

#if defined (_MPI)      
      call MPI_barrier(MPI_comm_world, ierror)
#endif 

      call gather_um  

      if (root) then

         if (mod(Tn,10)==0) write(*,*) Tn
         call set_up_tree
         call eval_Bernstein         
         call eval_lnGam_tree

      end if         

      !--------------------------------------------------------------
      do m=1,MC ! the ranks have to come here for scatter_Csim below

         if (root) then
            !if (mod(m,100)==0) write(*,*) m

            if (m>1) call update_sst
            call update_nv    
            call update_S 
            call update_H
            call update_pi

            if (eval_grd) then

               !--Write out um
               if ((Tn==T1).and.(m==1)) then 
                  open(18,file=dirw//'sc1_um.out',action='write',status='replace')
                  do i=1,Tn
                     write(18,'(7F10.5)') um(i,:)
                  end do
                  close(18)                  
               end if

               !--Evaluate copula density over a grid
               if ((m>burnin).and.(mod(m,100)==0)) then
                  call set_up_grid
                  call eval_cop_dens
                  dev = dev + 1._wp
               end if
            end if

         end if

         if (m>burnin) then ! the ranks have to come here for later scatter_Csim
            if (root) then

               !--Obtain # of active nodes in the tree for s>0
               Anodes = sum(dt) - 1 ! exclude level s=0
               !--multinomial and cumulative probabilities for sampling
               allocate(PrMn(Anodes), sPrMn(Anodes), NCd(Anodes)) 
               call copula_simulate
               deallocate(PrMn, sPrMn, NCd)

            end if

            !--root sends Csim to ranks (need to account for dim > Nranks or Nranks=1)
            call scatter_Csim ! the ranks have to get here
       
            !--Each rank inverts Csim through marginal = Msim   
            do itr=1,Niter ! distributed among ranks with MPI, otherwise root
               cd = (itr-1)*Nranks + rank + 1 ! each rank gets its current dimension
               if (cd <= dim) then ! in case Niter*Nranks > dim (with MPI)
            
                  do i=1,Nsim
                     call rootfind(mcdf_pred, -10._wp, 10._wp, 0.000001_wp, & 
                                  & Msim(i,cd), Csim(i,cd), cd, i)                 
                  end do

               end if
            end do

            !--ranks send Msim to root
            call gather_Msim
            if ((write_sim).and.(root).and.(Tn==T1)) call write_CMsim

            if (root) then
               call eval_lprsim 
               call quicksort_nr(lprsim)
               VaR(m-burnin) = lprsim(qid)
               ES(m-burnin)  = sum(lprsim(1:qid))/real(qid)
            end if

         end if
      end do

      if (root) then
         open(17,file=dirw//'sc1_VE.out',action='write',position='append')    
         write(17,'(I6,2F15.7)') Tn, sum(VaR)/real(MC-burnin), sum(ES)/real(MC-burnin)
         close(17)
      end if

      if (root) then
         write(*,*) 'pish:'
         do s=0,smax
            write(*,'(A5,I3,A12,F10.5)') 'level', s, ' sum(pish) =', sum(pish(s,1:dt(s)))
         end do
         write(*,'(A11,F10.5)') 'sum(pish) =', sum(sum(pish,dim=1))
      end if

   end do

   if ((eval_grd).and.(root)) then
      cdens = cdens/dev ! take an average
      cdens2 = cdens2/dev
      call write_cop_dens
   end if

   if (root) then
      call timing_finish
      call timing_report
      if (write_sim) close(16)
      close(17)
   end if

   !--Cleanup
   call arrays_deallocate 
   if (eval_grd) call grid_arrays_deallocate

#if defined (_MPI)
   call MPI_barrier(MPI_comm_world, ierror)
   call MPI_finalize(ierror)
#endif 

contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE calfun
 !--------------------------------------------------------------------   
 SUBROUTINE calfun(thin, f)
   implicit none
   real(dp), intent(IN)   :: thin(:)
   real(dp), intent(OUT)  :: f
   integer :: k, q

   omega(cd) = thin(1) 
   alpha(:,cd) = thin(1+1:1+Pdim)
   beta(:,cd) = thin(1+Pdim+1:1+Pdim+Qdim)

   f = -ln_mkernel(cd)

   if (.not.((f<1.e10).and.(f>-1.e10))) then
      f = 10000._wp
   end if

 END SUBROUTINE calfun

 !--------------------------------------------------------------------
 !                    SUBROUTINE write_cop_dens
 !--------------------------------------------------------------------   
 SUBROUTINE write_cop_dens

   open(11,file=dirw//'sc1_cdens.out',action='write',status='replace') 
   ct = 0
   do i=1,gsz
      do j=1,gsz
         ct = ct + 1
         write(11,'(2F10.3,F15.7)') gr(ct,1), gr(ct,2), cdens(ct)
      end do
   end do
   close(11) 
  
   open(11,file=dirw//'sc1_cdens2.out',action='write',status='replace') 
   ct = 0
   do i=1,gsz
      do j=1,gsz
         ct = ct + 1
         write(11,'(2F10.3,F15.7)') gr(ct,1), gr(ct,2), cdens2(ct)
      end do
   end do
   close(11)   

 END SUBROUTINE write_cop_dens

 !--------------------------------------------------------------------
 !                    SUBROUTINE write_CMsim
 !--------------------------------------------------------------------   
 SUBROUTINE write_CMsim
  integer :: i, d

   do i=1,Nsim
      do d=1,dim
         write(16,'(2I6,I3,2F15.7)') m, i, d, Csim(i,d), Msim(i,d)
      end do
   end do

 END SUBROUTINE write_CMsim 

end program SBP
