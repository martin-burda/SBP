MODULE global_data

#if defined (_MPI)
 use MPI_f08
#endif
 implicit none

 integer, parameter :: dp = selected_real_kind(p=12,r=60) 
 integer, parameter :: sp = selected_real_kind(p=6,r=37)  
 integer, parameter :: wp = dp

 !--MPI info
 integer :: Nranks, rank, ierror
 real(wp):: Nranks_wp
 logical :: root
#if defined (_MPI) 
 type(MPI_status) :: status 
 type(MPI_Datatype) :: MPI_type  
#endif 

 !--Constants
 character(len=2), parameter :: dirr = './' 
 character(len=2), parameter :: dirw = './' 
 real(wp), parameter :: pi = 3.1415926535897_wp
 real(wp), parameter :: srtp_inv = 1._wp/sqrt(2._wp*pi)
 integer,  parameter :: TT = 606 
 integer,  parameter :: T1 = TT-252 
 integer             :: Tn       
 real(wp)            :: Tn_wp
 integer,  parameter :: dim = 7 
 real(wp), parameter :: dim_wp = real(dim)
 integer,  parameter :: smax = 10
 integer,  parameter :: MC = 1100
 integer,  parameter :: burnin = 100
 real(wp), parameter :: prior_a = 10._wp
 real(wp), parameter :: prior_b = 1._wp

 integer  :: ttps(0:smax)  
 real(wp) :: ttpsr(0:smax) 

 real(wp), allocatable :: lret(:,:), um(:,:), ln_um(:,:), ln_1_um(:,:), & 
 & Bp(:,:), hir(:,:,:), htreer(:,:,:), nsh(:,:), vsh(:,:), v_nsh(:,:),  & 
 & Ssh(:,:), Hsh(:,:), pish(:,:), prsu(:,:), prs(:,:), un(:),           &
 & lnGam1(:,:,:), lnGam2(:,:,:), lnGam3(:,:,:), Bdn(:,:,:),             &
 & lnGam_tree1(:,:,:), lnGam_tree2(:,:,:), lnGam_tree3(:,:,:),          &
 & lnBe(:,:), lnBep(:), lnBep2(:)

 integer,  allocatable :: hi(:,:,:), htree(:,:,:), hp(:,:), dt(:),      &
 & pn(:,:), pc(:), ones(:), sst(:) 
 integer :: Niter, itr, cd

 !--Density over 2D grid
 logical, parameter    :: eval_grd = .false.
 integer, parameter    :: d1 = 6 
 integer, parameter    :: d2 = 7 
 integer, parameter    :: gsz = 99 
 real(wp), allocatable :: gr(:,:), grdat(:,:), ln_grdat(:,:), ln_1_grdat(:,:)
 real(wp), allocatable :: cdens(:), cdens2(:)

 !--For marginals sampling
 integer, parameter :: Pdim = 1   
 integer, parameter :: Qdim = 1   
 integer, parameter :: init = min(Pdim,Qdim) 
 integer, parameter :: Kmax = TT  
 real(wp), parameter:: prior_lama = 1._wp
 real(wp), parameter:: prior_lamb = 1._wp 
 real(wp), parameter:: ones_Kmax(Kmax) = 1._wp
 integer  :: zt(TT,dim)
 real(wp) :: nu(dim), eta(Kmax,dim), evec(Kmax,dim), rho(Kmax,dim), ut(TT,dim), & 
 & xi(TT,dim), mu(dim), lambda(Kmax,dim), lambda_inv(Kmax,dim), &
 & omega(dim), alpha(Pdim,dim), beta(Qdim,dim), ns(Kmax), ht(TT,dim), &
 & ac_mu(dim), tm_mu(dim), ac_ab(dim), tm_ab(dim), acr_mu(dim), acr_ab(dim), &
 & lret_mu(TT,dim), dtheta(1+Pdim+Qdim)
 real(wp), allocatable :: mPrMn(:,:), smPrMn(:,:), munif(:), Ncdf(:,:,:)

 !--For simulations for VaR and ES
 logical, parameter    :: write_sim = .false. 
 integer, parameter    :: NSim = 100
 real(wp), parameter   :: NSim_wp = real(NSim)
 real(wp), allocatable :: PrMn(:), sPrMn(:)
 real(wp) :: Csim(Nsim,dim), Msim(Nsim,dim)
 real(wp) :: VaR(MC-burnin), ES(MC-burnin)
 real(sp) :: lprsim(Nsim)
 integer, allocatable  :: NCd(:) 
 integer  :: Anodes 

 !--Index of the VaR quantile
 real(wp), parameter :: alp = 0.05_wp
 real(wp), parameter :: qid_wp = NSim_wp*alp 
 integer, parameter  :: qid = nint(qid_wp) 

 END MODULE global_data
