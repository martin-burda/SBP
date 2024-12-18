MODULE arrays

 use global_data
 
 implicit none
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE arrays_allocate
 !-------------------------------------------------------------------- 
 SUBROUTINE arrays_allocate

  allocate(lret(TT,dim), um(TT,dim), ln_um(TT,dim), ln_1_um(TT,dim),  & 
  & Bp(0:smax,TT), hir(0:smax,TT,dim), hi(0:smax,TT,dim),             & 
  & htree(0:smax,TT,dim), htreer(0:smax,TT,dim), hp(0:smax,TT),       &
  & dt(0:smax), pn(smax,TT), pc(smax), ones(TT), sst(TT),             & 
  & nsh(0:smax,TT), vsh(0:smax,TT), v_nsh(0:smax,TT), Ssh(0:smax,TT), &
  & Hsh(smax,TT), pish(0:smax,TT), prsu(0:Smax,TT), prs(0:smax,TT),   &
  & un(TT))

  allocate(mPrMn(TT,Kmax), smPrMn(TT,Kmax), munif(TT), Ncdf(TT,TT,dim))

  allocate(lnGam1(0:smax,TT,dim), lnGam2(0:smax,TT,dim),             & 
  & lnGam3(0:smax,TT,dim), Bdn(0:smax,TT,dim),                       &
  & lnGam_tree1(0:smax,TT,dim), lnGam_tree2(0:smax,TT,dim),          &
  & lnGam_tree3(0:smax,TT,dim), lnBe(gsz**2,dim),                    &
  & lnBep(gsz**2), lnBep2(gsz**2), cdens(gsz**2), cdens2(gsz**2))

 END SUBROUTINE arrays_allocate

 !--------------------------------------------------------------------
 !                    SUBROUTINE grid_arrays_allocate
 !-------------------------------------------------------------------- 
 SUBROUTINE grid_arrays_allocate

  allocate(gr(gsz**2,2), grdat(gsz**2,dim), ln_grdat(gsz**2,dim), & 
  & ln_1_grdat(gsz**2,dim))

 END SUBROUTINE grid_arrays_allocate

 !--------------------------------------------------------------------
 !                    SUBROUTINE grid_arrays_deallocate
 !-------------------------------------------------------------------- 
 SUBROUTINE grid_arrays_deallocate

  deallocate(gr, grdat, ln_grdat, ln_1_grdat)
 
 END SUBROUTINE grid_arrays_deallocate 

 !--------------------------------------------------------------------
 !                    SUBROUTINE Arrays_deallocate
 !-------------------------------------------------------------------- 
 SUBROUTINE arrays_deallocate

  deallocate(lret, um, ln_um, ln_1_um, Bp, hir, hi, htree, htreer,  & 
  & hp, dt, pn, pc, ones, sst, nsh, vsh, v_nsh, Ssh, Hsh, pish,     &
  & prsu, prs, un)

  deallocate(mPrMn, smPrMn, munif)

  deallocate(lnGam1, lnGam2, lnGam3, Bdn, lnGam_tree1, lnGam_tree2, & 
  & lnGam_tree3, lnBe, lnBep, lnBep2, cdens, cdens2)

 END SUBROUTINE arrays_deallocate

END MODULE arrays
