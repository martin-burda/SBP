MODULE rootfindmod

 use global_data
 use marginals

 implicit none
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE rootfind
 !--------------------------------------------------------------------
 !--Adapted from Numerical Recipes 
 !--Using Brent's method, find the root of a function func
 !--known to lie between x1 and x2
 SUBROUTINE rootfind(func, x1, x2, tol, zbrent, y, cd, i)
  real(wp), intent(in) :: x1, x2, tol, y
  integer, intent(in)  :: cd, i
  real(wp), intent(out):: zbrent
  interface
  function func(x, y, cd)
  use global_data
  implicit none
  real(wp), intent(in) :: x, y
  integer, intent(in)  :: cd
  real(wp)             :: func
  end function func
  end interface
  integer, parameter  :: itmax = 1000
  real(wp), parameter :: eps = epsilon(x1)
  real(wp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
  integer  :: iter   
  logical, parameter  :: debug_flag = .false.

   a  = x1
   b  = x2
   fa = func(a, y, cd)
   fb = func(b, y, cd)

   if ((fa>0._wp .and. fb>0._wp) .or. (fa<0._wp .and. fb<0._wp)) then
      write(*,*)'root must be bracketed for zbrent'
      STOP
   end if

   c  = b
   fc = fb
   do iter=1,itmax
      if ((fb>0.0 .and. fc>0.0) .or. (fb<0.0 .and. fc< 0.0)) then
         c  = a
         fc = fa
         d  = b - a
         e  = d
      end if
      if (abs(fc) < abs(fb)) then
         a  = b
         b  = c
         c  = a
         fa = fb
         fb = fc
         fc = fa
      end if

      tol1 = 2._wp*eps*abs(b) + 0.5_wp*tol
      xm = 0.5_wp*(c-b)

      if (abs(xm) <= tol1 .or. fb == 0._wp) then
         zbrent = b
         return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
         s = fb/fa
         if (a == c) then
            p = 2._wp*xm*s
            q = 1._wp - s
         else
            q = fa/fc
            r = fb/fc
            p = s*(2._wp*xm*q*(q-r)-(b-a)*(r-1._wp))
            q = (q-1._wp)*(r-1._wp)*(s-1._wp)
         end if
         if (p > 0.0) q=-q
         p=abs(p)
         if (2._wp*p<min(3._wp*xm*q-abs(tol1*q),abs(e*q))) then
            e = d
            d = p/q
         else
            d = xm
            e = d
         end if
      else
         d = xm
         e = d
      end if
      a  = b
      fa = fb
      b  = b + merge(d,sign(tol1,xm), abs(d) > tol1)
      fb = func(b, y, cd)

      if (debug_flag) then
         write(*,*)'arg in brent =', b
         write(*,*)'fn  in brent =', fb
      end if

   end do

   write(*,*)'zbrent: exceeded maximum iterations'
   zbrent = b

   END SUBROUTINE rootfind 

END MODULE rootfindmod
