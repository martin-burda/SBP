MODULE timers

 use global_data

 implicit none
  real     :: start_time, finish_time, time_cpu ! for function cpu_time
  real(dp) :: time_wallclock                    ! for function system_clock
  integer  :: start_time2, finish_time2, crate  ! for function system_clock
  real(wp) :: MPI_time_wallclock                ! for function MPI_Wtime
  real(wp) :: start_time3, finish_time3         ! for function MPI_Wtime
 contains
 !--------------------------------------------------------------------
 !                    SUBROUTINE timing_start
 !--------------------------------------------------------------------  
 SUBROUTINE timing_start

   call cpu_time(start_time)
   call system_clock(count_rate=crate)
   call system_clock(count=start_time2)
#if defined (_MPI)   
   start_time3 = MPI_Wtime() 
#endif      

 END SUBROUTINE timing_start

 !--------------------------------------------------------------------
 !                    SUBROUTINE timing_finish
 !--------------------------------------------------------------------  
 SUBROUTINE timing_finish

   call cpu_time(finish_time)
   call system_clock(count=finish_time2)
#if defined (_MPI)   
   finish_time3 = MPI_Wtime()
#endif     
   time_cpu = finish_time - start_time ! in seconds
   time_wallclock = real(finish_time2 - start_time2)/real(crate) ! in seconds
#if defined (_MPI)   
   MPI_time_wallclock = finish_time3 - start_time3
#endif 

 END SUBROUTINE timing_finish

 !--------------------------------------------------------------------
 !                    SUBROUTINE timing_report
 !--------------------------------------------------------------------  
 SUBROUTINE timing_report

   write(*,'(A20,F10.3)') "time_cpu =          ", time_cpu
   write(*,'(A20,F10.3)') "time_wallclock =    ", time_wallclock
#if defined (_MPI)    
   write(*,'(A20,ES15.5)') "MPI_time_wallclock =", MPI_time_wallclock
#endif 

 END SUBROUTINE timing_report

 END MODULE timers
