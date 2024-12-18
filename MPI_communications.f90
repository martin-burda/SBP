MODULE MPI_communications

 use global_data

 implicit none 
 contains 
 !--------------------------------------------------------------------
 !                    SUBROUTINE gather_um
 !--------------------------------------------------------------------  
 SUBROUTINE gather_um
  integer :: rk

#if defined (_MPI)
   if (Nranks>1) then
      call MPI_barrier(MPI_comm_world, ierror)

      do itr=1,Niter 
         cd = (itr-1)*Nranks + rank + 1 

         do rk=1,Nranks-1
            if (root) cd = (itr-1)*Nranks + rk + 1 
            if (cd <= dim) then

               if (rank==rk) then 
                  call MPI_send(buf = um(1:Tn,cd), count = Tn, & 
                  & datatype = MPI_type, dest = 0, tag = 1, comm = MPI_comm_world, & 
                  & ierror = ierror)
               end if
         
               if (rank==0) then                  
                  call MPI_recv(buf = um(1:Tn,cd), count = Tn, &
                  & datatype = MPI_type, source = rk, tag = 1, comm = MPI_comm_world, &
                  & status = status, ierror = ierror)
               end if
            end if

         end do
      end do

   end if
#endif 

 END SUBROUTINE gather_um

 !--------------------------------------------------------------------
 !                    SUBROUTINE scatter_Csim
 !--------------------------------------------------------------------  
 SUBROUTINE scatter_Csim
  integer :: rk

#if defined (_MPI)
   if (Nranks>1) then
      call MPI_barrier(MPI_comm_world, ierror)

      do itr=1,Niter 
         cd = (itr-1)*Nranks + rank + 1 

         do rk=1,Nranks-1
            if (root) cd = (itr-1)*Nranks + rk + 1 
            if (cd <= dim) then 

               if (rank==0) then 
                  call MPI_send(buf = Csim(:,cd), count = Nsim, & 
                  & datatype = MPI_type, dest = rk, tag = 2, comm = MPI_comm_world, & 
                  & ierror = ierror)
               end if
         
               if (rank==rk) then                  
                  call MPI_recv(buf = Csim(:,cd), count = Nsim, &
                  & datatype = MPI_type, source = 0, tag = 2, comm = MPI_comm_world, &
                  & status = status, ierror = ierror)
               end if
            end if

         end do
      end do

   end if
#endif 

 END SUBROUTINE scatter_Csim

 !--------------------------------------------------------------------
 !                    SUBROUTINE gather_Msim
 !--------------------------------------------------------------------  
 SUBROUTINE gather_Msim
   integer :: rk
 
#if defined (_MPI)
    if (Nranks>1) then
       call MPI_barrier(MPI_comm_world, ierror)
 
       do itr=1,Niter 
          cd = (itr-1)*Nranks + rank + 1 
 
          do rk=1,Nranks-1
             if (root) cd = (itr-1)*Nranks + rk + 1 
             if (cd <= dim) then 
 
                if (rank==rk) then 
                   call MPI_send(buf = Msim(:,cd), count = Nsim, & 
                   & datatype = MPI_type, dest = 0, tag = 3, comm = MPI_comm_world, & 
                   & ierror = ierror)
                end if
          
                if (rank==0) then                  
                   call MPI_recv(buf = Msim(:,cd), count = Nsim, &
                   & datatype = MPI_type, source = rk, tag = 3, comm = MPI_comm_world, &
                   & status = status, ierror = ierror)
                end if
             end if
 
          end do
       end do
 
    end if
#endif 
 
  END SUBROUTINE gather_Msim 

END MODULE MPI_communications
