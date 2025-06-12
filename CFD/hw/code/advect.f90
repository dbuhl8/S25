!Author: Dante Buhl
! CFD HW 4

program advect
  use LinAl
  use NumDE
  use CFD

  implicit none

  ! Global vars
  integer :: i, j, fn

  ! Simvars
  real(kind=kr) :: sf = 0.9
  real(kind=kr) :: xstart=0, xstop=1, tstart=0, tstop=0.3, dx, dt
  integer :: Nx, Nt, ic, ngc = 2, lv1, lv2, file, counter = 20
  character :: method, submethod, bc_type
  character, dimension(7) :: method_array = (/'f','p','p','p','p','p','p'/)
  character, dimension(7) :: submethod_array = (/'u','u', 'd', 'c', 'm', 'o', 'v'/)
  real(kind=kr), allocatable :: X(:), U(:,:), cell_bound(:), U_grid(:,:), t(:,:)

  ! in order to run all (77 base cases this needs to be looped efficiently)
  ! i.e. vary the method and submethod variables from a masked array and then
  ! rerun the following code block in each loop
 
  ! global vars 
  Nx = 32
  Nt = 1000

  print *, " ----------------------------------------------------------- "
  print *, " "
  
  ! method
  do lv1 = 1, 11
    !submethod
    do lv2 = 1, 7
      submethod = submethod_array(lv2)
      method = method_array(lv2)

      ! run code with selected submethod
      allocate(U(Nx+2*ngc,Nt), U_grid(Nx+1+2*ngc,Nt), t(1,Nt))
      t(1,1) = tstart

      ! this specifies which IC to use (1-11)
      ! 1 corresponds to IC 7.38
      ! 11 corresponds to IC 7.48
      ic = lv1
      if (lv1 .eq. 11) then
        bc_type = 'p'
      else 
        bc_type = 'o'
      end if
      
      call grid_init(xstart, xstop, nx, ngc, x, cell_bound, dx)
      call advect_init(U(:,1),U_grid(:,1),cell_bound,nx,ngc,dx,ic)
      call bc(U(:,1),nx,ngc,bc_type)
      call cfl(U(:,1),dx,sf,dt)

      do i = 1, Nt-1
        U(ngc+1:nx+ngc,i+1) = advect_update(U(:,i),nx,ngc,dx,dt,method,submethod)
        t(1,i+1) = t(1,i) + dt
        call bc(U(:,i+1),nx,ngc,bc_type)
        call cfl(U(:,i+1),dx,sf,dt)
        if (mod(i,20) .eq. 0) then
          print *, 'Completed Timestep :', i,', current time: ', t(1,i+1)
        end if
        if (t(1,i+1) .ge. tstop) then
          exit
        end if
      end do
  
      fn = counter
      open(fn,file="U_IC"//trim(str(ic))//"Method"//method//submethod//".dat")

      ! open output file, write to file
      call write_data(U,nx+2*ngc,nt,fn)

      close(fn)
      counter = counter + 1
      fn = counter
      open(fn,file="t_IC"//trim(str(ic))//"Method"//method//submethod//".dat")

      ! open output file, write to file
      call write_data(t,1,nt,fn)

      close(fn)
      counter = counter + 1

      ! deallocate memory for next loop
      deallocate(U, U_grid, x, cell_bound,t)
      print *, " "
      print *, " ----------------------------------------------------------- "
      print *, " "
    end do 
  end do 
  

  contains
    include 'advect_update.f90'  
end program advect
