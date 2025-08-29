subroutine soln_PPM(dt)

#include "definition.h"  

  use grid_data
  use sim_data
  use slopeLimiter
  use eigensystem

  implicit none
  real, intent(IN) :: dt
  integer :: i

  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NSYS_VAR) :: vecL,vecR,sigL,sigR,sig2L,sig2R,vec2L,vec2R
  integer :: kWaveNum
  real :: lambdaDtDx, delC1, delC2
  real, dimension(NUMB_VAR)  :: delV,delL,delR,aLR,delLL,dVL,dVC,dVR
  real, dimension(NUMB_VAR,3) :: C
  real, dimension(NUMB_WAVE) :: delW
  integer :: nVar
  

  ! we need conservative eigenvectors
  conservative = .false.

  do i = gr_ibeg-1, gr_iend+1
    C = 0.
    delLL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i-1)-gr_V(DENS_VAR:PRES_VAR,i-2)
    delL(DENS_VAR:PRES_VAR)  = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
    delR(DENS_VAR:PRES_VAR)  = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
    delRR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+2)-gr_V(DENS_VAR:PRES_VAR,i+1)

    do nVar = DENS_VAR, PRES_VAR
      if (sim_limiter == 'minmod') then
         call minmod(delR(nVar),delRR(nVar),dvR(nVar))
         call minmod(delL(nVar),delR(nVar),dvC(nVar))
         call minmod(delLL(nVar),delL(nVar),dvL(nVar))
      elseif (sim_limiter == 'vanLeer') then
         !call vanLeer(delL(nVar),delR(nVar),delV(nVar))
         call vanLeer(delR(nVar),delRR(nVar),dvR(nVar))
         call vanLeer(delL(nVar),delR(nVar),dvC(nVar))
         call vanLeer(delLL(nVar),delL(nVar),dvL(nVar))

      elseif (sim_limiter == 'mc') then
         !call mc(delL(nVar),delR(nVar),delV(nVar))
         call mc(delR(nVar),delRR(nVar),dvR(nVar))
         call mc(delL(nVar),delR(nVar),dvC(nVar))
         call mc(delLL(nVar),delL(nVar),dvL(nVar))
      endif
    end do 
     
    vecL = 0.5*(gr_V(DENS_VAR:PRES_VAR,i-1)+gr_V(DENS_VAR:PRES_VAR,i)) - (1./6)*&
      (dvC - dvL)
    vecR = 0.5*(gr_V(DENS_VAR:PRES_VAR,i)+gr_V(DENS_VAR:PRES_VAR,i+1)) - (1./6)*&
      (dvR - dvC)

    C(:,3) = (6./(gr_dx**2))*(0.5*(vecL+vecR) - gr_V(DENS_VAR:PRES_VAR))
    C(:,2) = (vecR-vecL)/gr_dx
    C(:,1) = gr_V(DENS_VAR:PRES_VAR,i) - C(:,3)*(gr_dx**2)/12.

    ! NEED TO IMPLEMENT CONDITION 1 on Page 159 of the lecture note. 
    ! it will probably be easier to just return the FOG method rather than shift
    ! the parameters vecL, vecR, C
    do nVar = DENS_VAR, PRES_VAR 
      if ((vecR(nvar) - gr_V(nvar,i))*(-vecL(nvar) + gr_V(nvar,i)) .le. 0) then
        C(nVar,2:3) = 0
      else
        if (-(vecR(nVar)-vecL(nVar))**2 .gt. 6*(vecR(nVar)-vecL(nVar))*&
          (gr_V(nVar,i)-(vecR(nVar)+vecL(nVar))/2)) then
          vecR(nVar) = 3*gr_V(nVar,i) - 2*vecL(nVar)
        elseif (-(vecR(nVar)-vecL(nVar))**2 .gt. 6*(vecR(nVar)-vecL(nVar))*&
          (gr_V(nVar,i)-(vecR(nVar)+vecL(nVar))/2)) then
          vecL(nVar) = 3*gr_V(nVar,i) - 2*vecR(nVar)
        endif
      endif
    end do 
    ! update the coefficients after the correction
    C(:,3) = (6./(gr_dx**2))*(0.5*(vecL+vecR) - gr_V(DENS_VAR:PRES_VAR))
    C(:,2) = (vecR-vecL)/gr_dx
    C(:,1) = gr_V(DENS_VAR:PRES_VAR,i) - C(:,3)*(gr_dx**2)/12.


    call eigenvalues(gr_V(DENS_VAR:GAME_VAR,i),lambda)
    call left_eigenvectors (gr_V(DENS_VAR:GAME_VAR,i),conservative,leig)
    call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR,i),conservative,reig)

    ! primitive limiting
    if (.not. sim_charLimiting) then
       do kWaveNum = 1, NUMB_WAVE
          ! slope limiting
          ! deltas in primitive vars
          !delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
          !delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
          !do nVar = DENS_VAR,PRES_VAR
             !if (sim_limiter == 'minmod') then
                !call minmod(delL(nVar),delR(nVar),delV(nVar))
             !elseif (sim_limiter == 'vanLeer') then
                !call vanLeer(delL(nVar),delR(nVar),delV(nVar))
             !elseif (sim_limiter == 'mc') then
                !call mc(delL(nVar),delR(nVar),delV(nVar))
             !endif
          !enddo
          ! project primitive delta to characteristic vars
          delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),dvC(DENS_VAR:PRES_VAR))
       enddo
    elseif (sim_charLimiting) then
        !stop
        !STUDENTS: PLEASE FINISH THIS CHARACTERISTIC LIMITING
        !(THE IMPLEMENTATION SHOULD NOT BE LONGER THAN THE PRIMITIVE LIMITING CASE)
        do kWaveNum = 1, NUMB_WAVE
           ! slope limiting
           ! deltas in primitive vars
           !delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
           !delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
           !do nVar = DENS_VAR,PRES_VAR
           if (sim_limiter == 'minmod') then
              !call minmod(delL(nVar),delR(nVar),delV(nVar))
              call minmod(&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delL(DENS_VAR:PRES_VAR)),&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delR(DENS_VAR:PRES_VAR)),&
                delW(kWaveNum))
           elseif (sim_limiter == 'vanLeer') then
              !call vanLeer(delL(nVar),delR(nVar),delV(nVar))
              call vanLeer(&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delL(DENS_VAR:PRES_VAR)),&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delR(DENS_VAR:PRES_VAR)),&
                delW(kWaveNum))
           elseif (sim_limiter == 'mc') then
              !call mc(delL(nVar),delR(nVar),delV(nVar))
              call mc(&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delL(DENS_VAR:PRES_VAR)),&
                dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),&
                  delR(DENS_VAR:PRES_VAR)),&
                delW(kWaveNum))
           endif
           !enddo
           !delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),delV(DENS_VAR:PRES_VAR))
           ! project primitive delta to characteristic vars
        enddo
    endif



     ! set the initial sum to be zero
     sigL(DENS_VAR:ENER_VAR) = 0.
     sigR(DENS_VAR:ENER_VAR) = 0.
     vecL(DENS_VAR:ENER_VAR) = 0.
     vecR(DENS_VAR:ENER_VAR) = 0.
     sig2L(DENS_VAR:ENER_VAR) = 0.
     sig2R(DENS_VAR:ENER_VAR) = 0.
     vec2L(DENS_VAR:ENER_VAR) = 0.
     vec2R(DENS_VAR:ENER_VAR) = 0.
    
     do kWaveNum = 1, NUMB_WAVE
        ! lambdaDtDx = lambda*dt/dx
        lambdaDtDx = lambda(kWaveNum)*dt/gr_dx
        delC1 = gr_dx*dot_product(leig(DENS_VAR:PRES_VAR,nWaveNum),C(:,2))
        delC2 = (gr_dx**2)*dot_product(leig(DENS_VAR:PRES_VAR,nWaveNum),C(:,3))
        !if (sim_riemann == 'roe') then
        ! pretty sure that this solver only uses ROE
          if (lambdaDtDx .gt. 0) then
            vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delC1
            vec2R(DENS_VAR:PRES_VAR) = 0.25*(1.0 - 2*lambdaDtDx+(4./3)*lambdaDtDx**2)*&
              reig(DENS_VAR:PRES_VAR,kWaveNum)*delC2
            sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)
            sig2R(DENSE_VAR:PRES_VAR) = sig2R(DENS_VAR:PRES_VAR) + vec2R(DENS_VAR:PRES_VAR)
          else
            vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delC1
            vec2L(DENS_VAR:PRES_VAR) = 0.25*(1.0 + 2*lambdaDtDx+(4./3)*lambdaDtDx**2)*&
              reig(DENS_VAR:PRES_VAR,kWaveNum)*delC2
            sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
            sig2L(DENSE_VAR:PRES_VAR) = sig2L(DENS_VAR:PRES_VAR) + vec2L(DENS_VAR:PRES_VAR)
          end if
        !elseif (sim_riemann == 'hll') then
          !vecR(DENS_VAR:PRES_VAR) = 0.5*(1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
          !sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)

          !vecL(DENS_VAR:PRES_VAR) = 0.5*(-1.0 - lambdaDtDx)*reig(DENS_VAR:PRES_VAR,kWaveNum)*delW(kWaveNum)
          !sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
        !endif

        ! Let's make sure we copy all the cell-centered values to left and right states
        ! this will be just FOG
        gr_vL(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
        gr_vR(DENS_VAR:NUMB_VAR,i) = gr_V(DENS_VAR:NUMB_VAR,i)
        
        ! Now PLM reconstruction for dens, velx, and pres
        gr_vL(DENS_VAR:PRES_VAR,i) = gr_V(DENS_VAR:PRES_VAR,i) + sigL(DENS_VAR:PRES_VAR) + sig2L(DENS_VAR:PRES_VAR)
        gr_vR(DENS_VAR:PRES_VAR,i) = gr_V(DENS_VAR:PRES_VAR,i) + sigR(DENS_VAR:PRES_VAR) + sig2R(DENS_VAR:PRES_VAR)
     end do
  end do
 
  return
end subroutine soln_PPM
