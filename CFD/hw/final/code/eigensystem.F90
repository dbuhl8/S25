module eigensystem

#include "definition.h"
  
  use grid_data

contains

  subroutine eigenvalues(V,lambda)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda

    real :: a, u

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))!;print*,a,V(GAMC_VAR),V(PRES_VAR),V(DENS_VAR)
    u = V(VELX_VAR)
    
    lambda(SHOCKLEFT) = u - a
    lambda(CTENTROPY) = u
    lambda(SHOCKRGHT) = u + a
    
    return
  end subroutine eigenvalues


  
  subroutine right_eigenvectors(V,conservative,reig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig

    real :: a, u, d, g, ekin, hdai, hda
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)
    d = V(DENS_VAR)
    g = V(GAMC_VAR) - 1.
    ekin = 0.5*u**2
    hdai = 0.5*d/a
    hda  = 0.5*d*a
    
    if (conservative) then
       !! Conservative eigenvector
       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VELX_VAR,SHOCKLEFT) = u - a
       reig(PRES_VAR,SHOCKLEFT) = ekin + a**2/g - a*u
       reig(:,SHOCKLEFT) = -hdai*reig(:,SHOCKLEFT)

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VELX_VAR,CTENTROPY) = u
       reig(PRES_VAR,CTENTROPY) = ekin
       
       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VELX_VAR,SHOCKRGHT) = u + a
       reig(PRES_VAR,SHOCKRGHT) = ekin + a**2/g + a*u
       reig(:,SHOCKRGHT) = hdai*reig(:,SHOCKRGHT)
       
    else
       !! Primitive eigenvector
       !! STUDENTS: PLEASE FINISH THIS PRIMITIVE RIGHT EIGEN VECTORS
       print*,'eigeysystem.F90: right eigenvectors'
       stop
    endif
    
    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: a, u, d, g, ekin, hdai, hda
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VELX_VAR)
    d = V(DENS_VAR)
    g = V(GAMC_VAR) - 1.
    ekin = 0.5*u**2
    hdai = 0.5*d/a
    hda  = 0.5*d*a
    
    if (conservative) then
       !! Conservative eigenvector
       !! STUDENTS: PLEASE FINISH THIS CONSERVATIVE LEFT EIGEN VECTORS
       print*,'eigeysystem.F90: left conservative eigenvectors'
       stop
       
    else
       !! Primitive eigenvector
       !! STUDENTS: PLEASE FINISH THIS PRIMITIVE LEFT EIGEN VECTORS
       print*,'eigeysystem.F90: left prim eigenvectors'
       stop

    endif
    
    return
  end subroutine left_eigenvectors


  
end module eigensystem
