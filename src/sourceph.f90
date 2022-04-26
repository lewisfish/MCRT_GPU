module sourcephMOD

    use randomMod
    use photonMod
    use iso_fortran_env, only : real32

    implicit none

    contains
        subroutine sourceph(packet,xmax,ymax,zmax,twopi,nxg,nyg,nzg,h)

            implicit none
            !$acc routine nohost

            real(kind=real32), intent(IN) :: xmax,ymax,zmax,twopi
            integer, intent(IN)           :: nxg,nyg,nzg
            integer, intent(INOUT)        :: h
            type(photon), intent(OUT)     :: packet

            ! emit photon isotropically from origin
            packet%xp=0._real32
            packet%yp=0._real32
            packet%zp=0._real32

            packet%cost=2. * ran2(h) - 1.
            packet%sint=(1. - packet%cost**2)
            packet%sint=sqrt(packet%sint)

            packet%phi=twopi*ran2(h)
            packet%cosp=cos(packet%phi)
            packet%sinp=sin(packet%phi)

            ! Set photon direction cosines for direction of travel
            packet%nxp=packet%sint*packet%cosp  
            packet%nyp=packet%sint*packet%sinp
            packet%nzp=packet%cost
            ! Linear Grid 
            packet%xcell=int(nxg*(packet%xp+xmax)/(2.*xmax))+1
            packet%ycell=int(nyg*(packet%yp+ymax)/(2.*ymax))+1
            packet%zcell=int(nzg*(packet%zp+zmax)/(2.*zmax))+1
        end subroutine sourceph
end module sourcephMOD