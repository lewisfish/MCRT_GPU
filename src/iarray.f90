module iarrayMOD

    implicit none

    contains
        subroutine iarray(xface,yface,zface,rhokap,jmean,nxg,nyg,nzg)

            use iso_fortran_env, only : real32

            implicit none

            integer :: nxg,nyg,nzg

            real(kind=real32), allocatable, intent(OUT) :: xface(:),yface(:),zface(:)
            real(kind=real32), allocatable, intent(OUT) :: rhokap(:, :, :),jmean(:, :, :)

            ! Initialize array values to be zero
            allocate(rhokap(nxg, nyg, nzg))
            allocate(jmean(nxg, nyg, nzg))
            allocate(xface(nxg+1))
            allocate(yface(nyg+1))
            allocate(zface(nzg+1))
            xface= 0.
            yface= 0.
            zface= 0.

            rhokap= 0.
            jmean = 0.

        end subroutine iarray
end module iarrayMOD