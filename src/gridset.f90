module gridsetMOD
    implicit none

    contains
        subroutine gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,kappa,nxg,nyg,nzg)

            use densityMOD
            use iso_fortran_env, only : real32

            implicit none

            integer, intent(IN) :: nxg,nyg,nzg

            real(kind=real32), intent(INOUT) :: xface(:),yface(:),zface(:), rhokap(:,:,:)
            real(kind=real32), intent(IN) :: xmax,ymax,zmax,kappa

            integer :: i,j,k
            real(kind=real32) :: x,y,z,rho

            print *, 'Setting up density grid....'

            ! Linear Cartesian grid. Set up grid faces
            do i=1,nxg+1
            xface(i)=(i-1)*2.*xmax/nxg
            end do
            do i=1,nyg+1
            yface(i)=(i-1)*2.*ymax/nyg
            end do
            do i=1,nzg+1
            zface(i)=(i-1)*2.*zmax/nzg
            end do

            ! Loop through x, y, and z to set up grid density.
            do i=1,nxg
                do j=1,nyg
                    do k=1,nzg
                        x=xface(i)-xmax+xmax/nxg
                        y=yface(j)-ymax+ymax/nyg
                        z=zface(k)-zmax+zmax/nzg

                        !Call density setup subroutine 
                        call density(x,y,z,rho)
                        rhokap(i,j,k)=rho*kappa
                    end do
                end do
            end do
        end subroutine gridset
end module gridsetMOD