module densityMOD
    implicit none
    contains
        subroutine density(x,y,z,rho)

            use iso_fortran_env, only : real32

            implicit none

            real(kind=real32) :: x,y,z,rho
            real(kind=real32) :: r

            ! calculate some distances for use in setting up density 
            ! structure. Note that distances are in units of xmax, ymax, and zmax 
            ! as called from the loop over cells in gridset.f
            r=sqrt(x**2 + y**2 + z**2)

            ! Set up uniform density sphere within the grid
            if(r > 1.) then
                rho=0.
            else
                rho=1.
            endif
        end subroutine density
end module densityMOD