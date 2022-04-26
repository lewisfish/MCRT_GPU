module photonMod

    implicit none

    type photon
        real :: xp, yp, zp
        real :: nxp, nyp, nzp
        real :: cosp, sinp, sint, cost, phi
        integer :: xcell, ycell, zcell
    end type photon

end module photonMod