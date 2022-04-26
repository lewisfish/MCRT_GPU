program mcpolar
    use stokesMOD
    use tauint2MOD
    use sourcephMOD
    use iarrayMOD
    use gridsetMOd
    use openacc
    use randomMod   
    use iso_fortran_env, only : real32
    use photonMod

    implicit none

    integer, parameter :: nxg=201,nyg=201,nzg=201

    type(photon) :: packet
    real(kind=real32), allocatable :: rhokap(:, :, :), xface(:), yface(:), zface(:),jmean(:, :, :)
    real :: ran

    ! Parameter declarations
    integer :: nphotons, iseed, j, tflag, i, seed, seq, offset
    real(kind=real32) :: nscatt, kappa,albedo,hgg,xmax,ymax,zmax, pi,twopi,g2,delta

    logical :: rflag
    ! Read in parameters from the file input.params
    open(newunit=j,file='res/input.params',status='old')
    read(j,*) nphotons
    read(j,*) kappa
    read(j,*) albedo
    read(j,*) hgg
    read(j,*) xmax
    read(j,*) ymax
    read(j,*) zmax
    close(j)

    pi=4.*atan(1.)
    twopi=2.*pi
    g2=hgg*hgg  ! Henyey-Greenstein parameter, hgg^2

    ! Initialize arrays to zero
    call iarray(xface,yface,zface,rhokap,jmean,nxg,nyg,nzg)

    ! Set up density grid
    call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,kappa,nxg,nyg,nzg)

    ! Set small distance for use in optical depth integration routines 
    ! for roundoff effects when crossing cell walls
    delta=1.e-5*(2.*xmax/nxg)

    nscatt=0
    ! Setup OpenACC
    !$acc parallel private(rflag) copy(jmean,rhokap,xface,yface,zface)
    rflag = .false.
    !$acc loop private(packet,tflag, iseed) independent reduction(+:nscatt)
    
    ! Loop over nph photons from each source
    do j=1,nphotons
        !Set seed for OpenACC
        if(.not.rflag)then
            iseed = abs(1000 * (__pgi_gangidx()+__pgi_vectoridx()+__pgi_workeridx()))
            rflag = .true.
        end if

        if(mod(j,10000) == 0)print *, j,' scattered photons completed'

        ! Release photon from point source
        call sourceph(packet,xmax,ymax,zmax,twopi,nxg,nyg,nzg,iseed)

        ! Find scattering location
        call tauint2(packet,xmax,ymax,zmax,xface,yface,zface,rhokap,jmean,tflag,iseed,delta,nxg,nyg,nzg)

        ! Photon scatters in grid until it exits (tflag=1) 
        do while(tflag == 0)
            if(ran2(iseed) < albedo) then
                ! Scatter photon into new direction and update Stokes parameters
                call stokes(packet,hgg,g2,pi,twopi, iseed)
                nscatt=nscatt+1
            else
                exit
            endif

            ! Find next scattering location
            call tauint2(packet,xmax,ymax,zmax,xface,yface,zface,rhokap,jmean,tflag,iseed,delta,nxg,nyg,nzg)
        end do
    end do
    !$acc end parallel

    print*,'Avereage number of scatterings = ',(nscatt/nphotons)
    open(newunit=j,file="data/jmean.dat", access="stream", form="unformatted")
    write(j)jmean * ((2.*xmax*2.*ymax)/(nphotons * (2. * xmax / nxg) * (2. * ymax / nyg) * (2. * zmax / nzg)))
    close(j)
end program mcpolar