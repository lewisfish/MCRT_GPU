module tauint2MOD

    use randomMod
    use photonMod
    use iso_fortran_env, only : real32

    implicit none

    contains

        subroutine tauint2(packet,xmax,ymax,zmax,xface,yface,zface,rhokap,jmean,tflag,h,delta,nxg,nyg,nzg)

            implicit none
            !$acc routine nohost

            type(photon) :: packet
            integer, intent(IN) :: nxg,nyg,nzg

            real(kind=real32), intent(IN) :: xface(:),yface(:),zface(:), rhokap(:,:,:)
            real(kind=real32), intent(INOUT) :: jmean(201,201,201)

            integer, intent(INOUT) :: tflag, h
            real(kind=real32), intent(IN) :: xmax,ymax,zmax,delta

            integer :: celli,cellj,cellk
            real(kind=real32) :: tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
            real(kind=real32) :: dx,dy,dz,smax

            ! tflag=0 means photon is in envelope
            tflag=0

            ! generate random optical depth tau
            tau=-log(ran2(h))

            ! set the cumulative distance and optical depth (d and taurun) 
            ! along the photon path to zero.  set the current photon coordinates.
            ! note that the origin of the (xcur,ycur,zcur) system is at the 
            ! bottom corner of the grid.
            taurun=0.
            d=0.
            xcur=packet%xp+xmax
            ycur=packet%yp+ymax
            zcur=packet%zp+zmax

            celli=packet%xcell
            cellj=packet%ycell
            cellk=packet%zcell

            ! calculate smax -- maximum distance photon can travel
            if(packet%nxp > 0.) then
                dsx=(2.*xmax-xcur)/packet%nxp
            elseif(packet%nxp < 0.) then
                dsx=-xcur/packet%nxp
            elseif(packet%nxp == 0.) then
                dsx=1.e2*xmax
            endif

            if(packet%nyp > 0.) then
                dsy=(2.*ymax-ycur)/packet%nyp
            elseif(packet%nyp < 0.) then
                dsy=-ycur/packet%nyp
            elseif(packet%nyp == 0.) then
                dsy=1.e2*ymax
            endif

            if(packet%nzp > 0.) then
                dsz=(2.*zmax-zcur)/packet%nzp
            elseif(packet%nzp < 0.) then
                dsz=-zcur/packet%nzp
            elseif(packet%nzp == 0.) then
                dsz=1.e2*zmax
            endif

            smax=min(dsx,dsy,dsz)
            if(smax < delta) then
                tflag=1
                return
            endif

            ! integrate through grid
            do while((taurun < tau).and.(d < (.999*smax)))

                ! find distance to next x, y, and z cell walls.  
                ! note that dx is not the x-distance, but the actual distance along 
                !the direction of travel to the next x-face, and likewise for dy and dz.
                if(packet%nxp > 0.) then
                    dx=(xface(celli+1)-xcur)/packet%nxp
                    if(dx < delta) then
                        xcur=xface(celli+1)
                        celli=celli+1
                        dx=(xface(celli+1)-xcur)/packet%nxp
                    endif
                elseif(packet%nxp < 0.) then
                    dx=(xface(celli)-xcur)/packet%nxp
                    if(dx < delta) then
                        xcur=xface(celli)
                        dx=(xface(celli-1)-xcur)/packet%nxp
                        celli=celli-1
                    endif
                elseif(packet%nxp == 0.) then
                    dx=1.e2*xmax
                endif

                if(packet%nyp > 0.) then
                    dy=(yface(cellj+1)-ycur)/packet%nyp
                    if(dy < delta) then
                        ycur=yface(cellj+1)
                        cellj=cellj+1
                        dy=(yface(cellj+1)-ycur)/packet%nyp
                    endif
                elseif(packet%nyp < 0.) then
                    dy=(yface(cellj)-ycur)/packet%nyp
                    if(dy < delta) then
                        ycur=yface(cellj)
                        dy=(yface(cellj-1)-ycur)/packet%nyp
                        cellj=cellj-1
                    endif
                elseif(packet%nyp == 0.) then
                    dy=1.e2*ymax
                endif

                if(packet%nzp > 0.) then
                    dz=(zface(cellk+1)-zcur)/packet%nzp
                    if(dz < delta) then
                        zcur=zface(cellk+1)
                        cellk=cellk+1
                        dz=(zface(cellk+1)-zcur)/packet%nzp
                    endif
                elseif(packet%nzp < 0.) then
                    dz=(zface(cellk)-zcur)/packet%nzp
                    if(dz < delta) then
                        zcur=zface(cellk)
                        dz=(zface(cellk-1)-zcur)/packet%nzp
                        cellk=cellk-1
                    endif
                elseif(packet%nzp == 0.) then
                    dz=1.e2*zmax
                endif

                if((celli > nxg).or.(celli < 0)) print*,'celli',celli,xcur
                if((cellj > nyg).or.(cellj < 0)) print*,'celli',cellj,ycur
                if((cellk > nyg).or.(cellk < 0)) print*,'celli',cellk,zcur


                ! distances are only zero if photon is on cell wall.  if it is 
                ! on cell wall then set to arbitrary large distance, since we will
                ! in fact hit another wall

                if( (dx == 0.) .or. ((abs(dx)) < (delta)) ) dx=1.e2*xmax
                if( (dy == 0.) .or. ((abs(dy)) < (delta)) ) dy=1.e2*ymax
                if( (dz == 0.) .or. ((abs(dz)) < (delta)) ) dz=1.e2*zmax

                ! find distance to next cell wall -- minimum of dx, dy, and dz
                dcell=min(dx,dy,dz)
                if(dcell <= 0.) then
                    print *,'tauint2: dcell < 0',dx,dy,dz,packet%nxp,packet%nyp,packet%nzp
                    print *,xcur,ycur,zcur,celli,cellj,cellk
                endif
                if(dx < 0.) dcell=min(dy,dz)
                if(dy < 0.) dcell=min(dx,dz)
                if(dz < 0.) dcell=min(dx,dy)

                ! optical depth to next cell wall is 
                ! taucell= (distance to cell)*(opacity of current cell)
                taucell=dcell*rhokap(celli,cellj,cellk)

                ! if taurun+taucell>tau then scatter at distance d+d1.  
                ! update photon position and cell.  
                ! if taurun+taucell<tau then photon moves distance dcell 
                ! (i.e. ends up on next cell wall) and update photon position
                ! and cell.
                if((taurun+taucell).ge.tau) then
                    d1=(tau-taurun)/rhokap(celli,cellj,cellk)
                    d=d+d1
                    !$acc atomic
                    jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + d1
                    taurun=taurun+taucell
                    xcur=xcur+d1*packet%nxp
                    ycur=ycur+d1*packet%nyp
                    zcur=zcur+d1*packet%nzp

                    ! Linear Grid
                    celli=int(nxg*xcur/(2.*xmax))+1
                    cellj=int(nyg*ycur/(2.*ymax))+1
                    cellk=int(nzg*zcur/(2.*zmax))+1
                else
                    !$acc atomic
                    jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + dcell
                    d=d+dcell
                    taurun=taurun+taucell
                    xcur=xcur+dcell*packet%nxp
                    ycur=ycur+dcell*packet%nyp
                    zcur=zcur+dcell*packet%nzp

                    ! Linear Grid 
                    celli=int(nxg*xcur/(2.*xmax))+1
                    cellj=int(nyg*ycur/(2.*ymax))+1
                    cellk=int(nzg*zcur/(2.*zmax))+1
                endif
            end do

            ! calculate photon final position.  if it escapes envelope then
            ! set tflag=1.  if photon doesn't escape leave tflag=0 and update 
            ! photon position.
            if((d.ge.(.999*smax))) then
                tflag=1
            else
                packet%xp=packet%xp+d*packet%nxp
                packet%yp=packet%yp+d*packet%nyp
                packet%zp=packet%zp+d*packet%nzp
                packet%xcell=int(nxg*(packet%xp+xmax)/(2.*xmax))+1
                packet%ycell=int(nyg*(packet%yp+ymax)/(2.*ymax))+1
                packet%zcell=int(nzg*(packet%zp+zmax)/(2.*zmax))+1
            endif
        end subroutine tauint2
end module tauint2MOD