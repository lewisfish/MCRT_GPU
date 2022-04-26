module stokesMOD

   use randomMod
   use photonMod

   implicit none

   contains

      subroutine stokes(packet,hgg,g2,pi,twopi,iseed)

         use iso_fortran_env, only : real32

         implicit none
         !$acc routine nohost

         integer, intent(INOUT) :: iseed
         type(photon), intent(INOUT) :: packet
         real(kind=real32), intent(IN) :: hgg, g2, pi, twopi
         real(kind=real32) :: costp, sintp, phip, bmu, b, ri1, ri3, cosi3, sini3
         real(kind=real32) :: cosb2, sinbt, cosi2, sini1, cosi1, sini2, bott, cosdph
      
         !***** isotropic scattering if g = 0.0 ******************************
         if(hgg == 0.0) then
         packet%cost=2. * ran2(iseed) - 1.
         packet%sint=(1. -packet%cost**2)
            if(packet%sint <= 0.)then
            packet%sint = 0.
            else
            packet%sint=sqrt(packet%sint)
            endif
      
         packet%phi=TWOPI*ran2(iseed)
         packet%sinp=sin(packet%phi)
         packet%cosp=cos(packet%phi)
      
         packet%nxp=packet%sint*packet%cosp
         packet%nyp=packet%sint*packet%sinp
         packet%nzp=packet%cost
      
         else
      
         !***** heyney greenstein scattering ********************************
      
            costp=packet%cost
            sintp=packet%sint
            phip=packet%phi
      
            bmu=((1.+g2)-((1.-g2)/(1.-hgg+2.*hgg*ran2(iseed)))**2)/(2.*hgg)
            cosb2=bmu**2
            b=cosb2-1.
      
            if(abs(bmu) > 1.) then
               if(bmu > 1.) then
                  bmu=1.
                  cosb2=1.
                  b=0.
               else
                  bmu=-1.
                  cosb2=1.
                  b=0.
               end if
            end if
            sinbt=sqrt(1.-cosb2)
            ri1=TWOPI*ran2(iseed)
      
            if(ri1 > PI) then
               ri3=TWOPI-ri1
               cosi3=cos(ri3)
               sini3=sin(ri3)
      
               if(bmu == 1.) then
                  return
               else
                  if(bmu == -1.) then
                     return
                  end if
               end if
      
            packet%cost=costp*bmu+sintp*sinbt*cosi3
               if(abs(packet%cost) < 1.) then
               packet%sint=abs(sqrt(1. -packet%cost**2))
                  sini2=sini3*sintp/packet%sint
                  bott=packet%sint*sinbt
                  cosi2=costp/bott-packet%cost*bmu/bott
               else
               packet%sint=0.
                  sini2=0.
                  if(packet%cost >= 1.)  cosi2=-1.
                  if(packet%cost <= -1.) cosi2=1.
               end if
      
               cosdph=-cosi2*cosi3+sini2*sini3*bmu
               if(abs(cosdph) > 1.) then
                  if(cosdph > 1.) then
                     cosdph=1.
                  else
                     cosdph=-1.
                  end if
               end if
               
            packet%phi=phip+acos(cosdph)
               if(packet%phi > TWOPI)packet%phi=packet%phi-TWOPI
               if(packet%phi.lt.0.)   packet%phi=packet%phi+TWOPI
      
            else  
               cosi1=cos(ri1)
               sini1=sin(ri1)
               if(bmu == 1.) then
                  return
               else
                  if(bmu == -1.) then
                     return
               end if
            end if
      
            packet%cost=costp*bmu+sintp*sinbt*cosi1
               if(abs(packet%cost).lt.1.) then
               packet%sint=abs(sqrt(1. -packet%cost**2))
                  sini2=sini1*sintp/packet%sint
                  bott=packet%sint*sinbt
                  cosi2=costp/bott-packet%cost*bmu/bott
               else
               packet%sint=0.
                  sini2=0.
                  if(packet%cost >= 1.)  cosi2=-1.
                  if(packet%cost <= -1.) cosi2=1.
               end if
      
               cosdph=-cosi1*cosi2+sini1*sini2*bmu
               if(abs(cosdph) > 1.) then
                  if(cosdph > 1.) then
                     cosdph=1.
                  else
                     cosdph=-1.
                  end if
               end if
            packet%phi=phip-acos(cosdph)
               if(packet%phi > TWOPI)packet%phi=packet%phi-TWOPI
               if(packet%phi < 0.)   packet%phi=packet%phi+TWOPI
            end if
      
         packet%cosp=cos(packet%phi)
         packet%sinp=sin(packet%phi)
      
         packet%nxp=packet%sint*packet%cosp
         packet%nyp=packet%sint*packet%sinp
         packet%nzp=packet%cost
      
         end if
   end subroutine stokes
end module stokesMOD