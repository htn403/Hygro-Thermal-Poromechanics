      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME
      FLUX(1) = -1.0d0*10.0d0**(-8.0d0)*(SOL-0.68d0)
      RETURN
      END
      
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
      FIELD(1) = STATEV(7)
      RETURN
      END
      
C ABAQUS format user material subroutine UMAT for simplified 
C nonlinear diffusion model
C  SVD1-6  = shrinkage strain (d) IV=0
C  SVD7  = Cond  (d) 33.72e-11
C  SVD8  = Kp  (d) IV=m1=1.5d0 + ((wc-0.18d0)/0.15)**2
C  SVD9  = Alpha (d) IV=0.05      
C  temp  = h       IV = 1.0d0      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      double precision :: young, nuy, lambda, muy, alp, wc, hf, au, da,
     1 s, Kp, m1, m2, c1, hc, cf_0, bet, h_p, cf, r, k0, TEMPNEW, ca, dh
      double precision, dimension(NTENS,NTENS) :: ones, diag
      double precision, dimension(NTENS,NTENS) :: DDSDE
      double precision, dimension(NTENS,NTENS) :: DDSDS
      double precision, dimension(NTENS,NTENS) :: DSH
      double precision, dimension(NTENS) :: SIGold


      integer :: j

C-----------------------------------------------------------------------------
C                    Start new variable definitions here
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
C              Read in material properties from the input file
C-----------------------------------------------------------------------------
      young = PROPS(1)
      nuy = PROPS(2)
      wc  = 0.3
C     Calculate useful variables
      lambda  = young*nuy/((1+nuy)*(1-2*nuy))
      muy     = young/(2*(1+nuy))
      hf      = 0.75d0
      m1      = 1.5d0 + ((wc-0.18d0)/0.15d0)**2  !m1 for sorption isotherm slope
      m2      = 0.73d0                            !m2 for sorption isotherm slope
      h_p     = 1.0d0-0.13d0/(m1**1.3d0)         !h* for isotherm slope
      r       = 2.0d0
      k0      = 130.0d0*10.0d0**(-5.0d0)                  !shrinkage strain and humidity constant
      au      = 0.46d0 + 0.95d0*(wc-0.17d0)**(0.6d0) !Ultimate shrinkage strain
      ca      = 10.0d0**(-4.0d0)
!Initialize
      DDSDDE = 0.d0
      DDSDE  = 0.d0
      DDSDS  = 0.d0
      DSH    = 0.d0
      SIGold = 0.d0
      ones   = 0.d0
      diag   = 0.d0
      DDSDDT  = 0.d0
!Define Tensor coefficients
      ones(1:NTENS,1:NTENS) = 1.d0
      diag(1:NTENS,1:NTENS) = 0.d0
      forall(K1=1:NTENS) diag(K1,K1) = 1.d0

!Define DDSDDE, DDSDE, DDSDS related to Normal Stress
      DDSDE(1:NDI,1:NDI) = lambda*
     1 ones(1:NDI,1:NDI) + 2.d0*muy*diag(1:NDI,1:NDI)
      
      DDSDDE(1:NDI,1:NDI) = lambda*
     1 ones(1:NDI,1:NDI) + 2.d0*muy*diag(1:NDI,1:NDI)
      
      DDSDS(1:NDI,1:NDI) = -diag(1:NDI,1:NDI)
      DSH(1:NDI,1:NDI)   = ones(1:NDI,1:NDI)      

 !Define DDSDDE, DDSDE, DDSDS related to Shear Stress
      DDSDE(NDI+1:NTENS,NDI+1:NTENS) = muy
     1 *diag(NDI+1:NTENS,NDI+1:NTENS)
      DDSDDE(NDI+1:NTENS,NDI+1:NTENS) = muy
     1 *diag(NDI+1:NTENS,NDI+1:NTENS)
      DDSDS(NDI+1:NTENS,NDI+1:NTENS) = -diag(NDI+1:NTENS,NDI+1:NTENS)
      DSH(NDI+1:NTENS,NDI+1:NTENS)   = 0.0d0
 !Store present stress state ub array SIGold
      do i = 1,ntens
      SIGold(i) = stress(i)
      end do

 !Update stresses
C      DESH = k0*(au-0.9d0*0.05d0)*DTEMP/(STATEV(9)-0.9d0*0.05d0)!Shrinkage strain increment
      DESH = k0*DTEMP/STATEV(8)                             !Shrinkage strain increment
      do i = 1,ntens
      do j = 1,ntens
      stress(i) = stress(i) + DDSDDE(i,j)*(DSTRAN(j)-
     1 DSH(i,j)*DESH) + DDSDE(i,j)*(STRAN(j)-DSH(i,j)*STATEV(j)) 
     2 + DDSDS(i,j)*(SIGold(j)) + diag(i,j)*DDSDDT(i)*TEMP 
      end do
      end do
      do i = 1,ndi
      STATEV(i) = STATEV(i)+DESH
      end do
C-----------------------------------------------------------------------------
C       Update hydration degree
C-----------------------------------------------------------------------------
      if (TEMP .gt. hf) then
          da = ca*(au-STATEV(9))**2.0d0*(TEMP-hf)*DTIME
      else
          da = 0.0d0
      end if
      STATEV(9) = STATEV(9) + da
C-----------------------------------------------------------------------------
C       Calculate the desorption slope and current saturation degree
C-----------------------------------------------------------------------------
      TEMPNEW       = TEMP+DTEMP
      !Isotherm slope
      Kp            = 1.0d0/(m2+(m1-m2)/(1.0d0
     1 +((1.0d0-TEMPNEW)/(1.0d0-h_p))**3.0d0))
      STATEV(8) = Kp
C-----------------------------------------------------------------------------
C       Update the conductivity value using field variable 1
C-----------------------------------------------------------------------------
      c1             = 60.0d0*(1.0d0+12.0d0*(wc-0.17d0)**2.0d0)
     1 *STATEV(9)/au
      hc             = 0.77d0+0.22d0*(wc-0.17d0)**(0.5d0)
     1 + 0.15d0*(au/STATEV(9)-1)
      if (hc .gt. 0.99d0) then
          hc = 0.99d0
      end if
      cf_0           = 60.0d0*(1.0d0+12.0d0*(wc-0.17d0)**2.0d0)
     1 *STATEV(9)/au
      if (TEMPNEW .gt. 0.95d0) then
          cf = cf_0
      else
          cf = 0.1d0*cf_0+0.9d0*cf_0*(TEMPNEW/0.95d0)**4.0d0
      end if
      bet           = cf/c1
C-----------------------------------------------------------------------------
C       Update the conductivity value using field variable 1
C----------------------------------------------------------------------------- 
      STATEV(7) = 10.0d0**(-11.0d0)*c1*Kp*(bet+(1.0d0-bet)
     1 /(1.0d0+((1.0d0-TEMPNEW)/(1.0d0-hc))**r))
C-----------------------------------------------------------------------------
      return
      RETURN
      END SUBROUTINE UMAT
      


C #============================================================================
C #============================================================================


