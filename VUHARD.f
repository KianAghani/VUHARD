      subroutine vuhard(
C Read only -
     *     nblock, 
     *     jElem, kIntPt, kLayer, kSecPt, 
     *     lAnneal, stepTime, totalTime, dt, cmname,
     *     nstatev, nfieldv, nprops, 
     *     props, tempOld, tempNew, fieldOld, fieldNew,
     *     stateOld,
     *     eqps, eqpsRate,
C Write only -
     *     yield, dyieldDtemp, dyieldDeqps,
     *     stateNew )
      include 'vaba_param.inc'

      dimension props(nprops), tempOld(nblock), tempNew(nblock),
     1   fieldOld(nblock,nfieldv), fieldNew(nblock,nfieldv),
     2   stateOld(nblock,nstatev), eqps(nblock), eqpsRate(nblock),
     3   yield(nblock), dyieldDtemp(nblock), dyieldDeqps(nblock,2),
     4   stateNew(nblock,nstatev), jElem(nblock)
      
	  character*80 cmname

      real*8 SIGMA0,T1,Q1,T2,Q2
      real*8 N,PDOT0
      real*8 A,Tr,Tt
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 vp,dvp
      real*8 Teff,tp,dtp
!-----------------------------------------------------------------------
!     Read material properties
!-----------------------------------------------------------------------
!Isotropic work-hardening
      SIGMA0 = props(1)
      T1     = props(2)
      Q1     = props(3)
      T2     = props(4)
      Q2     = props(5)
!Strain rate sensitivity
      N      = props(6)
      PDOT0  = props(7)
!Temperature sensitivity
      A      = props(8)
      Tr     = props(9)
      Tt     = props(10)

!-----------------------------------------------------------------------
!     Compute yield stress and its derivatives
!-----------------------------------------------------------------------
      do i=1,nblock
         yield(i) =(SIGMA0+Q1*(1.0-exp(-(T1/Q1)*eqps(i)))
     +                    +Q2*(1.0-exp(-(T2/Q2)*eqps(i))))
     +                    
!Derivative with respect to equivalent plastic strain   
         dyieldDeqps(i,1) =(T1*exp(-(T1/Q1)*eqps(i))
     +                     +T2*exp(-(T2/Q2)*eqps(i)))
!Derivative with respect to equivalent plastic strain rate
         dyieldDeqps(i,2) = 0.0
!Derivative with respect to temperature
         dyieldDtemp(i)   = 0.0
      end do
!-----------------------------------------------------------------------
!     Apply rate-sensitivity if needed
!-----------------------------------------------------------------------
      if((N.gt.0.0).and.(PDOT0.gt.0.0))then
         do i=1,nblock
!Compute rate effects and its derivative    
            vp  = (1.0+eqpsRate(i)/PDOT0)**N
            dvp = (1.0+eqpsRate(i)/PDOT0)**(N-1.0)*(N/PDOT0)
!Apply rate effects to the required quantities
			dyieldDeqps(i,1) = dyieldDeqps(i,1)*vp
            dyieldDeqps(i,2) = yield(i)*dvp
            yield(i) = yield(i)*vp
         end do
      end if
!-----------------------------------------------------------------------
!     Apply temperature-sensitivity if needed
!-----------------------------------------------------------------------
      if((A.gt.0.0).and.(Tt.gt.0.0))then
         do i=1,nblock
!Compute rate effects and its derivative
            Teff = ((min(TEMPNEW(i),0.9995*Tt)-Tr)/(Tt-Tr))
            tp   = (1.0-Teff**A)
            dtp  =-(A/(Tt-Tr))*(Teff**(A-1.0))
!Apply rate effects to the required quantities
            dyieldDtemp(i)   = yield(i)*dtp
			dyieldDeqps(i,1) = dyieldDeqps(i,1)*tp
            dyieldDeqps(i,2) = dyieldDeqps(i,2)*tp 
            yield(i) = yield(i)*tp
         end do
      end if

      return
      end
