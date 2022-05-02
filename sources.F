!------------------------------------------------------------------------------------
!
!      FILE sources.F
!
!      This file is part of the FUNWAVE-TVD program under the Simplified BSD license
!
!-------------------------------------------------------------------------------------
!
!    Copyright (c) 2016, FUNWAVE Development Team
!
!    (See http://www.udel.edu/kirby/programs/funwave/funwave.html
!     for Development Team membership)
!
!    All rights reserved.
!
!    FUNWAVE_TVD is free software: you can redistribute it and/or modify
!    it under the terms of the Simplified BSD License as released by
!    the Berkeley Software Distribution (BSD).
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!
!    1. Redistributions of source code must retain the above copyright notice, this
!       list of conditions and the following disclaimer.
!    2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!    The views and conclusions contained in the software and documentation are those
!    of the authors and should not be interpreted as representing official policies,
!    either expressed or implied, of the FreeBSD Project.
!
!-------------------------------------------------------------------------------------
!
!    SourceTerms is subroutine for all source terms
!
!    HISTORY:
!       05/01/2010 Fengyan Shi
!       09/26/2013 Babak Tehranirad, added 2D Cd
!       08/18/2015 YoungKwang Choi, modified viscosity breaking
!       02/08/2016 Fengyan Shi, corrected wavemaker corresponding to
!                               conservative form of momentum equations
!
! --------------------------------------------------
SUBROUTINE SourceTerms
  USE GLOBAL

# if defined (VESSEL)
  USE VESSEL_MODULE
# endif

# if defined (METEO)
  USE METEO_MODULE
# endif

  IMPLICIT NONE
  REAL,DIMENSION(Mloc,Nloc) :: nu_vis
  LOGICAL :: PQ_scheme = .FALSE.
  REAL(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
  INTEGER :: kd,kf
  REAL(SP) :: AA                          ! for all WK's (optimized)
  REAL(SP), DIMENSION(Nfreq) :: BB,CC     ! for WK_IRR (optimized)
  REAL(SP), DIMENSION(NumWaveComp,1) :: BB1 ! for WK_TIME (optimized)
  REAL(SP), DIMENSION(NumFreq) :: BB2,CC2 ! for WK_DATA (optimized)
  INTEGER :: kf_index ! Chen, Salatin
  ! eddy viscosity, Salatin
  REAL(SP),DIMENSION(Mloc,Nloc) :: temp_smg
  ! temp files for sourceX and sourceY
  REAL(SP),DIMENSION(Mloc,Nloc) :: tempX,tempY,tempX_vis,tempY_vis
  INTEGER :: VTYPE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# if defined (CARTESIAN)
  DXg=DX
  DYg=DY
# else
  ! only for wavemaker
  DXg=DX(1,1)
  DYg=DY(1,1)
# endif

  ! wavemaker fyshi 02/08/2016
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the following code is optimized for accelerated simulation
  ! Chen, Salatin

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Regular waves
  IF(WAVEMAKER(1:6)=='WK_REG')THEN
    AA = TANH(PI/(Time_ramp*Tperiod)*TIME)*D_gen
    DO J=jlo,jhi
      DO I=ilo,ihi
        WaveMaker_Mass(I,J)= AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
        *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0_SP*PI/Tperiod*TIME)
      ENDDO
    ENDDO
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Regular,
  ! Irregular, TMA 1D, JONSWAP 1D, JONSWAP 2D, New Irregular WK
  IF(WaveMaker(1:6)=='WK_IRR'.OR.WaveMaker(1:6)=='TMA_1D'  &
  .OR.WaveMaker(1:6)=='JON_1D'.OR.WaveMaker(1:6)=='JON_2D' &
  .OR.WaveMaker(1:10)=='WK_NEW_IRR')THEN
    AA = TANH(PI/(Time_ramp/FreqPeak)*TIME)
    DO kf=1,Nfreq
      BB(kf) = COS(OMGN_IR(KF)*TIME)
      CC(kf) = SIN(OMGN_IR(KF)*TIME)
    ENDDO
    DO J=jlo,jhi
      DO I=ilo,ihi
        WK_Source=ZERO
        DO kf=1,Nfreq
          WK_Source=WK_Source+AA*(Cm(I,J,kf)*BB(kf) &
          +Sm(I,J,kf)*CC(kf))
        ENDDO
        WaveMaker_Mass(I,J)=WK_Source
      ENDDO
    ENDDO
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Time Series WK
  IF(WAVEMAKER(1:7)=='WK_TIME')THEN
    AA = TANH(PI/(Time_ramp*PeakPeriod)*TIME)
    DO kf=1,NumWaveComp
      BB1(kf,1)=COS(2.0_SP*PI/WAVE_COMP(kf,1)*TIME-WAVE_COMP(kf,3))
    ENDDO
    DO J=jlo,jhi
      DO I=ilo,ihi
        WK_Source=ZERO
        DO kf=1,NumWaveComp
          WK_Source=WK_Source+AA*D_genS(kf) &
          *EXP(-Beta_genS(kf)*(xmk_wk(I)-Xc_WK)**2)*BB1(kf,1)
        ENDDO
        WaveMaker_Mass(I,J)=WK_Source
      ENDDO
    ENDDO
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Measured Spectra WK's
  IF(WAVEMAKER(1:9)=='WK_DATA2D')THEN
    AA = TANH(PI/(Time_ramp/FreqPeak)*TIME)
    DO KF=1,NumFreq
      BB2(KF) = COS(OMGN2D(KF)*TIME)
      CC2(KF) = SIN(OMGN2D(KF)*TIME)
    ENDDO
    DO J=jlo,jhi
      DO I=ilo,ihi
        WK_Source=ZERO
        DO kf=1,NumFreq
          WK_Source=WK_Source+AA*(Cm(I,J,kf)*BB2(kf)+ &
          Sm(I,J,kf)*CC2(kf))
        ENDDO
        WaveMaker_Mass(I,J)=WK_Source
      ENDDO
    ENDDO
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Measured Spectra WK's
  IF(WAVEMAKER(1:13)=='WK_NEW_DATA2D')THEN
    AA = TANH(PI/(Time_ramp/FreqPeak)*TIME)
    ! Chen, Salatin (Optimzied)
    DO kf_index=1,size(loop_index)
      KF = loop_index(kf_index)
      BB2(KF) = COS(OMGN2D(KF)*TIME)
      CC2(KF) = SIN(OMGN2D(KF)*TIME)
    ENDDO
    DO J=jlo,jhi
      DO I=ilo,ihi
        WK_Source=ZERO
        ! Chen, Salatin (Optimzied)
        DO kf_index = 1,size(loop_index)
          kf = loop_index(kf_index)
          WK_Source=WK_Source+AA*(Cm(I,J,kf)*BB2(kf)+ &
          Sm(I,J,kf)*CC2(kf))
        ENDDO
        WaveMaker_Mass(I,J)=WK_Source
      ENDDO
    ENDDO
  ENDIF

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  nu_vis=ZERO
  ! ykchoi(08.18.2015)
  ! new variable :: WAVEMAKER_VIS
  ! should include ghost cells 02/28/2017 fyshi

  IF(VISCOSITY_BREAKING.OR.WAVEMAKER_VIS)THEN
    nu_vis = nu_break
  ENDIF

  IF(DIFFUSION_SPONGE)THEN
    ! should include ghost cells 02/28/2017, fyshi
    nu_vis = nu_vis + nu_sponge
  ENDIF

107   format(500f12.6)

# if defined (CARTESIAN)
# if defined (MIXING)

  ! Sub-grid-scale Turbulence
  ! Horizontal Eddy viscosity
  temp_smg = SQRT(2.0*uMAx**2+2.0*vMAy**2+(uMAy+vMAx)*(uMAy+vMAx))
  nu_smg = C_smg**2*DX*DY*temp_smg

  ! depth gradient term
  DO J=Jbeg-1,Jend+1
    DO I=Ibeg-1,Iend+1

      ! Smagorinsky
      tempX(I,J) = 1.0_SP/DX*( &
      0.5_SP*(nu_smg(I,J)+nu_smg(I+1,J)) &
      *1.0_SP/DX*(P(I+1,J)-P(I,J))  &
      - 0.5_SP*(nu_smg(I,J)+nu_smg(I-1,J)) &
      *1.0_SP/DX*(P(I,J)-P(I-1,J))) &
      + 0.5_SP/DY*( &
      0.5_SP*(nu_smg(I,J+1)+nu_smg(I,J)) &
      *1.0_SP/DY*(P(I,J+1)-P(I,J)) &
      -0.5_SP*(nu_smg(I,J-1)+nu_smg(I,J)) &
      *1.0_SP/DY*(P(I,J)-P(I,J-1))) &
      + 1.0_SP/DY*( &
      nu_smg(I,J+1)*0.5_SP/DX*(Q(I+1,J+1)-Q(I-1,J+1)) &
      -nu_smg(I,J-1)*0.5_SP/DX*(Q(I+1,J-1)-Q(I-1,J-1)))

      tempY(I,J) = 1.0_SP/DY*( &
      0.5_SP*(nu_smg(I,J)+nu_smg(I,J+1)) &
      *1.0_SP/DY*(Q(I,J+1)-Q(I,J))  &
      - 0.5_SP*(nu_smg(I,J)+nu_smg(I,J-1)) &
      *1.0_SP/DY*(Q(I,J)-Q(I,J-1))) &
      + 0.5_SP/DX*( &
      0.5_SP*(nu_smg(I+1,J)+nu_smg(I,J)) &
      *1.0_SP/DX*(Q(I+1,J)-Q(I,J)) &
      -0.5_SP*(nu_smg(I-1,J)+nu_smg(I,J)) &
      *1.0_SP/DX*(Q(I,J)-Q(I-1,J))) &
      + 1.0_SP/DX*( &
      nu_smg(I+1,J)*0.5_SP/DY*(P(I+1,J+1)-P(I+1,J-1)) &
      -nu_smg(I-1,J)*0.5_SP/DY*(P(I-1,J+1)-P(I-1,J-1)))

      ! Viscosity Breaking
      IF(PQ_scheme)THEN
        ! it turns out P and Q are not exchanged at processor
        ! interface. it affects edges, make PQ_scheme=false
        tempX_vis(I,J) = 0.5_SP/DX*( &
        nu_vis(I+1,J)* &
        1.0_SP/DX*(P(I+2,J)-P(I+1,J)) &
        -nu_vis(I-1,J)* &
        1.0_SP/DX*(P(I,J)-P(I-1,J)) ) &
        + 1.0_SP/DY*( &
        0.5_SP*(nu_vis(I,J+1)+nu_vis(I,J))* &
        0.5_SP/DY*(P(I,J+1)+P(I+1,J+1)-P(I,J)-P(I+1,J)) &
        -0.5_SP*(nu_vis(I,J-1)+nu_vis(I,J))* &
        0.5_SP/DY*(P(I,J)+P(I+1,J)-P(I,J-1)-P(I+1,J-1)) )

        tempY_vis(I,J) = 0.5_SP/DY*( &
        nu_vis(I,J+1)* &
        1.0_SP/DY*(Q(I,J+2)-Q(I,J+1)) &
        -nu_vis(I,J-1)* &
        1.0_SP/DY*(Q(I,J)-Q(I,J-1)) ) &
        + 1.0_SP/DX*( &
        0.5_SP*(nu_vis(I+1,J)+nu_vis(I,J))* &
        0.5_SP/DX*(Q(I+1,J)+Q(I+1,J+1)-Q(I,J)-Q(I,J+1)) &
        -0.5_SP*(nu_vis(I-1,J)+nu_vis(I,J))* &
        0.5_SP/DX*(Q(I,J)+Q(I,J+1)-Q(I-1,J)-Q(I-1,J+1)) )

      ELSE

        tempX_vis(I,J) = 0.5_SP/DX*( &
        (nu_vis(I+1,J)+nu_vis(I,J)) &
        *1.0_SP/DX*(HU(I+1,J)-HU(I,J)) &
        -(nu_vis(I-1,J)+nu_vis(I,J)) &
        *1.0_SP/DX*(HU(I,J)-HU(I-1,J)) ) &
        + 0.5_SP/DY*( &
        (nu_vis(I,J+1)+nu_vis(I,J)) &
        *1.0_SP/DY*(HU(I,J+1)-HU(I,J)) &
        -(nu_vis(I,J-1)+nu_vis(I,J)) &
        *1.0_SP/DY*(HU(I,J)-HU(I,J-1)) )

        tempY_vis(I,J) = 0.5_SP/DX*( &
        (nu_vis(I+1,J)+nu_vis(I,J)) &
        *1.0_SP/DX*(HV(I+1,J)-HV(I,J)) &
        -(nu_vis(I-1,J)+nu_vis(I,J)) &
        *1.0_SP/DX*(HV(I,J)-HV(I-1,J)) ) &
        + 0.5_SP/DY*( &
        (nu_vis(I,J+1)+nu_vis(I,J)) &
        *1.0_SP/DY*(HV(I,J+1)-HV(I,J)) &
        -(nu_vis(I,J-1)+nu_vis(I,J)) &
        *1.0_SP/DY*(HV(I,J)-HV(I,J-1)) )
      ENDIF ! end pq_scheme
    ENDDO
  ENDDO

  ! Mixing
# endif
  ! Cartesian
# endif


  ! depth gradient term
  DO J=Jbeg,Jend
    DO I=Ibeg,Iend
# if defined (CARTESIAN)
      ! second order, move the second term to left-hand side
      SourceX(I,J)=GRAV*(Eta(I,J))/DX*(Depthx(I+1,J)-Depthx(I,J)) &
      *MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))   &
      *U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
      ! end Manning
# endif
      ! dispersion
      ! (h+eta)(u*\nabla V4 + V4 * \nabla u - v1pp-v2-v3)
      + Gamma1*MASK9(I,J)*Max(H(I,J),MinDepthFrc)*(  &
      U(I,J)*0.5_SP*(U4(I+1,J)-U4(I-1,J))/DX+V(I,J)  &
      *0.5_SP*(U4(I,J+1)-U4(I,J-1))/DY &
      +U4(I,J)*0.5_SP*(U(I+1,J)-U(I-1,J))/DX+V4(I,J) &
      *0.5_SP*(U(I,J+1)-U(I,J-1))/DY &
      -Gamma2*MASK9(I,J)*(U1pp(I,J)+U2(I,J)+U3(I,J)) &
      )&
      ! Ht(-V4+V1p) = div(M)*(U4-U1p)
      +Gamma1*MASK9(I,J)* &
      ((P(I+1,J)-P(I,J))/DX+(Q(I,J+1)-Q(I,J))/DY)*(U4(I,J)-U1p(I,J)) &
      ! wavemaker
      +WaveMaker_Mass(I,J)*U(I,J)

      IF(FRICTION_SPONGE) THEN
        ! note that, compared to wei et al, we used flux.
        ! So need multiply D
        SourceX(I,J) = SourceX(I,J) &
        -CD_4_SPONGE(I,J)*U(I,J) &
        *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*Depth(I,J)
      ENDIF

      IF(WaveMakerCurrentBalance)THEN
# if defined (PARALLEL)
        xmk=(I-Ibeg)*DXg + (iista-1)*DXg
        ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        xmk=(I-Ibeg)*DXg
        ymk=(J-Jbeg)*DYg
# endif
        IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
        ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          SourceX(I,J) = SourceX(I,J) &
          -WaveMakerCd*U(I,J) &
          *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        ENDIF
      ENDIF ! current balance

      IF(BREAKWATER) THEN
        ! note that, compared to wei et al, we used flux.
        ! So need multiply D
        SourceX(I,J) = SourceX(I,J) &
        - CD_breakwater(I,J)*U(I,J) &
        *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*Depth(I,J)
      ENDIF

      SourceY(I,J)=GRAV*(Eta(I,J))/DY*(Depthy(I,J+1)-Depthy(I,J)) &
      *MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))  &
      *V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
      ! end Manning
# endif
      ! dispersion
      ! (h+eta)(u*\nabla V4 + V4 * \nabla u -v1pp-v2-v3)
      + Gamma1*MASK9(I,J)*Max(H(I,J),MinDepthFrc)*( &
      U(I,J)*0.5_SP*(V4(I+1,J)-V4(I-1,J))/DX+V(I,J) &
      *0.5_SP*(V4(I,J+1)-V4(I,J-1))/DY &
      +U4(I,J)*0.5_SP*(V(I+1,J)-V(I-1,J))/DX+V4(I,J) &
      *0.5_SP*(V(I,J+1)-V(I,J-1))/DY  &
      -Gamma2*MASK9(I,J)*(V1pp(I,J)+V2(I,J)+V3(I,J)) &
      )&
      ! Ht(-V4+V1p) = div(Q)*(V4-V1p)
      +Gamma1*MASK9(I,J) &
      *((P(I+1,J)-P(I,J))/DX+(Q(I,J+1)-Q(I,J))/DY) &
      *(V4(I,J)-V1p(I,J))  &
      ! wavemaker
      +WaveMaker_Mass(I,J)*V(I,J)

      IF(FRICTION_SPONGE) THEN
        ! note that, compared to wei et al, we used flux.
        ! So need multiply D
        SourceY(I,J) = SourceY(I,J) &
        -CD_4_SPONGE(I,J)*V(I,J) &
        *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*Depth(I,J)
      ENDIF

      IF(WaveMakerCurrentBalance)THEN
# if defined (PARALLEL)
        xmk=(I-Ibeg)*DXg + (iista-1)*DXg
        ymk=(J-Jbeg)*DYg + (jjsta-1)*DYg
# else
        xmk=(I-Ibeg)*DXg
        ymk=(J-Jbeg)*DYg
# endif
        IF(ABS(xmk-Xc_WK)<Width_WK.AND. &
          ABS(ymk-Yc_WK)<Ywidth_WK/2.0_SP)THEN
          SourceY(I,J) = SourceY(I,J) &
          -WaveMakerCd*V(I,J) &
          *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        ENDIF
      ENDIF ! current balance

      IF(BREAKWATER) THEN
        ! note that, compared to wei et al, we used flux.
        ! So need multiply D
        SourceY(I,J) = SourceY(I,J) &
        -CD_BREAKWATER(I,J)*V(I,J) &
        *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*Depth(I,J)
      ENDIF

# if defined (METEO)
      IF(WindForce)THEN
        SourceX(I,J)=SourceX(I,J) &
        +MASK_WIND(I,J)*RHO_AW*Cdw*WindU2D(I,J) &
        *SQRT(WindU2D(I,J)*WindU2D(I,J) &
        + WindV2D(I,J)*WindV2D(I,J))

        SourceY(I,J)=SourceY(I,J) &
        +MASK_WIND(I,J)*RHO_AW*Cdw*WindV2D(I,J) &
        *SQRT(WindU2D(I,J)*WindU2D(I,J) &
        + WindV2D(I,J)*WindV2D(I,J))
      ENDIF   ! end wind
      ! end METEO
# endif

# if defined (MIXING)
      SourceX(I,J) = SourceX(I,J) ! + tempX(I,J)
      SourceY(I,J) = SourceY(I,J) ! + tempY(I,J)
# endif

      ! eddy viscosity breaking
      ! ykchoi(08.18.2015)
      ! new variable :: WAVEMAKER_VIS
      !IF(VISCOSITY_BREAKING.OR.DIFFUSION_SPONGE)THEN
      IF(VISCOSITY_BREAKING.OR.DIFFUSION_SPONGE.OR.WAVEMAKER_VIS)THEN
        SourceX(I,J) = SourceX(I,J) + tempX_vis(I,J)
        SourceY(I,J) = SourceY(I,J) + tempY_vis(I,J)
      ENDIF  ! end eddy viscosity breaking
      ! else Cartesian
# else

      ! second order, move the second term to left-hand side
# if defined (ZALPHA)
      SourceX(I,J)=GRAV*(Eta(I,J))*SlopeX(I,J)*MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))  &
      *U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# endif
      ! dispersion
      ! Ht(+V1p) = div(M)*(-U1p)
      +Gamma1*MASK9(I,J) &
      *((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
      *(U4(I,J)-U1p(I,J)) &
      ! Coriolis
      +Coriolis(I,J)*0.5_SP*(Q(I,J)+Q(I,J+1))

      SourceY(I,J)=GRAV*(Eta(I,J))*SlopeY(I,J)*MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))  &
      *V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# endif
      ! dispersion
      ! Ht(+V1p) = div(Q)*(-V1p)
      +Gamma1*MASK9(I,J) &
      *((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
      *(V4(I,J)-V1p(I,J)) &
      ! Coriolis
      -Coriolis(I,J)*0.5_SP*(P(I,J)+P(I+1,J))
      ! else ZALPHA
# else
      SourceX(I,J)=GRAV*(Eta(I,J))*SlopeX(I,J)*MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))  &
      *U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*U(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# endif
      ! dispersion
      ! Ht(+V1p) = div(M)*(-U1p)
      +Gamma1*MASK9(I,J) &
      *((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
      *(-U1p(I,J)) &
      ! Coriolis
      +Coriolis(I,J)*0.5_SP*(Q(I,J)+Q(I,J+1))

      SourceY(I,J)=GRAV*(Eta(I,J))*SlopeY(I,J)*MASK(I,J) &
      ! friction
# if defined (MANNING)
      -(GRAV*Cd(I,J)**2/  &
      Max(H(I,J),MinDepthFrc)**(0.333333_SP))  &
      *V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# else
      -Cd(I,J)*V(I,J)*SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)) &
# endif
      ! dispersion
      ! Ht(+V1p) = div(Q)*(-V1p)
      +Gamma1*MASK9(I,J) &
      *((P(I+1,J)-P(I,J))/DX(I,J)+(Q(I,J+1)-Q(I,J))/DY(I,J)) &
      *(-V1p(I,J)) &
      ! Coriolis
      -Coriolis(I,J)*0.5_SP*(P(I,J)+P(I+1,J))
      ! end if ZALPHA
# endif
      ! end if Cartesian
# endif

    ENDDO
  ENDDO

# if defined (VESSEL)
  SourceX = SourceX + VesselPressureX
  SourceY = SourceY + VesselPressureY
# endif

# if defined (METEO)
  IF(AirPressure)THEN
    SourceX = SourceX + StormPressureX
    SourceY = SourceY + StormPressureY
  ENDIF
# endif

END SUBROUTINE SourceTerms
