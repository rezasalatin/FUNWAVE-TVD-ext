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
  LOGICAL :: PQ_scheme = .TRUE.
  REAL(SP) :: xmk,ymk,Dxg,Dyg,WK_Source
  INTEGER :: kd,kf
  REAL(SP) :: AA                          ! for all WK's (optimized)
  REAL(SP), DIMENSION(Nfreq) :: BB,CC     ! for WK_IRR (optimized)
  REAL(SP), DIMENSION(NumWaveComp,1) :: BB1 ! for WK_TIME (optimized)
  REAL(SP), DIMENSION(NumFreq) :: BB2,CC2 ! for WK_DATA (optimized)
  INTEGER :: kf_index ! Chen, Salatin
  ! variables eddy viscosity terms (horizontal, vertical, breaking)
  REAL(SP),DIMENSION(Mloc,Nloc) :: Temp1,Temp2,Temp3,Temp4
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
  ! two regular waves out of phase
  IF(WAVEMAKER(1:6)=='WK_REG')THEN
    AA = TANH(PI/(Time_ramp*Tperiod)*TIME)*D_gen

    IF(10.EQ.0)THEN
      ! do not execute this part (original code)
      DO J=jlo,jhi
        DO I=ilo,ihi
          WaveMaker_Mass(I,J)= AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
          *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME)
        ENDDO
      ENDDO
    ELSEIF(0.EQ.0)THEN
    ! SPLIT WAVEMAKER INTO TWO with 0.5 m gap in between according to
    ! Choi 2021 experimental design
    ! first part has phase of pi and second has zero deg phase
      DO J=Jbeg,Jend
        IF(ymk_wk(J).LE.(Yc_WK-2.0*DY)) THEN
          DO I=ilo,ihi
            WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
            *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME+PI)
          ENDDO
        ELSEIF(ymk_wk(J).GE.(Yc_WK+2.0*DY)) THEN
          DO I=ilo,ihi
            WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
            *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME)
          ENDDO
        ENDIF
      ENDDO
    ELSEIF(30.EQ.0)THEN
      ! split into four parts
      DO J=Jbeg,Jend
        IF(ymk_wk(J).GT.(Yc_WK/3.0+1.0*DY).AND.ymk_wk(J).LT.(Yc_WK-1.0*DY)) THEN
          DO I=ilo,ihi
          WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
          *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME+PI)
          ENDDO
        ELSEIF(ymk_wk(J).GT.(Yc_WK+1.0*DY).AND.ymk_wk(J).LT.(Yc_WK*5.0/3.0-1.0*DY)) THEN
          DO I=ilo,ihi
            WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
            *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME)
          ENDDO
        ELSEIF(ymk_wk(J).LT.(Yc_WK/3.0-1.0*DY)) THEN
          DO I=ilo,ihi
            WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
            *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME)
          ENDDO
        ELSEIF(ymk_wk(J).GT.(Yc_WK*5.0/3.0+1.0*DY)) THEN
          DO I=ilo,ihi
            WaveMaker_Mass(I,J)=AA*EXP(-Beta_gen*(xmk_wk(I)-Xc_WK)**2) &
            *SIN(rlamda*(ymk_wk(J)-ZERO)-2.0*PI/Tperiod*TIME+PI)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
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
  ! Measured Spectra WKs
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
  ! Measured Spectra WKs
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
    nu_vis = nu_vis + nu_break
  ENDIF

  IF(DIFFUSION_SPONGE)THEN
    ! should include ghost cells 02/28/2017, fyshi
    nu_vis = nu_vis + nu_sponge
  ENDIF

107   format(500f12.6)

# if defined (CARTESIAN)
  DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      ! depth gradient term
      ! second order, move the second term to left-hand side

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Source term in X direction
      SourceX(I,J)=GRAV*Eta(I,J)*(Depthx(I+1,J)-Depthx(I,J))/DX &
      *MASK(I,J) &
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Source term in Y direction
      ! Equation 22 in Shi et al. 2012
      SourceY(I,J)=GRAV*Eta(I,J)*(Depthy(I,J+1)-Depthy(I,J))/DY &
      *MASK(I,J) &
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! CURRENT BALANCE
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

          SourceY(I,J) = SourceY(I,J) &
          -WaveMakerCd*V(I,J) &
          *SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        ENDIF
      ENDIF ! current balance
    ENDDO
  ENDDO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Bottom Friction
# if defined (MANNING)
  Temp1 = -(GRAV*Cd**2/Max(H,MinDepthFrc)**(0.333333_SP))*U*SQRT(U*U+V*V)
  Temp2 = -(GRAV*Cd**2/Max(H,MinDepthFrc)**(0.333333_SP))*V*SQRT(U*U+V*V)
# else
  Temp1 = -Cd*U*SQRT(U*U+V*V)
  Temp2 = -Cd*V*SQRT(U*U+V*V)
  ! end Manning
# endif
  SourceX = SourceX + Temp1
  SourceY = SourceY + Temp2
  VTYPE=1
  CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp1,VTYPE,PERIODIC)
  CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp2,VTYPE,PERIODIC)
  ! bottom friction (vorticity source/sink)
  temp_output1 = Temp1  ! R_f with units m2/s2
  temp_output2 = Temp2
  Temp1 = ZERO
  Temp2 = ZERO


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! FRICTION SPONGE

  IF(FRICTION_SPONGE) THEN
    ! note that, compared to wei et al, we used flux.
    ! So need multiply D
    SourceX = SourceX &
    -CD_4_SPONGE*U*SQRT(U*U+V*V)*Depth
    SourceY = SourceY &
    -CD_4_SPONGE*V*SQRT(U*U+V*V)*Depth
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! BREAKWATER

  IF(BREAKWATER) THEN
    ! note that, compared to wei et al, we used flux.
    ! So need multiply D
    SourceX = SourceX &
    -CD_breakwater*U*SQRT(U*U+V*V)*Depth
    SourceY = SourceY &
    -CD_BREAKWATER*V*SQRT(U*U+V*V)*Depth
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! METEO MODULE

# if defined (METEO)
  IF(WindForce)THEN
    SourceX=SourceX &
    +MASK_WIND*RHO_AW &
    *Cdw*WindU2D*SQRT(WindU2D*WindU2D + WindV2D*WindV2D)

    SourceY=SourceY &
    +MASK_WIND*RHO_AW &
    *Cdw*WindV2D*SQRT(WindU2D*WindU2D + WindV2D*WindV2D)
  ENDIF
# endif

  ! breaking
  DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      Temp1(I,J) = 1.0_SP * (1.0_SP/DX*( &
      0.5_SP*(nu_vis(I,J)+nu_vis(I+1,J)) &
      *1.0_SP/DX*(P(I+1,J)-P(I,J))  &
      -0.5_SP*(nu_vis(I-1,J)+nu_vis(I,J)) &
      *1.0_SP/DX*(P(I,J)-P(I-1,J))) ) &
      + 0.5_SP * ( 1.0_SP/DY*( &
      0.5_SP*(nu_vis(I,J)+nu_vis(I,J+1)) &
      *1.0_SP/DY*(P(I,J+1)-P(I,J)) &
      -0.5_SP*(nu_vis(I,J-1)+nu_vis(I,J)) &
      *1.0_SP/DY*(P(I,J)-P(I,J-1))) &
      + 0.5_SP/DY*( &
      nu_vis(I,J+1) &
      *0.5_SP/DX*(Q(I+1,J+1)-Q(I-1,J+1)) &
      -nu_vis(I,J-1) &
      *0.5_SP/DX*(Q(I+1,J-1)-Q(I-1,J-1))) )

      Temp2(I,J) = 1.0_SP * (1.0_SP/DY*( &
      0.5_SP*(nu_vis(I,J)+nu_vis(I,J+1)) &
      *1.0_SP/DY*(Q(I,J+1)-Q(I,J))  &
      - 0.5_SP*(nu_vis(I,J-1)+nu_vis(I,J)) &
      *1.0_SP/DY*(Q(I,J)-Q(I,J-1))) ) &
      + 0.5_SP * ( 1.0_SP/DX*( &
      0.5_SP*(nu_vis(I,J)+nu_vis(I+1,J)) &
      *1.0_SP/DX*(Q(I+1,J)-Q(I,J)) &
      -0.5_SP*(nu_vis(I-1,J)+nu_vis(I,J)) &
      *1.0_SP/DX*(Q(I,J)-Q(I-1,J))) &
      + 0.5_SP/DX*( &
      nu_vis(I+1,J) &
      *0.5_SP/DY*(P(I+1,J+1)-P(I+1,J-1)) &
      -nu_vis(I-1,J)* &
      0.5_SP/DY*(P(I-1,J+1)-P(I-1,J-1))) )
    ENDDO
  ENDDO

  SourceX = SourceX + Temp1
  SourceY = SourceY + Temp2

  brksrcX_Window = CSHIFT(brksrcX_Window,-1,3)
  brksrcY_Window = CSHIFT(brksrcY_Window,-1,3)
  ! replace oldest variable by the newest
  brksrcX_Window(:,:,1) = Temp1
  brksrcY_Window(:,:,1) = Temp2
  ! do averaging
  Temp1 = SUM(brksrcX_Window,3)/MA_Window
  Temp2 = SUM(brksrcY_Window,3)/MA_Window

  VTYPE=1
  CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp1,VTYPE,PERIODIC)
  CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp2,VTYPE,PERIODIC)
  ! breaking (vorticity source)
  temp_output3 = Temp1
  temp_output4 = Temp2
  Temp1 = ZERO
  Temp2 = ZERO

  ! Horizontal Mixing
  IF(H_MIXING)THEN
    IF(MIX_MODEL.EQ.0)THEN
      ! default
      ! nu_smg = C_smg**2*DX*DY*SQRT(2.0*uMAx**2+2.0*vMAy**2+(uMAy+vMAx)**2)
      nu_smg = C_smg*DX*DY*SQRT(uMAx**2+vMAy**2+((uMAy+vMAx)**2)/2.0)
    ELSEIF(MIX_MODEL.EQ.1)THEN
    ! Leith
      VTYPE=1
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,vorticity,VTYPE,PERIODIC)
      CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,vorticity,Temp1)
      CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,vorticity,Temp2)
      Temp3 = Temp1**2+Temp2**2
      Temp1 = ZERO
      Temp2 = ZERO
      ! Lambda (C_smg) is 1 for Leith et al
      nu_smg = (C_smg*SQRT(DX*DY)/pi)**3*SQRT(Temp3)
    ELSEIF(MIX_MODEL.EQ.6)THEN
    ! Fox-Kemper
      VTYPE=1
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,vorticity,VTYPE,PERIODIC)
      CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,vorticity,Temp1)
      CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,vorticity,Temp2)
      Temp3 = Temp1**2+Temp2**2
      Temp1 = ZERO
      Temp2 = ZERO
      CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,uMAx+vMAy,Temp1)
      CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,uMAx+vMAy,Temp2)
      !temp_1_output = Temp3
      !temp_1_output = Temp1**2+Temp2**2
      ! Lambda (C_smg) is 1 for Fox-Kemper
      nu_smg = (C_smg*SQRT(DX*DY)/pi)**3*SQRT(Temp3+Temp1**2+Temp2**2)
      Temp1 = ZERO
      Temp2 = ZERO
    ELSEIF(MIX_MODEL.EQ.61)THEN
      VTYPE=1
      CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,vorticity,VTYPE,PERIODIC)
      CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,vorticity,Temp1)
      CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,vorticity,Temp2)
      Temp3 = Temp1**2+Temp2**2
      Temp1 = ZERO
      Temp2 = ZERO
      CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DX,Ux+Vy,Temp1)
      CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK,DY,Ux+Vy,Temp2)
      !temp_1_output = Temp3
      !temp_1_output = Temp1**2+Temp2**2
      ! Lambda (C_smg) is 1 for Fox-Kemper
      nu_smg = (C_smg*SQRT(DX*DY)/pi)**3*SQRT(Temp3+Temp1**2+Temp2**2)
      Temp1 = ZERO
      Temp2 = ZERO
    ENDIF

    ! when there is breaking, deactivate the horizontal mixing
    nu_smg = nu_smg*(1-CEILING(nu_vis))

    DO J=Jbeg,Jend
      DO I=Ibeg,Iend
        Temp1(I,J) = 1.0_SP * (1.0_SP/DX*( &
        0.5_SP*(nu_smg(I,J)+nu_smg(I+1,J)) &
        *1.0_SP/DX*(P(I+1,J)-P(I,J))  &
        -0.5_SP*(nu_smg(I-1,J)+nu_smg(I,J)) &
        *1.0_SP/DX*(P(I,J)-P(I-1,J))) ) &
        + 0.5_SP * ( 1.0_SP/DY*( &
        0.5_SP*(nu_smg(I,J)+nu_smg(I,J+1)) &
        *1.0_SP/DY*(P(I,J+1)-P(I,J)) &
        -0.5_SP*(nu_smg(I,J-1)+nu_smg(I,J)) &
        *1.0_SP/DY*(P(I,J)-P(I,J-1))) &
        + 0.5_SP/DY*( &
        nu_smg(I,J+1) &
        *0.5_SP/DX*(Q(I+1,J+1)-Q(I-1,J+1)) &
        -nu_smg(I,J-1) &
        *0.5_SP/DX*(Q(I+1,J-1)-Q(I-1,J-1))) )

        Temp2(I,J) = 1.0_SP * (1.0_SP/DY*( &
        0.5_SP*(nu_smg(I,J)+nu_smg(I,J+1)) &
        *1.0_SP/DY*(Q(I,J+1)-Q(I,J))  &
        - 0.5_SP*(nu_smg(I,J-1)+nu_smg(I,J)) &
        *1.0_SP/DY*(Q(I,J)-Q(I,J-1))) ) &
        + 0.5_SP * ( 1.0_SP/DX*( &
        0.5_SP*(nu_smg(I,J)+nu_smg(I+1,J)) &
        *1.0_SP/DX*(Q(I+1,J)-Q(I,J)) &
        -0.5_SP*(nu_smg(I-1,J)+nu_smg(I,J)) &
        *1.0_SP/DX*(Q(I,J)-Q(I-1,J))) &
        + 0.5_SP/DX*( &
        nu_smg(I+1,J) &
        *0.5_SP/DY*(P(I+1,J+1)-P(I+1,J-1)) &
        -nu_smg(I-1,J)* &
        0.5_SP/DY*(P(I-1,J+1)-P(I-1,J-1))) )
      ENDDO
    ENDDO

    SourceX = SourceX + Temp1
    SourceY = SourceY + Temp2
    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp1,VTYPE,PERIODIC)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Temp2,VTYPE,PERIODIC)
    ! mixing (vorticity source/sink)
    temp_output5 = Temp1
    temp_output6 = Temp2
    Temp1 = ZERO
    Temp2 = ZERO

  ELSE
    nu_smg = ZERO
  ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! else Cartesian
# else

  DO J=Jbeg,Jend
    DO I=Ibeg,Iend

      ! second order, move the second term to left-hand side
# if defined (ZALPHA)
      SourceX(I,J)=GRAV*Eta(I,J)*SlopeX(I,J)*MASK(I,J) &
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

      SourceY(I,J)=GRAV*Eta(I,J)*SlopeY(I,J)*MASK(I,J) &
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
    ENDDO
  ENDDO
  ! end if Cartesian
# endif


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

! calculate U1 and V1 
! Chen, Salatin
U1 = (U1p - U1p_Old)/DT + U1pp
V1 = (V1p - V1p_Old)/DT + V1pp
U1p_Old = U1p
V1p_Old = V1p
! calculate VorT
VorT = (vorticity - Vor_Old)/DT
Vor_Old = vorticity

END SUBROUTINE SourceTerms

! filtering based on FUNWAVE 1.0 code to smooth out nu_vis
SUBROUTINE fltr(var_fltr)
  USE GLOBAL

  IMPLICIT NONE
  REAL,DIMENSION(Mloc,Nloc),INTENT(INOUT) :: var_fltr
  REAL,DIMENSION(Mloc,Nloc) :: var_fltr_tmp
  REAL :: c1, c2, c3, c4, c5

  c1 = 186./256.
  c2 =  56./256.
  c3 = -28./256.
  c4 =   8./256.
  c5 =  -1./256.

  ! x direction

  DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      var_fltr_tmp(i,j)=c1*var_fltr(i,j) &
      +c2*(var_fltr(i+1,j)+var_fltr(i-1,j)) &
      +c3*(var_fltr(i+2,j)+var_fltr(i-2,j)) &
      +c4*(var_fltr(i+3,j)+var_fltr(i-3,j)) &
      +c5*(var_fltr(i+4,j)+var_fltr(i-4,j))
    ENDDO
  ENDDO

  DO J=Jbeg,Jend
    DO I=Ibeg,Iend
      var_fltr_tmp(i,j)=c1*var_fltr_tmp(i,j) &
      +c2*(var_fltr_tmp(i,j+1)+var_fltr_tmp(i,j-1)) &
      +c3*(var_fltr_tmp(i,j+2)+var_fltr_tmp(i,j-2)) &
      +c4*(var_fltr_tmp(i,j+3)+var_fltr_tmp(i,j-3)) &
      +c5*(var_fltr_tmp(i,j+4)+var_fltr_tmp(i,j-4))
    ENDDO
  ENDDO

  var_fltr = var_fltr_tmp

END SUBROUTINE fltr
