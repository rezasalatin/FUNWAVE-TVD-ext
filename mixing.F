!------------------------------------------------------------------------------------
!
!      FILE mixing.F
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
!-------------------------------------------------------------------------------------
!
!    MIXING_STUFF is subroutine to calculate mixing related, time-averaged properties
!    mean eta is also calculated.
!
!    HISTORY: 05/02/2011 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE MIXING_STUFF
  USE GLOBAL
  IMPLICIT NONE

  ! calculate mean for smagorinsky s mixing and wave height
  !ykchoi (be careful)
  !I think Umean, Vmean is not using in other routine.
  IF( time >= STEADY_TIME )THEN
    CALL CALCULATE_MEAN(T_INTV_mean,T_INTV_UVmean,T_sum,T_sumUVmean,&
    DT,Mloc,Nloc,U,V,ETA,ETA0,&
    Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,WaveHeightRMS, &
    WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
    HrmsSum,HavgSum,ETA2,ETA2sum,SigWaveHeight,ETAvariance, &
    Vor_Instant, Vor_Mean, Vor_Intensity, NumDT1, NumDT2,Vor2,&
    Vor2sum,Vorsum)
  ENDIF    !ykchoi

  ! Chen, Salatin, moving average of some properties (U,V,Eta,Vor,Syx)
  CALL CALCULATE_MovingAverage(DT,Mloc,Nloc,U,V,ETA,MASK9,DX,DY,Depth,&
  Ibeg,Iend,Jbeg,Jend,&
  OUT_uvMA,OUT_VorMA,OUT_RadStr,OUT_etaMA,&
  MA_Window,T_MA_Sum,T_MA_STEP,&
  Syx_temp,Syx_x,U_Window,V_Window,uMA,vMA,eta_MA_temp,eta_MA,Vor_MA,&
  ! Smagorinsky subgrid scale turbulence calculated here
  PERIODIC,Nghost,uMAx,uMAy,vMAx,vMAy,vis_smg_xMA_temp,vis_smg_yMA_temp,&
  vis_smg_xMA,vis_smg_yMA,vis_smg_x,vis_smg_y,nu_smg,nu_smg_MA,nu_smg_MA_temp)

END SUBROUTINE MIXING_STUFF


!-------------------------------------------------------------------------------------
!
!    CALCULATE_MEAN is subroutine to calculate mean u v required by
!      smagorinsky mixing and wave height
!      mean eta is also calculated.
!
!    HISTORY:
!      05/02/2011 Fengyan Shi
!                 Young-Kwang Choi added some time-averaging stuff
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_MEAN(T_INTV_mean,T_INTV_UVmean,T_sum,T_sumUVmean,&
  DT,M,N,U,V,ETA,ETA0,&
  Umean,Vmean,ETAmean,Usum,Vsum,ETAsum,&
  WaveHeightRMS, &
  WaveHeightAve,Emax,Emin,Num_Zero_Up,Ibeg,Iend,Jbeg,Jend, &
  HrmsSum,HavgSum,ETA2,ETA2sum,SigWaveHeight,ETAvariance, &
  Vor_Instant, Vor_Mean, Vor_Intensity, NumDT1, NumDT2,Vor2, &
  Vor2sum,Vorsum)

  ! calculate mean for smagorinsky mixing and wave height
  USE PARAM
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
  REAL(SP),DIMENSION(M,N),INTENT(IN)::U,V,ETA,ETA0
  REAL(SP),INTENT(IN) :: T_INTV_mean,T_INTV_UVmean,DT
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Umean,Vmean
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: WaveHeightRMS,WaveHeightAve
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Usum,Vsum
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: HrmsSum,HavgSum
  REAL(SP),INTENT(OUT) :: T_sum,T_sumUVmean
  REAL(SP)::Tmpe,Tmp_0
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Emax,Emin
  INTEGER,DIMENSION(M,N),INTENT(INOUT) :: Num_Zero_Up
  !!! Chen, Salatin
  INTEGER,INTENT(INOUT) :: NumDT1, NumDT2
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Vor_Instant,Vor_Mean,Vor_Intensity
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: Vor2,Vor2sum,Vorsum
  REAL(SP),DIMENSION(M,N) :: Vor_Variance
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: ETAmean,ETAsum,ETA2,ETA2sum
  REAL(SP),DIMENSION(M,N) :: ETAvariance
  REAL(SP),DIMENSION(M,N),INTENT(INOUT) :: SigWaveHeight

  T_sum=T_sum+DT ! for wave height and MSL
  T_sumUVmean=T_sumUVmean+DT ! Chen, Salatin, for current and vortcty field

  IF(T_sum.GE.T_INTV_mean)THEN
    ! Significant Wave Height
    ! Chen, Salatin
    ! variance equation with initially unknown etamean
    ETAsum =ETA*DT+ETAsum
    ETAmean=ETAsum/T_sum
    DO J=1,N
      DO I=1,M
        ETAvariance(I,J) = 1/real(NumDT1-1)*ETA2sum(I,J)&
        -ETAmean(I,J)*ETAmean(I,J)*real(NumDT1)/real(NumDT1-1)
      ENDDO
    ENDDO
    SigWaveHeight = 4.0*SQRT(ETAvariance)

    ETAsum=ZERO
    ETA2sum=ZERO
    NumDT1=ZERO
    T_sum=T_sum-T_INTV_mean

    ! Average and RMS wave height
    DO J=1,N
      DO I=1,M
        IF(Num_Zero_Up(I,J)>=2)THEN
          WaveHeightAve(I,J)=HavgSum(I,J)/Num_Zero_Up(I,J)
          WaveHeightRMS(I,J)=SQRT(HrmsSum(I,J)/Num_Zero_Up(I,J))
        ENDIF
        !        Num_Zero_Up(I,J)=ZERO
        !        HavgSum(I,J)=ZERO
        !        HrmsSum(I,J)=ZERO
      ENDDO
    ENDDO

    CALL PREVIEW_MEAN
  ELSE

    NumDT1 = NumDT1+1
    ETAsum=ETA*DT+ETAsum
    DO J=1,N
      DO I=1,M
        ETA2(I,J) = ETA(I,J)*ETA(I,J)
      ENDDO
    ENDDO
    ETA2sum=ETA2+ETA2sum

    ! wave height
    ! Chen, Salatin
    ETAmean = ETAsum/T_sum ! effective for long simulations,
    ! for short simulation just comment it
    DO J=1,N
      DO I=1,M
        if(Eta(i,j)>Emax(i,j)) Emax(i,j) = Eta(i,j)
        if(Eta(i,j)<Emin(i,j)) Emin(i,j) = Eta(i,j)
        Tmpe = Eta(i,j)-ETAmean(i,j)
        Tmp_0 = Eta0(i,j)-ETAmean(i,j)
        if(Tmpe>Tmp_0.and.Tmpe*Tmp_0<=Zero) then
          Num_Zero_Up(i,j) = Num_Zero_Up(i,j)+1
          if(Num_Zero_Up(i,j)>=2) then
            HavgSum(i,j) = HavgSum(i,j)+Emax(i,j)-Emin(i,j)
            HrmsSum(i,j) = HrmsSum(i,j)+(Emax(i,j)-Emin(i,j))**2
          endif
          ! reset Emax and Emin to find next wave
          Emax(i,j) = -1000.
          Emin(i,j) =  1000.
        endif
      ENDDO
    ENDDO

  ENDIF  ! end average time

  !-------------------------------------------------------------------------------------
  !
  !    PREVIEW_UVMEAN is subroutine for print-out of mean U and V field data
  !
  !  HISTORY:
  !    05/13/2020  Chen, Salatin UVmean
  !    01/25/2020  Chen, Salatin Vorticity
  !
  ! since U and V are preferably averaged over 2Tp,
  ! a separate function is assigned for that
  !-------------------------------------------------------------------------------------

  IF(T_sumUVmean.GE.T_INTV_UVmean)THEN

    Usum=U*DT+Usum
    Vsum=V*DT+Vsum
    Umean=Usum/T_sumUVmean
    Vmean=Vsum/T_sumUVmean

    !Vorsum =Vor_Instant+Vorsum
    Vorsum =Vor_Instant*DT+Vorsum
    !Vor_Mean=Vorsum/NumDT2
    Vor_Mean=Vorsum/T_sumUVmean
    DO J=1,N
      DO I=1,M
        Vor_Variance(I,J) = 1/real(NumDT2-1)*Vor2sum(I,J)&
        -Vor_Mean(I,J)*Vor_Mean(I,J)*real(NumDT2)/real(NumDT2-1)
        Vor_Intensity(I,J) = SQRT(Vor_Variance(I,J))
      ENDDO
    ENDDO
    T_sumUVmean=T_sumUVmean-T_INTV_UVmean
    Usum=ZERO
    Vsum=ZERO
    Vorsum=ZERO
    Vor2sum=ZERO
    NumDT2=ZERO
    CALL PREVIEW_MEAN_UV

  ELSE

    Usum=U*DT+Usum
    Vsum=V*DT+Vsum
    !Vorticity Intensity
    NumDT2 = NumDT2+1
    !Vorsum =Vor_Instant+Vorsum
    Vorsum =Vor_Instant*DT+Vorsum
    DO J=1,N
      DO I=1,M
        Vor2(I,J) = Vor_Instant(I,J)*Vor_Instant(I,J)
      ENDDO
    ENDDO
    Vor2sum=Vor2+Vor2sum
  ENDIF

END SUBROUTINE CALCULATE_MEAN

!-------------------------------------------------------------------------------------
!
! Chen, Salatin, Oct 2021, Moving Average of some properties
! a moving average is carried out on instant velocity
! the window size is 2Tp for regular waves or 5Tp for irregular waves
! this window moves at each specified timestep
!
! Moving average of U, V, Vorticity, Eta, Radiation Shear Stress
!
!-------------------------------------------------------------------------------------
SUBROUTINE CALCULATE_MovingAverage(DT,M,N,U,V,ETA,MASK,DX,DY,Depth,&
  Ibeg,Iend,Jbeg,Jend, &
  OUT_uvMA,OUT_VorMA,OUT_RadStr,OUT_etaMA,&
  MA_Window,T_MA_Sum,T_MA_STEP,&
  Syx_temp,Syx_x,U_Window,V_Window,uMA,vMA,eta_MA_temp,eta_MA,Vor_MA,&
  PERIODIC,Nghost,uMAx,uMAy,vMAx,vMAy,vis_smg_xMA_temp,vis_smg_yMA_temp,&
  vis_smg_xMA,vis_smg_yMA,vis_smg_x,vis_smg_y,nu_smg,nu_smg_MA,nu_smg_MA_temp)

  USE PARAM
  IMPLICIT NONE
  REAL(SP),INTENT(IN) :: DT
  INTEGER, INTENT(IN) :: M,N,Ibeg,Iend,Jbeg,Jend
  REAL(SP),DIMENSION(M,N),INTENT(IN)::U,V,ETA,Depth,MASK
  REAL(SP),INTENT(IN) :: DX,DY
  LOGICAL :: OUT_uvMA,OUT_VorMA,OUT_RadStr,OUT_etaMA
  ! moving average, temporal properties
  INTEGER, INTENT(IN) :: MA_Window
  REAL(SP),INTENT(IN) :: T_MA_STEP
  REAL(SP),INTENT(OUT) :: T_MA_Sum
  ! moving average of current speed
  REAL(SP),DIMENSION(M,N) :: U_temp1,V_temp1
  REAL(SP),DIMENSION(M,N,MA_Window),INTENT(INOUT) :: U_Window,V_Window
  REAL(SP),DIMENSION(M,N) :: U_Window_Sum,V_Window_Sum
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: uMA,vMA
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: uMAx, uMAy, vMAx, vMAy
  ! moving averaged of radiation shear stress
  REAL(SP),DIMENSION(M,N,MA_Window),INTENT(INOUT) :: Syx_temp
  REAL(SP),DIMENSION(M,N) :: Syx_temp_sum,Syx
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Syx_x
  ! moving average of surface elevation
  REAL(SP),DIMENSION(M,N,MA_Window),INTENT(INOUT) :: eta_MA_temp
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: eta_MA
  REAL(SP),DIMENSION(M,N) :: eta_MA_Sum
  ! moving average of vorticity
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Vor_MA
  ! eddy viscosity
  REAL(SP),DIMENSION(M,N,MA_Window),INTENT(INOUT) :: vis_smg_xMA_temp,vis_smg_yMA_temp,nu_smg_MA_temp
  REAL(SP),DIMENSION(M,N),INTENT(IN) :: vis_smg_x,vis_smg_y,nu_smg
  REAL(SP),DIMENSION(M,N),INTENT(OUT) :: vis_smg_xMA,vis_smg_yMA,nu_smg_MA
  REAL(SP),DIMENSION(M,N) :: vis_smg_xMA_SUM,vis_smg_yMA_SUM,nu_smg_MA_SUM

  INTEGER :: VTYPE,Nghost
  LOGICAL :: PERIODIC

  T_MA_Sum=T_MA_Sum+DT

  IF(T_MA_Sum.GE.T_MA_STEP)THEN
    ! shift the window and make space for new variables
    IF(OUT_uvMA.OR.OUT_RadStr.OR.OUT_VorMA)THEN
      U_Window = CSHIFT(U_Window,-1,3)
      V_Window = CSHIFT(V_Window,-1,3)
      ! replace oldest variable by the newest
      U_Window(:,:,1) = U
      V_Window(:,:,1) = V
      ! only sum up the values within the window
      U_Window_Sum = SUM(U_Window,3)
      V_Window_Sum = SUM(V_Window,3)
      ! do averaging
      uMA = U_Window_Sum/MA_Window
      vMA = V_Window_Sum/MA_Window
      CALL DERIVATIVE_X(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,vMA,vMAx)
      CALL DERIVATIVE_Y(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,vMA,vMAy)
      CALL DERIVATIVE_X(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,uMA,uMAx)
      CALL DERIVATIVE_Y(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DY,uMA,uMAy)
      VTYPE=1
      CALL PHI_COLL(M,N,Ibeg,Iend,Jbeg,Jend,Nghost,uMAx,VTYPE,PERIODIC)
      CALL PHI_COLL(M,N,Ibeg,Iend,Jbeg,Jend,Nghost,vMAy,VTYPE,PERIODIC)
      VTYPE=3
      CALL PHI_COLL(M,N,Ibeg,Iend,Jbeg,Jend,Nghost,uMAy,VTYPE,PERIODIC)
      Vtype=2
      CALL PHI_COLL(M,N,Ibeg,Iend,Jbeg,Jend,Nghost,vMAx,VTYPE,PERIODIC)

      ! subtracting average velocity from the instant velocity
      U_temp1 = U - uMA
      V_temp1 = V - vMA
    ENDIF
    IF(OUT_RadStr)THEN
      ! Calculate the radiation shear stress
      Syx_temp = CSHIFT(Syx_temp,-1,3)
      Syx_temp(:,:,1) = RHO_WATER*U_temp1*V_temp1*(Depth+Eta)
      Syx_temp_sum = SUM(Syx_temp,3)
      Syx = Syx_temp_sum/MA_Window
      CALL DERIVATIVE_X(M,N,Ibeg,Iend,Jbeg,Jend,MASK,DX,Syx,Syx_x)
    ENDIF
    IF(OUT_VorMA)THEN
      ! moving average of vorticity
      Vor_MA=vMAx-uMAy
    ENDIF
    IF(OUT_etaMA)THEN
      ! moving average of the surface elevation
      eta_MA_temp = CSHIFT(eta_MA_temp,-1,3)
      eta_MA_temp(:,:,1) = ETA
      eta_MA_Sum = SUM(eta_MA_temp,3)
      eta_MA = eta_MA_Sum/MA_Window
    ENDIF

    ! moving average of the horizontal viscosity
    vis_smg_xMA_temp = CSHIFT(vis_smg_xMA_temp,-1,3)
    vis_smg_yMA_temp = CSHIFT(vis_smg_yMA_temp,-1,3)
    nu_smg_MA_temp = CSHIFT(nu_smg_MA_temp,-1,3)

    vis_smg_xMA_temp(:,:,1) = vis_smg_x
    vis_smg_yMA_temp(:,:,1) = vis_smg_y
    nu_smg_MA_temp(:,:,1) = nu_smg

    vis_smg_xMA_SUM = SUM(vis_smg_xMA_temp,3)
    vis_smg_yMA_SUM = SUM(vis_smg_yMA_temp,3)
    nu_smg_MA_SUM = SUM(nu_smg_MA_temp,3)

    vis_smg_xMA = vis_smg_xMA_SUM/MA_Window
    vis_smg_yMA = vis_smg_yMA_SUM/MA_Window
    nu_smg_MA = nu_smg_MA_SUM/MA_Window

    ! go to next time step
    T_MA_Sum=T_MA_Sum-T_MA_STEP

  ENDIF
END SUBROUTINE CALCULATE_MovingAverage
