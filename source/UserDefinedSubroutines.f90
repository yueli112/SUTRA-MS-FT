!     SUBROUTINE        B  C  T  I  M  E       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE WHICH ALLOWS THE USER TO SPECIFY:     
! ***   (1) TIME-DEPENDENT SPECIFIED PRESSURES AND TIME-DEPENDENT       
! ***       CONCENTRATIONS OR TEMPERATURES OF INFLOWS AT THESE POINTS   
! ***   (2) TIME-DEPENDENT SPECIFIED CONCENTRATIONS OR TEMPERATURES     
! ***   (3) TIME-DEPENDENT FLUID SOURCES AND CONCENTRATIONS             
! ***       OR TEMPERATURES OF INFLOWS AT THESE POINTS                  
! ***   (4) TIME-DEPENDENT ENERGY OR SOLUTE MASS SOURCES                
!                                                                       
      SUBROUTINE BCTIME (IPBCT, IUBCT, IQSOPT, IQSOUT, GNUP, PM1,UM1,UM2, GNUPCJ, TDLEVEL, OUTSEEP)  ! CHENGJI 2013-09-02
      USE PARAMS 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE GRAVEC
      USE SutraStorage, ONLY : IPBC, PBC, IUBC, UBC, &
                               QIN, UIN, QUIN, IQSOP, IQSOU, &
                               X, Y, Z, &
                               SpecifiedPBC, &
                               MultiSpeciesBC,SLL1,SLL,SW1,SW,SII1,SII2 ,CFREEZ,SII ,SIOUT1,SIOUT2,SIOUT,SPACE,SPA ,EXTRA,VOL,RHOCHANGE,UCHANGE,QIN1,UCHANGE1,CTRNQ1,VOL1, &               !!!
                               IFREEZE1,  IFREEZE2,UMMIN
       USE SutraZoneModule  ,ONLY:NODEDATA !!为了nodedata%por
      USE M_TIDE
      USE M_SEEPAGE
      USE M_CONTROL
      USE SutraMSPrecision
      USE BCname
	  USE M_OUT
       USE M_SWCC
      USE surface
      IMPLICIT NONE

      integer (I4B) :: &
        IPBCT, IUBCT, IQSOPT, IQSOUT
		
!============================================= CHENGJI 2013-09-02	
      REAL(DP):: GNUPCJ, GNUP, PM1,UM1,UM2
      REAL(DP):: TDLEVEL
      REAL(DP):: SEEPX,SEEPZ
      INTEGER(I4B) :: OUTSEEP ! CHENGJI 2015-08-03
      DIMENSION GNUPCJ(NBCN) 
      DIMENSION PM1(NN)
      DIMENSION UM1 (NN, NSPE)
      DIMENSION UM2 (NN, NSPE)      
!================================================================
!================================================= !MT - 04112017	
      REAL(DP):: DRDTSW, SWT
      REAL(DP):: TIMEIND,TIMEIND2
      REAL(DP):: RHOSWB
      INTEGER(I4B):: NOMTH, NOYEAR

!============================================= yu  蒸发
      REAL (DP) :: &
           SHCAP,DENAIR,ETrate1,ETes, ETlv,ETra,ETrs, &
           RHOFWB,rainhk,DELTULRAIN,FLOWGREEM,  &
           ETa1,ETqa,ETQS,ETQg,ETes2, ETQs2,    &
           Rneng,Heng,LEeng,heatback,CHI,ETa2,QSED   !!,VNYU,VNFYU,SWRM1YU,SWRESYU,AAYU,AAPVNYU,AAPVNNYU
              
      REAL(DP),dimension(:),allocatable :: ETrate2,     &
           ETFLUXBACKSALT,  FLOWGREEMin

     REAL (DP) :: DELTASI
!================================================================	  
    
    
!     LOCAL VARIABLES
      INTEGER (I4B) :: &
        I, IP, IU, IUP, IQP, IQU, &
        K, &
        NSOPI, NSOUI
     REAL(DP):: THD          ,YUZHI, CC1,CC2                   !!
      
      IF(SEEP.EQ.1) THEN
        OPEN(8,FILE='SEEPAGE_FACE.DAT')
        SEEPX=XSTART
        SEEPZ=ZSTART
      ENDIF

         
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     PBC(IP)   = SPECIFIED PRESSURE VALUE AT IP(TH) SPECIFIED          
!                 PRESSURE NODE                                         
!     UBC(IP,K) = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY   
!                 INFLOW OCCURRING AT IP(TH) SPECIFIED PRESSURE NODE    
!                 FOR EACH SPECIES                                      
!     IPBC(IP)  = ACTUAL NODE NUMBER OF IP(TH) SPECIFIED PRESSURE NODE  
!                 {WHEN NODE NUMBER I=IPBC(IP) IS NEGATIVE (I<0),       
!                 VALUES MUST BE SPECIFIED FOR PBC AND UBC.}            
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     UBC(IUP,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE AT     
!                   IU(TH) SPECIFIED CONCENTRATION OR TEMPERATURE NODE  
!                   (WHERE IUP=IU+NPBC)                                 
!                   FOR EACH SPECIES                                    
!     IUBC(IUP,K) = ACTUAL NODE NUMBER OF IU(TH) SPECIFIED              
!                   CONCENTRATION OR TEMPERATURE NODE                   
!                   (WHERE IUP=IU+NPBC)                                 
!                   {WHEN NODE NUMBER I=IUBC(IU) IS NEGATIVE (I<0),     
!                   A VALUE MUST BE SPECIFIED FOR UBC.}                 
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOP(IQP) = NODE NUMBER OF IQP(TH) FLUID SOURCE NODE.            
!                  {WHEN NODE NUMBER I=IQSOP(IQP) IS NEGATIVE (I<0),    
!                  VALUES MUST BE SPECIFIED FOR QIN AND UIN.}           
!     QIN(-I)    = SPECIFIED FLUID SOURCE VALUE AT NODE (-I)            
!     UIN(-I,K)  = SPECIFIED CONCENTRATION OR TEMPERATURE VALUE OF ANY  
!                  INFLOW OCCURRING AT FLUID SOURCE NODE (-I)           
!                  FOR EACH SPECIES                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     IQSOU(IQU,K) = NODE NUMBER OF IQU(TH) ENERGY OR                   
!                    SOLUTE MASS SOURCE NODE                            
!                    {WHEN NODE NUMBER I=IQSOU(IQU) IS NEGATIVE (I<0),  
!                    A VALUE MUST BE SPECIFIED FOR QUIN.}               
!                    FOR EACH SPECIES                                   
!     QUIN(-I,K)  = SPECIFIED ENERGY OR SOLUTE MASS SOURCE VALUE        
!                   AT NODE (-I)                                        
!                   FOR EACH SPECIES                                    
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!.....ADDITIONAL USEFUL VARIABLES                                       
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!     "FUNITS" ARE UNIT NUMBERS FOR INPUT AND OUTPUT FILES              
!         AS ASSIGNED IN THE INPUT FILE, "SUTRA.FIL"                    
!                                                                       
!     X(I), Y(I), AND Z(I) ARE THE X-, Y-, AND Z-COORDINATES OF NODE I  
!     (FOR 2-D PROBLEMS, Z(I) IS THE THICKNESS AT NODE I)               
!                                                                       
!     GRAVX, GRAVY AND GRAVZ ARE THE X-, Y-, AND Z-COMPONENTS OF THE    
!     GRAVITY VECTOR                                                    
!     (FOR 2-D PROBLEMS, GRAVZ = 0)                                     
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
!                                                                       
!                                                                       
!.....NSOPI IS ACTUAL NUMBER OF FLUID SOURCE NODES                      
      NSOPI = NSOP - 1 
!.....NSOUI IS ACTUAL NUMBER OF ENERGY OR SOLUTE MASS SOURCE NODES      
      NSOUI = MNSOU - 1 
!                                                                       
!                                                                      
!      READ(5,*) TDLEVEL
!      IF(MOD(IT,354).EQ.0) THEN
!	   REWIND(5)
!      ENDIF
!
IF(IFEVAP.EQ.1)THEN
     allocate (ETFLUXBACKSALT(NSOPI)) 
      allocate (ETrate2(NSOPI))               
      allocate (FLOWGREEMin(NSOPI))                                                                      
      CALL ZERO (ETrate2, NSOUI, 0.0D0)
      CALL ZERO (FLOWGREEMin, NSOUI, 0.0D0) 
      CALL ZERO (ETFLUXBACKSALT, NSOUI, 0.0D0)  
ENDIF      
!                                                                       
!                                                                       
     IF (IPBCT) 50, 240, 240 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (1):  SET TIME-DEPENDENT SPECIFIED PRESSURES OR           
!     CONCENTRATIONS (TEMPERATURES) OF INFLOWS AT SPECIFIED             
!     PRESSURE NODES                                                    
!                                                                       
   50 CONTINUE
      DO 200 IP = 1, NPBC 
         I = SpecifiedPBC(IP)%node
         IF (I) 100, 200, 200 
  100    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH PBC( ) CHANGES.                         
!     SpecifiedPBC(IP)%P =  ((          ))                                         
!     DO 150 K=1,NSPE                                                   
! 150 SpecifiedPBC(IP)%U(K) =  ((          ))  

!================================================= CHENGJI 2013-09-02
TDLEVEL=TM+TASP*SIN(2*3.1415926*TSEC/TPSP)+TANE*SIN(2*3.1415926*TSEC/TPNP)
 SpecifiedPBC(IP)%P = 9.81*(1000+SC*714.3)*(TDLEVEL-Y(IABS(I))) 
!	  WRITE(*,*) TDLEVEL, SWT, SC
     IF (IABS(I).EQ.1) THEN        !!!!把表面节点参数输出到屏幕；此处1代表平台第一个点（与编号有关）
      write (*,*) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^'
        write (*,*) 'Tide Level =', TDLEVEL ,'Tide Temperature =', SpecifiedPBC(IP)%U(1)  
      write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
     END IF
     
    DO 150 K=1,NSPE
      IF(K.EQ.1) THEN
        IF(Y(IABS(I)).LE.TDLEVEL) THEN              
            SpecifiedPBC(IP)%U(K) = BNDSWT		!!!equal to tidal temperature
        ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
            SpecifiedPBC(IP)%U(K) = FWT		!!!equal to ambient temperature
        ENDIF
      ELSEIF(K.EQ.2) THEN
			IF(Y(IABS(I)).LE.TDLEVEL) THEN              
                SpecifiedPBC(IP)%U(K) = SC
            ELSEIF(Y(IABS(I)).GT.TDLEVEL) THEN
                SpecifiedPBC(IP)%U(K) = 0
            ENDIF										 
      ENDIF
  150 CONTINUE	
  
     ! IF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).GT.THD)THEN
	    !SpecifiedPBC(IP)%P = 0
	    !GNUPCJ(IP)=GNUP  
	    !  IF(SEEP.EQ.1) THEN
	    !    IF(Y(IABS(I))+PM1(IABS(I))/(9.81*(1000+SC*714.3)).GT.SEEPZ) THEN
	    !     SEEPZ=Y(IABS(I))+PM1(IABS(I))/(9.81*(1000+SC*714.3))
     !        SEEPX=X(IABS(I))
     !       ENDIF
     !     ENDIF
     ! ELSEIF(PM1(IABS(I)).GT.0.AND.Y(IABS(I)).LE.THD)THEN
	    !GNUPCJ(IP) = GNUP
     ! ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).GT.THD)THEN
    	!GNUPCJ(IP)= 0
     ! ELSEIF(PM1(IABS(I)).LT.0.AND.Y(IABS(I)).LE.THD)THEN
    	!GNUPCJ(IP) = GNUP
     ! ENDIF
!=====================================================================                 
  200 END DO
  
!      IF(SEEP.EQ.1) THEN
!        WRITE(8,'(I10,5E16.7)') IT,TDLEVEL,SEEPX,SEEPZ,BNDSWT,SC
!      ENDIF
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  240 IF (IUBCT) 250, 440, 440 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (2):  SET TIME-DEPENDENT SPECIFIED                        
!     CONCENTRATIONS (TEMPERATURES)                                     
!                                                                       
250 CONTINUE 
    
    
    
    
!    IF (IT .eq. 1) then
!     ALLOCATE (SURTEM1(4320)) 
!  		!write (*,*) 'ITMAX ='  , ITMAX ,'it=',IT
!       OPEN(2024327,FILE='surfacetem.txt',STATUS='UNKNOWN')   
!    DO 1029 I=1,4320
!        READ(2024327,*) SURTEM1(I)
!1029 CONTINUE
!     CLOSE(2024327)
!    EndIF  
!    
    
    
    
    
      DO 400 K = 1, NSPE 
         DO 350 IU = 1, NUBC (K) 
            IUP = NPBC + IU 
            I = MultiSpeciesBC(K)%SpecifiedU(IU)%node

            IF (I) 300, 400, 400 
300         CONTINUE 
            
       IF(K.EQ.1) THEN 
       !   IF(IT.LT.2880)THEN    !!1440的步长30是12h
       !      MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(IT)-3
       !   !   MultiSpeciesBC(K)%SpecifiedU(IU)%U=-5
       !   ELSEIF(IT.GT.2880.AND.IT.LE.5760) THEN  
       !     MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(IT) 
       !    !  MultiSpeciesBC(K)%SpecifiedU(IU)%U=5
       !   
       !   ENDIF     
  IF(BOUNDSPE.EQ.1) THEN
            IF(IT.LE.4320) THEN
         MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(IT)   
         ELSEif (IT.ge.4320)  then
        MultiSpeciesBC(K)%SpecifiedU(IU)%U=-10
        ENDIF 
  ENDIF  
  IF(BOUNDSPE.EQ.2) THEN     
       IF(IT.LT.2160) THEN
        MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(IT)   
       ELSEif (IT.ge.2160)  then
       MultiSpeciesBC(K)%SpecifiedU(IU)%U=-5.5
       ENDIF
    ENDIF
    
   IF(BOUNDSPE.EQ.3) THEN     
         IF(IT.LT.6480)THEN    !!1440的步长30是12h
             MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(IT)
          ELSEIF(IT.Ge.6480) THEN  
            MultiSpeciesBC(K)%SpecifiedU(IU)%U=SURTEM1(6480) 
       ENDIF
    ENDIF  
       
       
       
       

        ELSEIF(K.EQ.2) THEN                                  
             MultiSpeciesBC(K)%SpecifiedU(IU)%U=SC
        ENDIF  
            
            
!            IF (I) 300, 400, 400 
!  300       CONTINUE 
!!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY TIME STEP        
!!              IN WHICH UBC( ) CHANGES.  IN ADDITION, IF FLUID          
!!              PROPERTIES ARE SENSITIVE TO 'U' THEN A FLOW SOLUTION     
!!              MUST OCCUR AS WELL
!        IF(K.EQ.1) THEN 
!           IF(Z(IABS(I)).LE.TDLEVEL) THEN              
!            MultiSpeciesBC(K)%SpecifiedU(IU)%U = BNDSWT    !!!equal to tidal temperature
!           ELSEIF(Z(IABS(I)).GT.TDLEVEL) THEN
!            MultiSpeciesBC(K)%SpecifiedU(IU)%U = FWT    !!! equal to ambient temperature
!           ENDIF
!            IF (IU.EQ.1) THEN        !!!!把表面节点参数输出到屏幕；此处1代表平台第一个点（与编号有关）
!!           write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!		    write (*,*) '   Soil Temperature =', UM1(IABS(I),1)
!		    write (*,*) '   Air Temperature =', FWT
!            write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!            ENDIF  
!        ELSEIF(K.EQ.2) THEN                                  
!             MultiSpeciesBC(K)%SpecifiedU(IU)%U=SC
!        ENDIF
  350    END DO 
400 END DO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
! 440 IF (IQSOPT) 450, 840, 840 
440 IF (IQSOPT) 450, 640, 640 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!.....SECTION (3):  SET TIME-DEPENDENT FLUID SOURCES/SINKS,             
!      OR CONCENTRATIONS (TEMPERATURES) OF SOURCE FLUID                 
!                                                                       
  450 CONTINUE 
      DO 600 IQP = 1, NSOPI 
         I = IQSOP (IQP) 
         IF (I) 500, 600, 600 
  500    CONTINUE 
!     NOTE : A FLOW AND TRANSPORT SOLUTION MUST OCCUR FOR ANY           
!            TIME STEP IN WHICH QIN( ) CHANGES.                         
!     QIN(-I) =   ((           ))                                       
!     NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                    
!            TIME STEP IN WHICH UIN( ) CHANGES.                         
!     DO 550 K=1,NSPE                                                   
!       UIN(-I,K) =   ((           ))                                   
! 550 CONTINUE


     IF(IT.LE.1) THEN
        SIOUT1(-I)    = 0
        SIOUT2(-I)=0
        SIOUT(-I)    = 0
        SPA(-I)      = 0
        SPACE(-I)    = 0
        SW1(-I)=0
        SLL1(-I)=0
        SII1(-I)=0
        SII2(-I)=0
        CFREEZ(-I)=UM1(-I,2)
        IFREEZE1(-I) = 0
        IFREEZE2(-I) = 0     
        UMMIN(-I)    = 1E-5
         RHOCHANGE(-I)=(1000+DRWDU(2)*UM1(-I,2)) 
         UCHANGE(-I)=UM1(-I,2)
         VOL1(-I)=0
         EXTRA(-I)=0
     ENDIF
                          !     IF(IFREEZE1(-I).EQ.1) EXTRA(-I)=1-SW1(-I)
QIN(-I) =0     
 IF(IFSUCTION.EQ.1)THEN    
     !VV=2.26/20000
     !IF(X(-I).EQ.0.OR.ABS(-I).GT.20201) THEN
     !  VV=VV/2 
     !  !WRITE(*,*)'X'
     !ENDIF
     !IF(Y(-I).EQ.0.OR.Y(-I).EQ.1) THEN
     !  VV=VV/2
     !  !WRITE(*,*)'Y'
     !ENDIF
     IF (UMMIN(-I) .GT. UM1(-I,1) )  UMMIN(-I)=UM1(-I,1)   !!UMMIN取到该节点降温过程的最小值
!     UMMIN(-I) = MIN(UM1(-I,1), UMMIN(-I))
     
     IF( UM1(-I,1) .LT. 0.0 .AND. UM2(-I,1) .GE. 0.0 )  THEN     !!!用IFREEZE1和IFREEZE2来判断进入冻结还是融化的外部冰程序。
         IFREEZE2(-I) = IFREEZE1(-I)
         IFREEZE1(-I) = IFREEZE1(-I)+1
      EXTRA(-I)=1-SW1(-I)
     !ELSEIF(UM1(-I,1) .GE. -0.5 .AND. UM2(-I,1) .LT. -0.5)    THEN
     !    IFREEZE2(-I) = IFREEZE1(-I)
     !    IFREEZE1(-I) = IFREEZE1(-I)-1     
     !ELSEIF(UM1(-I,1) .GE. 0.0 .AND. UM2(-I,1) .LT. 0.0 .AND. UMMIN(-I) .GT. -0.5 )    THEN  !!!加入UMMIN(-I)，用于冻融没结束的融化过程
     !    IFREEZE2(-I) = IFREEZE1(-I)
     !    IFREEZE1(-I) = IFREEZE1(-I)-1             
     ELSEIF(UM1(-I,1) .GT.0   .AND. UM2(-I,1) .GT. 0)    THEN  !!!加入UMMIN(-I)，用于冻融没结束的融化过程
         IFREEZE2(-I) = IFREEZE1(-I)
         IFREEZE1(-I) = IFREEZE1(-I)-1  
      !ELSEIF(UM1(-I,1) .LT. UM2(-I,1)  .AND. UMMIN(-I) .GT. -0.5.AND.UM2(-I,1).LT.0.AND.UM1(-I,1).LT.0)    THEN  !!!加入UMMIN(-I)，用于冻融没结束的融化过程
      !   IFREEZE2(-I) = IFREEZE1(-I)
      !   IFREEZE1(-I) = IFREEZE1(-I)+1      
     ENDIF
 QIN(-I)=0      

     IF(IT.GE.2) THEN
         !     IF(SIOUT1(-I).GT.1-SW1(-I))SIOUT1(-I)=1-SW1(-I)
    ! WRITE(*,*)'RGSJSYHJYTY'
        SPACE(-I)=1-SW1(-I) 
        SPA(-I)=SPACE(-I)-SIOUT1(-I)!!剩余外部空间    
       ! IF (IFREEZE1(-I) .GT. IFREEZE2(-I).AND.UM1(-I,1).GE.-1) THEN   
       !WRITE(*,*)'DEV1',DEV,DM1
         IF (IFREEZE1(-I) .GT. IFREEZE2(-I).AND.UM1(-I,1).GE.-1.and.UM1(-I,1).lt.(UM2(-I,1)-DEV/3600*delt)) THEN  !!变幅太慢不吸
           !  WRITE(*,*)'DEV2',DEV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      进入冻结的外部冰程序     
 !   IF(SII(-I).GT.0.AND.SPA(-I).GT.0.01  )  THEN
! IF(SII(-I).GT.0.AND.SPA(-I).GE.0.0  )  THEN
!IF ((PM1(-I) .GE. -100000).AND.(SW(-I).GT.0.15)) THEN

IF(SII1(-I).GT.0.0.AND.SPA(-I).GT.0.01  )  THEN   !!如果有冰(负温)且有空间（不饱和）才考虑吸水

!IF(SII1(-I).GT.0.AND.SPA(-I).GT.0.01 .AND. (SIOUT1(-I)+SII1(-I))/SLL1(-I) .LT. (1-SWRES1)/SWRES1  )  THEN   !!如果有冰(负温)且有空间（不饱和）才考虑吸水
          YUZHI=INT(SPA(-I) * 100.0) / 100.0 /DM1    ! INT(SPA(-I) * 1000.0) / 1000.0 ，减小阈值差异!!!
!!!!!!!!!!!!!  YUZHI=SPA(-I) /10，导致该层节点刚结冰时流速过大，应想办法限制第一次拿出去的量 
  QIN(-I) =0
IF(SIOUT1(-I).LE.EXTRA(-I))THEN
!IF(SLL1(-I).GT.0.2)THEN
!IF(SIOUT1(-I).LE.1)THEN
IF (SII1(-I) .LE. YUZHI) THEN
        SIOUT(-I) = SIOUT1(-I) + INT(SII1(-I) * 100.0) / 100.0                      
        QIN(-I) = - ( INT(SII1(-I) * 100.0) / 100.0  ) *VOL1(-I)*NodeData(1)%por*1000/DELT      !!  hardcoding, density = 1000 kg/m3, tstep = 30 s  
       
        !IF (SIOUT(-I) .GT. SPACE(-I) ) THEN
        !     SIOUT(-I) = SPACE(-I)
        !     QIN(-I) = -( SPA(-I) ) *VOL1(-I)*NodeData(1)%por*1000/DELT 
        !ENDIF    
    ELSEIF (SII1(-I) .GT. YUZHI) THEN
        SIOUT(-I) = SIOUT1(-I) + YUZHI                          
        QIN(-I) = -( YUZHI )*VOL1(-I)*NodeData(1)%por*1000/DELT        !!  hardcoding, density = 1000 kg/m3, tstep = 30 s
        !IF (SIOUT(-I) .GT. SPACE(-I) ) THEN
        !     SIOUT(-I) = SPACE(-I)
        !     QIN(-I) = -( SPA(-I) ) *VOL1(-I)*NodeData(1)%por*1000/DELT 
        !ENDIF 
    ENDIF 
 ENDIF  !IF(SIOUT1(-I).LE.0.45)THEN 
     ELSE
            QIN(-I) =0             
     ENDIF  
 IF(SII1(-I)+SLL1(-I)+SIOUT(-I).GT.0.99) THEN
     QIN(-I) =0 
     SIOUT(-I) = SIOUT1(-I)
 ENDIF 
 ! IF( QIN(-I).GT.-1E-4*VOL1(-I)*NodeData(1)%por*1000/DELT) THEN
 !    QIN(-I) =0 
 !    SIOUT(-I) = SIOUT1(-I)
 !ENDIF
 
! SII1(-I)+SLL1(-I)+SIOUT(-I).GT.0.97.AND.SIOUT(-I).GT.0.AND. QIN(-I).GT.-1E-4*VOL1(-I)*NodeData(1)%por*1000/DELT 
 !IF (SII1(-I)+SLL1(-I)+SIOUT1(-I).GT.0.97.AND.SIOUT(-I).GT.0.AND. QIN(-I).GT.-1E-4*VOL1(-I)*NodeData(1)%por*1000/DELT ) THEN
 ! QIN(-I) = SIOUT1(-I)/10*VOL1(-I)*NodeData(1)%por*1000/DELT
 !              UIN(-I,1)=UM1(-I,1)
 !            UIN(-I,2)=0
 !            SIOUT(-I)=SIOUT1(-I)*9/10
 ! IF ( SIOUT1(-I) .LT. 0.000001 ) SIOUT1(-I)=0
 !ENDIF
 !
 
 !    !   ELSEIF (IT.GT.4800.AND.Y(-I).GT.4.5) THEN
 ! ELSEIF (SIOUT(-I).GT.0.AND. QIN1(-I).GT.-1E-4*VOL1(-I)*NodeData(1)%por*1000/DELT) THEN          
 !     IF ( SIOUT1(-I) .GT. 0.00 ) THEN  !!外部储冰返还
 !  !                            
 !            QIN(-I) = SIOUT1(-I)/10*VOL1(-I)*NodeData(1)%por*1000/DELT
 !              UIN(-I,1)=UM1(-I,1)
 !            UIN(-I,2)=0
 !            SIOUT(-I)=SIOUT1(-I)*9/10
 !         
 !      IF ( SIOUT(-I) .LT. 0.000001 ) SIOUT(-I)=0  
 !      
 !      
 !       ELSE
 !           QIN(-I) =0             
 !      ENDIF
 !   IF(SII1(-I)+SLL1(-I)+SIOUT1(-I).GT.0.99) THEN
 !    QIN(-I) =0 
 !    SIOUT(-I) = SIOUT1(-I)
 !ENDIF 
 
 
!     ELSEIF (IFREEZE1(-I) .LT. IFREEZE2(-I) ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     进入融化的外部冰程序     
!QIN(-I) =0   
!!IF (PM1(-I) .LE. 10000) THEN
!!IF(IT.GT.10000)THEN
!IF(IT.GT.1)THEN
!     !  IF ( SIOUT1(-I) .GT. 0.05 ) THEN  !!外部储冰返还
!     IF ( SIOUT1(-I) .GT. 0.00 ) THEN  !!外部储冰返还
!   !          SII(-I)=SII1(-I)+SIOUT1(-I)/4      !???                        
!             QIN(-I) = SIOUT1(-I)/10*VOL1(-I)*NodeData(1)%por*1000/DELT
!               UIN(-I,1)=UM1(-I,1)
!             UIN(-I,2)=0
!             SIOUT(-I)=SIOUT1(-I)*9/10
!           !  SIOUTLEFT(-I)=SIOUTLEFT(-I)+SIOUT1(-I)*1/8
!       IF ( SIOUT(-I) .LT. 0.000001 ) SIOUT(-I)=0  
!       !ELSEIF (SIOUT1(-I) .GT. 0.0 .AND. SIOUT1(-I) .LE. 0.05 ) THEN
!       !       SII(-I)=SII1(-I)+SIOUT1(-I)               
!       !   !    QIN(-I) = SIOUT1(-I)*VOL1(-I)*NodeData(1)%por*RHOCHANGE(-I)/DELT  
!       !    QIN(-I) = SIOUT1(-I)*VOL1(-I)*NodeData(1)%por*1000/DELT  
!       !       UIN(-I,1)=UM1(-I,1)
!       !       UIN(-I,2)=UCHANGE(-I)
!       !       SIOUT(-I)=0
!    
!       !QIN(-I) = SIOUT1(-I)*VOL1(-I)*NodeData(1)%por*1000/DELT
!       !IF(SIOUT1(-I)+SW(-I).GT.1) QIN(-I)=(1-SW(-I))*VOL1(-I)*NodeData(1)%por*1000/DELT
!       !UIN(-I,1)=UM1(-I,1)
!       !UIN(-I,2)=UCHANGE(-I)
!       !SIOUT(-I)=0
!    
!       
!       
!       ELSE
!            QIN(-I) =0             
!       ENDIF
!  ENDIF   !IF(IT.GT10000 
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
  
 !      !   ELSEIF (IT.GT.4800.AND.Y(-I).GT.4.5) THEN
 ! IF (SIOUT1(-I).GT.0.01.AND. QIN1(-I).GT.-1E-6*VOL1(-I)*NodeData(1)%por*1000/DELT) THEN          
 !  IFREEZE2(-I) = IFREEZE1(-I)
 !  !                            
 !            QIN(-I) = SIOUT1(-I)/10*VOL1(-I)*NodeData(1)%por*1000/DELT
 !              UIN(-I,1)=UM1(-I,1)
 !            UIN(-I,2)=0
 !            SIOUT(-I)=SIOUT1(-I)*9/10
 !         
 !      IF ( SIOUT(-I) .LT. 0.00001 ) SIOUT(-I)=0  
 !      
 !      
 !  
 !   IF(SII1(-I)+SLL1(-I)+SIOUT1(-I).GT.0.99) THEN
 !    QIN(-I) =0 
 !    SIOUT(-I) = SIOUT1(-I)
 !ENDIF 
 ! 
 !ENDIF 
  IF(DM2.NE.9999)THEN
 IF (UM1(-I,1).LE.-1.AND.SIOUT1(-I).GT.0.0) THEN
  !    IF (IT.GT.2880.AND.SIOUT1(-I).GT.0.0) THEN
! IF(SLL(-I).EQ.0.1.AND.SIOUT1(-I).GT.0.0)THEN
     IFREEZE2(-I) = IFREEZE1(-I)
   !                            
             QIN(-I) = SIOUT1(-I)/DM2*VOL1(-I)*NodeData(1)%por*(1000+DRWDU(2)*UM1(-I,2))/DELT   
               UIN(-I,1)=UM1(-I,1)
             UIN(-I,2)=0
             SIOUT(-I)=SIOUT1(-I)*(DM2-1)/DM2
          
       IF ( SIOUT(-I) .LT. 0.000001 ) SIOUT(-I)=0  
       
       
   
    IF(SII1(-I)+SLL1(-I).GT.0.99) THEN
     QIN(-I) =0 
     SIOUT(-I) = SIOUT1(-I)
    ENDIF   
  
  ENDIF 
  !


! IF (UM1(-I,1).GT.UM2(-I,1).AND.UMMIN(-I).GT.-1.AND.UM1(-I,1).LE.-0.5 ) THEN
     IF (SUM(SII1(:)).LE.0 .and.IFREEZE1(-I) .LE. IFREEZE2(-I)) THEN                          !!!!!!换成IFREEZE2这个判别融化，别靠IT
!   IF (SII1(-I).LE.0.0 ) THEN   
!  IF (UM1(-I,1).GT.0 ) THEN   
!!  !IF (SII1(-I)+SIOUT1(-I) .LT. SII2(-I)+SIOUT2(-I) ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     进入融化的外部冰程序     
QIN(-I) =0   
     IF ( SIOUT1(-I) .GT. 0.00 ) THEN  !!外部储冰返还
   !          SII(-I)=SII1(-I)+SIOUT1(-I)/4      !???                        
             QIN(-I) = SIOUT1(-I)/DM2*VOL1(-I)*NodeData(1)%por*1000/DELT
               UIN(-I,1)=UM1(-I,1)
             UIN(-I,2)=0
             SIOUT(-I)=SIOUT1(-I)*(DM2-1)/DM2
           !  SIOUTLEFT(-I)=SIOUTLEFT(-I)+SIOUT1(-I)*1/8
       IF ( SIOUT(-I) .LT. 0.000001 ) SIOUT(-I)=0     
     ENDIF
!     
!     !IF(SII1(-I)+SLL1(-I).GT.0.99.AND.SIOUT1(-I) .GT. 0.00) THEN
!     !QIN(-I) =SIOUT1(-I)/200*VOL1(-I)*NodeData(1)%por*1000/DELT 
!     !           UIN(-I,1)=UM1(-I,1)
!     !        UIN(-I,2)=UCHANGE(-I)
!     !        SIOUT(-I)=SIOUT1(-I)*199/200
!     !        IF ( SIOUT(-I) .LT. 0.00001 ) SIOUT(-I)=0  
!     !ENDIF 
!   !    IF(SII1(-I)+SLL1(-I).GT.0.99.AND.SIOUT1(-I) .GT. 0.00) THEN
!   !        IF(SIOUT1(-I) .GT. 0.005)THEN
!   !  QIN(-I) =0.0005*VOL1(-I)*NodeData(1)%por*1000/DELT 
!   !             UIN(-I,1)=UM1(-I,1)
!   !          UIN(-I,2)=UCHANGE(-I)
!   !          SIOUT(-I)=SIOUT1(-I)-0.005
!   !        ENDIF
!   !          IF ( SIOUT1(-I) .LT. 0.005 )THEN 
!   !           QIN(-I) =SIOUT1(-I)*VOL1(-I)*NodeData(1)%por*1000/DELT    
!   !               UIN(-I,1)=UM1(-I,1)
!   !          UIN(-I,2)=UCHANGE(-I)
!   !          SIOUT(-I)=0  
!   !    ENDIF  
!   !ENDIF    
!       
  ENDIF
 ENDIF !IF(DM2.NE.9999)THEN   
     
  
     
     
     
  
 !ENDIF
     ENDIF  !IF(IT.GE.1)
 !IF(IT.EQ.400.AND.I==-302)THEN
 !QIN(302)=  0.2*VOL1(-I)*NodeData(1)%por*1000/DELT 
 !UIN(302,2)=UM1(302,2)
 !UIN(302,1)=UM1(302,1)
 !ENDIF
 
ENDIF!!!!IFSUCTION==1 
  ! IF(I==-302)WRITE(30701,'(2I9,20E15.6E3)')IT, 2,UM1(-I,1),UM1(-I,2),SLL(-I),SII(-I),SW(-I),QIN(-I) ,CFREEZ(-I)   
     
     
!IF(I==-302)WRITE(*,*)'12', QIN(302)       
     
!IF(IT.EQ.1.AND.ABS(I).EQ.1)THEN
!OPEN(240528,FILE='SPACE_AND_SIOUT.DAT')
!OPEN(240604,FILE='SUM.DAT')
!     WRITE(240528,575)   
!575 FORMAT(8X,'  IT          I           SPACE            SIOUT           SPA ')
!ENDIF 
! IF(ABS(I).EQ.20301) WRITE(240604,5751)IT, SUM(SIOUT),SUM(QIN),SUM(SW),SUM(SLL),SUM(SII)
!!WRITE(*,*)'SW1=',SW1(-I),SW(-I),'SLL1=',SLL1(-I),'SII1=',SII1(-I),'SPACE=',SPACE(-I),'SIOUT=',SIOUT(-I),'SPA=',SPA(-I)
!!WRITE(*,*)'QIN=',QIN(-I),'VOL=',VV
!IF(MOD(IT,10).EQ.0) THEN
!write (240528,574) IT,-I,SPACE(-I),SIOUT(-I),SPA(-I)
!ENDIF
!574 FORMAT(2I,12E15.5E3) 
!5751 FORMAT(1I,12E15.6) 
!  IF(IT.EQ.8640) THEN
!CLOSE(240528)
!  ENDIF
  
!  QIN(-I) = -0.01*X(-I)*Y(-I)*NodeData(1)%por*1000/300



  
IF(IFEVAP.EQ.1)THEN

!取上一步饱和度，温度，盐度

 SWYU   = SLL1(IABS(I))                   !!取上一步饱和度 SLL1(I)，认为冰没有蒸发，只有少量液态水蒸发
  TSSOIL(IQP)=UM1(IABS(I),1)                           !!温度和盐度
  SSSOIL(IQP)=UM1(IABS(I),2)

!!求蒸发率ETrate
por1=NodeData(1)%por  !!取了第一区域的孔隙度
  rainhk = 1.22D-6  !!saturated hydraulic conductivity (m/s)    !!目前没用到
 SHCAP = 1003    ! specific heat capacity of the air
     DENAIR = 1.205  ! air density

 ETlv = 2.5*1000000-2369.2*TSAIR(IT)  ! volumetric latent heat of vaporization
    ETra = 94.909 *  (Uwind(IT))**(-0.9036) !!! aerodynamic resistance  (s/m)
    ETrs = -8.05+41.4*(por1-por1*SWYU(IQP))               
    IF (ETrs .LE. 0) THEN
       ETrs=0
	  ENDIF

!    ETra = 230 /  Uwind(IT) !!! aerodynamic resistance  (s/m)  FAO
!    ETra = 132.5 / (1+0.54 * Uwind(IT)) !!! aerodynamic resistance  (s/m)  THom    
    ETes = 0.6108 * exp(17.27*(TSAIR(IT))/(TSAIR(IT)+237.3))   !!saturation vapor pressure (kg/m/s/s)
!    ETes2 = 0.6108 * exp(17.27*(TSSOIL(IQP))/(TSSOIL(IQP)+237.3))   !!saturation vapor pressure (kg/m/s/s)    
     ETQs=0.622*ETes/(ETP0(IT)/10-0.376*ETEs)  !!Saturated specific humidity
!     ETQs2=0.622*ETes2/(ETP0(IT)/10-0.378*ETEs2)   
!! #############################################################################################################
!    ETa1=EXP(0.018*PM1(IABS(I))/(1002-0.2205*TSSOIL(IQP))/8.31/(TSSOIL(IQP)+273.15)) !!relative humidity in soil
!    CHI=0.9D0+2.6D-3*(SSSOIL(IQP)/(1-SSSOIL(IQP)))*1.D3+2.6D-3*((SSSOIL(IQP)/(1-SSSOIL(IQP)))*1.D3)**5.3D-2
!    ETa2=EXP(-0.018*2*CHI*SSSOIL(IQP)/0.0585)
!    ETQg=ETa2*ETa1*ETQs2 
!! #############################################################################################################
       IF ((1.8*por1*SWYU(IQP)/(por1*SWYU(IQP)+0.30)).LE.1) THEN
       ETa1=(1.8*por1*SWYU(IQP)/(por1*SWYU(IQP)+0.30))
	  ELSE
      ETa1=1
	  ENDIF

        ETQg=ETa1*ETQs         
         ETQa=ETQs*Hraair(IT)/100 
           ETrate1=DENAIR*(ETQg-ETQa)/(ETra+ETrs)         !!!得到初步计算的蒸发率 ETrate1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!influx  还没有考虑降雨，所以目前FLOWGREEMin(IQP)=0 

      RHOFWB = 1002-0.2205*TSAIR(IT)  !!freshwater density(kg/m^3)
!       DELTULRAIN= Y(IABS(I))- Y(IABS(I)+ COLYU)
!	    FLOWGREEM= rainhk *(0.0 + DELTULRAIN*9.81*RHOFWB  -PM1(IABS(I)+ COLYU))/(DELTULRAIN*9.81*RHOFWB)   !!m/s
!        FLOWGREEMin(IQP)=FLOWGREEM
!      	IF (FLOWGREEM.GE.RAINYU(IT) ) THEN          !!hardcoding
!        FLOWGREEMin(IQP)=RAINYU(IT)
        FLOWGREEMin(IQP)=0        
!        ENDIF      

!         UIN(-I,1)=TSAIR(IT)
!         UIN(-I,2)=0.0001  
!	 QIN(-I) =   ((   Sarea(IQP) * RHOFWB *  FLOWGREEMin(IQP)  ))    !!!!!!!进流!!!!!!! (kg/s)
!     UIN(-I,1)=TSAIR(IT)
!     UIN(-I,2)=0.0     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!降雨↑↑↑↑


	 !!!! 判断是否发生蒸发    压力小于零   没有降雨 负出流
IF (IT .LT. 2) then         !!第一步默认没有
        QIN(IABS(I))=0  
        ETrate2(IQP) = 0
 ELSE        !!后面的步长时
    
       IF (ETrate1  .GE.0  ) THEN     !!蒸发率算的要是负值，就取0
            ETrate2(IQP)=ETrate1/1000  !     (m/s)       
       else
           ETrate2(IQP)=0                     !!ETrate2是最终蒸发率(m/s)
       ENDIF     
	   
		       IF (IQP.EQ.1) THEN    !!! 把表面节点蒸发率输出到屏幕  IF (IQP.EQ.561) THEN  
      write (*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      write (*,*) ' SATURATION =   ', SWYU(IQP)
	  write (*,*) ' *****    ET   =  ', ETrate2(IQP)*1000*3600*24,'MM/D'
               ENDIF 

       !!!算出蒸发率后赋值给出流节点  TDH和TDT是实测潮汐水位和水温
     IF ( TDH(IT)-Z(IABS(I)).LT. 0 .AND. PM1(IABS(I)) .GT. -1E10) THEN     !!PM1(IABS(I)) .GT. -1E10只要水不是太少，没淹潮汐就能蒸发

        QIN(-I) =  Sarea(IQP) * (1002-0.2205*TSSOIL(IQP))*  (-ETrate2(IQP))       !!!!!!!出流!!!!!!!  Sarea(IQP)也是输入条件（个数为NPBC?什么意思）
        UIN(-I,1)= 0.0    !!!
!        UIN(-I,2)= 0.0
    ELSE         
         QIN(-I) =  0.0       !!!!!!!出流!!!!!!! 
         UIN(-I,1)= 0.0   !!!
!         UIN(-I,2)= 0.0  
    endIF   
	
        IF (IQP.EQ.1) THEN        !!!!把表面节点参数输出到屏幕；此处1代表平台第一个点（与编号有关） IF (IQP.EQ.561)
      write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!        write (*,*) '   creekwaterlevel  ='  ,         TDH(IT)
!      write (*,*) DELTULRAIN,' I CAN INPUT',FLOWGREEM*1000*3600*24
!      write (*,*) '   RAINFALL  ='  ,       RAINYU(IT)
      write (*,*) ' RAINFALL-INFILTRATION OR ET =',  QIN(-I)/ Sarea(IQP)/1000 *1000*3600*24 ,'MM/D'
      write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
	      	 ENDIF  
			 
 ENDIF  !!结束后面步长的分支   
               
            ETFLUXBACKSALT(IQP)=QIN(-I)  !!!!   TO BE BACK  蒸发掉的流量？？  
            
END IF   !!!IF(IFEVAP.EQ.1)              
                
  550 CONTINUE

!  IF(I==-302) WRITE(*,*)'13', QIN(302)   !                                                         
  600 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!WRITE(*,*)'14', QIN(302)   !                                                                       
!                                                                       
  640 IF (IQSOUT) 650, 840, 840 
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!!.....SECTION (4):  SET TIME-DEPENDENT SOURCES/SINKS                    
!!     OF SOLUTE MASS OR ENERGY                                          
!!                                                                       
650 CONTINUE 
    
      DO 800 K = 1, NSPE 
         DO 750 IQU = 1, NSOU (K) 
            I = IQSOU (IQU, K) 
            IF (I) 700, 800, 800 
  700       CONTINUE 
!!       NOTE : A TRANSPORT SOLUTION MUST OCCUR FOR ANY                  
!!              TIME STEP IN WHICH QUIN( ) CHANGES.                      
!!       QUIN(-I,K) =   ((           ))

!VV=2.26/20000
!     IF(X(-I).EQ.0.OR.ABS(-I).GT.20201) THEN
!       VV=VV/2 
!       !WRITE(*,*)'X'
!     ENDIF
!     IF(Y(-I).EQ.0.OR.Y(-I).EQ.1) THEN
!       VV=VV/2
!       !WRITE(*,*)'Y'
!     ENDIF


IF(K.EQ.1)THEN
    QUIN(-I,K)=0
  
    IF(BOUNDSPE.EQ.1)  QUIN(-I,K)=0.0001*(-10-UM1(ABS(I),K))
   IF(BOUNDSPE.EQ.2)  QUIN(-I,K)=0 
   if (BOUNDSPE.EQ.3.and.it.gt.3720)  QUIN(-I,K)=0.001*(15-UM1(ABS(I),K))
    !      ! QUIN(-I,K)=0.003*(11-UM1(ABS(I),K))
    !      
    !     IF(ABS(I).LE.101.AND.ABS(I).GE.1) THEN
    !      QUIN(-I,K)= QUIN(-I,K)/2
    !     !IF(ABS(I).EQ.1) write (*,*) 'ZUO'
    !     ENDIF
    !     
    !      IF(ABS(I).GE.20201) THEN
    !      QUIN(-I,K)= QUIN(-I,K)/2
    !     !  IF(ABS(I).EQ.20301) write (*,*) 'YOU'
    !      ENDIF
! QUIN(-I,K)= QUIN(-I,K)-QIN(-I)*(334000+(4182 - 2108)*UM1(-I,1) ) 
ELSEIF (K.EQ.2)THEN           !!冻结冰的盐量返回至水中
   
  !           QUIN(-I,K)= 0.0  
  !           DELTASI = 0.0
  !   ! IF(I==-199)THEN
  !      !  IF (IFREEZE1(-I) .GT. IFREEZE2(-I) .AND. SII(-I).GT.0) THEN  
        IF(IT==1) THEN
            SII1(-I)=0
            SIOUT1(-I)=0
            SII2(-I)=0
            SIOUT2(-I)=0
            RHOCHANGE(-I)=(1000+DRWDU(2)*UM1(-I,2))  
            UCHANGE(-I)=UM1(-I,2) 
        ENDIF
        
  !          DELTASI= (SII(-I)-SII1(-I)+SIOUT(-I)-SIOUT1(-I))
  !          IF (IFREEZE1(-I) .LT. IFREEZE2(-I).AND.DELTASI .GT. 0) DELTASI = 0.0   !!! 冰       
  !         QUIN(-I,K)= DELTASI*VOL1(-I)*NodeData(1)%por*UCHANGE(-I)*RHOCHANGE(-I)/DELT
  !      ! QUIN(-I,K)= DELTASI*VOL1(-I)*NodeData(1)%por*UCHANGE(-I)*1000 /DELT
  !      !ENDIF     
  !     ! WRITE(*,*) SII(-I),SII1(-I),DELTASI,QUIN(-I,K)
  !  ! ENDIF
  !! IF(IT==200.AND. Y(-I) .EQ. 1) QUIN(-I,2)=1.3937E-5
  ! ! 
  !! IF(IT==233) QUIN(199,2)=-1.3937E-5 
QUIN(-I,K)=0
IF(IFSUCTION.EQ.1)THEN 
!IF(QIN(-I).LT.0)  QUIN(-I,K)=-QIN(-I)*UCHANGE(-I)- QIN1(-I)*  (UM1(-I,2)-UCHANGE1(-I)) - CTRNQ1(-I)*  (UM1(-I,2)-UCHANGE1(-I)) 
IF(QIN(-I).LT.0)  QUIN(-I,K)=-QIN(-I)*UCHANGE(-I)- QIN1(-I)*  (UM1(-I,2)-UCHANGE1(-I)) 
!IF(QIN(-I).GT.0)  QUIN(-I,K)=-QIN(-I)*UCHANGE(-I)
!IF(UM1(-I,1).GT.UM2(-I,1).AND.UMMIN(-I).GT.-1.AND.UM1(-I,1).LE.-0.5.AND.QIN(-I).GT.0)  QUIN(-I,K)=-QIN(-I)*UCHANGE(-I)
!IF(QIN(-I).GT.0)  QUIN(-I,K)=-QIN(-I)*(QIN(-I)/(VOL1(-I)*NodeData(1)%por*1000/DELT))/((QIN(-I)/(VOL1(-I)*NodeData(1)%por*1000/DELT))+SLL(-I))*UM1(-I,2)
  ! 
ENDIF!!IF(IFSUCTION.EQ.1)THEN 
! IF(I==-302)WRITE(30701,'(2I9,20E15.6E3)')IT, 2,UM1(-I,1),UM1(-I,2),SLL(-I),SII(-I),SIOUT(-I),QIN(-I)
    ! 
!IF(IT.EQ.1.AND.ABS(I).EQ.1)THEN
!OPEN(240528,FILE='SALT_TEST.DAT')
!ENDIF 
!IF(MOD(IT,1).EQ.0 .AND. ABS(I).EQ.10099) THEN
!write (240528,574) IT,-I,RHOCHANGE(-I),1000+710*UM1(-I,2),UM1(-I,1),UM1(-I,2)
!ENDIF
!574 FORMAT(2I,12E15.6)     
    ENDIF

  




!IF (K .EQ. 2) then
!          
!       IF (QIN(-I) .LT.0) THEN    !!如果有水蒸发走，就把带走的盐质量流回到表面
!           
!       QUIN(-I,k) =   (( QIN(-I)  *UM1(-I,K)        ))     !!!! 取通量负号  (kg/s)  重新获得水蒸发带走的盐
!	 ELSE
!	     QUIN(-I,k) = 0
!     ENDIF 
!  
!                
!ENDIF

!IF(IFEVAP.EQ.1)THEN
!
!      DO  IQU = 1, NSOU (K)
!               I = IQSOU (IQU, K)
!          if (IT .LT. 2) then
!        QUIN(-I,1)=0
!        QUIN(-I,2)=0
!          endif  
!     END DO                !!第一步取0
!
!       if (IT .GT. 1) then
!IF (K .EQ. 1) then 
!
!           DO  IQU = 1, NSOU (K)
!               I = IQSOU (IQU, K)
!            Heng= SHCAP*DENAIR/ETra*(TSSOIL(IQU)-TSAIR(IT))*( Sarea(IQU)  )        !!sensible heat(W/m/m)
!            Rneng= NETRn(IT) *( Sarea(IQU)  ) !净辐射
!            LEeng= ETrate2(IQU)*ETlv*1000 *( Sarea(IQU)  ) !潜热
!            heatback= -ETFLUXBACKSALT(IQU)*4182*TSSOIL(IQU)    !!离开的水的热量？？
!            QSED=0      !水传热????????不知道是什么
!	 	 IF(TDH(IT)-Z(IABS(I)).GT. 0) THEN    !!!! 如果潮汐淹没 
!              Heng  = 0
!              Rneng = 0
!              LEeng = 0
!              QSED = 1.5*(TDT(IT)-TSSOIL(IQU))/0.05 *( Sarea(IQU)  )        !!!! 网格垂向尺寸  潮汐传来的热？？0.05和1.5？？
!         endif
!           QUIN(-I,k)  =  Rneng-Heng-LEeng+heatback+QSED
!!         WRITE(20777,'(I10,5E16.6)')  IT,Rneng, -Heng,-LEeng, heatback, QSED
!            if (mod(IT,4).EQ.0.AND.K.EQ.1) then           !!! 输出显热和水传热的间隔
!                WRITE(20777,'(2I10,2E16.6)')  IT, -I,-Heng, QSED
!            endif				 
!        IF (IQU.EQ.1) THEN        !!!!把表面节点参数输出到屏幕；此处1代表平台第一个点（与编号有关）
!!      write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!        write (*,*) '   Heat_in ='  , QUIN(-I,1)
!		write (*,*) '   Soil Temperature ='  , TSSOIL(IQU)
!		write (*,*) '   Air Temperature ='  , TSAIR(IT)
!      write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
!        ENDIF  
!        enddo      
!           
! 
!ELSEIF (K .EQ. 2) then
!           DO  IQU = 1, NSOU (K)
!               I = IQSOU (IQU, K)
!       IF (ETFLUXBACKSALT(IQU) .LT.0) THEN    !!如果有水蒸发走，就把带走的盐质量流回到表面
!           
!       QUIN(-I,k) =   (( -ETFLUXBACKSALT(IQU)  *SSSOIL(IQU)         ))     !!!! 取通量负号  (kg/s)  重新获得水蒸发带走的盐
!	 ELSE
!	     QUIN(-I,k) = 0
!     ENDIF 
!  
!           enddo         
!
!           
!ENDIF
!
!       endif  
!
!ENDIF     !!!!IF(IFEVAP.EQ.1)       
   
  750 END DO       
  800 END DO 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
  840 CONTINUE 
!                                                                       
      RETURN 
      END SUBROUTINE BCTIME 

!                                                                       
!     SUBROUTINE        U  N  S  A  T              SUTRA VERSION 2D3D.1 
!                                                                       
! *** PURPOSE :                                                         
! ***  USER-PROGRAMMED SUBROUTINE GIVING:                               
! ***  (1)  SATURATION AS A FUNCTION OF PRESSURE ( SW(PRES) )           
! ***  (2)  DERIVATIVE OF SATURATION WITH RESPECT TO PRESSURE           
! ***       AS A FUNCTION OF EITHER PRESSURE OR SATURATION              
! ***       ( DSWDP(PRES), OR DSWDP(SW) )                               
! ***  (3)  RELATIVE PERMEABILITY AS A FUNCTION OF EITHER               
! ***       PRESSURE OR SATURATION ( REL(PRES) OR RELK(SW) )            
! ***                                                                   
! ***  CODE BETWEEN DASHED LINES MUST BE REPLACED TO GIVE THE           
! ***  PARTICULAR UNSATURATED RELATIONSHIPS DESIRED.                    
! ***                                                                   
! ***  DIFFERENT FUNCTIONS MAY BE GIVEN FOR EACH REGION OF THE MESH.    
! ***  REGIONS ARE SPECIFIED BY BOTH NODE NUMBER AND ELEMENT NUMBER     
! ***  IN INPUT DATA FILE FOR UNIT fINP.                                  
!                                                                       
      SUBROUTINE UNSAT (SW, DSWDP, RELK, PRES, KREG,SIOUT) 
      
      USE CONTRL
      USE M_SWCC
      USE SutraMSPrecision
    USE FRCONTRL
      IMPLICIT NONE
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!     E X A M P L E   C O D I N G   FOR                                 
!     MESH WITH TWO REGIONS OF UNSATURATED PROPERTIES USING             
!     THREE PARAMETER-UNSATURATED FLOW RELATIONSHIPS OF                 
!     VAN GENUCHTEN(1980)                                               
!        RESIDUAL SATURATION, SWRES, GIVEN IN UNITS {L**0}              
!        PARAMETER, AA, GIVEN IN INVERSE PRESSURE UNITS {m*(s**2)/kg}   
!        PARAMETER, VN, GIVEN IN UNITS {L**0}                           
!                                                                       
      INTEGER (I4B) :: &
        KREG
      REAL (DP) :: &
        SW, DSWDP, RELK, PRES,SIOUT

      !LOCAL VARIABLES
      REAL (DP) :: &
        SWRES, AA, VN, SWRM1, AAPVN, VNF, AAPVNN, DNUM, DNOM, SWSTAR 
!      REAL (DP) :: &
!        SWRES1, SWRES2,SWRES3, AA1, AA2, VN1, VN2 ,AA3, VN3, VN3
!!                                                                       
!!     DATA FOR REGION 1:                                                
!      DATA SWRES1/0.10E0/,   AA1/0.0005/,   VN1/1.70E0/ 
!      SAVE SWRES1, AA1, VN1 
!!     DATA FOR REGION 2:                                                
!      DATA SWRES2/0.1E0/,   AA2/0.00075/,   VN2/1.97E0/ 
!      SAVE SWRES2, AA2, VN2 
!!     DATA FOR REGION 2:                                                
!      DATA SWRES3/0.01E0/,   AA3/0.0018/,   VN3/2.9E0/ 
!      SAVE SWRES3, AA3, VN3 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
!                                                                       
! *** BECAUSE THIS ROUTINE IS CALLED OFTEN FOR UNSATURATED FLOW RUNS,   
! *** EXECUTION TIME MAY BE SAVED BY CAREFUL CODING OF DESIRED          
! *** RELATIONSHIPS USING ONLY INTEGER AND SINGLE PRECISION VARIABLES!  
! *** RESULTS OF THE CALCULATIONS MUST THEN BE PLACED INTO DOUBLE       
! *** PRECISION VARIABLES SW, DSWDP AND RELK BEFORE LEAVING             
! *** THIS SUBROUTINE.                                                  
!                                                                       
!                                                                       
!***********************************************************************
      ! KONG
!      GOTO 1800


!***********************************************************************
!                                                                       
!     SET PARAMETERS FOR CURRENT REGION, KREG                           
      GOTO (10, 20,30), KREG 
   10 SWRES = REAL(SWRES1)
      AA = REAL(AA1)
      VN = REAL(VN1) 
      GOTO 100 
   20 SWRES = REAL(SWRES2)
      AA = REAL(AA2)
      VN = REAL(VN2) 
      GOTO 100	  
   30 SWRES = REAL(SWRES3)                                                  
      AA = REAL(AA3)                                 !! 去掉了除以9800                                                            
      VN = REAL(VN3)				 					
  100 CONTINUE
    !IF (KREG.EQ.3) THEN        !!!!把表面节点参数输出到屏幕；此处1代表平台第一个点（与编号有关）
    !  write (*,*) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^'
    !    write (*,*) '   AA  =',  AA ,'   VN  =', VN    
    !  write (*,*) ' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
    !  END IF
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (1):                                                      
!     SW VS. PRES   (VALUE CALCULATED ON EACH CALL TO UNSAT)            
!     CODING MUST GIVE A VALUE TO SATURATION, SW.                       
!                                                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     THREE PARAMETER MODEL OF VAN GENUCHTEN(1980)                      
      SWRM1 = 1.E0 - SWRES 
      AAPVN = 1.E0 + (AA * ( - PRES) ) **VN 
      VNF = (VN - 1.E0) / VN 
      AAPVNN = AAPVN**VNF 
      SW = DBLE (SWRES + SWRM1 / AAPVNN) 
          !IF(SW .GT. 1-SIOUT) THEN
          !  SW = 1-SIOUT
          !ENDIF      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
 IF (IALSAT.EQ.1) THEN      
!***********************************************************************
!***********************************************************************
!.....SECTION (2):                                                      
!     DSWDP VS. PRES, OR DSWDP VS. SW   (CALCULATED ONLY WHEN IUNSAT=1) 
!     CODING MUST GIVE A VALUE TO DERIVATIVE OF SATURATION WITH         
!     RESPECT TO PRESSURE, DSWDP.                                       
!                                                                       
!  600 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      DNUM = AA * (VN - 1.E0) * SWRM1 * (AA * ( - PRES) ) ** (VN - 1.E0) 
      DNOM = AAPVN * AAPVNN 
      DSWDP = DBLE (DNUM / DNOM) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 END IF
      !     GOTO 1800 
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!***********************************************************************
!***********************************************************************
!.....SECTION (3):                                                      
!     RELK VS. P, OR RELK VS. SW   (CALCULATED ONLY WHEN IUNSAT=2)      
!     CODING MUST GIVE A VALUE TO RELATIVE PERMEABILITY, RELK.          
!                                                                       
! 1200 CONTINUE 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     GENERAL RELATIVE PERMEABILITY MODEL FROM VAN GENUCHTEN(1980)      
 !    SWSTAR = (SW - SWRES) / SWRM1 
!      RELK = DBLE (SQRT (SWSTAR) * (1.E0 - (1.E0 - SWSTAR** (1.E0 / VNF)) ** (VNF) ) **2)                                                 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                                                       
!***********************************************************************
!***********************************************************************
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
 1800 RETURN
	  
      END
