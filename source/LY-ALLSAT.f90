
    
    !!!ALLSAT calculates ice-water saturation, LIQSAT calculates liquid-water saturation, RELPERM calculates relative permeability.
    
    
    SUBROUTINE ALLSAT(SW,DSWDP,SLL,DSLDP,DSLDT,SII,DSIDP,DSIDT,RELK,&
                        PRES,TEMP,KREG,SIOUT,CONC,SII1,CFREEZ)  !yu
    USE CONTRL
    USE FRCONTRL
    USE ALLARR
    USE TIMES 
    USE SutraMSPrecision
    USE M_SWCC
     USE M_SFCC
    implicit none
     REAL (DP) :: &
        SW, DSWDP,SLL,DSLDP,DSLDT,SII,DSIDP,DSIDT, RELK, PRES,TEMP,SIOUT ,CONC,SII1,CFREEZ
     !LOCAL
     REAL (DP) :: DSLSATDT,SLSAT,TFREEZ
         REAL (DP) ::SWRES, AA, VN
    INTEGER (I4B) :: &
        KREG

     
    GOTO (10, 20,30), KREG 
    10 TFREEZ = REAL(TFREEZ1) 
      GOTO 100 
    20 TFREEZ = REAL(TFREEZ2) 
      GOTO 100 
    30 TFREEZ = REAL(TFREEZ3) 
100    CONTINUE
    
         GOTO (16, 26,36), KREG 
   16 SWRES = REAL(SWRES1)
      AA = REAL(AA1)
      VN = REAL(VN1) 
      GOTO 106 
   26 SWRES = REAL(SWRES2)
      AA = REAL(AA2)
      VN = REAL(VN2) 
      GOTO 106	  
   36 SWRES = REAL(SWRES3)                                                  
      AA = REAL(AA3)                                                                        
      VN = REAL(VN3)				 					
  106 CONTINUE         
   !    write(*,*)'CONC=',CONC
       
       
      
      IF ((IUNSAT.NE.0).AND.(PRES.LT.0D0)) THEN               
         CALL UNSAT (SW, DSWDP, RELK, PRES, KREG,SIOUT)     
                                                           
      ELSE                                                              
         SW = 1.D0                                                      
         IF (IALSAT.EQ.1) DSWDP = 0.D0                                  
      END IF                                                             
!    
      IF (TEMP.LT.TFREEZ.AND.IFREEZ.EQ.1) THEN   !!!!!!!!             
!........COMPUTE LIQUID-WATER SATURATION UNDER FULLY SATURATED          
!          CONDITIONS, SLSAT, AND ITS DERIVATIVE DSLSATDT.             
         CALL LIQSAT(SLSAT,DSLSATDT,TEMP,KREG,CONC,SII1,CFREEZ,SW)                          
     
!!!ice first
if(FIRS.eq.2)then
         SLL=MAX(SWRES,SW+SLSAT-1)                                              
         SII = SW - SLL                                                  
         IF (IALSAT.EQ.1) THEN                                          
            IF (SW+SLSAT-1.LE.SWRES) THEN                                         
!..............THE SATURATION/FREEZING MODEL USED IN SUTRA              
!              IMPLIES THE FOLLOWING RELATIONSHIPS BETWEEN              
!              THE LIQUID AND ICE SATURATION DERIVATIVES                
               DSLDP =   0                                           
               DSLDT = 0                                       
               DSIDP =   DSWDP                                          
               DSIDT =  0                                      
            ELSE                                                        
               DSLDP = DSWDP                                            
               DSLDT =  DSLSATDT                                           
               DSIDP = 0.D0                                             
               DSIDT = -DSLSATDT                                        
            END IF                                                      
         END IF    
 endif 
!!water first
if(FIRS.eq.1)then
          SLL = MIN(SLSAT, SW)                                            
         SII = SW - SLL                                                  
         IF (IALSAT.EQ.1) THEN                                          
            IF (SII.GT.0D0) THEN                                         
!..............THE SATURATION/FREEZING MODEL USED IN SUTRA              
!              IMPLIES THE FOLLOWING RELATIONSHIPS BETWEEN              
!              THE LIQUID AND ICE SATURATION DERIVATIVES                
               DSLDP = 0.D0                                             
               DSLDT = DSLSATDT                                         
               DSIDP = DSWDP                                            
               DSIDT = -DSLSATDT                                        
            ELSE                                                        
               DSLDP = DSWDP                                            
               DSLDT = 0.D0                                             
               DSIDP = 0.D0                                             
               DSIDT = 0.D0                                             
            END IF                                                      
         END IF 
endif     
         
      ELSE                                                              
!........NO FREEZING.                                                   
         SLL = SW                                                        
         SII = 0.D0                                                      
         IF (IALSAT.EQ.1) THEN                                          
            DSLDP = DSWDP                                               
            DSLDT = 0.D0                                                
            DSIDP = 0.D0                                                
            DSIDT = 0.D0                                                
         END IF                                                         
      END IF                                                            
!                                                                       
!.....COMPUTE RELATIVE PERMEABILITY, RELK, AS A FUNCTION OF LIQUID WATER   
  !      SATURATION, SL (IF NECESSARY).                                 
                                                
         IF (SLL.LT.1D0) THEN                                            
            CALL RELPERM(RELK,SLL,PRES,KREG,SII,SIOUT)                             
         ELSE                                                           
            RELK = 1.D0                                                 
         END IF    
                                                               
    RETURN 
    END SUBROUTINE ALLSAT     
    
    
    
                                                                         
   SUBROUTINE LIQSAT(SLSAT,DSLSATDT,TEMP,KREG,CONC,SII1,CFREEZ,SW) 
   USE CONTRL  
   USE ALLARR
   USE FRCONTRL
   USE M_SWCC
    USE M_SFCC
    USE SutraMSPrecision
    implicit none 
     REAL (DP) :: &
        SWRES, AA, VN, SWRM1, AAPVN, VNF, AAPVNN, DNUM, DNOM, SWSTAR   
   REAL (DP) :: &
        SLSAT,DSLSATDT,TEMP,TEMPREL,&
       SLSATRES,TLRES,TFREEZ,W,AASFCC,CONC,SII1,CFREEZ,SW,VNSFCC  
    INTEGER (I4B) :: &
        KREG  
     
     !  SET PARAMETERS FOR CURRENT REGION, KREG                           
      GOTO (12, 22,32), KREG 
   12 SLSATRES = REAL(SLSATRES1)
      TLRES = REAL(TLRES1)
      TFREEZ = REAL(TFREEZ1) 
      GOTO 102 
   22 SLSATRES = REAL(SLSATRES2)
      TLRES = REAL(TLRES2)
      TFREEZ = REAL(TFREEZ2) 
      GOTO 102	  
   32 SLSATRES = REAL(SLSATRES3)
      TLRES = REAL(TLRES3)
      TFREEZ = REAL(TFREEZ3) 				 					
102   CONTINUE
      
        GOTO (13, 23,33), KREG 
   13 SWRES = REAL(SWRES1)
      AA = REAL(AA1)
      VN = REAL(VN1) 
      GOTO 103 
   23 SWRES = REAL(SWRES2)
      AA = REAL(AA2)
      VN = REAL(VN2) 
      GOTO 103	  
   33 SWRES = REAL(SWRES3)                                                  
      AA = REAL(AA3)                                                                                    
      VN = REAL(VN3)				 					
  103 CONTINUE    
  

    TEMPREL = TEMP - TFREEZ
   AASFCC=AASFCC1/9810 
    
 IF(SPEFR.EQ.1)THEN    !!Piecewise linear function, considering the effect of solute on freezing point
          !IF (TEMPREL.LE.TLRES) THEN                                    
          !  SLSAT=SLSATRES                                        
          !ELSE                                                    
          !  SLSAT=(1.D0-SLSATRES)*(1.D0-TEMPREL/TLRES)+SLSATRES   
          !END IF     
          !
          !IF (IALSAT.EQ.1) THEN                                   
          !
          !  IF (TEMPREL.LT.TLRES) THEN                            
          !    DSLSATDT=0.D0                                       
          !  ELSE                                                  
          !    DSLSATDT=-(1.D0-SLSATRES)/TLRES                     
          !  END IF                                                
          !END IF                                                   
 
         IF(SII1.EQ.0)THEN
        IF(TEMPREL.GE.-KSF*CONC*1000)THEN     
           SLSAT=SW  
      !     DSLSATDT=0
           CFREEZ=CONC  
        ELSEIF(TEMPREL.LT.-KSF*CONC*1000.AND.TEMPREL.GT.(TLRES-KSF*CONC*1000))THEN   
           SLSAT=(1.D0-SLSATRES)*(1.D0-(TEMPREL+KSF*CONC*1000)/TLRES)+SLSATRES 
           DSLSATDT=-(1.D0-SLSATRES)/TLRES  
           CFREEZ=CONC  
         ELSEIF(TEMPREL.LE.(TLRES-KSF*CONC*1000))THEN  
             SLSAT=SLSATRES
             DSLSATDT=0.D0
             CFREEZ=CONC
        ENDIF    
  ELSEIF(SII1.GT.0)THEN         
          IF(TEMPREL.GE.-KSF*CFREEZ*1000)THEN     
          SLSAT=SW 
        !  DSLSATDT=0
          ELSEIF(TEMPREL.LT.-KSF*CFREEZ*1000.AND.TEMPREL.GT.(TLRES-KSF*CFREEZ*1000))THEN  
           SLSAT=(1.D0-SLSATRES)*(1.D0-(TEMPREL+KSF*CFREEZ*1000)/TLRES)+SLSATRES 
           DSLSATDT=-(1.D0-SLSATRES)/TLRES      
          ELSEIF(TEMPREL.LE.(TLRES-KSF*CFREEZ*1000))THEN  
             SLSAT=SLSATRES
             DSLSATDT=0.D0
          ENDIF 
  ENDIF
 
 ENDIF   
 
 IF(SPEFR.EQ.2)THEN           !!Exponential function
          W=0.5
         ! SLEXP=DEXP(-(TEMPREL/W)**2)
          SLSAT=(1.D0-SLSATRES)*DEXP(-(TEMPREL/W)**2)+SLSATRES
          DSLSATDT=-(1.D0-SLSATRES)*(2.D0*(TEMPREL)/(W*W))*DEXP(-(TEMPREL/W)**2)
ENDIF 
         !! 
         
  IF(SPEFR.EQ.3)THEN     !!Clapeyron equation , considering the effect of solute on freezing point    
         
   VNSFCC=2   
         
  IF(SII1.EQ.0)THEN
        IF(TEMPREL.GE.-KSF*CONC*1000)THEN     
           SLSAT=SW  
           CFREEZ=CONC  
        ELSEIF(TEMPREL.LT.-KSF*CONC*1000)THEN   
           SWRM1 = 1.E0 - SLSATRES 
           AAPVN = 1.E0 + (AASFCC * ( -(TEMPREL+CONC*KSF*1000)*334000*1000/273.15) ) **VNSFCC 
           VNF = (VNSFCC - 1.E0) / VNSFCC 
           AAPVNN = AAPVN**VNF 
           SLSAT=  DBLE (SLSATRES + SWRM1 / AAPVNN) 
           DNUM = AASFCC * (VNSFCC - 1.E0) * SWRM1 * (AASFCC * ( -(TEMPREL+CONC*KSF*1000)*334000*1000/273.15) ) ** (VNSFCC - 1.E0) 
           DNOM = AAPVN * AAPVNN 
           DSLSATDT=  DBLE (DNUM / DNOM*334000*1000/273.15)   
           CFREEZ=CONC  
        ENDIF    
  ELSEIF(SII1.GT.0)THEN         
          IF(TEMPREL.GE.-KSF*CFREEZ*1000)THEN     
          SLSAT=SW   
          ELSEIF(TEMPREL.LT.-KSF*CFREEZ*1000)THEN  
           SWRM1 = 1.E0 - SLSATRES 
           AAPVN = 1.E0 + (AASFCC * ( -(TEMPREL+CFREEZ*KSF*1000)*334000*1000/273.15) ) **VNSFCC 
           VNF = (VNSFCC - 1.E0) / VNSFCC 
           AAPVNN = AAPVN**VNF   
           SLSAT=  DBLE (SLSATRES + SWRM1 / AAPVNN) 
           DNUM = AASFCC * (VNSFCC - 1.E0) * SWRM1 * (AASFCC * ( -(TEMPREL+CFREEZ*KSF*1000)*334000*1000/273.15) ) ** (VNSFCC - 1.E0) 
           DNOM = AAPVN * AAPVNN 
           DSLSATDT=  DBLE (DNUM / DNOM*334000*1000/273.15)              
          ENDIF 
  ENDIF          
      
  !  IF (TEMPREL.LT.TLRES) THEN     !Truncated Clapeyron Equation
  !      SWRM1 = 1.E0 - SLSATRES 
  !         AAPVN = 1.E0 + (AASFCC * ( -(TLRES+CONC*KSF*1000)*334000*1000/273.15) ) **VNSFCC 
  !         VNF = (VNSFCC - 1.E0) / VNSFCC 
  !         AAPVNN = AAPVN**VNF 
  !         SLSAT=  DBLE (SLSATRES + SWRM1 / AAPVNN) 
  !          DSLSATDT=0.D0 
  !          
  !ENDIF
  !
  
 endif 
  
         
    RETURN 
    END SUBROUTINE LIQSAT
    
    
   SUBROUTINE RELPERM(RELK,SLL,PRES,KREG,SII,SIOUT )!! yu  
   USE CONTRL  
   USE FRCONTRL
   USE ALLARR
   USE SutraMSPrecision
   USE M_SWCC
    implicit none 
    REAL (DP) :: RELK,SLL,PRES,SII,SIOUT  !! yu 
   REAL (DP) :: & SWRES,VN,RKMIN
   INTEGER (I4B) :: KREG   
   !LOCAL
   REAL (DP) ::SLRM1,VNF,SLSTAR
    
  
    GOTO (10, 20,30), KREG 
   10 RKMIN = REAL(RKMIN1)
      SWRES = REAL(SWRES1)
      VN=REAL(VN1)
      GOTO 100 
   20 RKMIN = REAL(RKMIN2)
      SWRES = REAL(SWRES2)
      VN=REAL(VN2)
      GOTO 100	  
   30 RKMIN = REAL(RKMIN3)
      SWRES = REAL(SWRES3)
      VN=REAL(VN3)			 					
  100   CONTINUE
   
 
      SLRM1 = 1.D0 - SWRES                                          
      VNF = (VN - 1.D0)/VN                                          
   !!RELATIVE PERMEABILITY AS A FUNCTION OF LIQUID SATURATION      
     ! (RELK VS. SL)                                                 
      SLSTAR = (SLL - SWRES)/SLRM1                                  
      RELK = DSQRT(SLSTAR)*    &                                     
         (1.D0 - (1.D0 - SLSTAR**(1.D0/VNF))**VNF)**2               
    !!  ENFORCE THE MINIMUM VALUE OF RELK USING SUTRA PARAMETER RKMIN 
      IF (RELK .LT. RKMIN) RELK=RKMIN                               
  

   
   RETURN
   END SUBROUTINE RELPERM