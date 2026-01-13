!     SUBROUTINE        B  U  D  G  E  T       SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO CALCULATE AND OUTPUT FLUID MASS AND SOLUTE MASS OR            
! ***  ENERGY BUDGETS.                                                  
!                                                                       
      SUBROUTINE BUDGET (ML, IBCT, VOL, SW, DSWDP,SLL,DSLDP,DSLDT,SII,DSIDP,DSIDT, RHO, QIN, PVEC, &           !!!
                         PM1, PBC, QPLITR, IPBC, IQSOP, UVEC,PITER,UITER, UM1, UM2, UIN, QUIN,    &
                         IQSOU, UBC, IUBC, CS1, CS2, CS3, SL, SR, GNUPCJ,SIOUT,SII1,CFREEZ) ! CHENGJI 2013-09-02
      USE PARAMS 
      USE CONTRL 
      USE FRCONTRL
      USE MODSOR 
      USE FUNITS 
      USE DIMS
      USE TIMES
      USE SutraStorage, ONLY : SpecifiedPBC, &
                               MultiSpeciesBC,Y,X,Z
      use SutraMSPrecision
      USE SutraZoneModule
      USE MSErrorHandler
	  USE M_OUT
      USE SOLVI ! Chengji 2015-04-14
      USE BCNAME
        USE SALTPARAMS      !!!
      implicit none
	  
!============================================= CHENGJI 2013-09-02	
      REAL(DP):: GNUPCJ  
      DIMENSION GNUPCJ(NBCN)
!================================================================  
	  
      CHARACTER(13) UNAME (2) 
      integer (I4B) :: &
        IBCT, ML
      integer (I4B) :: &
        IQSOP (NSOP)
      integer (I4B) :: &
        IQSOU (MNSOU, NSPE)                                               
      integer (I4B) :: &
        IPBC (NBCN), IUBC (NBCN, NSPE)
      real (kind=8) :: &
        QIN (NN), UIN (NN, NSPE)
      real (DP) :: &
        QUIN (NN, NSPE)
      real (DP) :: &
        UBC (NBCN, NSPE),       &
        QPLITR (NBCN), PBC (NBCN)                                         
      real (DP) :: &
        VOL (NN), PVEC (NN), UVEC (NN, NSPE), SW (NN),SLL(NN),DSLDP(NN),DSLDT(NN),SII(NN),DSIDP(NN),DSIDT(NN),&      !!!
        DSWDP (NN), RHO (NN), PM1 (NN), UM1 (NN, NSPE), UM2 (NN,NSPE), &
        CS1 (NN, NSPE), CS2 (NN), CS3 (NN), SL (NN), SR (NN),PITER(NN),UITER(NN,NSPE),SIOUT(NN),SII1(NN),CFREEZ(NN),RHO1(NN),RHOCHANGE(NN),UCHANGE(NN)
 !     real (DP) :: &    ! COMMENTED BY CHENGJI
 !       STUNT (NSPE), STUPT (NSPE), STUTOT (NSPE) 
      !locals 
      integer (I4B) :: &
        I, K, &
        IP, IQU, IU, IUP, &
        MN,PPP
      integer (I4B) :: &
        INEGCT, IQP, NSOPI, NSOUI
      integer (I4B) :: imap
      real (DP) :: &
        CST, CWT, &
        RELK, &
        DNSTOT, FLDTOT, P0FTOT, P1FTOT, P0STOT, P1STOT, &
        QIUTOT, QQUTOT, QPUTOT, QULTOT, SLDTOT, &
      FIUTOT,USCPTOT,USLPTOT,USDLTOT,USLUTOT,USCPITOT,FLPTOT,FLUTOT,STETOT,STETOTSL,QFSTOT,FLUTOT1,FLUTOT2,&         !!!!!!!
        DUDT, EPRSV, ESRV, &
        QU, QPU
	  
      REAL(DP):: STPPOS, STPNEG, STUPOS(NSPE), STUNEG(NSPE), QINPOS, QINNEG, &
				 STPTOT, STUTOT, STFPOS, STFNEG, STFTOT, QINTOT, &
				 QPLPOS, QPLNEG, QPLTOT,      &
				 QFFPOS, QFFNEG, QFFTOT,      &
				 TERM,                        & ! Chengji 2013-09-24
                 TERM11,TERMSPE1,TERMSPE2,    & ! Chengji 2015-04-14
                 ACTFMB,ERFMBA,ERFMBR,        & ! Chengji 2015-08-06
                 QPSOLUTE                       ! Chengji 2016-01-05

      
       REAL(DP)::SCLPPOS, SCLPNEG,SLPPOS,SLPNEG ,SCIPPOS,SCIPNEG,SIPPOS ,SIPNEG,SLUTOT(NSPE),SIUTOT(NSPE),&            !!!
              SDLUPOS(NSPE),SDLUNEG(NSPE),SDLUTOT(NSPE),SDIUTOT(NSPE),SLUPOS(NSPE),SLUNEG(NSPE) ,&
              SDIUPOS(NSPE),SDIUNEG(NSPE),SIUPOS(NSPE), SIUNEG(NSPE),&
              SCLPTOT,SLPTOT,SCIPTOT,SIPTOT,STWTOT,&
             DPDT,EICE,WATOT,SATOT,SLTOT,SITOT,SATOTSW,SMTOT,SIOUTTOT
      integer (I4B) :: ro 
  real (DP) ::kr,j,rows,cols,SALTHEIGHT1,SALTHEIGHT2,QW1,QW2,QW3,QW4,max_valks,min_valks,max_valk,min_valk,max_vals,min_vals
    
      DATA UNAME (1) / 'CONCENTRATION' / , UNAME (2) / ' TEMPERATURE ' / 
      SAVE UNAME 
!                                                                       
!.....SAVE VALUES FOR ENERGY TRANSPORT                                  
      CST = CS 
      CWT = CW 
!                                                                       
      MN = 2 
      IF(IALSAT.NE.0) IALSAT=1
      IF (ME.EQ. - 1) MN = 1 
      WRITE (fLST, 10) 
   10 FORMAT(1H1) 
!.....SET UNSATURATED FLOW PARAMETERS, SW(I) AND DSWDP(I)               
 
  PPP=48*1800/DELT
  QW1=0
   QW2=0
    QW3=0
     QW4=0
  SALTHEIGHT1=0
    SALTHEIGHT2=0
       DO 30 I=1,NN
        CALL ALLSAT(SW(I),DSWDP(I),SLL(I),DSLDP(I),DSLDT(I),SII(I),DSIDP(I),DSIDT(I),RELK,PITER(I),UITER(I,1),NodeMap(I),SIOUT(I),UITER(I,2),SII1(I),CFREEZ(I))    
        
   !     
30     CONTINUE
       write(012301,'(I,101E15.6e3)')it,(uvec(k,1),k=5011,5511,5)
     write(012302,'(I,101E15.6e3)')it,(sll(k),k=5011,5511,5)
         write(012303,'(I,101E15.6e3)')it,(uvec(k,2),k=5011,5511,5)      
  
!.....CALCULATE COMPONENTS OF FLUID MASS BUDGET                         
   40 IF (ML - 1) 50, 50, 1000 
50    CONTINUE 
 
      
      
!====================================== Added according to SUTRA, Chengji 2013-09-24
      STPPOS = 0.D0                                            
      STPNEG = 0.D0                                            
      DO K=1,NSPE
	  STUPOS(K) = 0.D0   ! Influx for Kth species                                        
      STUNEG(K) = 0.D0   ! Outflux for Kth species
      ENDDO
      QINPOS = 0.D0                                            
      QINNEG = 0.D0
      
      STPTOT = 0.D0
      STUTOT = 0.D0
      STFPOS = 0.D0
      STFNEG = 0.D0
      STFTOT = 0.D0
      QINTOT = 0.D0
      
  !!*************Freeze new parameters↓↓****************************************    
           SCLPPOS = 0D0                                                    
           SCLPNEG = 0D0                                                    
           SLPPOS = 0D0                                                     
           SLPNEG = 0D0                                                     
           SCIPPOS = 0D0                                                    
           SCIPNEG = 0D0                                                    
           SIPPOS = 0D0                                                     
           SIPNEG = 0D0                                                     
     !.....SIP IS CHANGE IN MASS OF ICE DUE TO P CHANGE;                    
     !        SIP IS NON-ZERO ONLY DURING UNSATURATED FREEZE/THAW           
          DO K=1,NSPE
           SDLUPOS(K) = 0D0                                                    
           SDLUNEG(K) = 0D0                                                    
           SLUPOS(K) = 0D0                                                     
           SLUNEG(K) = 0D0                                                     
     !.....SLU IS CHANGE IN MASS OF LIQUID DUE TO U (I.E. TEMPERATURE) CHANG
     !        SLU IS NON-ZERO ONLY DURING FREEZE/THAW                       
           SDIUPOS(K) = 0D0                                                    
           SDIUNEG(K) = 0D0                                                    
           SIUPOS(K) = 0D0                                                     
           SIUNEG(K) = 0D0                                                     
     !.....SIU IS CHANGE IN MASS OF ICE DUE TO U (I.E. TEMPERATURE)CHANGE    
     
           ENDDO
      
	  
!***************************************************************************************** 
!!***************************************************************************************** 

	  
      DO 100 I=1,NN  
   
      !!pressure change
      DPDT = (PVEC(I)-PM1(I))/DELTP                                              !!DPDT          
                                                                         
        TERM = (1-ISSFLO/2)*RHO(I)*VOL(I)*SLL(I)*NodeData(NodeMap(i))%sop*DPDT 
        SCLPPOS = SCLPPOS + MAX(0D0, TERM)                               
        SCLPNEG = SCLPNEG + MIN(0D0, TERM)  
        SCLPTOT=SCLPPOS+SCLPNEG
        TERM = (1-ISSFLO/2)*RHO(I)*VOL(I)*NodeData(NodeMap(i))%por*DSLDP(I)*DPDT       
        SLPPOS = SLPPOS + MAX(0D0, TERM)                                 
        SLPNEG = SLPNEG + MIN(0D0, TERM)    
        SLPTOT=SLPPOS+SLPNEG
        TERM = (1-ISSFLO/2)*RHO(I)*VOL(I)*SII(I)*NodeData(NodeMap(i))%sopi*DPDT        
        SCIPPOS = SCIPPOS + MAX(0D0, TERM)                               
        SCIPNEG = SCIPNEG + MIN(0D0, TERM)                               
        TERM = (1-ISSFLO/2)*VOL(I)*NodeData(NodeMap(i))%por*RHO(I)*DSIDP(I)*DPDT  !! 
        SIPPOS = SIPPOS + MAX(0D0, TERM)                                 
        SIPNEG = SIPNEG + MIN(0D0, TERM)                                 
         SIPTOT=SIPPOS+SIPNEG
         
        STPPOS = SCLPPOS + SLPPOS + SCIPPOS + SIPPOS        ! Influx caused by pressure change                 
        STPNEG = SCLPNEG + SLPNEG + SCIPNEG + SIPNEG        ! Outflux caused by pressure change             
      
      
  !!temperature or solute change    
 
      
    !!Assuming the first U must represent temperature  
       TERM = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SLL(I)*DRWDU(1)*VOL(I)*(UVEC(I,1)-UM1(I,1))/DELTU            
       SDLUPOS(1) = SDLUPOS(1) + MAX(0D0, TERM)                            
       SDLUNEG(1) = SDLUNEG(1) + MIN(0D0, TERM)                            
       SDLUTOT(1)=SDLUPOS(1)+SDLUNEG(1)
       TERM = (1-ISSFLO/2)*VOL(I)*NodeData(NodeMap(i))%por*RHO(I)*DSLDT(I)*(UVEC(I,1)-UM1(I,1))/DELTU    !    
       SLUPOS(1) = SLUPOS(1) + MAX(0D0, TERM)                              
       SLUNEG(1) = SLUNEG(1) + MIN(0D0, TERM)                              
       SLUTOT(1)= SLUPOS(1)+ SLUNEG(1)                                                            
      !.SET ICE DENSITY CHANGE WITH T EQUAL TO THAT OF LIQUID         
      !    DUE TO THE 'NO VOLUME CHANGE' ASSUMPTION FOR FREEZE/THAW.  
                                              
                                                                     
       TERM = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SII(I)*DRWDU(1)*VOL(I)*(UVEC(I,1)-UM1(I,1))/DELTU         
       SDIUPOS(1) = SDIUPOS(1) + MAX(0D0, TERM)                            
       SDIUNEG(1) = SDIUNEG(1) + MIN(0D0, TERM)               
       SDIUTOT(1)=SDIUNEG(1)+SDIUPOS(1)
       TERM = (1-ISSFLO/2)*VOL(I)*NodeData(NodeMap(i))%por                             
       TERM = TERM*RHO(I)*DSIDT(I)*(UVEC(I,1)-UM1(I,1))/DELTU       !  
       SIUPOS(1) = SIUPOS(1) + MAX(0D0, TERM)                              
       SIUNEG(1) = SIUNEG(1) + MIN(0D0, TERM)                              
        SIUTOT(1)= SIUPOS(1)+SIUNEG(1)
        
       STUPOS(1) = SDLUPOS(1) + SLUPOS(1)  + SIUPOS(1)                  
       STUNEG(1) = SDLUNEG(1) + SLUNEG(1)  + SIUNEG(1)    
       STWTOT=SCLPTOT+SLPTOT+SIPTOT+SDLUTOT(1)+SLUTOT(1)+SIUTOT(1)
      
     !!If a solute is present, it alters the quality of the water.
       IF (NSPE.GT.1) THEN
       
           
      DO K=2,NSPE
      TERM = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SLL(I)*DRWDU(K)*VOL(I)*(UVEC(I,K)-UM1(I,K))/DELTU  
       SDLUPOS(K) = SDLUPOS(K) + MAX(0D0, TERM)                            
       SDLUNEG(K) = SDLUNEG(K) + MIN(0D0, TERM)                            
       SDLUTOT(K)=SDLUPOS(K)+SDLUNEG(K)
       
       
         TERM = (1-ISSFLO/2)*NodeData(NodeMap(i))%por*SII(I)*DRWDU(K)*VOL(I)*(UVEC(I,K)-UM1(I,K))/DELTU    
       SDIUPOS(K) = SDIUPOS(K) + MAX(0D0, TERM)                            
       SDIUNEG(K) = SDIUNEG(K) + MIN(0D0, TERM)               
       SDIUTOT(K)=SDIUNEG(K)+SDIUPOS(K)
       
      ENDDO  
      
      
      
       ENDIF
      
       IF(NSPE==1) STWTOT=SCLPTOT+SLPTOT+SIPTOT+SDLUTOT(1)+SLUTOT(1)+SIUTOT(1)
       IF(NSPE==2)STWTOT=SCLPTOT+SLPTOT+SCIPTOT+SIPTOT+ SDLUTOT(1)+SLUTOT(1)+SDIUTOT(1)+SIUTOT(1) +SDIUTOT(2)+SDLUTOT(2)
     
       
       
       
	  TERM = QIN(I)          
      QINPOS = QINPOS + MAX(0D0, TERM)  ! Influx by flux boundary                       
      QINNEG = QINNEG + MIN(0D0, TERM)  ! Outflux by flux boundary
100   CONTINUE
  
      STPTOT = STPPOS + STPNEG  ! Net flux caused by pressure change                  
      
      
      DO K=1,NSPE
	  STUTOT = STUTOT + STUPOS(K) + STUNEG(K) ! Net flux caused by concentration change         
      STFPOS = STPPOS + STUPOS(K) ! Influx caused by pressure change and concentration change  
      STFNEG = STPNEG + STUNEG(K) ! Outflux caused by pressure change and concentration change
      ENDDO
      STFTOT = STPTOT + STUTOT  ! Net flux caused by pressure change and concentration change                               
      QINTOT = QINPOS + QINNEG  ! Net flux caused by flux boundary
!
  
      QPLPOS = 0.D0                       
      QPLNEG = 0.D0 
      QPLTOT = 0.D0
                       
      DO 200 IP=1,NPBC                   
      I=IABS(SpecifiedPBC(IP)%node)                   
      TERM = GNUPCJ(IP)*(SpecifiedPBC(IP)%P-PVEC(I))
      QPLITR(IP) = GNUPCJ(IP)*(SpecifiedPBC(IP)%P-PVEC(I))
      QPLPOS = QPLPOS + MAX(0D0, TERM)  ! Influx caused by pressure-precribed boundary nodes
      QPLNEG = QPLNEG + MIN(0D0, TERM)  ! Outflux caused by pressure-precribed boundary nodes
  200 CONTINUE
                       
      QPLTOT = QPLPOS + QPLNEG   ! Net flux caused by pressure-prescribed boundary nodes        
      QFFPOS = QINPOS + QPLPOS           
      QFFNEG = QINNEG + QPLNEG           
      QFFTOT = QINTOT + QPLTOT 
      
    
																			 																	 																																					   
      IF (IBCT.EQ.4) GOTO 600 
      NSOPI = NSOP - 1 
      INEGCT = 0 
      DO 500 IQP = 1, NSOPI 
         I = IQSOP (IQP) 
         IF (I) 325, 500, 500 
  325    INEGCT = INEGCT + 1 
         IF (INEGCT.EQ.1) WRITE (fLST, 350) 
  350 FORMAT(///22X,'TIME-DEPENDENT FLUID SOURCES OR SINKS'//22X,       &
     &   ' NODE',5X,'INFLOW(+)/OUTFLOW(-)'/37X,'  (MASS/SECOND)'//)     
!         WRITE (fLST, 450) - I, QIN ( - I) 
!         WRITE (3, '(4E15.6e3)')  QIN ( - I)    !!Yu																								 											 
!  450 FORMAT(18X,I9,10X,1PD15.7) 
500 END DO 
! 
    
    !!!!!!!!!!!!!!!!!!!!!!!! water Variations
    IF(NSPE.EQ.1)THEN
      IF(IT.EQ.1)THEN
   WRITE(1103,9877) 
 9877 FORMAT(8X,'IT     SCLPTOT         SLPTOT         SIPTOT         SDLUTOT         SLUTOT         SIUTOT&
&     QINTOT         QPLTOT         STWTOT         QFFTOT ')
    ENDIF 
    WRITE(1103,'(1I9,20E15.6E3)')IT, SCLPTOT,SLPTOT,SIPTOT,SDLUTOT(1),SLUTOT(1),SIUTOT(1),QINTOT,QPLTOT,STWTOT,QFFTOT
    ENDIF   
  
   IF(NSPE.EQ.2)THEN
      IF(IT.EQ.1)THEN
   WRITE(1103,9870) 
 9870 FORMAT(8X,'IT     SCLPTOT         SLPTOT         SCIPTOT         SIPTOT     SDLUTOT(1)      SLUTOT      SDIUTOT(1)      SIUTOT  &
&     SDLUTOT(2)      SDIUTOT(2)       QINTOT         QPLTOT         STWTOT         QFFTOT ')
    ENDIF 
    WRITE(1103,'(1I9,20E15.6E3)')IT, SCLPTOT,SLPTOT,SCIPTOT,SIPTOT,SDLUTOT(1),SLUTOT(1),SDIUTOT(1),SIUTOT(1),SDLUTOT(2),SDIUTOT(2),QINTOT,QPLTOT,STWTOT,QFFTOT
    ENDIF  
    
 
    
  600 IF (NPBC.EQ.0) GOTO 800 
!      WRITE (fLST, 650)   ! CHENGJI 2013-09-03
  650 FORMAT(///22X,'FLUID SOURCES OR SINKS DUE TO SPECIFIED PRESSURES',&
     &   //22X,' NODE',5X,'INFLOW(+)/OUTFLOW(-)'/37X,'  (MASS/SECOND)'/)
! **************************************** XXH *******************************************																							 
!	 WRITE(20231,*) 'FLUID DUE TO SPECIFIED PRESSURES ("NODE","Q kg/s","FLUX mm/d"), IT=',IT
! **************************************** XXH *******************************************																												  															 																						  
      DO 700 IP = 1, NPBC 
!         I = IABS ( IPBC (IP) ) 
         I = IABS ( SpecifiedPBC(IP)%node ) 
!!         WRITE (fLST, 450) I, QPLITR (IP)
!	 WRITE(20231,20093) I, QPLITR (IP), QPLITR (IP)/Sarea(IP)/1000*1000*3600*24                
! 20093	FORMAT (I8,2E17.8)						 
!         WRITE (2,'(4E15.6e3)')  QPLITR (IP)  ! CHENGJI 2013-09-13, write node-wise fluid mass to the '.dat' file
  700 END DO 
!********************************************************************************************
!********************************************************************************************











!********************************************************************************************
!********************************************************************************************
!.....CALCULATE COMPONENTS OF ENERGY OR SOLUTE MASS BUDGET              
  800 IF (ML - 1) 1000, 5500, 1000 
 1000 CONTINUE 
!                                                                       
      DO 5000 K = 1, NSPE                
!.....GLOBAL COUNTER FOR SPECIES - DEFINED IN MODULE PARAMS             
         KSP = K 
!                                                                       
!.....WRITE SPECIES INFORMATION TO OUTPUT FILE                          
         WRITE (fLST, 1010) K, trim(adjustl(SPNAME(K)))
1010     FORMAT(//11X,'[SPECIES ',I3,'] - ',A) 
!                                                                       

  IF(K.EQ.1)THEN       !!!能量平衡
!.....ZERO VARIABLE FOR MASS BALANCE CALCULATIONS                        
         FLDTOT = 0.D0 
         SLDTOT = 0.D0 
         DNSTOT = 0.D0 
         P1FTOT = 0.D0 
         P1STOT = 0.D0 
         P0FTOT = 0.D0 
         P0STOT = 0.D0 
         QQUTOT = 0.D0 
         QIUTOT = 0.D0 
         
         FIUTOT=0
         USCPTOT=0
         USLPTOT=0
         USDLTOT=0
         USLUTOT=0
         FLPTOT=0
         FLUTOT=0
         STETOT=0
         QFSTOT=0
         
         FLUTOT1=0
         FLUTOT2=0
         
!.....SET ADSORPTION PARAMETERS                                         
!      IF (ME.EQ. - 1.AND.ADSMOD (K) .NE.'NONE      ') CALL ADSORB (CS1, CS2, CS3, SL, SR, UVEC)                                           
      IF (ME.EQ. - 1.AND.ADSMOD (K) .NE.'NONE      ') CALL ADSORB ()
!                                                                       
!.....SET APPROPRIATE VALUES FOR CS AND CW                              
         IF (K.EQ.NESP) THEN                 
            CS = CST 
            CW = CWT  
         ELSE                        
            CS = 0.0D0 
            CW = 1.0D0 
         ENDIF 
!                                                                       
         DO 1300 I = 1, NN 
            imap = NodeMap(i)
            ESRV = NodeData(imap)%por * SLL(I) * RHO(I) * VOL(I) 
            EPRSV = (1.D0 - NodeData(imap)%por ) * NodeData(imap)%rhos * VOL(I) 
            EICE = -HTLAT + CI*UVEC(I,1)                                           !!EICE
            DUDT = (1 - ISSTRA) * ( UVEC(I, K) - UM1(I, K) ) / DELTU 
            FLDTOT = FLDTOT + ESRV * CW * DUDT 
            SLDTOT = SLDTOT + EPRSV * CS1(I, K) * DUDT          
!======= FLDTOT + SLDTOT = total rate of change in stored solute mass in the region due to change in concentration
!======= FLDTOT + SLDTOT equals to (5.15a) on page 112 of the manual		
 !!!
 FIUTOT=FIUTOT+ NodeData(imap)%por * SII(I) * RHOI * VOL(I) *CI* DUDT          
USCPTOT=USCPTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,1)*RHO(I)*SLL(I)*NodeData(NodeMap(i))%sopi*(PVEC(I)-PM1(I))/DELTP 
USLPTOT=USLPTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,1)*RHO(I)*NodeData(imap)%por*DSLDP(I)*(PVEC(I)-PM1(I))/DELTP
USDLTOT=USDLTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,1)*NodeData(imap)%por*SLL(I)*DRWDU(1)*(UVEC(I,1)-UM1(I,1))/DELTU
USLUTOT=USLUTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,1)*NodeData(imap)%por*RHO(I)*DSLDT(I)*(1-ISSTRA)*(UVEC(I,1)-UM1(I,1))/DELTU
FLPTOT= FLPTOT+  (1-ISSFLO/2)*VOL(I)*CI*UVEC(I,1)*NodeData(imap)%por*RHO(I)*DSIDP(I)*(PVEC(I)-PM1(I))/DELTP
FLUTOT =FLUTOT+ VOL(I)*EICE*NodeData(imap)%por*1000*DSIDT(I)*(1-ISSTRA)*(UVEC(I,1)-UM1(I,1))/DELTU
STETOT =FLDTOT+SLDTOT+FIUTOT+ USCPTOT+USLPTOT+USDLTOT+USLUTOT+ FLPTOT+FLUTOT

FLUTOT1 =FLUTOT1+ VOL(I)*CI*UVEC(I,1)*NodeData(imap)%por*1000*DSIDT(I)*(1-ISSTRA)*(UVEC(I,1)-UM1(I,1))/DELTU
FLUTOT2 =FLUTOT2+ VOL(I)*(-HTLAT)*NodeData(imap)%por*1000*DSIDT(I)*(1-ISSTRA)*(UVEC(I,1)-UM1(I,1))/DELTU

            DNSTOT = DNSTOT + CW * UVEC(I, K) * (1 - ISSFLO / 2) *      &
              VOL(I) * ( RHO(I) * (SLL(I) * NodeData(imap)%sop + NodeData(NodeMap(i))%por * DSLDP(I) ) * &
              ( PVEC(I) - PM1(I) ) / DELTP + NodeData(imap)%por * SLL(I) *  &
              DRWDU (K) * ( UM1(I, K) - UM2(I, K) ) / DLTUM1) 
!======= DNSTOT = total rate of change in stored solute mass in the region due to change in stored fluid mass
!======= DNSTOT equals to (5.15b) on page 112 of the manual          
!         IF(K.EQ.NSPE) THEN
!          IF(MOD(I,NN1).EQ.1) THEN
!             WRITE(20098,'(2I8,10E15.7)') K, I, DUDT, (PVEC(I)-PM1(I))/DELTP, (UM1(I,K)-UM2(I,K))/DLTUM1,          &         
!                                                ESRV*CW*DUDT, EPRSV*CS1(I,K)*DUDT, CW*UVEC(I,K)*(1-ISSFLO/2)*      &
!                                 VOL(I)*(RHO(I)*(SW(I)*NodeData(imap)%sop+NodeData(NodeMap(i))%por*DSWDP(I))*      &
!                                                             (PVEC(I)-PM1(I))/DELTP+NodeData(imap)%por*SW(I)*      &
!                                                                        DRWDU(K)*(UM1(I,K)-UM2(I,K))/DLTUM1)
!          ENDIF
!         ENDIF
 			  
!            P1FTOT = P1FTOT + ESRV * PRODF1 (K) * UVEC(I, K) 
!            P1STOT = P1STOT + EPRSV * PRODS1 (K) * ( SL(I) * UVEC(I, K) + SR(I) )                                                  
!            P0FTOT = P0FTOT + ESRV * PRODF0 (K) 
!            P0STOT = P0STOT + EPRSV * PRODS0 (K) 
            P1FTOT = P1FTOT + ESRV  * ProdSorp(imap)%prodf1(K) * UVEC(I, K) 
!======= P1FTOT is the total rate of first-order solute mass production in the fluid, (5.16a) on page 113			
            P1STOT = P1STOT + EPRSV * ProdSorp(imap)%prods1(K) * ( SL(I) * UVEC(I, K) + SR(I) )  
!======= P1STOT is the total rate of first-order adsorbate production in the fluid, (5.16b) on page 113			
            P0FTOT = P0FTOT + ESRV  * ProdSorp(imap)%prodf0(K) 
            P0STOT = P0STOT + EPRSV * ProdSorp(imap)%prods0(K) 
!======= P0FTOT + P0STOT is the zero-order production of solute and adsorbate mass production in the fluid and solid matrix
!======= (5.17) on page 113
            QQUTOT = QQUTOT + QUIN(I, K) 
!======= QQUTOT is the diffusive-dispersive source of solute mass, (5.20) on page 113
            IF ( QIN(I) ) 1200, 1200, 1250 
 1200       QIUTOT = QIUTOT + QIN(I) * Ci * UVEC(I, K) 
            GOTO 1300 
 1250       QIUTOT = QIUTOT + QIN(I) * Ci * UIN(I, K)
!======= QIUTOT is the solute mass change at fluid source node, (5.18) on page 113 
   1300     END DO 
!
   
      

         
         
         
!.....MASS CHANGES DUE TO SPECIFIED PRESSURES                           
         QPUTOT = 0.
         
!         DO 1500 IP = 1, NPBC 
!            IF (QPLITR (IP) ) 1400, 1400, 1450 
! 1400       I = IABS ( SpecifiedPBC(IP)%node ) 
!            QPUTOT = QPUTOT + QPLITR (IP) * CW * UVEC(I, K)
!            QPSOLUTE = QPLITR (IP) * CW * UVEC(I, K)
!            GOTO 1500 
! 1450       QPUTOT = QPUTOT + QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
!            QPSOLUTE = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
 !1500    END DO      
         DO 1500 IP = 1, NPBC 
          IF (QPLITR(IP).LE.0) THEN 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * UVEC(I, K)
            QPSOLUTE = QPLITR (IP) * CW * UVEC(I, K)
          ELSEIF (QPLITR(IP).GT.0) THEN 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
            QPSOLUTE = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
          ENDIF
    1500 END DO         
!======= QPUTOT is the solute mass change at specified pressure node, (5.19) on page 113  
 

!                                                                       
!.....MASS CHANGES DUE TO SPECIFIED CONCENTRATIONS                      
         QULTOT = 0.D0 
         IF (NUBC (K) .EQ.0) GOTO 1540 
         DO 1510 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
            QPLITR (IUP) = GNUU (K) * ( MultiSpeciesBC(K)%SpecifiedU(IU)%U - UVEC(I, K) ) 
            QULTOT = QULTOT + QPLITR (IUP) 
!======= QULTOT is the solute mass due to specified concentration conditions, (5.21) on page 114
 1510    END DO 
!======= FLDTOT+SLDTOT+DNSTOT should be close to P1FTOT+P1STOT+P0FTOT+P0STOT+QIUTOT+QPUTOT+QQUTOT+QULTOT                                                                      
IF(K.EQ.1)THEN

ENDIF

 1540    IF (K.EQ.NESP) GOTO 1615 
		
!======= CHENGJI 2013-09-03, output total solute mass to 'dat' file		
!		WRITE(3,'(I6,3E15.7)') IT,FLDTOT+SLDTOT+DNSTOT,P1FTOT+P1STOT+P0FTOT+P0STOT+QIUTOT+QPUTOT+QQUTOT+QULTOT
!                                                                       
!.....OUTPUT SOLUTE MASS BUDGET                                         
 1550    WRITE (fLST, 1600) IT, FLDTOT, SLDTOT, DNSTOT, P1FTOT, P1STOT,   &
         P0FTOT, P0STOT, QIUTOT, QPUTOT, QQUTOT, QULTOT                 
 1600 FORMAT(//11X,'S O L U T E   B U D G E T      AFTER TIME STEP ',I5,&
     &   ',   IN (SOLUTE MASS/SECOND)'///11X,1PD15.7,5X,'NET RATE OF ', &
     &   'INCREASE(+)/DECREASE(-) OF SOLUTE DUE TO CONCENTRATION CHANGE'&
     &   /11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/DECREASE(-) OF ',     &
     &   'ADSORBATE'/11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/',         &
     &   'DECREASE(-) OF SOLUTE DUE TO CHANGE IN MASS OF FLUID'//11X,   &
     &   1PD15.7,5X,'NET FIRST-ORDER PRODUCTION(+)/DECAY(-) OF SOLUTE'  &
     &   /11X,1PD15.7,5X,'NET FIRST-ORDER PRODUCTION(+)/DECAY(-) OF ',  &
     &   'ADSORBATE'/11X,1PD15.7,5X,'NET ZERO-ORDER PRODUCTION(+)/',    &
     &   'DECAY(-) OF SOLUTE'/11X,1PD15.7,5X,'NET ZERO-ORDER ',         &
     &   'PRODUCTION(+)/DECAY(-) OF ADSORBATE'/11X,1PD15.7,5X,          &
     &   'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH FLUID SOURCES AND SINKS'&
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH ',      &
     &   'INFLOWS OR OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'          &
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF SOLUTE THROUGH ',      &
     &   'SOLUTE SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-)'&
     &  ,' OF SOLUTE AT POINTS OF SPECIFIED CONCENTRATION')             
         GOTO 1645 
!                                                                       
!.....OUTPUT ENERGY BUDGET                                              
 1615    WRITE (fLST, 1635) IT, FLDTOT, SLDTOT, DNSTOT, P0FTOT, P0STOT,   &
         QIUTOT, QPUTOT, QQUTOT, QULTOT                                 
 !        WRITE(3,'(2I6,4E15.7)') IT,K, QIUTOT, QPUTOT,QQUTOT, QULTOT      !!Xu																							
 1635 FORMAT(//11X,'E N E R G Y   B U D G E T      AFTER TIME STEP ',I5,&
     &   ',   IN (ENERGY/SECOND)'///11X,1PD15.7,5X,'NET RATE OF ',      &
     &   'INCREASE(+)/DECREASE(-) OF ENERGY IN FLUID DUE TO TEMPERATURE'&
     & ,' CHANGE'/11X,1PD15.7,5X,'NET RATE OF INCREASE(+)/DECREASE(-) ',&
     &   'OF ENERGY IN SOLID GRAINS'/11X,1PD15.7,5X,'NET RATE OF ',     &
     &   'INCREASE(+)/DECREASE(-) OF ENERGY IN FLUID DUE TO CHANGE IN ',&
     &   'MASS OF FLUID'//11X,1PD15.7,5X,'NET ZERO-ORDER PRODUCTION(+)/'&
     &  ,'LOSS(-) OF ENERGY IN FLUID'/11X,1PD15.7,5X,'NET ZERO-ORDER ', &
     &   'PRODUCTION(+)/LOSS(-) OF ENERGY IN SOLID GRAINS'              &
     &   /11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ENERGY THROUGH FLUID ',&
     &   'SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ',  &
     &   'ENERGY THROUGH INFLOWS OR OUTFLOWS AT POINTS OF SPECIFIED ',  &
     &   'PRESSURE'/11X,1PD15.7,5X,'NET GAIN(+)/LOSS(-) OF ENERGY ',    &
     &   'THROUGH ENERGY SOURCES AND SINKS'/11X,1PD15.7,5X,'NET GAIN(+)'&
     &  ,'/LOSS(-) OF ENERGY THROUGH POINTS OF SPECIFIED TEMPERATURE')  
!                                                                       
 1645    NSOPI = NSOP - 1 
         IF (NSOPI.EQ.0) GOTO 2000 
         IF (ME) 1649, 1649, 1659 
 1649    WRITE (fLST, 1650) 
 1650 FORMAT(///22X,'SOLUTE SOURCES OR SINKS AT FLUID SOURCES AND ',    &
     &   'SINKS'//22X,' NODE',8X,'SOURCE(+)/SINK(-)'/32X,               &
     &   '(SOLUTE MASS/SECOND)'/)                                       
! **************************************** XXH *******************************************																							 
!	 WRITE(20233,*) 'SOLUTE AT FLUID SOURCES ("NODE","Q (kg/s)"), IT=',IT
! **************************************** XXH *******************************************																						  
         GOTO 1680 
 1659    WRITE (fLST, 1660) 
 1660 FORMAT(///22X,'ENERGY SOURCES OR SINKS AT FLUID SOURCES AND ',    &
     &   'SINKS'//22X,' NODE',8X,'SOURCE(+)/SINK(-)'/37X,               &
     &   '(ENERGY/SECOND)'/)                                            
! **************************************** XXH *******************************************																							 
!	 WRITE(20233,*) 'ENERGY AT FLUID SOURCES ("NODE","Q (J/s)"), IT=',IT
! **************************************** XXH *******************************************
 1680    DO 1900 IQP = 1, NSOPI 
            I = IABS (IQSOP (IQP) ) 
            IF ( QIN(I) ) 1700, 1700, 1750 
 1700       QU = QIN(I) * CW * UVEC(I, K) 
            GOTO 1800 
 1750       QU = QIN(I) * CW * UIN(I, K) 
! 1800       WRITE (fLST, 450) I, QU 
 1800       CONTINUE
!           WRITE (20233, '(1E15.6e3)') QU   !!Yu													 
 1900    END DO 
!                                                                       
 2000    IF (NPBC.EQ.0) GOTO 2500 
         IF (ME) 2090, 2090, 2150 
 2090    WRITE (fLST, 2100) 
 2100 FORMAT(///22X,'SOLUTE SOURCES OR SINKS DUE TO FLUID INFLOWS OR ', &
     &   'OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'//22X,' NODE',8X,    &
     &   'SOURCE(+)/SINK(-)'/32X,'(SOLUTE MASS/SECOND)'/)               
! **************************************** XXH *******************************************																							 
!	      WRITE(20232,*) 'SOLUTE OF SPECIFIED PRESSURE ("NODE","Q kg/s"), IT=',IT
! **************************************** XXH *******************************************																						  
         GOTO 2190 
 2150    WRITE (fLST, 2160) 
 2160 FORMAT(///22X,'ENERGY SOURCES OR SINKS DUE TO FLUID INFLOWS OR ', &
     &   'OUTFLOWS AT POINTS OF SPECIFIED PRESSURE'//22X,' NODE',8X,    &
     &   'SOURCE(+)/SINK(-)'/37X,'(ENERGY/SECOND)'/)
 2190    DO 2400 IP = 1, NPBC 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            IF (QPLITR (IP) ) 2200, 2200, 2250 
 2200       QPU = QPLITR (IP) * CW * UVEC(I, K) 
 !           GOTO 2300 
 2250       QPU = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K) 
! 2300       WRITE (fLST, 450) I, QPU 
 2300    CONTINUE
!             WRITE (20232, '(1E15.6e3)')  QPU   !!Yu													  
 2400     END DO 
!
 2500    CONTINUE     
         IF (IBCT.EQ.4) GOTO 4500 
         NSOUI = NSOU (K) 
         INEGCT = 0 
         DO 3500 IQU = 1, NSOUI 
            I = IQSOU (IQU, K) 
            IF (I) 3400, 3500, 3500 
 3400       INEGCT = INEGCT + 1 
            IF (ME) 3450, 3450, 3460 
 3450       IF (INEGCT.EQ.1) WRITE (fLST, 3455) 
 3455 FORMAT(///22X,'TIME-DEPENDENT SOLUTE SOURCES AND SINKS'//22X,     &
     &   ' NODE',10X,'GAIN(+)/LOSS(-)'/30X,'  (SOLUTE MASS/SECOND)'//)  
! **************************************** XXH *******************************************																							 
!    WRITE(20234,*) 'SOLUTE OF TIME-DEPENDENT SOLUTE SOURCES ("NODE","Q kg/s"), IT=',IT
! **************************************** XXH *******************************************	
            GOTO 3475 
 3460       IF (INEGCT.EQ.1) WRITE (fLST, 3465) 
 3465 FORMAT(///22X,'TIME-DEPENDENT ENERGY SOURCES AND SINKS'//22X,     &
     &   ' NODE',10X,'GAIN(+)/LOSS(-)'/35X,'  (ENERGY/SECOND)'//)       
! **************************************** XXH *******************************************																							 
!    WRITE(20234,*) 'ENERGY OF TIME-DEPENDENT ENERGY SOURCES ("NODE","U J/s"), IT=',IT
! **************************************** XXH *******************************************
 3475       CONTINUE 
!            WRITE (fLST, 3490) - I, QUIN ( - I, K) 
!            WRITE (20234, '(1E15.6e3)')  QUIN ( - I, K)  !!YU
 3490 FORMAT(22X,I9,10X,1PD15.7) 
 3500    END DO 
!                                                                       
 4500    IF (NUBC (K) .EQ.0) GOTO 4888
         IF (K.EQ.NESP) GOTO 4610 
 4600    WRITE (fLST, 4650) 
 4650 FORMAT(///22X,'SOLUTE SOURCES OR SINKS DUE TO SPECIFIED ',        &
     &   'CONCENTRATIONS'//22X,' NODE',10X,'GAIN(+)/LOSS(-)'/30X,       &
     &   '  (SOLUTE MASS/SECOND)'/)                                     
! **************************************** XXH *******************************************																							 
 !   WRITE(20235,*) 'SOLUTE OF SPECIFIED CONCENTRATIONS ("NODE","Q kg/s"), IT=',IT
! **************************************** XXH *******************************************	
         GOTO 4690 
 4610    WRITE (fLST, 4660) 
 4660 FORMAT(///22X,'ENERGY SOURCES OR SINKS DUE TO SPECIFIED ',        &
     &   'TEMPERATURES'//22X,' NODE',10X,'GAIN(+)/LOSS(-)'/35X,         &
     &   '  (ENERGY/SECOND)'/)                                          
! **************************************** XXH *******************************************																							 
 !   WRITE(20235,*) 'ENERGY OF SPECIFIED TEMPERATURES ("NODE","Q J/s"), IT=',IT
! **************************************** XXH *******************************************
 4690    CONTINUE 
         DO 4700 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
!            WRITE (fLST, 450) I, QPLITR (IUP)
!            WRITE (20235, '(1E15.6e3)') QPLITR (IUP)   !!YU
4700     END DO 
!  
4888 CONTINUE
  IF (K.EQ.1)THEN         !!!!
     QFSTOT=QIUTOT+QPUTOT+QQUTOT+QULTOT 
   !!!!!!!!!!!!!!!!!!!!!!!!!!! temperature change 
     IF(IT.EQ.1)THEN
   WRITE(1102,9876)
    
9876 FORMAT(8X,'IT     FLDTOT         SLDTOT         FIUTOT         USCPTOT         USLPTOT         USDLTOT&
&     USLUTOT         FLPTOT         FLUTOT         QIUTOT         QPUTOT         QQUTOT         QULTOT         STETOT         QFSTOT ')
    ENDIF 
    WRITE(1102,'(1I9,20E15.6E3)')IT, FLDTOT,SLDTOT,FIUTOT,USCPTOT,USLPTOT,USDLTOT,USLUTOT,FLPTOT,FLUTOT,QIUTOT,QPUTOT,QQUTOT,QULTOT,STETOT,QFSTOT
    

    ENDIF         
         
         
    ENDIF  !!! IF(K.EQ.1)THEN          
   
    
    IF(K.NE.1)  THEN   !!!solute change
   
         FLDTOT = 0.D0 !1
         SLDTOT = 0.D0 
         DNSTOT = 0.D0 
         P1FTOT = 0.D0 
         P1STOT = 0.D0 
         P0FTOT = 0.D0 
         P0STOT = 0.D0 
         QQUTOT = 0.D0 
         QIUTOT = 0.D0 
         
         FIUTOT=0 !7
         USCPTOT=0 !2
         USLPTOT=0 !3
         USDLTOT=0
         USLUTOT=0 !4
         FLPTOT=0 !5
         FLUTOT=0 !6
         STETOT=0
         STETOTSL=0
         QFSTOT=0
    
         WATOT=0!总水质量
         SATOT=0!总盐质量
         SATOTSW=0
         SLTOT=0
         SITOT=0
         SIOUTTOT=0
         SMTOT=0
     IF (ME.EQ. - 1.AND.ADSMOD (K) .NE.'NONE      ') CALL ADSORB ()
!                                                                       
!.....SET APPROPRIATE VALUES FOR CS AND CW                              
         IF (K.EQ.NESP) THEN                 
            CS = CST 
            CW = CWT  
         ELSE                        
            CS = 0.0D0 
            CW = 1.0D0 
         ENDIF 
    
         DO 1301 I = 1, NN 
            imap = NodeMap(i)
            ESRV = NodeData(imap)%por * SLL(I) * RHO(I) * VOL(I) 
            EPRSV = (1.D0 - NodeData(imap)%por ) * NodeData(imap)%rhos * VOL(I) 
            EICE = -HTLAT + CI*UVEC(I,1)                                           !
            DUDT = (1 - ISSTRA) * ( UVEC(I, K) - UM1(I, K) ) / DELTU 
            FLDTOT = FLDTOT + ESRV * CW * DUDT !1
         !   SLDTOT = SLDTOT + EPRSV * CS1(I, K) * DUDT          
!======= FLDTOT + SLDTOT = total rate of change in stored solute mass in the region due to change in concentration
!======= FLDTOT + SLDTOT equals to (5.15a) on page 112 of the manual		
 !!!!新增,最后算的是总存储量的变化
 FIUTOT=FIUTOT+ NodeData(imap)%por * SII(I) * RHOI * VOL(I) *CW* DUDT  !7        
USCPTOT=USCPTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,K)*RHO(I)*SLL(I)*NodeData(NodeMap(i))%sopi*(PVEC(I)-PM1(I))/DELTP !2
USLPTOT=USLPTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,K)*RHO(I)*NodeData(imap)%por*DSLDP(I)*(PVEC(I)-PM1(I))/DELTP !3
!USDLTOT=USDLTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,K)*NodeData(imap)%por*SLL(I)*DRWDU(1)*(UVEC(I,1)-UM1(I,1))/DELTU
USLUTOT=USLUTOT+ (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,K)*NodeData(imap)%por*RHO(I)*DSLDT(I)*(1-ISSTRA)*(UVEC(I,K)-UM1(I,K))/DELTU !4
FLPTOT= FLPTOT+  (1-ISSFLO/2)*VOL(I)*CW*UVEC(I,K)*NodeData(imap)%por*RHO(I)*DSIDP(I)*(PVEC(I)-PM1(I))/DELTP!5
FLUTOT =FLUTOT+ VOL(I)*CW*UVEC(I,K)*NodeData(imap)%por*1000*DSIDT(I)*(1-ISSTRA)*(UVEC(I,K)-UM1(I,K))/DELTU !6
STETOT =FLDTOT+FIUTOT+ USCPTOT+USLPTOT+USLUTOT+ FLPTOT+FLUTOT
STETOTSL=FLDTOT+USCPTOT+USLPTOT+USLUTOT
!!!总水总盐质量
SLTOT=SLTOT+VOL(I)*RHO(I)*NodeData(imap)%por*SLL(I)
SITOT=SITOT+VOL(I)*RHOI*NodeData(imap)%por*SII(I)
WATOT=SLTOT+SITOT
SATOT=SATOT+VOL(I)*RHO(I)*NodeData(imap)%por*SLL(I)*UVEC(I,2)
SATOTSW=SATOTSW+VOL(I)*RHO(I)*NodeData(imap)%por*SW(I)*UVEC(I,2)
SMTOT=SMTOT+SM(I)  
SIOUTTOT=SIOUTTOT+SIOUT(I)
!IF(I==299)WRITE(*,*)'RHO=',RHO(299)
            DNSTOT = DNSTOT + CW * UVEC(I, K) * (1 - ISSFLO / 2) *      &
              VOL(I) * ( RHO(I) * (SLL(I) * NodeData(imap)%sop + NodeData(NodeMap(i))%por * DSLDP(I) ) * &
              ( PVEC(I) - PM1(I) ) / DELTP + NodeData(imap)%por * SLL(I) *  &
              DRWDU (K) * ( UM1(I, K) - UM2(I, K) ) / DLTUM1) 
         
                   QQUTOT = QQUTOT + QUIN(I, K) 
!======= QQUTOT is the diffusive-dispersive source of solute mass, (5.20) on page 113
            IF ( QIN(I) ) 1201, 1201, 1251 
 1201       QIUTOT = QIUTOT + QIN(I) * CW * UVEC(I, K) 
            GOTO 1301 
 1251       QIUTOT = QIUTOT + QIN(I) * CW * UIN(I, K)
!======= QIUTOT is the solute mass change at fluid source node, (5.18) on page 113 
   1301     END DO   
         
   QPUTOT = 0.      
       DO 1501 IP = 1, NPBC 
          IF (QPLITR(IP).LE.0) THEN 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * UVEC(I, K)
            QPSOLUTE = QPLITR (IP) * CW * UVEC(I, K)
          ELSEIF (QPLITR(IP).GT.0) THEN 
            QPUTOT = QPUTOT + QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
            QPSOLUTE = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K)
          ENDIF
!            WRITE(4,'(2I8,4E15.7)') K, IABS(SpecifiedPBC(IP)%node), QPSOLUTE   !same to 20235
    1501 END DO        
         
     QULTOT = 0.D0 
         IF (NUBC (K) .EQ.0) GOTO 1541 
         DO 1511 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
            QPLITR (IUP) = GNUU (K) * ( MultiSpeciesBC(K)%SpecifiedU(IU)%U - UVEC(I, K) ) 
            QULTOT = QULTOT + QPLITR (IUP) 
!======= QULTOT is the solute mass due to specified concentration conditions, (5.21) on page 114
1511     END DO     
  
         
   !!!!底下的没在算通量了      
1541          GOTO 1646 
         
  1646    NSOPI = NSOP - 1 
         IF (NSOPI.EQ.0) GOTO 2001 
   DO 1901 IQP = 1, NSOPI 
            I = IABS (IQSOP (IQP) ) 
            IF ( QIN(I) ) 1701, 1701, 1751 
 1701       QU = QIN(I) * CW * UVEC(I, K) 
            GOTO 1801 
 1751       QU = QIN(I) * CW * UIN(I, K) 
! 1800       WRITE (fLST, 450) I, QU 
 1801       CONTINUE
!           WRITE (20233, '(1E15.6e3)') QU   !!Yu													 
 1901    END DO 
!                                                                       
 2001    IF (NPBC.EQ.0) GOTO 2501 


           DO 2401 IP = 1, NPBC 
            I = IABS ( SpecifiedPBC(IP)%node ) 
            IF (QPLITR (IP) ) 2201, 2201, 2251 
 2201       QPU = QPLITR (IP) * CW * UVEC(I, K) 
 !           GOTO 2300 
 2251       QPU = QPLITR (IP) * CW * SpecifiedPBC(IP)%U(K) 

!             WRITE (20232, '(1E15.6e3)')  QPU   !!Yu													  
 2401     END DO 
!
 2501    CONTINUE     
         IF (IBCT.EQ.4) GOTO 4501 
         NSOUI = NSOU (K) 
         INEGCT = 0 
         DO 3501 IQU = 1, NSOUI 
            I = IQSOU (IQU, K) 
            IF (I) 3401, 3501, 3501 
 3401       INEGCT = INEGCT + 1 

!            WRITE (20234, '(1E15.6e3)')  QUIN ( - I, K)  !!YU
 3491 FORMAT(22X,I9,10X,1PD15.7) 
 3501    END DO 
!                                                                       
 4501    IF (NUBC (K) .EQ.0) GOTO 4889
         
         DO 4701 IU = 1, NUBC (K) 
            IUP = IU + NPBC 
            I = IABS ( MultiSpeciesBC(K)%SpecifiedU(IU)%node ) 
!            WRITE (fLST, 450) I, QPLITR (IUP)
!            WRITE (20235, '(1E15.6e3)') QPLITR (IUP)   !!YU
4701     END DO 
!  
4889 CONTINUE         
         

     QFSTOT=QIUTOT+QPUTOT+QQUTOT+QULTOT 
   !!!!!!!!!!!!!!!!!!!!!!!!!!! solute change  
     IF(IT.EQ.1)THEN
   WRITE(1108,98762)
    
98762 FORMAT(8X,'IT     FLDTOT         USCPTOT         USLPTOT  &
&     USLUTOT         FLPTOT         FLUTOT         FIUTOT         QIUTOT         QPUTOT         QQUTOT         QULTOT         STETOT         QFSTOT         STETOTSL ')
    ENDIF 
    WRITE(1108,'(1I9,20E15.6E3)')IT, FLDTOT,USCPTOT,USLPTOT,USLUTOT,FLPTOT,FLUTOT,FIUTOT,QIUTOT,QPUTOT,QQUTOT,QULTOT,STETOT,QFSTOT,STETOTSL
    

  IF(IT.EQ.1)THEN
   WRITE(1109,98764)
    
98764 FORMAT(8X,'IT     SLTOT         SITOT         WATOT         SATOT         SMTOT         SM+SA')
    ENDIF  
   WRITE(1109,'(1I9,20E15.6E3)')IT, SLTOT,SITOT,WATOT,SATOT,SMTOT ,SMTOT+SATOT ,SATOTSW,SIOUTTOT   
         
         
         
         
         
    
    ENDIF  !!! IF(K.EQ.2)THEN 
    
 5000 END DO 
!.....END OF SPECIES LOOP                                               
!                                                                       
5500 CONTINUE 
 
                                                                  
!.....RESET VALUES FOR ENERGY TRANSPORT FOR SUBSEQUENT CALCULATIONS     
      CS = CST 
      CW = CWT 
!                                                                             
      RETURN 
      END SUBROUTINE BUDGET                         
