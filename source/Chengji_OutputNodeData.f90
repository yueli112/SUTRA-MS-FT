!     SUBROUTINE        O  U  T  N  O  D  E    SUTRA-MS VERSION 2004.1
!                                                                       
! *** PURPOSE :                                                         
! ***  TO PRINT NODE COORDINATES, PRESSURES OR HEADS, CONCENTRATIONS OR 
! ***  TEMPERATURES, AND SATURATIONS IN A FLEXIBLE, COLUMNWISE FORMAT.  
! ***  OUTPUT IS TO UNIT fNOD.                                            
!                                                                       
      SUBROUTINE OUTNODE2(PVEC, UVEC, SW, SLL,SII,IN, X, Y, Z, TITLE1, TITLE2) 
      USE ITERAT 
      USE CONTRL 
      USE SOLVI 
      USE JCOLS 
      USE PARAMS
      
      USE FRPARAMS         !!
      USE FRCONTRL
      
      USE FUNITS 
      USE DIMS
      USE DIMX
      USE TIMES
      USE KPRINT
      USE GRAVEC
      USE PLT1
      USE SutraMSPrecision
      USE M_CONTROL
      IMPLICIT NONE
      real (DP) :: &
        PVEC (NN), UVEC (NN, NSPE), SW (NN) ,SLL(NN), SII(NN)
      real (DP) :: &
        X (NN), Y (NN), Z (NN) 
      real (DP) :: &
        VCOL (NCOLMX), VVAR (9 + NSPE-1)               !!VVAR∂‡¡À2
      CHARACTER(1) TITLE1 (80), TITLE2 (80) 
      CHARACTER(8) HORP 
      CHARACTER(13) TORC (NSPE) 
      CHARACTER(1) CPHORP, CPTORC, CPSATU 
      integer (I4B) :: &
        IN (NIN), IIN (8)
      real (I4B) :: &
        TT (99999)
      integer (I4B) :: &
        ITT (99999), ISTORC (99999), ISHORP (99999), &
        ISSATU (99999) 
      LOGICAL PRINTN 
      !LOCAL VARIABLES
      CHARACTER (LEN=15) &
        TCHAR
      INTEGER (I4B) :: &
        LCHAR
      INTEGER (I4B) :: &
        TS, &
        LCHORP, LCTORC, &
        I, &
        JT, JTMAX, &
        K, KK, KT, KTMAX, &
        M
      REAL (DP) :: &
        DKTM2, DELTK
      
      INTEGER (I4B) :: A
      !LOCAL VARIABLES Chengji 2015-08-10
      !REAL (DP):: XNT(NN),YNT(NN),ZNT(NN)
      !REAL (DP):: MEANP(NN),MEANT(NN),MEANC(NN),MEANS(NN)    
!.....Calculate headers on time step 1                                  
!.....and create output on time steps greater than or equal to 1      
 !     IF (IT.EQ.1) THEN
 !          XNT(I)=0.
 !          YNT(I)=0.
 !          ZNT(I)=0.
 !          MEANP(I)=0.
	!       MEANT(I)=0.
 !          MEANC(I)=0.
	!       MEANS(I)=0.
 !      ENDIF 
 !       IF(MOD(IT,NSTEP).EQ.1) THEN          !    print interval
 !         
 !      DO 35760 I=1,NN
 !          XNT(I)=0.
 !          YNT(I)=0.
 !          ZNT(I)=0.
 !          MEANP(I)=0.
	!       MEANT(I)=0.
 !          MEANC(I)=0.
	!       MEANS(I)=0.
 !35760      CONTINUE
 !       Endif                      

!.....  The nodewise data for this time step  
      DO 980 I = 1, NN    
            VVAR (1) = DBLE (I) 
            VVAR (2) = X (I) 
            VVAR (3) = Y (I) 
            VVAR (4) = Z (I) 
            VVAR (5) = PVEC (I) 
            KK = 0 
            DO K = 6, (5 + NSPE) 
            KK = KK + 1 
            VVAR (K) = UVEC (I, KK) 
            ENDDO 
            VVAR (6 + NSPE) = SW (I)
            VVAR (7 + NSPE) = SLL (I)                !!
            VVAR (8 + NSPE) = SII (I)
            DO 972 M = 1, NCOLS5 
               VCOL (M) = VVAR (J5COL (M) ) 
  972       END DO 

 !! ******** SUM UP THE VELOCITY OVER ONE TIDAL CYCLE ********
 !           XNT(I) = XNT(I) + X (I)
 !           YNT(I) = YNT(I) + Y (I)
 !           ZNT(I) = ZNT(I) + Z (I)
 !           MEANP(I)= MEANP(I)+PVEC(I)
 !           MEANT(I)= MEANT(I)+UVEC(I,1)
	!		MEANS(I)= MEANS(I)+SW(I)
 !           MEANC(I)= MEANC(I)+UVEC(I,2)
 !! ******** SUM UP THE VELOCITY OVER ONE TIDAL CYCLE ******** 
  980    END DO
!        if (mod(IT,NSTEP)  .EQ. 0) then           !    steps over a tidal cycle
!	  WRITE (200209,*) 'VARIABLES = "X","Y", "Z", "P", "T","C", "S"'
!      WRITE (200209,*)' Time step =',IT

        	!DO 89689 I=1,NN
!	         write (200209,'(7E15.7)') XNT(I)/NSTEP,YNT(I)/NSTEP,ZNT(I)/NSTEP,MEANP(I)/NSTEP, MEANT(I)/NSTEP,MEANC(I)/NSTEP,MEANS(I)/NSTEP
!             write (200209,'(7E15.7)') XNT(I)/4,YNT(I)/4,ZNT(I)/4,MEANP(I)/4, MEANT(I)/4,MEANC(I)/4,MEANS(I)/4
!89689       CONTINUE
!        Endif 
        
      RETURN 
!                                                                       
      END SUBROUTINE OUTNODE2
   
         
                         