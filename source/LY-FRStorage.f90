module FRStorage                                    
 USE DIMS    
 USE DIMX   
 USE MSErrorHandler
       use SutraMSPrecision 
       implicit none   
    
  !  real (DP), allocatable :: & 
  !  SLL(:),  DSLDP(:),DSLDT(:),SII(:),  DSIDP(:),DSIDT(:)
  !  
  !REAL (DP), ALLOCATABLE :: &  
  !  SLL_(:),SII_(:)
  !  
  ! public & 
  !SLL,  DSLDP,DSLDT,SII,  DSIDP,DSIDT,&  
  ! SLL_,SII_ 
  !  
    end module  FRStorage
 
    
   MODULE ALLARR               !Related to the regional parameters                                  
   IMPLICIT NONE
   DOUBLE PRECISION :: FIRS,SPEFR,SLSATRES1,TLRES1,TFREEZ1,SLSATRES2,TLRES2,TFREEZ2,&
                       SLSATRES3,TLRES3,TFREEZ3,RKMIN1,RKMIN2,RKMIN3
 
    
    
    
    END MODULE ALLARR
 
    
    MODULE FRPARAMS 
        USE SutraMSPrecision
     REAL (DP) ::CI,COMPI,RHOI   ,&
                       HTLAT,SIGMAI 
      
         
     PUBLIC COMPI,CI,SIGMAI, RHOI, HTLAT                    !
      END MODULE FRPARAMS    
        
        
        
  MODULE FRCONTRL 
        USE SutraMSPrecision       
       INTEGER (I4B) :: IALSAT,IFREEZ
       PUBLIC  IALSAT,IFREEZ
      END MODULE FRCONTRL    
        
        
        
        
  MODULE SALTPARAMS         !!
     USE SutraMSPrecision   
     REAL (DP) ::SOLUBILITY,KSALT
     real (DP), allocatable ::TOTWA(:),SM(:),TOTSA(:)
      PUBLIC SOLUBILITY,KSALT,TOTWA,SM,TOTSA
  END MODULE SALTPARAMS     