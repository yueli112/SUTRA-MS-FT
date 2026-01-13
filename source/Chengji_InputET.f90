SUBROUTINE INPUTET ()
! THE SUBROUTINE IS USED TO READ THE "ET.INP" FILE
! DATA CONTAINED IN THE FILE INCLUDE:
! 1 -- TIDAL PARAMETERS
! 2 -- PERTURBATION OF PERMEABILITY FIELD
! 3 -- SOIL RETENTION CURVE PARAMETERS

  USE M_SWCC
  USE M_TIDE
  USE M_CONTROL
  USE M_SEEPAGE
  USE ALLARR                !
  USE FRPARAMS               !
  USE SALTPARAMS            !
  USE FRCONTRL
  USE M_OUT
  USE M_SFCC
  IMPLICIT NONE
  INTEGER(4) :: fET
  INTEGER(4) :: NLSKIP
  INTEGER(4):: IOS
  
 
  fET=1221
  OPEN(UNIT=fET,FILE='ET.INP',STATUS='OLD')
  
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) TASP,TANE,TPSP,TPNP,TM,RHOST,SC,BNDSWT,FWT,FWH
  if(ios/=0) call ErrorIO('Error reading tide parameters')

  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) SWRES1,AA1,VN1,SWRES2,AA2,VN2,SWRES3,AA3,VN3
  if(ios/=0) call ErrorIO('Error reading SWCC parameters')
  
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) AVGVEL,NSTEP,TLAB,HOMO,AMP
  if(ios/=0) call ErrorIO('Error reading CONTROL parameters')
  
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) SEEP,XSTART,ZSTART
  if(ios/=0) call ErrorIO('Error reading SEEPAGE parameters')

  !!!*************************************************************************freeze
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) IFREEZ,COMPI,CI,SIGMAI,RHOI,HTLAT                                         !6
  if(ios/=0) call ErrorIO('Error reading ice property parameters')
  
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios)FIRS,SPEFR, SLSATRES1,TLRES1,TFREEZ1,RKMIN1,SLSATRES2,TLRES2,TFREEZ2,&            !14
                       RKMIN2,SLSATRES3,TLRES3,TFREEZ3,RKMIN3
  if(ios/=0) call ErrorIO('Error reading liquid saturation parameters')

    
  !!!************************************************************************Salt precipitation
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) SOLUBILITY ,KSALT                                       
  if(ios/=0) call ErrorIO('Error reading salt parameters')
  
  
  
    CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) COLYU, S_OUT, OUT_INTE, ITtot
  if(ios/=0) call ErrorIO('Error reading OUT parameters')  
  
  
  CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) IFEVAP,IFSUCTION,BOUNDSPE,DM1,DM2,DEV
  if(ios/=0) call ErrorIO('Error reading evaporation and crysuction  parameters')  
  
   CALL SKPCOM (fET, NLSKIP)
  READ (fET,*,iostat=ios) AASFCC1,KSF
  if(ios/=0) call ErrorIO('Error reading SFCC  parameters')  
  
  
RETURN

END SUBROUTINE