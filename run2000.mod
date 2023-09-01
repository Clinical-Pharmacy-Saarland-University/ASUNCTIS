;; 1. Based on: run2000
;; 2. Description: PK Model+Friberg model, Circ0(Platelets,LDH,BL_Lymphos,IMMUNO), QC Dataset, -AE, PAPER

;; x1. Author: dings


;; 3. Label:

$PROBLEM APOGENIX
$INPUT ID TIME Screening DROP Randomization Site DROP CMT AMT DV MDV EVID RATE FLAGDUR FLAGTIME OBS ALT_U_L AST_U_L Bilirubin_mg_dL Creatinine_mg_dL CRP_mg_L DIME_D_ng_mL Glucose_mg_dL IL_6_corrected_pg_mL IL_6_pg_mL LDH_U_L Lymphocytes_10_9_L Platelets_10_9_L PROCAL_corrected_ng_mL PROCAL_ng_mL Uric_acid_umol_L WBC_10_9_L AGE WEIGHT HEIGHT SEX DOSEMG CORTICO IMMUNO DROP FLAGTREAT FLAGSTATUS TIMEsinceHosp BMI BSA eGFR BL_Lymphocytes BL_CRP BL_WHO CORTICO_YN IMMUNO_YN FLAG
$DATA     ..\DATASET\APG101_NONMEM_DATASET_V04_QC.csv IGNORE=@ IGNORE(CMT.GE.3) IGNORE(FLAG.EQ.1)

$SUBROUTINES ADVAN13 TOL=5

$MODEL
COMP(CENT)          ;central compartment        CMT 1
COMP(Peri1)         ;Peripheral                 CMT 2
COMP(Peri2)         ;Peripheral 2               CMT 3
COMP(STEM)			; Stem Cells		    	CMT 4
COMP(TRN1)			; Transit 1			        CMT 5
COMP(TRN2)			; Transit 2			        CMT 6
COMP(TRN3)			; Transit 3			        CMT 7
COMP(CIRC)	        ; Circ Cells	            CMT 8

$PK
;-------------------PK----------------------------------
FACT = 343 ;µmol/L
HILL = 31.6
MAX = (Uric_acid_umol_L**HILL/(Uric_acid_umol_L**HILL+FACT**HILL))
V1 = 3.94 *(WEIGHT/79.6)**(0.704)*EXP(ETA(2))
RATIO = 0.62
SIALICEC50 = 0.378
SIALICHILL = 8.03
CLMAX = 0.0147
MAX2 = CLMAX-(CLMAX*RATIO**SIALICHILL/(RATIO**SIALICHILL+SIALICEC50**SIALICHILL))
CL = (0.00773 * (1-MAX) *EXP(ETA(1)) + (MAX2))*24;L/h -> L/day
Q2 = 0.199*24 ;L/h -> L/day
V2 = 6.96 *EXP(ETA(3))  ;L ;*(PROTEIN/70)**FACT6
Q3 = 0.00985 *(Bilirubin_mg_dL/8.6*0.0585)**0.472*EXP(ETA(4))*24 ;L/h -> L/day; Bilirubin µmol/L -> mg/dl *(ALBUMIN/43)**(-1.89) 
V3 = 106 ;L
K10 = CL/V1
K12 = Q2/V1
K21 = Q2/V2
K13 = Q3/V1
K31 = Q3/V3
S1 = V1/1000

;-------------------PD----------------------------------
Circ0 = THETA(1)*EXP(ETA(5))*(Platelets_10_9_L/259)**THETA(7)*(LDH_U_L/323)**THETA(8)*(BL_Lymphocytes/0.95)**THETA(9)*(1+IMMUNO*THETA(10))
Baseline = BL_Lymphocytes
A_0(4) = Baseline
A_0(5) = Baseline
A_0(6) = Baseline
A_0(7) = Baseline
A_0(8) = Baseline

MTT = THETA(2)*EXP(ETA(6))

KTR = 4/MTT
KSYN = KTR
KCIRC = KTR

GA = THETA(3)

Emax = THETA(4)
EC50 = THETA(5)

CRP = CRP_mg_L
IF(CRP_mg_L.LT.10) CRP = 10
CRP_Effect = (CRP/10)**THETA(6)

$DES
;-------------------PK----------------------------------
DADT(1) = K21*A(2) - K12*A(1) + K31*A(3) - K13*A(1) - K10*A(1)
DADT(2) = -K21*A(2) + K12*A(1)
DADT(3) = -K31*A(3) + K13*A(1)
CONC = A(1)/V1
;-------------------PD----------------------------------
EFFECT = 1+Emax*CONC/(CONC+EC50)
RBD = (Circ0/A(8))**(GA*EFFECT)
DADT(4) = KSYN*A(4)*RBD - KTR*A(4)
DADT(5) = KTR*A(4) - KTR*A(5) 
DADT(6) = KTR*A(5) - KTR*A(6) 
DADT(7) = KTR*A(6) - KTR*A(7) 
DADT(8) = KTR*A(7) - KCIRC*A(8)*CRP_Effect

$ERROR
IPRED = A(1)/V1
IF(CMT.EQ.2) IPRED = A(8)
W = IPRED
DEL = 0
IF(IPRED.EQ.0) DEL = 0.001
IRES	= DV -  IPRED
IWRES   = IRES/(W+DEL)
Y  = IPRED *(1+EPS(1)) + EPS(2)

$THETA
(0, 1) ;1 Circ0 (10**9/L)
(0.1, 1,10) ;2 MTT (days)
(0, 0.02,1) ;3 GA
(-1, 5) ;4 emax
(0.5, 2) ;5 EC50
(-10, 0.1) ;6 CRP-effect
(-10, 0.5, 10);7 Circ0(Platelets)
(-10, -0.5, 10);8 Circ0(LDH)
(-10, 0.1, 19);9 Circ0(BL_Lymphos)
(-1, -0.5, 19);10 Circ0(IMMUNO)

$OMEGA  
 0 FIX  ;0.764 FIX    
 0 FIX  ;0.0368 FIX    
 0 FIX  ;0.279 FIX    
 0 FIX  ;0.638 FIX    
 0.1 ; IIV_Circ0
 0.1 ; IIV_MTT

$SIGMA
 0.1 ; PD PE
 0 FIX ; PD AE

$EST METHOD=1 INTER MAXEVAl=99999 NOABORT PRINT=1 SIG=3 POSTHOC
$COV
$TABLE ID TIME ETAS(1:LAST) CMT MDV EVID IPRED IWRES CWRES NOPRINT ONEHEADER FILE=sdtab2000
