# 2nd edition of constants (included NFkB noncanonical pathway) w/ simplified Apoptosis module
# Define global constants for model
#---------------------------------------------
const ANTIGEN_DOSE = 0.03;   # 0.03 for high dose of BCR, 0.01nM and 0.003nM for medium & low dose of a-BCR
const CD40L_DOSE = 100;   # 100nM for high dose and 10nM, 3nM for medium and low dose
const CD40L_DELAY = 0;   # 0, 1, 3, 5, or 8hr delay


const CONVERSION = 1;   # unit conversion (default in hr-1). Set to 1/60 or 1/3600 to convert to min-1 or sec-1.
const APOP_CONVERSION = 0.0078;   # apoptosis #molecule to nM conversion: 10e9 * molecule / [(6.02 * 10e23) * (212 * 10e-18 * 1000)] = 0.78 * 10e-2 = 0.0078

const EPS = 1;    # translation modifier
const IKK_MOD = 1.5;    # IKK modifer (scale IKK from 0-100 AU to 0-150 nM)
# const BASAL_IKK = 0.01;
# const TOTAL_IKK = 170;    # unit: nM (170 in RoyMitchell 2019, 288.48 in Mitchell 2018, 140 in Mitchell 2022 review)
# const BASAL_NIK = 0.1;    # ~ 5nM at baseline in Mitchell 2022 review)
# const TOTAL_NIK = 50;    # unit: nM (~50 at 48h in Mitchell 2022 review)
const NIK_KM = 10;     # Km for NIK-mediated processing of IkBd and p100

const BASAL_TBCL2 = 277.2133 * 60 * CONVERSION;    # basal Bcl2t transcription, unit: molecules / hr
# const cRelFeedback = 10;    # cRel induction modifier
# const constituitiveMultiplier = (1+(6/150))/(1+cRelFeedback*(6/150));
const P100_MOD = 9e-4;    # default: 9e-4

# thresholds for tMyC, tBcl2, tCycD transcription
const MYCTHR = 40;     # unit: nM (40 in Mitchell 2018)
const BCL2THR = 60;
const CYCDTHR = 40;
const GROWTHR = 40;

const SCALE_CELLULAR2MEDIA = 0.001;   # scale for external ligands like CD40L, a-IgM, LPS, etc.
const SCALE_MITO2CELLULAR = 1/0.07;   # scale for mitochondrial species


# Parameter distributions constants
#---------------------------------------------
const FOUNDER_RATE_CV = 0.112; # a CV of 11.2% in rate constants result in 32% CV in concentrations in steady state
const FOUNDER_CONC_CV = 0.00; # initial concentrations don't get distributed directly -- they are distributed due to distribution of rate constants after the cell reach steady states

const DAUGHTER_RATE_CV = 0.00;
const DAUGHTER_CONC_CV = 0.00;
const DAUGHTER_PARTITION_CV = 0.072; # how asymmetric each division needs to be

const FOUNDER_CELL_NUM = 8; # 125 in Mitchell et al. 2018
const MAX_GEN = 7; # maximum generation
const BURN_IN_PERIOD = 240; # hrs, pre-simulation phase
const GLOBAL_END_TIME = 140; # hrs



# Set up the ODE parameters
#---------------------------------------------
const method = Tsit5(); # Integration algorithm
const abserr = 1.0e-5; # Absolute error (og value: 1.0e-8)
const relerr = 1.0e-3; # Relative error (og value: 1.0e-6)
const saveat = 0.01; # deprecated by save_everystep = false, but in case we want to save timepoints for plots



# Define indices for rates & reactionFlux
#---------------------------------------------
const TOTAL_PARAMS = 16;   # number of params for each species

# (1) Initial concentration
const INITIALCONC = 1;
# (2) Scaling factor (for unit conversion)
const SCALE = 2;
# (3) Basal synethesis rate
const BASALSYNTHESIS = 3;
# (4) Induced synthesis rate
const SYNTHESIS = 4;
# (5) Translation rate
const TRANSLATION = 5;
# (6) Basal degradation rate
const BASALDECAY = 6;
# (7) Induced degradation (or repression) rate
const DECAY = 7;
# (7) Repression (irrelevant to rate, only to reactionFlux)
const REPRESSION = 8;
# (8) Complex association rate
const ASSOCIATION = 9;
# (9) Complex dissociation rate
const DISSOCIATION = 10;
# (10) Catalytic / induced phosphorylation rate
const PHOSPHORYLATION = 11;
# (11) Catalytic / induced dephosphorylation rate
const DEPHOSPHORYLATION = 12;
# (10) Catalytic / induced activation rate
const ACTIVATION = 13;
# (11) Catalytic / induced deactivation rate
const DEACTIVATION = 14;
# (12) Transportation rate into Mitochondria / Nucleus
const TRANSPORTIN = 15;
# (13) Transportation rate out of Mitochondria / Nucleus
const TRANSPORTOUT = 16;



# Define number of species for each module
#---------------------------------------------
const RECEPTOR_SPECIES = 20-2; #-2 because NFKB_SPECIES start from 3
const NFKB_SPECIES = 106;
const DIFF_SPECIES = 4;
const APOPTOSIS_SPECIES = 12; # 59 in original apoptosis module
const PROLIF_SPECIES = 25;
const TOTAL_SPECIES = RECEPTOR_SPECIES + NFKB_SPECIES + DIFF_SPECIES + APOPTOSIS_SPECIES + PROLIF_SPECIES + 1;


# Define indices for Receptor species
#-------------------------------------------------------------------------------------------------------------------
# 1 : ANTIGEN
const ANTIGEN = 1;
# 2 : BCR
const BCR = 2;
# 3 : Antigen-BCR
const ABCR = 3;
# 4 : inactivated CBM
const CBM = 4;
# 5 : activated CBM-complex
const ACBM = 5;
# 6 : inhibited CBM
const ICBM = 6;

# 7 : CD40L
const CD40L = 7;
# 8 : CD40R
const CD40R = 8;
# 9 : CD40L-R complex
const CD40LR = 9;
# 10 : TRAF6_off
const TRAF6_OFF = 10;
# 11 : TRAF6
const TRAF6 = 11;
# 12 : TRAF3-cIAP1/2 complex
const TRAF3 = 12;

# 13 : TAK1
const TAK1 = 13;
# 14 : activated TAK1
const ATAK1 = 14;
# 15 : IKK (IKK1)
const IKK_OFF = 15;
# 16 : pIKK (IKK2)
const IKK2 = 16;
# 17 : ppIKK (IKK3)
const IKK3 = 17;
# 18 : pppIKK (IKK4)
const IIKK = 18;
# 19 : total activated IKK (IKK2+IKK3)
const IKK = 19;
# 20 : NIK
const NIK = 20;

# Define indices for NFkB species
#-------------------------------------------------------------------------------------------------------------------
# 3 : tIkBa
const TIKBA = 3+RECEPTOR_SPECIES;
# 4 : tIkBb
const TIKBB = 4+RECEPTOR_SPECIES;
# 5 : tIkBe
const TIKBE = 5+RECEPTOR_SPECIES;
# 6 : IkBa
const IKBA = 6+RECEPTOR_SPECIES;
# 7 : IkBb
const IKBB = 7+RECEPTOR_SPECIES;
# 8 : IkBe
const IKBE = 8+RECEPTOR_SPECIES;
# 9 : IkBd
const IKBD = 9+RECEPTOR_SPECIES;
# 10 : IkBan
const NIKBA = 10+RECEPTOR_SPECIES;
# 11 : IkBbn
const NIKBB = 11+RECEPTOR_SPECIES;
# 12 : IkBen
const NIKBE = 12+RECEPTOR_SPECIES;
# 13 : IkBdn
const NIKBD = 13+RECEPTOR_SPECIES;
# 14 : tRelA
const TRELA = 14+RECEPTOR_SPECIES;
# 15 : tp50
const TP50 = 15+RECEPTOR_SPECIES;
# 16 : tRelB
const TRELB = 16+RECEPTOR_SPECIES;
# 17 : tP100
const TP100 = 17+RECEPTOR_SPECIES;
# 18 : tcRel
const TCREL = 18+RECEPTOR_SPECIES;
# 19 : RelA
const RELA = 19+RECEPTOR_SPECIES;
# 20 : p50
const P50 = 20+RECEPTOR_SPECIES;
# 21 : RelB
const RELB = 21+RECEPTOR_SPECIES;
# 22 : P100
const P100 = 22+RECEPTOR_SPECIES;
# 23 : cRel
const CREL = 23+RECEPTOR_SPECIES;
# 24 : P52
const P52 = 24+RECEPTOR_SPECIES;
# 25 : RelAn
const NRELA = 25+RECEPTOR_SPECIES;
# 26 : p50n
const NP50 = 26+RECEPTOR_SPECIES;
# 27 : nRelB
const NRELB = 27+RECEPTOR_SPECIES;
# 28 : nP100
const NP100 = 28+RECEPTOR_SPECIES;
# 29 : cReln
const NCREL = 29+RECEPTOR_SPECIES;
# 30 : nP52
const NP52 = 30+RECEPTOR_SPECIES;
# 31 : RelA:RelA
const AA = 31+RECEPTOR_SPECIES;
# 32 : RelA:p50
const A50 = 32+RECEPTOR_SPECIES;
# 33 : RelA:p52
const A52 = 33+RECEPTOR_SPECIES;
# 34 : RelB:p50
const B50 = 34+RECEPTOR_SPECIES;
# 35 : RelB:p52
const B52 = 35+RECEPTOR_SPECIES;
# 36 : cRel:p50
const C50 = 36+RECEPTOR_SPECIES;
# 37 : cRel:p52
const C52 = 37+RECEPTOR_SPECIES;
# 38 : cRel:p100
const C100 = 38+RECEPTOR_SPECIES;
# 39 : p50:p50
const P50P50 = 39+RECEPTOR_SPECIES;
# 40 : p52:p52
const P52P52 = 40+RECEPTOR_SPECIES;
# 41 : RelA:RelAn
const NAA = 41+RECEPTOR_SPECIES;
# 42 : RelA:p50n
const NA50 = 42+RECEPTOR_SPECIES;
# 43 : RelA:p52n
const NA52 = 43+RECEPTOR_SPECIES;
# 44 : RelB:p50n
const NB50 = 44+RECEPTOR_SPECIES;
# 45 : RelB:p52n
const NB52 = 45+RECEPTOR_SPECIES;
# 46 : cRel:p50n
const NC50 = 46+RECEPTOR_SPECIES;
# 47 : cRel:p52n
const NC52 = 47+RECEPTOR_SPECIES;
# 48 : cRel:p100n
const NC100 = 48+RECEPTOR_SPECIES;
# 49 : p50:p50n
const NP50P50 = 49+RECEPTOR_SPECIES;
# 50 : p52:p52n
const NP52P52 = 50+RECEPTOR_SPECIES;
# 51 : IkBa-RelA:RelA
const IKBAAA = 51+RECEPTOR_SPECIES;
# 52 : IkBb-RelA:RelA
const IKBBAA = 52+RECEPTOR_SPECIES;
# 53 : IkBe-RelA:RelA
const IKBEAA = 53+RECEPTOR_SPECIES;
# 54 : IkBd-RelA:RelA
const IKBDAA = 54+RECEPTOR_SPECIES;
# 55 : IkBa-RelA:RelAn
const NIKBAAA = 55+RECEPTOR_SPECIES;
# 56 : IkBb-RelA:RelAn
const NIKBBAA = 56+RECEPTOR_SPECIES;
# 57 : IkBe-RelA:RelAn
const NIKBEAA = 57+RECEPTOR_SPECIES;
# 58 : IkBd-RelA:RelAn
const NIKBDAA = 58+RECEPTOR_SPECIES;
# 59 : IkBa-RelA:p50
const IKBAA50 = 59+RECEPTOR_SPECIES;
# 60 : IkBb-RelA:p50
const IKBBA50 = 60+RECEPTOR_SPECIES;
# 61 : IkBe-RelA:p50
const IKBEA50 = 61+RECEPTOR_SPECIES;
# 62 : IkBd-RelA:p50
const IKBDA50 = 62+RECEPTOR_SPECIES;
# 63 : IkBa-RelA:p50n
const NIKBAA50 = 63+RECEPTOR_SPECIES;
# 64 : IkBb-RelA:p50n
const NIKBBA50 = 64+RECEPTOR_SPECIES;
# 65 : IkBe-RelA:p50n
const NIKBEA50 = 65+RECEPTOR_SPECIES;
# 66 : IkBd-RelA:p50n
const NIKBDA50 = 66+RECEPTOR_SPECIES;
# 67 : IkBa-RelA:p52
const IKBAA52 = 67+RECEPTOR_SPECIES;
# 68 : IkBb-RelA:p52
const IKBBA52 = 68+RECEPTOR_SPECIES;
# 69 : IkBe-RelA:p52
const IKBEA52 = 69+RECEPTOR_SPECIES;
# 70 : IkBd-RelA:p52
const IKBDA52 = 70+RECEPTOR_SPECIES;
# 71 : IkBa-RelA:p52n
const NIKBAA52 = 71+RECEPTOR_SPECIES;
# 72 : IkBb-RelA:p52n
const NIKBBA52 = 72+RECEPTOR_SPECIES;
# 73 : IkBe-RelA:p52n
const NIKBEA52 = 73+RECEPTOR_SPECIES;
# 74 : IkBd-RelA:p52n
const NIKBDA52 = 74+RECEPTOR_SPECIES;
# 75 : IkBa-RelB:p50
const IKBAB50 = 75+RECEPTOR_SPECIES;
# 76 : IkBb-RelB:p50
const IKBBB50 = 76+RECEPTOR_SPECIES;
# 77 : IkBe-RelB:p50
const IKBEB50 = 77+RECEPTOR_SPECIES;
# 78 : IkBd-RelB:p50
const IKBDB50 = 78+RECEPTOR_SPECIES;
# 79 : IkBa-RelB:p50n
const NIKBAB50 = 79+RECEPTOR_SPECIES;
# 80 : IkBb-RelB:p50n
const NIKBBB50 = 80+RECEPTOR_SPECIES;
# 81 : IkBe-RelB:p50n
const NIKBEB50 = 81+RECEPTOR_SPECIES;
# 82 : IkBd-RelB:p50n
const NIKBDB50 = 82+RECEPTOR_SPECIES;
# 83 : IkBa-RelB:p52
const IKBAB52 = 83+RECEPTOR_SPECIES;
# 84 : IkBb-RelB:p52
const IKBBB52 = 84+RECEPTOR_SPECIES;
# 85 : IkBe-RelB:p52
const IKBEB52 = 85+RECEPTOR_SPECIES;
# 86 : IkBd-RelB:p52
const IKBDB52 = 86+RECEPTOR_SPECIES;
# 87 : IkBa-RelB:p52n
const NIKBAB52 = 87+RECEPTOR_SPECIES;
# 88 : IkBb-RelB:p52n
const NIKBBB52 = 88+RECEPTOR_SPECIES;
# 89 : IkBe-RelB:p52n
const NIKBEB52 = 89+RECEPTOR_SPECIES;
# 90 : IkBd-RelB:p52n
const NIKBDB52 = 90+RECEPTOR_SPECIES;
# 91 : IkBa-cRel:p50
const IKBAC50 = 91+RECEPTOR_SPECIES;
# 92 : IkBb-cRel:p50
const IKBBC50 = 92+RECEPTOR_SPECIES;
# 93 : IkBe-cRel:p50
const IKBEC50 = 93+RECEPTOR_SPECIES;
# 94 : IkBd-cRel:p50
const IKBDC50 = 94+RECEPTOR_SPECIES;
# 95 : IkBa-cRel:p50n
const NIKBAC50 = 95+RECEPTOR_SPECIES;
# 96 : IkBb-cRel:p50n
const NIKBBC50 = 96+RECEPTOR_SPECIES;
# 97 : IkBe-cRel:p50n
const NIKBEC50 = 97+RECEPTOR_SPECIES;
# 98 : IkBd-cRel:p50n
const NIKBDC50 = 98+RECEPTOR_SPECIES;
# 99 : IkBa-cRel:p52
const IKBAC52 = 99+RECEPTOR_SPECIES;
# 100 : IkBb-cRel:p52
const IKBBC52 = 100+RECEPTOR_SPECIES;
# 101 : IkBe-cRel:p52
const IKBEC52 = 101+RECEPTOR_SPECIES;
# 102 : IkBd-cRel:p52
const IKBDC52 = 102+RECEPTOR_SPECIES;
# 103 : IkBa-cRel:p52n
const NIKBAC52 = 103+RECEPTOR_SPECIES;
# 104 : IkBb-cRel:p52n
const NIKBBC52 = 104+RECEPTOR_SPECIES;
# 105 : IkBe-cRel:p52n
const NIKBEC52 = 105+RECEPTOR_SPECIES;
# 106 : IkBd-cRel:p52n
const NIKBDC52 = 106+RECEPTOR_SPECIES;



# Define indices for Differentiation species
#--------------------------------------------------------------------------------------------------------------------
# 1 : Pax-5
const PAX5 = 1+NFKB_SPECIES+RECEPTOR_SPECIES;
# 2 : Bcl-6
const BCL6 = 2+NFKB_SPECIES+RECEPTOR_SPECIES;
# 3 : Blimp-1
const BLIMP1 = 3+NFKB_SPECIES+RECEPTOR_SPECIES;
# 4 : IRF-4
const IRF4 = 4+NFKB_SPECIES+RECEPTOR_SPECIES;



# MODULE 1: Define indices for NFkB regulatory species
#--------------------------------------------
# 01 : BCL2T
const BCL2T = 1+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 02 : Bcl2
const BCL2 = 2+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;

# MODULE 2: Define indices for CASPASE 8 module (initiator caspase) species
#--------------------------------------------
# 03 : pC8
const PC8 = 3+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 04 : C8
const C8 = 4+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;

# MODULE 3: Define indices for MOMP (pore-forming and transporting) species
#--------------------------------------------
# 05 : Mito
const MITO = 5+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 06 : MOMP (Mitochondrial Outerlayer Permeabilization)
const MOMP = 6+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;

# MODULE 4: Define indices for CASPASE 3 module (effector caspase) species
#--------------------------------------------
# 07 : pC3
const PC3 = 7+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 08 : C3
const C3 = 8+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;

# MODULE 5: Define indices for CASPASE 6 feedback module species
#--------------------------------------------
# 09 : pC6
const PC6 = 9+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 10 : C6
const C6 = 10+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;

# MODULE 6: Define indices for cell death species
#--------------------------------------------
# 11 : PARP
const PARP = 11+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 12 : cPARP
const CPARP = 12+NFKB_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;



# Define indices for Cell Cycle species
#--------------------------------------------------------------------------------------------------------------------
# 1 : tMyc
const TMYC = 1+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 2 : Myc
const MYC = 2+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 3 : tE2F
const TE2F = 3+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 4 : E2F
const E2F = 4+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 5 : Rb
const RB = 5+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 6 : ppRB
const PPRB = 6+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 7 : E2F:Rb
const E2FRB = 7+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 8* : tCycD (deprecated in new model)
const TCYCD = 8+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 9 : CycD
const CYCD = 9+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 10 : CycE
const CYCE = 10+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 11 : CycA
const CYCA = 11+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 12 : CycB
const CYCB = 12+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 13 : p27
const P27 = 13+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 14 : CycD:p27
const CYCDP27 = 14+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 15 : CycE:p27
const CYCEP27 = 15+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 16 : CycA:p27
const CYCAP27 = 16+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 17 : IEP
const IEP = 17+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 18 : PPX
const PPX = 18+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 19 : Cdc20 (total)
const CDC20 = 19+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 20 : Cdc20P (active)
const CDC20P = 20+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 21 : Cdh1
const CDH1 = 21+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
#----------------------------------------
# 22 : GM
const GM = 22+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 23 : Mass
const MASS = 23+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 24 : Rb growth switch
const RBGS = 24+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;
# 25 : Generations
const GEN = 25+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES+RECEPTOR_SPECIES;




# Define species names
#--------------------------------------------------------------------------------------------------------------------
speciesNames = Array{String,1}(undef, TOTAL_SPECIES);

# Receptor layer :
#--------------------------------------------------------------------------------------------------------------------
# 1 : Antigen
speciesNames[ANTIGEN] = "Antigen";
# 2 : B-cell Receptor
speciesNames[BCR] = "B-cell Receptor";
# 3 : Antigen-BCR complex
speciesNames[ABCR] = "Antigen-BCR";
# 4 : inactivated CBM complex
speciesNames[CBM] = "inactivated CBM";
# 5 : activated CBM complex
speciesNames[ACBM] = "activated CBM";
# 6 : inhibited CBM complex
speciesNames[ICBM] = "inhibited CBM";

# 7 : CD40 ligand
speciesNames[CD40L] = "CD40 ligand";
# 8 : CD40 receptor
speciesNames[CD40R] = "CD40 receptor";
# 9 : CD40L-R complex
speciesNames[CD40LR] = "CD40L-R complex";
# 10 : inactivated TRAF6
speciesNames[TRAF6_OFF] = "inactivated TRAF6";
# 11 : activated TRAF6
speciesNames[TRAF6] = "activated TRAF6";

# 12 : inactivated TAK1
speciesNames[TAK1] = "inactivated TAK1";
# 13 : activated TAK1
speciesNames[ATAK1] = "activated TAK1";
# 14 : inactivated IKK
speciesNames[IKK_OFF] = "inactivated IKK";
# 15 : activated pIKK
speciesNames[IKK2] = "activated pIKK";
# 16 : activated ppIKK
speciesNames[IKK3] = "activated ppIKK";
# 17 : inhibited IKK
speciesNames[IIKK] = "inhibited pppIKK";
# 18 : activated IKK
speciesNames[IKK] = "total activated IKK";
# 19 : TRAF3-cIAP1/2 complex
speciesNames[TRAF3] = "TRAF3-cIAP1/2 complex";
# 20 : NIK
speciesNames[NIK] = "NIK";


# NFkB layer :
#--------------------------------------------------------------------------------------------------------------------
# 3 : IkBa transcript
speciesNames[TIKBA] = "IkBa transcript";
# 4 : IkBb transcript
speciesNames[TIKBB] = "IkBb transcript";
# 5 : IkBe transcript
speciesNames[TIKBE] = "IkBe transcript";
# 6 : IkBa protein
speciesNames[IKBA] = "IkBa protein";
# 7 : IkBb protein
speciesNames[IKBB] = "IkBb protein";
# 8 : IkBe protein
speciesNames[IKBE] = "IkBe protein";
# 9 : IkBd protein
speciesNames[IKBD] = "IkBd protein";
# 10 : nuclear IkBa
speciesNames[NIKBA] = "nuclear IkBa";
# 11 : nuclear IkBb
speciesNames[NIKBB] = "nuclear IkBb";
# 12 : nuclear IkBe
speciesNames[NIKBE] = "nuclear IkBe";
# 13 : nuclear IkBd
speciesNames[NIKBD] = "nuclear IkBd";
# 14 : RelA transcript
speciesNames[TRELA] = "RelA transcript";
# 15 : p50 transcript
speciesNames[TP50] = "p50 transcript";
# 16 : RelB transcript
speciesNames[TRELB] = "RelB transcript";
# 17 : P100 transcript
speciesNames[TP100] = "P100 transcript";
# 18 : cRel transcript
speciesNames[TCREL] = "cRel transcript";
# 19 : RelA protein
speciesNames[RELA] = "RelA protein";
# 20 : p50 protein
speciesNames[P50] = "p50 protein";
# 21 : RelB protein
speciesNames[RELB] = "RelB protein";
# 22 : P100 protein
speciesNames[P100] = "P100 protein";
# 23 : cRel protein
speciesNames[CREL] = "cRel protein";
# 24 : P52 protein
speciesNames[P52] = "P52 protein";
# 25 : nuclear RelA
speciesNames[NRELA] = "nuclear RelA";
# 26 : nuclear p50
speciesNames[NP50] = "nuclear p50";
# 27 : nuclear RelB
speciesNames[NRELB] = "nuclear RelB";
# 28 : nuclear P100
speciesNames[NP100] = "nuclear P100";
# 29 : nuclear cRel
speciesNames[NCREL] = "nuclear cRel";
# 30 : nuclear P52
speciesNames[NP52] = "nuclear P52";
# 31 : RelA:RelA protein
speciesNames[AA] = "RelA:RelA protein";
# 32 : RelA:p50 protein
speciesNames[A50] = "RelA:p50 protein";
# 33 : RelA:p52 protein
speciesNames[A52] = "RelA:p52 protein";
# 34 : RelB:p50 protein
speciesNames[B50] = "RelB:p50 protein";
# 35 : RelB:p52 protein
speciesNames[B52] = "RelB:p52 protein";
# 36 : cRel:p50 protein
speciesNames[C50] = "cRel:p50 protein";
# 37 : cRel:p52 protein
speciesNames[C52] = "cRel:p52 protein";
# 38 : cRel:p100 protein
speciesNames[C100] = "cRel:p100 protein";
# 39 : p50:p50 protein
speciesNames[P50P50] = "p50:p50 protein";
# 40 : p52:p52 protein
speciesNames[P52P52] = "p52:p52 protein";
# 41 : nuclear RelA:RelA
speciesNames[NAA] = "nuclear RelA:RelA";
# 42 : nuclear RelA:p50
speciesNames[NA50] = "nuclear RelA:p50";
# 43 : nuclear RelA:p52
speciesNames[NA52] = "nuclear RelA:p52";
# 44 : nuclear RelB:p50
speciesNames[NB50] = "nuclear RelB:p50";
# 45 : nuclear RelB:p52
speciesNames[NB52] = "nuclear RelB:p52";
# 46 : nuclear cRel:p50
speciesNames[NC50] = "nuclear cRel:p50";
# 47 : nuclear cRel:p52
speciesNames[NC52] = "nuclear cRel:p52";
# 48 : nuclear cRel:p100
speciesNames[NC100] = "nuclear cRel:p100";
# 49 : nuclear p50:p50
speciesNames[NP50P50] = "nuclear p50:p50";
# 50 : nuclear p52:p52
speciesNames[NP52P52] = "nuclear p52:p52";
# 51 : IkBa-RelA:RelA protein
speciesNames[IKBAAA] = "IkBa-RelA:RelA protein";
# 52 : IkBb-RelA:RelA protein
speciesNames[IKBBAA] = "IkBb-RelA:RelA protein";
# 53 : IkBe-RelA:RelA protein
speciesNames[IKBEAA] = "IkBe-RelA:RelA protein";
# 54 : IkBd-RelA:RelA protein
speciesNames[IKBDAA] = "IkBd-RelA:RelA protein";
# 55 : nuclear IkBa-RelA:RelA
speciesNames[NIKBAAA] = "nuclear IkBa-RelA:RelA";
# 56 : nuclear IkBb-RelA:RelA
speciesNames[NIKBBAA] = "nuclear IkBb-RelA:RelA";
# 57 : nuclear IkBe-RelA:RelA
speciesNames[NIKBEAA] = "nuclear IkBe-RelA:RelA";
# 58 : nuclear IkBd-RelA:RelA
speciesNames[NIKBDAA] = "nuclear IkBd-RelA:RelA";
# 59 : IkBa-RelA:p50 protein
speciesNames[IKBAA50] = "IkBa-RelA:p50 protein";
# 60 : IkBb-RelA:p50 protein
speciesNames[IKBBA50] = "IkBb-RelA:p50 protein";
# 61 : IkBe-RelA:p50 protein
speciesNames[IKBEA50] = "IkBe-RelA:p50 protein";
# 62 : IkBd-RelA:p50 protein
speciesNames[IKBDA50] = "IkBd-RelA:p50 protein";
# 63 : nuclear IkBa-RelA:p50
speciesNames[NIKBAA50] = "nuclear IkBa-RelA:p50";
# 64 : nuclear IkBb-RelA:p50
speciesNames[NIKBBA50] = "nuclear IkBb-RelA:p50";
# 65 : nuclear IkBe-RelA:p50
speciesNames[NIKBEA50] = "nuclear IkBe-RelA:p50";
# 66 : nuclear IkBd-RelA:p50
speciesNames[NIKBDA50] = "nuclear IkBd-RelA:p50";
# 67 : IkBa-RelA:p52 protein
speciesNames[IKBAA52] = "IkBa-RelA:p52 protein";
# 68 : IkBb-RelA:p52 protein
speciesNames[IKBBA52] = "IkBb-RelA:p52 protein";
# 69 : IkBe-RelA:p52 protein
speciesNames[IKBEA52] = "IkBe-RelA:p52 protein";
# 70 : IkBd-RelA:p52 protein
speciesNames[IKBDA52] = "IkBd-RelA:p52 protein";
# 71 : nuclear IkBa-RelA:p52
speciesNames[NIKBAA52] = "nuclear IkBa-RelA:p52";
# 72 : nuclear IkBb-RelA:p52
speciesNames[NIKBBA52] = "nuclear IkBb-RelA:p52";
# 73 : nuclear IkBe-RelA:p52
speciesNames[NIKBEA52] = "nuclear IkBe-RelA:p52";
# 74 : nuclear IkBd-RelA:p52
speciesNames[NIKBDA52] = "nuclear IkBd-RelA:p52";
# 75 : IkBa-RelB:p50 protein
speciesNames[IKBAB50] = "IkBa-RelB:p50 protein";
# 76 : IkBb-RelB:p50 protein
speciesNames[IKBBB50] = "IkBb-RelB:p50 protein";
# 77 : IkBe-RelB:p50 protein
speciesNames[IKBEB50] = "IkBe-RelB:p50 protein";
# 78 : IkBd-RelB:p50 protein
speciesNames[IKBDB50] = "IkBd-RelB:p50 protein";
# 79 : nuclear IkBa-RelB:p50
speciesNames[NIKBAB50] = "nuclear IkBa-RelB:p50";
# 80 : nuclear IkBb-RelB:p50
speciesNames[NIKBBB50] = "nuclear IkBb-RelB:p50";
# 81 : nuclear IkBe-RelB:p50
speciesNames[NIKBEB50] = "nuclear IkBe-RelB:p50";
# 82 : nuclear IkBd-RelB:p50
speciesNames[NIKBDB50] = "nuclear IkBd-RelB:p50";
# 83 : IkBa-RelB:p52 protein
speciesNames[IKBAB52] = "IkBa-RelB:p52 protein";
# 84 : IkBb-RelB:p52 protein
speciesNames[IKBBB52] = "IkBb-RelB:p52 protein";
# 85 : IkBe-RelB:p52 protein
speciesNames[IKBEB52] = "IkBe-RelB:p52 protein";
# 86 : IkBd-RelB:p52 protein
speciesNames[IKBDB52] = "IkBd-RelB:p52 protein";
# 87 : nuclear IkBa-RelB:p52
speciesNames[NIKBAB52] = "nuclear IkBa-RelB:p52";
# 88 : nuclear IkBb-RelB:p52
speciesNames[NIKBBB52] = "nuclear IkBb-RelB:p52";
# 89 : nuclear IkBe-RelB:p52
speciesNames[NIKBEB52] = "nuclear IkBe-RelB:p52";
# 90 : nuclear IkBd-RelB:p52
speciesNames[NIKBDB52] = "nuclear IkBd-RelB:p52";
# 91 : IkBa-cRel:p50 protein
speciesNames[IKBAC50] = "IkBa-cRel:p50 protein";
# 92 : IkBb-cRel:p50 protein
speciesNames[IKBBC50] = "IkBb-cRel:p50 protein";
# 93 : IkBe-cRel:p50 protein
speciesNames[IKBEC50] = "IkBe-cRel:p50 protein";
# 94 : IkBd-cRel:p50 protein
speciesNames[IKBDC50] = "IkBd-cRel:p50 protein";
# 95 : nuclear IkBa-cRel:p50
speciesNames[NIKBAC50] = "nuclear IkBa-cRel:p50";
# 96 : nuclear IkBb-cRel:p50
speciesNames[NIKBBC50] = "nuclear IkBb-cRel:p50";
# 97 : nuclear IkBe-cRel:p50
speciesNames[NIKBEC50] = "nuclear IkBe-cRel:p50";
# 98 : nuclear IkBd-cRel:p50
speciesNames[NIKBDC50] = "nuclear IkBd-cRel:p50";
# 99 : IkBa-cRel:p52 protein
speciesNames[IKBAC52] = "IkBa-cRel:p52 protein";
# 100 : IkBb-cRel:p52 protein
speciesNames[IKBBC52] = "IkBb-cRel:p52 protein";
# 101 : IkBe-cRel:p52 protein
speciesNames[IKBEC52] = "IkBe-cRel:p52 protein";
# 102 : IkBd-cRel:p52 protein
speciesNames[IKBDC52] = "IkBd-cRel:p52 protein";
# 103 : nuclear IkBa-cRel:p52
speciesNames[NIKBAC52] = "nuclear IkBa-cRel:p52";
# 104 : nuclear IkBb-cRel:p52
speciesNames[NIKBBC52] = "nuclear IkBb-cRel:p52";
# 105 : nuclear IkBe-cRel:p52
speciesNames[NIKBEC52] = "nuclear IkBe-cRel:p52";
# 106 : nuclear IkBd-cRel:p52
speciesNames[NIKBDC52] = "nuclear IkBd-cRel:p52";


# Apoptosis layer :
#--------------------------------------------------------------------------------------------------------------------
# 1 : RelA
speciesNames[RELA] = "RelA";
# 2 : cRel
speciesNames[CREL] = "cRel";
# 3 : B-cell Receptor
speciesNames[BCR] = "B-cell Receptor";
# 4 : Bcl2t
speciesNames[BCL2T] = "Bcl2 transcript";
# 5 : Bcl2
speciesNames[BCL2] = "Bcl2 protein";
# 6 : pC8
speciesNames[PC8] = "Procaspase 8";
# 7 : C8
speciesNames[C8] = "Caspase 8";
# 8 : Mitochondria
speciesNames[MITO] = "Mitochondria";
# 9 : Mitochondrial outer membrane permeabilization
speciesNames[MOMP] = "MOMP";
# 10 : pC3
speciesNames[PC3] = "Procaspase 3";
# 11 : C3
speciesNames[C3] = "Caspase 3";
# 12 : pC6
speciesNames[PC6] = "Procaspase 6";
# 9 : C6
speciesNames[C6] = "Caspase 6";
# 10 : PARP
speciesNames[PARP] = "PARP";
# 11 : cPARP
speciesNames[CPARP] = "cleaved PARP";

#------------------------------------------------
# Differentiation layer :
# 62 : Pax-5
speciesNames[PAX5] = "Pax-5";
# 63 : Bcl-6
speciesNames[BCL6] = "Bcl-6";
# 64 : Blimp-1
speciesNames[BLIMP1] = "Blimp-1";
# 65 : IRF-4
speciesNames[IRF4] = "IRF-4";
#------------------------------------------------
# Cell cycle layer :
# 66 : tMyc
speciesNames[TMYC] = "Myc transcript";
# 67 : Myc
speciesNames[MYC] = "Myc protein";
# 68 : tE2F
speciesNames[TE2F] = "E2F transcript";
# 69 : E2F
speciesNames[E2F] = "E2F protein";
# 70 : Rb
speciesNames[RB] = "Retinoblastoma";
# 71 : ppRB
speciesNames[PPRB] = "Phosphorylated Rb";
# 72 : E2F:Rb
speciesNames[E2FRB] = "E2F:Rb complex";
# 73 : CycD
speciesNames[CYCD] = "Cyclin D protein";
# 73 : tCycD
speciesNames[TCYCD] = "Cyclin D transcript";
# 74 : CycE
speciesNames[CYCE] = "Cyclin E protein";
# 75 : CycA
speciesNames[CYCA] = "Cyclin A protein";
# 76 : CycB
speciesNames[CYCB] = "Cyclin B protein";
# 77 : p27
speciesNames[P27] = "p27 protein";
# 78 : CycD:p27
speciesNames[CYCDP27] = "Cyclin D:p27 complex";
# 79 : CycE:p27
speciesNames[CYCEP27] = "Cyclin E:p27 complex";
# 80 : CycA:p27
speciesNames[CYCAP27] = "Cyclin A:p27 complex";
# 81 : IEP
speciesNames[IEP] = "IE Phosphorylated";
# 82 : PPX
speciesNames[PPX] = "Generic phosphatase";
# 83 : Cdc20 (total)
speciesNames[CDC20] = "Cdc20 protein";
# 84 : Cdc20P (active)
speciesNames[CDC20P] = "Active form Cdc20 protein";
# 85 : Cdh1
speciesNames[CDH1] = "Cdh1 protein";
#------------------------------------------------
# 86 : GM
speciesNames[GM] = "General Machinery of the cell";
# 87 : Mass
speciesNames[MASS] = "Mass of cell";
# 88 : RBGS
speciesNames[RBGS] = "[Rb] growth switch";
# 89 : Generation number
speciesNames[GEN] = "Generation";
