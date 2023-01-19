# Define global constants for model
#---------------------------------------------
const CONVERSION = 1;   # unit conversion (default in hr-1). Set to 1/60 or 1/3600 to convert to min-1 or sec-1.

const EPS = 1;    # translation modifier
const BASAL_IKK = 0.001;
const TOTAL_IKK = 288.48;    # unit: nM (170 in RoyMitchell 2019, 288.48 in Mitchell 2018)
const BASAL_TBCL2 = 277.2133 * 60 * CONVERSION;    # basal Bcl2t transcription, unit: molecules / hr
const cRelFeedback = 10;    # cRel induction modifier
const constituitiveMultiplier = (1+(6/150))/(1+cRelFeedback*(6/150));
const IKK_TYPE = "CpG2";

# thresholds for tMyC, tBcl2, tCycD
const MYCTHR = 40;
const BCL2THR = 40;
const CYCDTHR = 40;
const GROWTHR = 40;



# Parameter distributions constants
#---------------------------------------------
const FOUNDER_RATE_CV = 0.112; # a CV of 11.2% in rate constants result in 32% CV in concentrations in steady state
const FOUNDER_CONC_CV = 0.00; # initial concentrations don't get distributed directly -- they are distributed due to distribution of rate constants after the cell reach steady states

const DAUGHTER_RATE_CV = 0.00;
const DAUGHTER_CONC_CV = 0.00;
const DAUGHTER_PARTITION_CV = 0.072; # how asymmetric each division needs to be

const FOUNDER_CELL_NUM = 8; # 125 in Mitchell et al. 2018?
const MAX_GEN = 7; # maximum generation
const BURN_IN_PERIOD = 24; # hrs, pre-simulation phase
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
const NFKB_SPECIES = 62;
const DIFF_SPECIES = 4;
const APOPTOSIS_SPECIES = 59;
const PROLIF_SPECIES = 25;
const TOTAL_SPECIES = NFKB_SPECIES + DIFF_SPECIES + APOPTOSIS_SPECIES + PROLIF_SPECIES;



# Define indices for NFkB species
#-------------------------------------------------------------------------------------------------------------------
# 1 : IKK
const IKK = 1;

# MODULE 1: Define indices for IkBs & their mRNA transcripts
#-----------------------------------------------------------------------
# 2 : tIkBa
const TIKBA = 2;
# 3 : tIkBb
const TIKBB = 3;
# 4 : tIkBe
const TIKBE = 4;
# 5 : tIkBd
const TIKBD = 5;
#---------------------
# 6 : IkBa
const IKBA = 6;
# 7 : IkBb
const IKBB = 7;
# 8 : IkBe
const IKBE = 8;
# 9 : IkBd
const IKBD = 9;
#---------------------
# 10 : nIkBa
const NIKBA = 10;
# 11 : nIkBb
const NIKBB = 11;
# 12 : nIkBe
const NIKBE = 12;
# 13 : nIkBd
const NIKBD = 13;

# MODULE 2: Define indices for NFkB monomers & their mRNA transcripts
#-----------------------------------------------------------------------
# 14 : tRelA
const TRELA = 14;
# 15 : tP50
const TP50 = 15;
# 16 : tcRel
const TCREL = 16;
#---------------------
# 17 : RelA
const RELA = 17;
# 18 : P50
const P50 = 18;
# 19 : cRel
const CREL = 19;
#---------------------
# 20 : nRelA
const NRELA = 20;
# 21 : nP50
const NP50 = 21;
# 22 : ncRel
const NCREL = 22;

# MODULE 3: Define indices for NFkB dimers & their IkB complexes
#-----------------------------------------------------------------------
# 23 : RelA:RelA
const AA = 23;
# 24 : RelA:p50
const A50 = 24;
# 25 : p50:p50
const P50P50 = 25;
# 26 : cRel:p50
const C50 = 26;
#---------------------
# 27 : nRelA:RelA
const NAA = 27;
# 28 : nRelA:p50
const NA50 = 28;
# 29 : np50:p50
const NP50P50 = 29;
# 30 : ncRel:p50
const NC50 = 30;
#---------------------
# 31 : IkBa-RelA:RelA
const IKBAAA = 31;
# 32 : IkBb-RelA:RelA
const IKBBAA = 32;
# 33 : IkBe-RelA:RelA
const IKBEAA = 33;
# 34 : IkBd-RelA:RelA
const IKBDAA = 34;

# 35 : nIkBa-RelA:RelA
const NIKBAAA = 35;
# 36 : nIkBb-RelA:RelA
const NIKBBAA = 36;
# 37 : nIkBe-RelA:RelA
const NIKBEAA = 37;
# 38 : nIkBd-RelA:RelA
const NIKBDAA = 38;
#---------------------
# 39 : IkBa-RelA:p50
const IKBAA50 = 39;
# 40 : IkBb-RelA:p50
const IKBBA50 = 40;
# 41 : IkBe-RelA:p50
const IKBEA50 = 41;
# 42 : IkBd-RelA:p50
const IKBDA50 = 42;

# 43 : nIkBa-RelA:p50
const NIKBAA50 = 43;
# 44 : nIkBb-RelA:p50
const NIKBBA50 = 44;
# 45 : nIkBe-RelA:p50
const NIKBEA50 = 45;
# 46 : nIkBd-RelA:p50
const NIKBDA50 = 46;
#---------------------
# 47 : IkBa-p50:p50
const IKBA5050 = 47;
# 48 : IkBb-p50:p50
const IKBB5050 = 48;
# 49 : IkBe-p50:p50
const IKBE5050 = 49;
# 50 : IkBd-p50:p50
const IKBD5050 = 50;

# 51 : nIkBa-p50:p50
const NIKBA5050 = 51;
# 52 : nIkBb-p50:p50
const NIKBB5050 = 52;
# 53 : nIkBe-p50:p50
const NIKBE5050 = 53;
# 54 : nIkBd-p50:p50
const NIKBD5050 = 54;
#---------------------
# 55 : IkBa-cRel:p50
const IKBAC50 = 55;
# 56 : IkBb-cRel:p50
const IKBBC50 = 56;
# 57 : IkBe-cRel:p50
const IKBEC50 = 57;
# 58 : IkBd-cRel:p50
const IKBDC50 = 58;

# 59 : nIkBa-cRel:p50
const NIKBAC50 = 59;
# 60 : nIkBb-cRel:p50
const NIKBBC50 = 60;
# 61 : nIkBe-cRel:p50
const NIKBEC50 = 61;
# 62 : nIkBd-cRel:p50
const NIKBDC50 = 62;



# Define indices for Differentiation species
#--------------------------------------------------------------------------------------------------------------------
# 1 : Pax-5
const PAX5 = 1+NFKB_SPECIES;
# 2 : Bcl-6
const BCL6 = 2+NFKB_SPECIES;
# 3 : Blimp-1
const BLIMP1 = 3+NFKB_SPECIES;
# 4 : IRF-4
const IRF4 = 4+NFKB_SPECIES;



# Define indices for Apoptosis species
#--------------------------------------------------------------------------------------------------------------------
# MODULE 1: Define indices for NFkB regulatory species
#--------------------------------------------
# 1 : tBCL2
const TBCL2 = 1+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 2: Define indices for Receptor & DISC system species
#--------------------------------------------
# 2 : L
const L = 2+NFKB_SPECIES+DIFF_SPECIES;
# 3 : R
const R = 3+NFKB_SPECIES+DIFF_SPECIES;
# 4 : L:R
const LR = 4+NFKB_SPECIES+DIFF_SPECIES;
# 5 : DISC
const DISC = 5+NFKB_SPECIES+DIFF_SPECIES;
# 6 : flip
const FLIP = 6+NFKB_SPECIES+DIFF_SPECIES;
# 7 : flip:DISC
const FLIPDISC = 7+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 3: Define indices for CASPASE 8 module (initiator caspase) species
#--------------------------------------------
# 8 : pC8
const PC8 = 8+NFKB_SPECIES+DIFF_SPECIES;
# 9 : DISC:pC8
const DISCPC8 = 9+NFKB_SPECIES+DIFF_SPECIES;
# 10 : C8
const C8 = 10+NFKB_SPECIES+DIFF_SPECIES;
# 11 : Bar
const BAR = 11+NFKB_SPECIES+DIFF_SPECIES;
# 12 : Bar:C8
const BARC8 = 12+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 4: Define indices for species before Bax enter mitochondria
#--------------------------------------------
# 13 : Bid
const BID = 13+NFKB_SPECIES+DIFF_SPECIES;
# 14 : C8:Bid
const C8BID = 14+NFKB_SPECIES+DIFF_SPECIES;
# 15 : tBid
const TBID = 15+NFKB_SPECIES+DIFF_SPECIES;
# 16 : cBcl2
const CBCL2 = 16+NFKB_SPECIES+DIFF_SPECIES;
# 17 : cBcl2:tBid
const CBCL2TBID = 17+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 5: Define indices for species involved in the oligomerization of Bax
#--------------------------------------------
# 18 : Bax
const BAX = 18+NFKB_SPECIES+DIFF_SPECIES;
# 19 : tBid:Bax
const TBIDBAX = 19+NFKB_SPECIES+DIFF_SPECIES;
# 20 : aBax
const ABAX = 20+NFKB_SPECIES+DIFF_SPECIES;
# 21 : mBax
const MBAX = 21+NFKB_SPECIES+DIFF_SPECIES;
# 22 : Bcl2
const BCL2 = 22+NFKB_SPECIES+DIFF_SPECIES;
# 23 : mBax:Bcl2
const MBAXBCL2 = 23+NFKB_SPECIES+DIFF_SPECIES;
# 24 : Bax2
const BAX2 = 24+NFKB_SPECIES+DIFF_SPECIES;
# 25 : Bax2:Bcl2
const BAX2BCL2 = 25+NFKB_SPECIES+DIFF_SPECIES;
# 26 : Bax4
const BAX4 = 26+NFKB_SPECIES+DIFF_SPECIES;
# 27 : Bax4:Bcl2
const BAX4BCL2 = 27+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 6: Define indices for MOMP (pore-forming and transporting) species
#--------------------------------------------
# 28 : Mito
const MITO = 28+NFKB_SPECIES+DIFF_SPECIES;
# 29 : Bax4:Mito
const BAX4MITO = 29+NFKB_SPECIES+DIFF_SPECIES;
# 30 : aMito
const AMITO = 30+NFKB_SPECIES+DIFF_SPECIES;
# 31 : mCytoC
const MCYTOC = 31+NFKB_SPECIES+DIFF_SPECIES;
# 32 : aMito:mCytoC
const AMITOMCYTOC = 32+NFKB_SPECIES+DIFF_SPECIES;
# 33 : aCytoC
const ACYTOC = 33+NFKB_SPECIES+DIFF_SPECIES;
# 34 : mSmac
const MSMAC = 34+NFKB_SPECIES+DIFF_SPECIES;
# 35 : aMito:mSmac
const AMITOMSMAC = 35+NFKB_SPECIES+DIFF_SPECIES;
# 36 : aSmac
const ASMAC = 36+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 7: Define indices for Feedforward pathway 1 (XIAP) species
#--------------------------------------------
# 37 : XIAP
const XIAP = 37+NFKB_SPECIES+DIFF_SPECIES;
# 38 : cSmac
const CSMAC = 38+NFKB_SPECIES+DIFF_SPECIES;
# 39 : cSmac:XIAP
const CSMACXIAP = 39+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 8: Define indices for Feedforward pathway 2 (Apoptosome) species
#--------------------------------------------
# 40 : cCytoC
const CCYTOC = 40+NFKB_SPECIES+DIFF_SPECIES;
# 41 : Apaf
const APAF = 41+NFKB_SPECIES+DIFF_SPECIES;
# 42 : Apaf:cCytoC
const APAFCCYTOC = 42+NFKB_SPECIES+DIFF_SPECIES;
# 43 : aApaf
const AAPAF = 43+NFKB_SPECIES+DIFF_SPECIES;
# 44 : pC9
const PC9 = 44+NFKB_SPECIES+DIFF_SPECIES;
# 45 : Apop
const APOP = 45+NFKB_SPECIES+DIFF_SPECIES;
# 46 : Apop:XIAP
const APOPXIAP = 46+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 9: Define indices for CASPASE 3 module (effector caspase) species
#--------------------------------------------
# 47 : pC3
const PC3 = 47+NFKB_SPECIES+DIFF_SPECIES;
# 48 : C8:pC3
const C8PC3 = 48+NFKB_SPECIES+DIFF_SPECIES;
# 49 : C3
const C3 = 49+NFKB_SPECIES+DIFF_SPECIES;
# 50 : XIAP:C3
const XIAPC3 = 50+NFKB_SPECIES+DIFF_SPECIES;
# 51 : Apop:pC3
const APOPPC3 = 51+NFKB_SPECIES+DIFF_SPECIES;
# 52 : uC3
const UC3 = 52+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 10: Define indices for CASPASE 6 feedback module species
#--------------------------------------------
# 53 : pC6
const PC6 = 53+NFKB_SPECIES+DIFF_SPECIES;
# 54 : C3:pC6
const C3PC6 = 54+NFKB_SPECIES+DIFF_SPECIES;
# 55 : C6
const C6 = 55+NFKB_SPECIES+DIFF_SPECIES;
# 56 : C6:pC8
const C6PC8 = 56+NFKB_SPECIES+DIFF_SPECIES;

# MODULE 11: Define indices for cell death species
#--------------------------------------------
# 57 : PARP
const PARP = 57+NFKB_SPECIES+DIFF_SPECIES;
# 58 : C3:PARP
const C3PARP = 58+NFKB_SPECIES+DIFF_SPECIES;
# 59 : tPARP
const TPARP = 59+NFKB_SPECIES+DIFF_SPECIES;



# Define indices for Cell Cycle species
#--------------------------------------------------------------------------------------------------------------------
# 1 : tMyc
const TMYC = 1+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 2 : Myc
const MYC = 2+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 3 : tE2F
const TE2F = 3+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 4 : E2F
const E2F = 4+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 5 : Rb
const RB = 5+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 6 : ppRB
const PPRB = 6+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 7 : E2F:Rb
const E2FRB = 7+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 8* : tCycD (deprecated in new model)
const TCYCD = 8+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 9 : CycD
const CYCD = 9+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 10 : CycE
const CYCE = 10+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 11 : CycA
const CYCA = 11+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 12 : CycB
const CYCB = 12+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 13 : p27
const P27 = 13+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 14 : CycD:p27
const CYCDP27 = 14+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 15 : CycE:p27
const CYCEP27 = 15+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 16 : CycA:p27
const CYCAP27 = 16+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 17 : IEP
const IEP = 17+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 18 : PPX
const PPX = 18+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 19 : Cdc20 (total)
const CDC20 = 19+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 20 : Cdc20P (active)
const CDC20P = 20+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 21 : Cdh1
const CDH1 = 21+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
#----------------------------------------
# 22 : GM
const GM = 22+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 23 : Mass
const MASS = 23+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 24 : Rb growth switch
const RBGS = 24+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
# 25 : Generations
const GEN = 25+NFKB_SPECIES+APOPTOSIS_SPECIES+DIFF_SPECIES;
