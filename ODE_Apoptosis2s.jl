# 2nd edition of Apoptosis model (added utilities to solve ODE separately from NFkB model), simplified module
# include("HelperFunctions.jl");
# Compute reaction fluxes for apoptosis module
#--------------------------------------------

function computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase, time, birthday, inputCurves)
    @inbounds begin
        #--------------------------------------------
        # MODULE 1: Bcl2 mRNA transcript
        #--------------------------------------------
        if phase == 2
            concentration[ABCR] = inputCurves(time + birthday, idxs=ABCR);
            concentration[NA50] = inputCurves(time + birthday, idxs=NA50);
            concentration[NC50] = inputCurves(time + birthday, idxs=NC50);
            concentration[IKK] = inputCurves(time + birthday, idxs=IKK);
        end
        # 01. 0 --nNFkB--> tBcl2 activity : Activation (((kA50*[A50n]+kC50*[C50n]+kIKK*[IKK])/Kd)^expn/(k1+((kA50*[A50n]+kC50*[C50n]+kIKK*[IKK])/Kd)^expn)*(1-k2)+k2 ----> IKK activity represents other transcription factor
        # reactionFlux[BCL2T, ACTIVATION] = hillActivation(Ka = 1.0, na = 2.0, X = (0.3*concentration[NA50] + 0.6*concentration[NC50] + 0.1*concentration[IKK])/BCL2THR) * (1-Srates[BCL2T, ACTIVATION]) + Srates[BCL2T, ACTIVATION];
        reactionFlux[BCL2T, ACTIVATION] = hillActivation(Ka = 1.0, na = 2.0, X = (0.3*(concentration[NA50]+concentration[NA52]) + 0.6*(concentration[NC50]+concentration[NC52]) + 0.1*concentration[IKK])/BCL2THR) * (1-Srates[TBCL2, ACTIVATION]) + Srates[TBCL2, ACTIVATION];
        # 02. tBcl2 Act ----> tBcl2 : Synthesis
        reactionFlux[BCL2T, SYNTHESIS] = Srates[BCL2T, SYNTHESIS] * max(reactionFlux[BCL2T, ACTIVATION], 0.01);
        # 03. tBcl2 ----> 0 : Basal decay
        reactionFlux[BCL2T, BASALDECAY] = Srates[BCL2T, BASALDECAY] * concentration[BCL2T];
        # 04. tBcl2 ----> Bcl2 : translation
        if phase == 1 # assume that during equilibration the cell is receiving a survival signal
            reactionFlux[BCL2, TRANSLATION] = EPS * BASAL_TBCL2;
        else
            reactionFlux[BCL2, TRANSLATION] = EPS * Srates[BCL2, TRANSLATION] * concentration[BCL2T];
        end
        # 05. Bcl2 ----> 0 : Basal decay
        reactionFlux[BCL2, BASALDECAY] = Srates[BCL2, BASALDECAY] * concentration[BCL2];


        #--------------------------------------------
        # MODULE 2: CASPASE 8 module (initiator caspase)
        #--------------------------------------------
        # 06. 0 ----> pC8 : Basal synthesis
        reactionFlux[PC8, BASALSYNTHESIS] = Srates[PC8, BASALSYNTHESIS];
        # 07. pC8 ----> 0 : Basal decay
        reactionFlux[PC8, BASALDECAY] = Srates[PC8, BASALDECAY] * concentration[PC8];
        # 08. pC8 ----> C8 : Basal synthesis (constitutive pC8 processing)
        reactionFlux[C8, BASALSYNTHESIS] = Srates[C8, BASALSYNTHESIS] * concentration[PC8];
        # 09. pC8 + BCR ----> C8 + BCR : BCR-mediated pC8 processing
        reactionFlux[C8, ACTIVATION] = Srates[C8, ACTIVATION] * concentration[PC8] * concentration[ABCR];
        # 10. C8 ----> 0 : Basal decay
        reactionFlux[C8, BASALDECAY] = Srates[C8, BASALDECAY] * concentration[C8];


        #----------------------------------------
        # MODULE 3: MOMP (pore-forming and transporting)
        #----------------------------------------
        # 11. Mito ----> MOMP : Caspase-activated MOMP
        reactionFlux[MOMP, ACTIVATION] = Srates[MOMP, ACTIVATION] * concentration[MITO] * concentration[C8];
        # 12. MOMP ----> Mito : Basal decay
        reactionFlux[MOMP, BASALDECAY] = Srates[MOMP, BASALDECAY] * concentration[MOMP];
        # 12. MOMP ----> Mito : Bcl2-activated MOMP deactivation
        reactionFlux[MOMP, DEACTIVATION] = Srates[MOMP, DEACTIVATION] * concentration[MOMP] * concentration[BCL2];

        #----------------------------------------
        # MODULE 4: CASPASE 3 module (effector caspase)
        #----------------------------------------
        # 13. 0 ----> pC3 : Basal synthesis
        reactionFlux[PC3, BASALSYNTHESIS] = Srates[PC3, BASALSYNTHESIS];
        # 14. pC3 ----> 0 : Basal decay
        reactionFlux[PC3, BASALDECAY] = Srates[PC3, BASALDECAY] * concentration[PC3];
        # 15. pC3 + C8 ----> C3 + C8 : Activation
        reactionFlux[C3, ACTIVATION] = Srates[C3, ACTIVATION] * concentration[PC3] * concentration[C8];
        # 16. pC3 + MOMP ----> C3 + MOMP : synthesis
        reactionFlux[C3, SYNTHESIS] = Srates[C3, SYNTHESIS] * concentration[PC3] * hillActivation(Ka = 500.0, na = 2.0, X = concentration[MOMP]);
        # 17. C3 ----> 0 : Basal decay
        reactionFlux[C3, BASALDECAY] = Srates[C3, BASALDECAY] * concentration[C3];

        #----------------------------------------
        # MODULE 5: CASPASE 6 feedback module
        #----------------------------------------
        # 18. 0 ----> pC6 : Basal synthesis
        reactionFlux[PC6, BASALSYNTHESIS] = Srates[PC6, BASALSYNTHESIS];
        # 19. pC6 ----> 0 : Basal decay
        reactionFlux[PC6, BASALDECAY] = Srates[PC6, BASALDECAY] * concentration[PC6];
        # 20. pC6 + C3 ----> C6 + C3 : C3-mediated pC6 processing (C6 ACTIVATION)
        reactionFlux[C6, ACTIVATION] = Srates[C6, ACTIVATION] * concentration[C3] * concentration[PC6];
        # 21. C6 ----> 0 : Basal decay
        reactionFlux[C6, BASALDECAY] = Srates[C6, BASALDECAY] * concentration[C6];

        # 22. pC8 + C6 ----> C8 + C6 : C6-mediated pC8 processing (C8 SYNTHESIS)
        reactionFlux[C8, SYNTHESIS] = Srates[C8, SYNTHESIS] * concentration[PC8] * concentration[C6];

        #----------------------------------------
        # MODULE 6: Cell death module
        #----------------------------------------
        # 23. 0 ----> PARP : Basal synthesis
        reactionFlux[PARP, BASALSYNTHESIS] = Srates[PARP, BASALSYNTHESIS];
        # 24. PARP ----> 0 : Basal decay
        reactionFlux[PARP, BASALDECAY] = Srates[PARP, BASALDECAY] * concentration[PARP];
        # 25. C3 + PARP ----> C3 + CPARP : PARP deactivation
        reactionFlux[PARP, DEACTIVATION] = Srates[PARP, DEACTIVATION] * concentration[PARP] * concentration[C3];
        # 26. cPARP ----> 0 : Basal decay
        reactionFlux[CPARP, BASALDECAY] = Srates[CPARP, BASALDECAY] * concentration[CPARP];

    end
    nothing
end

function ApoptosisNettFluxes!(nettFlux, reactionFlux)

    # Compute apoptosis net fluxes
    @inbounds begin
        #------------------------------------------------
        # MODULE 1
        #--------------------------------------------
        # 01 : Bcl2t
        nettFlux[BCL2T] = reactionFlux[BCL2T, SYNTHESIS] - reactionFlux[BCL2T, BASALDECAY];
        # 02 : Bcl2
        nettFlux[BCL2] = reactionFlux[BCL2, TRANSLATION] - reactionFlux[BCL2, BASALDECAY];

        # MODULE 3
        #--------------------------------------------
        # 03 : pC8
        nettFlux[PC8] = reactionFlux[PC8, BASALSYNTHESIS] - reactionFlux[PC8, BASALDECAY] - reactionFlux[C8, BASALSYNTHESIS] - reactionFlux[C8, ACTIVATION] - reactionFlux[C8, SYNTHESIS];
        # 04 : C8
        nettFlux[C8] = reactionFlux[C8, BASALSYNTHESIS] + reactionFlux[C8, ACTIVATION] + reactionFlux[C8, SYNTHESIS] - reactionFlux[C8, BASALDECAY];

        # MODULE 6
        #--------------------------------------------
        # 05 : Mito
        nettFlux[MITO] = reactionFlux[MOMP, DEACTIVATION] - reactionFlux[MOMP, ACTIVATION] + reactionFlux[MOMP, BASALDECAY];
        # 06 : MOMP
        nettFlux[MOMP] = reactionFlux[MOMP, ACTIVATION] - reactionFlux[MOMP, DEACTIVATION] - reactionFlux[MOMP, BASALDECAY];

        # MODULE 9
        #--------------------------------------------
        # 07 : pC3
        nettFlux[PC3] = reactionFlux[PC3, BASALSYNTHESIS] - reactionFlux[PC3, BASALDECAY] - reactionFlux[C3, ACTIVATION];
        # 08 : C3
        nettFlux[C3] = reactionFlux[C3, ACTIVATION] + reactionFlux[C3, SYNTHESIS] - reactionFlux[C3, BASALDECAY];

        # MODULE 10
        #--------------------------------------------
        # 09 : pC6
        nettFlux[PC6] = reactionFlux[PC6, BASALSYNTHESIS] - reactionFlux[PC6, BASALDECAY] - reactionFlux[C6, ACTIVATION];
        # 10 : C6
        nettFlux[C6] = reactionFlux[C6, ACTIVATION] - reactionFlux[C6, BASALDECAY];

        # MODULE 11
        #--------------------------------------------
        # 11 : PARP
        nettFlux[PARP] = reactionFlux[PARP, BASALSYNTHESIS] - reactionFlux[PARP, BASALDECAY] - reactionFlux[PARP, DEACTIVATION];
        # 12 : cPARP
        nettFlux[CPARP] = reactionFlux[PARP, DEACTIVATION] - reactionFlux[CPARP, BASALDECAY];
        #------------------------------------------------

    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute Apoptosis reaction fluxes
    computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase, nothing, nothing, nothing);
    # Compute Apoptosis net fluxes
    ApoptosisNettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, phase = 2)
function computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase, time, birthday, inputCurves)
    #------------------------------------------------
    # Compute Apoptosis reaction fluxes
    computeApoptosisFluxes!(concentration, reactionFlux, Srates, phase, time, birthday, inputCurves);
    # Compute Apoptosis net fluxes
    ApoptosisNettFluxes!(nettFlux, reactionFlux);
    nothing
end
