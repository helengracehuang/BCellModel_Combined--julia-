# 2nd edition of NFkB model (optimized with for loops instead of dot operations)
# include("HelperFunctions.jl");
# Compute reaction fluxes for NFkB module
#--------------------------------------------

function computeNFkBFluxes!(concentration, reactionFlux, Srates, phase; delay=nothing, time=nothing, birthday=nothing, IKKCurve=nothing, historicFlux=nothing)
    @inbounds begin
        # get historic values for nuclear NFkB concentrations (for delayed transcription)
        if phase == 2
            p = (birthday, IKKCurve, Srates, reactionFlux, historicFlux);
            # hist_NAA_075 = delay(p, time-0.75/CONVERSION; idxs=NAA);
            # hist_NA50_075 = delay(p, time-0.75/CONVERSION; idxs=NA50);
            # hist_NC50_075 = delay(p, time-0.75/CONVERSION; idxs=NC50);
            #
            # hist_NA50_1 = delay(p, time-1.0/CONVERSION; idxs=NA50);
            # hist_NC50_1 = delay(p, time-1.0/CONVERSION; idxs=NC50);
            #
            # hist_NAA_15 = delay(p, time-1.5/CONVERSION; idxs=NAA);
            # hist_NA50_15 = delay(p, time-1.5/CONVERSION; idxs=NA50);

            delay(historicFlux, p, time-0.75/CONVERSION);
            hist_NAA_075 = historicFlux[NAA];
            hist_NA50_075 = historicFlux[NA50];
            hist_NC50_075 = historicFlux[NC50];
            delay(historicFlux, p, time-1.0/CONVERSION);
            hist_NA50_1 = historicFlux[NA50];
            hist_NC50_1 = historicFlux[NC50];
            delay(historicFlux, p, time-1.5/CONVERSION);
            hist_NAA_15 = historicFlux[NAA];
            hist_NA50_15 = historicFlux[NA50];
        end

        #-----------------------------------------------
        #--------------------------------------------
        # MODULE 1: IKK input
        #--------------------------------------------
        # 01. Prescribe IKK concentration according to input functions
        if phase == 1
            concentration[IKK] = BASAL_IKK;
        else
            concentration[IKK] = IKKCurve(time + birthday);
        end
        # 02. IKK : Deactivation (basal IKK activity)
        reactionFlux[IKK, DEACTIVATION] = EPS * BASAL_IKK * TOTAL_IKK;
        # 03. IKK : Activation (IKK activity)
        reactionFlux[IKK, ACTIVATION] = EPS * concentration[IKK] * TOTAL_IKK;
        # line 83 in matlab ODE (not sure about the linear IKK activity?)

        #--------------------------------------------
        # MODULE 2: IkB Reactions
        #--------------------------------------------
        # IkB mRNA
        # 0 --nNFkB--> tIkB : nuclear NFkB dimer-induced synthesis
        # const rate * ( 1 + sum( w*([d]/Kd)^hill)) / ( 1 + sum( ([d]/Kd)^hill))
        reactionFlux[TIKBA, BASALSYNTHESIS] = Srates[TIKBA, BASALSYNTHESIS] * hillInduction(25, 200, 1, concentration[NAA], concentration[NA50], concentration[NC50]; Kd = 150.0, hill = 1.1);
        reactionFlux[TIKBB, BASALSYNTHESIS] = Srates[TIKBB, BASALSYNTHESIS];
        reactionFlux[TIKBE, BASALSYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS] * hillInduction(25, 25, 250, concentration[NAA], concentration[NA50], concentration[NC50]; Kd = 150.0, hill = 1.1);
        reactionFlux[TIKBD, BASALSYNTHESIS] = Srates[TIKBD, BASALSYNTHESIS] * hillInduction(200, 200, concentration[NAA], concentration[NA50]; Kd = 150.0, hill = 1.1);
        if phase == 1  # pre-simulation phase (no delay)
            for i in TIKBA:TIKBD
                reactionFlux[i, SYNTHESIS] = reactionFlux[i, BASALSYNTHESIS];
            end
        else
            reactionFlux[TIKBA, SYNTHESIS] = reactionFlux[TIKBA, BASALSYNTHESIS]; # no delay for tIkba/b
            reactionFlux[TIKBB, SYNTHESIS] = reactionFlux[TIKBB, BASALSYNTHESIS];
            if time > 0.75
                reactionFlux[TIKBE, SYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS] * hillInduction(25, 25, 250, hist_NAA_075, hist_NA50_075, hist_NC50_075; Kd = 150.0, hill = 1.1);
            else # transcription is at basal level
                reactionFlux[TIKBE, SYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS];
            end
            if time > 1.5
                reactionFlux[TIKBD, SYNTHESIS] = Srates[TIKBD, BASALSYNTHESIS] * hillInduction(200, 200, hist_NAA_15, hist_NA50_15; Kd = 150.0, hill = 1.1);
            else # transcription is at basal level
                reactionFlux[TIKBD, SYNTHESIS] = Srates[TIKBE, BASALSYNTHESIS];
            end
        end
        # tIkB ----> 0 : Basal decay
        for i in TIKBA:TIKBD
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # IkB free in cytoplasm
        for i in IKBA:IKBD
            # tIkB ----> IkB : translation
            reactionFlux[i, TRANSLATION] = EPS * Srates[i, TRANSLATION] * concentration[i-IKBA+TIKBA];
            # IkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
            # IkB + IKK ----> IKK : IKK-mediated free IkB decay (basal IKK only for IkBd)
            reactionFlux[i, DECAY] = reactionFlux[IKK, ACTIVATION] * Srates[i, DECAY] * concentration[i];
            # IkB ----> nIkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
        end
        reactionFlux[IKBD, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBD];

        # IkB free in nucleus
        for i in NIKBA:NIKBD
            # nIkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
            # nIkB ----> IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 2: NFkB monomer Reactions
        #--------------------------------------------
        # NFkB mRNA
        # 0 --nNFkB--> tNFkB : nuclear NFkB dimer-induced synthesis
        # const rate * multiplier * ( 1 + sum( w*([d]/Kd)^hill)) / ( 1 + sum( ([d]/Kd)^hill))
        reactionFlux[TRELA, BASALSYNTHESIS] = Srates[TRELA, BASALSYNTHESIS];
        reactionFlux[TP50, BASALSYNTHESIS] = constituitiveMultiplier * Srates[TP50, BASALSYNTHESIS] * hillInduction(10, 10, concentration[NA50], concentration[NC50]; Kd = 150.0, hill = 1.0);
        reactionFlux[TCREL, BASALSYNTHESIS] = constituitiveMultiplier * Srates[TCREL, BASALSYNTHESIS] * hillInduction(10, 10, concentration[NA50], concentration[NC50]; Kd = 150.0, hill = 1.0);
        if phase == 1  # pre-simulation phase (no delay)
            for i in TRELA:TCREL
                reactionFlux[i, SYNTHESIS] = reactionFlux[i, BASALSYNTHESIS];
            end
        else
            reactionFlux[TRELA, SYNTHESIS] = reactionFlux[TRELA, BASALSYNTHESIS]; # no delay for tRelA
            if time > 1
                reactionFlux[TP50, SYNTHESIS] = constituitiveMultiplier * Srates[TP50, BASALSYNTHESIS] * hillInduction(10, 10, hist_NA50_1, hist_NC50_1; Kd = 150.0, hill = 1.0);
                reactionFlux[TCREL, SYNTHESIS] = constituitiveMultiplier * Srates[TCREL, BASALSYNTHESIS] * hillInduction(10, 10, hist_NA50_1, hist_NC50_1; Kd = 150.0, hill = 1.0);
            else # transcription is at basal level
                reactionFlux[TP50, SYNTHESIS] = constituitiveMultiplier * Srates[TP50, BASALSYNTHESIS];
                reactionFlux[TCREL, SYNTHESIS] = constituitiveMultiplier * Srates[TCREL, BASALSYNTHESIS];
            end
        end
        # Blimp1 ----| tcRel : Repression (multiplicative factor)
        reactionFlux[TCREL, REPRESSION] = hillRepression(Kr = 1.0, nr = 2.0, X = concentration[BLIMP1]);
        reactionFlux[TCREL, SYNTHESIS] = reactionFlux[TCREL, REPRESSION] * reactionFlux[TCREL, SYNTHESIS];
        for i in TRELA:TCREL
            # tNFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # NFkB free in cytoplasm
        for i in RELA:CREL
            # tNFkB ----> NFkB : translation
            reactionFlux[i, TRANSLATION] = EPS * Srates[i, TRANSLATION] * concentration[i-(RELA-TRELA)];
            # NFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        # NFkB free in nucleus
        for i in NRELA:NCREL
            # nNFkB ----> 0 : Basal decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 3: NFkB monomer <-> dimer Reactions
        #--------------------------------------------
        # NFkB + NFkB ----> NFkB:NFkB : association (cytoplasm)
        reactionFlux[AA, ASSOCIATION] = Srates[AA, ASSOCIATION] * concentration[RELA] ^ 2;
        reactionFlux[A50, ASSOCIATION] = Srates[A50, ASSOCIATION] * concentration[RELA] * concentration[P50];
        reactionFlux[P50P50, ASSOCIATION] = Srates[P50P50, ASSOCIATION] * concentration[P50] ^ 2;
        reactionFlux[C50, ASSOCIATION] = Srates[C50, ASSOCIATION] * concentration[CREL] * concentration[P50];
        # nNFkB + nNFkB ----> nNFkB:NFkB : association (nucleus)
        reactionFlux[NAA, ASSOCIATION] = Srates[NAA, ASSOCIATION] * concentration[NRELA] ^ 2;
        reactionFlux[NA50, ASSOCIATION] = Srates[NA50, ASSOCIATION] * concentration[NRELA] * concentration[NP50];
        reactionFlux[NP50P50, ASSOCIATION] = Srates[NP50P50, ASSOCIATION] * concentration[NP50] ^ 2;
        reactionFlux[NC50, ASSOCIATION] = Srates[NC50, ASSOCIATION] * concentration[NCREL] * concentration[NP50];

        for i in AA:C50
            # NFkB:NFkB ----> NFkB + NFkB : dissociation (cytoplasm)
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # NFkB:NFkB ----> nNFkB:NFkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # NFkB:NFkB ----> 0 : Basal decay (cytoplasm)
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        for i in NAA:NC50
            # nNFkB:NFkB ----> nNFkB + nNFkB : dissociation (nucleus)
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nNFkB:NFkB ----> NFkB:NFkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nNFkB:NFkB ----> 0 : Basal decay (nucleus)
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 4: AA:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAAA:IKBDAA
            # AA + IkB ----> AA:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[AA] * concentration[i-(IKBAAA-IKBA)];
            # AA:IkB ----> AA + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # AA:IkB ----> nAA:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # AA:IkB + IKK ----> AA + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
            reactionFlux[i, DECAY] = reactionFlux[IKK, ACTIVATION] * Srates[i-(IKBAAA-IKBA), DECAY] * concentration[i];
            # AA:IkB ----> AA : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] = reactionFlux[i, DECAY] + Srates[i, DECAY] * concentration[i];
            # AA:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        reactionFlux[IKBDAA, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDAA];

        # in nucleus
        for i in NIKBAAA:NIKBDAA
            # nAA + nIkB ----> nAA:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NAA] * concentration[i-(NIKBAAA-NIKBA)];
            # nAA:IkB ----> nAA + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nAA:IkB ----> AA:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nAA:IkB ----> nAA : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nAA:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 5: A50:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAA50:IKBDA50
            # A50 + IkB ----> A50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[A50] * concentration[i-(IKBAA50-IKBA)];
            # A50:IkB ----> A50 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # A50:IkB ----> nA50:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # A50:IkB + IKK ----> A50 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
            reactionFlux[i, DECAY] = reactionFlux[IKK, ACTIVATION] * Srates[i-(IKBAA50-IKBA), DECAY] * concentration[i];
            # A50:IkB ----> A50 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] = reactionFlux[i, DECAY] + Srates[i, DECAY] * concentration[i];
            # A50:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        reactionFlux[IKBDA50, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDA50];

        # in nucleus
        for i in NIKBAA50:NIKBDA50
            # nA50 + nIkB ----> nA50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NA50] * concentration[i-(NIKBAA50-NIKBA)];
            # nA50:IkB ----> nA50 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nA50:IkB ----> A50:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nA50:IkB ----> nA50 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nA50:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 6: 5050:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBA5050:IKBD5050
            # 5050 + IkB ----> 5050:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[P50P50] * concentration[i-(IKBA5050-IKBA)];
            # 5050:IkB ----> 5050 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # 5050:IkB ----> n5050:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # 5050:IkB + IKK ----> 5050 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
            reactionFlux[i, DECAY] = reactionFlux[IKK, ACTIVATION] * Srates[i-(IKBA5050-IKBA), DECAY] * concentration[i];
            # 5050:IkB ----> 5050 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] = reactionFlux[i, DECAY] + Srates[i, DECAY] * concentration[i];
            # 5050:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        reactionFlux[IKBD5050, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBD5050];

        # in nucleus
        for i in NIKBA5050:NIKBD5050
            # n5050 + nIkB ----> n5050:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NP50P50] * concentration[i-(NIKBA5050-NIKBA)];
            # n5050:IkB ----> n5050 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # n5050:IkB ----> 5050:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # n5050:IkB ----> n5050 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # n5050:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end

        #--------------------------------------------
        # MODULE 7: C50:IkB Reactions
        #--------------------------------------------
        # in cytoplasm
        for i in IKBAC50:IKBDC50
            # C50 + IkB ----> C50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[C50] * concentration[i-(IKBAC50-IKBA)];
            # C50:IkB ----> C50 + IkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # C50:IkB ----> nC50:IkB : transport into nucleus
            reactionFlux[i, TRANSPORTIN] = Srates[i, TRANSPORTIN] * concentration[i];
            # C50:IkB + IKK ----> C50 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
            reactionFlux[i, DECAY] = reactionFlux[IKK, ACTIVATION] * Srates[i-(IKBAC50-IKBA), DECAY] * concentration[i];
            # C50:IkB ----> C50 : IKK-independent NFkB-bound IkB decay
            reactionFlux[i, DECAY] = reactionFlux[i, DECAY] + Srates[i, DECAY] * concentration[i];
            # C50:IkB ----> IkB : IkB-bound NFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
        reactionFlux[IKBDC50, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDC50];

        # in nucleus
        for i in NIKBAC50:NIKBDC50
            # nC50 + nIkB ----> nC50:IkB : association
            reactionFlux[i, ASSOCIATION] = Srates[i, ASSOCIATION] * concentration[NC50] * concentration[i-(NIKBAC50-NIKBA)];
            # nC50:IkB ----> nC50 + nIkB : dissociation
            reactionFlux[i, DISSOCIATION] = Srates[i, DISSOCIATION] * concentration[i];
            # nC50:IkB ----> C50:IkB : transport out of nucleus
            reactionFlux[i, TRANSPORTOUT] = Srates[i, TRANSPORTOUT] * concentration[i];
            # nC50:IkB ----> nC50 : IKK-independent nNFkB-bound nIkB decay
            reactionFlux[i, DECAY] = Srates[i, DECAY] * concentration[i];
            # nC50:IkB ----> nIkB : nIkB-bound nNFkB decay
            reactionFlux[i, BASALDECAY] = Srates[i, BASALDECAY] * concentration[i];
        end
    end
    nothing
end

function NFkBNettFluxes!(nettFlux, reactionFlux)
    @inbounds begin
        #------------------------------------------------
        # MODULE 1
        #--------------------------------------------
        # 2-5 : IkB mRNA
        for i in TIKBA:TIKBD
            nettFlux[i] = reactionFlux[i, SYNTHESIS] - reactionFlux[i, BASALDECAY];
        end
        # 6-9 : IkB free in cytoplasm
        for i in IKBA:IKBD
            nettFlux[i] = reactionFlux[i, TRANSLATION] - reactionFlux[i, BASALDECAY] - reactionFlux[i, DECAY] -
            reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBA+NIKBA, TRANSPORTOUT] -
            reactionFlux[i-IKBA+IKBAAA, ASSOCIATION] + reactionFlux[i-IKBA+IKBAAA, DISSOCIATION] + reactionFlux[i-IKBA+IKBAAA, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAA50, ASSOCIATION] + reactionFlux[i-IKBA+IKBAA50, DISSOCIATION] + reactionFlux[i-IKBA+IKBAA50, BASALDECAY] -
            reactionFlux[i-IKBA+IKBA5050, ASSOCIATION] + reactionFlux[i-IKBA+IKBA5050, DISSOCIATION] + reactionFlux[i-IKBA+IKBA5050, BASALDECAY] -
            reactionFlux[i-IKBA+IKBAC50, ASSOCIATION] + reactionFlux[i-IKBA+IKBAC50, DISSOCIATION] + reactionFlux[i-IKBA+IKBAC50, BASALDECAY];
        end
        # 10-13 : IkB free in nucleus
        for i in NIKBA:NIKBD
            nettFlux[i] = - reactionFlux[i, BASALDECAY] +
            reactionFlux[i-NIKBA+IKBA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] -
            reactionFlux[i-NIKBA+NIKBAAA, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAAA, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAAA, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAA50, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAA50, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAA50, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBA5050, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBA5050, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBA5050, BASALDECAY] -
            reactionFlux[i-NIKBA+NIKBAC50, ASSOCIATION] + reactionFlux[i-NIKBA+NIKBAC50, DISSOCIATION] + reactionFlux[i-NIKBA+NIKBAC50, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 2
        #--------------------------------------------
        # NFkB mRNA
        for i in TRELA:TCREL
            nettFlux[i] = reactionFlux[i, SYNTHESIS] - reactionFlux[i, BASALDECAY];
        end
        # NFkB free in cytoplasm
        for i in RELA:CREL
            nettFlux[i] = reactionFlux[i, TRANSLATION] - reactionFlux[i, BASALDECAY];
        end
        nettFlux[RELA] += - 2 * reactionFlux[AA, ASSOCIATION] + 2 * reactionFlux[AA, DISSOCIATION];
        nettFlux[RELA] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[P50] += - 2 * reactionFlux[P50P50, ASSOCIATION] + 2 * reactionFlux[P50P50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        nettFlux[CREL] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        # NFkB free in nucleus
        for i in NRELA:NCREL
            nettFlux[i] = - reactionFlux[i, BASALDECAY];
        end
        nettFlux[NRELA] += - 2 * reactionFlux[NAA, ASSOCIATION] + 2 * reactionFlux[NAA, DISSOCIATION];
        nettFlux[NRELA] += - reactionFlux[NA50, ASSOCIATION] + reactionFlux[NA50, DISSOCIATION];
        nettFlux[NP50] += - 2 * reactionFlux[NP50P50, ASSOCIATION] + 2 * reactionFlux[NP50P50, DISSOCIATION];
        nettFlux[NP50] += - reactionFlux[NA50, ASSOCIATION] + reactionFlux[NA50, DISSOCIATION];
        nettFlux[NP50] += - reactionFlux[NC50, ASSOCIATION] + reactionFlux[NC50, DISSOCIATION];
        nettFlux[NCREL] += - reactionFlux[NC50, ASSOCIATION] + reactionFlux[NC50, DISSOCIATION];

        #------------------------------------------------
        # MODULE 3
        #--------------------------------------------
        # NFkB free in cytoplasm
        for i in AA:C50
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-AA+NAA, TRANSPORTOUT] - reactionFlux[i, BASALDECAY];
        end
        # NFkB free in nucleus
        for i in NAA:NC50
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NAA+AA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 4
        #--------------------------------------------
        for i in IKBAAA:IKBDAA
            # AA free in cytoplasm
            nettFlux[AA] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # AA:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAAA+NIKBAAA, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAAA:NIKBDAA
            # AA free in nucleus
            nettFlux[NAA] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # AA:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAAA+IKBAAA, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 5
        #--------------------------------------------
        for i in IKBAA50:IKBDA50
            # A50 free in cytoplasm
            nettFlux[A50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A50:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAA50+NIKBAA50, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAA50:NIKBDA50
            # A50 free in nucleus
            nettFlux[NA50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # A50:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAA50+IKBAA50, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 6
        #--------------------------------------------
        for i in IKBA5050:IKBD5050
            # 5050 free in cytoplasm
            nettFlux[P50P50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # 5050:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBA5050+NIKBA5050, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBA5050:NIKBD5050
            # 5050 free in nucleus
            nettFlux[NP50P50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # 5050:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBA5050+IKBA5050, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        #------------------------------------------------
        # MODULE 7
        #--------------------------------------------
        for i in IKBAC50:IKBDC50
            # C50 free in cytoplasm
            nettFlux[C50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C50:IkB in cytoplasm
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] - reactionFlux[i, TRANSPORTIN] + reactionFlux[i-IKBAC50+NIKBAC50, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

        for i in NIKBAC50:NIKBDC50
            # C50 free in nucleus
            nettFlux[NC50] += - reactionFlux[i, ASSOCIATION] + reactionFlux[i, DISSOCIATION] + reactionFlux[i, DECAY];
            # C50:IkB in nucleus
            nettFlux[i] = reactionFlux[i, ASSOCIATION] - reactionFlux[i, DISSOCIATION] + reactionFlux[i-NIKBAC50+IKBAC50, TRANSPORTIN] - reactionFlux[i, TRANSPORTOUT] - reactionFlux[i, DECAY] - reactionFlux[i, BASALDECAY];
        end

    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeNFkBFluxes!(concentration, reactionFlux, Srates, phase);
    # Compute NFkB net fluxes
    NFkBNettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNFkBNettFluxes!(nettFlux, concentration, delay, reactionFlux, Srates, phase, time, birthday, IKKCurve, historicFlux)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeNFkBFluxes!(concentration, reactionFlux, Srates, phase; delay=delay, time=time, birthday=birthday, IKKCurve=IKKCurve, historicFlux=historicFlux);
    # Compute NFkB net fluxes
    NFkBNettFluxes!(nettFlux, reactionFlux);
    nothing
end
