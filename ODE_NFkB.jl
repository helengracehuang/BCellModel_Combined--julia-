# include("HelperFunctions.jl");
# Compute reaction fluxes for NFkB module
#--------------------------------------------

function computeNFkBFluxes!(concentration, reactionFlux, Srates, phase; delay=nothing, time=nothing, birthday=nothing, IKKCurve=nothing, historicFlux=nothing)
    @inbounds begin
        # get historic values for nuclear NFkB concentrations (for delayed transcription)
        if phase == 2
            p = (birthday, IKKCurve, Srates, reactionFlux, phase, historicFlux);
            hist_NAA_075 = delay(p, time-0.75/CONVERSION; idxs=NAA);
            hist_NA50_075 = delay(p, time-0.75/CONVERSION; idxs=NA50);
            hist_NC50_075 = delay(p, time-0.75/CONVERSION; idxs=NC50);

            hist_NA50_1 = delay(p, time-1.0/CONVERSION; idxs=NA50);
            hist_NC50_1 = delay(p, time-1.0/CONVERSION; idxs=NC50);

            hist_NAA_15 = delay(p, time-1.5/CONVERSION; idxs=NAA);
            hist_NA50_15 = delay(p, time-1.5/CONVERSION; idxs=NA50);
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
            @. reactionFlux[TIKBA:TIKBD, SYNTHESIS] = (@view reactionFlux[TIKBA:TIKBD, BASALSYNTHESIS]);
        else
            @. reactionFlux[TIKBA:TIKBB, SYNTHESIS] = (@view reactionFlux[TIKBA:TIKBB, BASALSYNTHESIS]); # no delay for tIkba/b
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
        @. reactionFlux[TIKBA:TIKBD, BASALDECAY] = (@view Srates[TIKBA:TIKBD, BASALDECAY]) * (@view concentration[TIKBA:TIKBD]);

        # IkB free in cytoplasm
        # tIkB ----> IkB : translation
        @. reactionFlux[IKBA:IKBD, TRANSLATION] = EPS * (@view Srates[IKBA:IKBD, TRANSLATION]) * (@view concentration[TIKBA:TIKBD]);
        # IkB ----> 0 : Basal decay
        @. reactionFlux[IKBA:IKBD, BASALDECAY] = (@view Srates[IKBA:IKBD, BASALDECAY]) * (@view concentration[IKBA:IKBD]);
        # IkB + IKK ----> IKK : IKK-mediated free IkB decay (basal IKK only for IkBd)
        @. reactionFlux[IKBA:IKBE, DECAY] = reactionFlux[IKK, ACTIVATION] * (@view Srates[IKBA:IKBE, DECAY]) * (@view concentration[IKBA:IKBE]);
        reactionFlux[IKBD, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBD];
        # IkB ----> nIkB : transport into nucleus
        @. reactionFlux[IKBA:IKBD, TRANSPORTIN] = (@view Srates[IKBA:IKBD, TRANSPORTIN]) * (@view concentration[IKBA:IKBD]);

        # IkB free in nucleus
        # nIkB ----> 0 : Basal decay
        @. reactionFlux[NIKBA:NIKBD, BASALDECAY] = (@view Srates[NIKBA:NIKBD, BASALDECAY]) * (@view concentration[NIKBA:NIKBD]);
        # nIkB ----> IkB : transport out of nucleus
        @. reactionFlux[NIKBA:NIKBD, TRANSPORTOUT] = (@view Srates[NIKBA:NIKBD, TRANSPORTOUT]) * (@view concentration[NIKBA:NIKBD]);

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
            @. reactionFlux[TRELA:TCREL, SYNTHESIS] = (@view reactionFlux[TRELA:TCREL, BASALSYNTHESIS]);
        else
            reactionFlux[TRELA, SYNTHESIS] = reactionFlux[TRELA, BASALSYNTHESIS]; # no delay for tRelA
            if time > 1
                reactionFlux[TP50, SYNTHESIS] = constituitiveMultiplier * Srates[TP50, BASALSYNTHESIS] * hillInduction(10, 10, hist_NA50_1, hist_NC50_1; Kd = 150.0, hill = 1.0);
                reactionFlux[TCREL, SYNTHESIS] = constituitiveMultiplier * Srates[TCREL, BASALSYNTHESIS] * hillInduction(10, 10, hist_NA50_1, hist_NC50_1; Kd = 150.0, hill = 1.0);
            else # transcription is at basal level
                @. reactionFlux[TP50:TCREL, SYNTHESIS] = constituitiveMultiplier * (@view Srates[TP50:TCREL, BASALSYNTHESIS]);
            end
        end
        # Blimp1 ----| tcRel : Repression (multiplicative factor)
        reactionFlux[TCREL, REPRESSION] = hillRepression(Kr = 1.0, nr = 2.0, X = concentration[BLIMP1]);
        reactionFlux[TCREL, SYNTHESIS] = reactionFlux[TCREL, REPRESSION] * reactionFlux[TCREL, SYNTHESIS];
        # tNFkB ----> 0 : Basal decay
        @. reactionFlux[TRELA:TCREL, BASALDECAY] = (@view Srates[TRELA:TCREL, BASALDECAY]) * (@view concentration[TRELA:TCREL]);

        # NFkB free in cytoplasm
        # tNFkB ----> NFkB : translation
        @. reactionFlux[RELA:CREL, TRANSLATION] = EPS * (@view Srates[RELA:CREL, TRANSLATION]) * (@view concentration[TRELA:TCREL]);
        # NFkB ----> 0 : Basal decay
        @. reactionFlux[RELA:CREL, BASALDECAY] = (@view Srates[RELA:CREL, BASALDECAY]) * (@view concentration[RELA:CREL]);

        # NFkB free in nucleus
        # nNFkB ----> 0 : Basal decay
        @. reactionFlux[NRELA:NCREL, BASALDECAY] = (@view Srates[NRELA:NCREL, BASALDECAY]) * (@view concentration[NRELA:NCREL]);

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

        # NFkB:NFkB ----> NFkB + NFkB : dissociation (cytoplasm)
        @. reactionFlux[AA:C50, DISSOCIATION] = (@view Srates[AA:C50, DISSOCIATION]) * (@view concentration[AA:C50]);
        # nNFkB:NFkB ----> nNFkB + nNFkB : dissociation (nucleus)
        @. reactionFlux[NAA:NC50, DISSOCIATION] = (@view Srates[NAA:NC50, DISSOCIATION]) * (@view concentration[NAA:NC50]);

        # NFkB:NFkB ----> nNFkB:NFkB : transport into nucleus
        @. reactionFlux[AA:C50, TRANSPORTIN] = (@view Srates[AA:C50, TRANSPORTIN]) * (@view concentration[AA:C50]);
        # nNFkB:NFkB ----> NFkB:NFkB : transport out of nucleus
        @. reactionFlux[NAA:NC50, TRANSPORTOUT] = (@view Srates[NAA:NC50, TRANSPORTOUT]) * (@view concentration[NAA:NC50]);

        # NFkB:NFkB ----> 0 : Basal decay (cytoplasm)
        @. reactionFlux[AA:C50, BASALDECAY] = (@view Srates[AA:C50, BASALDECAY]) * (@view concentration[AA:C50]);
        # nNFkB:NFkB ----> 0 : Basal decay (nucleus)
        @. reactionFlux[NAA:NC50, BASALDECAY] = (@view Srates[NAA:NC50, BASALDECAY]) * (@view concentration[NAA:NC50]);

        #--------------------------------------------
        # MODULE 4: AA:IkB Reactions
        #--------------------------------------------
        # AA + IkB ----> AA:IkB : association (cytoplasm)
        @. reactionFlux[IKBAAA:IKBDAA, ASSOCIATION] = (@view Srates[IKBAAA:IKBDAA, ASSOCIATION]) * concentration[AA] * (@view concentration[IKBA:IKBD]);
        # nAA + nIkB ----> nAA:IkB : association (nucleus)
        @. reactionFlux[NIKBAAA:NIKBDAA, ASSOCIATION] = (@view Srates[NIKBAAA:NIKBDAA, ASSOCIATION]) * concentration[NAA] * (@view concentration[NIKBA:NIKBD]);

        # AA:IkB ----> AA + IkB : dissociation (cytoplasm)
        @. reactionFlux[IKBAAA:IKBDAA, DISSOCIATION] = (@view Srates[IKBAAA:IKBDAA, DISSOCIATION]) * (@view concentration[IKBAAA:IKBDAA]);
        # nAA:IkB ----> nAA + nIkB : dissociation (nucleus)
        @. reactionFlux[NIKBAAA:NIKBDAA, DISSOCIATION] = (@view Srates[NIKBAAA:NIKBDAA, DISSOCIATION]) * (@view concentration[NIKBAAA:NIKBDAA]);

        # AA:IkB ----> nAA:IkB : transport into nucleus
        @. reactionFlux[IKBAAA:IKBDAA, TRANSPORTIN] = (@view Srates[IKBAAA:IKBDAA, TRANSPORTIN]) * (@view concentration[IKBAAA:IKBDAA]);
        # nAA:IkB ----> AA:IkB : transport out of nucleus
        @. reactionFlux[NIKBAAA:NIKBDAA, TRANSPORTOUT] = (@view Srates[NIKBAAA:NIKBDAA, TRANSPORTOUT]) * (@view concentration[NIKBAAA:NIKBDAA]);

        # AA:IkB + IKK ----> AA + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
        @. reactionFlux[IKBAAA:IKBEAA, DECAY] = reactionFlux[IKK, ACTIVATION] * (@view Srates[IKBA:IKBE, DECAY]) * (@view concentration[IKBAAA:IKBEAA]);
        reactionFlux[IKBDAA, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDAA];

        # AA:IkB ----> AA : IKK-independent NFkB-bound IkB decay
        @. reactionFlux[IKBAAA:IKBDAA, DECAY] = (@view reactionFlux[IKBAAA:IKBDAA, DECAY]) + (@view Srates[IKBAAA:IKBDAA, DECAY]) * (@view concentration[IKBAAA:IKBDAA]);
        # nAA:IkB ----> nAA : IKK-independent nNFkB-bound nIkB decay
        @. reactionFlux[NIKBAAA:NIKBDAA, DECAY] = (@view Srates[NIKBAAA:NIKBDAA, DECAY]) * (@view concentration[NIKBAAA:NIKBDAA]);

        # AA:IkB ----> IkB : IkB-bound NFkB decay
        @. reactionFlux[IKBAAA:IKBDAA, BASALDECAY] = (@view Srates[IKBAAA:IKBDAA, BASALDECAY]) * (@view concentration[IKBAAA:IKBDAA]);
        # nAA:IkB ----> nIkB : nIkB-bound nNFkB decay
        @. reactionFlux[NIKBAAA:NIKBDAA, BASALDECAY] = (@view Srates[NIKBAAA:NIKBDAA, BASALDECAY]) * (@view concentration[NIKBAAA:NIKBDAA]);

        #--------------------------------------------
        # MODULE 5: A50:IkB Reactions
        #--------------------------------------------
        # A50 + IkB ----> A50:IkB : association (cytoplasm)
        @. reactionFlux[IKBAA50:IKBDA50, ASSOCIATION] = (@view Srates[IKBAA50:IKBDA50, ASSOCIATION]) * concentration[A50] * (@view concentration[IKBA:IKBD]);
        # nA50 + nIkB ----> nA50:IkB : association (nucleus)
        @. reactionFlux[NIKBAA50:NIKBDA50, ASSOCIATION] = (@view Srates[NIKBAA50:NIKBDA50, ASSOCIATION]) * concentration[NA50] * (@view concentration[NIKBA:NIKBD]);

        # A50:IkB ----> A50 + IkB : dissociation (cytoplasm)
        @. reactionFlux[IKBAA50:IKBDA50, DISSOCIATION] = (@view Srates[IKBAA50:IKBDA50, DISSOCIATION]) * (@view concentration[IKBAA50:IKBDA50]);
        # nA50:IkB ----> nA50 + nIkB : dissociation (nucleus)
        @. reactionFlux[NIKBAA50:NIKBDA50, DISSOCIATION] = (@view Srates[NIKBAA50:NIKBDA50, DISSOCIATION]) * (@view concentration[NIKBAA50:NIKBDA50]);

        # A50:IkB ----> nA50:IkB : transport into nucleus
        @. reactionFlux[IKBAA50:IKBDA50, TRANSPORTIN] = (@view Srates[IKBAA50:IKBDA50, TRANSPORTIN]) * (@view concentration[IKBAA50:IKBDA50]);
        # nA50:IkB ----> A50:IkB : transport out of nucleus
        @. reactionFlux[NIKBAA50:NIKBDA50, TRANSPORTOUT] = (@view Srates[NIKBAA50:NIKBDA50, TRANSPORTOUT]) * (@view concentration[NIKBAA50:NIKBDA50]);

        # A50:IkB + IKK ----> A50 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
        @. reactionFlux[IKBAA50:IKBEA50, DECAY] = reactionFlux[IKK, ACTIVATION] * (@view Srates[IKBA:IKBE, DECAY]) * (@view concentration[IKBAA50:IKBEA50]);
        reactionFlux[IKBDA50, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDA50];

        # A50:IkB ----> A50 : IKK-independent NFkB-bound IkB decay
        @. reactionFlux[IKBAA50:IKBDA50, DECAY] = (@view reactionFlux[IKBAA50:IKBDA50, DECAY]) + (@view Srates[IKBAA50:IKBDA50, DECAY]) * (@view concentration[IKBAA50:IKBDA50]);
        # nA50:IkB ----> nA50 : IKK-independent nNFkB-bound nIkB decay
        @. reactionFlux[NIKBAA50:NIKBDA50, DECAY] = (@view Srates[NIKBAA50:NIKBDA50, DECAY]) * (@view concentration[NIKBAA50:NIKBDA50]);

        # A50:IkB ----> IkB : IkB-bound NFkB decay
        @. reactionFlux[IKBAA50:IKBDA50, BASALDECAY] = (@view Srates[IKBAA50:IKBDA50, BASALDECAY]) * (@view concentration[IKBAA50:IKBDA50]);
        # nA50:IkB ----> nIkB : nIkB-bound nNFkB decay
        @. reactionFlux[NIKBAA50:NIKBDA50, BASALDECAY] = (@view Srates[NIKBAA50:NIKBDA50, BASALDECAY]) * (@view concentration[NIKBAA50:NIKBDA50]);

        #--------------------------------------------
        # MODULE 6: 5050:IkB Reactions
        #--------------------------------------------
        # 5050 + IkB ----> 5050:IkB : association (cytoplasm)
        @. reactionFlux[IKBA5050:IKBD5050, ASSOCIATION] = (@view Srates[IKBA5050:IKBD5050, ASSOCIATION]) * concentration[P50P50] * (@view concentration[IKBA:IKBD]);
        # n5050 + nIkB ----> n5050:IkB : association (nucleus)
        @. reactionFlux[NIKBA5050:NIKBD5050, ASSOCIATION] = (@view Srates[NIKBA5050:NIKBD5050, ASSOCIATION]) * concentration[NP50P50] * (@view concentration[NIKBA:NIKBD]);

        # 5050:IkB ----> 5050 + IkB : dissociation (cytoplasm)
        @. reactionFlux[IKBA5050:IKBD5050, DISSOCIATION] = (@view Srates[IKBA5050:IKBD5050, DISSOCIATION]) * (@view concentration[IKBA5050:IKBD5050]);
        # n5050:IkB ----> n5050 + nIkB : dissociation (nucleus)
        @. reactionFlux[NIKBA5050:NIKBD5050, DISSOCIATION] = (@view Srates[NIKBA5050:NIKBD5050, DISSOCIATION]) * (@view concentration[NIKBA5050:NIKBD5050]);

        # 5050:IkB ----> n5050:IkB : transport into nucleus
        @. reactionFlux[IKBA5050:IKBD5050, TRANSPORTIN] = (@view Srates[IKBA5050:IKBD5050, TRANSPORTIN]) * (@view concentration[IKBA5050:IKBD5050]);
        # n5050:IkB ----> 5050:IkB : transport out of nucleus
        @. reactionFlux[NIKBA5050:NIKBD5050, TRANSPORTOUT] = (@view Srates[NIKBA5050:NIKBD5050, TRANSPORTOUT]) * (@view concentration[NIKBA5050:NIKBD5050]);

        # 5050:IkB + IKK ----> 5050 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
        @. reactionFlux[IKBA5050:IKBE5050, DECAY] = reactionFlux[IKK, ACTIVATION] * (@view Srates[IKBA:IKBE, DECAY]) * (@view concentration[IKBA5050:IKBE5050]);
        reactionFlux[IKBD5050, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBD5050];

        # 5050:IkB ----> 5050 : IKK-independent NFkB-bound IkB decay
        @. reactionFlux[IKBA5050:IKBD5050, DECAY] = (@view reactionFlux[IKBA5050:IKBD5050, DECAY]) + (@view Srates[IKBA5050:IKBD5050, DECAY]) * (@view concentration[IKBA5050:IKBD5050]);
        # n5050:IkB ----> n5050 : IKK-independent nNFkB-bound nIkB decay
        @. reactionFlux[NIKBA5050:NIKBD5050, DECAY] = (@view Srates[NIKBA5050:NIKBD5050, DECAY]) * (@view concentration[NIKBA5050:NIKBD5050]);

        # 5050:IkB ----> IkB : IkB-bound NFkB decay
        @. reactionFlux[IKBA5050:IKBD5050, BASALDECAY] = (@view Srates[IKBA5050:IKBD5050, BASALDECAY]) * (@view concentration[IKBA5050:IKBD5050]);
        # n5050:IkB ----> nIkB : nIkB-bound nNFkB decay
        @. reactionFlux[NIKBA5050:NIKBD5050, BASALDECAY] = (@view Srates[NIKBA5050:NIKBD5050, BASALDECAY]) * (@view concentration[NIKBA5050:NIKBD5050]);

        #--------------------------------------------
        # MODULE 7: C50:IkB Reactions
        #--------------------------------------------
        # C50 + IkB ----> C50:IkB : association (cytoplasm)
        @. reactionFlux[IKBAC50:IKBDC50, ASSOCIATION] = (@view Srates[IKBAC50:IKBDC50, ASSOCIATION]) * concentration[C50] * (@view concentration[IKBA:IKBD]);
        # nC50 + nIkB ----> nC50:IkB : association (nucleus)
        @. reactionFlux[NIKBAC50:NIKBDC50, ASSOCIATION] = (@view Srates[NIKBAC50:NIKBDC50, ASSOCIATION]) * concentration[NC50] * (@view concentration[NIKBA:NIKBD]);

        # C50:IkB ----> C50 + IkB : dissociation (cytoplasm)
        @. reactionFlux[IKBAC50:IKBDC50, DISSOCIATION] = (@view Srates[IKBAC50:IKBDC50, DISSOCIATION]) * (@view concentration[IKBAC50:IKBDC50]);
        # nC50:IkB ----> nC50 + nIkB : dissociation (nucleus)
        @. reactionFlux[NIKBAC50:NIKBDC50, DISSOCIATION] = (@view Srates[NIKBAC50:NIKBDC50, DISSOCIATION]) * (@view concentration[NIKBAC50:NIKBDC50]);

        # C50:IkB ----> nC50:IkB : transport into nucleus
        @. reactionFlux[IKBAC50:IKBDC50, TRANSPORTIN] = (@view Srates[IKBAC50:IKBDC50, TRANSPORTIN]) * (@view concentration[IKBAC50:IKBDC50]);
        # nC50:IkB ----> C50:IkB : transport out of nucleus
        @. reactionFlux[NIKBAC50:NIKBDC50, TRANSPORTOUT] = (@view Srates[NIKBAC50:NIKBDC50, TRANSPORTOUT]) * (@view concentration[NIKBAC50:NIKBDC50]);

        # C50:IkB + IKK ----> C50 + IKK : IKK-mediated NFkB-bound IkB decay (basal IKK only for IkBd)
        @. reactionFlux[IKBAC50:IKBEC50, DECAY] = reactionFlux[IKK, ACTIVATION] * (@view Srates[IKBA:IKBE, DECAY]) * (@view concentration[IKBAC50:IKBEC50]);
        reactionFlux[IKBDC50, DECAY] = reactionFlux[IKK, DEACTIVATION] * Srates[IKBD, DECAY] * concentration[IKBDC50];

        # C50:IkB ----> C50 : IKK-independent NFkB-bound IkB decay
        @. reactionFlux[IKBAC50:IKBDC50, DECAY] = (@view reactionFlux[IKBAC50:IKBDC50, DECAY]) + (@view Srates[IKBAC50:IKBDC50, DECAY]) * (@view concentration[IKBAC50:IKBDC50]);
        # nC50:IkB ----> nC50 : IKK-independent nNFkB-bound nIkB decay
        @. reactionFlux[NIKBAC50:NIKBDC50, DECAY] = (@view Srates[NIKBAC50:NIKBDC50, DECAY]) * (@view concentration[NIKBAC50:NIKBDC50]);

        # C50:IkB ----> IkB : IkB-bound NFkB decay
        @. reactionFlux[IKBAC50:IKBDC50, BASALDECAY] = (@view Srates[IKBAC50:IKBDC50, BASALDECAY]) * (@view concentration[IKBAC50:IKBDC50]);
        # nC50:IkB ----> nIkB : nIkB-bound nNFkB decay
        @. reactionFlux[NIKBAC50:NIKBDC50, BASALDECAY] = (@view Srates[NIKBAC50:NIKBDC50, BASALDECAY]) * (@view concentration[NIKBAC50:NIKBDC50]);
    end
    nothing
end

function NFkBNettFluxes!(nettFlux, reactionFlux)
    @inbounds begin
        #------------------------------------------------
        # MODULE 1
        #--------------------------------------------
        # 2-5 : IkB mRNA
        @. nettFlux[TIKBA:TIKBD] = (@view reactionFlux[TIKBA:TIKBD, SYNTHESIS]) - (@view reactionFlux[TIKBA:TIKBD, BASALDECAY]);
        # 6-9 : IkB free in cytoplasm
        @. nettFlux[IKBA:IKBD] = (@view reactionFlux[IKBA:IKBD, TRANSLATION]) - (@view reactionFlux[IKBA:IKBD, BASALDECAY]) - (@view reactionFlux[IKBA:IKBD, DECAY]) - (@view reactionFlux[IKBA:IKBD, TRANSPORTIN]) + (@view reactionFlux[NIKBA:NIKBD, TRANSPORTOUT]) - (@view reactionFlux[IKBAAA:IKBDAA, ASSOCIATION]) + (@view reactionFlux[IKBAAA:IKBDAA, DISSOCIATION]) + (@view reactionFlux[IKBAAA:IKBDAA, BASALDECAY]) - (@view reactionFlux[IKBAA50:IKBDA50, ASSOCIATION]) + (@view reactionFlux[IKBAA50:IKBDA50, DISSOCIATION]) + (@view reactionFlux[IKBAA50:IKBDA50, BASALDECAY]) - (@view reactionFlux[IKBA5050:IKBD5050, ASSOCIATION]) + (@view reactionFlux[IKBA5050:IKBD5050, DISSOCIATION]) + (@view reactionFlux[IKBA5050:IKBD5050, BASALDECAY]) - (@view reactionFlux[IKBAC50:IKBDC50, ASSOCIATION]) + (@view reactionFlux[IKBAC50:IKBDC50, DISSOCIATION]) + (@view reactionFlux[IKBAC50:IKBDC50, BASALDECAY]);
        # 10-13 : IkB free in nucleus
        @. nettFlux[NIKBA:NIKBD] = - (@view reactionFlux[NIKBA:NIKBD, BASALDECAY]) + (@view reactionFlux[IKBA:IKBD, TRANSPORTIN]) - (@view reactionFlux[NIKBA:NIKBD, TRANSPORTOUT]) - (@view reactionFlux[NIKBAAA:NIKBDAA, ASSOCIATION]) + (@view reactionFlux[NIKBAAA:NIKBDAA, DISSOCIATION]) + (@view reactionFlux[NIKBAAA:NIKBDAA, BASALDECAY]) - (@view reactionFlux[NIKBAA50:NIKBDA50, ASSOCIATION]) + (@view reactionFlux[NIKBAA50:NIKBDA50, DISSOCIATION]) + (@view reactionFlux[NIKBAA50:NIKBDA50, BASALDECAY]) - (@view reactionFlux[NIKBA5050:NIKBD5050, ASSOCIATION]) + (@view reactionFlux[NIKBA5050:NIKBD5050, DISSOCIATION]) + (@view reactionFlux[NIKBA5050:NIKBD5050, BASALDECAY]) - (@view reactionFlux[NIKBAC50:NIKBDC50, ASSOCIATION]) + (@view reactionFlux[NIKBAC50:NIKBDC50, DISSOCIATION]) + (@view reactionFlux[NIKBAC50:NIKBDC50, BASALDECAY]);

        #------------------------------------------------
        # MODULE 2
        #--------------------------------------------
        # NFkB mRNA
        @. nettFlux[TRELA:TCREL] = (@view reactionFlux[TRELA:TCREL, SYNTHESIS]) - (@view reactionFlux[TRELA:TCREL, BASALDECAY]);
        # NFkB free in cytoplasm
        @. nettFlux[RELA:CREL] = (@view reactionFlux[RELA:CREL, TRANSLATION]) - (@view reactionFlux[RELA:CREL, BASALDECAY]);
        nettFlux[RELA] += - 2 * reactionFlux[AA, ASSOCIATION] + 2 * reactionFlux[AA, DISSOCIATION];
        nettFlux[RELA] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[P50] += - 2 * reactionFlux[P50P50, ASSOCIATION] + 2 * reactionFlux[P50P50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[A50, ASSOCIATION] + reactionFlux[A50, DISSOCIATION];
        nettFlux[P50] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        nettFlux[CREL] += - reactionFlux[C50, ASSOCIATION] + reactionFlux[C50, DISSOCIATION];
        # NFkB free in nucleus
        @. nettFlux[NRELA:NCREL] = - (@view reactionFlux[NRELA:NCREL, BASALDECAY]);
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
        @. nettFlux[AA:C50] = (@view reactionFlux[AA:C50, ASSOCIATION]) - (@view reactionFlux[AA:C50, DISSOCIATION]) - (@view reactionFlux[AA:C50, TRANSPORTIN]) + (@view reactionFlux[NAA:NC50, TRANSPORTOUT]) - (@view reactionFlux[AA:C50, BASALDECAY]);
        # NFkB free in nucleus
        @. nettFlux[NAA:NC50] = (@view reactionFlux[NAA:NC50, ASSOCIATION]) - (@view reactionFlux[NAA:NC50, DISSOCIATION]) + (@view reactionFlux[AA:C50, TRANSPORTIN]) - (@view reactionFlux[NAA:NC50, TRANSPORTOUT]) - (@view reactionFlux[NAA:NC50, BASALDECAY]);

        #------------------------------------------------
        # MODULE 4
        #--------------------------------------------
        # AA free in cytoplasm
        nettFlux[AA] += - sum(@view reactionFlux[IKBAAA:IKBDAA, ASSOCIATION]) + sum(@view reactionFlux[IKBAAA:IKBDAA, DISSOCIATION]) + sum(@view reactionFlux[IKBAAA:IKBDAA, DECAY]);
        # AA free in nucleus
        nettFlux[NAA] += - sum(@view reactionFlux[NIKBAAA:NIKBDAA, ASSOCIATION]) + sum(@view reactionFlux[NIKBAAA:NIKBDAA, DISSOCIATION]) + sum(@view reactionFlux[NIKBAAA:NIKBDAA, DECAY]);
        # AA:IkB in cytoplasm
        @. nettFlux[IKBAAA:IKBDAA] = (@view reactionFlux[IKBAAA:IKBDAA, ASSOCIATION]) - (@view reactionFlux[IKBAAA:IKBDAA, DISSOCIATION]) - (@view reactionFlux[IKBAAA:IKBDAA, TRANSPORTIN]) + (@view reactionFlux[NIKBAAA:NIKBDAA, TRANSPORTOUT]) - (@view reactionFlux[IKBAAA:IKBDAA, DECAY]) - (@view reactionFlux[IKBAAA:IKBDAA, BASALDECAY]);
        # AA:IkB in nucleus
        @. nettFlux[NIKBAAA:NIKBDAA] = (@view reactionFlux[NIKBAAA:NIKBDAA, ASSOCIATION]) - (@view reactionFlux[NIKBAAA:NIKBDAA, DISSOCIATION]) + (@view reactionFlux[IKBAAA:IKBDAA, TRANSPORTIN]) - (@view reactionFlux[NIKBAAA:NIKBDAA, TRANSPORTOUT]) - (@view reactionFlux[NIKBAAA:NIKBDAA, DECAY]) - (@view reactionFlux[NIKBAAA:NIKBDAA, BASALDECAY]);

        #------------------------------------------------
        # MODULE 5
        #--------------------------------------------
        # A50 free in cytoplasm
        nettFlux[A50] += - sum(@view reactionFlux[IKBAA50:IKBDA50, ASSOCIATION]) + sum(@view reactionFlux[IKBAA50:IKBDA50, DISSOCIATION]) + sum(@view reactionFlux[IKBAA50:IKBDA50, DECAY]);
        # A50 free in nucleus
        nettFlux[NA50] += - sum(@view reactionFlux[NIKBAA50:NIKBDA50, ASSOCIATION]) + sum(@view reactionFlux[NIKBAA50:NIKBDA50, DISSOCIATION]) + sum(@view reactionFlux[NIKBAA50:NIKBDA50, DECAY]);
        # A50:IkB in cytoplasm
        @. nettFlux[IKBAA50:IKBDA50] = (@view reactionFlux[IKBAA50:IKBDA50, ASSOCIATION]) - (@view reactionFlux[IKBAA50:IKBDA50, DISSOCIATION]) - (@view reactionFlux[IKBAA50:IKBDA50, TRANSPORTIN]) + (@view reactionFlux[NIKBAA50:NIKBDA50, TRANSPORTOUT]) - (@view reactionFlux[IKBAA50:IKBDA50, DECAY]) - (@view reactionFlux[IKBAA50:IKBDA50, BASALDECAY]);
        # A50:IkB in nucleus
        @. nettFlux[NIKBAA50:NIKBDA50] = (@view reactionFlux[NIKBAA50:NIKBDA50, ASSOCIATION]) - (@view reactionFlux[NIKBAA50:NIKBDA50, DISSOCIATION]) + (@view reactionFlux[IKBAA50:IKBDA50, TRANSPORTIN]) - (@view reactionFlux[NIKBAA50:NIKBDA50, TRANSPORTOUT]) - (@view reactionFlux[NIKBAA50:NIKBDA50, DECAY]) - (@view reactionFlux[NIKBAA50:NIKBDA50, BASALDECAY]);

        #------------------------------------------------
        # MODULE 6
        #--------------------------------------------
        # 5050 free in cytoplasm
        nettFlux[P50P50] += - sum(@view reactionFlux[IKBA5050:IKBD5050, ASSOCIATION]) + sum(@view reactionFlux[IKBA5050:IKBD5050, DISSOCIATION]) + sum(@view reactionFlux[IKBA5050:IKBD5050, DECAY]);
        # 5050 free in nucleus
        nettFlux[NP50P50] += - sum(@view reactionFlux[NIKBA5050:NIKBD5050, ASSOCIATION]) + sum(@view reactionFlux[NIKBA5050:NIKBD5050, DISSOCIATION]) + sum(@view reactionFlux[NIKBA5050:NIKBD5050, DECAY]);
        # 5050:IkB in cytoplasm
        @. nettFlux[IKBA5050:IKBD5050] = (@view reactionFlux[IKBA5050:IKBD5050, ASSOCIATION]) - (@view reactionFlux[IKBA5050:IKBD5050, DISSOCIATION]) - (@view reactionFlux[IKBA5050:IKBD5050, TRANSPORTIN]) + (@view reactionFlux[NIKBA5050:NIKBD5050, TRANSPORTOUT]) - (@view reactionFlux[IKBA5050:IKBD5050, DECAY]) - (@view reactionFlux[IKBA5050:IKBD5050, BASALDECAY]);
        # 5050:IkB in nucleus
        @. nettFlux[NIKBA5050:NIKBD5050] = (@view reactionFlux[NIKBA5050:NIKBD5050, ASSOCIATION]) - (@view reactionFlux[NIKBA5050:NIKBD5050, DISSOCIATION]) + (@view reactionFlux[IKBA5050:IKBD5050, TRANSPORTIN]) - (@view reactionFlux[NIKBA5050:NIKBD5050, TRANSPORTOUT]) - (@view reactionFlux[NIKBA5050:NIKBD5050, DECAY]) - (@view reactionFlux[NIKBA5050:NIKBD5050, BASALDECAY]);

        #------------------------------------------------
        # MODULE 7
        #--------------------------------------------
        # C50 free in cytoplasm
        nettFlux[C50] += - sum(@view reactionFlux[IKBAC50:IKBDC50, ASSOCIATION]) + sum(@view reactionFlux[IKBAC50:IKBDC50, DISSOCIATION]) + sum(@view reactionFlux[IKBAC50:IKBDC50, DECAY]);
        # C50 free in nucleus
        nettFlux[NC50] += - sum(@view reactionFlux[NIKBAC50:NIKBDC50, ASSOCIATION]) + sum(@view reactionFlux[NIKBAC50:NIKBDC50, DISSOCIATION]) + sum(@view reactionFlux[NIKBAC50:NIKBDC50, DECAY]);
        # C50:IkB in cytoplasm
        @. nettFlux[IKBAC50:IKBDC50] = (@view reactionFlux[IKBAC50:IKBDC50, ASSOCIATION]) - (@view reactionFlux[IKBAC50:IKBDC50, DISSOCIATION]) - (@view reactionFlux[IKBAC50:IKBDC50, TRANSPORTIN]) + (@view reactionFlux[NIKBAC50:NIKBDC50, TRANSPORTOUT]) - (@view reactionFlux[IKBAC50:IKBDC50, DECAY]) - (@view reactionFlux[IKBAC50:IKBDC50, BASALDECAY]);
        # C50:IkB in nucleus
        @. nettFlux[NIKBAC50:NIKBDC50] = (@view reactionFlux[NIKBAC50:NIKBDC50, ASSOCIATION]) - (@view reactionFlux[NIKBAC50:NIKBDC50, DISSOCIATION]) + (@view reactionFlux[IKBAC50:IKBDC50, TRANSPORTIN]) - (@view reactionFlux[NIKBAC50:NIKBDC50, TRANSPORTOUT]) - (@view reactionFlux[NIKBAC50:NIKBDC50, DECAY]) - (@view reactionFlux[NIKBAC50:NIKBDC50, BASALDECAY]);
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
    computeNFkBFluxes!(concentration, reactionFlux, Srates, phase; delay, time, birthday, IKKCurve, historicFlux);
    # Compute NFkB net fluxes
    NFkBNettFluxes!(nettFlux, reactionFlux);
    nothing
end
