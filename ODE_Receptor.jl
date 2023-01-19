# include("HelperFunctions.jl");
# Compute reaction fluxes for CD40 module
#--------------------------------------------
function computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, time)
    @inbounds begin
        # MODULE 1 : CD40 ligand-receptor interaction
        # 1 : CD40L basal decay
        reactionFlux[CD40L, BASALDECAY] = Srates[CD40L, BASALDECAY] * concentration[CD40L];
        # 2 : CD40R basal synthesis
        reactionFlux[CD40R, BASALSYNTHESIS] = Srates[CD40R, BASALSYNTHESIS];
        # 3 : CD40R basal decay
        reactionFlux[CD40R, BASALDECAY] = Srates[CD40R, BASALDECAY] * concentration[CD40R];
        # 4 : CD40L-R association
        reactionFlux[CD40LR, ASSOCIATION] = Srates[CD40LR, ASSOCIATION] * concentration[CD40L] * concentration[CD40R];
        # 5 : CD40L-R dissociation
        reactionFlux[CD40LR, DISSOCIATION] = Srates[CD40LR, DISSOCIATION] * concentration[CD40LR];
        # 6 : CD40L-R basal decay
        reactionFlux[CD40LR, BASALDECAY] = Srates[CD40LR, BASALDECAY] * concentration[CD40LR];

        # MODULE 2 : Canonical pathway activation
        # 7 : TRAF6 activation (induced by CD40L-R complex)
        reactionFlux[TRAF6, ACTIVATION] = Srates[TRAF6, ACTIVATION] * concentration[CD40LR] * concentration[TRAF6_OFF];
        # 8 : TRAF6 deactivation
        reactionFlux[TRAF6, DEACTIVATION] = Srates[TRAF6, DEACTIVATION] * concentration[TRAF6];
        # 9 : TAK1 activation (induced by TRAF6)
        reactionFlux[TAK1, ACTIVATION] = Srates[TAK1, ACTIVATION] * concentration[TRAF6] * concentration[TAK1_OFF];
        # 10 : TAK1 deactivation
        reactionFlux[TAK1, DEACTIVATION] = Srates[TAK1, DEACTIVATION] * concentration[TAK1];
        # 11 : IKK activation
        reactionFlux[IKK, ACTIVATION] = Srates[IKK, ACTIVATION] * concentration[IKK_OFF] * hillActivation(Ka = 3.0, na = 2.0, X = concentration[TAK1]);
        # 12 : IKK basal activation
        reactionFlux[IKK, BASALSYNTHESIS] = Srates[IKK, BASALSYNTHESIS] * concentration[IKK_OFF];
        # 13 : IKK deactivation (IKK inhibition)
        reactionFlux[IKK, DEACTIVATION] = Srates[IKK, DEACTIVATION] * concentration[IKK];
        # 14 : IKK renewal
        reactionFlux[IKK_OFF, SYNTHESIS] = Srates[IKK_OFF, SYNTHESIS] * concentration[IKK_I];

        # MODULE 3 : Noncanonical pathway activation
        # 15 : TRAF3 complex basal synthesis
        reactionFlux[TRAF3, BASALSYNTHESIS] = Srates[TRAF3, BASALSYNTHESIS];
        # 16 : TRAF3 complex basal decay
        reactionFlux[TRAF3, BASALDECAY] = Srates[TRAF3, BASALDECAY] * concentration[TRAF3];
        # 17 : TRAF3 complex decay (induced by CD40L-R complex)
        reactionFlux[TRAF3, DECAY] = Srates[TRAF3, DECAY] * concentration[CD40LR] * concentration[TRAF3];
        # 18 : NIK basal synthesis
        reactionFlux[NIK, BASALSYNTHESIS] = Srates[NIK, BASALSYNTHESIS];
        # 19 : NIK basal decay
        reactionFlux[NIK, BASALDECAY] = Srates[NIK, BASALDECAY] * concentration[NIK];
        # 20 : NIK decay (induced by TRAF3 complex)
        reactionFlux[NIK, DECAY] = Srates[NIK, DECAY] * concentration[NIK] * hillActivation(Ka = 10.0, na = 2.0, X = concentration[TRAF3]);

    end
    nothing
end

function CD40NettFluxes!(nettFlux, reactionFlux)
    @inbounds begin
        # 1 : CD40 ligand
        nettFlux[CD40L] = - reactionFlux[CD40L, BASALDECAY] - reactionFlux[CD40LR, ASSOCIATION] * SCALE_CELLULAR2MEDIA + reactionFlux[CD40LR, DISSOCIATION] * SCALE_CELLULAR2MEDIA;
        # 2 : CD40 receptor
        nettFlux[CD40R] = reactionFlux[CD40R, BASALSYNTHESIS] - reactionFlux[CD40R, BASALDECAY] - reactionFlux[CD40LR, ASSOCIATION] + reactionFlux[CD40LR, DISSOCIATION];
        #------------------------------------------------
        # 3 : CD40L-R complex
        nettFlux[CD40LR] = reactionFlux[CD40LR, ASSOCIATION] - reactionFlux[CD40LR, DISSOCIATION] - reactionFlux[CD40LR, BASALDECAY];

        # 4 : inactivated TRAF6
        nettFlux[TRAF6_OFF] = - reactionFlux[TRAF6, ACTIVATION] + reactionFlux[TRAF6, DEACTIVATION];
        # 5 : activated TRAF6
        nettFlux[TRAF6] = reactionFlux[TRAF6, ACTIVATION] - reactionFlux[TRAF6, DEACTIVATION];
        # 6 : inactivated TAK1
        nettFlux[TAK1_OFF] = - reactionFlux[TAK1, ACTIVATION] + reactionFlux[TAK1, DEACTIVATION];
        # 7 : activated TAK1
        nettFlux[TAK1] = reactionFlux[TAK1, ACTIVATION] - reactionFlux[TAK1, DEACTIVATION];
        # 8 : inactivated IKK
        nettFlux[IKK_OFF] = - reactionFlux[IKK, ACTIVATION] - reactionFlux[IKK, BASALSYNTHESIS] + reactionFlux[IKK_OFF, SYNTHESIS];
        # 9 : inhibited IKK
        nettFlux[IKK_I] = reactionFlux[IKK, DEACTIVATION] - reactionFlux[IKK_OFF, SYNTHESIS];
        # 10 : IKK (activated IKK)
        nettFlux[IKK] = reactionFlux[IKK, ACTIVATION] + reactionFlux[IKK, BASALSYNTHESIS] - reactionFlux[IKK, DEACTIVATION];

        # 11 : TRAF3-cIAP1/2 complex
        nettFlux[TRAF3] = reactionFlux[TRAF3, BASALSYNTHESIS] - reactionFlux[TRAF3, BASALDECAY] - reactionFlux[TRAF3, DECAY];
        # 12 : NIK
        nettFlux[NIK] = reactionFlux[NIK, BASALSYNTHESIS] - reactionFlux[NIK, BASALDECAY] - reactionFlux[NIK, DECAY];
    end
    nothing
end

# time-independent ODE (pre-simulation: phase = 1)
function computeCD40NettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, nothing);
    # Compute NFkB net fluxes
    CD40NettFluxes!(nettFlux, reactionFlux);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeCD40NettFluxes!(nettFlux, concentration, reactionFlux, Srates, phase, time)
    #------------------------------------------------
    # Compute NFkB reaction fluxes
    computeCD40Fluxes!(concentration, reactionFlux, Srates, phase, time);
    # Compute NFkB net fluxes
    CD40NettFluxes!(nettFlux, reactionFlux);
    nothing
end
