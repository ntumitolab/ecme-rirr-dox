_inhibit_c4(DOX=0μM, KI_DOX_C4=165μM) = hil(KI_DOX_C4, DOX, 3)

function get_c4_eqs(;
    dpsi, cytc_ox, cytc_rd,
    O2=6μM, MT_PROT=1,
    C4_INHIB_DOX=_inhibit_c4(),
    CYANIDE_BLOCK=0
)
    @parameters begin
        rhoC4 = 325μM
        δ₅_C4 = 0.5
        K34_C4 = 2.9445e10Hz / mM^3 ## @pH7
        K43_C4 = 2.9E-6Hz / mM^3
        K35_C4 = 750Hz / mM
        K36_C4 = 4.826e11Hz / mM
        K63_C4 = 4.826Hz / mM
        K37_C4 = 2.92367e6Hz
        K73_C4 = 0.029236Hz ## @pH7
    end

    @variables begin
        vO2(t)
        vHresC4(t)
        C4_Y(t)
        C4_Yr(t)
        C4_YO(t)
        C4_YOH(t)
    end

    C4_CONC = rhoC4 * MT_PROT
    C4_INHIB = C4_INHIB_DOX * (1 - CYANIDE_BLOCK) ## Cyanide is a competitive inhibitor of O2 binding to complex IV
    aδ = exp(-iVT * δ₅_C4 * dpsi)
    a1mδ = exp(iVT * (1 - δ₅_C4) * dpsi)
    fHm = h_m * inv(1E-7Molar)
    fHi = h_i * inv(1E-7Molar)
    f_hm = fHm * aδ
    f_hi = fHi * a1mδ
    f_cr = cytc_rd
    f_co = cytc_ox * a1mδ
    a12 = K34_C4 * f_cr^3 * f_hm^4  # a12 = K34 * exp(-δ₅_C4 * 4 * vfrt) * cytc_rd^3 * h_m^4
    a21 = K43_C4 * f_co^3 * f_hi    # K43 * exp((1 - δ₅_C4) * 4 * vfrt) * cytc_ox^3 * h_i
    a23 = C4_INHIB * K35_C4 * O2
    a32 = 0
    a34 = K36_C4 * f_cr * f_hm^3    # K36 * exp(-δ₅_C4 * 3 * vfrt) * cytc_rd * h_m^3
    a43 = K63_C4 * f_co * f_hi^2    # K63 * exp((1 - δ₅_C4) * 3 * vfrt) * cytc_ox * h_i^2
    a41 = K37_C4 * f_hm  # K37 * exp(-δ₅_C4 * vfrt) * h_m
    a14 = K73_C4 * f_hi  # K73 * exp((1 - δ₅_C4) * vfrt) * h_i

    # Weight of each state (from KA pattern)
    C4_e1 = a41 * a34 * (a21 + a23)
    C4_e2 = a12 * a41 * a34
    C4_e3 = a23 * a12 * (a41 + a43) + a43 * a14 * (a21 + a23)
    C4_e4 = a34 * (a14 * (a21 + a23) + a23 * a12)
    denC4 = C4_e1 + C4_e2 + C4_e3 + C4_e4
    # Reaction rates
    # v34 = C4_CONC * (C4_Y * a12 - C4_Yr * a21)
    v35 = C4_CONC * C4_Yr * a23
    # v36 = C4_CONC * (a34 * C4_YO - a43 * C4_YOH)
    # v37 = C4_CONC * (a41 * C4_YOH - a14 * C4_Y)

    c4eqs = [
        ## C4_CONC ~ cytc_rd + cytc_ox,
        C4_Y ~ C4_e1 / denC4,
        C4_Yr ~ C4_e2 / denC4,
        C4_YO ~ C4_e3 / denC4,
        C4_YOH ~ C4_e4 / denC4,
        vO2 ~ v35,
        vHresC4 ~ 8 * vO2,  ## charge movement
    ]
    return (; c4eqs, vO2, vHresC4)
end
