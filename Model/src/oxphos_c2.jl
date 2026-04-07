_inhibit_c2(DOX=0μM, KI_DOX_C2=2000μM) = hil(KI_DOX_C2, DOX, 3)

# Complex II (SDH)
# Reversible rapid equlibrium random Bi-Bi enzyme catalytic mechanism
function get_c2_eqs(;
    vSDH, Q_n, QH2_n,
    succinate, fumarate, oaa,
    C2_INHIB_DOX=_inhibit_c2(),     ## Doxorubicin concentration
    MT_PROT=1,                      ## OXPHOS protein content scale factor
)

    @parameters begin
        ## Reaction rate constant of SDH (complex II)
        VF_C2 = 250mM / minute
        ## Inhibition constant for OAA
        KI_OAA_C2 = 150μM
        KM_SUC_C2 = 30μM
        KM_Q_C2 = 0.3μM
        KM_FUM_C2 = 25μM
        KM_QH2_C2 = 1.5μM
        ## midpoint potential of fumarate -> succinate
        Em_FUM_SUC = 40mV
        ## midpoint potential of Q -> QH2
        Em_Q_QH2 = 100mV
        ## equlibrium constant of SDH
        KEQ_C2 = exp(2iVT * (Em_Q_QH2 - Em_FUM_SUC))
        ## Reverse rate following Haldane relationship
        VR_C2 = VF_C2 * KM_FUM_C2 * KM_QH2_C2 / (KEQ_C2 * KM_SUC_C2 * KM_Q_C2)
    end

    C2_INHIB = hil(KI_OAA_C2, oaa) * C2_INHIB_DOX * MT_PROT
    A = succinate / KM_SUC_C2
    B = Q_n / KM_Q_C2
    P = fumarate / KM_FUM_C2
    Q = QH2_n / KM_QH2_C2
    eqs_c2 = [vSDH ~ C2_INHIB * (VF_C2 * A * B - VR_C2 * P * Q) / (1 + A + B + A * B + P + Q + P * Q)]
    return eqs_c2
end
