# 40-state local control CICR (Gauthier 2012), including L-type calcium channel (LCC) and Ryanodine receptor (RyR)
function get_cicr40_sys(ca_i, ca_jsr, ca_o, vm, A_CAP, V_SS_SINGLE=0.203E-6pL; name=:cicrsys)
    @parameters begin
        R_RYR = 19.6 / ms
        R_XFER = 220 / ms
        P_LCC = 9.13E-13 * cm^3 * Hz # Single LCC conductance
        A_LCC = 5
        B_LCC = 7
        γ_LCC = 7.5 / (ms * mM)
        ω_LCC = 0.05 / ms
        f_LCC = 0.85 / ms
        g_LCC = 2 / ms
        KRYR_12 = 5265 / mM^2 / ms
        KRYR_21 = 1500 / ms
        KRYR_23 = 2.358E8 / mM^2 / ms
        KRYR_32 = 9.6 / ms
        KRYR_34 = 1.415E6 / mM^2 / ms
        KRYR_43 = 13.65 / ms
        KRYR_45 = 0.07 / ms
        KRYR_54 = 93.385 / mM^2 / ms
        KRYR_56 = 1.887E7 / mM^2 / ms
        KRYR_65 = 30 / ms
        KRYR_25 = 2.358E6 / mM^2 / ms
        KRYR_52 = 0.001235 / ms
        NCaRU = 50000 # Number of LCC/RyR pairs
    end

    @variables begin
        xca(t)[1:4, 1:10]
        ca_ss(t)[1:4]
        ca_ssavg(t)
        pcicr(t)[1:4]
        α_lcc(t)
        β_lcc(t)
        y_inf(t)
        τ_yca(t)
        ICaL(t)
        Jlcc(t)
        Jryr(t)
        Jxfer(t)
    end

    # Utility functions
    "Subspace Ca macro state"
    _ca(i, j) = ca_ss[(j==3)*2+(i==3)+1]

    "Add LCC transition rates"
    function _lcc!(rs, kf, kb, i, j1, j2, symm=true)
        v = kf * xca[i, j1] - kb * xca[i, j2]
        rs[i, j1] -= v
        rs[i, j2] += v
        if symm
            _lcc!(rs, kf, kb, i, j1 + 5, j2 + 5, false)
        end
        return rs
    end

    "Add RyR transition rates"
    function _ryr!(rs, kf, kb, j, i1, i2)
        v = kf * xca[i1, j] - kb * xca[i2, j]
        rs[i1, j] -= v
        rs[i2, j] += v
    end

    v = vm / mV
    wlcc = P_LCC / V_SS_SINGLE * exprel(-2iVT * vm)
    wryr = R_RYR
    wcai = R_XFER
    α′ = α_lcc * A_LCC
    β′ = β_lcc / B_LCC
    kb = y_inf / τ_yca
    kf = (1 - y_inf) / τ_yca

    rates = fill(Num(0), 4, 10)

    # LCC transition
    for i in 1:4
        _lcc!(rates, α_lcc, β_lcc, i, 1, 2)
        _lcc!(rates, f_LCC, g_LCC, i, 2, 3)
        _lcc!(rates, γ_LCC * _ca(i, 1), ω_LCC, i, 1, 4, false)
        _lcc!(rates, γ_LCC * _ca(i, 6), ω_LCC, i, 6, 9, false)
        _lcc!(rates, α′, β′, i, 4, 5)
        _lcc!(rates, A_LCC * γ_LCC * _ca(i, 2), ω_LCC / B_LCC, i, 2, 5, false)
        _lcc!(rates, A_LCC * γ_LCC * _ca(i, 7), ω_LCC / B_LCC, i, 7, 10, false)
        for j in 1:5
            _lcc!(rates, kf, kb, i, j, j + 5, false)
        end
    end

    # RyR transition
    for j in 1:10
        _ryr!(rates, KRYR_12 * _ca(1, j)^2, KRYR_21, j, 1, 2)
        k32 = KRYR_32 * hil(KRYR_43, KRYR_34 * _ca(3, j)^2)
        _ryr!(rates, KRYR_23 * _ca(2, j)^2, k32, j, 2, 3)
        k42 = KRYR_52 * hil(KRYR_65, KRYR_56 * _ca(4, j)^2)
        _ryr!(rates, KRYR_25 * _ca(2, j)^2, k42, j, 2, 4)
        k34 = KRYR_45 * hil(KRYR_34 * _ca(3, j)^2, KRYR_43)
        k43 = KRYR_54 * _ca(4, j)^2 * hil(KRYR_65, KRYR_56 * _ca(4, j)^2)
        _ryr!(rates, k34, k43, j, 3, 4)
    end

    eqs = [
        ca_ssavg ~ pcicr[1] * ca_ss[1] + pcicr[2] * ca_ss[2] + pcicr[3] * ca_ss[3] + pcicr[4] * ca_ss[4],
        1 ~ pcicr[1] + pcicr[2] + pcicr[3] + pcicr[4],
        pcicr[2] ~ xca[1, 3] + xca[2, 3] + xca[4, 3],
        pcicr[3] ~ xca[3, 1] + xca[3, 2] + xca[3, 4] + xca[3, 5] + xca[3, 6] + xca[3, 7] + xca[3, 8] + xca[3, 9] + xca[3, 10],
        pcicr[4] ~ xca[3, 3],
        ca_ss[1] ~ ca_i,
        ca_ss[2] ~ (wlcc * 0.341ca_o + wcai * ca_i) / (wlcc + wcai),
        ca_ss[3] ~ (wryr * ca_jsr + wcai * ca_i) / (wryr + wcai),
        ca_ss[4] ~ (wlcc * 0.341ca_o + wcai * ca_i + wryr * ca_jsr) / (wlcc + wcai + wryr),
        α_lcc ~ 0.78133 / ms * exp(0.08 * (v - 8)),
        β_lcc ~ 0.34319 / ms * exp(-0.086676 * (v - 8)),
        y_inf ~ 0.3 + 0.7 * expit(-0.2 * (v + 25)),
        τ_yca ~ 20ms + 500ms * expit(-0.1 * (v + 28)) * expit(inv(14) * (v + 42)),
        Jxfer ~ R_XFER * (ca_ssavg - ca_i),
        ICaL ~ NCaRU / A_CAP * P_LCC * (pcicr[2] * ghk(1, vm, ca_ss[2], 0.341 * ca_o, 2) + pcicr[4] * ghk(1, vm, ca_ss[4], 0.341 * ca_o, 2)),
        Jlcc ~ - ICaL * A_CAP / (2 * V_SS_SINGLE * Faraday),
        Jryr ~ R_RYR * (pcicr[3] * (ca_jsr - ca_ss[3]) + pcicr[4] * (ca_jsr - ca_ss[4]))
    ]

    xca1 = Num(1)
    for i in 1:4, j in 1:10
        if i == 1 && j == 1
            continue
        end
        push!(eqs, D(xca[i, j]) ~ rates[i, j])
        xca1 -= xca[i, j]
    end
    push!(eqs, xca[1, 1] ~ xca1)

    return System(eqs, t; name)
end
