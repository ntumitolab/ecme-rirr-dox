# mtDNA model by de Oliveira et al.
@with_kw struct MTDNAParams
    MT_DNA = 0.75  # Default Mitocondrial DNA content
    α = 5.625e-11  # Rate constant of mtDNA reapir
    β = 6.940e-11  # reaction rate constant of hydroxyl radical and mtDNA
    γ = 8.100e-9  # reaction rate constant of DOX and mtDNA
    κ = 0.045  # half saturation concentration of mtDNA repair
    K_B = 1.42   # rate constant of mt protein synthesis
    K_PROT = 0.32  # half saturation concentration of mt protein synthesis
    MT_PROT = K_B * hil(MT_DNA, K_PROT)
end
