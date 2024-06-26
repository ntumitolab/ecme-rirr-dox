#=
ROS scavenging system
=#

"""
Rate of superoxide dismutase. Based on (McAdam, 1977)
"""
function _vsod(sox, h2o2, K1, K3, K5, KI_H2O2, E0)
	f_h2o2 = K3 * (1 + h2o2 / KI_H2O2)
	f_sox  = K1 * sox
	denom = K5 * (2 * K1 + f_h2o2) + f_h2o2 * f_sox
	return 2 * E0 * K5 * f_sox * (K1 + f_h2o2) / denom
end

"ROS diffusion and detox system"
function get_ros_sys(dpsi, sox_m, nadph_i, V_MITO_V_MYO=0.615; name=:rossys)
	@parameters begin
		# SOD params
		K1_SOD = 1200/mM/ms  	# 2nd order rate constant of SOD
		K3_SOD = 24/mM/ms  		# 2nd order rate constant of SOD
		K5_SOD = 2.4E-4/ms  	# 1st order rate constant of SOD
		KI_H2O2_SOD = 0.5mM  	# Inhibition constant of H2O2
		# ET_SOD_I = 1.43μM  	# Cytosolic SOD concentration (Zhou, 2009)
		ET_SOD_I = 3μM			# Cytosolic SOD concentration
		ET_SOD_M = 0.3μM  		# Mitochondrial SOD concentration
		# glutathione peroxidase (GPX)
		𝚽1_GPX = 5E-3mM*ms  # Rate constant of GPX
		𝚽2_GPX = 0.75mM*ms  # Rate constant of GPX
		ET_GPX = 10μM  # GPX concentration ()
		# thioredoxin peroxidase (TPX)
		𝚽1_TPX = 3.83mM*ms  # Rate constant of TPX
		𝚽2_TPX = 1.85mM*ms  # Rate constant of TPX
		ET_TPX = 3μM  # TPX concentration
		# GR (glutathion reductace) parameters
		K1_GR = 5Hz  # Catalytic constant of GR
		ET_GR = 10μM  # Mitochondrial concentration of GR
		KM_NADPH_GR = 15μM  # Michaelis constant for NADPH of GR
		KM_GSSG_GR = 60μM  # Michaelis constant for oxidized GSH of GR
		ΣGSH_i = 1mM  # Glutathione pool (mM)
		# Thioredoxin reductase (TR) parameters
		K1_TR = 22.75Hz  # Catalytic constant
		ET_TR = 3.5E-4mM  # Enzyme concentration
		KM_NADPH_TR = 0.065mM  # Michaelis constant for NADPH
		KM_TRXSS_TR = 0.035mM  # Michaelis constant for oxidized thioredoxin of TR
		TRX_T_TR = 0.025mM  # Thioredoxin pool(mM)
		# Catalase
		K1_CAT = 17.0/(mM*ms)  # Rate constant of CAT
		ET_CAT = 0.01mM  # Total pool of CAT
		FR_CAT = 0.05/mM  # H2O2 inhibition factor of CAT
		# IMAC (Inner mitochondrial anion channel) from Cortassa et al. (2004)
		A_IMAC = 1E-3  # Basal IMAC conductance factor
		B_IMAC = 1E4   # Activation IMAC conductance factor by cytoplasmic superoxide
		KCC_SOX_IMAC = 10μM  # Michaelis constant for cytoplasmic superoxide of IMAC
		GL_IMAC = 3.5E-8mM/ms/mV  		# Leak conductance of IMAC (Zhou, 2009)
		G_MAX_IMAC = 3.9085E-6mM/ms/mV  # Maximal conductance of IMAC (Zhou, 2009)
		κ_IMAC = 0.07/mV  # Steepness factor
		DPSI_OFFSET_IMAC = 4mV  # Potential at half saturation
		J_IMAC = 0.1  # Fraction of IMAC conductance
	end

	@variables begin
		sox_i(t) = 8.57e-9mM
		h2o2_i(t) = 2.062e-9mM
		# h2o2_m(t)
		gsh_i(t) # Conserved
		# gsh_m(t)
		gssg_i(t) = 5e-6mM
		# gssg_m(t)
		# nadph_m(t)
		vSOD_i(t)
		vGPX_i(t)
		vGR_i(t)
		vCAT(t)   # Catalase flux
		vTrROS(t) # SOX flux via IMAC
		vIMAC(t)  # IMAC ion flux
		ΔVROS(t)  # Reversal potential of ROS
	end

	fv = GL_IMAC + G_MAX_IMAC * expit(κ_IMAC * (dpsi - DPSI_OFFSET_IMAC))
	gimac = (A_IMAC + B_IMAC * hil(sox_i, KCC_SOX_IMAC)) * fv
	vsod_i = _vsod(sox_i, h2o2_i, K1_SOD, K3_SOD, K5_SOD, KI_H2O2_SOD, ET_SOD_I)

	eqs = [
		ΔVROS ~ nernst(sox_i, sox_m, -1),
		vTrROS ~ J_IMAC * gimac * (dpsi + ΔVROS),
		vIMAC ~ gimac * dpsi,
		vGR_i ~ ET_GR * K1_GR * hil(nadph_i, KM_NADPH_GR) * hil(gssg_i, KM_GSSG_GR),
		vGPX_i ~ ET_GPX * h2o2_i * gsh_i / (𝚽1_GPX * gsh_i + 𝚽2_GPX * h2o2_i),
		vCAT ~ 2 * K1_CAT * ET_CAT * h2o2_i * exp(-FR_CAT * h2o2_i),
		vSOD_i ~ vsod_i,
		ΣGSH_i ~ gsh_i + 2 * gssg_i,
		D(sox_i) ~ V_MITO_V_MYO * vTrROS - vSOD_i,
		D(h2o2_i) ~ 0.5 * vSOD_i - vGPX_i - vCAT,
		D(gssg_i) ~ -0.5 * (vGR_i - vGPX_i)
	]
	return ODESystem(eqs, t; name)
end
