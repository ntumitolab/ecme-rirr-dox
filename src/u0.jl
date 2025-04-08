"Initial conditions after 1000 seconds of 1Hz pacing"
function build_u0(sys)
    return [
        sys.sox_i => 0.00011060634476071733,
        sys.h2o2_i => 8.545630293939543e-5,
        sys.gssg_i => 0.23982197972811886,
        sys.Q_n => 1942.499980908948,
        sys.Qdot_n => 84.2161187612923,
        sys.QH2_n => 6.588434982312799,
        sys.QH2_p => 6.550045756223868,
        sys.Qdot_p => 17.607057135837874,
        sys.cytb_1 => 239.11244394703508,
        sys.cytb_2 => 59.08056925495558,
        sys.cytb_3 => 26.111675130317995,
        sys.fes_ox => 284.85387674076117,
        sys.cytc1_ox => 316.5300293836441,
        sys.cytc_ox => 308.61843332069236,
        sys.isoc => 866.3867962737681,
        sys.akg => 0.4875990352014912,
        sys.scoa => 0.5786894917060189,
        sys.suc => 14.928386898447156,
        sys.fum => 14.727965021074025,
        sys.mal => 0.6567134464755138,
        sys.oaa => 5.338029993456324,
        sys.m_na => 0.0014401253920532193,
        sys.h_na => 0.9704371103255879,
        sys.j_na => 0.9787974017848147,
        sys.x_k => 0.002792628904582374,
        sys.ltr_ca => 23.13094176854962,
        sys.htr_ca => 137.62547836417795,
        sys.x_p0 => 0.006610808892067512,
        sys.x_p1 => 0.0058827289555559254,
        sys.x_p2 => 0.01104412905023128,
        sys.x_p3 => 0.009638651292251943,
        sys.x_n1 => 0.018919021226543722,
        sys.adp_ic => 53.54112971739686,
        sys.crp_i => 14700.132866399936,
        sys.crp_ic => 14626.430741004118,
        sys.po1_ryr => 0.00018838992135496773,
        sys.po2_ryr => 2.8561798974549356e-9,
        sys.pc2_ryr => 0.14220903491784556,
        sys.c1_lcc => 1.2401864898506483e-5,
        sys.c2_lcc => 5.7900430310736987e-11,
        sys.c3_lcc => 1.2014179601620536e-16,
        sys.c4_lcc => 9.340060381023616e-23,
        sys.o_lcc => 1.3786274285014654e-23,
        sys.cca0_lcc => 0.0038402708457007267,
        sys.cca1_lcc => 1.912470610924201e-7,
        sys.cca2_lcc => 3.571552246615562e-12,
        sys.cca3_lcc => 2.964394177375789e-17,
        sys.cca4_lcc => 9.226707387348273e-23,
        sys.x_yca => 0.7249578396553411,
        sys.vm => -85.42349534410306,
        sys.na_i => 10190.85030380128,
        sys.k_i => 144723.8904783004,
        sys.ca_i => 0.19040191947361992,
        sys.ca_nsr => 1201.452975673413,
        sys.ca_jsr => 1197.191386641691,
        sys.ca_ss => 0.19329079871841176,
        sys.ca_m => 0.9195634800668644,
        sys.adp_i => 52.70890499537141,
        sys.adp_m => 18.973530551446338,
        sys.nadh_m => 498.4156187334286,
        sys.dpsi => 161.20903046356412,
        sys.sox_m => 0.06725807034930333,
    ]
end
