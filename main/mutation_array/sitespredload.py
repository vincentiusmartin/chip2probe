from chip2probe.sitespredict.kompas import Kompas
from chip2probe.sitespredict.pwm import PWM
from chip2probe.sitespredict.imads import iMADS
from chip2probe.sitespredict.imadsmodel import iMADSModel
from chip2probe.sitespredict.kompaspwm import KompasPWM

# ========= iMADS =========
basesp = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe"
imads_ets_paths = ["%s/input/sitemodels/imads_models/Ets1_w12_GGAA.model" % basesp, "%s/input/sitemodels/imads_models/Ets1_w12_GGAT.model" % basesp]
imads_ets_cores = ["GGAA", "GGAT"]
imads_ets_models = [iMADSModel(path, core, 12, [1, 2, 3]) for path, core in zip(imads_ets_paths, imads_ets_cores)]
imads_ets = iMADS(imads_ets_models, 0.19) # 0.2128

imads_runx_paths = ["%s/input/sitemodels/imads_models/Runx1_w20_GAGGT.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GCGGC.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GCGGG.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GCGGT.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GTGGC.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GTGGG.model" % basesp,
                    "%s/input/sitemodels/imads_models/Runx1_w20_GTGGT.model" % basesp]
imads_runx_cores = ["GAGGT", "GCGGC", "GCGGG", "GCGGT", "GTGGC", "GTGGG", "GTGGT"]
imads_runx_models = [iMADSModel(path, core, 20, [1, 2, 3]) for path, core in zip(imads_runx_paths, imads_runx_cores)]
imads_runx = iMADS(imads_runx_models, 0.25) # 0.2128

# ========= PWM =========
pwm_runx = PWM("%s/input/sitemodels/pwm/runx1.txt"%basesp, 8, 17, log=True, reverse=True)
pwm_ets = PWM("%s/input/sitemodels/pwm/ets1.txt"% basesp, log=True, reverse=False)

# ========= KOMPAS =========
kompas_ets = Kompas("%s/input/sitemodels/kompas/Ets1_kmer_alignment.txt"% basesp, core_start = 11, core_end = 15, core_center = 12)
kompas_runx = Kompas("%s/input/sitemodels/kompas/Runx1_kmer_alignment.txt"% basesp, core_start = 12, core_end = 17, core_center = 14)

# ========= KompasPWM =========
kp_ets = KompasPWM(kompas_ets, pwm_ets)
kp_runx = KompasPWM(kompas_runx, pwm_runx)
