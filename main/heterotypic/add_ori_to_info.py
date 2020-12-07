
import pandas as pd

if __name__ == "__main__":
    basepath = "/Users/vincentiusmartin/Research/chip2gcPBM/chip2probe/output/heterotypic/EtsRunx_v1/customseq"
    info = pd.read_csv("%s/custom_info_fz.csv" % basepath)
    origori = info[info["sequence"] == info["ori_seq"]][["ori_seq","ets_pos","runx_pos"]]
    origori["orig_ori"] = info.apply(lambda x: "er" if x["ets_pos"] < x["runx_pos"] else "re",axis=1)
    origori = origori[["ori_seq","orig_ori"]]
    info.merge(origori, on="ori_seq").to_csv("custom_info_fz_worigori.csv", index=False)
