

def pwm_score(dfpos, pwm, startcol, corelen, flanklen, seqcol="Sequence"):
    """
    e.g. pwm_score(dft, pwm_ets, "ets_start", 4, 3)
    """
    res = []
    ori = []
    core = []
    for idx, row in dfpos.iterrows():
        seq = row[seqcol][row[startcol]-flanklen:row[startcol]+corelen+flanklen]
        core.append(row[seqcol][row[startcol]:row[startcol]+corelen])
        if len(seq) < (corelen + 2*flanklen): # error
            res.append(-999)
            ori.append(-999)
        else:
            pred = pwm.predict_sequence(seq,zero_thres=False)[0]
            res.append(pred["score"])
            cur_ori = 0 if pred["orientation"] == -1 else 1
            ori.append(cur_ori)
    return res, ori, core
