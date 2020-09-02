
if __name__ == "__main__":
    # ------- labeling sequences, for now only take wt -------
    keyword = "all_clean_seqs"
    df_genomics = []
    for df in dflist:
        df_gen = df[df["Name"].str.contains(keyword,na=False)].sort_values(by=["Name"])
        df_gen[["Name","type","rep","ori"]] = df_gen["Name"].str.rsplit("_", n = 3, expand = True)
        df_gen = df_gen[df_gen["type"] == "wt"] #only take wt
        df_genomics.append(df_gen)

    # we filter each group if median intensity is lower than cutoff
    oris_df = []

    for ori in ["o1","o2"]:
        df1 = df_genomics[0][df_genomics[0]["ori"] == ori][["Name","rep","Alexa488Adjusted"]]
        df2 = df_genomics[1][df_genomics[1]["ori"] == ori][["Name","rep","Alexa488Adjusted"]]
        df_comb = df1.merge(df2, on=["Name","rep"])
        print("Number of distinct names before filtering %d" % df_comb["Name"].nunique())
        median_df = df_comb.groupby("Name").median().reset_index()
        above_cut = median_df[(median_df["Alexa488Adjusted_x"] > cutoff) & (median_df["Alexa488Adjusted_y"] > cutoff)]
        print("Number of distinct names after filtering %d" % above_cut["Name"].nunique())

        df_comb = df_comb.merge(above_cut[["Name"]], on=["Name"])
        lbled = df_comb.groupby("Name").apply(lambda x :
            pd.Series({
                "p_coop": st.wilcox(
                        x["Alexa488Adjusted_y"].tolist(),
                        x["Alexa488Adjusted_x"].tolist(),
                        "greater"
                    ),
                "p_anti": st.wilcox(
                        x["Alexa488Adjusted_y"].tolist(),
                        x["Alexa488Adjusted_x"].tolist(),
                        "less"
                    ),
            })).reset_index()
        # too low statistical power, can't do fdr correction
        # lbled["p_coop"] = sm.fdrcorrection(lbled["p_coop"])[1]
        # lbled["p_anti"] = sm.fdrcorrection(lbled["p_anti"])[1]
        lbled["label"] = lbled.apply(
            lambda x:
                "cooperative" if x["p_coop"] < 0.06 else
                "anticoop" if x["p_anti"] < 0.06 else
                "additive",
                axis=1
            )

        median_df.merge(lbled, on="Name").to_csv("test_%s.csv"%ori)
        df_lbled = median_df.merge(lbled, on="Name")[["Name","Alexa488Adjusted_x","Alexa488Adjusted_y","label"]]
        arr.plot_classified_labels(df_lbled, filepath="%s.png"%ori)

    o1df = pd.read_csv("test_o1.csv")
    o2df = pd.read_csv("test_o2.csv")
    olist = [o1df,o2df]
    both_ori = olist[0].merge(olist[1], on=["Name"])
    both_ori["Alexa488Adjusted_x"] = (both_ori["Alexa488Adjusted_x_x"] + both_ori["Alexa488Adjusted_x_y"]) / 2
    both_ori["Alexa488Adjusted_y"] = (both_ori["Alexa488Adjusted_y_x"] + both_ori["Alexa488Adjusted_y_y"]) / 2
    both_ori["label"] = both_ori.apply(lambda x:
        "cooperative" if x["label_x"] == "cooperative" and x["label_y"] == "cooperative" else
        "anticoop" if x["label_x"] == "anticoop" and x["label_y"] == "anticoop" else
        "additive",
        axis=1
    )
    both_ori = both_ori[["Name","Alexa488Adjusted_x","Alexa488Adjusted_y","label"]]
    arr.plot_classified_labels(both_ori, filepath="both.png")

    print(o1df[["label","Name"]].groupby("label").count())
