import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages

class BindingSite:

    def __init__(self,bpos1,bpos2,distance):
        if bpos1 < bpos2:
            self.bpos1 = bpos1
            self.bpos2 = bpos2
        else:
            self.bpos1 = bpos2
            self.bpos2 = bpos1
        self.distance = distance

class SitesFinder:

    def __init__(self, pwmpath, escorepath,
                 pwm_startidx=1,pwm_endidx=-1,pwm_log=True):

        self.pwmdict = self.read_dna_pwm(pwmpath,pwm_startidx,pwm_endidx,pwm_log)
        self.escoredict = self.read_escore(escorepath)

    # ======== PWM RELATED ========

    def read_dna_pwm(self, pwmfile, startidx=1, endidx=-1, log=True):
        start = startidx-1
        with open(pwmfile,'r') as f:
            end = len(f.readline().split(":")[1].split("\t")) if endidx == -1 else endidx
        pwm = dict({})
        with open(pwmfile,'r') as f:
            for line in f:
                base,scores = line.strip().split(":")
                base = base.strip()
                scores = scores.strip().split("\t")[start:end]
                if log:
                    with np.errstate(divide='ignore'): # ignore divide by zero warning
                        pwm[base] = [np.log2(float(score)/0.25) for score in scores]

                else:
                    pwm[base] = [float(score) for score in line.split("\t")[1:]]
        return pwm

    def get_binding_logpwm(self, seq, rectangle=False, rectangle_color='aqua', startrectangleoffset=0, endrectangleoffset=0):
        '''
        rectangle: if true, return a rectangle object for plotting: (x,y),width,height with color as rectangle_color
        '''
        lenpwm = len(self.pwmdict['A'])
        pwmscores = []
        for i in range(0,len(seq)-lenpwm+1):
            score = sum([self.pwmdict[seq[j]][j-i] for j in range(i,i+lenpwm)])
            # below: lenpwm+1 to handle for odd case
            infodict = {"position":i+(lenpwm+1)//2,"seq":seq[i:i+lenpwm],"score":score} #position,sequence,score
            if rectangle:
                # lenpwm-2: because the rectangle width should only cover the binding sites
                infodict["rectangle"]=patches.Rectangle((i+1+startrectangleoffset,0),lenpwm-2-endrectangleoffset,score,
                                  facecolor=rectangle_color,alpha=0.9,edgecolor='black')
            pwmscores.append(infodict)
        return pwmscores

    def lineplot_pwm_escore(self,probes, indexes=[],
                            start_pwm_offset=0, end_pwm_offset=0, bottom_cutoff=0,
                            plot_mut_pos=True,
                            filepath="plot_all.pdf"):
        '''
        offset -> how much differs we want from the beginning and end of the pwm
        matrix
        '''
        if not indexes:
            indexes = probes.table["wt"].index.values

        seqdict =  probes.get_seq("wt",indexes)

        n = 0
        numcol = 4 # adjust the number of columns in the plot
        numrow = 4 # adjust the number of rows in the plot
        with PdfPages(filepath) as pdf:
            for key in seqdict: #seqdict
                if n == 0:
                    fig = plt.figure(figsize=(12,12))
                    fig.subplots_adjust(hspace=0.4,wspace=0.5)
                n+=1

                seq = seqdict[key]
                ax = fig.add_subplot(numcol,numrow,n) # how many plots are in a col,row
                # ================ PWM PLOTTING ================
                pwmsites = self.get_binding_logpwm(seq,True,startrectangleoffset=start_pwm_offset,
                                            endrectangleoffset=end_pwm_offset)
                y_pwm = [x["score"] for x in pwmsites]
                x_pwm = [x["position"] for x in pwmsites]
                ax.plot(x_pwm,y_pwm,linewidth=2.5)
                for site in pwmsites:
                    ax.add_patch(site["rectangle"])

                # ================ ESCORE PLOTTING ================
                escoresites = self.get_binding_escore(seq)
                y_escore = [x["score"] * 10 for x in escoresites]
                x_escore = [x["position"] for x in escoresites]
                ax.plot(x_escore, y_escore, linewidth=2.5, color="darkorange")

                ax.axhline(4,color='red',linestyle='dashed',linewidth=1)

                # ================ MUTATION POSITION PLOTTING ================
                if plot_mut_pos:
                    mutpos = probes.get_mutpos(indexes)
                    for pos in mutpos[key]:
                        ax.axvline(pos,color='purple',linestyle='dashed',linewidth=1)

                # ================ SET PLOTTING LIMIT ================

                ax.axhline(y=0,color='gray')
                ax.set_ylim(bottom=bottom_cutoff,top=max(y_pwm+y_escore)+0.5) # this is based on the pwm
                ax.set_xlim(0.5,)

                # custom the first x axis
                xi = [seq[i] for i in range(0,len(seq))]
                ax.set_xticks(range(1,len(seq)+1))
                ax.set_xticklabels(xi,fontdict={'fontsize':5})
                ax.xaxis.set_label_text('sequence')

                # Make the second axis:
                ax2 = ax.twiny()
                ax2.set_xlim(ax.get_xlim())

                # ================ FINALIZE ================
                ax.yaxis.set_label_text('PWM log score')
                ax.tick_params(direction="in")
                ax2.tick_params(direction="in")

                ax.set_title("Sequence %d"%key,pad=20)
                if n == numcol*numrow:
                    pdf.savefig(fig)
                    plt.close()
                    n = 0
            pdf.savefig(fig)
            plt.close()

    # ======== ESCORE RELATED ========

    def read_escore(self,escorefile):
        df = pd.read_csv(escorefile,sep="\t")
        escoredict = {}
        for index, row in df.iterrows():
            escoredict[row["8-mer"]] = row["E-score"]
            escoredict[row["8-mer.1"]] = row["E-score"] # for the reverse complement
        return escoredict

    def get_binding_escore(self,seq):
        seq_escore = []
        kmer = 8
        for i in range(0,len(seq)-kmer+1):
            score = self.escoredict[seq[i:i+kmer]]
            seq_escore.append({"position":i+(kmer+1)//2,"seq":seq[i:i+kmer],"score":score,"listindex":i})
        return seq_escore

    # ======== FILTER RELATED ========
    def binding_sites(self,pwmsites,escoresites,seqid=-1):
        E_CUTOFF = 0.4
        PWMRANGE = 3
        # ^ hard coded just because our pwmseqlen is 4 while the motif is len 6

        signifcount = 0
        startidx = -1
        bindingsites = []
        for i in range (0,len(escoresites)):
            escoresite = escoresites[i]
            pwmsite = pwmsites[i]
            #pwmsite: {'position': 24, 'seq': 'CCACATGA', 'score': 2.937143374527686}
            if escoresite["score"] >= E_CUTOFF:
                if signifcount == 0:
                    startidx = i
                signifcount += 1
            if escoresite["score"] < 0.4 or i == len(escoresites)-1:
                if signifcount > 0:
                    if signifcount >= 2:
                        startpwm = max(startidx-PWMRANGE,0)
                        endpwm = i+PWMRANGE
                        maxpwm = max(pwmsites[startpwm:endpwm],key=lambda x:x['score'])

                        # startpos: the start of binding
                        bind = {"startpos":escoresites[startidx]['position'],  "escorelength":signifcount, "escorelistidx":escoresites[startidx]['listindex'],
                        "position":maxpwm["position"],"pwmscore":maxpwm["score"]}
                        bindingsites.append(bind)
                    startidx = -1
                    signifcount = 0
        return bindingsites

    def validate_bindsites(self, bindsites, escoresites, nonspecific_ecutoff=0.35):
        if len(bindsites) != 2:
            return False,"Error: Found %d binding sites\n" % len(bindsites),0

        # since we only have sites with size 2 and we know the distance between sites
        # should be between 10-20 bp
        distance = abs(bindsites[1]["position"]-bindsites[0]["position"])
        if distance < 10 or distance > 20:
            return False,"Error: Incorrect distance between sites: %d\n" % distance,distance

        flag = False
        for i in range(bindsites[0]["escorelistidx"]+bindsites[0]["escorelength"],bindsites[1]["escorelistidx"]):
            if escoresites[i]["score"] < nonspecific_ecutoff:
                flag = True
                break

        if not flag:
            return False,"Error: Position between binding sites don't have site with escore < 0.35: %d\n" % escoresites[i]["score"],distance

        bsite = BindingSite(bindsites[0]["position"],bindsites[1]["position"],distance)
        return True,str(bindsites)+"\n",bsite

    def filter_sequences(self,probes,indexes=[],tofile=False,label="sequence"):
        """
        returntype: linear or table
        """
        table = []
        log = ""

        # when getting binding site, always use "wt"
        seqdict =  probes.get_seq("wt",indexes)

        for key in seqdict:
            if indexes and key not in indexes:
                continue
            seq = seqdict[key]
            log += "--Sequence %d\n"%key
            pwmsites = self.get_binding_logpwm(seq)
            escoresites = self.get_binding_escore(seq)
            bsites = self.binding_sites(pwmsites,escoresites)
            validated,msg,bsite = self.validate_bindsites(bsites,escoresites)
            log += msg
            if validated:
                seq = probes.get_seq("wt",[key])[key]
                table.append({"index":key,
                              "sequence":seq,
                              "bpos1":bsite.bpos1,
                              "bpos2":bsite.bpos2,
                              "distance":bsite.distance,
                              "label":label})

        if tofile:
            '''
            with open("%s-idxs.txt"%label,'w') as f:
                f.write(">%s_filtered\n"%label)
                f.write(",".join(str(idx) for idx in filtered_idxs))
                f.write("\n")
                f.write(">%s_distance\n"%label)
                f.write(",".join(str(idx) for idx in distances))
            '''
            with open("%s-log.txt"%label,'w') as f:
                f.write(log)

        return pd.DataFrame(table)

    # TODO: make this more general
    def filtertrain_to_csv(self,probes,categories,filename="training.csv"):
        coop_tbl = self.filter_sequences(probes,indexes=categories["coop_overlap"],label="cooperative")
        additive_tbl = self.filter_sequences(probes,indexes=categories["additive_overlap"],label="additive")

        print("Number filtered cooperative %d"%len(coop_tbl))
        print("Number filtered additive %d"%len(additive_tbl))

        combined = pd.concat([coop_tbl,additive_tbl])
        combined[["sequence","bpos1","bpos2","distance","label"]].to_csv(filename,index=False)
