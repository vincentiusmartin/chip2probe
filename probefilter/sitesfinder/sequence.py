'''
Created on Jul 25, 2019

@author: vincentiusmartin
'''
import collections
import copy

class BindingSite(object):
    '''
    classdocs
    '''


    def __init__(self, site_pos, imads_score, site_start, site_end, core_start, core_end, sequence_in_site, barrier = 1):
        '''
        barrier: #bp from the core that shouldn't be mutated
        '''
        self.site_pos = site_pos
        self.imads_score = imads_score
        self.site_start = site_start
        self.site_end = site_end 
        self.core_start = core_start
        self.core_end = core_end
        self.site_sequence = sequence_in_site
        self.barrier = barrier
        
    def __str__(self):
        return "site_pos: {}, imads_score: {}, site_start: {}, site_end: {}, core_start: {}, core_end: {}, site_sequence {}".format(self.site_pos, self.imads_score, self.site_start, self.site_end, self.core_start, self.core_end, self.site_sequence)

MutatedSequence = collections.namedtuple('MutatedSequence', 'sequence, mutpos, plot_functions')

class Sequence(object):
    '''
    classdocs
    '''

    def __init__(self, escore_preds, imads_preds, escore_cutoff=0.4):
        # we need to check the order between escore and imads
        self.sequence = escore_preds.sequence
        self.bsites = self.get_bsite_escore_imads(escore_preds, imads_preds, escore_cutoff)
        
    def __str__(self):
        return "Sequence object: {}\nsites {}".format(self.sequence,str(self.bsites))

    def mutate_escore_seq_at_pos(self, seq, pos, pbmescore):
        if len(seq) != 8:
            raise Exception("sequence must be at length of 8")
            # TODO: make length can be anything
        nucleotides = "ACGT"
        all_muts = {}
        for n in nucleotides:
            seqlist = list(seq)
            seqlist[pos] = n
            mutseq = "".join(seqlist)
            all_muts[mutseq] = pbmescore.predict_sequence(mutseq).predictions[0]
        minseq = min(all_muts, key=lambda x:all_muts[x]['score'])
        return minseq
    
    def get_non_core_intersecting_maxescore(self, full_seq, site_index, pbmescore):
        """
        we need to include bsite_to_mutate as parameter since we will keep mutating this
        """
        s_start = self.bsites[site_index].site_start
        s_end = self.bsites[site_index].site_end
        site_seq = full_seq[s_start:s_end]
        epreds = pbmescore.predict_sequence(site_seq)
        
        # use large value
        max_escore = {"score":float("inf")}
        while max_escore["score"] == float("inf") and epreds.predictions:
            max_val = max(epreds.predictions, key=lambda x:x['score'])
            max_val_seqpos = self.bsites[site_index].site_start + max_val['position']
            for i in range(len(self.bsites)):
                if i == site_index:
                    continue
                # if the maximum escore is within core of another binding site, then don't mutate
                unmutated_start = self.bsites[i].core_start - self.bsites[i].barrier
                unmutated_end = self.bsites[i].core_end + self.bsites[i].barrier
                if max_val_seqpos >= unmutated_start and max_val_seqpos < unmutated_end:
                    epreds.predictions.remove(max_val)
                else:
                    max_escore = copy.copy(max_val)
        return max_escore
        
    
    def eliminate_site(self, sequence, site_index, pbmescore, escore_threshold = 0.3):
        """
        This assumes that input is a sequence of a site
        """
        maxescore = 1 # just to initialize with a large value
        prevmax = -1 # just to avoid infinite loop
        site_mutpos = []
        full_sequence = str(sequence)
        
        flag = True
        while flag and prevmax != maxescore:    
            # INITIALIZE max escore calculation
            # if we don't need non-core-intersecting: 
            #epreds = pbmescore.predict_sequence(full_sequence)
            #maxepreds = max(epreds.predictions, key=lambda x:x['score'])
            prevmax = float(maxescore)
            maxepreds = self.get_non_core_intersecting_maxescore(full_sequence, site_index, pbmescore)
            maxescore = maxepreds['score'] 
            if maxescore < escore_threshold: # our conidition is met
                flag = False
            else:
                if maxescore == float("inf"): # no e-score that can be chosen
                    print("No e-score site can be mutated for sequence %s" % full_sequence)
                    # since there is no site to be mutated, then just use empty list
                    return full_sequence, []
                seq_tomutate = maxepreds["escore_seq"]
                midpos = len(seq_tomutate) // 2
                mutated_escore_seq = self.mutate_escore_seq_at_pos(seq_tomutate, midpos, pbmescore)
                
                # mutate the sequence
                mut_start = self.bsites[site_index].site_start + maxepreds["start_idx"]
                mut_end = mut_start + len(mutated_escore_seq)
                full_sequence = full_sequence[:mut_start] + mutated_escore_seq + full_sequence[mut_end:]
                site_mutpos.append(mut_start + midpos)
        return full_sequence, site_mutpos
        
    def abolish_sites(self, sites, pbmescore, mode = "to_eliminate"):
        """
        type can either be to_eliminate or to_keep
        """
        sites_to_mutate = []
        if mode == "to_eliminate":
            sites_to_mutate = sites
        elif mode == "to_keep":
            sites_to_mutate = [i for i in range(len(self.bsites)) if i not in sites]
        else:
            raise Exception("type should either be to_eliminate or to_keep")
        
        mutpos = []
        mutated = str(self.sequence)
        for i in sites_to_mutate:
            mutated, site_mutpos = self.eliminate_site(mutated, i, pbmescore)
            mutpos.extend(site_mutpos)
        
        functions = []
        for pos in mutpos:
            functions.append({"func":"axvline","args":[pos],"kwargs":{"color":"purple", "linestyle":"dashed", "linewidth":1}})

        return MutatedSequence(mutated, mutpos, functions)
                
    """  
    prev_idx = idx - 1
    # we don't want to mutate if the position if within the barrier
    mutable_prev = 0 if prev_idx == -1 else bsites[prev_idx].core_end + bsites[prev_idx].barrier
    mutable_cur_start = bsites[idx].site_start
    start_mutseq = max(mutable_prev, mutable_cur_start)
    
    next_idx = idx + 1
    mutable_next = len(seq) if next_idx >= len(bsites) else bsites[next_idx].core_start + bsites[next_idx].barrier
    mutable_cur_end = bsites[idx].site_end
    end_mutseq = min(mutable_next, mutable_cur_end)
    
    to_mutate = seq[start_mutseq:end_mutseq]
    epreds = pbmescore.predict_sequence(to_mutate)
    print(epreds)
    """
    
    def get_bsite_escore_imads(self, escore_preds, imads_preds, escore_cutoff):
        sequence = escore_preds.sequence
        escores = escore_preds.predictions
        signifcount = 0
        startidx = -1
        bindingsites = []
        for i in range(0, len(escores)):
            escoresite = escores[i]
            if escoresite["score"] > escore_cutoff:
                if signifcount == 0:
                    startidx = i
                signifcount += 1
            elif escoresite["score"] < escore_cutoff or i == len(escores)-1:
                if signifcount > 0:
                    # if we have found sufficient e-scores above the cutoff then get the binding sites
                    if signifcount >= 2:
                        # startpos: the start of binding
                        escore_bind = {"startpos":escores[startidx]['position'],  "escorelength":signifcount, 
                                "escore_startidx":escores[startidx]['start_idx']}
                        # we get e-scores now we look for its imads binding site
                        for imads_pred in imads_preds.predictions:
                            escore_start = escore_bind['startpos']
                            escore_end = escore_start + escore_bind["escorelength"]
                            core_start = imads_pred["core_start"]
                            core_end = core_start + imads_pred["core_width"]
                            if escore_start <= core_end and core_start <= escore_end: #overlap
                                # a site is found,
                                # also sequence can be negative or have length above sequence length since we appended
                                # flanking regions, so we need to define site start and site end
                                site_start = max(0,imads_pred["site_start"])
                                site_end = min(imads_pred["site_start"] + imads_pred["site_width"], len(sequence))
                                bsite = BindingSite(site_pos = (core_start + core_end) // 2,
                                                    imads_score = imads_pred["score"],
                                                    site_start = site_start,
                                                    site_end = site_end,
                                                    core_start = imads_pred["core_start"],
                                                    core_end = imads_pred["core_start"] + imads_pred["core_width"],
                                                    sequence_in_site = sequence[site_start:site_end])
                                bindingsites.append(bsite)
                    startidx = -1
                    signifcount = 0        
        return bindingsites
    
    def sites_to_dict(self, bindingsites):
        out_dict = {}
        for i in range(0,len(bindingsites)):
            attrs = [attr for attr in dir(bindingsites[i]) if not callable(getattr(bindingsites[i], attr)) and not attr.startswith("__")]
            for attr in attrs:
                out_dict["%s_%d" % (attr, i+1)] = getattr(bindingsites[i], attr)
        return out_dict
    
    def get_sites_dict(self):
        return self.sites_to_dict(self.bsites)
        
    def site_exist(self):
        return len(self.bsites) != 0
    
    def site_count(self):
        return len(self.bsites)
    