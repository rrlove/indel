import allel

class AffectedTranscript():
    
    def __init__(self,chrom,name,feature_table_entry):
        self.chrom = chrom
        self.name = name
        self.feature = feature_table_entry
        self.variants = []
        self.length_changes = []
        self.vtbl_indices = []
        
    def add_variant(self,variant):
        self.variants.append(variant)
        self.vtbl_indices.append(variant.vtbl_index)
        
    def count_n_variants(self):
        self.n_variants = len(self.variants)
        
    def count_length_change(self):
        for variant in self.variants:
            self.length_changes.append(set(variant.length_change))
        
        self.length_change_combos = set([sum(list(itertools.product(*self.length_changes))[i]) for i in range(len(list(itertools.product(*self.length_changes))))])
        
    def phase(self,phased_genotype_array_by_chrom):
        self.genos_phased = genosDict[chrom][testDict[name].vtbl_indices]
        ##self.genos_phased = allel.phase_by_transmission(self.genos_chunk, window_size=len(self.genos_chunk))
        
    def extract_haplotypes(self):
        self.genotypes_mother = self.genos_phased[:,0]
        self.haplotypes_mother = self.genotypes_mother.to_haplotypes()
        self.haplotypes_progeny_maternal = allel.HaplotypeArray(self.genos_phased[:, 2:, 0])
        
        self.genotypes_father = self.genos_phased[:,1]
        self.haplotypes_father = self.genotypes_father.to_haplotypes()
        self.haplotypes_progeny_paternal = allel.HaplotypeArray(self.genos_phased[:, 2:, 1])
        
        self.haplotypes = self.haplotypes_mother.concatenate(self.haplotypes_father, axis=1).concatenate(self.haplotypes_progeny_maternal, axis=1).concatenate(self.haplotypes_progeny_paternal, axis=1)
        
        self.length_change_list = []
        self.unique_haplotypes = {}

        for c in range(self.haplotypes.shape[1]):
            changes = []
            haplotype = []
            for a,b in enumerate(self.haplotypes[:,c]):
                haplotype.append(b)
                if b == 0:
                    changes.append(0)
                else:
                    change = len(self.variants[a].record["ALT"][0][b-1]) - len(self.variants[a].record["REF"][0])
                    changes.append(change)
    
            self.length_change_list.append(sum(changes))
            self.unique_haplotypes[tuple(haplotype)] = sum(changes)
    
        self.length_changes = set(self.length_change_list)
