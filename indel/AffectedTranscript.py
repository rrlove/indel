import allel
import itertools

class AffectedTranscript():
    
    def __init__(self,chrom,name):
        
        self.chrom = chrom
        self.name = name
        self.intersecting_variants = []
        self.vtbl_indices = []
        self.variant_length_changes = []
        self.n_variants = 0
    
    '''a collection of variant records, the index to access the associated 
    genotypes, and the length change(s) of each'''
    
    def add_variant(self,variant):
        
        self.intersecting_variants.append(variant)
        self.vtbl_indices.append(variant.vt_index)
        self.variant_length_changes.append(set(variant.length_change))
        self.n_variants += 1
    
    ##run this method after you have added all variants
    def extract_vtbl(self,vtbl):

        assert len(self.vtbl_indices) > 0,\
        "This transcript has no associated variants"
        self.vtbl_indices.sort()
        self.vtbl = vtbl[vtbl["CHROM"] ==\
                         self.chrom.encode()][[self.vtbl_indices]]

    ##for each variant, keep track of all possible length changes
    def calculate_all_possible_length_changes(self):
                        
        ##calculate all possible "paths" 
        ##through the collection of possible length changes
        ##first, assemble that collection of possible length changes
        all_possible_changes =\
        list(itertools.product(*self.variant_length_changes))
        
        ##now, sum all of these and keep the unique set
        self.possible_length_changes =\
        set([sum(i) for i in all_possible_changes])
        
    def extract_haplotypes(self,phased_genotype_array_by_chrom):
        
        self.genos_phased =\
        phased_genotype_array_by_chrom[self.chrom][[self.vtbl_indices]]
        
        haplotypes_mother = self.genos_phased[:,0].to_haplotypes()
        haplotypes_father = self.genos_phased[:,1].to_haplotypes()
        
        self.haplotypes =\
        haplotypes_mother.concatenate(haplotypes_father, axis=1)
             
    def calculate_haplotype_based_length_changes(self):
        
        self.haplotype_length_changes = {}

        ##for all the different haplotypes for this transcript, 
        ##calculate net change per haplotype
        ##I DON'T want to remove variants that can't be phased, I think; 
        ##right now, I only want
        ##to narrow down the possible combinations of length changes 
        ##to make them computationally tactable
        
        for haplotype_index in range(self.haplotypes.shape[1]):
            changes = []
            haplotype = self.haplotypes[:, haplotype_index]
            
            ##get the index of each variant in the haplotype, 
            ##and the variant itself
            
            for variant_index,genotype_index in enumerate(haplotype):
                
                if genotype_index == 0:
                    ##keep track of the cumulative length change of the 
                    ##haplotype, which is 0 if the variant is REF
                    changes.append(0)
                
                else:
                    ##if the variant is not REF, calculate the length change
                    change =\
                    len(self.intersecting_variants[variant_index]\
                        .record["ALT"][0][genotype_index-1]) -\
                    len(self.intersecting_variants[variant_index]\
                        .record["REF"][0])
                    changes.append(change)
    
            self.haplotype_length_changes[tuple(haplotype)] = sum(changes)
