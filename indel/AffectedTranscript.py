import numpy as np

class AffectedTranscript():
    
    def __init__(self, chrom, name):
        
        self.chrom = chrom
        self.name = name
        self.n_exons = 0
        self.ranges = []
        self.intersecting_variants = []
        self.vtbl_indices = []
        self.n_variants = 0
        
    def add_variant(self,variant):
        
        self.intersecting_variants.append(variant)
        self.vtbl_indices.append(variant.vt_index)
        self.n_variants += 1
    
    ##run this method after you have added all variants
    def extract_vtbl(self,vtbl):
        
        if not len(self.vtbl_indices) > 0:
            
            raise ValueError("This transcript has no associated variants")
            
        else:

            self.vtbl_indices.sort()
            self.vtbl = vtbl[vtbl["CHROM"] ==\
                         self.chrom.encode('ascii')][[self.vtbl_indices]]
        
    def extract_haplotypes(self,phased_genotype_array_by_chrom):
        
        self.genos_phased =\
        phased_genotype_array_by_chrom[self.chrom][self.vtbl_indices]
        
        haplotypes_mother = self.genos_phased[:,0].to_haplotypes()
        haplotypes_father = self.genos_phased[:,1].to_haplotypes()
        
        self.haplotypes =\
        haplotypes_mother.concatenate(haplotypes_father, axis=1)
        
        self.variants_by_haplotype = {}
    
        for i in range(4):
        
            hap = self.haplotypes[:,i]
        
            if np.sum(hap) > 0:
        
                self.variants_by_haplotype[i] = hap
             
    def calculate_haplotype_based_length_changes(self):
        
        self.haplotype_length_changes = {}

        ##calculate net length change for each haplotype 
        
        for haplotype_index in range(self.haplotypes.shape[1]):
            changes = []
            haplotype = self.haplotypes[:, haplotype_index]
            
            ##first, get the index and identity of each variant
            
            for variant_index,genotype_index in enumerate(haplotype):
                
                if genotype_index == 0:
                    ##if variant is REF, length change is 0
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
