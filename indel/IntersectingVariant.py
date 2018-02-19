import numpy as np

class IntersectingVariant():
    
    def __init__(self,chrom_name,position):
        self.chrom = chrom_name
        self.position = position
        self.expression =\
        "(POS == {position}) & (CHROM == {name})".format\
        (position=self.position, name=self.chrom.encode('ascii'))
        self.length_change = [0]
        
    def add_vtbl_record(self,variant_table):
        ##pull out the variant table row describing this position
        self.record = variant_table.query(self.expression)
        
        ##record its position in the variant table-- crucial for filtering genotypes
        index = np.where(variant_table["POS"] == self.position)[0][0]
        
        assert index >= 0, "Can't have a negative variant table index"
        
        self.vt_index = index
        
        ##for every ALT allele, record the length change
        self.num_alleles = self.record["num_alleles"][0]
        
        for i in range(self.num_alleles-1):
            length_change = len(self.record["ALT"][0][i]) - len(self.record["REF"][0])
            self.length_change.append(length_change)
        
    def add_genos(self,genotypes,variant_table):
        ##attach the corresponding genotype vector
        self.genos = genotypes[variant_table.eval(self.expression)]               
