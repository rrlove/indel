import numpy as np


class AffectedTranscript():

    def __init__(self, chrom, name):

        self.chrom = chrom
        self.genos = None
        self.haplotypes = None
        self.indices = None
        self.name = name
        self.positions = None
        self.ranges = None
        self.vtbl = None

    @property
    def n_exons(self):

        return len(self.ranges)
    
    @property
    def n_indels(self):
        
        return np.sum(~ self.vtbl["is_snp"])
    
    def add_indices(self, indices):
        
        self.indices = indices
        self.indices.sort()

    def extract_vtbl(self, vtbl):

        if self.indices is None:

            raise ValueError("This transcript has no associated variants")

        self.indices.sort()

        self.vtbl = vtbl[vtbl["CHROM"] ==\
                         self.chrom.encode('ascii')][[self.indices]]

    def extract_genos(self, genos):

        if self.indices is None:

            raise ValueError("This transcript has no associated variants")

        self.indices.sort()

        self.genos = genos[[self.indices]]

    def extract_haplotypes(self):

        # Check that parental genotypes are phased

        if not (np.sum(np.sum(self.genos.is_phased[:, 0:2], axis=1) == 2) ==
                len(self.genos)):

            raise ValueError("Genotypes are not completely phased")

        haplotypes_mother = self.genos[:, 0].to_haplotypes()
        haplotypes_father = self.genos[:, 1].to_haplotypes()

        self.haplotypes =\
        haplotypes_mother.concatenate(haplotypes_father, axis=1)
