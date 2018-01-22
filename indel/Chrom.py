import allel
import h5py
import numpy as np
import pandas as pd

class Chrom():
    
    def __init__(self,name,path_to_hdf5):
        self.name = str(name)
        self.path = path_to_hdf5
        
    def read_calls(self):
        ##because a callset always has the same name, and is processed chrom by chrom, I can simplify the syntax
        ##and change all references to self.callset[self.name] to self.callset
        self.raw_callset = h5py.File(self.path,mode='r')
        self.callset = self.raw_callset[self.name]
        
    def read_metadata(self,path_to_metadata):
        self.metadata = pd.read_csv(path_to_metadata, sep='\t')
        
        assert len(self.metadata) == self.callset["samples"].shape[0] \
        , "Metadata and callset sample lists are different lengths"
        
    def check_metadata_parental(self):
        ##for any task requiring parental samples to be at the beginning of the array
        
        assert self.metadata["role"][0] == \
        "parent", "non-parental sample at beginning, index 0"
        assert self.metadata["role"][1] == \
        "parent", "non-parental sample at beginning, index 1"
        
        for i in range(2,len(self.metadata)):
            
            assert self.metadata["role"][i] == \
            "progeny", "non-progeny sample with index " + str(i)
        
    def ID_positions(self,feature_table):
        
        '''take a feature table filtered down to entries of interest, 
        return a boolean of positions in the callset that overlap those features.
        also return a boolean of features in the feature table overlapped by the callset.
        This latter object will be useful for transcript identification (or will it?
        I'll have filtered out some of the positions...)'''
        
        try:
            self.allPositions = allel.SortedIndex(self.callset["variants"]["POS"])
        except ValueError:
            print("Load one chromosome at a time")
            
        self.exonicPositions, self.overlappedFeatures = \
        self.allPositions.locate_intersection_ranges(
            feature_table[feature_table.seqid == self.name].start,
            feature_table[feature_table.seqid == self.name].end)
        
    def extract_exonic(self):
        
        '''return genotypes for all exonic positions, associated genotype qualities,
        and a variant table of the exonic positions'''        
        
        self.genotypes = allel.GenotypeArray(self.callset["calldata/genotype"]).\
        subset(sel0=self.exonicPositions)
        
        self.GQ = self.callset["calldata/GQ"][self.exonicPositions,:]

        self.vt = allel.VariantChunkedTable(self.callset["variants"])[:]\
        [self.exonicPositions]
        
        assert len(self.genotypes) == len(self.GQ) == len(self.vt)\
        , "Genotypes, genotype qualities, and variant table do not have the same length"
        
        ##write test to make sure vt doesn't return any of the various indexing errors above/below
        
    def filter_GQ_MQ_QD(self,min_GQ=20,min_QD=15,min_MQ=40,allowed_missing=0):
        '''return a set of genotypes and variants filtered for minimum # of present
        genotypes of a certain quality (GQ, allowed_missing), minimum map quality (MQ), and 
        minimum depth-adjusted quality score (QD)
        '''
        
        self.num_present = self.GQ.shape[1]
        
        self.gq_filter = self.GQ >= min_GQ
        self.gq_bool = np.sum(self.gq_filter, axis=1) >= (self.num_present - allowed_missing)
        
        ##filtering by GQ, which evaluates based on the genotype quality object
        self.genotypes_GQ_filtered = self.genotypes.subset(sel0=self.gq_bool)
        self.vtbl_GQ_filtered = self.vt[self.gq_bool]
        
        ##filtering by MQ and GQ, which evaluates based on the variant table object
        
        self.filter_expression = '( (MQ >= {MQ}) & (QD >= {QD}) )'.format(
            MQ = min_MQ, QD = min_QD)
        
        self.variant_qual_bool = self.vtbl_GQ_filtered.eval(self.filter_expression)
        
        self.gt_qual_filtered = self.genotypes_GQ_filtered.subset(sel0=self.variant_qual_bool)
        
        self.vt_qual_filtered = self.vtbl_GQ_filtered[self.variant_qual_bool]
        
    def filter_on_parents(self):
        ##remove any sites with missing or non-segregating parents
        
        self.check_metadata_parental()
        
        ##remove sites with missing parents
        self.parents_called = np.all(self.gt_qual_filtered[:,:2].is_called(), axis=1)
        self.gt_parents_called = self.gt_qual_filtered.subset(sel0=self.parents_called)
        self.vt_parents_called = self.vt_qual_filtered[self.parents_called]
        
        ##remove sites with non-segregating variants
        self.segregating = self.gt_parents_called[:,:2].count_alleles().is_segregating()
        self.gt_filtered = self.gt_parents_called.subset(sel0=self.segregating)
        self.vt_filtered = self.vt_parents_called[self.segregating]
        ##I could reduce intermediate objects by stacking the vt filtering here...
        
    def ID_autosomal_Mendelian_violations(self):
        
        self.check_metadata_parental()
        
        ##now that we've checked that, actually find the violations
        self.auto_gt_violations = \
        allel.mendel_errors(self.gt_filtered[:,:2],self.gt_filtered[:,2:])
        self.auto_site_violations = np.sum(self.auto_gt_violations, axis=1) 
        
    def filter_heterozygous_heterogametes(self,heterogametic,error_tolerance=0):
        
        '''This method should only be run on the sex chromosome. It is part of the filtering
        process with ID_Mendelian_violations_sex. 
        
        Because the heterogametic sex only has one
        copy of each sex chromosome, any heterozygous genotypes on the sex chromosome must be
        sequencing errors (in this representation, the genotype is "doubled" on the sex 
        chromosome for the heterogametic sex.)
        
        This method takes a list of heterozygous samples and returns filtered genotypes and 
        variant tables where sites with more errors than the threshold have been removed. By
        default, the threshold is 0. I strongly discourage changing this threshold, especially
        without specifically filtering out errors in the heterogametic PARENT. Omitting that
        could throw off the identification of Mendelian violations on that chromosome.
        
        By using this method together with ID_Mendelian_violations_sex, you should be able to 
        check every sample on the sex chromosome, whether heterogametic or not.'''
        
        self.het_errors = np.sum(self.gt_filtered[:,heterogametic].is_het(), axis=1)
        self.het_errors_bool = self.het_errors <= error_tolerance
        
        self.gt_het_error_filtered = self.gt_filtered.subset(sel0=self.het_errors_bool)
        self.vt_het_error_filtered = self.vt_filtered[self.het_errors_bool]
    
    def ID_Mendelian_violations_sex(self,parents_homo_progeny): 
        '''parents_homo_progeny is an array holding the indices of the parental samples plus 
        the samples homogametic for the sex chromosome (in the case of Anopheles, females)
        heterogametic is an array holding the indices of all heterogametic samples
        '''
        ##parental samples have to be together at the beginning of the array, so check
        self.check_metadata_parental()
        
        ##to force running filter_heterozygous_heterogametes first, use those output objects as input
        
        self.no_het_progeny = self.gt_het_error_filtered.subset(sel1=parents_homo_progeny)
        
        self.sex_gt_violations = \
        allel.mendel_errors(self.no_het_progeny[:,:2],self.no_het_progeny[:,2:])
        
        self.sex_site_violations = np.sum(self.sex_gt_violations, axis=1)
        
    def remove_Mendelian_violations(self,permitted_violations=0,sex=False):
        
        if sex is True:
            
            self.Mendel_bool = self.sex_site_violations <= permitted_violations
            self.vt_Mendel_filtered = self.vt_het_error_filtered[self.Mendel_bool]
            self.gt_Mendel_filtered = self.gt_het_error_filtered.subset(sel0=self.Mendel_bool)
            
        elif sex is False:
            
            self.Mendel_bool = self.auto_site_violations <= permitted_violations
            self.vt_Mendel_filtered = self.vt_filtered[self.Mendel_bool]
            self.gt_Mendel_filtered = self.gt_filtered.subset(sel0=self.Mendel_bool)
    
    def phase_and_filter(self,permitted_nonphased=0,window=25):
        
        self.check_metadata_parental()
        
        self.phased = allel.phase_by_transmission(self.gt_Mendel_filtered, window_size=window)
        assert self.num_present == self.gt_Mendel_filtered.shape[1], \
        "You have a different number of samples from when you started!"
        
        self.phased_Bool = np.sum(self.phased.is_phased, axis=1) >= \
        (self.num_present - permitted_nonphased)
        
        self.phased_genos = self.phased.subset(sel0=self.phased_Bool)
        self.unphased_genos_at_phased_site = \
        self.gt_Mendel_filtered.subset(sel0=self.phased_Bool)
        self.vt_phased = self.vt_Mendel_filtered[self.phased_Bool]

    
