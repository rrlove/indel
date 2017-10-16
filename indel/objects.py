class Chrom():
    
    def __init__(self,name,path_to_hdf5):
        self.name = str(name)
        self.path = path_to_hdf5
        
    def read_calls(self):
        self.callset = h5py.File(self.path,mode='r')
        
    def read_metadata(self,path_to_metadata):
        self.metadata = pd.read_csv(path_to_metadata, sep='\t')
        
        assert len(self.metadata) == self.callset[self.name]["samples"].shape[0] \
        , "Metadata and callset sample lists are different lengths"
        
    def check_metadata_parental(self):
        ##for any task requiring parental samples to be at the beginning of the array
        
        assert self.metadata["role"][0] == "parent", "non-parental sample at beginning, index 0"
        assert self.metadata["role"][1] == "parent", "non-parental sample at beginning, index 1"
        
        for i in range(2,len(self.metadata)):
            
            assert self.metadata["role"][i] == "progeny", "non-progeny sample with index " + str(i)
        
    def ID_positions(self,feature_table):
        
        '''take a feature table filtered down to entries of interest, 
        return a boolean of positions in the callset that overlap those features.
        also return a boolean of features in the feature table overlapped by the callset.
        This latter object will be useful for transcript identification (or will it?
        I'll have filtered out some of the positions...)'''
        
        try:
            self.allPositions = allel.SortedIndex(self.callset[self.name]["variants"]["POS"])
        except ValueError:
            print("Load one chromosome at a time")
            
        self.exonicPositions, self.overlappedFeatures = \
        self.allPositions.locate_intersection_ranges(
            feature_table[feature_table.seqid == self.name.encode('ascii')].start,
            feature_table[feature_table.seqid == self.name.encode('ascii')].end)
        
    def extract_exonic(self):
        
        '''return genotypes for all exonic positions, associated genotype qualities,
        and a variant table of the exonic positions'''
        
        self.genotypes = allel.GenotypeArray(self.callset[self.name]["calldata/genotype"]).\
        subset(sel0=self.exonicPositions)
        self.GQ = self.callset[self.name]["calldata/GQ"][self.exonicPositions,:]
        
        self.vt = allel.VariantChunkedTable(self.callset[self.name]["variants"])[:]\
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
        self.gt_violations = allel.mendel_errors(self.gt_filtered[:,:2],self.gt_filtered[:,2:])
        self.site_violations = np.sum(self.gt_violations, axis=1)
    
        
