import allel
import h5py
import numpy as np
import pandas as pd


class Chrom():

    def __init__(self, name, path_to_hdf5):
        self.path = path_to_hdf5
        self.name = name
        self.callset = None
        self.genotypes = None
        self.gq = None

    def read_calls(self):
        self.callset = h5py.File(self.path, mode='r')[self.name]

        self.genotypes =\
        allel.GenotypeArray(self.callset["calldata/genotype"])

        self.gq = self.callset["calldata/GQ"][:]

        self.vt = allel.VariantChunkedTable(self.callset["variants"])[:]

    def read_metadata(self, path_to_metadata):
        self.metadata = pd.read_csv(path_to_metadata, sep="\t")

        if not len(self.metadata) == self.callset["samples"].shape[0]:

            raise ValueError(("Metadata and callset sample lists "
                              + "are different lengths"))

    def check_metadata_parental(self):
        '''for any task requiring parental samples
        to be at the beginning of the array
        '''

        for i in range(2):

            if not self.metadata["role"][i] == "parent":

                raise ValueError(
                    "non-parental sample at the start, index {i}".format(
                        i=i))

        for i in range(2, len(self.metadata)):

            if not self.metadata["role"][i] == "progeny":

                raise ValueError("non-progeny sample with index {i}".format(
                    i=i))

    def id_positions(self, feature_table):

        '''take a feature table filtered down to entries of interest,
        return a boolean of positions in the callset that overlap those
        features. also return a boolean of features in the feature table
        overlapped by the callset. This latter object will be useful for
        transcript identification (or will it?
        I'll have filtered out some of the positions...)
        '''

        try:
            self.all_positions =\
            allel.SortedIndex(self.callset["variants"]["POS"])

        except ValueError:
            print("Load one chromosome at a time")

        region_feature = feature_table[feature_table.seqid == self.name]

        self.exonic_positions, self.overlapped_features = \
        self.all_positions.locate_intersection_ranges(
            region_feature.start,
            region_feature.end)

    def extract_exonic(self):

        '''return genotypes for all exonic positions,
        associated genotype qualities, and a variant table
        of the exonic positions
        '''

        self.genotypes = self.genotypes.subset(sel0=self.exonic_positions)

        self.gq = self.gq[self.exonic_positions, :]

        self.vt = self.vt[self.exonic_positions]

        if not len(self.genotypes) == len(self.gq) == len(self.vt):

            raise ValueError(("Genotypes, genotype qualities, and "
                              + "variant table do not have the same length"))

    def filter_gq_mq_qd(self,
                        min_gq=20,
                        min_qd=15,
                        min_mq=40,
                        allowed_missing=0):
        '''return a set of genotypes and variants where genotypes below the
        minimum genotype quality (GQ) are masked, and sites below the minimum
        map quality (MQ), depth-adjusted quality score (QD), and missingness
        (allowed_missing) are filtered out
        '''

        # Hide genotypes of poor quality, so they register as missing
        self.genotypes.mask = self.gq < min_gq

        # How many samples have to be called for a site to pass?
        self.num_present = self.gq.shape[1]

        min_called = self.num_present - allowed_missing

        # Get the list of sites with sufficient # of samples called
        self.missingness_filter =\
        np.sum(self.genotypes.is_called(), axis=1) >= min_called

        # Filter the genotypes and variant table to keep only these sites
        self.genotypes_gq_filtered =\
        self.genotypes.subset(sel0=self.missingness_filter)
        self.vtbl_gq_filtered = self.vt[self.missingness_filter]

        # Filtering by MQ and QD, which evaluates on the variant table

        self.filter_expression = '( (MQ >= {MQ}) & (QD >= {QD}) )'.format(
            MQ=min_mq, QD=min_qd)

        self.variant_qual_bool =\
        self.vtbl_gq_filtered.eval(self.filter_expression)

        self.gt_qual_filtered =\
        self.genotypes_gq_filtered.subset(sel0=self.variant_qual_bool)

        self.vt_qual_filtered = self.vtbl_gq_filtered[self.variant_qual_bool]

    def filter_on_parents(self):

        self.check_metadata_parental()

        # Remove sites with missing parents
        self.parents_called =\
        np.all(self.gt_qual_filtered[:, :2].is_called(), axis=1)

        self.gt_parents_called =\
        self.gt_qual_filtered.subset(sel0=self.parents_called)

        self.vt_parents_called = self.vt_qual_filtered[self.parents_called]

        # Remove sites with non-segregating variants
        self.segregating =\
        self.gt_parents_called[:, :2].count_alleles().is_segregating()

        self.gt_filtered = self.gt_parents_called.subset(sel0=self.segregating)

        self.vt_filtered = self.vt_parents_called[self.segregating]

    def id_autosomal_mendelian_violations(self):

        self.check_metadata_parental()

        self.auto_gt_violations = \
        allel.mendel_errors(self.gt_filtered[:, :2], self.gt_filtered[:, 2:])
        self.auto_site_violations = np.sum(self.auto_gt_violations, axis=1)

    def filter_heterozygous_heterogametes(self, heterogametic, error_tol=0):

        '''This method should only be run on the sex chromosome. It is part of
        the filtering process with id_mendelian_violations_sex.

        Because the heterogametic sex only has one
        copy of each sex chromosome, any heterozygous genotypes on the sex
        chromosome must be sequencing errors (in this representation, the
        genotype is "doubled" on the sex chromosome for the heterogametic sex.)

        This method takes a list of heterozygous samples and returns filtered
        genotypes and variant tables where sites with more errors than the
        threshold have been removed. By default, the threshold is 0. I
        strongly discourage changing this threshold, especially without
        specifically filtering out errors in the heterogametic PARENT.
        Omitting that could throw off the identification of Mendelian
        violations on that chromosome.

        By using this method together with ID_Mendelian_violations_sex,
        you should be able to check every sample on the sex chromosome,
        whether heterogametic or not.
        '''

        self.het_errors =\
        np.sum(self.gt_filtered[:, heterogametic].is_het(), axis=1)

        self.het_errors_bool = self.het_errors <= error_tol

        self.gt_het_error_filtered =\
        self.gt_filtered.subset(sel0=self.het_errors_bool)

        self.vt_het_error_filtered = self.vt_filtered[self.het_errors_bool]

    def id_mendelian_violations_sex(self, parents_homo_progeny):
        '''parents_homo_progeny is an array holding the indices of the parental
        samples plus the samples homogametic for the sex chromosome (in the
        case of Anopheles, females)
        heterogametic is an array holding the indices of all heterogametic
        samples
        '''
        # Check that parental samples are together at beginning of metadata
        self.check_metadata_parental()

        '''To force running filter_heterozygous_heterogametes first,
        use those output objects as input
        '''

        self.no_het_progeny =\
        self.gt_het_error_filtered.subset(sel1=parents_homo_progeny)

        self.sex_gt_violations =\
        allel.mendel_errors(
            self.no_het_progeny[:, :2], self.no_het_progeny[:, 2:])

        self.sex_site_violations = np.sum(self.sex_gt_violations, axis=1)

    def remove_mendelian_violations(self, permitted_violations=0, sex=False):

        if sex:

            self.mendel_bool =\
            self.sex_site_violations <= permitted_violations

            self.vt_mendel_filtered =\
            self.vt_het_error_filtered[self.mendel_bool]

            self.gt_mendel_filtered =\
            self.gt_het_error_filtered.subset(sel0=self.mendel_bool)

        elif not sex:

            self.mendel_bool =\
            self.auto_site_violations <= permitted_violations

            self.vt_mendel_filtered = self.vt_filtered[self.mendel_bool]

            self.gt_mendel_filtered =\
            self.gt_filtered.subset(sel0=self.mendel_bool)

    def phase_and_filter_parents(self, permitted_nonphased=0, window=25):

        self.check_metadata_parental()

        self.phased =\
        allel.phase_by_transmission(
            self.gt_mendel_filtered, window_size=window)

        if not self.num_present == self.gt_mendel_filtered.shape[1]:

            raise ValueError(
                ("You have a different number of samples"
                 + " from when you started!"))

        self.parents_phased_bool =\
        np.sum(self.phased.is_phased[:, 0:2], axis=1) == 2

        self.phased_genos = self.phased.subset(sel0=self.parents_phased_bool)

        self.unphased_genos_at_phased_site = \
        self.gt_mendel_filtered.subset(sel0=self.parents_phased_bool)

        self.vt_phased = self.vt_mendel_filtered[self.parents_phased_bool]

        if permitted_nonphased < len(self.metadata):

            self.phased_bool = (np.sum(self.phased_genos.is_phased, axis=1) >=
                                (self.num_present - permitted_nonphased))

            self.phased_genos = self.phased_genos.subset(sel0=self.phased_bool)

            self.unphased_genos_at_phased_site = \
            self.gt_mendel_filtered.subset(sel0=self.phased_bool)

            self.vt_phased = self.vt_mendel_filtered[self.phased_bool]
    