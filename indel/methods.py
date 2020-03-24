'''This module contains the indel methods.'''

import allel
import h5py
import numpy as np
import pandas as pd

#from indel import AffectedTranscript

class ZeroVariantError(Exception):
    pass

def check_shape(genotypes, vt):
    
    if genotypes.shape[0] == 0 or vt.shape[0] == 0:
        
        raise ZeroVariantError("all variants have been filtered out")

def read_calls(path, name):
    
    callset = h5py.File(path, mode='r')[name]

    genotypes = allel.GenotypeArray(callset["calldata/genotype"])

    gq = callset["calldata/GQ"][:]

    vt = allel.VariantChunkedTable(callset["variants"])[:]
    
    return(callset, genotypes, gq, vt)
    
def read_metadata(path_to_metadata, check_parental=False, sep="\t"):

    metadata = pd.read_csv(path_to_metadata, sep=sep)

    if check_parental:

        for i in range(2):

            if not metadata["role"][i] == "parent":

                raise ValueError(
                    "non-parental sample at the start, index {i}".format(
                        i=i))

        for i in range(2, len(metadata)):

            if not metadata["role"][i] == "progeny":

                raise ValueError("non-progeny sample with index {i}".format(
                    i=i))

    return metadata

def id_positions(name, vt, feature_table):
    
    try:
        all_positions = allel.SortedIndex(vt["POS"])

    except ValueError:
        print("Load one chromosome at a time")

    region_feature = feature_table[feature_table.seqid == name]

    positions, overlapped_features = all_positions.locate_intersection_ranges(
        region_feature.start,
        region_feature.end)
    
    if not np.sum(positions) > 0:
        
        raise ZeroVariantError("no positions overlap")
    
    return(positions, overlapped_features)
    
def extract_positions(genotypes, gq, vt, positions_bool):

    '''return genotypes for all exonic positions,
    associated genotype qualities, and a variant table
    of the exonic positions
    '''
    
    check_shape(genotypes, vt)

    genotypes = genotypes.subset(sel0=positions_bool)

    gq = gq[positions_bool, :]

    vt = vt[positions_bool]

    if not len(genotypes) == len(gq) == len(vt):

        raise ValueError(("Genotypes, genotype qualities, and "
                          + "variant table do not have the same length"))
        
    return(genotypes, gq, vt)
##is this really worth having as a separate function? I'll see.
    
def filter_quality(genotypes, gq, vt, min_gq=20, min_qd=15, min_mq=40, 
                   allowed_missing=0):

    '''return a set of genotypes and variants where genotypes below the
    minimum genotype quality (GQ) are masked, and sites below the minimum
    map quality (MQ), depth-adjusted quality score (QD), and missingness
    (allowed_missing) are filtered out
    '''

    check_shape(genotypes, vt)

    # Hide genotypes of poor quality, so they register as missing
    genotypes.mask = gq < min_gq

    # How many samples have to be called for a site to pass?
    num_present = gq.shape[1]

    min_called = num_present - allowed_missing

    # Get the list of sites with sufficient # of samples called
    missingness_filter = np.sum(genotypes.is_called(), axis=1) >= min_called

    # Filter the genotypes and variant table to keep only these sites
    genotypes_gq_filtered = genotypes.subset(sel0=missingness_filter)
    vt_gq_filtered = vt[missingness_filter]

    # Filtering by MQ and QD, which evaluates on the variant table

    filter_expression = '( (MQ >= {MQ}) & (QD >= {QD}) )'.format(
        MQ=min_mq, QD=min_qd)

    variant_qual_bool = vt_gq_filtered.eval(filter_expression)

    genotypes_qual_filtered = \
    genotypes_gq_filtered.subset(sel0=variant_qual_bool)

    vt_qual_filtered = vt_gq_filtered[variant_qual_bool]

    return(genotypes_qual_filtered, vt_qual_filtered)
    
def filter_on_parents(genotypes, vt):
    
    check_shape(genotypes, vt)
    
    # Remove sites with missing parents
    parents_called = np.all(genotypes[:, :2].is_called(), axis=1)

    genotypes_parents_called = genotypes.subset(sel0=parents_called)

    vt_parents_called = vt[parents_called]

    # Remove sites with non-segregating variants
    segregating = \
    genotypes_parents_called[:, :2].count_alleles().is_segregating()

    gt_filtered = genotypes_parents_called.subset(sel0=segregating)

    vt_filtered = vt_parents_called[segregating]

    return(gt_filtered, vt_filtered)
    
def filter_heterozygous_heterogametes(genotypes, vt, heterogametic, 
                                      error_tol=0):

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
    
    check_shape(genotypes, vt)

    het_errors = np.sum(genotypes[:, heterogametic].is_het(), axis=1)

    het_errors_bool = het_errors <= error_tol

    gt_het_error_filtered = genotypes.subset(sel0=het_errors_bool)

    vt_het_error_filtered = vt[het_errors_bool]
    
    return(gt_het_error_filtered, vt_het_error_filtered)
    
def id_mendelian_violations(genotypes, sex=False, parents_homo_progeny=None):
    
    if genotypes.shape[0] == 0:
        
        raise ZeroVariantError("all variants have been filtered out")
    
    if sex:
        
        if parents_homo_progeny is None:
            
            raise ValueError(
                    "must pass indices of parental samples " + \
                    "and homogametic samples")
            
        no_het_progeny = genotypes.subset(sel1=parents_homo_progeny)
        
        gt_violations = allel.mendel_errors(no_het_progeny[:, :2], 
                                            no_het_progeny[:, 2:])
        
    else:
        
        gt_violations = allel.mendel_errors(genotypes[:, :2], genotypes[:, 2:])
    
    site_violations = np.sum(gt_violations, axis=1)

    return site_violations

def remove_mendelian_violations(genotypes, vt, site_violations, 
                                permitted_violations=0):

    check_shape(genotypes, vt)
    
    mendel_bool = site_violations <= permitted_violations

    vt_mendel_filtered = vt[mendel_bool]

    genotypes_mendel_filtered = genotypes.subset(sel0=mendel_bool)
    
    return(genotypes_mendel_filtered, vt_mendel_filtered)

def phase_and_filter_parents(genotypes, vt, window=25):
    
    check_shape(genotypes, vt)

    phased = allel.phase_by_transmission(genotypes, window_size=window)

    parents_phased_bool = np.sum(phased.is_phased[:, 0:2], axis=1) == 2

    genotypes_phased = phased.subset(sel0 = parents_phased_bool)

    vt_phased = vt[parents_phased_bool]

    ##for now, I will not handle extracting the sites that could be phased 
    ##from the unphased genotypes
    
    return(genotypes_phased, vt_phased)
    
def filter_on_type(genotypes, vt, variant_type):
    
    check_shape(genotypes, vt)

    if variant_type == "SNP":

        variant_type_bool = vt["is_snp"]

    elif variant_type == "indel":

        variant_type_bool = (~ vt["is_snp"])

    else:

        raise ValueError("variant_type must be 'SNP' or 'indel'" +
                             variant_type)

    vt_filtered = vt[variant_type_bool]
    genotypes_filtered = genotypes.subset(sel0 = variant_type_bool)
    
    return(genotypes_filtered, vt_filtered)
    
def cross_workflow(name, callset_path, features, filter_on="indel", 
                   phased=False, sex=False,
                   hets=None, X_phasing=None, min_gq=20, min_qd=15, min_mq=40, 
                   allowed_missing=0):
    
    callset, genotypes, gq, vt = read_calls(callset_path, name)
    
    positions, overlapped_features = id_positions(name, vt, features)

    genotypes, gq, vt = extract_positions(genotypes, gq, vt, positions)
    
    genotypes, vt = filter_quality(genotypes, gq, vt, min_gq=min_gq, 
                                   min_qd=min_qd, min_mq=min_mq, 
                                   allowed_missing=allowed_missing)
    
    genotypes, vt = filter_on_parents(genotypes, vt)
    
    if sex:
        
        if hets is None:
            
            raise ValueError("must pass indices of heterogametic samples")
            
        else:
            
            genotypes, vt = filter_heterozygous_heterogametes(genotypes, vt, 
                                                              hets)
    
        mendelian_violations = id_mendelian_violations(genotypes, sex=True, 
                                                       parents_homo_progeny=\
                                                       X_phasing)

    else:
        
        mendelian_violations = id_mendelian_violations(genotypes)
        
    genotypes, vt = remove_mendelian_violations(genotypes, vt, 
                                                mendelian_violations)
    
    if phased:
        
        genotypes, vt = phase_and_filter_parents(genotypes, vt)
    
    genotypes, vt = filter_on_type(genotypes, vt, variant_type=filter_on)
    
    return(genotypes, vt)

def wild_workflow(name, callset_path, features, filter_on="indel"):
    
    callset, genotypes, gq, vt = read_calls(callset_path, name)
    
    positions, overlapped_features = id_positions(name, vt, features)

    genotypes, gq, vt = extract_positions(genotypes, gq, vt, positions)
        
    genotypes, vt = filter_on_type(genotypes, vt, variant_type=filter_on)
    
    return(genotypes, vt)

'''def run_workflow(chrom, name, variant_table, positions,
                 genotypes, feature_table, feature_type):
    
    if feature_type == "CDS":
        
        field = "Parent"
        
    elif feature_type == "mRNA":
        
        field = "ID"

    transcript = AffectedTranscript(chrom=chrom,
                                    name=name)

    cds_frame = feature_table[feature_table[field] == name]

    transcript.ranges = np.stack((cds_frame["start"], cds_frame["end"]),
                                 axis=1)

    transcript.positions =\
        positions.intersect_ranges(transcript.ranges[:, 0],
                                   transcript.ranges[:, 1]).values
                                   
    indices =\
    np.array([positions.locate_key(key) for key in transcript.positions])

    transcript.add_indices(indices)

    if len(transcript.indices) > 0:

        transcript.extract_vtbl(variant_table)

        transcript.extract_genos(genotypes)

    return transcript'''
