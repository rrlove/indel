'''This module contains the indel methods.'''

import numpy as np

from indel import AffectedTranscript

def run_workflow(chrom, name, variant_table, positions,
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

    return transcript
