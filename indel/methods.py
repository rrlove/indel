'''This module contains the indel methods.'''

import numpy as np

from indel import AffectedTranscript

def run_workflow(chrom, name, variant_table, positions,
                 genotypes, feature_table):

    transcript = AffectedTranscript(chrom=chrom,
                                    name=name)

    CDS_frame = feature_table[feature_table["Parent"] == name]

    transcript.ranges = np.stack((CDS_frame["start"], CDS_frame["end"]),
                                 axis=1)

    transcript.positions =\
        positions.intersect_ranges(transcript.ranges[:, 0],
                                   transcript.ranges[:, 1]).values

    transcript.indices =\
        np.array([positions.locate_key(key) for key in transcript.positions])

    if len(transcript.indices) > 0:

        transcript.extract_vtbl(variant_table)

        transcript.extract_genos(genotypes)

    return transcript
