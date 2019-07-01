#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:58:25 2019

@author: beccalove
"""

import unittest

import numpy as np
import numpy.testing as npt

import allel
import indel

class TestWorkflow(unittest.TestCase):

    def setUp(self):

        ##mock variant table

        records = [[b'3R', 101, b'CG', [b'C', '', ''], 2],
                   [b'3R', 180, b'A', [b'AA', '', ''], 2],
                   [b'3R', 207, b'A', [b'AT', b'ATT', b'ATTT'], 4],
                   [b'3R', 220, b'T', [b'TC', '', ''], 2],
                   [b'3R', 284, b'A', [b'ACAG', b'ACAGCAG', ''], 3],
                   [b'3R', 316, b'GAAA', [b'G', '', ''], 2],
                   [b'3R', 347, b'TTGTGTG', [b'T', '', ''], 2],
                   [b'3R', 1091, b'G', [b'GA', b'GAAGA', ''], 3],
                   [b'3R', 74807, b'ATATA', [b'A', b'ATA', ''], 3],
                   [b'3R', 16045837, b'CA', [b'C', b'CAA', b'CAAA'], 4],
                   [b'3R', 28475632, b'T', [b'TT', '', ''], 2]]

        dtype = [('CHROM', 'S4'),
                 ('POS', 'u4'),
                 ('REF', 'S100'),
                 ('ALT', ('S100', 3)),
                 ('num_alleles', 'uint8')]

        self.variant_table = allel.VariantTable(
            records, dtype=dtype, index=('CHROM', 'POS'))

        ##mock genotypes

        self.genotypes = allel.GenotypeArray([
            [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[1, 1], [0, 1], [0, 1], [0, 1], [0, 0], [0, 0]],
            [[0, 1], [2, 2], [0, 2], [0, 2], [0, 2], [0, 2]],
            [[1, 1], [0, 0], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[2, 2], [0, 1], [0, 2], [1, 2], [0, 2], [1, 2]],
            [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[1, 1], [0, 0], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
            [[0, 0], [1, 2], [0, 1], [0, 2], [0, 2], [0, 2]],
            [[0, 1], [0, 1], [1, 1], [0, 0], [0, 1], [0, 0]]
            ], dtype='i1')

        ##mock feature table
        feature_table_data = [
            ['3R', 'DB', 'CDS', 95, 120, -1, '+', 0,
             'gene3-PA', 'gene3-RA'],
            ['3R', 'DB', 'CDS', 206, 296, -1, '+', 1,
             'gene3-PA', 'gene3-RA'],
            ['3R', 'DB', 'CDS', 310, 363, -1, '+', 0,
             'gene3-PA', 'gene3-RA'],
            ['3R', 'DB', 'CDS', 15780, 15809, -1, '+',
             0, 'gene24-PA', 'gene24-RA'],
            ['3R', 'DB', 'CDS', 21075362, 21075460, -1,
             '+', 0, 'gene436-PA', 'gene436-RA']
            ]

        feature_table_dtype = [
            ('seqid', 'S4'),
            ('source', 'S2'),
            ('type', 'S15'),
            ('start', int),
            ('end', int),
            ('score', '<f8'),
            ('strand', 'S1'),
            ('phase', '<i8'),
            ('ID', 'S15'),
            ('Parent', 'S15'),
            ]

        self.feature_table = \
        allel.FeatureTable(feature_table_data, dtype=feature_table_dtype)

        self.chrom = '3R'
        self.name = b'gene3-RA'

        self.positions = allel.SortedIndex(self.variant_table["POS"])

        self.affected_transcript = indel.run_workflow(self.chrom,
                                                      self.name,
                                                      self.variant_table,
                                                      self.positions,
                                                      self.genotypes,
                                                      self.feature_table)

    def test_ranges(self):

        test_ranges = np.array([[95, 120], [206, 296], [310, 363]])

        npt.assert_array_equal(test_ranges, self.affected_transcript.ranges)

    def test_positions(self):

        test_positions = np.array([101, 207, 220, 284, 316, 347])

        npt.assert_array_equal(test_positions,
                               self.affected_transcript.positions)

    def test_indices(self):

        test_indices = np.array([0, 2, 3, 4, 5, 6])

        npt.assert_array_equal(test_indices, self.affected_transcript.indices)

    def test_n_exons(self):

        self.assertEqual(3, self.affected_transcript.n_exons)
