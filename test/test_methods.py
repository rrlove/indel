#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:58:25 2019

@author: beccalove
"""

import unittest

import numpy as np
import numpy.testing as npt
import pandas as pd

import allel
import indel

class MethodsTestCase(unittest.TestCase):

    def setUp(self):

        ##mock genotype qualities
        self.gq = np.array([[27, 23, 19, 21, 15, 24],
                       [28, 21, 23, 20, 29, 25],
                       [32, 20, 29, 31, 30, 20],
                       [21, 27, 29, 32, 30, 33],
                       [22, 23, 29, 30, 27, 21],
                       [21, 24, 26, 27, 26, 25],
                       [20, 32, 29, 31, 30, 28],
                       [27, 26, 25, 24, 30, 31],
                       [25, 24, 32, 31, 28, 25]])

        ##one variant to get filtered out by GQ
        ##one variant to get filtered out as non-CDS
        ##one variant to get filtered out as Mendelian error
        ##one variant to get filtered out by MQ
        ##one variant to get filtered out by QD
        ##one variant to get filtered out as non-seg
        ##one variant to get filtered out with a missing parent
        ##two variants to rule them all, and in the darkness bind them

        ##mock genotypes
        self.genotypes =\
        allel.GenotypeArray(
            [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
             [[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
             ##the variant below is filtered out as non-segregating
             [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
             [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
             [[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             ##the variant below is filtered out as a Mendelian error
             [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]],
             ##the variant below is filtered out as missing parent
             [[-1, -1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]
             ], dtype='i1')

        ##mock variant table
        ##default min QD 15
        ##default min MQ 40
        ##variants filtered out:1 (GQ),
        ##2 (QD),
        ##3 (non-CDS),
        ##4 (non-seg),
        ##5 (MQ),
        ##8 (Mendelian error),
        ##9 (missing parent)
        records = [(b'2L', 20, b'A', b'ATC', 15, 35.7, 45, (6, 6)),
                   ##the variant below gets filtered out for QD
                   (b'2L', 53, b'CG', b'C', 18, 14, 45, (9, 3)),
                   ##the variant below is filtered out as non-CDS
                   (b'2L', 207, b'ATATA', b'A', 13, 29.2, 41, (3, 9)),
                   (b'2L', 602, b'TCTCT', b'T', 27, 16, 41, (0, 12)),
                   ##the variant below is filtered out for MQ
                   (b'2L', 754, b'TCTGT', b'TC', 18, 33.2, 29, (6, 6)),
                   (b'2L', 1090, b'G', b'GTTT', 17, 31.4, 43, (6, 6)),
                   (b'2L', 1114, b'GCAT', b'G', 20, 37.2, 44, (6, 6)),
                   (b'2L', 1129, b'T', b'TTTT', 32, 15, 40, (3, 9)),
                   (b'2L', 1150, b'AAA', b'A', 30, 20, 45, (0, 10))
                  ]

        dtype = [('CHROM', 'S4'),
                 ('POS', 'u4'),
                 ('REF', 'S100'),
                 ('ALT', 'S100'),
                 ('DP', int),
                 ('QD', float),
                 ('MQ', int),
                 ('AC', (int, 2))]

        self.vt = allel.VariantTable(records,
                                dtype=dtype, index=('CHROM', 'POS'))

        self.callset = {"variants": self.vt, 
                        "calldata/genotype": self.genotypes,
                        "calldata/GQ": self.gq}

        ##mock feature table
        feature_table_data = [
            ["2L", "DB", b'gene', 15, 2000, -1, "+", -1, "gene1", "."],
            ["2L", "DB", b'CDS', 17, 77, -1, "+", 0, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 569, 803, -1, "+", 2, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 1001, 1166, -1, "+", 2, "gene1-PA", 
             "gene1-RA"],
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

        feature_table_names = tuple(t[0] for t in feature_table_dtype)

        self.test_feature_table = \
        allel.FeatureTable(feature_table_data, names=feature_table_names)
        
        self.test_cds = \
        self.test_feature_table[self.test_feature_table.type == b'CDS']

        ##mock metadata
        sample_names = pd.Series(["A", "B", "C", "D", "E", "F"])
        sample_roles = \
        pd.Series(["parent", "parent", "progeny",
                   "progeny", "progeny", "progeny"])
        sample_sexes = pd.Series(["F", "M", "F", "M", "F", "M"])
        
        metadata = \
        pd.DataFrame({'samples': sample_names, 'role': sample_roles,\
                      'sex': sample_sexes})
        
        self.name = "2L"
        
    def test_id_positions(self):

        actual_positions, actual_features = indel.id_positions(self.name, 
                                                               self.vt, 
                                                               self.test_cds)

        expected_positions = [True, True, False, True, True,
                              True, True, True, True]

        expected_features = [True, True, True]

        npt.assert_equal(expected_features, actual_features)
        npt.assert_equal(expected_positions, actual_positions)
        
    def test_extract_positions(self):

        ##these differ from the original test data
        ##because one variant doesn't overlap a CDS
        expected_exonic_gt_shape = (8, 6, 2)
        expected_exonic_gq_shape = (8, 6)
        expected_exonic_vt_shape = (8,)
        
        positions, _ = indel.id_positions(self.name, self.vt, self.test_cds)
        
        genotypes, gq, vt = indel.extract_positions(self.genotypes, self.gq, 
                                              self.vt, positions)

        self.assertEqual(expected_exonic_gt_shape,
                         genotypes.shape)
        
        self.assertEqual(expected_exonic_gq_shape, gq.shape)
        
        self.assertEqual(expected_exonic_vt_shape, vt.shape)
    
    def test_quality_filtering_basic(self):
        #this test and the ones below use allel.GenotypeArray(self.genotypes)
        #to create a copy. otherwise, the original genotypes are masked
        #which has inappropriate downstream effects
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt)
        
        self.assertEqual(genotypes.shape, (5, 6, 2))
        self.assertEqual(vt.shape, (5,))

    def test_quality_filtering_gq(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, min_gq=15)
        
        self.assertEqual(genotypes.shape, (6, 6, 2))
        self.assertEqual(vt.shape, (6,))
        
    def test_quality_filtering_mq(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, min_mq=25)
        
        self.assertEqual(genotypes.shape, (6, 6, 2))
        self.assertEqual(vt.shape, (6,))
        
    def test_quality_filtering_qd(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, min_qd=10)
        
        self.assertEqual(genotypes.shape, (6, 6, 2))
        self.assertEqual(vt.shape, (6,))
        
    def test_quality_filtering_missingness(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, allowed_missing=2)
        
        ##it's now seven instead of six because low quality genotypes are
        ##masked and therefore appear as missing
        ##so the low-gq variant is now not removed
        self.assertEqual(genotypes.shape, (7, 6, 2))
        self.assertEqual(vt.shape, (7,))
        
    def test_quality_filtering_combo_one(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, min_gq=15, min_mq=20)
        
        self.assertEqual(genotypes.shape, (7, 6, 2))
        self.assertEqual(vt.shape, (7,))
        
    def test_quality_filtering_combo_two(self):
        
        genotypes, vt = indel.filter_quality(
                allel.GenotypeArray(self.genotypes), 
                self.gq, self.vt, min_gq=25, 
                                       allowed_missing=2)
        
        self.assertEqual(genotypes.shape, (6, 6, 2))
        self.assertEqual(vt.shape, (6,))
        
    def test_filter_on_parents(self):
        
        genotypes, vt = indel.filter_on_parents(self.genotypes, self.vt)

        self.assertEqual(genotypes.shape, (7, 6, 2))
        self.assertEqual(vt.shape, (7, ))

        expected_parental_filter_genos = \
        [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
         [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
         [[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
         [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
         [[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
         [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
         [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]]
         ]

        npt.assert_array_equal(expected_parental_filter_genos,
                               genotypes)
        
    def test_filter_heterozygous_heterogametes(self):

        '''I reuse the original mock genotypes and variant filters and
        "pretend" 2L is a sex chromosome'''
        
        heterogametic = [2, 3]
        
        genotypes, vt = indel.filter_heterozygous_heterogametes(self.genotypes,
                                                          self.vt,
                                                          heterogametic)

        expected_gt_output =\
        [[[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
         [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]],
         [[-1, -1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]]

        expected_vt_shape = (3,)

        npt.assert_array_equal(expected_gt_output,
                               genotypes)
        
        self.assertEqual(expected_vt_shape,
                         vt.shape)

    def test_filter_heterozygous_heterogametes_error_tol(self):

        '''I reuse the original mock genotypes and variant filters and
        "pretend" 2L is a sex chromosome'''
        
        heterogametic = [2, 3]
        
        genotypes, vt = indel.filter_heterozygous_heterogametes(self.genotypes,
                                                          self.vt,
                                                          heterogametic, 
                                                          error_tol=1)

        expected_gt_output = \
             [[[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
             [[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
             [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
             [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
             [[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]],
             [[-1, -1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]
             ]

        expected_vt_shape = (7,)

        npt.assert_array_equal(expected_gt_output,
                               genotypes)
        
        self.assertEqual(expected_vt_shape,
                         vt.shape)

    def test_autosomal_mendelian_filtering(self):

        violations = indel.id_mendelian_violations(self.genotypes[:-1])

        expected_violations = np.array([0,0,0,0,0,0,0,1])

        npt.assert_array_equal(violations, expected_violations)

    def test_sex_mendelian_filtering(self):
        
        violations = indel.id_mendelian_violations(self.genotypes[:-1], 
                                                   sex=True,
                                                   parents_homo_progeny=\
                                                   [0,1,3,4,5])

        expected_violations = np.array([0,0,0,0,0,0,0,0])

        npt.assert_array_equal(violations, expected_violations)

    def test_raises_error(self):
        
        self.assertRaises(ValueError, indel.id_mendelian_violations, 
                          self.genotypes[:-1], True)
        
    def test_remove_mendelian_violations_one(self):
    
        site_violations = np.array([0,0,1,0,1,4,0,0])

        genotypes, vt = indel.remove_mendelian_violations(self.genotypes[:-1], 
                                                          self.vt[:-1], 
                                                          site_violations)

        expected_gt_output = [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
                 [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
                 [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
                 [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
                 [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]]
                 ]

        npt.assert_array_equal(expected_gt_output, genotypes)

        self.assertEqual(vt.shape, (5,))
    
    def test_remove_mendelian_violations_two(self):
    
        site_violations = np.array([0,0,1,0,1,4,0,0])

        genotypes, vt = indel.remove_mendelian_violations(self.genotypes[:-1], 
                                                    self.vt[:-1], 
                                                    site_violations, 
                                                    permitted_violations=2)

        expected_gt_output = [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
                 [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
                 [[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
                 [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
                 [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
                 [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
                 [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]]
                 ]

        npt.assert_array_equal(expected_gt_output, genotypes)

        self.assertEqual(vt.shape, (7,))
        
    def test_filter_on_type_snp(self):

        records = [(b'2L', 10, 'A', 'T', True),
                   (b'2L', 20, 'GGG', 'G', False),
                   (b'2L', 30, 'C', 'T', True),
                   (b'2L', 49, 'CTGCTGC', 'C', False),
                   (b'2L', 110, 'G', 'A', True),
                   (b'2L', 300, 'TT', 'T', False),
                   (b'2L', 351, 'T', 'TTT', False)]

        dtype = [('CHROM', 'S4'),
                 ('POS', 'u4'),
                 ('REF', 'S100'),
                 ('ALT', 'S100'),
                 ('is_snp', '?')]

        test_vt = allel.VariantTable(records, dtype=dtype, 
                                                  index=('CHROM','POS'))

        test_genotypes = allel.GenotypeArray([[[0, 0], [1, 1], [0, 1], [0, 1]],
                                              [[0, 0], [1, 1], [0, 1], [1, 1]],
                                              [[1, 1], [0, 0], [1, 1], [1, 1]],
                                              [[0, 0], [1, 1], [0, 0], [0, 1]],
                                              [[0, 0], [1, 1], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]]
                                             ], dtype='i1')

        genotypes, vt = indel.filter_on_type(test_genotypes, test_vt, 
                                             variant_type="SNP")

        expected_genos = allel.GenotypeArray([[[0, 0], [1, 1], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [1, 1], [1, 1]],
                                              [[0, 0], [1, 1], [0, 1], [0, 1]]
                                              ], dtype='i1')

        npt.assert_array_equal(genotypes, expected_genos)

        self.assertEqual(vt.shape, (3,))
        
    def test_filter_on_type_indel(self):

        records = [(b'2L', 10, 'A', 'T', True),
                   (b'2L', 20, 'GGG', 'G', False),
                   (b'2L', 30, 'C', 'T', True),
                   (b'2L', 49, 'CTGCTGC', 'C', False),
                   (b'2L', 110, 'G', 'A', True),
                   (b'2L', 300, 'TT', 'T', False),
                   (b'2L', 351, 'T', 'TTT', False)]

        dtype = [('CHROM', 'S4'),
                 ('POS', 'u4'),
                 ('REF', 'S100'),
                 ('ALT', 'S100'),
                 ('is_snp', '?')]

        test_vt = allel.VariantTable(records, dtype=dtype, 
                                                  index=('CHROM','POS'))

        test_genotypes = allel.GenotypeArray([[[0, 0], [1, 1], [0, 1], [0, 1]],
                                              [[0, 0], [1, 1], [0, 1], [1, 1]],
                                              [[1, 1], [0, 0], [1, 1], [1, 1]],
                                              [[0, 0], [1, 1], [0, 0], [0, 1]],
                                              [[0, 0], [1, 1], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]]
                                             ], dtype='i1')

        genotypes, vt = indel.filter_on_type(test_genotypes, test_vt, 
                                             variant_type="indel")

        expected_genos = allel.GenotypeArray([[[0, 0], [1, 1], [0, 1], [1, 1]],
                                              [[0, 0], [1, 1], [0, 0], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]],
                                              [[1, 1], [0, 0], [0, 1], [0, 1]]
                                             ], dtype='i1')

        npt.assert_array_equal(genotypes, expected_genos)

        self.assertEqual(vt.shape, (4,))
        
    def test_filter_on_type_catches_bad_type(self):
        
        self.assertRaises(ValueError,
                          indel.filter_on_type, self.genotypes, self.vt, 
                          variant_type = "duck")

class IntegrationTest(unittest.TestCase):

    def setUp(self):

        ##mock genotype qualities
        self.gq = np.array([[27, 23, 19, 21, 15, 24],
                       [28, 21, 23, 20, 29, 25],
                       [32, 20, 29, 31, 30, 20],
                       [21, 27, 29, 32, 30, 33],
                       [22, 23, 29, 30, 27, 21],
                       [21, 24, 26, 27, 26, 25],
                       [20, 32, 29, 31, 30, 28],
                       [27, 26, 25, 24, 30, 31],
                       [25, 24, 32, 31, 28, 25]])

        ##one variant to get filtered out by GQ
        ##one variant to get filtered out as non-CDS
        ##one variant to get filtered out as Mendelian error
        ##one variant to get filtered out by MQ
        ##one variant to get filtered out by QD
        ##one variant to get filtered out as non-seg
        ##one variant to get filtered out with a missing parent
        ##two variants to rule them all, and in the darkness bind them

        ##mock genotypes
        self.genotypes =\
        allel.GenotypeArray(
            [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
             [[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
             ##the variant below is filtered out as non-segregating
             [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
             [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
             [[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             ##the variant below is filtered out as a Mendelian error
             [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]],
             ##the variant below is filtered out as missing parent
             [[-1, -1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]
             ], dtype='i1')

        ##mock variant table
        ##default min QD 15
        ##default min MQ 40
        ##variants filtered out:1 (GQ),
        ##2 (QD),
        ##3 (non-CDS),
        ##4 (non-seg),
        ##5 (MQ),
        ##8 (Mendelian error),
        ##9 (missing parent)
        records = [(b'2L', 20, b'A', b'ATC', 15, 35.7, 45, (6, 6)),
                   ##the variant below gets filtered out for QD
                   (b'2L', 53, b'CG', b'C', 18, 14, 45, (9, 3)),
                   ##the variant below is filtered out as non-CDS
                   (b'2L', 207, b'ATATA', b'A', 13, 29.2, 41, (3, 9)),
                   (b'2L', 602, b'TCTCT', b'T', 27, 16, 41, (0, 12)),
                   ##the variant below is filtered out for MQ
                   (b'2L', 754, b'TCTGT', b'TC', 18, 33.2, 29, (6, 6)),
                   (b'2L', 1090, b'G', b'GTTT', 17, 31.4, 43, (6, 6)),
                   (b'2L', 1114, b'GCAT', b'G', 20, 37.2, 44, (6, 6)),
                   (b'2L', 1129, b'T', b'TTTT', 32, 15, 40, (3, 9)),
                   (b'2L', 1150, b'AAA', b'A', 30, 20, 45, (0, 10))
                  ]

        dtype = [('CHROM', 'S4'),
                 ('POS', 'u4'),
                 ('REF', 'S100'),
                 ('ALT', 'S100'),
                 ('DP', int),
                 ('QD', float),
                 ('MQ', int),
                 ('AC', (int, 2))]

        self.vt = allel.VariantTable(records,
                                dtype=dtype, index=('CHROM', 'POS'))

        self.callset = {"variants": self.vt, 
                        "calldata/genotype": self.genotypes,
                        "calldata/GQ": self.gq}

        ##mock feature table
        feature_table_data = [
            ["2L", "DB", b'gene', 15, 2000, -1, "+", -1, "gene1", "."],
            ["2L", "DB", b'CDS', 17, 77, -1, "+", 0, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 569, 803, -1, "+", 2, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 1001, 1166, -1, "+", 2, "gene1-PA", 
             "gene1-RA"],
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

        feature_table_names = tuple(t[0] for t in feature_table_dtype)

        self.test_feature_table = \
        allel.FeatureTable(feature_table_data, names=feature_table_names)
        
        self.test_cds = \
        self.test_feature_table[self.test_feature_table.type == b'CDS']

        ##mock metadata
        sample_names = pd.Series(["A", "B", "C", "D", "E", "F"])
        sample_roles = \
        pd.Series(["parent", "parent", "progeny",
                   "progeny", "progeny", "progeny"])
        sample_sexes = pd.Series(["F", "M", "F", "M", "F", "M"])
        
        metadata = \
        pd.DataFrame({'samples': sample_names, 'role': sample_roles,\
                      'sex': sample_sexes})
        
        self.name = "2L"



'''class TestWorkflow(unittest.TestCase):

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
                                                      self.feature_table,
                                                      "CDS")

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
'''