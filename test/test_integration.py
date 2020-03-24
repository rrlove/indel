#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 18:17:21 2020

@author: beccalove
"""

import unittest

import numpy as np
import numpy.testing as npt
import pandas as pd

import allel
import indel

class IntegrationTestsAuto(unittest.TestCase):

    def setUp(self):

        ##feature table
        feature_table_data = [
            ["2L", "DB", b'gene', 15, 2000, -1, "+", -1, "gene1", "."],
            ["2L", "DB", b'CDS', 17, 77, -1,\
             "+", 0, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 569, 803, -1,\
             "+", 2, "gene1-PA", "gene1-RA"],
            ["2L", "DB", b'CDS', 1001, 1166, -1,\
             "+", 2, "gene1-PA", "gene1-RA"],
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
        
        self.test_cds =\
        self.test_feature_table[self.test_feature_table.type == b'CDS']
        
        self.name = "2L"
        
        self.path = "int_test_auto.h5"
        
        self.dtypes = [('CHROM', 'S4'),
                       ('POS', 'u4'),
                       ('REF', 'S100'),
                       ('ALT', 'S100'),
                       ('DP', int),
                       ('QD', float),
                       ('MQ', int),
                       ('AC', (int, 2)),
                       ('is_snp', '?')]
        
    def test_1(self):
        
        ##this tests the cross workflow with default quality filters, 
        ##retaining indels
        
        genotypes, vt = indel.cross_workflow(self.name, self.path, 
                                             self.test_cds)
        
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
            [[[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]]
            ], dtype='i1')
        
        ##create the expected output vtbl
        test_1_records = \
        [(b'2L', 1090, b'G', b'GTTT', 17, 31.4, 43, (6, 6), False),
         (b'2L', 1114, b'GCAT', b'G', 20, 37.2, 44, (6, 6), False)]

        expected_vt = allel.VariantTable(test_1_records,
                                       dtype=self.dtypes, 
                                       index=('CHROM', 'POS'))

        
        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)
        
    def test_2(self):
        
        ##this tests the cross workflow with default quality filters, 
        ##retaining SNPs
        
        genotypes, vt = indel.cross_workflow(self.name, self.path, 
                                             self.test_cds, filter_on="SNP")
        
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
            [[[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
             [[0, 0], [0, 1], [0, 0], [0, 1], [0, 0], [0, 1]]], dtype='i1')
        
        ##create the expected output vtbl
        test_2_records = [(b'2L', 65, b'T', b'G', 21, 21, 52, (3, 9), True),
                          (b'2L', 1160, b'C', b'T', 25, 22, 51, (9, 3), True)
                         ]

        expected_vt = allel.VariantTable(test_2_records,
                        dtype=self.dtypes, index=('CHROM', 'POS'))

        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)

    def test_3(self):
        
        ##this tests the cross workflow with modified quality filters, 
        ##retaining indels
        
        self.assertRaises(indel.ZeroVariantError, indel.cross_workflow, 
                          self.name, self.path, self.test_cds, min_gq=25, 
                          min_qd=25, allowed_missing=0)
        
    def test_4(self):
        
        ##this tests the wild-caught workflow, retaining indels
        
        genotypes, vt = indel.wild_workflow(self.name, self.path, 
                                            self.test_cds)
        
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
            [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             [[0, 0], [0, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
             [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
             [[0, 1], [0, 1], [0, 0], [0, 1], [1, 1], [0, 1]],
             [[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
             [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]],
             [[-1, -1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]],
            ], dtype='i1')

        
        ##create the expected output vtbl
        test_4_records = \
        [(b'2L', 20, b'A', b'ATC', 15, 35.7, 45, (6, 6), False),
           (b'2L', 53, b'CG', b'C', 18, 14, 45, (9, 3), False),
           (b'2L', 602, b'TCTCT', b'T', 27, 16, 41, (0, 12), False),
           (b'2L', 754, b'TCTGT', b'TC', 18, 33.2, 29, (6, 6), False),
           (b'2L', 1090, b'G', b'GTTT', 17, 31.4, 43, (6, 6), False),
           (b'2L', 1114, b'GCAT', b'G', 20, 37.2, 44, (6, 6), False),
           (b'2L', 1129, b'T', b'TTTT', 32, 15, 40, (3, 9), False),
           (b'2L', 1150, b'AAA', b'A', 30, 20, 45, (0, 10), False),
          ]

        expected_vt = allel.VariantTable(test_4_records,
                        dtype=self.dtypes, index=('CHROM', 'POS'))

        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)
    
    def test_5(self):
        
        ##this tests the wild-caught workflow, retaining SNPs
        
        genotypes, vt = indel.wild_workflow(self.name, self.path, 
                                            self.test_cds, filter_on="SNP")
        
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
            [[[1, 1], [0, 1], [0, 1], [1, 1], [0, 1], [1, 1]],
             [[0, 0], [0, 1], [0, 0], [0, 1], [0, 0], [0, 1]]
            ], dtype='i1')
        
        ##create the expected output vtbl
        test_5_records = [(b'2L', 65, b'T', b'G', 21, 21, 52, (3, 9), True),
               (b'2L', 1160, b'C', b'T', 25, 22, 51, (9, 3), True)
              ]

        expected_vt = allel.VariantTable(test_5_records,
                        dtype=self.dtypes, index=('CHROM', 'POS'))

        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)

class IntegrationTestsSex(unittest.TestCase):

    def setUp(self):

        ##feature table
        feature_table_data = [
                    ["X", "DB", b'gene', 15, 2000, -1, "+", -1, "gene1", "."],
                    ["X", "DB", b'CDS', 17, 47, -1,\
                     "+", 0, "gene3-PA", "gene3-RA"],
                    ["X", "DB", b'CDS', 569, 803, -1,\
                     "+", 2, "gene3-PA", "gene3-RA"],
                    ["X", "DB", b'CDS', 1001, 1166, -1,\
                     "+", 2, "gene3-PA", "gene3-RA"],
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
        
        self.name = "X"
        
        self.path = "./int_test_X.h5"
        
        self.dtypes = [('CHROM', 'S4'),
                       ('POS', 'u4'),
                       ('REF', 'S100'),
                       ('ALT', 'S100'),
                       ('DP', int),
                       ('QD', float),
                       ('MQ', int),
                       ('AC', (int, 2)),
                       ('is_snp', '?')]
        
        ##create the hets and the parents & homo progeny indices lists
        self.het_phasing = [1, 3, 5]
        
        self.parents_homo_progeny = [0, 1, 2, 4]
        
    def test_1(self):
        
        ##this tests the cross workflow with default quality filters, retaining indels
        
        genotypes, vt = indel.cross_workflow(self.name, self.path, 
                                             self.test_cds, sex=True,
                                             hets=self.het_phasing, 
                                             X_phasing=\
                                             self.parents_homo_progeny)
        
        ##create the expected output genotypes
        expected_genotypes =\
        allel.GenotypeArray(
            [[[1, 1], [0, 0], [0, 1], [1, 1], [0, 1], [1, 1]],
             [[0, 1], [1, 1], [0, 1], [0, 0], [1, 1], [0, 0]],
             [[0, 1], [1, 1], [0, 1], [1, 1], [1, 1], [0, 0]],
             ], dtype='i1')

        ##create the expected output vtbl
        test_1_records = \
        [(b'X', 582, b'ATATA', b'A', 13, 29.2, 41, (4, 8), False),
         (b'X', 754, b'TCTGT', b'TC', 18, 33.2, 51, (6, 6), False),
         (b'X', 1129, b'T', b'TTTT', 32, 15, 40, (5, 7), False)]

        expected_vt = allel.VariantTable(test_1_records,
                                dtype=self.dtypes, index=('CHROM', 'POS'))
        
        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)
        
    def test_2(self):
        
        ##this tests the cross workflow with default quality filters, retaining SNPs
        
        genotypes, vt = indel.cross_workflow(self.name, self.path, 
                                             self.test_cds, filter_on="SNP",
                                             sex=True, hets=self.het_phasing,
                                             X_phasing=\
                                             self.parents_homo_progeny)
                
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
            [[[0, 0], [1, 1], [0, 1], [0, 0], [0, 1], [0, 0]]], dtype='i1')
        
        ##create the expected output vtbl
        test_2_records = [(b'X', 1160, b'C', b'T', 25, 22, 51, (8, 4), True)]

        expected_vt = allel.VariantTable(test_2_records,
                        dtype=self.dtypes, index=('CHROM', 'POS'))

        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)

    def test_3(self):
        
        ##this tests the cross workflow with modified quality filters, retaining indels
        
        genotypes, vt = indel.cross_workflow(self.name, self.path, 
                                             self.test_cds,
                                             sex=True, hets=self.het_phasing,
                                             X_phasing=\
                                             self.parents_homo_progeny,
                                             min_gq=22, min_qd=10, min_mq=30, 
                                             allowed_missing=3)
                
        ##create the expected output genotypes
        expected_genotypes = allel.GenotypeArray(
        [[[0, 0], [1, 1], [0, 1], [0, 0], [0, 1], [0, 0]],
         [[1, 1], [0, 0], [0, 1], [1, 1], [0, 1], [1, 1]],
         [[0, 1], [1, 1], [0, 1], [0, 0], [1, 1], [0, 0]],
         [[0, 1], [0, 0], [0, 1], [1, 1], [0, 0], [0, 0]],
         [[0, 1], [1, 1], [0, 1], [1, 1], [1, 1], [0, 0]],
         ], dtype='i1')

        ##create the expected output vtbl
        test_3_records = \
        [(b'X', 20, b'A', b'ATC', 15, 35.7, 45, (8, 4), False),
         (b'X', 582, b'ATATA', b'A', 13, 29.2, 41, (4, 8), False),
         (b'X', 754, b'TCTGT', b'TC', 18, 33.2, 51, (6, 6), False),
         (b'X', 1090, b'G', b'GTTT', 17, 31.4, 38, (8, 4), False),
         (b'X', 1129, b'T', b'TTTT', 32, 15, 40, (5, 7), False)]

        expected_vt = allel.VariantTable(test_3_records,
                        dtype=self.dtypes, index=('CHROM', 'POS'))

        ##compare expected and actual outputs
        npt.assert_array_equal(genotypes, expected_genotypes)
        npt.assert_array_equal(vt, expected_vt)
