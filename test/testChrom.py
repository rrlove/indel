##TODO 10-11-17 track down unexpected behavior in importing multi-chrom files-- assert in Chrom and test case in ChromTestCase may not be needed

import allel
import indel
import numpy as np
import numpy.testing as npt
import pandas as pd
import unittest


class ChromTestCase(unittest.TestCase):

        def setUp(self):
            
            ##mock genotype qualities
            gq = np.array([[27, 23, 19, 21, 15, 24],
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
            gt =\
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
                                     ],dtype='i1')
           
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
            records = [[b'2L', 20, b'A',b'ATC', 15, 35.7, 45, (6, 6)],
                        ##the variant below gets filtered out for QD
                    [b'2L', 53, b'CG', b'C',18, 14, 45, (9, 3)],
                    ##the variant below is filtered out as non-CDS
                    [b'2L', 207, b'ATATA', b'A', 13, 29.2, 41, (3, 9)],
                    [b'2L', 602, b'TCTCT', b'T', 27, 16, 41, (0, 12)],
                    ##the variant below is filtered out for MQ
                    [b'2L', 754, b'TCTGT', b'TC', 18, 33.2, 29, (6, 6)],
                    [b'2L', 1090, b'G', b'GTTT', 17, 31.4, 43, (6, 6)],
                    [b'2L', 1114, b'GCAT', b'G', 20, 37.2, 44, (6, 6)],
                    [b'2L', 1129, b'T', b'TTTT', 32, 15, 40, (3, 9)],
                    [b'2L', 1150, b'AAA', b'A', 30, 20, 45, (0, 10)]
                      ]
            dtype = [('CHROM', 'S4'),
                ('POS', 'u4'),
                ('REF', 'S100'),
                ('ALT', 'S100'),
                ('DP', int),
                ('QD', float),
                ('MQ', int),
                ('AC', (int, 2))]
            
            vt = allel.VariantTable(records, 
                                    dtype=dtype, index=('CHROM','POS'))
            
            callset = {"variants": vt, "calldata/genotype": gt,\
                       "calldata/GQ": gq}
            
            self.chrom = indel.Chrom("2L",'fake_path')
            self.chrom.callset = callset
            self.chrom.genotypes = gt
            self.chrom.vt = vt
            self.chrom.gq = gq
            
            ##mock feature table
            feature_table_data = [
                ["2L", "DB", b'gene', 15, 2000, -1, "+", -1, "gene1" , "."],
                ["2L", "DB", b'CDS', 17, 77, -1,\
                 "+", 0, "gene1-PA", "gene1-RA"],
                ["2L", "DB", b'CDS', 569, 803, -1,\
                 "+", 2, "gene1-PA", "gene1-RA"],
                ["2L", "DB", b'CDS', 1001, 1166, -1,\
                 "+", 2, "gene1-PA", "gene1-RA"],
                ]
            
            feature_table_dtype = [
                ('seqid','S4'),
                ('source','S2'),
                ('type','S15'),
                ('start',int),
                ('end',int),
                ('score','<f8'),
                ('strand','S1'),
                ('phase','<i8'),
                ('ID','S15'),
                ('Parent','S15'),
                ]
            
            feature_table_names = tuple(t[0] for t in feature_table_dtype)
            
            
            self.test_feature_table = \
            allel.FeatureTable(feature_table_data, names=feature_table_names)
            ##self.chrom.read_calls()
            
            ##mock metadata
            sample_names = pd.Series(["A","B","C","D","E","F"])
            sample_roles = \
            pd.Series(["parent", "parent", "progeny",
                       "progeny", "progeny", "progeny"])
            sample_sexes = pd.Series(["F", "M", "F", "M", "F", "M"])
            self.chrom.metadata = \
            pd.DataFrame({'samples': sample_names, 'role': sample_roles,\
                          'sex': sample_sexes})
            
            self.chrom.ID_positions(
                    self.test_feature_table[
                            self.test_feature_table["type"] == b'CDS'])
            
        def test_check_metadata_parental(self):
            
            ##make some ill-formatted metadata
            ##replace the real metadata with it
            bad_names = pd.Series(["A","B","C","D"])
            bad_roles_1 = pd.Series(["parent","progeny","progeny","progeny"])
            bad_roles_2 = pd.Series(["progeny","parent","progeny","progeny"])
            bad_roles_3 = pd.Series(["parent","parent","parent","parent"])
            
            self.md_test_chrom = indel.Chrom("3L","fake_path")
            
            self.md_test_chrom.metadata = \
            pd.DataFrame({'samples': bad_names, 'role': bad_roles_1})
            self.assertRaises(AssertionError,
                              self.md_test_chrom.check_metadata_parental)
            
            self.md_test_chrom.metadata = \
            pd.DataFrame({'samples': bad_names, 'role': bad_roles_2})
            self.assertRaises(AssertionError,
                              self.md_test_chrom.check_metadata_parental)
            
            self.md_test_chrom.metadata = \
            pd.DataFrame({'samples': bad_names, 'role': bad_roles_3})
            self.assertRaises(AssertionError,
                              self.md_test_chrom.check_metadata_parental)
        
        def test_ID_positions_right_type(self):
            
            self.assertEqual(type(self.chrom.allPositions),
                             allel.model.ndarray.SortedIndex)
           
        def test_ID_positions_right_content(self):

            expected_features = [True,True,True]
            expected_positions = [True,True,False,True,True,
                                  True,True,True,True]
            
            npt.assert_equal(expected_features, self.chrom.overlappedFeatures)
            npt.assert_equal(expected_positions, self.chrom.exonicPositions)
            
        def test_extract_exonic(self):

            self.chrom.extract_exonic()
            
            ##these differ from the original test data 
            ##because one variant doesn't overlap a CDS
            expected_exonic_gt_shape = (8, 6, 2)
            expected_exonic_gq_shape = (8, 6)
            expected_exonic_vt_shape = (8,)
            
            self.assertEqual(expected_exonic_gt_shape, 
                             self.chrom.genotypes.shape)
            self.assertEqual(expected_exonic_gq_shape, self.chrom.gq.shape)
            self.assertEqual(expected_exonic_vt_shape, self.chrom.vt.shape)
            
        def test_quality_filtering(self):
            
            self.chrom.extract_exonic()
            self.chrom.filter_GQ_MQ_QD(allowed_missing=1)
            
            self.assertEqual(self.chrom.num_present, 6)
            self.assertEqual(self.chrom.genotypes_GQ_filtered.shape, (7, 6, 2))
            self.assertEqual(self.chrom.vtbl_GQ_filtered.shape, (7,))
            
            self.assertEqual(self.chrom.gt_qual_filtered.shape, (5, 6, 2))
            self.assertEqual(self.chrom.vt_qual_filtered.shape, (5,))
            
        def test_filter_on_parents(self):
            
            self.chrom.extract_exonic()
            self.chrom.filter_GQ_MQ_QD()
            self.chrom.filter_on_parents()
            
            self.assertEqual(self.chrom.gt_filtered.shape, (3, 6, 2))
            self.assertEqual(self.chrom.vt_filtered.shape, (3, ))
            
            expected_parental_filter_genos = \
            [[[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
            [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
            [[0, 1], [1, 1], [0, 0], [1, 1], [1, 1], [1, 1]]]
            
            npt.assert_array_equal(expected_parental_filter_genos, 
                                   self.chrom.gt_filtered)
            
        def test_autosomal_Mendelian_filtering(self):
            
            self.chrom.extract_exonic()
            self.chrom.filter_GQ_MQ_QD()
            self.chrom.filter_on_parents()
            self.chrom.ID_autosomal_Mendelian_violations()
            
            expected_autosomal_violations = [0, 0, 1]
            npt.assert_equal(expected_autosomal_violations, 
                             self.chrom.auto_site_violations)
            
        def test_filter_heterozygous_heterogametes(self):
            
            '''I reuse the original mock genotypes and variant filters and
            "pretend" 2L is a sex chromosome'''
            self.chrom.gt_filtered = self.chrom.callset["calldata/genotype"]
            self.chrom.vt_filtered = self.chrom.callset["variants"]
            
            heterogametic = [2,3]
            self.chrom.filter_heterozygous_heterogametes(heterogametic)
            
            expected_gt_output = [[[1, 1], [1, 1], [1, 1],
                                   [1, 1], [1, 1], [1, 1]],
                                  [[0, 1], [1, 1], [0, 0], 
                                   [1, 1], [1, 1], [1, 1]],
                                  [[-1, -1], [1, 1], [1, 1], 
                                   [1, 1], [1, 1], [1, 1]]]
            
            expected_vt_shape = (3,)
            
            self.assertEqual(expected_vt_shape, 
                             self.chrom.vt_het_error_filtered.shape)
            npt.assert_array_equal(expected_gt_output, 
                                   self.chrom.gt_het_error_filtered)
         
        def test_sex_Mendelian_filtering(self):
            
            self.chrom.gt_het_error_filtered = \
            allel.GenotypeArray([
                                [[0, 0], [1, 1], [0, 1],
                                 [0, 1], [0, 1], [0, 1]],
                                [[0, 0], [1, 1], [0, 1],
                                 [1, 1], [0, 1], [1, 1]],
                                [[1, 1], [0, 0], [1, 1],
                                 [1, 1], [1, 1], [1, 1]],
                                [[0, 0], [1, 1], [0, 0],
                                 [0, 1], [1, 1], [0, 1]],
                                     ],dtype='i1')
            
            parents_homo_progeny = [0,1,4,5]
        
            expected_gt_violations = [[0,0],[0,1],[1,1],[1,0]]
            
            self.chrom.ID_Mendelian_violations_sex(parents_homo_progeny)
            
            npt.assert_array_equal(expected_gt_violations, 
                                   self.chrom.sex_gt_violations)
            
            expected_site_violations = np.array([0,1,2,1])
            
            npt.assert_array_equal(expected_site_violations, 
                                   self.chrom.sex_site_violations)
            
        def test_remove_Mendelian_violations_auto(self):
            
            self.chrom.extract_exonic()
            self.chrom.filter_GQ_MQ_QD()
            self.chrom.filter_on_parents()
            self.chrom.ID_autosomal_Mendelian_violations()
            self.chrom.remove_Mendelian_violations()
            
            expected_Mendel_error_filtered_genos = \
            [[[0, 1], [0, 1], [0, 1], [1, 1], [0, 1], [0, 0]],
             [[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]]]
            
            npt.assert_array_equal(expected_Mendel_error_filtered_genos,
             self.chrom.gt_Mendel_filtered)
            
        def test_remove_Mendelian_violations_sex(self):
            
            self.chrom.gt_het_error_filtered = \
            allel.GenotypeArray([
                                [[0, 0], [1, 1], [0, 1],
                                 [0, 1], [0, 1], [0, 1]],
                                [[0, 0], [1, 1], [0, 1],
                                 [1, 1], [0, 1], [1, 1]],
                                [[1, 1], [0, 0], [1, 1], 
                                 [1, 1], [1, 1], [1, 1]],
                                [[0, 0], [1, 1], [0, 0], 
                                 [0, 1], [1, 1], [0, 1]],
                                     ],dtype='i1')
            
            records = [[b'2L', 20, 'A', 'ATC'],
                      [b'2L', 53, 'CG', 'C'],
                      [b'2L', 207, 'ATATA', 'A'],
                      [b'2L', 602, 'TCTCT', 'T']]
        
            dtype = [('CHROM', 'S4'),
                ('POS', 'u4'),
                ('REF', 'S100'),
                ('ALT', 'S100')]
            
            self.chrom.vt_het_error_filtered = \
            allel.VariantTable(records, dtype=dtype, index=('CHROM','POS'))
            
            self.chrom.sex_site_violations = np.array([0,1,2,1])
            
            self.chrom.remove_Mendelian_violations(sex=True)
            
            expected_Mendel_error_filtered_genos = \
            [[[0, 0], [1, 1], [0, 1], [0, 1], [0, 1], [0, 1]]]
            
            npt.assert_array_equal\
            (expected_Mendel_error_filtered_genos, 
             self.chrom.gt_Mendel_filtered)
    
        def test_phase_catches_misshapen_genos(self):
        
            self.phase_test_chrom = indel.Chrom("3L","fake_path")
            
            sample_names = pd.Series(["A","B","C","D"])
            sample_roles = pd.Series(["parent", "parent", 
                                      "progeny", "progeny"])
            sample_sexes = pd.Series(["F", "M", "F", "M"])
            self.phase_test_chrom.metadata = \
            pd.DataFrame({'samples': sample_names, 'role': sample_roles,\
                          'sex': sample_sexes})
                    
            self.phase_test_chrom.num_present = 11
            self.phase_test_chrom.gt_Mendel_filtered = \
            allel.GenotypeArray([[[0,0],[0,1],[0,1],[0,0]],
                                 [[1,1],[0,1],[0,1],[1,1]]], dtype='i1')
        
            self.assertRaises(AssertionError, 
                              self.phase_test_chrom.phase_and_filter)
        
        @unittest.expectedFailure
        def test_phase_output(self):
            
            self.chrom.gt_Mendel_filtered = \
            allel.GenotypeArray([[[0,0],[1,1],[0,1],[0,1]],
                                 [[0,1],[0,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,1],[0,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,1],[0,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,1],[0,1],[0,0],[1,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,0],[1,1],[0,1],[0,1]],
                                 [[0,1],[0,1],[0,0],[1,1]]], dtype='i1')
            
            
            records = [[b'2L', 1, 'A', 'AT'],
                       [b'2L', 2, 'C', 'CG'],
                       [b'2L', 3, 'TC', 'T'],
                       [b'2L', 4, 'GG', 'G'],
                       [b'2L', 5, 'T', 'TA'],
                       [b'2L', 6, 'G', 'GC'],
                       [b'2L', 7, 'ATA', 'A'],
                       [b'2L', 8, 'C', 'CC'],
                       [b'2L', 9, 'CCC', 'C'],
                       [b'2L', 10, 'G', 'GA'],
                       [b'2L', 11, 'GT', 'G'],
                       [b'2L', 12, 'T', 'TTTT'],
                       [b'2L', 13, 'ATAA', 'A'],
                       [b'2L', 14, 'C', 'CT'],
                       [b'2L', 15, 'CCC', 'C']]
        
            dtype = [('CHROM', 'S4'),
                ('POS', 'u4'),
                ('REF', 'S100'),
                ('ALT', 'S100')]
            
            self.chrom.vt_Mendel_filtered = \
            allel.VariantTable(records, dtype=dtype, index=('CHROM','POS'))
            
            self.chrom.num_present = 4
            
            self.chrom.phase_and_filter(window=10)
            
            expected_phased_Bool = \
            [True,False,True,False,True,False,True,True,False,
             True,True,True,True,True,False]
            
            expected_genos_shape = [10,4]
            expected_vt_shape = [10,]
            
            npt.assert_array_equal(expected_phased_Bool,self.chrom.phased_Bool)
            self.assertEqual(expected_genos_shape,self.chrom.phased_genos)
            self.assertEqual(expected_vt_shape,self.chrom.vt_phased)        
        
        @unittest.expectedFailure
        def test_setup(self):
            
            print(self.md_test_chrom)
            
        @unittest.expectedFailure            
        def test_use_of_npt(self):
            
            test_array_1 = [[5,6],[6,5]]
            test_array_2 = [[5,6],[6,6]]
            npt.assert_array_equal(test_array_1,test_array_2)
                                    
if __name__ =='__main__':
	unittest.main()
