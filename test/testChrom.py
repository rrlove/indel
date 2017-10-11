##TODO 10-11-17 track down unexpected behavior in importing multi-chrom files-- assert in Chrom and test case in ChromTestCase may not be needed

import allel
import indel
import numpy.testing as npt
import unittest


class ChromTestCase(unittest.TestCase):

        def setUp(self):
            self.chrom = Chrom("3L","/Users/beccalove 1/Desktop/Indels/code/indel/test/data/test.h5")
            self.test_feature_table = allel.FeatureTable.from_gff3("/Users/beccalove 1/Desktop/Indels/code/indel/test/data/test.gff3", attributes=['Parent'])
            ##self.badchrom = Chrom("3L","/Users/beccalove 1/Desktop/Indels/code/indel/test/data/badtest.h5")
                
        def testReadCalls(self):
            self.chrom.read_calls()
            self.assertEqual(len(self.chrom.callset), 1)
                
        def test_ID_positions_right_type(self):
            self.chrom.ID_positions(self.test_feature_table)
            self.assertEqual(type(self.chrom.allPositions), allel.model.ndarray.SortedIndex)
           
        def test_ID_positions_right_content(self):
            self.chrom.ID_positions(self.test_feature_table)
            expected_features = [True,False,False,False,True,False,True,True,False,False]
            expected_positions = [False,True,True,True,True,True,False,False,False,False,False,False,False,False,False]
            
            npt.assert_equal(expected_features, self.chrom.overlappedFeatures)
            npt.assert_equal(expected_positions, self.chrom.exonicPositions)
            
        '''def test_multichrom_file_breaks(self):
            self.badchrom.read_calls()
            self.assertRaises(ValueError, self.badchrom.ID_positions(self.test_feature_table))
            self.assertRaises(SyntaxError, self.badchrom.ID_positions(self.test_feature_table))'''
