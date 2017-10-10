import allel
import unittest
import indel

class ChromTestCase(unittest.TestCase):

        def setUp(self):
            self.chrom = Chrom("3L","/Users/beccalove 1/Desktop/Indels/code/indel/test/test.h5")
            self.test_feature_table = allel.FeatureTable.from_gff3("/Users/beccalove 1/Desktop/Indels/code/indel/test/test.gff3", attributes=['Parent'])
                
        def testReadCalls(self):
            self.chrom.read_calls()
            self.assertEqual(len(self.chrom.callset), 1)
                
        def test_ID_positions(self):
            
            self.chrom.ID_positions(self.test_feature_table)
            self.assertEqual(type(self.chrom.allPositions), allel.model.ndarray.SortedIndex)
