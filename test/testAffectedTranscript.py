import allel
import indel
import numpy as np
import numpy.testing as npt
import unittest

class AffectedTranscriptTestCase(unittest.TestCase):

        def setUp(self):
            
            ##mock variant table
            
            records = [[b'3R',101,b'CG',[b'C','',''],2],
                       [b'3R',207,b'A',[b'AT',b'ATT',b'ATTT'],4],
                       [b'3R',220,b'T',[b'TC','',''],2],
                       [b'3R',284,b'A',[b'ACAG',b'ACAGCAG',''],3],
                       [b'3R',316,b'GAAA',[b'G','',''],2],
                       [b'3R',347,b'TTGTGTG',[b'T','',''],2],
                       [b'3R',1091,b'G',[b'GA',b'GAAGA',''],3],
                       [b'3R',74807,b'ATATA',[b'A',b'ATA',''],3],
                       [b'3R',16045837,b'CA',[b'C',b'CAA',b'CAAA'],4],
                       [b'3R',28475632,b'T',[b'TT','',''],2]]
            
            dtype = [('CHROM', 'S4'),
                ('POS', 'u4'),
                ('REF', 'S100'),
                ('ALT', ('S100',3),),
                ('num_alleles', 'uint8')]
            
            self.variant_table = allel.VariantTable(
                records, dtype=dtype, index=('CHROM','POS'))
            
            ##mock genotypes
            
            self.genotypes = allel.GenotypeArray([
              [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
              [[0,1],[2,2],[0,2],[0,2],[0,2],[0,2]],
              [[1,1],[0,0],[0,1],[0,1],[0,1],[0,1]],
              [[2,2],[0,1],[0,2],[1,2],[0,2],[1,2]],
              [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
              [[1,1],[0,0],[0,1],[0,1],[0,1],[0,1]],
                [[0,0],[1,1],[0,1],[0,1],[0,1],[0,1]],
                [[1,1],[1,1],[1,1],[1,1],[1,1],[1,1]],
                [[0,0],[1,2],[0,1],[0,2],[0,2],[0,2]],
                [[0,1],[0,1],[1,1],[0,0],[0,1],[0,0]]
            ], dtype='i1')
            
            ##mock feature table
            feature_table_data = [
                ['3R', 'DB', 'CDS', 95, 120, -1, '+', 0, 'gene3-PA', 'gene3-RA'],
                ['3R', 'DB', 'CDS', 206, 296, -1, '+', 1, 'gene3-PA', 'gene3-RA'],
                ['3R', 'DB', 'CDS', 310, 363, -1, '+', 0, 'gene3-PA', 'gene3-RA'],
                ['3R', 'DB', 'CDS', 15780, 15809, -1, '+', 0, 'gene24-PA', 'gene24-RA'],
                ['3R', 'DB', 'CDS', 21075362, 21075460, -1, '+', 0, 'gene436-PA', 'gene436-RA']
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
            allel.FeatureTable(feature_table_data, dtype=feature_table_dtype)
            
            ##mock IntersectingVariants
            
            #self.variants = 

            self.variants = [indel.IntersectingVariant('3R',101),
                            indel.IntersectingVariant('3R',207),
                            indel.IntersectingVariant('3R',220),
                            indel.IntersectingVariant('3R',284),
                            indel.IntersectingVariant('3R',316),
                            indel.IntersectingVariant('3R',347)]
            
            for variant in self.variants:
                
                variant.add_vtbl_record(self.variant_table)
                variant.add_genos(self.genotypes, self.variant_table)
            
            self.affected_transcript = indel.AffectedTranscript(b'3R','gene3-RA')
            
        def add_variants_to_transcript(self):
            
            for variant in self.variants:
                
                self.affected_transcript.add_variant(variant)
        
        def test_add_variant(self):
                
            self.add_variants_to_transcript()
                    
            ##test the vtbl_indices
                
            test_vtbl_indices = [0,1,2,3,4,5]
            npt.assert_array_equal(test_vtbl_indices, 
                                   self.affected_transcript.vtbl_indices)
                
            ##test the n_variants
                
            self.assertEqual(self.affected_transcript.n_variants, 6)
                
            ##test the variant_length_changes
                
            test_variant_length_changes = [
                {0,-1},{0,1,2,3},{0,1},{0,3,6},{0,-3},{0,-6}]
                
            npt.assert_array_equal(test_variant_length_changes, 
                                    self.affected_transcript.variant_length_changes)

        def test_extract_vtbl(self):
            
            self.add_variants_to_transcript()
                
            self.affected_transcript.extract_vtbl(self.variant_table)
                
            test_vtbl = self.variant_table[0:6]
                
            npt.assert_array_equal(test_vtbl, self.affected_transcript.vtbl)
                
        def test_calculate_all_possible_length_changes(self):
            
            self.add_variants_to_transcript()
                
            self.affected_transcript.extract_vtbl(self.variant_table)
            
            self.affected_transcript.calculate_all_possible_length_changes()
                
            test_possible_length_changes = {
                -10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10
            }
                
            npt.assert_array_equal(test_possible_length_changes, 
                                    self.affected_transcript.possible_length_changes)
                
                
        def test_extract_haplotypes(self):
                        
            ##TO-DO: what happens if the vtbl indices are out of order?
            ##TO-DO: what happens if the vtbl indices are negative?
            ##TO-DO: what happens if some of the genotypes can't be phased?
            
            self.affected_transcript.vtbl_indices = [0,3,5,6,7]
            
            mock_phased_genotype_array = {self.affected_transcript.chrom:
                                          allel.GenotypeArray([[[0,1],[0,0],[1,0],[0,0],[1,0]],
                                                   [[1,1],[0,0],[1,0],[1,0],[1,0]],
                                                   [[0,0],[0,1],[0,1],[0,0],[0,0]],
                                                   [[1,0],[0,0],[0,0],[1,0],[0,0]],
                                                   [[0,1],[0,0],[1,0],[0,0],[1,0]],
                                                   [[1,1],[0,0],[1,0],[1,0],[1,0]],
                                                   [[1,1],[0,0],[1,0],[1,0],[1,0]],
                                                   [[1,1],[0,0],[1,0],[1,0],[1,0]]], dtype='i1')}
            
            
            expected_genos_phased = np.array([[[0,1],[0,0],[1,0],[0,0],[1,0]],
                                             [[1,0],[0,0],[0,0],[1,0],[0,0]],
                                             [[1,1],[0,0],[1,0],[1,0],[1,0]],
                                             [[1,1],[0,0],[1,0],[1,0],[1,0]],
                                             [[1,1],[0,0],[1,0],[1,0],[1,0]]])
            
            expected_haplotypes = np.array([[0,1,0,0],
                                           [1,0,0,0],
                                           [1,1,0,0],
                                           [1,1,0,0],
                                           [1,1,0,0]])

            self.affected_transcript.extract_haplotypes(mock_phased_genotype_array)
            
            npt.assert_array_equal(self.affected_transcript.genos_phased, expected_genos_phased)
            
            npt.assert_array_equal(self.affected_transcript.haplotypes, expected_haplotypes)                
      

if __name__ =='__main__':
	unittest.main()
