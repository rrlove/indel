import allel
import indel
import numpy as np
import numpy.testing as npt
import unittest

class IntersectingVariantTestCase(unittest.TestCase):

        def setUp(self):
            
            ##mock genotypes
            self.genotypes = allel.GenotypeArray([[[0, 2], [1, 1], [0, 1], [1, 2]],
                                      [[0, 0], [0, 1], [0, 1], [0, 0]],
                                      [[1, 2], [0, 1], [0, 2], [1, 2]],
                                      [[0, 1], [2, 3], [0, 3], [0, 1]]],
                                     dtype='i1')
           
            ##mock variant table
            
            records = [[b'2L', 20, b'A', [b'ATC',b'ATCTC',b''], 3],
                    [b'2L', 53, b'CG', [b'C',b'',b''], 2],
                    [b'2L', 207, b'ATATA', [b'A',b'AT',b''], 3],
                    [b'2L', 1090, b'G', [b'GTTT',b'GT',b'GTT'], 4]
                      ]
            
            dtype = [('CHROM', 'S4'),
                ('POS', 'u4'),
                ('REF', 'S100'),
                ('ALT', ('S100',3),),
                ('num_alleles', 'uint8')]
            
            self.variant_table = allel.VariantTable(records, dtype=dtype, index=('CHROM','POS'))
            
            self.variant = indel.IntersectingVariant('2L',207)
            
        def test_initialization(self):
            
            expected_expression = "(POS == 207) & (CHROM == b'2L')"

            self.assertEqual(expected_expression, self.variant.expression)
            
        def test_add_vtbl_record(self):
            
            self.variant.add_vtbl_record(self.variant_table)
            
            self.assertEqual([3], self.variant.num_alleles)
            
            ##list is initialized with 0 to allow for haplotypes where this position has REF
            expected_length_changes = [0, -4, -3]
            
            self.assertEqual(expected_length_changes, self.variant.length_change)
            
            self.assertEqual(2, self.variant.vt_index)
            
        def test_add_genos(self):
            
            self.variant.add_genos(self.genotypes,self.variant_table)
            
            expected_genos = [[[1, 2], [0, 1], [0, 2], [1, 2]]]
            
            npt.assert_array_equal(expected_genos, self.variant.genos)

if __name__ =='__main__':
        unittest.main()

