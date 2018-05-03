#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:46:33 2018

@author: beccalove
"""

import indel
import indel.variants
import unittest

class TranscriptTestCase(unittest.TestCase):
    
    def setUp(self):
        
        test_seq = "ATGCTTCATCAGCAGCAGCCCGGGGCCTAG"
        self.transcript = indel.variants.Transcript("2L", test_seq)
        
    def test_translate_seq(self):
        
        ##test the translated versions
        
        self.transcript.seq_changed = "ATGCGTCATCAGCAGCAGCAGTCCGGGGCCTAG"
        
        unchanged_translation = "MLHQQQPGA*"
        
        changed_translation = "MRHQQQQSGA*"
        
        self.transcript.translate_seqs()
        
        self.assertEqual(self.transcript.seq_translated, 
                         unchanged_translation)
        
        self.assertEqual(self.transcript.seq_changed_translated,
                         changed_translation)
        
    def test_populate_exon_seq_length_check(self):
        ##pass a set of exons whose collective length is less than the length
        ##of the sequence, and make sure an AssertionError is raised

        exon1 = indel.ExonSequence("2L", "foo", 1, "+", 5, 21)
        exon2 = indel.ExonSequence("2L", "foo", 2, "+", 77, 90)
        
        self.transcript.add_exon(exon1)
        self.transcript.add_exon(exon2)
        
        self.assertRaises(AssertionError, self.transcript.populate_exon_seq)
        
    def test_populate_exon_seq_correct_seqs(self):
        
        exon1 = indel.ExonSequence("2L", "foo", 1, "+", 5, 21)#17
        exon2 = indel.ExonSequence("2L", "foo", 2, "+", 50, 55)#6
        exon3 = indel.ExonSequence("2L", "foo", 3, "+", 70, 74)#5
        exon4 = indel.ExonSequence("2L", "foo", 4, "+", 103, 107)#5
        
        for exon in [exon1, exon2, exon3, exon4]:
            self.transcript.add_exon(exon)
        
        self.transcript.populate_exon_seq()
        
        self.assertEqual(self.transcript.exons[0].sequence,
                         "ATGCGTCATCAGCAGCA")
        
        self.assertEqual(self.transcript.exons[1].sequence,
                         "GCAGTC")
        
        self.assertEqual(self.transcript.exons[1].sequence,
                         "CGGGG")
                
        self.assertEqual(self.transcript.exons[1].sequence,
                         "CCTAG")
