#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:02 2018

@author: beccalove
"""

class ExonSequence:
    
    def __init__(self, chrom, transcript, exon_number, strand, start, end):
        
        self.chrom = chrom
        self.transcript = transcript
        self.exon_number = exon_number
        self.strand = strand
        self.start = start
        self.end = end
        self.length = abs(self.end - self.start) + 1
        self.variants = []
        
    def add_sequence(self, sequence):
        
        self.sequence = sequence
        
    def change(self, variants):
        
        mutable_sequence_list = list(self.sequence)
        
        for variant in variants:
            
            if not variant.chrom == self.chrom:
                return("Variant chrom and exon chrom don't match")
            
            if not self.start <= variant.pos <= self.end:
                return("Variant is out of bounds of exon")
            
            index = variant.pos - self.start
            
            if len(variant.alt) >= len(variant.ref):
                ##variant is insertion or SNP
                
                if not variant.ref == self.sequence[index]:
                    return("Reference alleles don't match at" + variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
            elif len(variant.ref) > len(variant.alt):##variant is a deletion
                
                test_string = self.sequence[(index):(index + len(variant.ref))]
                
                if not variant.ref == test_string:
                    return("Reference alleles don't match at" + variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
                del \
                mutable_sequence_list[(index + 1):(index + len(variant.ref))]
                
        self.changed_sequence = ''.join(mutable_sequence_list)
