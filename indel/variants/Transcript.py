#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:36:43 2018

@author: beccalove
"""

from Bio import Seq
import itertools

class Transcript:

    def __init__(self, name, seq):

        self.name = name
        self.seq = seq
        self.exons = []
        self.num_exons = 0
        self.seq_translated = ''
        self.seq_changed = ''
        self.seq_changed_translated = ''
        self.fluffy_bunnies = True

    def add_exon(self, exon):

        self.exons.append(exon)
        self.num_exons += 1

    def translate_seqs(self, which="both"):

        if which == "both" or which == "original":

            self.seq_translated = str(Seq.Seq(self.seq).translate())

        if which == "both" or which == "changed":

            self.seq_changed_translated =\
            str(Seq.Seq(self.seq_changed).translate())

    def populate_exon_seq(self):

        ##first, check to make sure that the length of the transcript
        ##is the same as the collective length of the exons

        assert len(self.exons) > 0, "Must add exons to transcript first"

        assert len(self.seq) == \
        sum([exon.length for exon in self.exons]), \
        "Length of sequence does not match total length of exons"

        for exon in self.exons:

            if exon.exon_number == 1:

                sequence = self.seq[:exon.length]
                exon.add_sequence(sequence)

                start = exon.length

            elif exon.exon_number > 1 and exon.exon_number < self.num_exons:

                end = start + exon.length

                sequence = self.seq[start:end]
                exon.add_sequence(sequence)

                start = end

            elif exon.exon_number == self.num_exons:

                sequence = self.seq[start:]

                exon.add_sequence(sequence)

    def parse_variants_list(self, variants):

        for exon in self.exons:

            variant_bool = \
            [exon.start <= variant.pos <= exon.end for variant in variants]

            exon.variants = itertools.compress(variants, variant_bool)

    def assemble_changed_seq(self):

        self.seq_changed =\
        ''.join([exon.changed_sequence for exon in self.exons])
