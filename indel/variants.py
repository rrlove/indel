#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 16:46:44 2018

@author: beccalove
"""

from Bio import Seq
#from Bio import SeqIO
from collections import namedtuple
import itertools
#import numpy as np
import requests

VariantRecord = namedtuple(
    'VariantRecord',['chrom','pos','ref','alt'])

class RequestReturnError(Exception):
    pass

class ExonSequence:
    
    def __init__(self, chrom, transcript, exon_number, strand, start, end):
        
        self.chrom = chrom
        self.transcript = transcript
        self.exon_number = exon_number
        self.strand = strand
        self.start = start
        self.end = end
        self.length = abs(self.end - self.start)
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
            
            if len(variant.alt) >= len(variant.ref):##variant is insertion or SNP
                
                if not variant.ref == self.sequence[index]:
                    return("Reference alleles don't match at" + variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
            elif len(variant.ref) > len(variant.alt):##variant is a deletion
                
                test_string = self.sequence[(index):(index + len(variant.ref))]
                
                if not variant.ref == test_string:
                    return("Reference alleles don't match at" + variant.pos)
                
                mutable_sequence_list[index] = variant.alt
                
                del mutable_sequence_list[(index + 1): (index + len(variant.ref))]
                
        self.changed_sequence = ''.join(mutable_sequence_list)
        
class Transcript:
    
    def __init__(self, name, seq):
        
        self.name = name
        self.seq = seq
        self.exons = []
        self.seq_translated = ''
        self.seq_changed = ''
        self.seq_changed_translated = ''
        
    def add_exon(self, exon):
        
        self.exons.append(exon)
        self.num_exons = len(self.exons)
        
    def translate_seqs(self, which = "both"):
        
        if which == "both" or which == "original":
            
            self.seq_translated = str(Seq.Seq(self.seq).translate())
            
        if which == "both" or which == "changed":
            
            self.seq_changed_translated = str(Seq.Seq(self.seq_changed).translate())
            
    def populate_exon_seq(self):

        for exon in self.exons:

            if exon.exon_number == 1:

                sequence = self.seq[:exon.length + 1]
                exon.add_sequence(sequence)

                start = exon.length + 1

            elif exon.exon_number > 1 and exon.exon_number < self.num_exons:

                end = start + exon.length + 1

                sequence = self.seq[start:end]
                exon.add_sequence(sequence)

                start = end

            elif exon.exon_number == self.num_exons:

                sequence = self.seq[start:]

                exon.add_sequence(sequence)
    
    def parse_variants_list(self, variants):
        
        for exon in self.exons:
    
        #print([exon.start <= variant.pos <= exon.end for variant in test_variants])
            variant_bool = \
            [exon.start <= variant.pos <= exon.end for variant in variants]
            
            exon.variants = itertools.compress(variants,variant_bool)
    
    def assemble_changed_seq(self):
        
        self.seq_changed = ''.join([exon.changed_sequence for exon in self.exons])
        
def make_POST_request(gene, 
                      post_server = "https://www.vectorbase.org/rest", 
                      post_ext = "/sequence/id/",
                     feature = "cds"):

    headers = {'Content-type' : 'application/json', 'Accept' : 'application/json'}
    
    post_string = post_server + post_ext
    
    feature_seq_payload = '{"ids" : ["' + gene + '"], "type" : "' + feature + '"}'
    
    feature_seq = requests.post(post_string, data = feature_seq_payload, headers = headers)
    
    if not feature_seq.status_code == requests.codes.ok:
        
        raise RequestReturnError("POST request status: ", feature_seq.status_code)
    
    #transcript = Transcript(name = feature_seq.json()[0]["id"], 
     #                       seq = feature_seq.json()[0]["seq"])
    
    return feature_seq.json()

def make_transcript(feature_seq_json):
    
    transcript = Transcript(name = feature_seq_json[0]["id"], 
                            seq = feature_seq_json[0]["seq"])
    
    return transcript

def make_GET_request(gene,
                    get_server = "https://www.vectorbase.org/rest",
                    feature_types = ["transcript","exon","cds"]):
    
    headers = {'Content-type' : 'application/json', 'Accept' : 'application/json'}
    
    get_ext = "/overlap/id/" + gene + "?"
    
    get_string = get_server + get_ext + "".join(\
            ["feature=" + feature + ";" for feature in feature_types]).rstrip(";")
    
    feature_coords = requests.get(get_string, headers = headers)
    
    if not feature_coords.status_code == requests.codes.ok:
        
        raise RequestReturnError("GET request status: ", feature_coords.status_code)
        
    return feature_coords.json()

def make_exons(feature_coords_json, gene):
    
    cds_list = []
    exon_list = []
    exon_object_list = []
    
    for item in feature_coords_json:
    
        if item["feature_type"] == "cds":

            if item["Parent"] == gene:

                cds_list.append(item)

        elif item["feature_type"] == "exon":

            if item["Parent"] == gene:

                exon_list.append(item)
            
    if not len(cds_list) == len(exon_list):
        
        return("Must be same # of exons and CDSes")
        
    for i in range(len(exon_list)):
    
        exon = ExonSequence(chrom = exon_list[i]["seq_region_name"],
                        transcript = gene,
                       exon_number = exon_list[i]["rank"],
                       strand = exon_list[i]["strand"],
                       start = cds_list[i]["start"],
                       end = cds_list[i]["end"])
        
        exon_object_list.append(exon)
        
    return exon_object_list

def run_workflow(gene_name, variants):
    
    transcript = make_transcript(make_POST_request(gene = gene_name))
    
    for exon in make_exons(make_GET_request(gene = gene_name), gene_name):
        transcript.add_exon(exon)
        
    transcript.populate_exon_seq()
    
    transcript.parse_variants_list(variants)
    
    for exon in transcript.exons:
    
        exon.change(exon.variants)
        
    transcript.assemble_changed_seq()
    
    transcript.translate_seqs()
    
    return transcript
        
        