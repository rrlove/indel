#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:37:17 2018

@author: beccalove
"""

from collections import namedtuple
from . import Transcript
from . import ExonSequence
import requests

VariantRecord = namedtuple(
    'VariantRecord',['chrom','pos','ref','alt'])

class RequestReturnError(Exception):
    pass

def make_POST_request(gene, 
                      post_server = "https://www.vectorbase.org/rest", 
                      post_ext = "/sequence/id/",
                     feature = "cds"):

    headers = {'Content-type' : 'application/json',
               'Accept' : 'application/json'}
    
    post_string = post_server + post_ext
    
    feature_seq_payload = '{"ids" : ["' + gene + '"],\
    "type" : "' + feature + '"}'
    
    feature_seq = requests.post(post_string, 
                                data = feature_seq_payload, headers = headers)
    
    if not feature_seq.status_code == requests.codes.ok:
        
        raise RequestReturnError("POST request status: ", 
                                 feature_seq.status_code)
    
    return feature_seq.json()

def make_transcript(feature_seq_json):
    
    transcript = Transcript(name = feature_seq_json[0]["id"], 
                            seq = feature_seq_json[0]["seq"])
    
    return transcript

def make_GET_request(gene,
                    get_server = "https://www.vectorbase.org/rest",
                    feature_types = ["transcript","exon","cds"]):
    
    headers = {'Content-type' : 'application/json', \
               'Accept' : 'application/json'}
    
    get_ext = "/overlap/id/" + gene + "?"
    
    get_string =\
    get_server + get_ext +\
    "".join(["feature=" + feature +\
             ";" for feature in feature_types]).rstrip(";")
    
    feature_coords = requests.get(get_string, headers = headers)
    
    if not feature_coords.status_code == requests.codes.ok:
        
        raise RequestReturnError("GET request status: ",
                                 feature_coords.status_code)
        
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