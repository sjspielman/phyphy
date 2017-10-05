#!/usr/bin/env python

##############################################################################
##  phyhy: Python HyPhy: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################
"""
    Unit tests for the phyphy_parser.py module.
"""

import unittest
import os
from phyphy import *
ZERO=1e-8
DECIMAL=8


"""

 _unpack_json(self):
 _determine_analysis_from_json(self):
 _count_partitions(self):
 _extract_slac_sitetable(self, raw, slac_by, slac_ancestral_type):
 _parse_sitemethod_to_csv(self, delim, slac_by = "by-site", slac_ancestral_type = "AVERAGED"):
 _parse_absrel_to_csv(self, delim):
 _parse_relrates_to_csv(self, delim):
 _reform_rate_phrase(self, phrase):
 _replace_tree_info(self, tree, node, newvalue):
 extract_model_names(self):
 extract_model_component(self, model_name, component):
 extract_model_logl(self, model_name):
 extract_model_estimated_parameters(self, model_name):
 extract_model_aicc(self, model_name):
 extract_model_rate_distributions(self, model_name):
 extract_model_frequencies(self, model_name, as_dict = False):
 extract_branch_sets(self, by_set = False):
 extract_input_tree(self, partition = None):
 reveal_branch_attributes(self):
 extract_branch_attribute(self, attribute_name, partition = None, map = False):
 map_branch_attribute(self, attribute_name, partition = None):
 extract_model_tree(self, model, partition = None):
 reveal_fields(self):
 extract_csv(self, csv, delim = ",", slac_by = "by-site", slac_ancestral_type = "AVERAGED"):
 extract_timers(self):
 extract_site_logl(self):
 extract_evidence_ratios(self):

"""

### copied from a pyvolve unittest so i can format knowledgably once i feel like writing tests.
# 
# class phyphy_parser_(unittest.TestCase):
# 
#     
#     def setUp(self):
#         ''' 
#             Tree and frequency set-up.
#         '''
#         dummy_json = 
#         self.myparser = HyPhyParser() 
#         
#     def test_evolver_mrca_nucleotide(self):
#         rootseq = "AAATTTCCCGGG"
#         m = Model("nucleotide")
#         p = Partition(root_sequence = rootseq, models = m)
#         evolve = Evolver(partitions = p, tree = self.tree)
#         evolve(ratefile = False, infofile=False, seqfile=False)
#         seqdict = evolve.get_sequences(anc = True)
#         self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for nucleotide evolution.")
# 
# 
# 
#     def test_evolver_mrca_aminoacid(self):
#         rootseq = "ARGYMMKLPQ"
#         m = Model("WAG")
#         p = Partition(root_sequence = rootseq, models = m)
#         evolve = Evolver(partitions = p, tree = self.tree)
#         evolve(ratefile = False, infofile=False, seqfile=False)
#         seqdict = evolve.get_sequences(anc = True)
#         self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for amino acid evolution.")
# 
# 
#     def test_evolver_mrca_codon(self):
#         rootseq = "AACCGATTTGGCCAT"
#         m = Model("codon", {"beta":0.5})
#         p = Partition(root_sequence = rootseq, models = m)
#         evolve = Evolver(partitions = p, tree = self.tree)
#         evolve(ratefile = False, infofile=False, seqfile=False)
#         seqdict = evolve.get_sequences(anc = True)
#         self.assertTrue(seqdict["root"] == rootseq, msg = "MRCA not preserved for codon evolution.")
# 

            
# def run_evolver_test():
# 
#     run_tests = unittest.TextTestRunner()
# 
#     print "Testing evolver no het, one partition"
#     test_suite0 = unittest.TestLoader().loadTestsFromTestCase(evolver_singlepart_nohet_tests)
#     run_tests.run(test_suite0)

            
            
            
            
            
            
            