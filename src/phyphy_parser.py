#!/usr/bin/env python

##############################################################################
##  phyhy: *P*ython *HyPhy*: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################
"""
    Parse JSON output from a standard HyPhy analysis.
"""

import os
import re
import json
import pyvolve
from .phyphy_runner import *
from copy import deepcopy

class JSONFields():
    """
        This class defines the strings of relevant JSON keys. 
        Note that these strings correspond precisely to those in the HyPhy distribution. See file: `TemplateBatchFiles/libv3/all-terms.bf` in the `terms.json` namespace.
    """
    
    def __init__(self):
        
        self.input                = "input"
        self.input_info           = "info"
        self.input_filename       = "file name"
        self.input_sites          = "number of sites"
        self.input_sequences      = "number of sequences"
        self.input_npartitions    = "partition count"
        self.input_trees          = "trees"
        
        
        self.analysis_description      = "analysis"
        self.analysis_description_info = "info"
        
        self.substitution_rate    = re.compile(r"Substitution rate from [\w-]+ (\w) to [\w-]+ (\w)")
        
        self.model_fits           = "fits"
        self.log_likelihood       =  "Log Likelihood"
        self.aicc                 = "AIC-c"
        self.estimated_parameters = "estimated parameters"
        self.frequencies          = "Equilibrium frequencies"
        self.nucleotide_gtr       = "Nucleotide GTR"    
        self.generic_mg94xrev     = "MG94xREV"      
        self.per_branch_omega     = "Per-branch omega"
        self.omega                = "omega"
        self.proportion           = "proportion"
        self.rate_distributions   = "Rate Distributions"
        self.nonsyn_syn_ratio_for = "non-synonymous/synonymous rate ratio for"

        self.MLE                       = "MLE"
        self.MLE_headers               = "headers"
        self.MLE_content               = "content"
        self.relative_site_rates       = "Relative site rate estimates"
        self.UB                        = "UB"
        self.LB                        = "LB"

        self.tested          = "tested"
        
        self.site_logl       = "Site Log Likelihood"
        self.evidence_ratios = "Evidence Ratios"        
        
        self.LRT            = "LRT"
        self.uncorrected_p  = "Uncorrected P-value"
        self.corrected_p    = "Corrected P-value"
        self.baseline_omega = "Baseline MG94xREV omega ratio"
        self.rate_classes   = "Rate classes"
   
        self.branch_attributes   = "branch attributes"
        self.attributes          = "attributes"
        self.attribute_type      = "attribute type"
        
        self.timers        = "timers"
        self.order         = "order"
        self.display_order = "display order"

        




class AnalysisNames():
    """
        This class defines the names of analyses which we can parse.
    """
    def __init__(self):
        self.absrel   = "ABSREL"
        self.busted   = "BUSTED"
        self.fel      = "FEL"
        self.fubar    = "FUBAR"
        self.leisr    = "LEISR"
        self.meme     = "MEME"
        self.relax    = "RELAX"
        self.slac     = "SLAC"
        
        self.all_analyses              = [self.absrel, self.busted, self.fel, self.fubar, self.leisr, self.meme, self.relax, self.slac]
        self.site_analyses             = [self.fel, self.fubar, self.meme, self.slac]
        self.single_partition_analyses = [self.absrel, self.leisr, self.relax]

        self.slac_by = ["by-site", "by-branch"]
        self.slac_ancestral_type = ["AVERAGED", "RESOLVED"]


class Genetics():
    """
        Class to define codes used. 
        Primarily (only?) used to extract frequencies as dictionaries.
    """
    def __init__(self):
        self.nucleotides  = ["A", "C", "G", "T"]
        self.amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        self.codons       = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.genetics     = {4: self.nucleotides, 20: self.amino_acids, 61: self.codons}
 
    

class HyPhyParser():
    """
        This class parses JSON output and contains a variety of methods for pulling out various pieces of information.
    """    
    
    def __init__(self, content):
        """
            Initialize a HyPhyParser instance.
            
            Required arguments:
                1. **content**, The input content to parse. Two types of input may be provided here:
                    + The path to a JSON file to parse, provided as a string
                    + A phyphy `Analysis` (i.e. `BUSTED`, `SLAC`, `FEL`, etc.) object which has been used to execute a HyPhy analysis through the phyphy interface
        """
        self.fields = JSONFields()
        self.genetics = Genetics()
        
        self.analysis_names = AnalysisNames()
        self.allowed_analyses = self.analysis_names.all_analyses
        
        ### Input ###
        if type(content) == str:
            self.json_path = content
        elif isinstance(content, Analysis):
            self.json_path = content.final_path
        else:
            raise AssertionError("\n[ERROR]: Expected a single argument. Provide either the path to the JSON to parse, or an `Analysis` object which has been executed.")
        assert(os.path.exists(self.json_path)), "\n[ERROR]: JSON file does not exist."
        
        self._unpack_json()
        self._determine_analysis_from_json()
        self._count_partitions()             ## This seems to be generally useful


    ############################## PRIVATE FUNCTIONS #################################### 
    def _unpack_json(self):
        """
            Unpack JSON into dictionary
        """ 
        self.json = None
        with open (self.json_path, "rU") as f:
            self.json = json.load(f)
        assert(self.json is not None and len(self.json)!=0), "\n[ERROR]: Unable to obtain JSON contents."

    
    def _determine_analysis_from_json(self):
        """
            Determine the relevant analysis name directly from the JSON description field.
        """

        json_info = self.json[ self.fields.analysis_description ][ self.fields.analysis_description_info ].upper()
        
        for name in self.allowed_analyses:
        
            find_analysis = re.search(name.upper(), json_info)
            if find_analysis is not None:
                self.analysis = name
                break
        
        assert(self.analysis is not None), "\n[ERROR]: Could not determine analysis from JSON. Please ensure that the JSON is correctly formatted and created with HyPhy version >=2.3.4."



    def _count_partitions(self):
        """
            Define self.npartitions, the number of partitions in analysis.
        """
        if self.analysis in self.analysis_names.single_partition_analyses:
            self.npartitions = 1
        else:
            self.npartitions = int(self.json[ self.fields.input ][ self.fields.input_npartitions ])
   


    def _extract_slac_sitetable(self, raw, slac_by, slac_ancestral_type):
        """
            Extract the specific SLAC tables of interest for parsing to CSV.
        """
        final = {}
        for x in range(self.npartitions):
            part = raw[str(x)]
            subset = part[slac_by][slac_ancestral_type]
            final[str(x)] = subset
        return final          
            
        


    def _parse_sitemethod_to_csv(self, delim, slac_by = "by-site", slac_ancestral_type = "AVERAGED"):
        """
            Extract a CSV from a **site-level** method JSON, including FEL, SLAC, MEME, FUBAR.
        """
        site_block =  self.json[ self.fields.MLE ]
        raw_header = site_block[ self.fields.MLE_headers ]
        raw_content = site_block[ self.fields.MLE_content]
        if self.analysis == self.analysis_names.slac:
            raw_content = self._extract_slac_sitetable(raw_content, slac_by, slac_ancestral_type)
            
        final_header = "site,"+delim.join( [x[0].replace(" ","_") for x in raw_header] )
        if self.npartitions > 1:
            final_header = "partition," + final_header
        
        site_count = 1
        final_content = ""
        for part in raw_content:
            for row in raw_content[part]:
                row = "\n" + str(site_count) + delim + delim.join(str(x) for x in row)
                if self.npartitions > 1:
                    row = str(part) + delim + row
                final_content += row
                site_count += 1
        
        with open(self.csv, "w") as f:
            f.write(final_header + final_content)



    def _parse_absrel_to_csv(self, delim):
        """
            Extract a CSV from an aBSREL JSON. 
            CSV contents:
                Node name, Baseline MG94 omega, Number of inferred rate classes, Tested (bool), Proportion of selected sites, LRT, uncorrected P, bonferroni-holm P
        """
        
        header = delim.join( ["node", "baseline_omega", "number_rate_classes", "tested", "prop_sites_selected", "LRT", "uncorrected_P", "corrected_P"] )
        attr = self.json[ self.fields.branch_attributes ]["0"] ## Only allowed single partition for ABSREL
        node_names = list( self.json[ self.fields.branch_attributes ]["0"].keys())  
        
        full_rows = ""
        for node in node_names:
            
            try:
                d = attr[str(node)]
            except:
                raise KeyError("\n[ERROR]: Unable to parse JSON.")
            rates = d[self.fields.rate_distributions]
            if len(rates) > 1:
                for pair in rates:
                    if pair[0] > 1.:
                        prop = str(pair[1])
                        break
            else:
                prop = "0"
            
            if d[ self.fields.LRT ] == 1 and d[ self.fields.uncorrected_p ] ==  1 and d[ self.fields.corrected_p ] == 1:
                run = "0"
            else:
                run  = "1"
                
            row = "\n" + delim.join([node, 
                                    str(d[self.fields.baseline_omega]), 
                                    str(d[self.fields.rate_classes]), 
                                    run,
                                    prop,
                                    str(d[self.fields.LRT]), 
                                    str(d[self.fields.uncorrected_p]), 
                                    str(d[self.fields.corrected_p]) ])
            full_rows += row
        
        with open(self.csv, "w") as f:
            f.write(header + full_rows)


    def _parse_leisr_to_csv(self, delim):
        """
            Extract a CSV from a LEISR analysis 
            CSV contents:
                site, rate, lower95, upper95
                
            Required arguments:
                1. **delim**, the delimitor for output file. Default: "," (i.e., CSV).
        """
        
        header = delim.join( ["site", "rate", "lower_95_bound", "upper_95_bound"] )
        attr = self.json[ self.fields.relative_site_rates ]
        nsites = self.json[ self.fields.input ][ self.fields.input_sites ]
        
        full_rows = ""
        for i in range(nsites):
            
            site = str(i+1)
            row = "\n" + delim.join( [ site, str(attr[site][self.fields.MLE]), str(attr[site][self.fields.LB]), str(attr[site][self.fields.UB]) ])
            full_rows += row
         
        with open(self.csv, "w") as f:
            f.write(header + full_rows)
       
        
    
    def _reform_rate_phrase(self, phrase):
        """
            Convert rate phrase to simpler key, i.e. "Substitution rate from nucleotide A to nucleotide C" returns simply "AC"
                
            Required arguments:
                1. **phrase**, the key to reform
        """
        find = re.search(self.fields.substitution_rate, phrase)
        if find:
            source = find.group(1).upper()
            target = find.group(2).upper()
        
            return str(source + target)
        else:
            raise AssertionError("\n[ERROR]: Bad rate reform.")



    def _replace_tree_info(self, tree, node, newvalue):
        """
            Faciliate attribute mapping onto a phylogeny by replacing the branch length for a given node with the provided value (ie, replace <stuff> in :Node<stuff>)
                            
            Required arguments:
                1. **tree**, the newick tree string of interest
                2. **node**, the node node whose branch length should be replaced
                3. **newvalue**, the new value to replace the old branch length
        """
        return re.sub(node + ":[\.\de-]+?([,\)])", node + ":" + newvalue + "\\1", tree)
        


    ############################################## PUBLIC FUNCTIONS ################################################### 


    ################################################ MODEL FITS #######################################################
    def extract_model_names(self):
        """
            Return a list of all model names in the `fits` JSON field.
        """

        self.fitted_models = list( self.json[ self.fields.model_fits ].keys())
        return self.fitted_models


    def extract_model_component(self, model_name, component):
        """
            Return a model component for a given model name found in the `fits` JSON field. 
            
            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
                2. **component**, the component of the model to return. 

            Note there are a variety of convenience methods which wrap this function to extract all components (note that not all analyses will have all of these components):
                + .extract_model_logl(model_name) returns the log likelihood of a given model fit
                + .extract_model_estimated_parameters(model_name) returns the number of estimated parameters in a given model fit
                + .extract_model_aicc(model_name) returns the small-sample AIC (AIC-c) for a given model fit
                + .extract_model_rate_distributions(model_name) returns rate distributions for a given model fit 
                + .extract_model_frequencies(model_name) returns the equilibrium frequencies for the given model fit
        """            
        try:
            model_fit = self.json[ self.fields.model_fits ][ model_name ]
        except:
            raise KeyError("\n[ERROR]: Invalid model name.")
            
        try:
            component = model_fit[component]
        except: 
            raise KeyError("\n[ERROR]: Invalid model component.")
            
        return component


    def extract_model_logl(self, model_name):
        """
            Return log likelihood (as a float) for a given model that appears in the the `fits` field.

            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
        """
        return float( self.extract_model_component(model_name, self.fields.log_likelihood) )


    def extract_model_estimated_parameters(self, model_name):
        """
            Return estimated parameters (as an int) for a given model that appears in the `fits` field.

            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
        """
        return int( self.extract_model_component(model_name, self.fields.estimated_parameters) )


    def extract_model_aicc(self, model_name):
        """
            Return AICc (as a float) for a given model that appears in the `fits` field.
            
            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
        """
        return float( self.extract_model_component(model_name, self.fields.aicc) )


    def extract_model_rate_distributions(self, model_name):
        """
            Return rate distributions, as a reformatted dictionary, for a given model that appears in the `fits` field.            
            NOTE: Currently assumes dS = 1 for all initial MG94xREV fits, as in the current HyPhy implementation (True in <=2.3.4).

            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
        """
        rawrates = self.extract_model_component(model_name, self.fields.rate_distributions)
        rates = {}
        
        if model_name == self.fields.nucleotide_gtr: 
            for k,v in rawrates.items():
                rates[ self._reform_rate_phrase(k) ] = v       
        
        elif self.fields.generic_mg94xrev in model_name:
            if self.analysis == self.analysis_names.absrel:
                rates = rawrates[self.fields.per_branch_omega]   
            else:
                rates = {}
                for k,v in rawrates.items():
                    find = re.search(r""+self.fields.nonsyn_syn_ratio_for +" \*(\w+)\*", k)
                    if find:
                        rates[find.group(1)] = {self.fields.omega: v[0][0], self.fields.proportion: 1.0}
                    else:
                        rates[self.fields.omega] = v[0][0]  
        else:
            rates = rawrates
        return rates



    def extract_model_frequencies(self, model_name, as_dict = False):
        """
            Return a list of equilibrium frequencies (in alphabetical order) for a given model that appears in the field `fits`.
            
            Required arguments:
                1. **model_name**, the name of the model of interest. Note that all model names can be revealed with the method `.extract_model_names()`
            Optional keyword arguments:
                1. **as_dict**, Boolean to indicate if the frequencies should be returned as a dictionary. Default: False.
        """
        fraw = self.extract_model_component(model_name, self.fields.frequencies)
        f = [float(x[0]) for x in fraw]
        
        if as_dict:
            try:
                fdict = dict(zip( self.genetics.genetics[len(f)], f))
            except:
                raise AssertionError("\n[ERROR]: Unknown frequencies found in JSON. Please report bug to github.com/veg/hyphy/issues.")            
            return fdict
        else:
            return f   
    ###################################################################################################################
    

     
    ################################################ BRANCH SETS ######################################################
  
    def extract_branch_sets(self, by_set = False):
        """
            Return branch set designations as a dictionary for all nodes. 
            By default, this function will return the branch sets "as is" is the JSON field `tested`, where keys are node and values are the branch set to which the given node belongs
            NOTE: Assumes that all partitions share the same branch sets.
            
            Optional keyword arguments:
                1. **by_set**, Boolean to indicate if the returned dictionary should use *branch sets* as keys, and values are a *list* of nodes in that branch set. Default: False.
        """
        try:
            branch_sets = self.json[ self.fields.tested ]["0"]
        except:
            raise KeyError("\n[ERROR]: Provided JSON has no branch set designations")
            
        if not by_set:
            return branch_sets
        else:
            branch_sets2 = {}
            for k,v in branch_sets.items():
                if v in branch_sets2:
                    branch_sets2[v].append(k)
                else:
                    branch_sets2[v] = [k]
            return branch_sets2
     ###################################################################################################################



    ######################################## TREE WITH VARIOUS ATTRIBUTES ##############################################
    
    def extract_input_tree(self, partition = None):
        """
            Return the inputted newick phylogeny, whose nodes have been labeled by HyPhy (if node labels were not present).
            For analyses with a single partition, returns a string.
            For analyses with multiple partitions (and hence multiple trees), returns a *list* of trees. 
            
            
            Return branch set designations as a dictionary for all nodes. 
            By default, this function will return the branch sets "as is" is the JSON field `tested`, where keys are node and values are the branch set to which the given node belongs
            NOTE: Assumes that all partitions share the same branch sets.
            
            Optional keyword arguments:
                1. **partition**, Integer indicating which partition's tree to return (as a string) if multiple partitions exist. NOTE: PARTITIONS ARE ORDERED FROM 0.
        """
        tree_field = self.json[ self.fields.input ][ self.fields.input_trees ]
        if self.npartitions == 1:
            self.input_tree = str(tree_field["0"]) + ";"
            return self.input_tree
        
        else:
            self.input_tree = []
            for i in range(len(tree_field)):
                self.input_tree.append( str(tree_field[str(i)]) + ";" )
        
            if partition is not None:
                try:
                    return self.input_tree[partition]
                except:
                    raise KeyError("\n[ERROR]: Partition not found. Note that partitions are enumerated starting from 0.")
            else:
                return self.input_tree

    
    def reveal_branch_attributes(self):
        """
            Return a dictionary of all the attributes in the `branch attributes` field and their attribute type (node label or branch label).
        """
        attributes_info = self.json[ self.fields.branch_attributes ][ self.fields.attributes ]
        self.attribute_names = {}
        for x in attributes_info:
            if x == self.fields.display_order:
                continue
            else:
                self.attribute_names[str(x)] = str(attributes_info[x][self.fields.attribute_type])
        return self.attribute_names
        

        
    
    def extract_branch_attribute(self, attribute_name, partition = None):
        """
            Return dictionary of attributes for given attribute, where keys are nodes and values are attributes.
            If there are multiple partitions, default returns a dictionary with all partitions. 
            If partition = [some integer], only the attribute for the given partition will be returned. NOTE: PARTITION STARTS FROM 0. 
            
            Importantly, the values for all returned dictionaries will be **strings**, except for the extraction of rate distributions .
                      
        """
        self.reveal_branch_attributes() ## Needed to create self.attribute_names
        assert(attribute_name in self.attribute_names), "\n[ERROR]: Specified attribute does not exist in JSON."

        total_attr_dict = {}
        for x in range(self.npartitions):
            attr_dict = {}
            partition_attributes = self.json[ self.fields.branch_attributes ][str(x)]      
            for node in partition_attributes:
                attribute_value = str( partition_attributes[node][attribute_name] )
                attr_dict[str(node)] = attribute_value   
            total_attr_dict[str(x)] = attr_dict   

        if partition is not None:
            try:
                return total_attr_dict[str(partition)]
            except:
                raise KeyError("\n[ERROR]: Partition not found. Note that partitions are enumerated starting from 0.")
        else:
            ## Ignore a partition argument if there is only 1 partition.
            if self.npartitions == 1:
                return attr_dict
            else:
                return total_attr_dict
        
        
        
        
        
    def map_branch_attribute(self, attribute_name, partition = None):
        """
            Return an attribute mapped onto the newick phylogeny.
            If there are multiple partitions, default returns a list of mapped for all partitions. 
            If partition = [some integer], only the attribute for the given partition will be returned. NOTE: PARTITION STARTS FROM 0.            
        """
        assert(attribute_name != self.fields.rate_distributions), "\n[ERROR]: Cannot map rate distributions onto a tree."
                
        attr_dict = self.extract_branch_attribute(attribute_name)
        self.extract_input_tree()       ## Needed to grab input tree (grab all)
        ptree = deepcopy( self.input_tree )
        
        if self.npartitions == 1:
            for node in attr_dict:
                ptree = self._replace_tree_info( ptree, node, str(attr_dict[node]) )
            return ptree
        
        else:
            many_trees = []
            for ptree in self.input_tree:
                for node in attr_dict:
                    ptree = self._replace_tree_info( ptree, node, str(attr_dict[node]) )
                many_trees.append(ptree)
        
        if partition is not None:
            try:
                return many_trees[partition]
            except:
                raise KeyError("\n[ERROR]: Partition not found. Note that partitions are enumerated starting from 0.")
        else:
            return many_trees
    
 
    def extract_model_tree(self, model, partition = None):
        """
            Return newick phylogeny fitted to a certain model, i.e. with branch lengths optimized for specified model.
            This is just a special case of map_branch_attribute.
        """
        return self.map_branch_attribute(model, partition = partition)
    ############################################################################################################################



    
    ################################################### MISCELLANEOUS ##########################################################
 

    def reveal_fields(self):
        """
            Return list of top-level JSON fields.
        """
        return [str(x) for x in list( self.json.keys() )]
        
        

    def extract_csv(self, csv, delim = ",", slac_by = "by-site", slac_ancestral_type = "AVERAGED"):
        """
            Extract results to a CSV. Currently only for SLAC, MEME, FEL, FUBAR, aBSREL, and relative rates. Other analyses do not lend well to CSV.
            Note that we don't export BUSTED (could be site logl or ERs) because these quantities are not actual tests for selection and users should not be readily able to treat as such. They can grab on their own if that interested.
        """       
        
        self.csv = csv

        ### Standard site output ###
        if self.analysis in self.analysis_names.site_analyses:
            assert(slac_by in self.analysis_names.slac_by), "\n[ERROR]: Argument `slac_by` must be either 'by-site' or 'by-branch'."
            assert(slac_ancestral_type in self.analysis_names.slac_ancestral_type), "\n[ERROR]: Argument `slac_ancestral_type` must be either 'AVERAGED' or 'RESOLVED'."
            self._parse_sitemethod_to_csv(delim)
       

        ### aBSREL ###
        elif self.analysis == self.analysis_names.absrel:
            self._parse_absrel_to_csv(delim)
        
        ## LEISR ##
        elif self.analysis == self.analysis_names.leisr:
            self._parse_leisr_to_csv(delim)
            
        else:
            print("\nContent from provided analysis is not convertable to CSV.")
 
 
    def extract_timers(self):
        """
            Extract dictionary of timers, with display order removed
        """
        raw = self.json[self.fields.timers]
        final = {}
        for step in raw:
            del raw[step][self.fields.order]
            for k,v in raw[step].items():
                final[str(step)] = float(v)
        return final
        
 

    def extract_site_logl(self):
        """
            Extract BUSTED site log likelihoods, as dictionary
        """
        assert(self.analysis == self.analysis_names.busted), "\n[ERROR]: Site Log Likelihoods are specific to BUSTED."
        
        raw = self.json[self.fields.site_logl]
        site_logl = {}
        for k,v in raw.items():
            site_logl[str(k)] = v[0]
        
        return site_logl
    
    
    def extract_evidence_ratios(self):
        """
            Extract BUSTED ERs, as dictionary
        """                    
        assert(self.analysis == self.analysis_names.busted), "\n[ERROR]: Site Log Likelihoods are specific to BUSTED."
        raw = self.json[self.fields.evidence_ratios]
        if len(raw) == 0:
            print("\n[Warning] Evidence ratios are only computed for BUSTED models with significant tests for selection. Note further that they should be interpretted only as **descriptive** measures of selection, NOT statistical tests.")
            return None
        else:
            ev_ratios = {}
            for k,v in raw.items():
                ev_ratios[str(k)] = v[0]
        return ev_ratios
    ###################################################################################################################
    
