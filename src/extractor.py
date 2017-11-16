#!/usr/bin/env python

##############################################################################
##  phyhy: *P*ython *HyPhy*: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################
"""
    Parse (Extract!) JSON output from a standard HyPhy analysis.
"""

import os
import re
import json
import dendropy
from .analysis import *
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
        
        
        self.analysis_description         = "analysis"
        self.analysis_description_info    = "info"
        self.analysis_description_version = "version"
        
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
        self.original_name       = "original name"
        
        self.timers        = "timers"
        self.order         = "order"
        self.display_order = "display order"
        
        self.slac_by_site = "by-site"

        




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
        self.site_analyses             = [self.fel, self.fubar, self.meme, self.slac, self.leisr]
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
 
    

class Extractor():
    """
        This class parses JSON output and contains a variety of methods for pulling out various pieces of information.
    """    
    
    def __init__(self, content):
        """
            Initialize a Extractor instance.
            
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
        
        ## These are generally useful to have around ##
        self._count_partitions()          ### ---> self.npartitions        
        self._obtain_input_tree()         ### ---> self.input_tree, self.input_tree_dendropy
        self._obtain_fitted_models()      ### ---> self.fitted_models
        self._obtain_branch_attributes()  ### ---> self.branch_attributes, self.attribute_names
        self._obtain_original_names()     ### ---> self.original_names
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

        ### LEISR version error out ###
        version_field = self.json[ self.fields.analysis_description ][ self.fields.analysis_description_version ]
        assert( str(version_field) != "0.1alpha" ), "\n[ERROR]: LEISR analysis to parse was produced with HyPhy 2.3.6, which is not supported. Please re-analyze with version >=2.3.7 to use with phyphy."




    def _count_partitions(self):
        """
            Define self.npartitions, the number of partitions in analysis.
        """
        if self.analysis in self.analysis_names.single_partition_analyses:
            self.npartitions = 1
        else:
            self.npartitions = int(self.json[ self.fields.input ][ self.fields.input_npartitions ])
   

    def _obtain_input_tree(self):
        """
            Save the input tree(s) as either a string (single partition analysis), or as a dictionary (multiple partition analysis).
        """
        tree_field = self.json[ self.fields.input ][ self.fields.input_trees ]
        self.input_tree = {}
        self.input_tree_dendropy = {}
        for i in range(len(tree_field)):
            self.input_tree[i] = str(tree_field[str(i)]) + ";"
            self.input_tree_dendropy[i] = dendropy.Tree.get(data = self.input_tree[i], schema = "newick", preserve_underscores=True)     


    def _obtain_fitted_models(self):
        """
            Obtain list of all models in fits/attributes.
        """
        self.fitted_models = list( self.json[ self.fields.model_fits ].keys())      


    def _obtain_original_names(self):
        """
            Obtain original names dictionary, selecting only 0th partition.
        """
        self.original_names = self.extract_branch_attribute(self.fields.original_name, partition = 0)



    def _obtain_branch_attributes(self):
        """
            Obtain two things:
                - the full branch attributes dictionary (sans attributes part), self.branch_attributes
                - dictionary of attribute names, as attributes:attribute_type, self.attribute_names
        """
        raw = self.json[ self.fields.branch_attributes ]
        self.branch_attributes = {}
        for key in raw:
            try:
                self.branch_attributes[int(key)] = raw[key]
            except:
                pass
                
        self.attribute_names = {}
        for x in raw[ self.fields.attributes ]:
            if x == self.fields.display_order:
                continue
            else:
                self.attribute_names[x] = str(raw[ self.fields.attributes ][x][self.fields.attribute_type])      
        

    def _extract_slac_sitetable(self, raw):
        """
            Extract the specific SLAC tables of interest for parsing to CSV.
        """
        final = {}
        for x in range(self.npartitions):
            part = raw[str(x)]
            subset = part[self.fields.slac_by_site][self.slac_ancestral_type]
            final[str(x)] = subset
        return final          
       
       
    
    def _clean_meme_html_header(self, raw_header):
        """
            MEME has html tags all over it and this has to go.
            THIS IS VERY VERY HARDCODED, but flexible enough to not hurt much if someone changes the headers one day.
        """
        if raw_header[0] == "alpha;":
            raw_header[0] = "alpha"
        if raw_header[1] == "&beta;<sup>-</sup>":
            raw_header[1] = "beta_neg"
        if raw_header[2] == "p<sup>-</sup>":
            raw_header[2] = "prop_beta_neg"
        if raw_header[3] == "&beta;<sup>+</sup>":
            raw_header[3] = "beta_pos"
        if raw_header[4] == "p<sup>+</sup>":
            raw_header[4] = "prop_beta_pos"
        if raw_header[7] == "# branches under selection":
            raw_header[7] = "num branches under selection"
        
        return raw_header

        


    def _parse_sitemethod_to_csv(self, delim):
        """
            Extract a CSV from a **site-level** method JSON, including FEL, SLAC, MEME, FUBAR, LEISR.
        """
        site_block =  self.json[ self.fields.MLE ]
        raw_header = site_block[ self.fields.MLE_headers ]
        raw_header = [str(x[0]) for x in raw_header]
        raw_content = site_block[ self.fields.MLE_content]
        
        if self.analysis == self.analysis_names.slac:
            raw_content = self._extract_slac_sitetable(raw_content)
        if self.analysis == self.analysis_names.meme:
            raw_header = self._clean_meme_html_header(raw_header)

        final_header = "site,"+delim.join( [x.replace(" ","_") for x in raw_header] )
            
        if self.npartitions > 1:
            final_header = "partition," + final_header
        
        site_count = 1
        final_content = ""
        for i in range(self.npartitions):
            for row in raw_content[str(i)]:
                outrow = str(site_count) + delim + delim.join(str(x) for x in row)
                if self.npartitions > 1:
                    outrow = "\n" + str(i) + delim + outrow
                else:
                    outrow = "\n" + outrow
                final_content += outrow
                site_count += 1
        
        with open(self.csv, "w") as f:
            f.write(final_header + final_content)



    def _parse_absrel_to_csv(self, delim, original_names):
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
                
            if original_names is True:
                try:
                    outnode = str( d[self.fields.original_name] )
                except:
                    outnode = node
            else:
                outnode = node
                
                
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



    def _replace_tree_branch_length(self, dtree, bl_dict):
        """
            Replacing the branch length for a given node with the provided value (ie, replace <stuff> in :Node<stuff>)
                            
            Required arguments:
                1. **dtree**, the single dendropy tree to manipulate
                2. **bl_dict**, dictionary of new branch lengths
        """
        for node in dtree.preorder_node_iter():
            if node.parent_node is None:
                continue
            if node.is_internal():
                node.edge.length = bl_dict[node.label]
            else:
                node.edge.length = bl_dict[node.taxon.label]
        return dtree
        

       

    def _replace_node_labels(self, dtree, label_dict):
        """
            For a leaf, tack _<label> onto the name
            For an internal node, replace the node name with the label entirely
                            
            Required arguments:
                1. **dtree**, the single dendropy tree to manipulate
                2. **label_dict**, the dictionary of node:newthing
        """
        for node in dtree.preorder_node_iter():
            if node.parent_node is None:
                continue
            if node.is_internal():
                node.label = label_dict[node.label]
        return dtree
                
                
                

    def _tree_to_original_names(self, dtree):
        """
            Convert node names in a tree to original names
            
            Required arguments:
                1. **dtree**, the single dendropy tree to manipulate

        """
        for node in self.original_names:
            itsme = dtree.find_node_with_taxon_label(node)
            itsme.taxon.label = self.original_names[node]
        return dtree
    ############################################## PUBLIC FUNCTIONS ################################################### 


    ################################################ MODEL FITS #######################################################
    def reveal_model_names(self):
        """
            Reveal a list of all model names in the `fits` JSON field.
        """
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
        assert(model_name in self.fitted_models), "\n[ERROR]: Invalid model name."
        model_fit = self.json[ self.fields.model_fits ][ model_name ]
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
                raise AssertionError("\n[ERROR]: Unknown frequencies found in JSON.")            
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
    def reveal_branch_attributes(self):
        """
            Return a dictionary of all the attributes in the `branch attributes` field and their attribute type (node label or branch label).
        """
        return self.attribute_names

        

    def extract_input_tree(self, partition = None, original_names = False):
        """
            Return the inputted newick phylogeny, whose nodes have been labeled by HyPhy (if node labels were not present).
            For analyses with a single partition OR for a request for a specific partition's tree, returns a string.
            For analyses with multiple partitions (and hence multiple trees), returns a *dictionary* of trees. 
            
            Optional keyword arguments:
                1. **partition**, Integer indicating which partition's tree to return (as a string) if multiple partitions exist. NOTE: PARTITIONS ARE ORDERED FROM 0. This argument is ignored for single-partitioned analyses.
                2. **original_names**, Boolean (Default: False) if should update with original names before returning
        """
        
        if original_names is True:
            original_tree = {}
            for key in self.input_tree_dendropy:
                t = deepcopy(self.input_tree_dendropy[key])
                original_tree[key] = self._tree_to_original_names(t).as_string("newick", unquoted_underscores=True).strip()
        else:
            original_tree = self.input_tree
        if partition is None:
            if self.npartitions == 1:
                return original_tree[0]
            else:
                return original_tree
        else:
            return original_tree[int(partition)] 
    
        
    def extract_branch_attribute(self, attribute_name, partition = None):
        """

            Return dictionary of attributes for given attribute, where keys are nodes and values are attributes.
            If there are multiple partitions, default returns a dictionary with all partitions. 
            If partition = [some integer], only the attribute for the given partition will be returned. NOTE: PARTITION STARTS FROM 0. 
            
            Importantly, the values for all returned dictionaries will be **strings**, except for the extraction of rate distributions .
            
            Required positional arguments:
                1. **attribute_name**, the name of the attribute to obtain. Attribute names available can be revealed with the method `.reveal_branch_attributes()`.
                
            Optional keyword arguments:
                1. **partition**, Integer indicating which partition's tree to return (as a string) if multiple partitions exist. NOTE: PARTITIONS ARE ORDERED FROM 0. This argument is ignored for single-partitioned analyses.      
        """        
        assert(attribute_name in self.attribute_names), "\n[ERROR]: Specified attribute does not exist in JSON."
        if self.npartitions == 1:
            partition = 0
            
        attr_dict = {}
        for x in range(self.npartitions):
            partition_attr = {}
            for node in self.branch_attributes[x]:
                try:
                    attribute_value = str( self.branch_attributes[x][node][attribute_name] )
                    partition_attr[str(node)] = attribute_value
                except:
                    assert(attribute_name == self.fields.original_name), "\n[ERROR] Could not extract branch attribute."
                    pass  
            attr_dict[x] = partition_attr   
        if self.npartitions == 1:
            return attr_dict[0]
        else:
            if partition is None:
                return attr_dict
            else:
                return attr_dict[int(partition)]
        
        
        
        
    def map_branch_attribute(self, attribute_name, partition = None, original_names = False, as_label = False, update_branch_lengths = None):
        """
            Return the newick phylogeny with specified attribute mapped into the phylogeny.
            If there are multiple partitions, default returns a dictionary of mapped trees for all partitions. 
            If partition is specified, only the attribute for the given partition will be returned. NOTE: PARTITION STARTS FROM 0.            

            Required positional arguments:
                1. **attribute_name**, the name of the attribute to obtain. Attribute names available can be revealed with the method `.reveal_branch_attributes()`.
                
            Optional keyword arguments:
                1. **partition**, Integer indicating which partition's tree to return (as a string) if multiple partitions exist. NOTE: PARTITIONS ARE ORDERED FROM 0. This argument is ignored for single-partitioned analyses.      
                2. **original_names**, reformat the tree with the original names (as opposed to hyphy-friendly names with forbidden characters replaced). In most cases hyphy and original names are identical. Default: False.
                3. **as_label**, attributes should be added as labels onto the tree, as opposed to replacing the branch lengths. **This means that tips will have no information**  Default: False
                4. **update_branch_lengths**, string model name, indicting that branch lengths should be replaced with the given model fit's optimized lengths. Default: None. Note that this argument is ONLY USED when as_label = True.

        """
        assert(attribute_name != self.fields.rate_distributions), "\n[ERROR]: Cannot map rate distributions onto a tree."
        
        bl_dict = None
        if update_branch_lengths is not None:
            assert(update_branch_lengths in self.fitted_models and update_branch_lengths in self.reveal_branch_attributes()), "\n [ERROR]: Model with which to update branch length is not available."
            bl_dict = self.extract_branch_attribute(update_branch_lengths, partition = partition)
        dtree = deepcopy( self.input_tree_dendropy )
        
        mapped_trees = {}
        for key in dtree:
            t = dtree[key]
            attr_dict = self.extract_branch_attribute(attribute_name, partition = key)
            if as_label:
                if bl_dict is not None:   
                    t = self._replace_tree_branch_length( t, bl_dict )
                t = self._replace_node_labels( t, attr_dict )
            else: 
                t = self._replace_tree_branch_length( t, attr_dict )
            if original_names is True:
                t = self._tree_to_original_names(t)
            mapped_trees[key] = t.as_string("newick", unquoted_underscores=True).strip()    
        if self.npartitions == 1:
            return mapped_trees[0]
        else:
            if partition is None:
                return mapped_trees
            else:
                return mapped_trees[int(partition)]

 
 
    def extract_model_tree(self, model, partition = None, original_names = False):
        """
            Return newick phylogeny fitted to a certain model, i.e. with branch lengths optimized for specified model.
            This is just a special case of map_branch_attribute.

            Required positional arguments:
                1. **model**, the name of the model whose optimized tree you wish to obtain. Models names available can be revealed with the method `.reveal_model_names()`.
                
            Optional keyword arguments:
                1. **partition**, Integer indicating which partition's tree to return (as a string) if multiple partitions exist. NOTE: PARTITIONS ARE ORDERED FROM 0. This argument is ignored for single-partitioned analyses.      
                2. **original_names**, reformat the tree with the original names (as opposed to hyphy-friendly names with forbidden characters replaced). In most cases hyphy and original names are identical. Default: False.

        """
        return self.map_branch_attribute(model, partition = partition, original_names = original_names)
    ############################################################################################################################



    
    ################################################### MISCELLANEOUS ##########################################################
 

    def reveal_fields(self):
        """
            Return list of top-level JSON fields.
        """
        return [str(x) for x in list( self.json.keys() )]
        
        

    def extract_csv(self, csv, delim = ",", original_names = True, slac_ancestral_type = "AVERAGED"):
        """
            
            Extract a CSV from JSON, for certain methods:
                + FEL
                + SLAC
                + MEME
                + FUBAR
                + LEISR
                + aBSREL
                
            Required positional arguments:
                1. **csv**, File name for output CSV
                
            Optional keyword arguments:
                1. **delim**, A different delimitor for the output, e.g. "\t" for tab
                2. **original_names**, An *ABSREL* specific boolean argument to indicate whether HyPhy-reformatted branch should be used in output csv (False), or original names as present in the input data alignment should be used (True). Default: True
                3. **slac_ancestral_type**, A *SLAC* specific argument, either "AVERAGED" (Default) or "RESOLVED" (case insensitive) to indicate whether reported results should be from calculations done on either type of ancestral counting.

        """       
        
        self.csv = csv
        
        ### FEL, MEME, SLAC, FUBAR, LEISR ###
        if self.analysis in self.analysis_names.site_analyses:
            slac_ancestral_type = slac_ancestral_type.upper() 
            assert(slac_ancestral_type in self.analysis_names.slac_ancestral_type), "\n[ERROR]: Argument `slac_ancestral_type` must be either 'AVERAGED' or 'RESOLVED' (case insensitive)."
            self.slac_ancestral_type = slac_ancestral_type
            self._parse_sitemethod_to_csv(delim)
       
        ### aBSREL ###
        elif self.analysis == self.analysis_names.absrel:
            assert(type(original_names) == bool), "\n[ERROR]: Argument `original_names` must be boolean."
            self._parse_absrel_to_csv(delim, original_names)

            
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
    
