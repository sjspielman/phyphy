#!/usr/bin/env python

##############################################################################
##  phyhy: *P*ython *HyPhy*: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################

	
	
"""
    Execute a standard HyPhy analysis.
"""
    
import sys
import subprocess
import os
import shutil
import re
from Bio import Phylo
from copy import deepcopy
try:
    from io import StringIO
except:
    from StringIO import StringIO
from math import ceil

if __name__ == "__main__":
    print("\nThis is the Analysis module in `phyphy`. Please consult docs for `phyphy` usage." )
    sys.exit()

from .hyphy import *


_GENETIC_CODE = {
                  1: "Universal",
                  2: "Vertebrate mtDNA",
                  3: "Yeast mtDNA",
                  4: "Mold/Protozoan mtDNA",
                  5: "Invertebrate mtDNA",
                  6: "Ciliate Nuclear",
                  9: "Echinoderm mtDNA",
                  10: "Euplotid Nuclear",
                  12: "Alt. Yeast Nuclear",
                  13: "Ascidian mtDNA",
                  14: "Flatworm mtDNA",
                  15: "Blepharisma Nuclear",
                  16: "Chlorophycean mtDNA",
                  21: "Trematode mtDNA",
                  22: "Scenedesmus obliquus mtDNA",
                  23: "Thraustochytrium mtDNA",
                  24: "Pterobranchia mtDNA",
                  25: "SR1 and Gracilibacteria",
                  26: "Pachysolen Nuclear"
                }


        

class Analysis(object):
       
            
    def __init__(self, **kwargs):
        """
            Parent class for all analysis methods, which include the following children:
        
            + ABSREL
            + BUSTED
            + FEL
            + FUBAR
            + LEISR
            + MEME
            + RELAX
            + SLAC                        

            **Do not use this parent class. Instead, see child classes for analysis-specific arguments and examples.**

        """
    
        self.hyphy = kwargs.get("hyphy", None)
        if self.hyphy is None:
            self.hyphy = HyPhy()
        
        
        self.alignment = kwargs.get("alignment", None) ### alignment only
        self.tree      = kwargs.get("tree", None)    ### tree only
        self.data      = kwargs.get("data", None)    ### combined alignment and tree or NEXUS
        self._check_files()

        self.user_json_path = kwargs.get("output", None)
        if self.user_json_path is not None:
            dirname = os.path.dirname(os.path.abspath(self.user_json_path))
            assert( os.path.exists(dirname) ),"\n[ERROR]: Provided output path does not exist."
                    
        self.genetic_code = kwargs.get("genetic_code", "Universal")
        assert(self.genetic_code in list(_GENETIC_CODE.values()) or self.genetic_code in list(_GENETIC_CODE.keys())), "\n[ERROR] Incorrect genetic code specified. Consult NCBI for options: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. \nNOTE that HyPhy supports only options up to and including 26."
        if self.genetic_code in list(_GENETIC_CODE.keys()):
            for k,v in _GENETIC_CODE.items():
                if k == self.genetic_code:
                    self.genetic_code = str(v)
                    break
        self.genetic_code = self.genetic_code.replace(" ", "\ ") # Sigh.
        
        
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/SelectionAnalyses/"
        self.shared_branch_choices = ("All", "Internal", "Leaves", "Unlabeled branches")
        
        ######## 2.3.7 models #######
        self.available_protein_models = ("JC69", "WAG", "LG", "JTT", "mtMAM", "cpREV", "HIVBm", "HIVWm", "AB")
        self.available_nucleotide_models = ("GTR", "HKY85", "JC69")
    
    
    
    
    def _format_yesno(self, argument):
        """
            Format argument to be Yes/No from True/False
        """
        self.yesno_truefalse = {True: "Yes", False: "No"}
        
        if type(argument) == str:
            argument = argument.capitalize()
            assert(argument in list(self.yesno_truefalse.values())),"\n[ERROR]: Incorrect Yes/No argument."
        elif type(argument) is bool:
            argument = self.yesno_truefalse[argument]
        else:
            raise TypeError("\n[ERROR]: Incorrect Yes/No argument.")
        return argument            
    
    
    def _check_files(self):
        """
            Check provided paths for alignment+tree or data. Assign input hyphy variables accordingly.
            Additionally extract the tree string for use in label finding.
        """
        
        assert(self.data is not None or (self.alignment is not None and self.tree is not None)), "\n[ERROR]: You must supply argument `data` (file with alignment and tree) OR arguments `alignment` and `tree` (separate files containing respective contents)."
        if self.alignment is not None:
            assert(os.path.exists(self.alignment)), "\n[ERROR] Provided alignment not found, check path?"
            assert(os.path.exists(self.tree)), "\n[ERROR] A tree must be provided. As needed, check path?"
            self.hyphy_alignment = os.path.abspath(self.alignment)
            self.hyphy_tree      = os.path.abspath(self.tree)
            with open(self.hyphy_tree, "r") as f:
                self.tree_string = f.read().strip()
        else:
            assert(os.path.exists(self.data)), "\n[ERROR] Provided data not found, check path?"
            self.hyphy_alignment = os.path.abspath(self.data)
            try:
                ## NOTE: We must use biopython, *not* dendropy here, because HyPhy uses curly braces for labels and dendropy will not parse a tree with labels formatted in this manner that HyPhy requires
                t = Phylo.read(self.hyphy_alignment, "nexus") ## If no error, tree is there and there will be no prompt 
                tree_handle = StringIO()
                Phylo.write(t, tree_handle, "newick")
                self.tree_string = tree_handle.getvalue().strip()
                self.hyphy_tree = ""            
            
            except: 
                self.hyphy_tree = "Y" # Use the tree found in the file
                with open(self.hyphy_alignment, "r") as f:
                    alnstring = f.read()
                    find_tree = re.search(r"(\(.+\);)", alnstring)
                    if find_tree:
                        self.tree_string = find_tree.group(1)
                    else:
                        raise AssertionError("\n[ERROR] Malformed or missing tree in input data.")


   
    def _sanity_branch_selection(self):
        """
            Ensure appropriate value provided for branch selection.
        """
        self._find_all_labels()
        allowed = list(self.shared_branch_choices) + self._all_labels
        assertion = "\n[ERROR]: Bad branch selection. Must be one of: " + ", ".join(["'"+str(x)+"'" for x in allowed]) + "."
        assert(self.branches in allowed), assertion
            
        

    def _find_all_labels(self):
        """
            Parse the tree string to find all the labels.
        """
        ## since all characters accepted inside {}, march along the tree to grab each one
        self._all_labels = []
        label = ""
        curly = False
        for i in range(len(self.tree_string)):
            if self.tree_string[i] == "{":
                curly = True
                continue 
            if self.tree_string[i] == "}":
                curly = False
                if label not in self._all_labels:
                    self._all_labels.append(label)
                label = ""
            if curly:
                label += self.tree_string[i]
            
    def _build_analysis_command(self):
        print("Parent method. Not run.")


    def _build_full_command(self):
        """
            Construct full command to execute analysis.
        """    
        self._build_analysis_command()
        self.run_command = " ".join([self.hyphy.hyphy_call, self.analysis_command])


    def _execute(self):
        """
            Execute HyPhy
        """

        if self.hyphy.quiet:
            quiet = open("/dev/null", "w")
            check = subprocess.call(self.run_command, shell = True, stdout = quiet, stderr = subprocess.STDOUT)
        else:    
            check = subprocess.call(self.run_command, shell = True)
        assert(check == 0), "\n[ERROR] HyPhy failed to run."


    def run_analysis(self):
        """
            Execute an Analysis and save output.

            **Examples:**
               
               >>> ### Execute a default FEL analysis
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat")
               >>> myfel.run_analysis()         
        """
        self._execute()
        self._save_output() 





    def _save_output(self):
        """
            Move JSON to final location. 
        """        
        if self.user_json_path is not None:       
            self.final_path = self.user_json_path
        else:
            self.final_path = self.default_json_path   
        shutil.move(self.default_json_path, self.final_path)




class FEL(Analysis):

    def __init__(self, **kwargs):
        """
            
            Initialize and execute a FEL analysis.
            
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **srv**, Employ synonymous rate variation in inference (i.e. allow dS to vary across sites?). Values "Yes"/"No" or True/False accepted. Default: True.
                3. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** in your tree are accepted
                4. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                5. **alpha**, The p-value threshold for calling sites as positively or negatively selected. Note that this argument has 0 bearing on JSON output. Default: 0.1
                6. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.


            **Examples:**
               
               >>> ### Define a default FEL analysis, where data is contained in a single file
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default FEL analysis, where alignment and tree are in separate files 
               >>> myfel = FEL(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")

               >>> ### Define a FEL analysis, with a specified path to output JSON
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat", output="/path/to/json/output.json")
               
               >>> ### Define a FEL analysis, with a one-rate approach (i.e. synonymous rate variation turned off) 
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat", srv=False)

               >>> ### Define FEL analysis, specifying only to use internal branches to test for selection
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat", branches="Internal")
               
               >>> ### Define FEL analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myfel = FEL(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)
               
               >>> ### Execute a defined FEL instance
               >>> myfel.run_analysis()
        """                
        
        super(FEL, self).__init__(**kwargs)
        
        self.batchfile = "FEL.bf"
        self.default_json_path = self.hyphy_alignment + ".FEL.json"

        self.srv   = kwargs.get("srv", True) ## They can provide T/F or Yes/No
        self.srv   = self._format_yesno(self.srv)
        self.alpha = str( kwargs.get("alpha", 0.1) )

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
        self._build_full_command()

        
    def _build_analysis_command(self):
        """
            Construct the FEL command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree ,
                                           self.branches , 
                                           self.srv , 
                                           self.alpha ])

        
      
class FUBAR(Analysis):       
        
    def __init__(self, **kwargs):
        """

            Initialize and execute a FUBAR analysis.

            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data. 
                3. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
                4. **grid_size**, Number of grid points per rate grid dimension (Default: 20, allowed [5,50])
                5. **nchains**, Number of MCMC chains to run (Default: 5, allowed [2,20])
                6. **chain_length**, The length of each chain (Default: 2e6, allowed [5e5,5e7])
                7. **burnin**, Number of samples to consider as burn-in (Default 1e6, allowed [ceil(chain_length/20),ceil(95*chain_length/100)])
                8. **samples_per_chain**, Number of samples to draw per chain (Default 100, allowed [50,chain_length-burnin])
                9. **alpha**, The concentration parameter of the Dirichlet prior (Default 0.5, allowed[0.001,1])
                10. **cache**, Name (and path to) output FUBAR cache. Default: goes to same directory as provided data. Provide the argument **False** to not save the cache (this argument simply sends it to /dev/null)


            **Examples:**
               
               >>> ### Define a default FUBAR analysis, where data is contained in a single file
               >>> myfubar = FUBAR(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default FUBAR analysis, where alignment and tree are in separate files 
               >>> myfubar = FUBAR(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")
               
               >>> ### Define a FUBAR analysis using a 10x10 grid and alpha parameter of 0.75
               >>> myfubar = FUBAR(data = "/path/to/data_with_tree.dat", grid_size = 10, alpha = 0.75 )

               >>> ### Define FUBAR analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myfubar = FUBAR(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)

               >>> ### Execute a defined FUBAR instance
               >>> myfubar.run_analysis()
        """                



        super(FUBAR, self).__init__(**kwargs)
        
        self.batchfile = "FUBAR.bf"
        self.default_json_path = self.hyphy_alignment + ".FUBAR.json"
        self.cache             = kwargs.get("cache", None)

        self.grid_size         = kwargs.get("grid_size", 20)
        self.nchains           = kwargs.get("nchains", 5)
        self.chain_length      = kwargs.get("chain_length", 2e6)
        self.burnin            = kwargs.get("burnin", 1e6)
        self.samples_per_chain = kwargs.get("samples_per_chain", 100)
        self.alpha             = kwargs.get("alpha", 0.5)
        
        self._sanity_fubar()
        self._build_full_command()
    
    
    
    def _sanity_fubar(self):
        """
            Sanity check FUBAR arguments
        """
        assert(self.grid_size >=5 and self.grid_size <=50), "\n[ERROR]: FUBAR grid size must be in range [5,50]."
        assert(self.nchains >=2 and self.nchains <=20), "\n[ERROR]: FUBAR nchains must be in range [2,20]."
        assert(self.chain_length >=5e5 and self.chain_length <=5e7), "\n[ERROR]: FUBAR chain length must be in range [5e5,5e7]."
        assert(self.burnin >=ceil(self.chain_length/20) and self.burnin <= ceil(95*self.chain_length/100)), "\n[ERROR]: FUBAR burnin size out of range."
        assert(self.samples_per_chain >=50 and self.samples_per_chain <= (self.chain_length - self.burnin)), "\n[ERROR]: FUBAR samples_per_chain out of range."
        assert(self.alpha >=0.001 and self.alpha <= 1), "\n[ERROR]: FUBAR Dirichlet prior parameter alpha must in be in range [0.001,1]."
        
        self.default_cache_path = self.default_json_path.replace(".json", ".cache")
        
        if self.cache is not None:
            if self.cache is False:
                self.cache_path = "/dev/null/"
            else:
                dirname = os.path.dirname(os.path.abspath(self.cache))
                assert( os.path.exists(dirname) ),"\n[ERROR]: Provided path to output cache does not exist."  
                self.cache_path = self.cache
        else:
            self.cache_path = self.default_cache_path
        
        
    def _build_analysis_command(self):
        """
            Construct the MEME command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree ,
                                           str(self.grid_size),
                                           str(self.nchains), 
                                           str(self.chain_length),
                                           str(self.burnin),
                                           str(self.samples_per_chain),
                                           str(self.alpha) ])


        
    def _save_output(self):
        """
            Move JSON  and cache to final location. 
        """        
        if self.user_json_path is not None:       
            self.final_path = self.user_json_path
        else:
            self.final_path = self.default_json_path   
        shutil.move(self.default_json_path, self.final_path)
        
        if self.cache_path != self.default_cache_path:
            shutil.move(self.default_cache_path, self.cache_path)
       

class MEME(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **alpha**, The p-value threshold for calling sites as positively or negatively selected. Note that this argument has 0 bearing on JSON output. Default: 0.1
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.


            **Examples:**
               
               >>> ### Define a default MEME analysis, where data is contained in a single file
               >>> mymeme = MEME(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default MEME analysis, where alignment and tree are in separate files 
               >>> mymeme = MEME(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")
               
               >>> ### Define a MEME analysis, specifying that selection be tested only on leaves
               >>> mymeme = MEME(data = "/path/to/data_with_tree.dat", branches = "Leaves" )

               >>> ### Define a MEME analysis, specifying a custom JSON output file
               >>> mymeme = MEME(data = "/path/to/data_with_tree.dat", output = "meme.json" )

               >>> ### Define MEME analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> mymeme = MEME(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)

               >>> ### Execute a defined MEME instance
               >>> mymeme.run_analysis()
        """                

        super(MEME, self).__init__(**kwargs)
        
        self.batchfile = "MEME.bf"
        self.default_json_path = self.hyphy_alignment + ".MEME.json"
        
        self.alpha = str( kwargs.get("alpha", 0.1) )
        
        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
        self._build_full_command()
        
    def _build_analysis_command(self):
        """
            Construct the MEME command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree ,
                                           self.branches , 
                                           self.alpha ])



class SLAC(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **alpha**, The p-value threshold for calling sites as positively or negatively selected. Note that this argument has 0 bearing on JSON output. Default: 0.1
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
                6. **bootstrap**, The number of samples used to assess ancestral reconstruction uncertainty, in [0,100000]. Default:100.
            
            **Examples:**
               
               >>> ### Define a default SLAC analysis, where data is contained in a single file
               >>> myslac = SLAC(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default SLAC analysis, where alignment and tree are in separate files 
               >>> myslac = SLAC(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")
               
               >>> ### Define a SLAC analysis, specifying that selection be tested only on leaves
               >>> myslac = SLAC(data = "/path/to/data_with_tree.dat", branches = "Leaves" )

               >>> ### Define a SLAC analysis, specifying 150 bootstrap replicates be used for ASR uncertainty
               >>> myslac = SLAC(data = "/path/to/data_with_tree.dat", bootstrap = 150 )
               
               >>> ### Define a SLAC analysis, specifying a custom JSON output file
               >>> myslac = SLAC(data = "/path/to/data_with_tree.dat", output = "slac.json" )

               >>> ### Define SLAC analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myslac = SLAC(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)

               >>> ### Execute a defined SLAC instance
               >>> myslac.run_analysis()
        """                

        super(SLAC, self).__init__(**kwargs)
        
        self.batchfile = "SLAC.bf"
        self.default_json_path = self.hyphy_alignment + ".SLAC.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
    
        self.alpha = str( kwargs.get("alpha", 0.1) )
    
        self.bootstrap_samples = kwargs.get("bootstrap", 100)
        self.range_bootstrap_samples = [0,100000]
        assert(self.bootstrap_samples >= self.range_bootstrap_samples[0] and self.bootstrap_samples <= self.range_bootstrap_samples[1]), "\n [ERROR] Number of samples to assess ASR uncertainty must be in range [0,100000]."
        self.bootstrap_samples = str(self.bootstrap_samples)
        self._build_full_command()
        
    def _build_analysis_command(self):
        """
            Construct the SLAC command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree ,
                                           self.branches , 
                                           self.bootstrap_samples, 
                                           self.alpha ])





class ABSREL(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.

            **Examples:**
               
               >>> ### Define a default ABSREL analysis, where data is contained in a single file
               >>> myabsrel = ABSREL(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default ABSREL analysis, where alignment and tree are in separate files 
               >>> myabsrel = ABSREL(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")
               
               >>> ### Define a ABSREL analysis, specifying that selection be tested only on leaves
               >>> myabsrel = ABSREL(data = "/path/to/data_with_tree.dat", branches = "Leaves" )
 
               >>> ### Define ABSREL analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myabsrel = ABSREL(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)
                             
               >>> ### Execute a defined ABSREL instance
               >>> myabsrel.run_analysis()

        """                
                
        super(ABSREL, self).__init__(**kwargs)
        
        self.batchfile = "aBSREL.bf"
        self.default_json_path = self.hyphy_alignment + ".ABSREL.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
        self._build_full_command()

    def _build_analysis_command(self):
        """
            Construct the aBSREL command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree,
                                           self.branches
                                         ])



class BUSTED(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.

            **Examples:**
               
               >>> ### Define a default BUSTED analysis, where data is contained in a single file
               >>> mybusted = BUSTED(data = "/path/to/data_with_tree.dat")

               >>> ### Define a default BUSTED analysis, where alignment and tree are in separate files 
               >>> mybusted = BUSTED(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre")
               
               >>> ### Define a BUSTED analysis, specifying that selection be tested only on leaves
               >>> mybusted = BUSTED(data = "/path/to/data_with_tree.dat", branches = "Leaves" )
               
               >>> ### Define BUSTED analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> mybusted = BUSTED(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)
               
               >>> ### Execute a defined BUSTED instance
               >>> mybusted.run_analysis()
        """                
                
        super(BUSTED, self).__init__(**kwargs)
        
        self.batchfile = "BUSTED.bf"
        self.default_json_path = self.hyphy_alignment + ".BUSTED.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection() 
        self._build_full_command()

    def _build_analysis_command(self):
        """
            Construct the BUSTED command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.genetic_code ,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree,
                                           self.branches
                                         ])
       



class RELAX(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **test_label**, The label (must be found in your tree) corresponding to the **test** branch set
                3. **reference_label**, The label f(must be found in your tree) corresponding to the **reference** branch set. **Only provide this argument if your tree has multiple labels in it.**
                4. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                5. **analysis_type**, "All" (run hypothesis test and fit descriptive models) or "Minimal" (only run hypothesis test). Default: "All".
                6. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.


            **Examples:**
               
               >>> ### Define a default RELAX analysis, where data is contained in a single file and test branches are labeled "test"
               >>> myrelax = RELAX(data = "/path/to/data_with_tree.dat", test_label = "test")

               >>> ### Define a default RELAX analysis, where alignment and tree are in separate files and test branches are labeled "test"
               >>> myrelax = RELAX(alignment = "/path/to/alignment.fasta", tree = "/path/to/tree.tre", test_label = "test")
               
               >>> ### Define a default RELAX analysis, where data is contained in a single file, test branches are labeled "test", and reference branches are labeled "ref"
               >>> myrelax = RELAX(data = "/path/to/data_with_tree.dat", test_label = "test", reference_label = "ref")

               >>> ### Define a default RELAX analysis, where data is contained in a single file, test branches are labeled "test", and the Minimal analysis is run
               >>> myrelax = RELAX(data = "/path/to/data_with_tree.dat", test_label = "test", analysis_type = "Minimal")
 
               >>> ### Define RELAX analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myrelax = RELAX(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)
              
               >>> ### Execute a defined RELAX instance
               >>> myrelax.run_analysis()
        """                
                
        super(RELAX, self).__init__(**kwargs)
        
        self.batchfile = "RELAX.bf"
        self.default_json_path = self.hyphy_alignment + ".RELAX.json"

        self._find_all_labels()
        if len(self._all_labels) == 0:
            raise AssertionError("\n[ERROR] RELAX requires at least one label in the tree. Visit http://veg.github.io/phylotree.js/ for assistance labeling your tree.")
        
        self.test_label = kwargs.get("test_label", None)
        assert(self.test_label in self._all_labels), "\n [ERROR] You must provide a `test_label` arguement that corresponds to a label in your **labeled tree**. Visit http://veg.github.io/phylotree.js/ for assistance labeling your tree."
        
        self.reference_label = kwargs.get("reference_label", None)
        if len(self._all_labels) > 1:
            if self.reference_label is None:
                print("\nWARNING: No branches were selected as 'reference' even though multiple labels exist in the tree. Defaulting to using all non-test branches.")
                self.reference_label = self.shared_branch_choices[-1]
            else:
                assert(self.reference_label in self._all_labels), "\n [ERROR] The value for `reference_label` must correspond to a label in your tree. To simply use all non-test branches as reference, do not provide the argument `reference_label`."
        assert(self.test_label != self.reference_label), "\n[ERROR] Must use different test and reference labels."

        self.allowed_types = ("All", "Minimal")
        self.analysis_type = kwargs.get("analysis_type", self.allowed_types[0]).capitalize()
        assert(self.analysis_type in self.allowed_types), "\n[ERROR] Incorrect analysis type specified. Provide either `All` or `Minimal`."
        self._build_full_command()

    def _build_analysis_command(self):
        """
            Construct the BUSTED command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        if self.reference_label is None:
            self.analysis_command = " ".join([  self.batchfile_with_path , 
                                                self.genetic_code ,
                                                self.hyphy_alignment ,
                                                self.hyphy_tree,
                                                self.test_label, 
                                                self.analysis_type
                                             ])
        else:
            self.analysis_command = " ".join([  self.batchfile_with_path , 
                                                self.genetic_code ,
                                                self.hyphy_alignment ,
                                                self.hyphy_tree,
                                                self.test_label, 
                                                self.reference_label,
                                                self.analysis_type
                                             ])


class LEISR(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)
                2. **type**, either "nucleotide" or "protein" indicating the type of data being analyzed

            Optional keyword arguments:
                1. **hyphy**, a :code:`HyPhy()` instance. Default: Assumes canonical HyPhy install.
                2. **model**, The model to use to fit relative rates, i.e. GTR for nucleotides or LG for amino acids. For full options, please see HyPhy. Default: JC69.
                3. **rate_variation**, Whether to apply rate variation to branch length optimization. Options include No, Gamma, GDD. Note that Gamma and GDD will use four categories each. Default: No
            

            **Examples:**
               
               >>> ### Define a LEISR Protein analysis, where data is contained in a single file and the WAG model with no rate variation is used
               >>> myleisr = LEISR(data = "/path/to/data_with_tree.dat", type = "protein", model = "WAG")

               >>> ### Define a LEISR Protein analysis, where data is contained in a single file and the WAG model with Gamma rate variation is used
               >>> myleisr = LEISR(data = "/path/to/data_with_tree.dat", type = "protein", model = "WAG", rate_variation = "Gamma")

               >>> ### Define a LEISR Nucleotide analysis, where data is contained in a single file and the HKY85 model with GDD rate variation is used
               >>> myleisr = LEISR(data = "/path/to/data_with_tree.dat", type = "nucleotide", model = "HKY85", rate_variation = "GDD")

               >>> ### Define LEISR analysis with a custom Hyphy, which is also defined here:
               >>> myhyphy = HyPhy(suppress_log = True, quiet = True) ## HyPhy will use default canonical install but run in full quiet mode
               >>> myleisr = LEISR(data = "/path/to/data_with_tree.dat", hyphy=myhyphy)

               >>> ### Execute a defined LEISR instance
               >>> myleisr.run_analysis()
         """                

        super(LEISR, self).__init__(**kwargs)
        
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/"
        self.batchfile = "LEISR.bf"
                
        self.default_json_path_choices = self.hyphy_alignment + ".LEISR.json"
        self.type_nucleotide = "Nucleotide"
        self.type_protein    = "Protein"
        
        self.type = kwargs.get("type", None)
        assert(self.type is not None),"\n[ERROR]: Must specify either 'nucleotide' or 'protein' for keyword argument `type` (case insensitive)."
        self.type = self.type.capitalize()
        assert(self.type == self.type_nucleotide or self.type == self.type_protein), "\n[ERROR]: Must specify either 'nucleotide' or 'protein' for keyword argument `type` (case insensitive)."
        
        self.rv   = kwargs.get("rate_variation", "No").capitalize()
        assert(self.rv in ["No", "Gamma", "Gdd"]), "\n[ERROR] Provided rate variation is unavailable. Use either `No`, `Gamma`, `GDD`."
        if self.rv == "Gdd":
            self.rv = "GDD"
            
        if self.type == self.type_nucleotide:
            self.model = kwargs.get("model", "GTR")
            assert(self.model in self.available_nucleotide_models), "\n[ERROR] Provided nucleotide model is unavailable."
            
        elif self.type == self.type_protein:
            self.model = kwargs.get("model", "JC69")
            assert(self.model in self.available_protein_models), "\n [ERROR] Provided protein model is unavailable."
        
        self._build_full_command()    
         
    def _build_analysis_command(self):
        """
            Construct the LEISR command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.hyphy_alignment ,
                                           self.hyphy_tree,
                                           self.type,
                                           self.model,
                                           self.rv
                                         ])

    


        
        
                