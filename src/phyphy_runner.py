#!/usr/bin/env python

##############################################################################
##  phyhy: *P*ython *HyPhy*: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################

	
	
"""
    Execute a standard HyPhy analysis.
"""
    
import subprocess
import os
import shutil
import re
from Bio import Phylo
from copy import deepcopy
from StringIO import StringIO
from math import ceil
import json

_DEFAULT_PATH = "/usr/local/lib/hyphy/"
_GENETIC_CODE = {
                  1: "Universal",
                  2: "Vertebrate mtDNA",
                  3: "Yeast mtDNA",
                  4: "Mold/Protozoan mtDNA",
                  5: "Invertebrate mtDNA",
                  6: "Ciliate Nuclear",
                  7: "Echinoderm mtDNA",
                  8: "Euplotid Nuclear",
                  9: "Alt. Yeast Nuclear",
                  10: "Ascidian mtDNA",
                  11: "Flatworm mtDNA",
                  12: "Blepharisma Nuclear",
                  13: "Chlorophycean mtDNA",
                  14: "Trematode mtDNA",
                  15: "Scenedesmus obliquus mtDNA",
                  16: "Thraustochytrium mtDNA",
                  17: "Pterobranchia mtDNA",
                  18: "SR1 and Gracilibacteria",
                  19: "Pachysolen Nuclear"
                }


class HyPhy():
    """
        This class creates a HyPhy instance. Generally this is only a necessary step if any of these applies:
            + You wish to use a local **build** of HyPhy (not a canonically installed build)
            + You wish to use a local **install** of HyPhy (installed elsewhere from /usr/local)
            + You wish to use a different HyPhy executable from the default, HYPHYMP
    
        Optional keyword arguments to __init__:
            1. **executable**, the desired executable to use (ie HYPHYMPI). Default: HYPHYMP
            2. **build_path**, the path to a **local hyphy build**. Use this argument if you have compiled hyphy in the downloaded hyphy/ directory and **did not run make install**
            3. **install_path**, the path to a **hyphy install**. Use this argument if you have specified a different installation path for hyphy, i.e. you provided `-DINSTALL_PREFIX=/other/path/` to cmake.
            4. **cpu**, the maximum number of CPUs per analysis. By default, HyPhy will take as many CPUs as it can/requires. This argument will limit the maximum.
            5. **quiet**, suppress screen output (Note, HyPhy will still creates messages.log and errors.log files, when applicable). Default: False
            6. **suppress_log**, suppress messages.log and errors.log files. Default: False. (If True, redirects to /dev/null/)
    """

    def __init__(self, **kwargs):


        self.executable    = kwargs.get("executable", "HYPHYMP")
        self.build_path    = kwargs.get("build_path", None)  
        self.install_path  = kwargs.get("install_path", None) 
        self.cpu           = kwargs.get("cpu", None)       
        self.quiet         = kwargs.get("quiet", False) ### run hyphy quietly
        self.suppress_log   = kwargs.get("suppres_stdout", False) ### send messages.log, errors.log to /dev/null
        
        
        
        ### Sanity checks for a local install ###
        if self.build_path is not None: 
            assert(os.path.exists(self.build_path)), "\n[ERROR] Build path does not exist."
            self.build_path = os.path.abspath(self.build_path) + "/" ## os.path.abspath will strip any trailing "/"
            self.libpath = self.build_path + "res/"
            assert(os.path.exists(self.libpath)), "\n[ERROR]: Build path does not contain a correctly built HyPhy."
            self.executable = self.build_path + self.executable
            self.hyphy_call = self.executable + " LIBPATH=" + self.libpath
        
        else: 
            ## Installed in non-default path
            if self.install_path is not None:
                assert(os.path.exists(self.install_path)), "\n[ERROR]: Install path does not exist."
                self.libpath = self.install_path + "lib/hyphy/"
                assert(os.path.exists(self.libpath)), "\n[ERROR]: Install path does not contain a correctly built HyPhy."               
                self.executable = self.build_path + self.executable
                self.hyphy_call = self.executable + " LIBPATH=" + self.libpath
            ## Installed in default path
            else:
                self.libpath = _DEFAULT_PATH
                self.hyphy_call = deepcopy(self.executable)
            
        ## Ensure executable exists somewhere
        with open("/dev/null", "w") as hushpuppies:
            exit_code = subprocess.call(["which", self.executable], stdout = hushpuppies, stderr = hushpuppies) # If you're reading this, I hope you enjoy reading hushpuppies as much as I enjoyed writing it. --SJS
            if exit_code == 1:
                raise AssertionError("\n[ERROR]: HyPhy executable not found. Please ensure it is properly installed or in your provided path.")

        if self.cpu is not None:
            self.hyphy_call += " CPU=" + str(self.cpu)
        
        if self.suppress_log is True:
            self.hyphy_call += "USEPATH=/dev/null/"

        
        
        

class Analysis(object):
    """
        Parent class for all analysis methods. 
        Child classes:
            + ABSREL
            + BUSTED
            + FEL
            + FUBAR
            + LEISR
            + MEME
            + RELAX
            + SLAC                        
    """
        
            
    def __init__(self, **kwargs):
        """
            Initialize a HyPhy analysis. 
            
            Required arguments:
                1. **alignment** _and_ **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)
                           
            Optional keyword arguments:
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **output**, name (and path to) to final output JSON file. Default: Goes to same directory as provided data
            
            See children classes for analysis-specific arguments.
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
                    
        ### Unused in AA analyses 
        self.genetic_code = kwargs.get("genetic_code", "Universal")
        assert(self.genetic_code in list(_GENETIC_CODE.values()) or self.genetic_code in list(_GENETIC_CODE.keys())), "\n[ERROR] Incorrect genetic code specified."
        if self.genetic_code in list(_GENETIC_CODE.keys()):
            for k,v in _GENETIC_CODE.items():
                if v == self.genetic_code:
                    self.genetic_code = str(k)
                    break
        self.genetic_code = self.genetic_code.replace(" ", "\ ") # Sigh.

        
        
        ### Will be overriden for SelectionAnalysis methods
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/SelectionAnalyses/"

        self.shared_branch_choices = ("All", "Internal", "Leaves", "Unlabeled branches")
        
        self.available_protein_models = ("JC69", "WAG", "LG", "JTT")
        self.available_nucleotide_models = ("GTR", "HKY85", "JC69")
    
    
    
    
    def _format_yesno(self, argument):
        """
            Format argument to be Yes/No from True/False
        """
        self.yesno_truefalse = {True: "Yes", False: "No"}
        if type(argument) == str:
            argument.capitalize()
            assert(argument in list(self.yesno_truefalse.keys())),"\n[ERROR]: Incorrect Yes/No argument."
        elif type(argument) is bool:
            argument = self.yesno_truefalse[argument]
        else:
            raise TypeError("\n[ERROR]: Incorrect Yes/No argument.")
        return argument            
    
    
    def _check_files(self):
        """
            Check provided paths for alignment+tree or data. Assign input hyphy variables accordingly.
            Additionally extract the tree string
        """
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
                        raise AssertionError("\n[ERROR] Malformed tree in input data.")


            
    def _build_command(self):
        print("Parent method. Not run.")

   
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
        
            
    def run_analysis(self):
        """
            Call HyPhy as a subprocess to run a given analysis. 
        """    
        self._build_analysis_command()
        full_command = " ".join([self.hyphy.hyphy_call, self.analysis_command])

        if self.hyphy.quiet:
            with open("/dev/null", "w") as quiet:
                check = subprocess.call(full_command, shell = True, stdout = quiet, stderr = quiet)
        else:    
            check = subprocess.call(full_command, shell = True)
        assert(check == 0), "\n[ERROR] HyPhy failed to run."

        self._save_output()


    def return_json(self):
        """
            Return **parsed*** JSON to user
            
            
        """


    def _save_output(self):
        """
            Move JSON to final location. 
        """        
        if self.user_json_path is None:
            final_path = self.default_json_path
        else:
            final_path = self.user_json_path            
        shutil.move(self.default_json_path, final_path)

 



class FEL(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **srv**, Employ synonymous rate variation in inference (i.e. allow dS to vary across sites?). Values "Yes"/"No" or True/False accepted. Default: True.
                3. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                4. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                5. **alpha**, The p-value threshold for calling sites as positively or negatively selected. Default: 0.1
                6. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                
        
        super(FEL, self).__init__(**kwargs)
        
        self.batchfile = "FEL.bf"
        self.default_json_path = self.hyphy_alignment + ".FEL.json"

        self.srv   = kwargs.get("srv", "Yes") ## They can provide T/F or Yes/No
        self.srv   = self._format_yesno(self.srv)
        self.alpha = str( kwargs.get("alpha", 0.1) )

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
        

        
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
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data. 
                3. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
                4. **grid_size**, Number of grid points per rate grid dimension (Default: 20, allowed [5,50])
                5. **nchains**, Number of MCMC chains to run (Default: 5, allowed [2,20])
                6. **chain_length**, The length of each chain (Default: 2e6, allowed [5e5,5e7])
                7. **burnin**, Number of samples to consider as burn-in (Default 1e6, allowed [ceil(chain_length/20),ceil(95*chain_length/100)])
                8. **samples_per_chain**, Number of samples to draw per chain (Default 100, allowed [50,chain_length-burnin])
                9. **alpha**, The concentration parameter of the Dirichlet prior (Default 0.5, allowed[0.001,1])
                10. **cache**, Name (and path to) output FUBAR cache. Default: goes to same directory as provided data. Provide the argument **False** to not save the cache (this argument simply sends it to /dev/null)
        """                



        super(FUBAR, self).__init__(**kwargs)
        
        self.batchfile = "FUBAR.bf"
        self.default_json_path = self.hyphy_alignment + ".MEME.json"
        self.cache             = kwargs.get("cache", None)

        self.grid_size         = kwargs.get("grid_size", 20)
        self.nchains           = kwargs.get("nchains", 5)
        self.chain_length      = kwargs.get("chain_length", 2e6)
        self.burnin            = kwargs.get("burnin", 1e6)
        self.samples_per_chain = kwargs.get("grid_size", 100)
        self.alpha             = kwargs.get("alpha", 0.5)
        
        self._sanity_fubar()
    
    
    
    
    def _sanity_fubar(self):
        """
            Sanity check FUBAR arguments
        """
        assert(self.grid_size >=5 and self.grid_size <=50), "\n[ERROR]: FUBAR grid size must be in range [5,50]."
        assert(self.nchains >=2 and self.nchains <=20), "\n[ERROR]: FUBAR nchains must be in range [2,20]."
        assert(self.chain_length >=5e5 and self.chain_length <=5e7), "\n[ERROR]: FUBAR chain length must be in range [5e5,5e7]."
        assert(self.burnin >=ceil(self.chain_length/20) and self.burnin <= ceil(95*self.chain_length/100)), "\n[ERROR]: FUBAR burnin size out of range."
        assert(self.samples_per_chain >=50 and self.samples_per_chain <= (self.chain_lenghth - self.burnin)), "\n[ERROR]: FUBAR samples_per_chain out of range."
        assert(self.alpha >=0.001 and self.alpha <= 1), "\n[ERROR]: FUBAR Dirichlet prior parameter alpha must in be in range [0.001,1]."
        
        self.default_cache_path = self.user_json_path.replace("json", "cache")
        if self.cache is False:
            self.cache_path = "/dev/null/"
        elif type(self.cache) == "str":
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
                                           self.grid_size,
                                           self.nchains, 
                                           self.chain_length,
                                           self.burnin,
                                           self.samples_per_chain,
                                           self.alpha ])


    def _save_output(self):
        """
            Move JSON  and cache to final location. 
        """        

        if self.user_json_path is None:
            final_json_path = self.default_json_path
        else:
            final_json_path = self.user_json_path            
        shutil.move(self.default_json_path, final_path)
        
        if self.cache_path != self.default_cache_path:
            shutil.move(self.default_cache_path, self.cache_path)
       

class MEME(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **alpha**, The p-value threshold for calling sites as positively selected. Default: 0.1
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                

        super(MEME, self).__init__(**kwargs)
        
        self.batchfile = "MEME.bf"
        self.default_json_path = self.hyphy_alignment + ".MEME.json"
        
        self.alpha = str( kwargs.get("alpha", 0.1) )
        
        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
    
        
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
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **alpha**, The p-value threshold for calling sites as positively selected. Default: 0.1
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
                6. **bootstrap_samples**, The number of samples used to assess ancestral reconstruction uncertainty, in [0,100000]. Default:100.
        """                

        super(SLAC, self).__init__(**kwargs)
        
        self.batchfile = "SLAC.bf"
        self.default_json_path = self.hyphy_alignment + ".SLAC.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
    
        self.alpha = str( kwargs.get("alpha", 0.1) )
    
        self.bootstrap_samples = kwargs.get("bootstrap_samples", 100)
        self.range_bootstrap_samples = [0,100000]
        assert(self.bootstrap_samples >= self.range_bootstrap_samples[0] and self.bootstrap_samples <= self.range_bootstrap_samples[1]), "\n [ERROR] Number of samples to assess ASR uncertainty must be in range [0,100000]."
        self.bootstrap_samples = str(self.bootstrap_samples)
        
        
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
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                
                
        super(ABSREL, self).__init__(**kwargs)
        
        self.batchfile = "aBSREL.bf"
        self.default_json_path = self.hyphy_alignment + ".json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()

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
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **branches**, Branches to consider in site-level selection inference. Values "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                4. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                
                
        super(BUSTED, self).__init__(**kwargs)
        
        self.batchfile = "BUSTED.bf"
        self.default_json_path = self.hyphy_alignment + ".BUSTED.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection() 


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
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **test_label**, The label (must be found in your tree) corresponding to the **test** branch set
                3. **reference_label**, The label f(must be found in your tree) corresponding to the **reference** branch set. **Only provide this argument if your tree has multiple labels in it.**
                4. **output**, Name (and path to) to final output JSON file. Default: Goes to same directory as provided data
                5. **analysis_type**, "All" (run hypothesis test and fit descriptive models) or "Minimal" (only run hypothesis test). Default: "All".
                6. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
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

        self.allowed_types = ("All", "Minimal")
        self.analysis_type = kwargs.get("analysis_type", self.allowed_types[0]).capitalize()
        assert(self.analysis_type in self.allowed_types), "\n[ERROR] Incorrect analysis type specified. Provide either `All` or `Minimal`."


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
                1. **hyphy**, a HyPhy() instance. Default: Assumes canonical HyPhy install.
                2. **model**, The nucleotide model to use to fit relative rates. Options include GTR, HKY85, or JC69 for Nucleotide (Default GTR), and LG, WAG, JTT, and JC69 for Protein (Default JC69).
                3. **rate_variation**, Whether to apply rate variation to branch length optimization. Options include No, Gamma, GDD (Default No). Note that Gamma and GDD will use four categories each.
                4. **plusF**, Only applicable to protein analyses, this option controls the protein model should use +F frequencies. +F means frequencies will be empirically read in from the provided data, in contrast to using the default model frequencies. Default: True.
            
         """                
        super(LEISR, self).__init__(**kwargs)
        
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/"
        self.batchfile = "LEISR.bf"
        self.default_json_path = self.hyphy_alignment + ".site-rates.json"
        self.type_nucleotide = "Nucleotide"
        self.type_protein    = "Protein"
        
        self.type = kwargs.get("type", None).capitalize()
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
            self.plus_f = kwargs.get("plusF", True)
            self.plus_f = self._format_yesno(self.plus_f)
            self.rv = self.rv + " " + self.plus_f ## To provide analysis command
        else:
            raise AssertionError("\n[ERROR]: Must specify either 'nucleotide' or 'protein' for keyword argument `type` (case insensitive).")
          
         
    def _build_analysis_command(self):
        """
            Construct the LEISR command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.type,
                                           self.model,
                                           self.rv,
                                           self.hyphy_alignment ,
                                           self.hyphy_tree
                                         ])



    

        
        
                