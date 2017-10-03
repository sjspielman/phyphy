"""
    Run a HyPhy analysis
"""
    
import subprocess
import os
import shutil
import re
from Bio import Phylo
from copy import deepcopy
from StringIO import StringIO

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

    def __init__(self, **kwargs):
        """
            Initialize a HyPhy instance. Generally this is only a necessary step if either applies:
                + You wish to use a local build of HyPhy (not an installed build)
                + You wish to use a different HyPhy executable (default is HYPHYMP)
        
            Optional keyword arguments:
                1) executable, the desired executable to use (ie HYPHYMPI). Default: HYPHYMP
                2) path, the path to a **local hyphy build**. Only use this argument if you **do not** want to use the installed hyphy in /usr/local.
                3) cpu, the maximum number of CPUs per analysis. 
                4) quiet, suppress screen output (still creates messages.log and errors.log, when applicable). Default: False
        """

        self.executable = kwargs.get("executable", "HYPHYMP")
        self.user_path  = kwargs.get("path", None)  ### only for local or non-standard installs
        self.cpu        = kwargs.get("cpu", None)       
        self.quiet      = kwargs.get("quiet", False) ### run hyphy quietly
        
        
        ### Sanity checks for a local/non-standard install ###
        if self.user_path is not None: 
            assert(os.path.exists(self.user_path)), "\n[ERROR] HyPhy not detected in provided path."
            self.user_path = os.path.abspath(self.user_path) + "/" ## os.path.abspath will strip any trailing "/"
            self.libpath = self.user_path + "res/"
            assert(os.path.exists(self.libpath)), "\n[ERROR]: User path does not contain a correctly built HyPhy."
            self.executable = self.user_path + self.executable
            self.hyphy_call = self.executable + " LIBPATH=" + self.libpath
        
        else: ## default install
            self.libpath = _DEFAULT_PATH
            self.hyphy_call = deepcopy(self.executable)
            
        ## Ensure executable exists somewhere
        with open("/dev/null", "w") as hushpuppies:
            exit_code = subprocess.call(["which", self.executable], stdout = hushpuppies, stderr = hushpuppies) # If you're reading this, I hope you enjoy reading hushpuppies as much as I enjoyed writing it. --SJS
            if exit_code == 1:
                raise AssertionError("\n[ERROR]: HyPhy executable not found. Please ensure it is properly installed or in your provided local install path.")

        if self.cpu is not None:
            self.hyphy_call += " CPU=" + str(self.cpu)

        
        
        

class Analysis(object):
    
    def __init__(self, **kwargs):
        """
            Parent class for all analysis methods.
        """

        self.hyphy = kwargs.get("hyphy", None)
        if self.hyphy is None:
            self.hyphy = HyPhy()
        
        
        self.alignment = kwargs.get("alignment", None) ### alignment only
        self.tree      = kwargs.get("tree", None)    ### tree only
        self.data      = kwargs.get("data", None)    ### combined alignment and tree or NEXUS
        self._check_files()

        self.alpha          = str( kwargs.get("alpha", 0.1) )   ### significance
        self.user_json_path = kwargs.get("output", None)
        
        
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
            Upon completion, move JSON to the user-specified location (if applicable).
        """    
        self._build_analysis_command()
        full_command = " ".join([self.hyphy.hyphy_call, self.analysis_command])
        

        if self.hyphy.quiet:
            with open("/dev/null", "w") as quiet:
                check = subprocess.call(full_command, shell = True, stdout = quiet, stderr = quiet)
        else:    
            check = subprocess.call(full_command, shell = True)
        assert(check == 0), "[ERROR] HyPhy failed to run."
        
        ### Move JSON to final resting place
        if self.user_json_path is not None:
            self.json_path = self.default_json_path
        else:
            self.json_path = self.user_json_path
            shutil.move(self.default_json_path, self.user_json_path)



 



class FEL(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, your HyPhy instance
                2. **two_rate**, "Yes"/"No" or True/False accepted
                3. **branches**, "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                4. **output**, New path and/or file for the JSON
                5. **alpha**, The p-value threshold for selection
                6. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                
                
        super(FEL, self).__init__(**kwargs)
        
        self.batchfile = "FEL.bf"
        self.default_json_path = self.hyphy_alignment + ".FEL.json"

        self.tworate = kwargs.get("two_rate", "Yes") ## They can provide T/F or Yes/No
        self.tworate = self._format_yesno(self.tworate)

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
                                           self.tworate , 
                                           self.alpha ])

        
        
        
        
        
        
       

class MEME(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, your HyPhy instance
                2. **branches**, "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, New path and/or file for the JSON
                4. **alpha**, The p-value threshold for selection
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                

        super(MEME, self).__init__(**kwargs)
        
        self.batchfile = "MEME.bf"
        self.default_json_path = self.hyphy_alignment + ".MEME.json"

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
                1. **hyphy**, your HyPhy instance
                2. **branches**, "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, New path and/or file for the JSON
                4. **bootstrap_samples**, The number of samples used to assess ancestral reconstruction uncertainty, in [0,100000]. Default:100.
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
        """                

        super(SLAC, self).__init__(**kwargs)
        
        self.batchfile = "SLAC.bf"
        self.default_json_path = self.hyphy_alignment + ".SLAC.json"

        self.branches = kwargs.get("branches", "All")
        self._sanity_branch_selection()
    
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
                1. **hyphy**, your HyPhy instance
                2. **branches**, "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, New path and/or file for the JSON
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
                1. **hyphy**, your HyPhy instance
                2. **branches**, "All", "Internal", "Leaves", "Unlabeled branches", or a **specific label** are accepted
                3. **output**, New path and/or file for the JSON
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
                1. **hyphy**, your HyPhy instance
                2. **test_label**, the label found in your tree(!) referring to the branches along which to test for relaxation of selection
                3. **reference_label**, the label found in your tree(!) referring to the reference branches used in the test for relaxation of selection. **Only provide this argument if your tree has multiple labels in it!**
                4. **analysis_type**, "All" (run hypothesis test and fit descriptive models) or "Minimal" (only run hypothesis test). Default: "All".
                5. **genetic_code**, the genetic code to use in codon analysis, Default: Universal. Consult NIH for details.
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
                print("WARNING: No branches were selected as 'reference' even though multiple labels exist in the tree. Defaulting to using all non-test branches.")
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


class RelativeProteinRates(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, your HyPhy instance
                2. **model**, the protein model to use to fit relative rates. Default: JC.
                3. **plusF**, should the model use +F frequencies? +F means frequencies will be empirically read in from alignment, in contrast to using the default model frequencies. Default: False.
         """                
        super(RelativeProteinRates, self).__init__(**kwargs)
        
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/ProteinAnalyses/"
        self.batchfile = "relative_prot_rates.bf"
        self.default_json_path = self.hyphy_alignment + ".site-rates.json"
        
        self.model = kwargs.get("model", "JC69")
        assert(self.model in self.available_protein_models), "\n [ERROR] Provided protein model is unavailable."
        
        self.plus_f = kwargs.get("plusF", "False")
        self.plus_f = self._format_yesno(self.plus_f)


    def _build_analysis_command(self):
        """
            Construct the relative_prot_rates command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.hyphy_alignment ,
                                           self.hyphy_tree,
                                           self.model,
                                           self.plus_f
                                         ])
   


class RelativeNucleotideRates(Analysis):

    def __init__(self, **kwargs):
        """
            Required arguments:
                1. **alignment** and **tree** OR **data**, either a file for alignment and tree separately, OR a file with both (combo FASTA/newick or nexus)

            Optional keyword arguments:
                1. **hyphy**, your HyPhy instance
                2. **model**, the nucleotide model to use to fit relative rates. Default: JC69.
         """                
        super(RelativeNucleotideRates, self).__init__(**kwargs)
        
        self.analysis_path = self.hyphy.libpath + "TemplateBatchFiles/"
        self.batchfile = "relative_nucleotide_rates.bf"
        self.default_json_path = self.hyphy_alignment + ".site-rates.json"
        
        self.model = kwargs.get("model", "JC69")
        assert(self.model in self.available_nucleotide_models), "\n[ERROR] Provided nucleotide model is unavailable."


    def _build_analysis_command(self):
        """
            Construct the relative_prot_rates command with all arguments to provide to the executable. 
        """
        self.batchfile_with_path = self.analysis_path + self.batchfile
        
        self.analysis_command = " ".join([ self.batchfile_with_path , 
                                           self.hyphy_alignment ,
                                           self.hyphy_tree,
                                           self.model,
                                         ])
  

    

def main():
    
    if __name__ == "__main__":
    
        ## Check out these sweet relative paths!!!! 
        ### when providing the data, give either alignment and tree OR data. 
        #aa_alignment = "data/aa.fasta"  ### file with AA alignment
        #codon_alignment = "data/seqs.fasta" ### file with codon alignment
        #tree      = "data/test.tre"  ### file with just tree
        #data      = "data/seqs.dat"    ### file with codon sequences *and* tree
    
        ### output file, hyphyhelper will move the output json to here for you
        #json = "out.json"
    
        ## FIRST, Create a HyPhy instance if you want to use a local (aka not installed into /usr/local) hyphy and/or specify other things. See __init__ docstring for the things.
        #hyphy = HyPhy(path = "/path/to/local/hyphy", quiet=False)
    
    
        #f = FEL(hyphy = hyphy, data = "original_part.nex", two_rate = False, output = json)     ## NOTE: This line could be used instead:   f = FEL(hyphy = hyphy, data = data, two_rate = False, output = json)
        #f.run_analysis()    

        ### FEL ###
        #f = FEL(hyphy = hyphy, alignment = codon_alignment, tree = tree, two_rate = False, output = json)     ## NOTE: This line could be used instead:   f = FEL(hyphy = hyphy, data = data, two_rate = False, output = json)
        #f.run_analysis()

        ### MEME ###
        #f = MEME(hyphy = hyphy, alignment = codon_alignment, tree = tree, output = json)
        #f.run_analysis()

        ### SLAC ###   ###### NOTE: CURRENT V2.3-DEV SLAC IMPLEMENTATION IS BROKEN. THIS WILL NOT RUN.
        #f = SLAC(hyphy = hyphy, alignment = codon_alignment, tree = tree, output = json)     
        #f.run_analysis()

        #### ABSREL ####
        #f = ABSREL(hyphy = hyphy, data = data, branches = "Internal", output = json)
        #f.run_analysis()
            
        #### RELAX ####
        #f = RELAX(hyphy = hyphy, data = data, test_label = "test", reference_label = "reference", analysis_type = "Minimal", output = json)
        #f.run_analysis()
 
        #### BUSTED ####
        #f = BUSTED(hyphy = hyphy, data = data, output = json, foreground = ["Node1", "Node2", "t3"])
        #f.run_analysis()        

        ### RelativeProteinRates ###
        #f = RelativeProteinRates(hyphy = hyphy, alignment = aa_alignment, tree = tree, model = "JC", plusF = True)
        #f.run_analysis()
    
        ### RelativeNucleotideRates ###
        f = RelativeNucleotideRates(hyphy = hyphy, alignment = codon_alignment, tree = tree, model = "GTR")
        f.run_analysis()        
main()
        
    
    

        
        
                