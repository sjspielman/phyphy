# phyphy [![DOI](http://joss.theoj.org/papers/10.21105/joss.00514/status.svg)](https://doi.org/10.21105/joss.00514)

**The current release is version 0.4.2**.

[`phyphy`](http://sjspielman.org/phyphy) aims to facilitate [HyPhy](http://hyphy.org) usage in two primary ways:
1) Execute standard analyses in a Python scripting environment
2) Conveniently parse various information from the resulting JSON output.

This functionality makes batch usage eminently more convenient. Never use the HyPhy prompt again!

`phyphy` is pronounced "feye-feye" and is so-named for "**P**ython **Hyphy**". Importantly, this name was chosen for maximal pronunciation enjoyment. `phyphy` is compatible with either Python2.7 or Python3.

`phyphy` is compatible with [**HyPhy version >=2.3.7**](https://github.com/veg/hyphy/releases) and supports the following analyses (in alphabetical order):

+ [`FEL` (Fixed Effects Likelihood)](http://hyphy.org/methods/selection-methods/#fel): Infer pervasive selection at individual sites using maximum likelihood
+ [`MEME` (Mixed Effects Model of Evolution)](http://hyphy.org/methods/selection-methods/#meme): Infer episodic selection at individual sites
+ [`SLAC` (Single-Likelihood Ancestor Counting)](http://hyphy.org/methods/selection-methods/#slac): Infer pervasive selection at individual sites using a counting-based approach
+ [`FUBAR` (Fast, Unconstrained Bayesian AppRoximation)](http://hyphy.org/methods/selection-methods/#fubar): Infer pervasive selection at individual sites using a Bayesian approach
+ [`aBSREL` (adaptive Branch-Site Random Effects Likelihood)](http://hyphy.org/methods/selection-methods/#absrel): Branch-site model of lineage-specific selection
+ [`BUSTED` (Branch-Site Unrestricted Statistical Test for Episodic Diversification)](http://hyphy.org/methods/selection-methods/#busted): Branch-site model of whole-gene selection
+ [`RELAX` (Test for selection RELAXation)](http://hyphy.org/methods/selection-methods/#relax): Test for relaxed or intensified selection on a specified branch set
+ [`LEISR` (Likelihood Estimation of Individual Site Rates)](https://doi.org/10.1101/206011): Infer relative evolutionary rates from protein or nucleotide data

Full API documentation, including some code examples, is available from [http://sjspielman.org/phyphy](http://sjspielman.org/phyphy).


# Table of contents

  * [Installation](#installation)
  * [Using phyphy](#using-phyphy)
  	* [Overview](#overview)
  	* [Defining HyPhy instances](#defining-hyphy-instances)
  	* [Executing HyPhy Analyses](#executing-hyphy-analyses)
	* [Parsing HyPhy output JSON](#parsing-hyphy-output-json)
	* [Extracting CSVs from HyPhy output JSON](#extracting-csvs-from-hyphy-output-json)
	* [Parsing annotated trees from HyPhy output JSON](#parsing-annotated-trees-from-hyphy-output-json)
  * [Get help](#get-help)
  * [A note for conda users](#conda-help)
  * [Citation](#citation)


## Installation

You can obtain `phyphy` from pip (or pip3!) with `pip install phyphy`.

Note `phyphy` has the following dependencies (pip will take care of these for you, if necessary):

+ `Biopython >= 1.67` [**ONLY** `phyphy <=0.4.1`, dependency removed in version `>=0.4.2`]
+ `ete3 >=3.1`

You can update your installed version with `pip install --upgrade phyphy`, when needed.

Alternatively, you can download from source, via the usual `setuptools` procedure. Briefly:

1. Download the current release [here](https://github.com/sjspielman/phyphy/releases), or clone the master branch, and `cd` in. 
2. **Build** the package with `python setup.py build` (with a `sudo` as needed)
3. Optionally, run tests with `python setup.py test`  (with a `sudo` as needed). 
	+ NOTE: the test suite is `python2.7` compatible only (do not attempt with a `python3` interpretter, or some tests will fail). 
4. **Install** per your own adventure:
	> To specify a different install directory, add the argument `--prefix=/path/to/my/favorite/directory` to the install line.

	+ To install as root: `sudo python setup.py install`
	+ To install for user only: `python setup.py install --user`


## Using phyphy

### Overview

`phyphy` has three primary modules:

+ [**`Hyphy`**](http://sjspielman.org/phyphy/hyphy.html)
	+ This module can be **optionally** used to specify that a *non-canonical* installation (i.e. not installed to `/usr/local/`) or a local build  (i.e. where `make install` was not run) of HyPhy. This module can additionally be used to specify that HyPhy be run quietly and/or without outputting its standard log files, `messages.log` and `errors.log`.
+ [**`Analysis`**](http://sjspielman.org/phyphy/analysis.html)
	+ This module contains the analysis methods to execute, named according to the analysis. For example, the `FEL` class within the `Analysis` module would be used to execute an FEL analysis.
	+ Unless a custom `HyPhy` object is given, assumes the executable **`HYPHYMP`** in `usr/local/bin` and a HyPhy library in `/usr/local/lib/hyphy` (these are the default `make install`) locations
+ [**`Extractor`**](http://sjspielman.org/phyphy/extractor.html)
	+ This module contains functionality to parse a given analysis output. `Extractor` makes it simple to extract information from a given analysis JSON output, including:
		+ Fitted model parameters
	   + Fitted phylogenies
	   + [Newick-extended format](https://home.cc.umanitoba.ca/~psgendb/doc/atv/NHX.pdf) phylogenies with branch features for downstream visualization in tools like the Python package [`ete3`](http://etetoolkit.org/) or the R package [`ggtree`](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
	   + CSV files, for methods FEL, MEME, SLAC, FUBAR, LEISR, aBSREL

### Examples

#### Defining HyPhy instances


Full API documentation is [here](http://sjspielman.org/phyphy/hyphy.html). 

Most use cases are shown here:

```python
import phyphy
	
## Use canonical install, but suppress creation of messages.log and errors.log
myhyphy = phyphy.HyPhy(suppress_log = True)

## Use canonical install, but suppress markdown output	
myhyphy = phyphy.HyPhy(quiet = True)	

## Use a **local build**
myhyphy = phyphy.HyPhy(build_path = /path/to/build/of/hyphy/)
		
## Use a **local install**
myhyphy = phyphy.HyPhy(build_path = /path/to/install/of/hyphy/)

## Specify that <=3 processes should be used by the default executable
myhyphy = phyphy.HyPhy(CPU=3)

## Specify that the MPI executable should be used, with the launcher mpirun and the given mpirun arguments (32 processes)
myhyphy = phyphy.HyPhy(executable   = "HYPHYMPI",
                       mpi_launcher = "mpirun",
                       mpi_options  = "-np 32")

## Once defined,the HyPhy instance can be passed to an analysis, for example this default FEL inference (see next section for details):
myfel = FEL(hyphy = myhyphy, data = "/path/to/data.nex", nexus=True)
myfel.run_analysis()
```	
		
#### Executing HyPhy Analyses

Possible analyses to define include the following:

+ `ABSREL`
+ `BUSTED`
+ `FEL`
+ `FUBAR`
+ `LEISR`
+ `MEME`
+ `RELAX`
+ `SLAC`

Running an analysis proceeds in two steps:
1. Define an analysis instance (e.g. `x=FEL(..args..)` would define an `FEL` analysis)
2. Execute the analysis with the method `.run_analysis()` (e.g. `x.run_analysis()`)

Full API documentation is [here](http://sjspielman.org/phyphy/analysis.html). Each analysis method is documented with examples and the available optional arguments. Importantly, while each analysis will have its own optional arguments, *at a mininum*, you must provide path(s) to input data. There are two mutually exclusive (but one is necessary) strategies for this:

+ When alignment and tree are in a single file, use the single keyword argument `data`. 
+ When alignment and tree are in separate files, use the two respective keyword arguments `alignment` and `tree`.

Note that `RELAX` requires one additional argument, at a minimum: the label in the tree corresponding to test branches.

	
Briefly, here is how one might define and run some `FEL` analyses:

```python
import phyphy

### Run a FEL analysis, for example, with data in a single file.
### When providing a NEXUS file to the data argument, include the argument `nexus=True` (in version >=0.4.2)
myfel = phyphy.FEL(data = "path/to/data.nex", nexus=True)
myfel.run_analysis()

### Run a FEL analysis, for example, with data in two files
myfel = phyphy.FEL(alignment = "path/to/aln.fasta", tree = "path/to/tree/tree.txt")
myfel.run_analysis()

### Run a FEL analysis, specifying the vertebrate mtDNA genetic code (NIH code #2)
myfel = phyphy.FEL(data = "path/to/data.nex", genetic_code = 2, nexus=True)
myfel.run_analysis()	

### Run a FEL analysis, specifying a custom-named output JSON. This time, not a nexus!	
myfel = phyphy.FEL(data = "path/to/data.dat", output = "customname.json")
myfel.run_analysis()		
```


#### Parsing HyPhy output JSON

Full API documentation is available [here](http://sjspielman.org/phyphy/extractor.html), and [this PDF](examples/json-fields.pdf) describes the contents JSON fields in standard analyses. 

An `Extractor` instance should be defined using a single argument, **either** an executed `Analysis` instance or a specific JSON file:

```python

import phyphy
	
### Run a FEL analysis, for example, with data in a single file
myfel = phyphy.FEL(data = "path/to/data.nex")
myfel.run_analysis()
	
### Create an Extractor instance to parse JSON produced by `myfel`
myext = phyphy.Extractor(myfel)
	
## Alternatively, create an Extractor instance with a JSON directly
myext = phyphy.Extractor("path/to/json.json")
```

There are several flavors of `Extractor` methods, all of which are detailed with examples in [the API](http://sjspielman.org/phyphy/extractor.html):

+ Reveal contents of the JSON
	+ `.reveal_fields()` returns a list of all top-level fields in the JSON, to see the overall structure of the file
	+ `.reveal_fitted_models()` returns a list names of models fitted to the data, from which various components can be extracted
	+ `.reveal_branch_attributes()`, returns a dictionary of branch attributes, i.e. branch-specific information, (generally, these correspond to names of the fitted models which represent the fitted branch lengths, and other analysis-specific elements). 

+ Extract input information
	+ `.extract_number_sequences()` returns the number of sequences in the input data
	+ `.extract_number_sites()` returns the number of sites in the input data (note, this will be length/3 for codon analyses)
	+ `.extract_partition_count()` returns the number of partitions in the analysis
	+ `.extract_input_tree()` returns the original inputted phylogeny, with HyPhy node annotations
	+ `.extract_input_file()` returns the provided file name for the analyzed dataset

+ Extract fitted model components
	+ `.extract_model_logl(<name of model>)` returns the Log Likelihood of the fitted model
	+ `.extract_model_estimated_parameters(<name of model>)` returns the number of estimated parameters in the fitted model
	+ `.extract_model_aicc(<name of model>)` returns the small-sample AIC (AICc) of the fitted model
	+ `.extract_model_rate_distributions(<name of model>)` returns a dictionary (when applicable) of fitted model rate distributions
	+ `.extract_model_frequenciesl(<name of model>)` returns the empirical equilibrium frequencies of the fitted model
	+ `.extract_model_logl(<name of model>)`

+ Extract miscallaneous information
	+ `.extract_branch_sets()` returns a dictionary of structure `<node name>:<branch set>`. Can be rearranged to a dictionary of structure `<branch set>:[list of nodes]` with the argument `by_set=True` (or dictionary, with argument `as_dict=True`)
	+ `.extract_branch_attribute(<attribute name)` returns a dictionary of the desired attribute values. 
	+ `.extract_site_logl()` and `.extract_evidence_ratios()` are BUSTED-specific methods to return these values, as dictionaries each
	+ `.extract_timers()` returns a dictionary of timers from the method (wall-clock time in seconds to complete each step in algorithm)
+ Extract a CSV, described in the next section.


#### Extracting CSVs from HyPhy output JSON

CSVs can be obtained for the methods FEL, SLAC, MEME, FUBAR, ABSREL, and LEISR, using the method `.extract_csv(<csv_file_name>)`. Consult [the API](http://sjspielman.org/phyphy/extractor.html#extractor.Extractor.extract_csv) for information on column headers, and more usage information. 

Briefly:

```python 

import phyphy

### Define a FEL Extractor, for example
e = Extractor("/path/to/FEL.json") 
e.extract_csv("fel.csv")  ## save to fel.csv

### tab-delimited output, as fel.tsv
e.extract_csv("fel.tsv", delim = "\t")
```


#### Parsing annotated trees from HyPhy output JSON

Of specific interest, `phyphy` uses the powerful Python package `ete3` to assist in tree manipulation, allowing for the extraction of specific trees that can be used for downstream processing or visualization in other tools:

+ The method `.extract_input_tree()` allows you to obtain the original inputted phylogeny, with HyPhy node annotations
+ The method `.extract_model_tree()` allows you to obtain the fitted phylogeny for a given model (i.e., branch lengths will be updated). This will be output in standard newick format
+ The method `.extract_feature_tree()` allows you to obtain an **annotated** tree in Newick eXtended format (NHX), where nodes are annotated with the provided feature (i.e., attribute). 
+ The method `.extract_absrel_tree()` is a special case of `.extract_feature_tree()` for specifically annotating branches based on whether an aBSREL analysis has found **evidence for selection**, at a given P-value threshold
+ Note, for multipartitioned analyses, you can specify a partition or obtain all partitions from either of these methods

Please consult [the documentation](http://sjspielman.org/phyphy/extractor.html) for examples and full usage information. Some brief examples follow here:

To determine which models can be used with `.extract_model_tree()`, use the method `. reveal_fitted_models()` to get a list of all models which were fitted:

```python

import phyphy

## Define a FEL Extractor, for example
e = phyphy.Extractor("/path/to/FEL.json") 
e.reveal_fitted_models()
  ['Nucleotide GTR', 'Global MG94xREV']
  
## Obtain tree fitted during the Global MG94xREV fit
## Consult the API for more arguments and customizations
e.extract_model_tree('Nucleotide GTR')
'((((Pig:0.192554792971,Cow:0.247996722936)Node3:0.101719189407,Horse:0.211310618381,Cat:0.273732369855)Node2:0.0644249932769,((RhMonkey:0.00372054481786,Baboon:0.0017701670358)Node9:0.0259206344918,(Human:0,Chimp:0.00182836999996)Node12:0.0178636195889)Node8:0.109431753602)Node1:0.284434196447,Rat:0.0670087588444,Mouse:0.120166947697);'
```


Feature trees in NHX format can be created using any of the available **branch attributes**. To determine which attributes are available, use the method `.reveal_branch_attributes()` to get a list of all dictionary attributes, where the keys are the attribute names and the values are the type of attribute (either a branch length, node label, or branch label). Any of these attributes can then be used as a feature:


```python

import phyphy

## Define an ABSREL Extractor, for example
e = phyphy.Extractor("/path/to/ABSREL.json") 
e.reveal_branch_attributes()
	{'Corrected P-value': 'branch label', 
	 'Rate classes': 'branch label', 
	 'original name': 'node label', 
	 'Full adaptive model': 'branch length', 
	 'Rate Distributions': 'branch label', 
	 'Uncorrected P-value': 'branch label', 
	 'Baseline MG94xREV': 'branch length', 
	 'Nucleotide GTR': 'branch length', 
	 'Baseline MG94xREV omega ratio': 'branch label', 
	 'LRT': 'branch label'}
 
  
## Any of those keys can then be used, or multiple in a list, to extract as features:
## Again, please consult the API for more arguments and customizations
e.extract_feature_tree('LRT')
'(0564_7:1[&&NHX:LRT=0],(((((0564_11:1[&&NHX:LRT=3.96048976951],0564_4:1[&&NHX:LRT=4.80587881259])Node20:1[&&NHX:LRT=3.30060030447],(0564_1:1[&&NHX:LRT=0.0105269546166],(0564_21:1[&&NHX:LRT=0],0564_5:1[&&NHX:LRT=4.51927751707])Node25:1[&&NHX:LRT=0])Node23:1[&&NHX:LRT=0])Node19:1[&&NHX:LRT=0.153859379099],0564_17:1[&&NHX:LRT=0])Node18:1[&&NHX:LRT=1.64667962972],((0564_13:1[&&NHX:LRT=0],(0564_15:1[&&NHX:LRT=4.97443221859])Node32:1[&&NHX:LRT=0])Node30:1[&&NHX:LRT=2.86518293899],((0564_22:1[&&NHX:LRT=0.114986865638],0564_6:1[&&NHX:LRT=0])Node36:1[&&NHX:LRT=0],0564_3:1[&&NHX:LRT=14.0568340492])Node35:1[&&NHX:LRT=22.65142315])Node29:1[&&NHX:LRT=1.50723222708])Node17:1[&&NHX:LRT=2.63431127725],0564_9:1[&&NHX:LRT=0])Node16:1[&&NHX:LRT=0],(((0557_24:1[&&NHX:LRT=0],0557_4:1[&&NHX:LRT=0],0557_2:1[&&NHX:LRT=0])Node9:1[&&NHX:LRT=0],0557_12:1[&&NHX:LRT=0])Node8:1[&&NHX:LRT=1.99177154715],((0557_21:1[&&NHX:LRT=0],0557_6:1[&&NHX:LRT=0.662153642024],0557_9:1[&&NHX:LRT=0],0557_11:1[&&NHX:LRT=2.65063418443],0557_13:1[&&NHX:LRT=0],0557_26:1[&&NHX:LRT=0],(0557_5:1[&&NHX:LRT=1.98893541781],0557_7:1[&&NHX:LRT=0])Node53:1[&&NHX:LRT=0.660753167251])Node6:1[&&NHX:LRT=0],0557_25:1[&&NHX:LRT=0])Node7:1[&&NHX:LRT=1.69045394756])Separator:1[&&NHX:LRT=14.127483568])[&&NHX:LRT=0];'

## With specific branch lengths for a fitted model:
e.extract_feature_tree('LRT', update_branch_lengths = "Full adaptive model")
'(0564_7:0.00708844[&&NHX:LRT=0],(((((0564_11:0.00527268[&&NHX:LRT=3.96048976951],0564_4:0.00714182[&&NHX:LRT=4.80587881259])Node20:0.0022574[&&NHX:LRT=3.30060030447],(0564_1:0.00583239[&&NHX:LRT=0.0105269546166],(0564_21:0.00121537[&&NHX:LRT=0],0564_5:0.00266921[&&NHX:LRT=4.51927751707])Node25:0.000797211[&&NHX:LRT=0])Node23:0.00142056[&&NHX:LRT=0])Node19:0.0019147[&&NHX:LRT=0.153859379099],0564_17:0.00605582[&&NHX:LRT=0])Node18:0.00100178[&&NHX:LRT=1.64667962972],((0564_13:0.0053066[&&NHX:LRT=0],(0564_15:0.00346989[&&NHX:LRT=4.97443221859])Node32:0.000752206[&&NHX:LRT=0])Node30:0.00188243[&&NHX:LRT=2.86518293899],((0564_22:0.00686981[&&NHX:LRT=0.114986865638],0564_6:0.00581523[&&NHX:LRT=0])Node36:0.00125905[&&NHX:LRT=0],0564_3:0.00791919[&&NHX:LRT=14.0568340492])Node35:0.0174886[&&NHX:LRT=22.65142315])Node29:0.0010489[&&NHX:LRT=1.50723222708])Node17:0.00156911[&&NHX:LRT=2.63431127725],0564_9:0.00551506[&&NHX:LRT=0])Node16:0.000783733[&&NHX:LRT=0],(((0557_24:0.00078793[&&NHX:LRT=0],0557_4:0.000787896[&&NHX:LRT=0],0557_2:0.000399166[&&NHX:LRT=0])Node9:0.00206483[&&NHX:LRT=0],0557_12:0.00267531[&&NHX:LRT=0])Node8:0.00118205[&&NHX:LRT=1.99177154715],((0557_21:0[&&NHX:LRT=0],0557_6:0.000391941[&&NHX:LRT=0.662153642024],0557_9:0.000402021[&&NHX:LRT=0],0557_11:0.00156985[&&NHX:LRT=2.65063418443],0557_13:0.000401742[&&NHX:LRT=0],0557_26:0.00079377[&&NHX:LRT=0],(0557_5:0.00117641[&&NHX:LRT=1.98893541781],0557_7:0[&&NHX:LRT=0])Node53:0.000391973[&&NHX:LRT=0.660753167251])Node6:0.00118062[&&NHX:LRT=0],0557_25:0.00220372[&&NHX:LRT=0])Node7:0.00103489[&&NHX:LRT=1.69045394756])Separator:0.00822051[&&NHX:LRT=14.127483568])[&&NHX:LRT=0];'
```

ABSREL selection trees can be created in much the same way, using the method `.extract_absrel_tree()`:

```python

import phyphy

## Define an ABSREL Extractor
e = phyphy.Extractor("/path/to/ABSREL.json") 
## Selected defined as be all P<=0.1, and use Full Adaptive Model branch lengths
e.extract_absrel_tree(p = 0.1, update_branch_lengths = "Full adaptive model")
'(0564_7:0.00708844[&&NHX:Selected=0],(((((0564_11:0.00527268[&&NHX:Selected=0],0564_4:0.00714182[&&NHX:Selected=0])Node20:0.0022574[&&NHX:Selected=0],(0564_1:0.00583239[&&NHX:Selected=0],(0564_21:0.00121537[&&NHX:Selected=0],0564_5:0.00266921[&&NHX:Selected=0])Node25:0.000797211[&&NHX:Selected=0])Node23:0.00142056[&&NHX:Selected=0])Node19:0.0019147[&&NHX:Selected=0],0564_17:0.00605582[&&NHX:Selected=0])Node18:0.00100178[&&NHX:Selected=0],((0564_13:0.0053066[&&NHX:Selected=0],(0564_15:0.00346989[&&NHX:Selected=0])Node32:0.000752206[&&NHX:Selected=0])Node30:0.00188243[&&NHX:Selected=0],((0564_22:0.00686981[&&NHX:Selected=0],0564_6:0.00581523[&&NHX:Selected=0])Node36:0.00125905[&&NHX:Selected=0],0564_3:0.00791919[&&NHX:Selected=1])Node35:0.0174886[&&NHX:Selected=1])Node29:0.0010489[&&NHX:Selected=0])Node17:0.00156911[&&NHX:Selected=0],0564_9:0.00551506[&&NHX:Selected=0])Node16:0.000783733[&&NHX:Selected=0],(((0557_24:0.00078793[&&NHX:Selected=0],0557_4:0.000787896[&&NHX:Selected=0],0557_2:0.000399166[&&NHX:Selected=0])Node9:0.00206483[&&NHX:Selected=0],0557_12:0.00267531[&&NHX:Selected=0])Node8:0.00118205[&&NHX:Selected=0],((0557_21:0[&&NHX:Selected=0],0557_6:0.000391941[&&NHX:Selected=0],0557_9:0.000402021[&&NHX:Selected=0],0557_11:0.00156985[&&NHX:Selected=0],0557_13:0.000401742[&&NHX:Selected=0],0557_26:0.00079377[&&NHX:Selected=0],(0557_5:0.00117641[&&NHX:Selected=0],0557_7:0[&&NHX:Selected=0])Node53:0.000391973[&&NHX:Selected=0])Node6:0.00118062[&&NHX:Selected=0],0557_25:0.00220372[&&NHX:Selected=0])Node7:0.00103489[&&NHX:Selected=0])Separator:0.00822051[&&NHX:Selected=1])[&&NHX:Selected=0];'
```

Note that, for all analysis, HyPhy will rename taxa when needed to "safe" names when certain forbidden characters are seen. In most cases, taxon names are safe. By default, extracted trees will contain these safe names. To force the output to use the original names as provided in the alignment/tree analyzed, use the argument `original_names = True` with any of these tree extraction method.

Finally, any NHX tree can be visualized with a variety of programmatic platforms, including [`ete3`](http://etetoolkit.org/) in Python3 or [`ggtree`](https://bioconductor.org/packages/release/bioc/html/ggtree.html) in R/Bioconductor. Examples of visualizing NHX trees with either of these two platforms are available in [`examples/visualize_feature_tree_ete3.py`](./examples/visualize_feature_tree_ete3.py) and [`examples/visualize_feature_tree_ggtree.R`](./examples/visualize_feature_tree_ggtree.R), respectively.  




## Get help

+ Full documentation is available from [http://sjspielman.org/phyphy](http://sjspielman.org/phyphy)

+ Please post questions and bugs to the [Issues page](https://github.com/sjspielman/phyphy/issues).

+ [This PDF](examples/json-fields.pdf) contains a reference describing the contents JSON fields in standard analyses.

## A note for conda users

If you are using the `conda` distribution (either via `miniconda` or `anaconda`) of HyPhy on Linux, you may need to issue the following command (with the appropriate path inserted) before launching your `phyphy` script: 

```
export LD_LIBRARY_PATH=/path/to/miniconda2/lib:$LD_LIBRARY_PATH
```

## Citation
Spielman, SJ (2018). *phyphy: Python package for facilitating the execution and parsing of HyPhy standard analyses.* Journal of Open Source Software, 3(21), 514, [https://doi.org/10.21105/joss.00514](https://doi.org/10.21105/joss.00514)
