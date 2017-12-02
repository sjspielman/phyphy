# phyphy

**IN DEVELOPMENT. NOT NECESSARILY EXPECTED TO WORK.**

`phyphy` aims to facilitate HyPhy usage in two primary ways:
1) Execute standard analyses in a Python scripting environment
2) Conveniently parse various information from the resulting JSON output.

`phyphy` is pronounced "feye-feye" and is so-named for "**P**ython **Hyphy**". Importantly, this name was chosen for increased pronunciation enjoyment. ``phyphy`` is compatible with either Python2.7 or Python3.


## Available Analyses

`phyphy` is compatible with HyPhy version >=2.3.7 and supports the following analyses (in alphabetical order):

+ aBSREL
+ BUSTED
+ FEL
+ FUBAR
+ LEISR
+ MEME
+ RELAX
+ SLAC


## Get help

+ Full documentation is available from [http://sjspielman.org/phyphy](http://sjspielman.org/phyphy)

+ Please post questions and bugs to the [Issues page](https://github.com/sjspielman/phyphy/issues) or contact `stephanie.spielman@temple.edu`.

+ [This PDF](data/json-fields.pdf) contains a reference describing the contents JSON fields in standard analyses.


## Installation

You can obtain `phyphy` from pip (or pip3!) with `pip install phyphy`.

Note `phyphy` has the following dependencies (pip will take care of these for you, if necessary):

+ `Biopython >= 1.67`
+ `ete3 >=3.1`

Alternatively, you can download from source, via the usual `setuptools` procedure. Briefly:

1. Download the current release [here](https://github.com/sjspielman/phyphy/releases), or clone the master branch, and `cd` in. 
2. **Build** the package with `python setup.py build` (with a `sudo` as needed)
3. Optionally, run tests with `python setup.py test`  (with a `sudo` as needed). 
	+ NOTE: the test suite is `python2.7` compatible only (do not attempt with a `python3` interpretter, or some tests will fail). 
4. **Install** per your own adventure:
	> To specify a different install directory, add the argument `--prefix=/path/to/my/favorite/directory` to the install line.

	+ To install as root: `sudo python setup.py install`
	+ To install for user only: `python setup.py install --user`


## Sparknotes usage

### Overview

`phyphy` has three primary modules:

+ `Hyphy`
	+ This module can be **optionally** used to specify that a *non-canonical* installation (i.e. not installed to `/usr/local/`) or a local build  (i.e. where `make install` was not run) of HyPhy. This module can additionally be used to specify that HyPhy be run quietly and/or without outputting its standard log files, `messages.log` and `errors.log`.
+ `Analysis`
	+ This module contains the analysis methods to execute, named according to the analysis. For example, the `FEL` class within the `Analysis` module would be used to execute an FEL analysis.
	+ Unless a custom `HyPhy` object is given, assumes the executable **`HYPHYMP`**
+ `Extractor`
	+ This module contains functionality to parse a given analysis output. `Extractor` makes it simple to extract information from a given analysis JSON output, including:
		+ Fitted model parameters
	   + Fitted phylogenies
	   + [Newick-extended format](https://home.cc.umanitoba.ca/~psgendb/doc/atv/NHX.pdf) phylogenies with branch features for downstream visualization in tools like the Python package [`ete3`](http://etetoolkit.org/) or the R package [`ggtree`](https://bioconductor.org/packages/release/bioc/html/ggtree.html)
	   + CSV files, for methods FEL, MEME, SLAC, FUBAR, LEISR, aBSREL

### Jupyter notebooks

The follow Jupyter notebooks show various example of `phyphy` usage:





### Examples

#### Defining HyPhy instances

Full API documentation is [here](http://sjspielman.org/phyphy/hyphy.html). Some brief examples are shown below.

```python
	import phyphy
	
	## Use canonical install, but suppress creation of messages.log and errors.log
	myhyphy = phyphy.HyPhy(suppress_log = True)
	
	## Use canonical install, but suppress markdown output	myhyphy = phyphy.HyPhy(quiet = True)	
	
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

	## Once defined,the HyPhy instance can be passed to an analysis, for example this default FEL inference:
	myfel = FEL(hyphy = myhyphy, data = "/path/to/data.nex")	myfel.run_analysis()
```	
		
#### Executing HyPhy Analyses

Full API documentation is [here](http://sjspielman.org/phyphy/analysis.html). Each analysis will have its own optional arguments, as detailed in the API. However, *at a mininum*, you must provide path(s) to input data. There are two mutually exclusive (but one is necessary) strategies for this:

+ When alignment and tree are in a single file, use the single keyword argument `data`. 
+ When alignment and tree are in separate files, use the two respective keyword arguments `alignment` and `tree`.
	
Possible analyses to define include the following:
	+ aBSREL
	+ BUSTED
	+ FEL
	+ FUBAR
	+ LEISR
	+ MEME
	+ RELAX
		+ Note that this analysis has one other required argument, the label in the tree corresponding to test branches.
	+ SLAC


```python
	import phyphy
	
	
	### Run a FEL analysis, for example, with data in a single file
	myfel = phyphy.FEL(data = "path/to/data.nex")
	myfel.run_analysis()
	
	### Run a FEL analysis, for example, with data in two files
	myfel = phyphy.FEL(alignment = "path/to/aln.fasta", tree = "path/to/tree/tree.txt")
	myfel.run_analysis()
	
	### Run a FEL analysis, specifying the vertebrate mtDNA genetic code (NIH code #2)
	myfel = phyphy.FEL(data = "path/to/data.nex", genetic_code = 2)
	myfel.run_analysis()	

	### Run a FEL analysis, specifying a custom-named output JSON
	myfel = phyphy.FEL(data = "path/to/data.nex", output = "customname.json")
	myfel.run_analysis()		
```


#### Parsing HyPhy output JSON

Full API documentation is [here](http://sjspielman.org/phyphy/extractor.html). Again, [this PDF](data/json-fields.pdf) describes the contents JSON fields in standard analyses available here. 

An `Extractor` instance should be defined using a single argument, **either** an executed `Analysis` instance or a specific JSON file:

```python
	import phyphy
	
	
	### Run a FEL analysis, for example, with data in a single file
	myfel = phyphy.FEL(data = "path/to/data.nex")
	myfel.run_analysis()
	
	### Create an Extractor instance to parse JSON produced by `myfel`
	myext = phyphy.Extractor(myfel)
	
	## Alternatively, create an Extractor instance with a JSON directly
	myext = Extractor("path/to/json.json")
```

A Jupyter notebook demonstrating usage of the `Extractor` class will be available shortly. 

Of interest, `phyphy` uses the powerful Python package `ete3` to assist in tree manipulation, allowing for the extraction of specific trees that can be used for downstream processing or visualization in other tools:
	+ The method `.extract_model_tree()` allows you to obtain the fitted phylogeny for a given model (i.e., branch lengths will be updated). This will be output in standard newick format
	+ The method `.extract_feature_tree()` allows you to obtain an **annotated** tree in Newick eXtended format (NHX), where nodes are annotated with the provided feature (i.e., attribute). 
	+ The method `.extract_absrel_tree()` is a special case of `.extract_feature_tree()` for specifically annotating branches based on whether an aBSREL analysis has found evidence for selection, at a given P-value threshold
	+ Note, for multipartitioned analyses, you can specify a partition or obtain all partitions from either of these methods

Any NHX tree can be visualized with a variety of programmatic platforms, including [`ete3`](http://etetoolkit.org/) in Python3 or [`ggtree`](https://bioconductor.org/packages/release/bioc/html/ggtree.html) in R/Bioconductor. Examples of creating such trees and visualizing them with either of these two platforms are available in `examples/` (FORTHCOMING).  

