# phyphy

**IN DEVELOPMENT. NOT NECESSARILY EXPECTED TO WORK.**

`phyphy` aims to facilitate HyPhy usage in two primary ways:
1) Execute standard analyses in a Python scripting environment
2) Conveniently parse various information from the resulting JSON output.

`phyphy` is pronounced "feye-feye" and is so-named for "**P**ython **Hyphy**". It can be used with either Python2.7 or Python3.

Future functionality will allow users to define and fit custom HyPhy models.

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

Note that the forthcoming `2.3.7` release lives in the `beta` branch in [the HyPhy repository](https://github.com/veg/hyphy).

## Get help

+ Please post questions and bugs to the [Issues page](https://github.com/sjspielman/phyphy/issues) or contact `stephanie.spielman@temple.edu`.

+ [This PDF](data/json-fields.pdf) contains a reference describing the contents JSON fields in most HyPhy methods.


## Install
<!--
You can obtain `phyphy` from pip (or pip3!) with `pip install phyphy`.
-->
Dependencies:

+ `Biopython`
+ `dendropy >=4.3`

Until release, you must download `phyphy` from source, with the usual `setuptools` procedure. Briefly:

1. Download the master branch and `cd` in. 
2. **Build** the package with `python setup.py build` (with a `sudo` as needed)
3. Optionally, run tests with `python setup.py test`  (with a `sudo` as needed). 
	+ NOTE: the test suite is `python2.7` compatible only (do not attempt with a `python3` interpretter, or tests will fail). 
4. **Install** per your own adventure:
	> To specify a different install directory, add the argument `--prefix=/path/to/my/favorite/directory` to the install line.

	+ To install as root: `sudo python setup.py install`
	+ To install for user only: `python setup.py install --user`


## Sparknotes usage



<!--
## Example script

```python
	from phyphy import *
	
	
	## Create a HyPhy instance. In general this is NOT NEEDED if the appropriate version is install in your system as the default HyPhy. If you are using a locally intalled version, you need to specify path. 
	# Possible arguments include:
	### 1) executable, the desired executable to use (ie HYPHYMPI). Default: HYPHYMP
	### 2) path, the path to a **local hyphy build**. Only use this argument if you **do not** want to use the installed hyphy in /usr/local.
	### 3) cpu, the maximum number of CPUs per analysis. 
   ### 4) quiet, suppress screen output (still creates messages.log and errors.log, when applicable). Default: False
	hyphy = HyPhy(path = /path/to/my/local/hyphy/, quiet=True)
	
	## Make some variables
	codon_alignment = "path/to/my/alignment.fasta"
	tree = "path/to/my/tree.tre"
	data = "path/to/file/with/both/alignment/and/tree/inside.dat" 
	json = "path/to/where/id/rather/save/the/json"
	## Create an analysis of your choice instance. The following runs a one-rate FEL
	f = FEL(hyphy = hyphy, alignment = codon_alignment, tree = tree, two_rate = False, output = json) ## Use help() to see available arguments    
	## NOTE: This line could be used instead (data rather than alignment and tree):   f = FEL(hyphy = hyphy, data = data, two_rate = False, output = json)
    
    
    ## Run the analysis
    f.run_analysis()
```

-->
	
	
