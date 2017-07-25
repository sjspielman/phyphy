# hyphyhelper

This repository contains various code bits for helping you to use HyPhy and analyze its output. 

**IN DEVELOPMENT. NOT NECESSARILY EXPECTED TO WORK.**

`hyphyhelper` provides functionality so that you can write and execute standard HyPhy analyses in a python environment. Never again worry about absolute paths or HyPhy idiosyncracies (mostly).

## Available Analyses

All implemented analyses are **only available** through the `v2.3-dev` branch in HyPhy. In other words, `hyphyhelper` will (almost entirely) **not work** with the current HyPhy release and/or `master` branch (as exists on 7/26/17).

+ RELAX
+ BUSTED
+ aBSREL
+ FEL
+ MEME
+ RelativeProteinRates
+ RelativeNucleotideRates

> + SLAC*
> 
> SLAC.bf is currently experiencing bugs so please don't try use SLAC in v2.3-dev yet.


## Example script

```{python}
	from hyphyhelper import *
	
	
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
	
	
	
	
	
	
	