---
title: 'phyphy: Python package for facilitating the execution and parsing of HyPhy standard analyses'
tags:
  - comparative sequence analysis
  - bioinformatics
  - viral evolution
  - phylogenetics
  - natural selection
  - molecular evolution
authors:
 - name: Stephanie J. Spielman
   orcid: 0000-0002-9090-4788
   affiliation: "1"
affiliations:
 - name: Institute for Genomics and Evolutionary Medicine, Temple Univeristy, Philadelphia, PA 19122
   index: 1
date: 4 December 2017
bibliography: paper.bib
---

# Summary

HyPhy [@HyPhy] is a powerful bioinformatics software platform for evolutionary comparative sequence analysis and testing hypotheses of natural selection from sequence data. Combined with its accompanying webserver Datamonkey [@datamonkey2], HyPhy has garnered over 3500 citations since its introduction in 2005 and has greatly accelerated the pace of biomedical and epidemiological research. I introduce phyphy (**P**ython **HyPhy**), a Python package aimed to i) faciliate the execution of standard HyPhy analyses, and ii) extract analysis information from the JSON-formatted HyPhy output into user-friendly formats. phyphy will greatly improve the batch-users' experience by allowing users to bypass the interactive HyPhy command-line prompt and execute hundreds or thousands of analyses directly from a Python script. In addition, phyphy makes it simple to obtain key information from an executed HyPhy analysis, including fitted model parameters, annotated phylogenies with analysis output formatted for downstream visualization, and CSV files containing the most relevant output for a given method. phyphy is compatible with Hyphy version >=2.3.7 and is freely available under a BSD-3 license from [https://github.com/sjspielman/phyphy](https://github.com/sjspielman/phyphy).

# References