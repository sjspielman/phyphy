from phyphy import *

"""
There are three main classes:

1) HyPhy
2) Analysis (parent, children are the analysis names)
3) HyPhyParser (needs a better name probably, but I didnt want to just do Parser. Maybe Extractor? we'll see)

The HyPhy class is **only** needed if you want to use a local or non-installed build, or use MPI. If you are using a hyphy version installed in the normal way, don't worry about this.
Otherwise, the docstring here (in src/phyphy_runner.py) should explain what you need. If you end up making a HyPhy instance, pass it into the Analysis:
However, it is USEFUL if you want hyphy to shut the hell up. Use quiet=T.

The Analysis class is the action. Options include:

1) data = ... OR alignment = .. and tree = .. . The data argument assumes concatenated data file, otherwise given alignment and tree separately. Relative paths.
2) If you have a HyPhy class, pass it in w/ keyword hyphy = ...
3) NOte that by default the json output goes to your alignment path, with the added suffix ".METHOD.json", ie "data.fasta.FEL.json"

Examples, with the arguments you'll want (mixing up data and tree/alignment to demonstrate):
x = FEL(data = /path/to/data.fna, srv = FALSE, output = /path/to/where/you/want/the/json) ## srv is synonymous rate varation, default is true
x.run_analysis()

x = SLAC(hyphy = theinstanceifyouhaveone, tree = /path/to/tree.nwk, alignment = /path/to/alignment.phy)
x.run_analysis()

## compatible with "beta" branch, not master or release, and actually I need to fix it for this too, soooo don't use this yet.
x = LEISR(data = /path/to/data.fna, type = "protein", model = "JC69"
x.run_analysis()



Then the parser class when you're done running. Takes either as an argument: the Analysis object (x) or path to a json.
It determines the analysis automatically. You'll mostly want the method .extract_csv()
p = HyPhyParser(x)
p.extract_csv("csvpath.csv")


Note, I've been using extension .fna to refer to file w/ fasta and tree.
"""



######## RUN LEISR: By far the fastest, but to use w/ phyphy you need the beta branch.
h = HyPhy(build_path = "/Users/sjspielman/hyphy/") ### build_path is where I have the beta branch locally built
x = LEISR(hyphy = h, data = "data/aa.fna", model = "WAG", type = "protein")
x.run_analysis()

y = HyPhyParser(x)
y.extract_csv("aa.csv")
assert 1==45
######## Run FEL, my /usr/local/lib install is >=2.3.6 so no need for HyPhy instance
x = FEL(data = "data/codon.fna", srv=False)
x.run_analysis()

y = HyPhyParser(x)
y.extract_csv("fel.csv")



######## Run SLAC, my /usr/local/lib install is >=2.3.6 so no need for HyPhy instance, except here I show you how to get hyphy to shut up with quiet (no stdout) and suppress_log (sends *log to /dev/null)
h = HyPhy(quiet=True, suppress_log=True)
x = SLAC(hyphy=h, data = "data/codon.fna")
x.run_analysis()

y = HyPhyParser(x)
y.extract_csv("slac.csv") ### ## Add argument: slac_ancestral_type = "RESOLVED" if you want resolved instead of averaged ancestral states, it matters about 0 but you never know..














