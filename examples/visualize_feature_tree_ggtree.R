#
#    SJS.
#    This example script provides a barebones example of how one might visualize a feature-annotated tree with ggtree.
#    
#    A wealth of ggtree documentation and examples are available from https://bioconductor.org/packages/release/bioc/html/ggtree.html
#

require(ggtree)

### Read in annotated tree to ggtree with read.nhx(), where tree is contained in file "absrel_selection.tre"
t <- read.nhx("absrel_selection.tre")

### Visualize tree with branches colored by `Selected` feature
ggtree(t, aes(color = Selected)) + 
    scale_color_continuous(low='black', high='red') + ## set low (0, not selected) to black, high (1, selected) to red
    geom_tiplab() +                                   ## show tip labels in tree
    geom_treescale()                                  ## show a scale bar
    
