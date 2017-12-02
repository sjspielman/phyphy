"""
    SJS.
    This example script provides a barebones example of how one might visualize a feature-annotated tree with ete3.
    ** Requires Python3 and the ete3 dependency PyQt5 (or PyQt4). **
    
    A wealth of ete3 documentation and examples are available from etetoolkit.org
"""
import ete3


### Example tree exported from phyphy, with .extract_absrel_tree(), with the feature `Selected` ###
mytree =  "(0564_7:0.00708844[&&NHX:Selected=0],(((((0564_11:0.00527268[&&NHX:Selected=0],0564_4:0.00714182[&&NHX:Selected=0])Node20:0.0022574[&&NHX:Selected=0],(0564_1:0.00583239[&&NHX:Selected=0],(0564_21:0.00121537[&&NHX:Selected=0],0564_5:0.00266921[&&NHX:Selected=0])Node25:0.000797211[&&NHX:Selected=0])Node23:0.00142056[&&NHX:Selected=0])Node19:0.0019147[&&NHX:Selected=0],0564_17:0.00605582[&&NHX:Selected=0])Node18:0.00100178[&&NHX:Selected=0],((0564_13:0.0053066[&&NHX:Selected=0],(0564_15:0.00346989[&&NHX:Selected=0])Node32:0.000752206[&&NHX:Selected=0])Node30:0.00188243[&&NHX:Selected=0],((0564_22:0.00686981[&&NHX:Selected=0],0564_6:0.00581523[&&NHX:Selected=0])Node36:0.00125905[&&NHX:Selected=0],0564_3:0.00791919[&&NHX:Selected=1])Node35:0.0174886[&&NHX:Selected=1])Node29:0.0010489[&&NHX:Selected=0])Node17:0.00156911[&&NHX:Selected=0],0564_9:0.00551506[&&NHX:Selected=0])Node16:0.000783733[&&NHX:Selected=0],(((0557_24:0.00078793[&&NHX:Selected=0],0557_4:0.000787896[&&NHX:Selected=0],0557_2:0.000399166[&&NHX:Selected=0])Node9:0.00206483[&&NHX:Selected=0],0557_12:0.00267531[&&NHX:Selected=0])Node8:0.00118205[&&NHX:Selected=0],((0557_21:0[&&NHX:Selected=0],0557_6:0.000391941[&&NHX:Selected=0],0557_9:0.000402021[&&NHX:Selected=0],0557_11:0.00156985[&&NHX:Selected=0],0557_13:0.000401742[&&NHX:Selected=0],0557_26:0.00079377[&&NHX:Selected=0],(0557_5:0.00117641[&&NHX:Selected=0],0557_7:0[&&NHX:Selected=0])Node53:0.000391973[&&NHX:Selected=0])Node6:0.00118062[&&NHX:Selected=0],0557_25:0.00220372[&&NHX:Selected=0])Node7:0.00103489[&&NHX:Selected=0])Separator:0.00822051[&&NHX:Selected=1])[&&NHX:Selected=0];"

### Read in to ete3, specifying format=1
t = ete3.Tree( mytree, format=1 )

## Define a treestyle, to show both leaf names and branch lengths in output
ts = ete3.TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True

## Define node styles, specifying that selected branches be colored red and non-selected branches be colored black.
style_selected = ete3.NodeStyle()
style_selected["vt_line_color"] = "red"
style_selected["hz_line_color"] = "red"

style_notselected = ete3.NodeStyle()
style_notselected["vt_line_color"] = "black"
style_notselected["hz_line_color"] = "black"

### Set style for nodes by traversing tree and applying styles where appropriate
### Note that all features are strings in the tree, and the `.Selected` feature comes directly from the feature string in the tree itself
for node in t.traverse("preorder"):
    if node.Selected=="1":
        node.set_style(style_selected)
    elif node.Selected=="0":
        node.set_style(style_notselected)

### Show in the interactive tree viewer
t.show(tree_style=ts)   
