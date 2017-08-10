"""
    Parse HyPhy JSON output
"""

import os
import re
import json

class JSONFields():
    """
        Define the fields used in JSONs. Essentially the terms.json namespace terms.
    """
    
    def __init__(self):
        
        self.input                = "input"
        self.input_info           = "info"
        self.input_filename       = "file name"
        self.input_sites          = "number of sites"
        self.input_sequences      = "number of sequences"
        self.input_npartitions    = "partition count"
        self.input_trees          = "trees"
        
        self.rate_distributions   = "rate distributions"
        
        self.analysis_description      = "analysis"
        self.analysis_description_info = "info"


        self.MLE                       = "MLE"
        self.MLE_headers               = "headers"
        self.MLE_content               = "content"
        

        

class HyPhyParser():
    
    def __init__(self, json_path, analysis = None):
        """
            Parse JSON output.
        """
        self.fields = JSONFields()
        self.json_path = json_path
        assert(os.path.exists(self.json_path)), "\nERROR: JSON file provided does not exist."
        
        ## Parse to dictionary
        self._unpack_json()
        
        ## Determine analysis method
        self.analysis = analysis
        if self.analysis is None:
            self._determine_analysis_from_json()

    ############################## PRIVATE FUNCTIONS #################################### 
    def _unpack_json(self):
        """
            Unpack JSON into dictionary and determine analysis type
        """ 
        with open (self.json_path, "rU") as f:
            self.json = json.load(f)
            
    
    def _determine_analysis_from_json(self):

        analyis_matches = ["FEL", "SLAC", "MEME", "FUBAR", "ABSREL", "RELAX", "BUSTED"]
        json_info = self.json[ self.fields.analysis_description ][ self.fields.analysis_description_info ].upper()
        
        for name in analyis_matches:
        
            find_analysis = re.search(name.upper(), json_info)
            if find_analysis is not None:
                self.analysis = name
                break
        
        assert(self.analysis is not None), "\nERROR: Could not determine analysis from JSON. Please ensure that the JSON is correctly formatted."
    
 
    def _extract_sitemethod_MLE_content(self, site_block):
        """
            Extract the CSV content component, RAW.
        """        
        raw_content = site_block[ self.fields.MLE_content]["0"]
        if self.analysis == "SLAC":
            raw_content = raw_content["by-site"]["AVERAGED"]
        return raw_content
        
        
                   
    def _parse_sitemethod_to_csv(self, delim):
        """
            Extract a CSV from a site-level method JSON, including FEL, SLAC, MEME
        """
        site_block =  self.json[ self.fields.MLE ]
        raw_header = site_block[ self.fields.MLE_headers ]
        raw_content = self._extract_sitemethod_MLE_content(site_block)

        final_header = delim.join( [x[0].replace(" ","_") for x in raw_header] ) + "\n"
        
        final_content = ""
        for row in raw_content:
            final_content += delim.join(str(x) for x in row) + "\n"
        
        with open(self.csv, "w") as f:
            f.write(final_header + final_content)
        
        


    ############################## PUBLIC FUNCTIONS #################################### 

    def extract_csv(self, csv, delim = ","):
        """
            Extract results to a CSV.
            
        """
        ###### TODO: THIS SHOULD ONLY BE ALLOWED FOR CERTAIN ANALYSES. #######
        ## assert(self.analysis in some_list_of_allowed_analyses) ##
        
        self.csv = csv
        
        if self.analysis in ["SLAC", "MEME", "FEL"]:
            self._parse_sitemethod_to_csv(delim)
            
    
    def detect_site_positive_selection(self, p):
        """
            Print out which sites are positively selected at given p-value (or posterior!) for a site method
        """
        
        self.p_column = {"FEL": 4, "MEME": 6, "SLAC": 8}
        self.FEL_columns = {"alpha": 0, "beta": 1}
        
        positive = {}
        
        content = self._extract_sitemethod_MLE_content( self.json[ self.fields.MLE ] )
        
        i = 0
        for row in content:
            i += 1
            p_value = float( row[ self.p_column[self.analysis] ] )
            # Don't assess FEL p-value if dS > dN
            if self.analysis == "FEL" and content[self.FEL_columns["alpha"]] > content[self.FEL_columns["beta"]]:
                continue
            if p_value <= p:
                positive[i] = p_value
        
        if len(positive) == 0:
            print ("\nNo sites were found to be under positive, diversifying selection by " +  self.analysis +  " at P <= " +  str(p))
        else:
            print("\nThe following sites were found to be under positive, diversifying selection by " +  self.analysis +  " at P <= " +  str(p) + " :")
            for site in positive:
                print(str(site) + "   (P=" + str(positive[site]) + ")")
            
                





# 
#     def 
#         """
#             Keys: [u'MLE', u'branch attributes', u'analysis', u'tested', u'data partitions', u'timers', u'fits', u'input']
#         """
#         site_field = self.json["MLE"]
#         
#         
#         
#  
# 

