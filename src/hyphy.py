#!/usr/bin/env python

##############################################################################
##  phyhy: *P*ython *HyPhy*: Facilitating the execution and parsing of standard HyPhy analyses.
##
##  Written by Stephanie J. Spielman (stephanie.spielman@temple.edu) 
##############################################################################

	
	
"""
    Define a custom HyPhy build/install.
"""


import sys
import os
import subprocess

if __name__ == "__main__":
    print("\nThis is the HyPhy module in `phyphy`. Please consult docs for `phyphy` usage." )
    sys.exit()

_DEFAULT_LIBPATH = "/usr/local/lib/hyphy/"


class HyPhy():
    """
        This class creates a HyPhy instance. Generally this is only necessary to use if any of these applies:
            + You wish to use a local **build** of HyPhy (not a canonically installed build)
            + You wish to use a local **install** of HyPhy (installed elsewhere from /usr/local)
            + You wish to use a different HyPhy executable from the default, which is HYPHYMP
    """

    def __init__(self, **kwargs):
        """
            Initiliaze a :code:`HyPhy()` instance.
        
            Optional keyword arguments:
                1. **executable**, the desired executable to use (ie HYPHYMPI). Default: HYPHYMP.
                2. **build_path**, the path to a **local hyphy build**. Use this argument if you have compiled hyphy in the downloaded hyphy/ directory and **did not run make install**
                3. **install_path**, the path to a **hyphy install**. Use this argument if you have specified a different installation path for hyphy, i.e. you provided `-DINSTALL_PREFIX=/other/path/` to cmake.
                4. **cpu**, the maximum number of processes per analysis. By default, HyPhy will take as many CPUs as it can/requires. This argument will limit the maximum. Default: None
                5. **quiet**, suppress screen output (Note, HyPhy will still creates messages.log and errors.log files, when applicable). Default: False
                6. **suppress_log**, suppress messages.log and errors.log files. Default: False. 
                7. **mpi_launcher**, mpi launcher. Default: :code:`mpirun`. Use this argument if are you specifying `HYPHYMPI` for executable.
                8. **mpi_options**, options to pass to the mpi launcher. Default: "".
                

            **Examples:**
               
               >>> ### Define a HyPhy() object, specifying a quiet analysis (no console output) with no log files output
               >>> my_hyphy = HyPhy(quiet = True, suppress_log = True)

               >>> ### Define a HyPhy() object, specifying that a maximum of 4 processes be used (for the default executable HYPHYMP)
               >>> my_hyphy = HyPhy(CPU = 4)
               
               >>> ### Define a HyPhy() object, specifying a specific *build path*
               >>> my_hyphy = HyPhy(build_path = "/path/to/local/build/hyphy/")

               >>> ### Define a HyPhy() object, specifying a specific *install path*
               >>> my_hyphy = HyPhy(install_path = "/path/to/local/install/hyphy/")               

               >>> ### Define a HyPhy() object, specifying that HYPHYMPI be used, with associated arguments for mpirun to use 32 processes
               >>> my_hyphy = HyPhy(executable = "HYPHYMPI", mpi_launcher = "mpirun", mpi_options = "-np 32")   
        """                                
                

        executable         = kwargs.get("executable", "HYPHYMP")
        self.build_path    = kwargs.get("build_path", None)         ### path if *built* locally (i.e. make install was NOT RUN)
        self.install_path  = kwargs.get("install_path", None)       ### path if installed locally (i.e. make install was run to a specified local path)
        self.cpu           = kwargs.get("cpu", None)                ### For use with MP
        self.mpi           = kwargs.get("mpi_launcher", "mpirun")   ### launcher for mpi
        self.mpiopts       = kwargs.get("mpi_options", "")          ### To pass to mpi environment   
        self.quiet         = kwargs.get("quiet", False)             ### If True, run hyphy quietly (no stdout/err)
        self.suppress_log  = kwargs.get("suppress_log", False)      ### If True, send messages.log, errors.log to /dev/null

      
        ### Checks for a local BUILD  ###
        if self.build_path is not None: 
            assert(os.path.exists(self.build_path)), "\n[ERROR] Build path does not exist."
            self.build_path = os.path.abspath(self.build_path) + "/" ## os.path.abspath will strip any trailing "/"
            self.libpath = self.build_path + "res/"
            assert(os.path.exists(self.libpath)), "\n[ERROR]: Build path does not contain a correctly built HyPhy."
            self.executable = self.build_path + executable
            self.hyphy_call = self.executable + " LIBPATH=" + self.libpath
        
        else: 
            ### Checks for a nonstandard (i.e. not in /usr/local/) INSTALL  ###
            if self.install_path is not None:
                assert(os.path.exists(self.install_path)), "\n[ERROR]: Install path does not exist."
                self.install_path = os.path.abspath(self.install_path) + "/"
                self.libpath = self.install_path + "lib/hyphy/"
                assert(os.path.exists(self.libpath)), "\n[ERROR]: Install path does not contain a correctly built HyPhy."               
                self.executable = self.install_path + "bin/" + executable
                self.hyphy_call = self.executable + " LIBPATH=" + self.libpath
            ## Installed in default path
            else:
                self.libpath = _DEFAULT_LIBPATH
                self.executable = executable
                self.hyphy_call = executable
        
        
        ## Ensure executable exists somewhere
        with open("/dev/null", "w") as hushpuppies:
            exit_code = subprocess.call(["which", self.executable], stdout = hushpuppies, stderr = hushpuppies) # If you're reading this, I hope you enjoy reading hushpuppies as much as I enjoyed writing it. --SJS
            if exit_code != 0:
                raise AssertionError("\n[ERROR]: HyPhy executable not found. Please ensure it is properly installed, or in your provided local path.")
        
        if executable == "HYPHYMPI":
            with open("/dev/null", "w") as hushpuppies:
                exit_code = subprocess.call(["which", self.mpi], stdout = hushpuppies, stderr = hushpuppies) 
                if exit_code != 0:
                    raise AssertionError("\n[ERROR]: MPI launcher not found (the default is `mpirun`).")    
            self.hyphy_call = self.mpi + " " + self.mpiopts + " " + self.hyphy_call
        else:
            if self.cpu is not None:
                self.hyphy_call += " CPU=" + str(self.cpu)
        
        if self.suppress_log is True:
            self.hyphy_call += " USEPATH=/dev/null/"

        
        