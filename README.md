# MockFBA

[![Coverage](https://codecov.io/gh/shadaba/MockFBA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/shadaba/MockFBA.jl)

# Goals:
  * Fast and efficient fiber-assignment mitigation for mocks
  * Easy to use, uses simple config file providing all information at one place
  * One combined fits file as final output, easy to understand final result
  * Efficiency as top priorities with additional features for mocks only
  * Follows the principle of minimal IO, and minimal memory allocation
  * Should be useful for almost any analysis using mocks and concerned about assignment
  * Can run on laptop as uses minimal memory, except in few preprocessing steps

# TO DO:
  * Capibility to assign STDs 
  * Capibility to account for SKY FIBERs
  * Additional focal plane constraints (example petal edges)
  * Field-Rotation
  * Parsing the evolving focal plane status


# Installation (USER)
  * Make sure you have Julia 1.7 installed: https://julialang.org/downloads/platform/
  * start julia by typing julia
  * use following command to install (This installs all julia depndencies)\\
      `using Pkg`
      `Pkg.add(url="https://github.com/shadaba/MockFBA.git")`
  * once installed you need to setup the path to the installed library to do so type following\\
      `julia -e "import MockFBA; println(pathof(MockFBA))"`
  * The path to source will be like: "path-to-dir/src/MockFBA.jl"
  * This will diplay the path to main src, define a shell variable to this path by adding following to bash
       export MOCKFBA_PATH="path-to-dir"
  * Note: If you want to update then simply re-run the Pkg.add step, but this may update the path as julia keeps older versions as it is. So to use the lates version after updating repeat the steps to update MOCKFBA_PATH variable otherwise it will continue to use old version of the installation.


# Installation (Developer)
  * Make sure you have Julia 1.7 installed: https://julialang.org/downloads/platform/
  * clone this repo
  * From the repo directory start julia
  * go to pkg mode by typing ]
  * Activate the project by : activate .
  * set MOCKFBA_PATH to the directory cloned  

# Python dependencies:
  * yaml
  * numpy
  * datetime
  * glob
  * fitsio
  * multiprocessing


# How to use:
  Once the config file is setup simply type (see the config.yaml in this repo for example)\\
   `python $MOCKFBA_PATH/MockFBA_run.py <config_file_name>`
  

# How to setup config file:
  An example config file is given in config.yaml with detailed comments
  Please go through it and make sure you provide all the appropriate paths and file names
  The master version of this config file can be found on this link:
  https://github.com/shadaba/MockFBA/blob/main/config.yaml



