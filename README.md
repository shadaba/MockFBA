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
  * use following command to install (This installs all julia depndencies)
      - `using Pkg`
      - `Pkg.add(url="https://github.com/shadaba/MockFBA.git")`
  * once installed you need to setup the path to the installed library to do so type following
      - `julia -e "import MockFBA; println(pathof(MockFBA))"`
  * The path to source will be like: "path-to-dir/src/MockFBA.jl"
  * This will display the path to main src, define a shell variable to this path by adding following to bash
       - `export MOCKFBA_PATH="path-to-dir"`
  * **Note:** If you want to update to newer versions then simply re-run the `Pkg.add` step, but this may update the path as julia keeps older versions as it is. So to use the latest version after updating repeat the steps to update `MOCKFBA_PATH` variable otherwise it will continue to use old version of the installation.


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
  Once the config file is setup simply type (see the config.yaml in this repo for example)
   * `python $MOCKFBA_PATH/MockFBA_run.py -config_file <config_file_name> -ncpu <number of cpu>`
   * You can also list all the options by typing
      - `python $MOCKFBA_PATH/MockFBA_run.py --h`



# How to setup config file:
  An example config file is given in config.yaml with detailed comments
  Please go through it and make sure you provide all the appropriate paths and file names
  The master version of this config file can be found on this link:
  https://github.com/shadaba/MockFBA/blob/main/config.yaml

  Currently this is setup for the FirstGen mock, more development is needed for user to apply general selection for various kind of mocks and data files.

# NERSC-job submission example
To submit a job on NERSC, the best strategy is to use shared queue as this will keep the io to lower level and effectively costing less computing time hours, though it uses larger wall-clock time.

```#!/bin/bash
#SBATCH -C haswell
#SBATCH -q shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH -J FBA1
#SBATCH -e FBA1

#Activate the python environment with necessary python libraries
source activate Abacus

#cd to your working directory for this code
cd /global/homes/s/shadaba/Projects/MockFBA

#Make sure the config_file full path is given, even if it is in local directory
python $MOCKFBA_PATH/MockFBA_run.py -config_file /global/homes/s/shadaba/Projects/MockFBA/config.yaml -ncpu 4
```

The above example uses 4 cpu but one can ask for as many cpus as available on a single node. This cannot work with multiple nodes with the current capability. Note that the pre-process and post-process steps can use maximum of as many cpu as many tracer in the input but individual tracers are worked on serially. Therefore, if these steps are significant part of wall-time then it is more optimal to use number of cpu close to the number of tracer.

Using larger number of process means each process needs to load and compile the code which also takes time so it is important to keep in mind that that the additional burden is small. Best efficiency will be for single cpu but if each cpu is performing significant of computation then this consideration becomes negligible.
