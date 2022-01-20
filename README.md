# MockFBA

[![Build Status](https://github.com/shadaba/MockFBA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shadaba/MockFBA.jl/actions/workflows/CI.yml?query=branch%3Amain)
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


# Installation
  * Make sure you have Julia 1.7 installed: https://julialang.org/downloads/platform/
  * clone this repo
  * From the repo directory start julia
  * go to pkg mode by typing ]
  * Activate the project by : activate .
  * Install dependencies by typing : Instantiate
  * Now it is ready to use

# Python dependencies (following python libraries are needed):
  * yaml
  * numpy
  * datetime
  * glob
  * fitsio
  * multiprocessing


# How to use:
  Once the config file is setup simply type
  python MockFBA_run.py <config_file_name>

# How to setup config file:
  An example config file is given in config.yaml with detailed comments
  Please go through it and make sure you provide all the appropriate paths and file names


