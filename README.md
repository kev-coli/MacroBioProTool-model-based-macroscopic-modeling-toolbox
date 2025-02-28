# MacroBioProTool : Virtual Bioprocessing Toolbox (Version 1.0)

## Authors
Mirko Pasquini (Simulations and Optimization), Kevin Colin (Modelling)

## Project Description
This library provides a comprehensive set of tools for model-based approaches in continuous bioprocessing. It enables the derivation of mathematical models of cell metabolism from raw data, the simulation of steady-state experiments using such models (i.e. in a digital twin fashion), and the optimization of feed media composition based on model-based optimization techniques. The goal is to support the continuous bioprocessing community in its transition towards digitalization, by offering a complete model-based pipeline (from raw data to optimized process design) to allow for an efficient and effective bioprocess development.

## Features
- **Macroscopic Model identification**
- **Steady-State Experiments Simulations**
- **Model-based Media Optimization**

## Technology
- **MATLAB** (version 2022a or later)

## Installation Guide

1. Clone the repository and add the whole folder to Matlab path
2. Install all the necessary (or recommended) dependencies below
3. Functions can directly be called in Matlab scripts (check ```usecase-examples``` folder to start)

### Dependencies
REQUIRED
  - Matlab version R2022a or newer
  - Matlab Optimization Toolbox
  - Matlab Parallel Computing Toolbox
(HIGHLY) RECOMMENDED
  - Mosek solver (https://www.mosek.com/)

Please use the script ```installation_test``` in the folder ```utils``` to check the dependencies requirements are met. If mosek is not installed, this should be specified in the macroscopic-modelling code (note that computations can be slow and results unexpected if mosek is not used). 

## Usage
Usecase examples are the best point to get started with the library. These can be found in the folder ```usecase-examples``` and present examples of model loading, kinetic parameters matrix/vector conversion and handling, model-based nominal optimization with and without ranking based on sensitivity and distance from available dataset.

A complete introductory guide to the library can be found in the repository.

The ```help``` feature in Matlab can be used for all functions to check specification. 

## Some troubleshooting and F.A.Q.

- **Q.1. My simulations/optimization returns very strange results (e.g. much different than what I expected). What is happening?**

It is difficult to answer this because many things could be happening: the model is not perfect, the optimization or non-linear equation solvers fail, etc. However, in our experience, in many cases the answer is simpler and usually related to external files not following the required format. For example the order of metabolites should be the same in all the model files (containing the stoichiometric matrix and the parameters matrix) and the data files. Most of the times this is not the case, and there are functions in the library that can pre-process these files so that compatibility is restored. These functions, to work properly, will read the metabolite labels in such files. _THESE LABELS SHOULD BE THE SAME IN EACH FILE._ ``Glc_ext``, ``Glc_ex``, ``Glc`` and ``Glc ext`` are all different labels and the above pre-processing functions will give you an error in the best-case, or produce unreliable results in the worst-case. Check that your external files follow the correct template/format for the library (you can find examples of models and data files in the folders ```models``` and ```data``` respectively).  

- **Q.2. I would like to use a function \texttt{xyz} but I do not know what are the input or the output parameters, and this document does not specify them. What can I do?**

You can access the documentation of each function by using the ```help``` command in Matlab, followed by the name of the function you want to use.

- **Q.3. I have no idea where to start. What is the best way to start using the library?**

A good $80\%$ of the library functionalities can be explored by simply taking a look at the usecase examples (in the ```usecase-examples``` library folder). These scripts will guide you through basic features of the library (loading a model, running a virtual experiment, optimizing using nominal optimization, etc.). Feel free to copy-paste these scripts and adjust them as you need.

- **Q.4. When I run my code I get an error saying ``Unrecognized function or variable [...]''. What should I do?**

Be sure to run the code in the correct folder (e.g. most of the usecase examples should be executed while being in the folder ```usecase-examples```) and be sure that all the library folders are visible to MATLAB (this can be done through Set Path in the Home tab of MATLAB). 

## How to Cite
If you use this library in your research or your work, please cite it as follows:

```bibtex
  @misc{virtual_bioprocessing,
    author = {Kévin Colin and Mirko Pasquini},
    title = {MacroBioProTool : Bioprocessing Modeling, Simulation, and Optimization Library},
    year = {2025},
    url = {https://github.com/kev-coli/MacroBioProTool-model-based-macroscopic-modeling-toolbox},
  }
```

