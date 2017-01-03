CMAQ_Dev
========

Community Multiscale Air Quality Model [US EPA Website](https://www.epa.gov/air-research/community-multi-scale-air-quality-cmaq-modeling-system-air-quality-management)

CMAQ is an active open-source development project of the U.S. EPA Computational Exposure Division
that consists of a suite of programs for conducting air quality model simulations.
CMAQ is supported by the CMAS Center (http://www.cmascenter.org).

CMAQ combines current knowledge in atmospheric science and air quality modeling with multi-processor
computing techniques in an open-source framework to deliver fast, technically sound estimates of ozone,
particulates, toxics, and acid deposition.

##Getting the CMAQ Repository
The publically available CMAQ Git repository is organized with each official public release stored as a branch in the USEPA/CMAQ repository.
To download code, scripts, and documentation from the CMAQ Git repository, specify the branch (i.e. version number) and issue the following command from within a working directory on your server. For example to get the 5.2Beta version in a local repository that you've named CMAQ_v5.2Beta:
```
git clone -b 5.2Beta https://github.com/USEPA/CMAQ.git CMAQ_v5.2Beta  
```

##CMAQ Repository Guide
Source code and scripts are organized as follows:
* **CCTM (CMAQ Chemical Transport Model):** code and scripts for running the 3D-CTM at the heart of CMAQ.
* **PREP:** Data preprocessing tools for important input files like initial and boundary conditions, etc.
* **POST:** Data postprocessing tools for aggregating and evaluating CMAQ output products (e.g. Combine, Site-Compare, etc)
* **UTIL:** Utilities for generating code and using CMAQ (e.g. chemical mechanism generation, IO-API, etc)

##Documentation
Release notes and Code documentation are included within this repository (they are version-controlled along with the code itself).  

[CMAQv5.2 Documentation](CCTM/docs/User_Manual/User_Manual.md)   
[CMAQv5.2 Release Notes](CCTM/docs/Release_Notes/README.md)   

##CMAQ Test Cases
Benchmark/tutorial input and output data for each CMAQ release version are available from the CMAS Center Software Clearinghouse. From http://www.cmascenter.org, select Download -> Software -> CMAQ and choose the version to get the tutorial data.

