# iSC.MEB
iSC.MEB: Integrate spatial clustering analysis with hidden Markov random field using empirical Bayes. iSC.MEB is a package for integrating and analyzing multiple spatially resolved transcriptomics (SRT) datasets, developed by the Jin Liu's lab. 

Check out our [Package Website](https://XiaoZhangryy.github.io/iSC.MEB/index.html) for a more complete description of the methods and analyses. 

iSC.MEB can be used to compare and contrast experimental datasets in a variety of contexts, for instance:

* Across experimental batches
* Across individuals
* Across different conditions (i.e., case and control)
* Across datasets with only partially shared cell/domain clusters

# Installation
To install the the packages "iSC.MEB", firstly, install the 'devtools' package. Besides, "iSC.MEB" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```{Rmd}
install.packages("devtools")
library(devtools)
install_github("XiaoZhangryy/iSC.MEB")
```

# Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [DLPFC data analysis](https://XiaoZhangryy.github.io/iSC.MEB/articles/iSC.MEB.DLPFC.html)
* [seqFISH data analysis](https://XiaoZhangryy.github.io/iSC.MEB/articles/iSC.MEB.seqFISH.html)


For the users that don't have set up system properly, the following setup on different systems can be referred.
## Setup on Windows system
First, download [Rtools](https://cran.r-project.org/bin/windows/Rtools/); second, add the Rtools directory to the environment variable. Users can follow [here](https://helpdeskgeek.com/windows-10/add-windows-path-environment-variable/#:~:text=Go%20ahead%20and%20click%20on%20the%20Environment%20Variables,you%20have%20to%20decide%20which%20one%20to%20edit) to add Windows PATH Environment Variable.


## Setup on MacOS system
First, install Xcode. Installation about Xcode can be referred [here](https://stackoverflow.com/questions/8291146/xcode-installation-on-mac#:~:text=You%20get%20it%20from%20the%20Mac%20App%20Store.,find%20the%20app%2C%20and%20click%20the%20install%20button).


Second, install "gfortran" for compiling C++ and Fortran at [here](https://github.com/fxcoudert/gfortran-for-macOS).


## Setup on Linux  system
For parallel computation on Linux, users must use the following system command to set the C_stack unlimited in case of the error `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

If you use conda environment on Linux system and some dependent packages (such as `scater`) can not normally installed, you can search R package at anaconda.org website. We take the `scater` package as example, and its search result is https://anaconda.org/bioconda/bioconductor-scater. Then you can install it in conda environment by following command.
```{Linux}

conda install -c bioconda bioconductor-scater
```
For the user not using conda environment, if  dependent packages (such as `scater`) not normally installed are in Bioconductor, then use the following command to install the dependent packages.
```{Linux}
# install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```
If  dependent packages (such as `DR.SC`) not normally installed are in CRAN, then use the following command to install the dependent packages.
```{Linux}
# install the package on CRAN
install.packages("DR.SC")
```


# Demonstration

For an example of typical iSC.MEB usage, please see our [Package Website](https://XiaoZhangryy.github.io/iSC.MEB/index.html) for a demonstration and overview of the functions included in iSC.MEB.

