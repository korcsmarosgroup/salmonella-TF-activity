## Tutorial

- [Usage - R](#usage)
  - [Install and load software](#install-and-load-software)
  - [Importing the PKN](#importing-the-PKN)
  - [Loading DE results](#importing-the-ontology)
  - [Importing the PKN](#importing-the-ontology)



### Usage - R

#### Install and load software

SalmonAct was tailored to be compatible with [decoupleR](https://saezlab.github.io/decoupleR/), an activity inference library capable of running multiple relevant methods.

decoupleR can be installed using the following commands:
`install.packages('BiocManager')`
`BiocManager::install('saezlab/decoupleR')`

Once the installation has finished, load the library using:
`library(decoupleR)`

#### Importing SalmonAct
Download the PKN hosted on this GitHub page (`data/db/salmonact.tsv`)
Read the file and store it as a variable, e.g.:
`net <- read.table(file = 'data/db/salmonact.tsv', sep='\t', header = T)`

#### Loading DE results
Load the differential expression results from the experimental data to be analysed
