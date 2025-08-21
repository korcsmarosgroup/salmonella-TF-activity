## Tutorial
This tutorial will guide you through a typical TF activity analysis. Essentially, the steps mirror a [normal decoupleR analysis](https://saezlab.github.io/decoupleR/articles/tf_bk.html), but instead of using the human PKN CollecTRI (Müller-Dott et al., 2023), we are using SalmonAct as the regulon resource, and naturally using corresponding *Salmonella* transcriptomics data. 


### Usage - R

#### Install and load software

SalmonAct was tailored to be compatible with [decoupleR](https://saezlab.github.io/decoupleR/), an activity inference library capable of running multiple relevant methods including their consensus.

decoupleR can be installed using the following commands:
```
install.packages('BiocManager')
BiocManager::install('saezlab/decoupleR')
```
Once the installation has finished, load the library using:
```
library(decoupleR)
```

#### Importing SalmonAct
Download the PKN hosted on this GitHub page (`data/db/salmonact.tsv`)
Read the file and store it as a variable, e.g.:
```
net <- read.table(file = 'salmonact.tsv', sep='\t', header = T)
```

#### Loading DE results
Load the differential expression results from the experimental data to be analysed, and map the identifiers to match the format of PKN. We used "Alternate SL1344 identifiers" (SL0000 format) for the resource. For S. Typhimurium strains we recommend using the supplementary table 2 published in `Canals et al., 2019`, uploaded to this repository under `data/analysis/d23480_annotation.tsv`. For further removed strains the protein orthology based mapping table may be used, available under `orthology_mapping/Mapping_table_salmonact.tsv`.

Below are the corresponding steps from the notebook replicating the analysis in the manuscript (`scripts/salmonact_notebook.Rmd`):

Loading mapping data from `Canals et al., 2019` for S. Typhimurium strains, selecting LT2 and SL1344 identifiers:
```
ann <- read_tsv("../data/analysis/d23580_annotation.tsv")
annotation_data <- read_tsv("../data/analysis/d23580_annotation.tsv") %>%
  select(`LT2 orthologue`, `Alternate SL1344 identifier`) %>%
  dplyr::rename(LT2 = `LT2 orthologue`, SL1344 = `Alternate SL1344 identifier`)
```
The resulting mapping table will look like this:
```
# A tibble: 5,126 × 2
   LT2     SL1344
   <chr>   <chr> 
 1 STM0001 SL0001
 2 STM0002 SL0002
 3 STM0003 SL0003
 4 STM0004 SL0004
 5 STM0005 SL0005
 6 STM0006 SL0006
 7 STM0007 SL0007
 8 STM0008 SL0008
 9 STM0009 SL0009
10 STM0010 SL0010
# ℹ 5,116 more rows
# ℹ Use `print(n = ...)` to see more rows
```
Loading the experiment data, in this case the ompR knockout results from `Perkins et al., 2013`. 


```
experiment <- read_tsv('../data/analysis/GSE35938.top.table.tsv')
```
To map the identifiers you will need to join the mapping table to the experimental results using the shared identifier from the mapping table. In this experiment the `ORF` column of the differential expression results contained LT2 locus tags, which we could map to using `left_join(annotation_data, by = c('ORF'='LT2'))`

Once complete, select the relevant columns and the mapped SL0000 identifier column (`select(SL1344, logFC, t, P.Value)`). 

We will need to make sure we drop all NA values, and do not have any duplicated gene names (`drop_na`, `distinct` functions). Finally, we will turn the object into a matrix, with the identifiers in the row names. 
```
# Map IDs, extract t-values per gene
deg <- experiment %>% left_join(annotation_data, by = c('ORF'='LT2')) %>% 
  ungroup() %>%  
  select(SL1344, logFC, t, P.Value) %>%
  dplyr::rename(ID=SL1344) %>% 
  drop_na() %>% 
  dplyr::distinct(ID, .keep_all=TRUE) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

# output:
> head(deg)
           logFC         t  P.Value
SL2237 -3.103376 -33.58213 5.48e-11
SL1384 -2.232289 -27.93897 2.96e-10
SL4312  1.888838  26.71750 4.46e-10
SL4311  2.542126  26.23583 5.27e-10
SL2253  2.096649  21.81188 2.84e-09
SL2252  1.687709  18.39595 1.33e-08
```
We are ready to run the activity analysis using decoupleR. We are using the `ulm` or Univariate Linear Model method, correcting for multiple testing, and filtering the results.
```
contrast_acts <- run_ulm(mat=deg[, 't', drop=FALSE], 
                         net=net, .source='source', .target='target',
                         .mor='mor', minsize = 5)

contrast_acts<-contrast_acts %>% 
  mutate(adjp = p.adjust(p_value, method =c('fdr'))) %>%
  dplyr::filter(adjp <= 0.05)

# results:
> contrast_acts %>%  arrange(score)
# A tibble: 20 × 6
   statistic source condition score  p_value     adjp
   <chr>     <chr>  <chr>     <dbl>    <dbl>    <dbl>
 1 ulm       ompR   t         -6.91 5.65e-12 2.85e-10
 2 ulm       phoP   t         -5.25 1.58e- 7 5.34e- 6
 3 ulm       malT   t         -4.82 1.49e- 6 3.01e- 5
 4 ulm       ttrR   t         -4.42 9.94e- 6 1.25e- 4
 5 ulm       nagC   t         -4.17 3.12e- 5 3.15e- 4
 6 ulm       arcA   t         -4.12 3.86e- 5 3.55e- 4
 7 ulm       yhcK   t         -3.53 4.17e- 4 3.01e- 3
 8 ulm       narL   t         -3.27 1.07e- 3 6.78e- 3
 9 ulm       gntR   t          2.67 7.52e- 3 3.80e- 2
10 ulm       yehT   t          3.00 2.73e- 3 1.45e- 2
11 ulm       modE   t          3.08 2.08e- 3 1.17e- 2
12 ulm       lrp    t          3.11 1.88e- 3 1.12e- 2
13 ulm       csgD   t          3.43 6.19e- 4 4.17e- 3
14 ulm       mraZ   t          3.76 1.73e- 4 1.34e- 3
15 ulm       mntR   t          3.91 9.27e- 5 7.81e- 4
16 ulm       narP   t          4.35 1.37e- 5 1.54e- 4
17 ulm       tdcA   t          4.47 8.15e- 6 1.18e- 4
18 ulm       fnr    t          4.67 3.04e- 6 5.12e- 5
19 ulm       dcuR   t          4.87 1.14e- 6 2.88e- 5
20 ulm       crp    t         12.5  1.89e-35 1.91e-33
```
