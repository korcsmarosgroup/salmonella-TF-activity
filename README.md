# `SalmonAct` - Transcription Factor Activity Inference for Salmonella Transcriptomics



## Repository Overview

This repository contains scripts, data, and tutorials to accompany the SalmonAct resource, a regulon tailored for inferring transcription factor activities in _Salmonella_. SalmonAct is a signed and directed prior knowledge network (PKN) made for _Salmonella_ transcriptomics, expanding on the SalmoNet2 interaction resource. Developed for compatibility with decoupleR, an R and Python library, this repository enables transcription factor activity estimation, aiding researchers in exploring *Salmonella* biology.

## Key Features:

**SalmonAct PKN**: A curated network of signed and directed interactions involving 191 transcription factors for _Salmonella_.

**Data Compatibility**: Harmonized interaction data from databases like PRODORIC, CollecTF, RegPrecise, SalComRegulon and RegulonDB, formatted for decoupleR.

**Analysis Scripts**: R scripts to replicate differential expression and transcription factor activity analyses shown in the manuscript.

**Mapping Table**: ID mapping table based on protein orthology established for SalmoNet2 to aid analysis in 20 strains

For a detailed tutorial please navigate to the `tutorial` folder, or follow the steps replicating the analysis shared in the manuscript in `scripts/salmonact_notebook.Rmd`.

