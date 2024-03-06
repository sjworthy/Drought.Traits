
# AusTraits v4.1.0


This is major/minor release of the AusTraits database.

  - austraits-4.1.0.zip: contains the compiled dataset and detailed of
    structure
  - austraits-4.1.0.rds: contains a version of the dataset for direct
    loading in R
  - source code v4.1.0.zip: contains the source materials used to build
    the compiled dataset

For details on access, structure and usage please visit
<https://doi.org/10.5281/zenodo.3568417>

This release was generated from source materials available at
<https://github.com/traitecoevo/austraits.build/releases/tag/v4.1.0> A
full set of changes in the source can be viewed at:
<https://github.com/traitecoevo/austraits.build/compare/v4.0.0>…v4.1.0

Compared to the last version, this release contains substantial
additions of new data and improvement of old data.

| version | dataset\_id |  taxa | locations | traits | records |
| :------ | ----------: | ----: | --------: | -----: | ------: |
| 4.0.0   |         296 | 34028 |      2697 |    470 | 1257443 |
| 4.1.0   |         296 | 34017 |      2697 |    464 | 1253250 |

This release contains:

  - Align units with UCUM standards
  - Correct small mistakes in plant woodiness and life history in datasets extracted from national and state floras (ABS_2022; WAH_2022_1; WAH_2022_2; NHNSW_2022; RBGV_2022; SAH_2022; NTH_2022)
  - Remove a small number of duplicate taxa from taxon.csv
  - Fix mistake in process.R script that was duplicating methods for some dataset by trait combinations
  - Further standardise trait names
  
