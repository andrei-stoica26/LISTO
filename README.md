# LISTO
`LISTO` is an evolving tool for performing comprehensive overlap assessments
on **lists comprising sets of strings** (such as lists of gene sets). 

While `LISTO` has been developed with scRNA-seq data analysis in mind, the 
methodology is fully applicable for the same problem arising in any other 
setting. Therefore, `LISTO` interacts with general R objects (list of vectors 
and data frames) rather than with scRNAseq objects.

The strings can be provided as **character vectors**, or as the 
**rownames of a data frame with a numeric column designated for ranking**. 
In the latter setting, the numeric columns will be used for generating cutoffs, 
p-values corresponding to each cutoff will be computed, and the median of these 
p-values will be taken as the p-value of the overlap

## Installation

To install `LISTO`, run the following R code:

```
devtools::install_github("andrei-stoica26/LISTO")
```

## Description and usage

This section will elaborate on the functionality and usage of `LISTO`. It 
discusses first the overlaps of individual **elements**, then the details of how
the **lists** of elements must be provided as input.

### Elements

Each **element** taking part in an individual overlap assessed by `LISTO` is a 
**set of strings**. Each overlap assessment of sets of strings answers the 
question of whether the sets intersect each other to a statistically 
significant extent.

### Number of elements in an overlap

`LISTO` currently supports assessments of overlaps of **two** sets of strings.
Partial support will be provided for assessments of overlaps of **three** 
sets of strings as well.

### Selection of elements in an overlap




### Lists




The elements included in the overlap assessments must be provided as **lists**.

Each individual element included in the list must have one of the following 
two input formats:

- A character vector.
    - In this setting, the elements of the vector will be used in the overlap
    assessments.
- A data frame with a numeric column used for ranking.
    - In this setting, the rownames of the data frame will be used in the
    overlap assessments.
    
The **number** of provided lists must be 2 or 3. The following settings will be
supported by `LISTO`:

- 2 lists comprising sets of strings selected from the same set. 
    - Example: Two lists of marker data frames (cluster markers and 
    experimental condition markers) from the same scRNA-seq dataset.
