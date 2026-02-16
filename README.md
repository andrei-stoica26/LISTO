# LISTO

`LISTO` is a tool for performing **comprehensive overlap assessments**
on lists comprising sets of strings, such as lists of gene sets. It can 
assess:

- Overlaps of pairs of sets of strings selected from the same
universe.
- Overlaps of pairs of sets of strings selected from different universes.
- Overlaps of triplets of sets of strings selected from the same universe. 

While `LISTO` has been developed with scRNA-seq data analysis in mind, 
the methodology is fully applicable for the same problem arising in any other 
setting. Thus, the implementation of `LISTO` uses general R objects 
(data frames, character vectors), rather than scRNA-seq-specific objects.

## Installation

To install `LISTO`, run the following R code:

```
devtools::install_github("andrei-stoica26/LISTO")
```

## Description and usage

This section will elaborate on the functionality and usage of `LISTO`. It 
discusses first the overlaps of individual **elements**, then the details of how
the **lists** of elements must be provided as input.


### Items

Each **item** taking part in an individual overlap assessed by `LISTO` is a 
**set of strings**. Each overlap assessment of sets of strings answers the 
question of whether the sets intersect each other to a statistically 
significant extent.

### Lists

The `runLISTO` function runs the entire LISTO pipeline. It requires two lists
as input. Each list can store two types of elements:

- Character vectors.
- Data frames with a numeric column specified by the `numCol` parameter. 

A third list, containing the same type of elements, can be optionally provided.

### Extracting items from lists

Items to be used in the overlap assessments are extracted from the input
lists as follows:

- **Character vectors**: They are used as such.

- **Data frames**: The rownames of the data frame are selected, and overlaps
are calculated based on cutoffs determined by the distinct values in the
column specified by `numCol`. The **median** of the resulting p-values 
is taken to be the p-value of the corresponding overlap.

