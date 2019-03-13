# Detection of Novel Exons

It is clear from the proteomic diversity and functional complexity of higher order organisms 
the outdated model that one mRNA codes for one protein is implausible [REF]. 
Alternative splicing allows for a single mRNA to code for a variety of protein. 
Currently, there are a handful of computational tools [REF] available to detect novel splice junctions, 
but a common limitation among some of these tools includes the uncertainty of 
parsing biologically relevant splice junctions against experimental error. 
Additionally, there are very few methods that address investigating alternative splicing on long reads.

Here we propose a method discover novel splice junctions, both on short and long reads, 
and validate them against known exons in GENCODE and our own annotated data.

## Getting Started

Before use, it is recommended to run your dataset through metadata scripts (? and if so add link here) to ensure it is complete and clean enough for analysis

### Prerequisites

input: GTF
dependencies: python 3.6.1, pandas, numpy, csv

```
add pipeline information here
```
