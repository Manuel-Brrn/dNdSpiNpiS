## dNdSpiNpiS - README

### 1. Overview

dNdSpiNpiS is a command-line tool designed to analyze genetic sequences, particularly focusing on calculating **dN/dS ratios**, **nucleotide diversity (piN, piS)**, and other **population genetics metrics**.

The program takes a **sequence alignment file** (in Fasta/PseudoFasta format) and computes various statistics. 
It allows you to specify an **ingroup** (focal species) and an **optional outgroup** (external species) for comparative analysis.

---

### 2. Usage

#### **Command Structure**
```bash
dNdSpiNpiS -alignment_file=seqfile.pfas -population_file=custom_pop.txt -ingroup=Species_name_1 -outgroup=Species_name_2 -gc_bin=8 -gapN_site=10 -gapN_seq=0.5 -min_nb_codon=10 -kappa=1 -remove_frameshift=no -compute_distances=no -tolerance_zone=20 -allow_internal_bc=no -print_header=yes -contig_usage=all -out=dNdSpiNpiS_output
```

---

### 3. Parameter Explanations

#### **Mandatory Parameters**

`-alignment_file`
- **Description**: Path to the input alignment file (Fasta/PseudoFasta format).
- **Format**: Sequence names must follow the structure:  
  `>Contig_name|Species_name|Individual_name|Allele_number`
- **Example**: `-alignment_file=seqfile.pfas`

`-ingroup`
- **Description**: Name of the ingroup (focal species) to analyze.
- **Example**: `-ingroup=Species_name_1`

---

#### **Optional Parameters**

`-population_file`
- **Description**: Path to a text file listing individuals from the focal species (one per line). If not provided, all individuals from the ingroup are used as a single population.
- **Example**: `-population_file=custom_pop.txt`

`-outgroup`
- **Description**: Name of the outgroup (external species) for comparative analysis.
- **Example**: `-outgroup=Species_name_2`

`-gc_bin`
- **Description**: Number of bins for GC allele frequency calculations. Must be > 0. It helps in categorizing the frequency of GC content (guanine and cytosine bases) across the sequences. Example: If -gc_bin=8, the tool will divide the GC frequencies into 8 bins (e.g., 0-12.5%, 12.5-25%, ..., 87.5-100%).
- **Default**: `8`
- **Example**: `-gc_bin=8`

`-gapN_site`
- **Description**: Specifies the minimum number of ingroup sequences that must not have gaps or undetermined characters (N) at a given site for the site to be retained in the analysis. Ensures that only sites with sufficient data (i.e., minimal missing or ambiguous data) are included in the calculations.
Example: If -gapN_site=10, at least 10 sequences in the ingroup must have valid data (no gaps or N) at a site for it to be considered.
- **Default**: `10`
- **Example**: `-gapN_site=10`

`-gapN_seq`
- **Description**: Sets the maximum proportion of gaps or undetermined characters (N) allowed in a sequence for it to be retained in the analysis. Filters out sequences with too much missing or ambiguous data. Example: If -gapN_seq=0.5, a sequence will be discarded if more than 50% of its sites contain gaps or N.
- **Default**: `0.5`
- **Example**: `-gapN_seq=0.5`

`-min_nb_codon`
- **Description**: Minimum number of complete codons required (>0).
- **Default**: `10`
- **Example**: `-min_nb_codon=10`

`-kappa`
- **Description**: Defines the transition/transversion ratio (Kappa) used in the analysis. Adjusts the weighting of transitions (e.g., A ↔ G or C ↔ T) versus transversions (e.g., A ↔ C or G ↔ T) in genetic distance calculations. Example: If -kappa=2, transitions are considered twice as likely as transversions.
- **Default**: `1`
- **Example**: `-kappa=1`

`-remove_frameshift`
- **Description**: This parameter determines whether the tool should reject or accept contigs that contain frameshift mutations. A frameshift is a genetic mutation caused by insertions or deletions that shift the reading frame of the sequence, often leading to incorrect protein translation.
Options:
    yes: Reject contigs with frameshifts (they will be excluded from the analysis).
    no: Accept contigs with frameshifts (they will be included in the analysis).
Example: If -remove_frameshift=yes, any contig with a frameshift mutation will be discarded.
- **Default**: `no`
- **Example**: `-remove_frameshift=no`

`-compute_distances`
- **Description**: Compute genetic distances (`yes`) or not (`no`). 
 Genetic distances measure how genetically different individuals are from each other, based on their sequences.
Options:
    yes: Compute genetic distances (this will add distance calculations to the output).
    no: Skip genetic distance calculations (this will make the analysis faster if distances are not needed).
- **Default**: `no`
- **Example**: `-compute_distances=no`

`-tolerance_zone`
- **Description**: Defines the number of positions at the beginning and end of a sequence where Stop codons or frameshifts are tolerated. Stop codons or frameshifts in the middle of a sequence are usually problematic, but they can sometimes occur at the edges due to sequencing errors or incomplete data. This parameter allows you to ignore such edge cases.
- **Default**: `20`
- **Example**: `-tolerance_zone=20`

`-allow_internal_bc`
- **Description**: Determines how the tool handles internal "bad codons" (e.g., Stop codons or frameshifts that occur outside the tolerance zone).
What are "bad codons"?: These are codons that disrupt the sequence, such as Stop codons (which prematurely end translation) or frameshifts (which disrupt the reading frame).
Options:
    yes: Tolerate internal bad codons by replacing them with NNN (an undetermined codon). The contig will still be used in the analysis.
    no: Reject the entire contig if it contains internal bad codons.
- **Default**: `no`
- **Example**: `-allow_internal_bc=no`

`-print_header`
- **Description**: Print a header line (`yes` or `no`).
- **Default**: `yes`
- **Example**: `-print_header=yes`

`-contig_usage`
- **Description**: Specifies which contigs to use (`all`, `even`, or `odd`).
- **Default**: `all`
- **Example**: `-contig_usage=all`

`-out`
- **Description**: Name of the output file.
- **Example**: `-out=dNdSpiNpiS_output`

---

#### **Bootstrap-Related Parameters**

`-bs_replicates`
- **Description**: Sets the number of bootstrap resampling replicates to perform. What is bootstrap resampling?: It's a statistical method where the tool randomly resamples the data (with replacement) to estimate confidence intervals and robustness of the results.
- More replicates increase the accuracy of the bootstrap analysis but also increase computation time.
- **Default**: `10000`
- **Example**: `-bs_replicates=10000`

`-bs_sort`
- **Description**: Sets the sorting precision for calculating confidence intervals during bootstrap analysis.
What is sorting precision?: It determines how finely the tool sorts the bootstrap results to calculate confidence intervals. Higher values provide more precise confidence intervals.
- **Default**: `1000000`
- **Example**: `-bs_sort=1000000`

`-bs_confidence`
- **Description**: Sets the confidence level for the bootstrap calculations.
What is a confidence level?: It represents the probability that the true value lies within the calculated confidence interval.
Range: Must be a value between 0 and 1 (e.g., 0.95 means 95% confidence).
- **Default**: `0.95`
- **Example**: `-bs_confidence=0.95`

---

### 4. Output

The tool generates an output file containing:
- **dN/dS ratios**
- **Nucleotide diversity (piN, piS)**
- **Allele frequencies**
- **Bootstrap confidence intervals (if enabled)**

---

