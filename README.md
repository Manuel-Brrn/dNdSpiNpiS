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

## 1. General Columns

### `Contig_name`
- **Description**: The name of the contig (sequence) being analyzed.

### `Status`
- **Description**: Indicates whether the contig passed or failed the filtering criteria (e.g., gaps, frameshifts, etc.).

---

## 2. SNP-Related Columns

### `BiAllelic_SNP`
- **Description**: Number of bi-allelic single nucleotide polymorphisms (SNPs) in the contig.
- **Calculation**: Counts SNPs where exactly two alleles are present in the population. A SNP is a variation at a single position in the DNA sequence among individuals.
At a specific position in the genome, one individual might have an A, while another individual might have a G.
An allele is one of the possible variants at a specific position in the genome, at a SNP site, the possible alleles could be A, T, C, or G.
When we say "two alleles are present in the population", it means that at that specific SNP site:
Only two different nucleotides (e.g., A and G) are observed across all individuals in the population.

### `TriAllelic_SNP`
- **Description**: Number of tri-allelic SNPs in the contig.
- **Calculation**: Counts SNPs where exactly three alleles are present in the population.

### `QuadriAllelic_SNP`
- **Description**: Number of quadri-allelic SNPs in the contig.
- **Calculation**: Counts SNPs where exactly four alleles are present in the population.

---

## 3. Codon and Fixation Columns

### `nN`
- **Description**: Number of non-synonymous sites in the contig. Non-synonymous sites (nN): Positions in a codon where a mutation changes the amino acid in the protein.
- **Calculation**: Estimated based on the codon sequence and the genetic code.

### `nS`
- **Description**: Number of synonymous sites in the contig. Synonymous sites (nS): Positions in a codon where a mutation does not change the amino acid (silent mutation).
- **Calculation**: Estimated based on the codon sequence and the genetic code.

-  **Codon Analysis**:
  The tool examines each codon in the sequence and determines which positions (1st, 2nd, or 3rd) are synonymous or non-synonymous based on the genetic code.
        For example: In the codon ATG (which codes for Methionine):
                The 3rd position (G) is synonymous because changing it to A, C, or T still codes for Methionine.
                The 1st and 2nd positions are non-synonymous because changes here can alter the amino acid.

    **Counting Sites**:
        For each codon, the tool counts how many of its positions are synonymous (nS) and how many are non-synonymous (nN).
        These counts are summed across all codons in the contig to get the total nN and nS.
   
    **Normalization**:
        The counts are often normalized by the total number of codons to account for differences in sequence length.

### `fixN`
- **Description**: Number of fixed non-synonymous differences between the ingroup and outgroup.
- **Calculation**: Compares the ingroup and outgroup sequences to identify fixed differences.

### `fixS`
- **Description**: Number of fixed synonymous differences between the ingroup and outgroup.
- **Calculation**: Compares the ingroup and outgroup sequences to identify fixed differences.

 ### `Fixed Differences`
    Fixed differences are mutations that are fixed (i.e., present in all individuals of the ingroup but absent in the outgroup, or vice versa).
    These represent genetic changes that have become permanent in one group relative to the other.

**Analysis**
    **Identifying Fixed Differences**:
        The tool compares the ingroup sequences (focal species) to the outgroup sequences (external species) at each codon position.
        For each codon:
            If the ingroup has a different amino acid than the outgroup, and this difference is consistent across all ingroup individuals, it is counted as a fixed non-synonymous difference (fixN).
            If the ingroup has a different nucleotide than the outgroup, but the amino acid remains the same, it is counted as a fixed synonymous difference (fixS).
    **Counting Differences**:
        The tool counts the total number of fixed non-synonymous (fixN) and synonymous (fixS) differences across all codons.

**Example**: 
Input Data
    Ingroup Sequence: ATG CTA GGA (Methionine, Leucine, Glycine)
    Outgroup Sequence: ATG CTG GGA (Methionine, Leucine, Glycine)

Step 1: Calculate nN and nS
    For each codon:
        ATG:
            Non-synonymous sites: 1st and 2nd positions.
            Synonymous sites: 3rd position.
        CTA:
            Non-synonymous sites: 1st and 2nd positions.
            Synonymous sites: 3rd position.
        GGA:
            Non-synonymous sites: 1st and 2nd positions.
            Synonymous sites: 3rd position.
    Total nN = 6 (2 per codon × 3 codons).
    Total nS = 3 (1 per codon × 3 codons).

Step 2: Calculate fixN and fixS
    Compare ingroup and outgroup sequences:
        ATG vs. ATG: No difference.
        CTA vs. CTG: Difference at the 3rd position (A → G).
            This is a synonymous difference because both codons code for Leucine.
        GGA vs. GGA: No difference.
    Total fixN = 0 (no non-synonymous differences).
    Total fixS = 1 (one synonymous difference).

---

## 4. Polymorphism Columns

### `polymN0.0`
- **Description**: Number of polymorphic non-synonymous sites at a frequency threshold of `0.0` (all polymorphic sites).
- **Calculation**: Counts non-synonymous sites with any level of polymorphism in the ingroup.

### `polymS0.0`
- **Description**: Number of polymorphic synonymous sites at a frequency threshold of `0.0` (all polymorphic sites).
- **Calculation**: Counts synonymous sites with any level of polymorphism in the ingroup.

### `polymN0.2`
- **Description**: Number of polymorphic non-synonymous sites at a frequency threshold of `0.2`.
- **Calculation**: Counts non-synonymous sites where the minor allele frequency (MAF) is ≥ `0.2`.

### `polymS0.2`
- **Description**: Number of polymorphic synonymous sites at a frequency threshold of `0.2`.
- **Calculation**: Counts synonymous sites where the minor allele frequency (MAF) is ≥ `0.2`.

---

## 5. McDonald-Kreitman (MK) Test Columns

### `MK_nN`
- **Description**: Number of non-synonymous sites used in the McDonald-Kreitman test.
- **Calculation**: Derived from the comparison of fixed and polymorphic sites between ingroup and outgroup, sum of non-synonymous polymorphic and fixed sites.

### `MK_nS`
- **Description**: Number of synonymous sites used in the McDonald-Kreitman test.
- **Calculation**: Derived from the comparison of fixed and polymorphic sites between ingroup and outgroup, sum of synonymous polymorphic and fixed sites.

---

## 6. GC Allele Frequency Columns

### `GCalfreq<=0.125` to `GCalfreq<=1`
- **Description**: Number of sites with GC allele frequencies falling into specific bins (e.g., `0-12.5%`, `12.5-25%`, ..., `87.5-100%`).
- **Calculation**: Based on the `-gc_bin` parameter, the tool divides GC frequencies into bins and counts the number of sites in each bin.

---

## 7. Species-Specific Columns

### `Ae_speltoides_nb_sequence`
- **Description**: Number of sequences from the species *Ae. speltoides* in the contig.
- **Calculation**: Counts sequences labeled with the species name *Ae. speltoides*.

### `Ae_speltoides_nb_complete_site`
- **Description**: Number of complete (non-gap, non-ambiguous) sites in *Ae. speltoides* sequences.
- **Calculation**: Counts sites without gaps or `N` in *Ae. speltoides* sequences.

### `Ae_speltoides_GC3`
- **Description**: GC content at the third codon position in *Ae. speltoides* sequences.
- **Calculation**: Calculates the proportion of G and C bases at the third codon position.

### `Ae_speltoides_piN`
- **Description**: Nucleotide diversity (π) at non-synonymous sites in *Ae. speltoides*.
- **Calculation**: Measures genetic diversity based on non-synonymous SNPs.

### `Ae_speltoides_piS`
- **Description**: Nucleotide diversity (π) at synonymous sites in *Ae. speltoides*.
- **Calculation**: Measures genetic diversity based on synonymous SNPs.

### `Ae_speltoides_Fit`
- **Description**: F-statistic for *Ae. speltoides*.
- **Calculation**: Measures genetic differentiation between populations.
The general formula for Fst​ is:
FST=(HT−HS)/HT
where:
    HT​ = Total genetic diversity in the entire population (expected heterozygosity across all populations).
    HS​ = Average genetic diversity within subpopulations.
Interpretation
    FST​=0 → No genetic differentiation (all populations are genetically identical).
    FST close to 1 → Complete genetic differentiation (populations share no alleles).
    FST between 0 and 1 → Some level of differentiation, with higher values indicating stronger genetic structuring.
  
### `Ae_speltoides_WeirCockerham84_Fit`
- **Description**: F-statistic calculated using the Weir and Cockerham (1984) method for *Ae. speltoides*.
- **Calculation**: A more refined method for estimating F<sub>ST</sub>.

---

## 8. 
---




