# Repeat Quinton et al., 1991, with genomic selection

<!-- `echo Repeat Quinton et al., 1991, with genomic selection | md5sum` => effb94... -->

## Summary

BLUP is better in short term selection.
But for long term selection, phenotypic selection is better.
Also, more genes are lost and inbreeding increases faster with BLUP.
Genomic selection (GS) has similar problems, at least in theory.

This simulation is to repeat Quinton et al., 1991.
I also add genomic selection here.

## Simulation details

500 sires and 500 dams are simulated with `MaCS`. `SLiM` can be used later.
They are randomly mated to have 21 sibs in each family to
produce 10,500 offspring.
500 of them are nucleus with half male and female calves.
Sires for the next generation are selected from nucleus.
Dams for the next generation are selected from all the offspring.

Selections are done with phenotypic values and EBV from BLUP and GS.
Inbreeding, $\Delta G$, gene loss are recorded.

## Progress

1. A starter population for selection
   1. [x] Formats of genotype and linkage map files
      - Genotypes
        - 24-byte header
        - genotype matrix
      - Linkage map
        - A DataFrame
          - Chromosome number
          - BP position
          - Allele frequency
        - Serialized with `Serialization`
      - Tree storage is currently quite beyond my capacity
        - The storage
        - SNP coding
   2. [x] Simulate a base population with `MaCS`.
      - Can link to `SLiM` later.
      - 2000 haplotypes.
      - chromosome lengths are from NCBI
      - ~19M SNP
      - Genotype and linkage map files
   3. [x] Sample 50k SNPs, and 1000 pairs of haplotypes from the base population.
   4. [ ] Mate 500 sires and 500 dams to produce 10,500 offspring.
      - Unique coding of the SNP alleles using `UInt16`.
      - odd for SNP allele 1, even for 0.
      - Can code at least $2^{14} = 16,384$ ID, or 32,768 haplotypes.
2. Selection simulation
   1. [ ] Selection with genomic selection, SNP-BLUP
   2. [ ] Selection with BLUP
   3. [ ] Selection with phenotypic values

## ToDo

1. [ ] Calculate LD to verify the sampling.
2. [ ] Make the `mate` function work.

!!! note
    The simulation is still in debugging stage.
    No result yet.
