# Workflow

An outline of the scripts used to identify CR1s in avian genomes and determine insertion timing as utilized in Galbraith et al., 2021

## Initial identification of CR1s

 - Repetitive sequences identified using CARP (Zeng et al. 2018) ran with default parameters

 - lineFinder.R used to identify LINEs by searching for conserved reverse transcriptase and endonuclease domains using RPSBLAST (Altschul et al. 1997) with the CDD database (Marchler-Bauer et al. 2017). Full length LINEs aligned with MAFFT (Katoh and Standley 2013) for manual curation.
 
 - BLASTN (Zhang et al. 2000) searches against Repbase (Bao et al. 2015) used to identify and LINEs which were not CR1s
 
## Reciprocal identification of CR1s

- initialPARepeatFinder used BLASTN (Zhang et al. 2000) to identify divergent CR1s using iterative BLASTN searches with CR1s identified above and avian CR1s extracted from Repbase (Bao et al., 2015). After each search full length CR1s were identified using the method in lineFinder.R and redundancy reduced by clustering using VSEARCH (Rognes et al. 2016)

## Insertion event timing

- best_small_hit_finder.R used to identify CR1s 100-600bp long and extend flanks for searches in related species

- Following BLASTN searches of related species from 100-600bp CR1s, species specific versions of smallOrthologueFinder.R (found in small_species_specific_ortho folder) were used to determine presence/absence of CR1s insertions in said related species

## References
Galbraith JD, Kortschak RD, Suh A, Adelson DL, 2021. Genome stability is in the eye of the beholder: recent retrotransposon activity varies significantly across avian diversity. preprint available at https://www.biorxiv.org/content/10.1101/2021.04.13.439746v1

Zeng L, Kortschak RD, Raison JM, Bertozzi T, Adelson DL. 2018. Superior ab initio identification, annotation and characterisation of TEs and segmental duplications from genome assemblies. PLoS One 13:e0193588.

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res. 25:3389–3402.

Marchler-Bauer A, Bo Y, Han L, He J, Lanczycki CJ, Lu S, Chitsaz F, Derbyshire MK, Geer RC, Gonzales NR, et al. 2017. CDD/SPARCLE: functional classification of proteins via subfamily domain architectures. Nucleic Acids Res. 45:D200–D203.

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol. Biol. Evol. 30:772–780.

Zhang Z, Schwartz S, Wagner L, Miller W. 2000. A greedy algorithm for aligning DNA sequences. J. Comput. Biol. 7:203–214.

Bao W, Kojima KK, Kohany O. 2015. Repbase Update, a database of repetitive elements in eukaryotic genomes. Mob. DNA 6:11.

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
