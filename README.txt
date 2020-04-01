HeterosisAI Project

1. SLiM simulation templates 
- Directory: slim/
- model0_neg.txt: Model_0 recessive deleterious model (adapted from Model3 in Kim et al. 2018 Plos Genetics)
- model0_neu.txt: Model_0 neutral model
- modelh_neg.txt: Model_h recessive deleterious model (human demography)
- modelh_neu.txt: Model_h neutral model 


2. Summary statistics calculation
- Directory: python/
- run_slim_get_stats.py: Python script that generates SLiM scripts based on parameters specified in command line, runs simulations, calculates summary statistics in 50kb windows, and writes output file for plotting

3. R plots (False positive rates, ecdf etc.)
- Directory: rplot/
- plot-FPR+ecdf.R: R script that reads simulation outputs (sample output for HYAL2 gene available under "stats/"), computes the critical values for each summary statistic, and plots FPR and ecdf for data visualization

4. Genomic structures for adaptive introgression candidates in human (used in SLiM)
- Directory: regions/
- Content: 26 candidate regions (5MB segment) + "chr11max" segment with the highest exon density in human genome