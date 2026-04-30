# MRMU Documentation

Welcome to the MRMU project documentation.

## Tutorial

- [MRMU R Package Tutorial](MRMU-Rpackage-Tutorial.md)
- [MRMU vs MRAPSS Comparison](MRMU-vs-MRAPSS-Comparison.md)

## Data Availability

The tutorial uses example GWAS summary statistics, a 1000 Genomes EUR PLINK
reference panel, and EUR LD score files. Download links will be added here after
the shared data links are available.

Expected folder structure after downloading data:

```text
MRMU_R_github/
  raw_gwas/
    Biomarker_Urate
    Metabolic_SBP
    Metabolic_DBP
    CAD_UKB
  plink_reference/
    all_1000G_EUR_Phase3.bed
    all_1000G_EUR_Phase3.bim
    all_1000G_EUR_Phase3.fam
  ldscore/
    eur_w_ld_chr/
      1.l2.ldscore.gz
      ...
      22.l2.ldscore.gz
      1.l2.M_5_50
      ...
      22.l2.M_5_50
  tools/
    plink_mac/
      plink
```

## GitHub Repository

Repository: [hxh0928/MRMU](https://github.com/hxh0928/MRMU)
