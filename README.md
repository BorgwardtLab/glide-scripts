# Glide companion scripts & tutorial

## Summary 
These Python scripts are provided for basic tasks around the usage of the [GLIDE package](https://github.com/BorgwardtLab/GLIDE) for the detection of epistasis in GWAS data. 

> T. Kam-Thong, C.-A. Azencott, L. Cayton, B. Puetz, A. Altmann, N. Karbalai, P. G. Saemann, B. Schoelkopf, B. Mueller-Myhsok, and K. M. Borgwardt. (2012) **GLIDE: GPU-Based Linear Regression for Detection of Epistasis**. _Human Heredity_, 73 (4), 220-236 [doi: 10.1159/000341885](http://www.karger.com/Article/FullText/341885)

They were mostly developed at the MLCB (Machine Learning for Computational Biology) research group of the Max Planck Institutes in Tuebingen (Germany). Thanks to Damian Roqueiro for stress-testing them and providing useful feedback!

## Documentation
The documentation, which contains additional CLI instructions and short scripts, is available under 
`documentation/HowToGLIDE.pdf`

Compile the documentation with 
`pdflatex -shell-escape documentation/HowToGLIDE.tex`

## Preprocessing scripts
* `qqplot.py`: Plot a Q-Q plot for (single SNP) p-values, from PLINK output.
* `naive_impute.py`: Impute missing values in GLIDE input file, using the most frequent value for this phenotype.	
* `plink2h5.py`: Convert a binary PLINK file into a `.h5` file.
* `transpose.py`: Transpose a PLINK `.raw` file into a `.glideIn` file.
* `split_glide.py`: Write the batch file that will run GLIDE sequentially on the "tiles" created beforehand by splitting the input file.

## Postprocessing scripts
* `compute_all_pvalues.py`:  Compute all p-values corresponding to the t-test scorse output by GLIDE. Run `compute_pvalues.py` on all outputs and match indices back to SNP names.
* `compute_pvalues.py`: Compute the p-values corresponding to the t-test scorse output by GLIDE.
* `compute_block_to_index.py`: Compute the dictionary (stored in `CONST.py`) used to convert between tile IDs and SNP indices.
* `CONST.py`: Provide the dictionary necessary to the conversion between tile IDs and indices.
* `rerun_confounder_h5.py`: Run linear regression again on a small number of SNP pairs. Include confounders. Read genotypes from `.h5` file (created with `plink2h5.py`).
* `rerun_confounder.py`: Run linear regression again on a small number of SNP pairs. Include confounders.
* `rerun_h5.py`: Run linear regression again on a small number of SNP pairs. Read genotypes from `.h5` file (created with `plink2h5.py`).
* `rerun.py`: Run linear regression again on a small number of SNP pairs.
* `utils.py`: Utilities that can be reused by several other scripts.
* `meta_analysis.py`: An example of running a meta-analysis of the pairs of SNPs selected by running GLIDE on 4 different cohorts.

## Contact
Any questions can be directed to Chloe-Agathe Azencott: `chloe-agathe [dot] azencott [at] mines-paristech [dot] fr`
