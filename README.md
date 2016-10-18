Author     Chloe-Agathe Azencott
Contact    chloe-agathe.azencott@mines-paristech.fr

These Python scripts are provided for basic tasks around the usage of the GLIDE package for the detection of epistasis in GWAS data. For usage, please refer to HowToGLIDE_public.pdf (which also contains useful command line instructions).

They were developed at the MLCB (Machine Learning for Computational Biology) research group of the Max Planck Institutes in Tuebingen (Germany).

FILES
	PREPROCESSING:
	qqplot.py
		Plot a Q-Q plot for (single SNP) p-values, from PLINK output.

	naive_impute.py
		Impute missing values in GLIDE input file,
		using the most frequent value for this phenotype.	
	plink2h5.py
		Convert binary PLINK files into h5 file.
	transpose.py
		Transpose PLINK .raw file into .glideIn file.

	split_glide.py
		Write the batch file that will run GLIDE sequentially
		on the "tiles" created beforehand by splitting the input file.

	POSTPROCESSING:
	compute_all_pvalues.py
		Compute all p-values corresponding to 
		the t-test scorse output by GLIDE:
		Run compute_pvalues.py on all outputs 
		and match indices back to SNP names.
	compute_pvalues.py
		Compute the p-values corresponding to 
		the t-test scorse output by GLIDE.
	compute_block_to_index.py
		Compute the dictionary (stored in CONST.py)
		used to convert between tile IDs and SNP indices.
	CONST.py
		Provide the dictionary necessary to the conversion 
		between tile IDs and indices.

	rerun_confounder_h5.py
		Run linear regression again on a small number of SNP pairs.
		Include confounders.
		Read genotypes from .h5 file (created with plink2h5.py).
	rerun_confounder.py
		Run linear regression again on a small number of SNP pairs.
		Include confounders.
	rerun_h5.py
		Run linear regression again on a small number of SNP pairs.
		Read genotypes from .h5 file (created with plink2h5.py).
	rerun.py
		Run linear regression again on a small number of SNP pairs.
	utils.py
		Utilities that can be reused by several other scripts.
	
	meta_analysis.py
		An example of running a meta-analysis of the pairs of SNPs
		selected by running GLIDE on 4 different cohorts.
	
	HowToGLIDE.pdf
		Documentation and additional command lines / shell scripts.
	
