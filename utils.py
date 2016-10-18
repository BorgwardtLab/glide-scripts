# utils.py --- Utilities that can be reused by several other scripts.

# contact: chloe-agathe.azencott@mines-paristech.fr

import numpy as np

from scipy import stats as st

import statsmodels.api as sm


def readSnp(snpIdx, glideInF):
    """
    Read the SNP at the given line index in glideInF.
    Return a masked array.
    """
    f = open(glideInF)
    for i, line in enumerate(f):
        if i == snpIdx:
            f.close()
            return np.ma.masked_values(\
                [int(val) if val != 'NA' else -9 for val in line.split()], -9)


def prepare_data(snp1x, snp2x, pheno, confounders=None):
    """
    Removed masked values and stack variables in matrix
    before fitting linear/logistic regression.
    """
    # Create mask
    fullMask = np.logical_or(snp1x.mask, snp2x.mask)
    snp1x.mask = fullMask
    snp2x.mask = fullMask
    pheno = np.ma.masked_array(pheno, mask=fullMask)
    if confounders:
        if len(confounders.shape) > 1:
            numCfdrs = confounders.shape[1]
            cfdrMask = np.repeat(fullMask[:, np.newaxis], numCfdrs, 1)
        else:
            cfdrMask = fullMask
    confounders = np.ma.masked_array(confounders, mask=cfdrMask)
    
    # Apply mask
    pheno = np.ma.compressed(pheno)
    snp1x = np.ma.compressed(snp1x)
    snp2x = np.ma.compressed(snp2x)
    confounders = np.ma.compressed(confounders)
    
    # Gather variables in X
    X = np.vstack((np.ones(snp1x.shape), snp1x))
    X = np.vstack((X, snp2x))
    if confounders:
        if len(confounders.shape) > 1:
            X = X.transpose()
            X = np.hstack((X, confounders))
            X = np.hstack((X, snp1x*snp2x))   
        else:
            X = np.vstack((X, confounders))
            X = np.vstack((X, snp1x*snp2x))   
            X = X.transpose()
    else:
        X = np.vstack((X, snp1x*snp2x))   
        X = X.transpose()

    return (X, pheno)


def testPair_logistic(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Use the logistic regression from statsmodels.
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)
    
    # Create model
    model = sm.Logit(pheno, X)
    try:
        # Fit model
        result = model.fit()
        print result.pvalues[-1]
        return (result.pvalues[-1], "\t".join(["%s" % pval for pval in result.pvalues]))
    except Exception:
        sys.stderr.write("Logit fail!\n")
        return (1., "N/A")


def testPair_linear(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Use the linear regression from statsmodels.
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)

    # Create model
    model = sm.OLS(pheno, X)
    try:
        # Fit model
        result = model.fit()
        print result.pvalues[-1]
        return (result.pvalues[-1], "\t".join(["%s" % pval for pval in result.pvalues]))
    except Exception:
        sys.stderr.write("Linear regression fail!\n")
        return (1., "N/A")
    

def testPair_linear_custom(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Custom linear regression implementation,
    as close as possible to that of GLIDE (bar the matrix inversion). 
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)
    m = len(pheno)
    nVars = X.shape[1]
    
    # Fit linear regression
    betas = np.dot(np.linalg.pinv(X), pheno)

    # t-scores
    y = np.dot(X, betas)
    v = np.sum((y-pheno)*(y-pheno)) / (m-4)

    varCovar = np.dot(np.transpose(X), X)
    varCovar = np.linalg.pinv(varCovar) * v
    
    t = betas / np.sqrt(np.diag(varCovar))

    # p-values
    p = [2 * st.t.cdf(-abs(t[i]), m-nVars) for i in range(nVars)]
    p = [pval if not np.isnan(pval) else 0. for pval in p]
    print p[-1]

    # string
    return (p[-1], "\t".join(["%s" % pval for pval in p]))
                                      






def testPair_logistic_intermediate(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Return effect size (beta) and squared standard error for snp1x:snp2x.
    Use the logistic regression from statsmodels.
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)
    
    # Create model
    model = sm.Logit(pheno, X)
    try:
        # Fit model
        result = model.fit()
        return [result.params[-1], np.square(result.bse[-1])]
    except Exception:
        sys.stderr.write("Logit fail!\n")
        return [np.nan, -np.inf]


def testPair_linear_intermediate(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Return effect size (beta) and squared standard error for snp1x:snp2x.
    Use the linear regression from statsmodels.
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)

    # Create model
    model = sm.OLS(pheno, X)
    try:
        # Fit model
        result = model.fit()
        return [result.params[-1], np.square(result.bse[-1])]
    except Exception:
        sys.stderr.write("Linear regression fail!\n")
        return [np.nan, -np.inf]
    

def testPair_linear_custom_intermediate(snp1x, snp2x, pheno, confounders=None):
    """
    Test the given pair of SNPs for association with the phenotype.
    Return effect size (beta) and squared standard error for snp1x:snp2x.
    Custom linear regression implementation,
    as close as possible to that of GLIDE (bar the matrix inversion). 
    """
    (X, pheno) = prepare_data(snp1x, snp2x, pheno, confounders)
    m = len(pheno)
    nVars = X.shape[1]
    
    # Fit linear regression
    betas = np.dot(np.linalg.pinv(X), pheno)

    # Compute standard error
    y = np.dot(X, betas)
    v = np.sum((y-pheno)*(y-pheno)) / (m-4)

    varCovar = np.dot(np.transpose(X), X)
    varCovar = np.linalg.pinv(varCovar) * v
    stderrsq = np.diag(varCovar)

    return [betas[-1], stderrsq[-1]]
    
                                      
