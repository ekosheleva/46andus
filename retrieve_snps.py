#Retrieves genotypes at relevant SNPs, given filename of 23andMe genotype raw data, and list of snps

import pandas as pd
import numpy as np

def retrieve_snps(f, snps):
    try:
        df = pd.read_csv(f, sep='\t', comment='#', error_bad_lines=False, warn_bad_lines=False) 
        df.columns = ["rsid", "chromosome", "position", "genotype"]
        df.index = df.rsid      

        snps.index = snps.ID
        geno = pd.merge(snps, df, how='left', left_index = True, right_index=True)
        geno = geno[["genotype"]]

        geno["ref_alleles"] = (geno["genotype"].str.get(0) == snps["Ref allele"]).astype(int) + (geno["genotype"].str.get(1) == snps["Ref allele"]).astype(int)
        geno["alt_alleles"] = (geno["genotype"].str.get(0) == snps["Alt allele"]).astype(int) + (geno["genotype"].str.get(1) == snps["Alt allele"]).astype(int)                   
        return geno
    except:
        print 'Error: please double check your file format; only 23andMe raw data accepted.'
    return 0

