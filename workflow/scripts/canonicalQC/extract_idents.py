#! /usr/bin/env python
import pandas as pd

inames  = pd.read_csv(
    snakemake.input.inames_tsv,
    sep='\t', header=None, names=['ident.name', 'id']
)
sel = list(inames[inames['id'] == snakemake.wildcards.lvl]['ident.name'])[0]

idents = pd.read_csv(snakemake.input.idents_tsv, sep='\t')
idents[sel].to_csv(snakemake.output[0], sep='\t', header=False)
