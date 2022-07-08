#!/usr/bin/env python

import pandas as pd


def load_genotypes():
    df = pd.read_excel("./genotypes.xls", sheet_name="Sheet 1")
    df['name'] = df['name'].apply(lambda x: x.split("|")[0])
    return pd.Series(df['phylo-major result'].values, df.name).to_dict()
