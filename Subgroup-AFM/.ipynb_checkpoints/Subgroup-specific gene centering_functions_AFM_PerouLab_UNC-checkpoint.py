# Written by: Aranzazu Fernandez-Martinez (afernan4@ad.unc.edu)
# 
# Copyright (c) 2020-2022 Aranzazu Fernandez-Martinez
# 
# Please cite:
# Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
# Hoadley KA, Shepherd J, Tolaney S, Henry NL, Dang C, Harris L, Berry D,
# Hahn O, Hudis C, Winer E, Partridge A, Perou CM, Carey LA. Survival,
# Pathologic Response, and Genomics in CALGB 40601 (Alliance),
# a Neoadjuvant Phase III Trial of Paclitaxel-Trastuzumab With or Without
# Lapatinib in HER2-Positive Breast Cancer. J Clin Oncol. 2020:JCO2001276.
# Epub 2020/10/24. doi: 10.1200/JCO.20.01276. PubMed PMID: 33095682.
# 
# These functions are adapted from "Zhao, Xi, Einar A. Rødland, Robert
# Tibshirani, and Sylvia Plevritis. "Molecular subtyping for clinically
# defined breast cancer subgroups." Breast Cancer Research 17, no. 1 (2015): 
# 29."
# 


def get_sigma(gene_expresion, groups):
    """Build quantile dataframe for each IHC group.
    
    For each gene it calculates the ECDF function within each IHC group,
    and calls it using the overall median of that gene.
    
    :param gene_expresion: Gene expresion dataframe, i.e. UNC232
    :param groups: IHC label dataframe for each sample on the gene expresion.
    
    :return: Quantile dataframe.
    """
    groups_cols = pd.Series(groups.values.ravel()).dropna().unique()
    res = pd.DataFrame({}, columns=groups_cols, index=gene_expresion.index)
    for name, values in gene_expresion.iterrows():
        for col in groups:
            unique = groups[col].dropna().unique()
            for u in unique:
                samples_from_group = groups.loc[groups[col] == u].index
                subset = values[samples_from_group].dropna()
                res[u][name] = ECDF(subset)(values.median())
    return res


def quantile_centering(expr_matrix, gene_quantile):
    """Do row centering based on the quantile and IHC group.
    
    :param expr_matrix: pandas.DataFrame where row are genes and columns are samples
    :param gene_quantile: pandas.DataFrame or Series containig all the genes in the
      first parameter and the value of the quantile to be used, i.e. .5 if one wants
      to do row centering using the mean.

    :return: Subgroup-specific Centered dataframe.
    """
    res = expr_matrix.copy()
    for name, values in expr_matrix.iterrows():
        q = gene_quantile.loc[name]
        q_value = expr_matrix.loc[name].quantile(q)
        res.loc[name] -= q_value
    return res


