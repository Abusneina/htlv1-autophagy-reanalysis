# Third-party data included or required

## Included in this repository

**MSigDB Hallmark collection** (`h.all.v7.0.symbols.gmt`)
Molecular Signatures Database, Broad Institute, MIT and the Regents of the
University of California. MSigDB v6.0 and above are distributed under a Creative
Commons style license. Please cite Liberzon A, Birger C, Thorvaldsdottir H,
Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database hallmark gene
set collection. Cell Syst. 2015;1:417-425.

**MCP-counter signatures** (`mcpcounter_signatures.txt`)
From the repository of Becht et al., Genome Biology 2016;17:218.

## Not included, download separately

The Reactome, KEGG legacy and KEGG MEDICUS collections are required to reproduce
the full 1,635-set index but are not redistributed here.

- **Reactome**: available from https://www.gsea-msigdb.org after free registration.
  Reactome content is CC BY 4.0 and could be redistributed, but it is omitted here
  so that every pathway collection is obtained from a single, versioned source and
  users are not left mixing a stale copy with current downloads.
- **KEGG MEDICUS**: available from the same source. KEGG MEDICUS derived sets are copyright Kanehisa
  Laboratories and are provided by MSigDB under a Creative Commons
  Attribution-ShareAlike 4.0 license; they are excluded here so that the ShareAlike
  obligation does not attach to this MIT-licensed repository.
- **KEGG legacy**: copyright 1995-2017 Kanehisa Laboratories, provided to MSigDB
  under a license granted to the Broad Institute for inclusion in that release.
  That permission does not extend to redistribution by third parties, so these
  sets must be obtained from MSigDB directly. See http://www.kegg.jp/kegg/docs/plea.html
  for the copyright holder's academic and commercial licensing policies.

To reproduce the full index, download the symbol GMT files for these collections
and place them in the repository root before running `csi_metric.py`, which reads
every `.gmt` file it finds there. The published index values for all 1,635 sets,
which contain set names and statistics but no gene memberships, are in
`results/results_CSI_index.csv`.
