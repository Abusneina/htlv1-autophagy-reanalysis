import pandas as pd, pyreadr
g = pyreadr.read_r('/tmp/grch38.rda')['grch38']
g = g[g['biotype'] == 'protein_coding'].dropna(subset=['symbol','chr','start','end','strand'])
g = g[g['chr'].astype(str).str.match(r'^(\d+|X|Y)$')]
g['chrom'] = 'chr' + g['chr'].astype(str)
g['start'] = g['start'].astype(int); g['end'] = g['end'].astype(int)
g['tss'] = g.apply(lambda r: r['start'] if r['strand'] > 0 else r['end'], axis=1)
g = g.sort_values('symbol').drop_duplicates('symbol')          # one locus per symbol
g[['symbol','chrom','start','end','tss','strand']].to_csv('regions_hg38.tsv', sep='\t', index=False)
print(g.shape[0], 'protein-coding loci with hg38 coordinates')
mods = pd.read_csv('code/modules_frozen_v1.tsv', sep='\t')
have = set(g['symbol'])
print('module genes located:', sum(x in have for x in mods.gene), '/', len(mods))
print('missing:', [x for x in mods.gene if x not in have])
