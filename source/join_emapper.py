
import polars as pl
seqmaps = value(oname_seqmaps[1])
emapper_out = f'{emapper_odir}/emapper.emapper.annotations',
emapper_out = pl.read_csv(emapper_out,
                          skip_lines=4,
                          separator='\t',
                          comment_prefix='##')
hash2original = []
for gs in seqmaps.values():
    for seqhash, header in gs:
        header, _ = header.split(' ', 1)
        hash2original.append((seqhash, header))

hash2original = pl.DataFrame(hash2original,
                             schema=['hash', 'original'],
                             orient='row')
emapper_out.join(hash2original, left_on='#query', right_on='hash')
