"""Write hg19 (GRCh37) versions of the hg38 precomputed score files.

    python make_hg19.py hg38ToHg19.over.chain.gz hg19.2bit OUTDIR file_hg38.csv.gz [...]

Every variant is lifted with the UCSC hg38ToHg19 chain (pyliftover). A variant is kept only
if its position maps to exactly one place on chr1-22, X or Y and the hg19 reference base
equals its REF allele (complemented when the chain maps it to the minus strand; ALT is
complemented too). Variants that fail are listed with the reason in
hg19_liftover_excluded.tsv; hg19 variants reached from two different hg38 variants are
excluded as well. Each output line is the hg19 chr, pos, ref, alt, then the hg38 chr, pos,
ref, alt (as chr_hg38, pos_hg38, ref_hg38, alt_hg38), then the rest of the input line
copied verbatim: the scores were computed on hg38 sequences and only the coordinates
change. Output is sorted by hg19 position and written as <name>_hg19.csv.gz (gzip with
mtime 0, so reruns are byte-identical). Requires pyliftover and py2bit.
"""
import gzip
import os
import sys
from collections import Counter

import py2bit
from pyliftover import LiftOver

chain, twobit, outdir, files = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
lo = LiftOver(chain)
genome = py2bit.open(twobit)
PRIMARY = [f'chr{c}' for c in list(range(1, 23)) + ['X', 'Y']]
RANK = {c: i for i, c in enumerate(PRIMARY)}
COMP = str.maketrans('ACGT', 'TGCA')
os.makedirs(outdir, exist_ok=True)
excluded = []

for path in files:
    name = os.path.basename(path).replace('.csv.gz', '')
    stats, kept, n = Counter(), [], 0
    with gzip.open(path, 'rt', newline='') as fh:
        header = fh.readline().rstrip('\n')
        assert header.split(',')[:4] == ['chr', 'pos', 'ref', 'alt'], header[:60]
        rest_header = header.split(',', 4)[4]
        for line in fh:
            n += 1
            c, p, r, a, rest = line.rstrip('\n').split(',', 4)
            hits = lo.convert_coordinate(c, int(p) - 1)
            reason = None
            if not hits:
                reason = 'unmapped'
            elif len(hits) > 1:
                reason = 'multiple_hg19_positions'
            else:
                c19, p0, strand, _ = hits[0]
                if c19 not in RANK:
                    reason = f'non_primary_contig:{c19}'
                else:
                    r19, a19 = (r, a) if strand == '+' else (r.translate(COMP), a.translate(COMP))
                    base = genome.sequence(c19, p0, p0 + 1).upper()
                    if base != r19:
                        reason = f'hg19_reference_is_{base}'
            if reason:
                stats[reason.split(':')[0]] += 1
                excluded.append((name, c, p, r, a, reason))
                continue
            stats['minus_strand' if strand == '-' else 'plus_strand'] += 1
            kept.append((RANK[c19], p0 + 1, r19, a19, f'{c19},{p0 + 1},{r19},{a19},{c},{p},{r},{a},{rest}'))
    keys = Counter(k[:4] for k in kept)
    shared = [k for k in kept if keys[k[:4]] > 1]
    for k in shared:
        c, p, r, a = k[4].split(',', 8)[4:8]
        excluded.append((name, c, p, r, a, 'hg19_variant_shared'))
    stats['hg19_variant_shared'] = len(shared)
    kept = [k for k in kept if keys[k[:4]] == 1]
    kept.sort(key=lambda k: k[:4])
    target = os.path.join(outdir, f'{name}_hg19.csv.gz')
    with open(target, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=6) as gz:
        gz.write(f'chr,pos,ref,alt,chr_hg38,pos_hg38,ref_hg38,alt_hg38,{rest_header}\n'.encode())
        for i in range(0, len(kept), 100000):
            gz.write(''.join(k[4] + '\n' for k in kept[i:i + 100000]).encode())
    print(f'{name}: {n:,} hg38 variants -> {len(kept):,} hg19 variants; ' +
          ', '.join(f'{k} {v:,}' for k, v in sorted(stats.items())), flush=True)
    del kept, keys, shared

with open(os.path.join(outdir, 'hg19_liftover_excluded.tsv'), 'w') as fh:
    fh.write('file\tchr_hg38\tpos_hg38\tref_hg38\talt_hg38\treason\n')
    for row in excluded:
        fh.write('\t'.join(row) + '\n')
print(f'excluded {len(excluded):,} variants in total', flush=True)
