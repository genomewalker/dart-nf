#!/usr/bin/env python3
"""Write report.md for one dart-nf sample from the published files. Standard library only."""
import argparse
import collections
import csv
import gzip
import re

ORF_SUFFIX = re.compile(r'_[+-]_\d+_\d+$')   # predict ORF id = <read>_<strand>_<frame>_<ordinal>
REJECT_BITS = {1: 'chance hit (decoy null)', 2: 'low EM margin', 4: 'high effective reference count'}


def opener(path):
    return gzip.open(path, 'rt', encoding='utf-8') if path.endswith('.gz') else open(path, 'rt', encoding='utf-8')


def count_fastx(path):
    reads = 0
    with opener(path) as f:
        first = f.readline()
        if not first:
            return 0
        if first.startswith('@'):
            lines = 1 + sum(1 for _ in f)
            return lines // 4
        reads = 1
        for line in f:
            if line.startswith('>'):
                reads += 1
    return reads


def count_lines(path, skip_prefix=None):
    n = 0
    with opener(path) as f:
        for line in f:
            if skip_prefix and line.startswith(skip_prefix):
                continue
            n += 1
    return n


def hits_stats(path):
    rows = 0
    queries = set()
    with opener(path) as f:
        for line in f:
            rows += 1
            queries.add(line.split('\t', 1)[0])
    reads = {ORF_SUFFIX.sub('', q) for q in queries}
    return rows, len(queries), len(reads)


def read_tsv(path):
    with opener(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def fnum(s):
    try:
        return float(s)
    except ValueError:
        return None


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def counter_table(rows, col, label):
    c = collections.Counter(r[col] for r in rows)
    return md_table([label, 'reads'], sorted(c.items(), key=lambda t: -t[1]))


def main():
    ap = argparse.ArgumentParser()
    for a in ('sample', 'reads', 'gff', 'decoys', 'hits', 'decoy_hits', 'annotated', 'gene_summary',
              'protein_summary', 'annotate_log', 'null_info', 'versions'):
        ap.add_argument('--' + a.replace('_', '-'), required=True)
    ap.add_argument('-o', '--out', required=True)
    args = ap.parse_args()

    n_reads = count_fastx(args.reads)
    n_genes = count_lines(args.gff, skip_prefix='#')
    with opener(args.decoys) as f:
        n_decoys = sum(1 for line in f if line.startswith('>'))
    hit_rows, hit_queries, hit_reads = hits_stats(args.hits)
    dhit_rows, dhit_queries, _ = hits_stats(args.decoy_hits)

    ann = read_tsv(args.annotated)
    genes = read_tsv(args.gene_summary)
    prot = read_tsv(args.protein_summary)
    with open(args.annotate_log, encoding='utf-8') as f:
        log = f.read().splitlines()
    with open(args.null_info, encoding='utf-8') as f:
        null_info = f.read().strip()
    with open(args.versions, encoding='utf-8') as f:
        versions = f.read().strip()

    doc = [f'# dart-nf report: {args.sample}', '', '```', versions, '```', '']

    doc += ['## Counts', '', md_table(['item', 'n'], [
        ('reads in', n_reads),
        ('ORFs predicted (GFF rows, one protein each)', n_genes),
        ('decoy peptides (reversed protein of every N-th read)', n_decoys),
        ('MMseqs2 alignments (hits.tsv rows)', hit_rows),
        ('ORFs with at least one hit', hit_queries),
        ('reads with at least one hit', hit_reads),
        ('decoy alignments / decoy peptides with a hit', f'{dhit_rows} / {dhit_queries}'),
        ('rows in the per-read table (reads whose best alignment has a mismatch)', len(ann)),
        ('genes in gene_summary', len(genes)),
        ('proteins in protein_summary', len(prot)),
    ]), '']

    doc += ['## Null model (decoy search)', '', 'Output of `dart damage-info` on the `.null` sidecar:', '',
            '```', null_info, '```', '']
    null_lines = []
    for i, l in enumerate(log):
        if l.startswith('Null model:'):
            null_lines.append(l)
        elif l.startswith('Warning: no null bin'):   # multi-line warning: continuation lines are indented
            block = [l]
            for cont in log[i + 1:]:
                if cont.startswith('  '):
                    block.append(cont.strip())
                else:
                    break
            null_lines.append(' '.join(block))
    reject = []
    total_line = ''
    rule, threshold = 'null_p', ''
    for l in log:
        m = re.match(r'null reject: bin (\S+) reads=(\d+) flagged=(\d+) \(([\d.]+)%\)(?: chance=(\S+))? '
                     r'bits at (null_p|null_fdr)=([\d.e-]+): (\S+)', l)
        if m:
            rule, threshold = m.group(6), m.group(7)
            reject.append((m.group(1), m.group(2), m.group(3), m.group(4) + '%', m.group(5) or '', m.group(8)))
        elif l.startswith('null reject:'):
            total_line = l
    for l in null_lines:
        doc.append(l)
    doc.append('')
    if rule == 'null_fdr':
        doc.append(f'Reads flagged as chance hits at a per-bin target-decoy FDR of {threshold} (`null reject:` lines of '
                   'damage-annotate). `expected chance` is the number of chance best hits expected in the bin among all '
                   'queries searched (queries times the share of decoys whose best hit fell in the bin). The last column '
                   'is the lowest kept bit score: reads below it are flagged; NA means nothing in the bin is flagged.')
        doc.append('')
        doc.append(md_table(['aligned length (aa)', 'reads with best hit', 'flagged', 'share', 'expected chance',
                             'lowest kept bits'], reject))
    else:
        doc.append(f'Reads flagged by `null_p > {threshold}` per aligned-length bin (`null reject:` lines of damage-annotate). '
                   'The last column is the bit score below which a read in that bin is flagged; NA means the bin holds too few decoys to flag anything.')
        doc.append('')
        doc.append(md_table(['aligned length (aa)', 'reads with best hit', 'flagged', 'share', 'bits at null_p threshold'],
                            [r[:4] + r[5:] for r in reject]))
    if total_line:
        doc += ['', total_line]
    doc.append('')
    qs = [fnum(r.get('null_q', 'NA')) for r in ann]
    qs = [q for q in qs if q is not None]
    if qs:
        doc.append('Other FDR cutoffs from the same run: rows of the per-read table (reads with a mismatch) whose '
                   '`null_q` exceeds each cutoff. Filter `null_q > x` on the table instead of rerunning.')
        doc.append('')
        doc.append(md_table(['FDR cutoff', 'flagged', 'kept', 'share flagged'],
                            [(c, sum(q > c for q in qs), sum(q <= c for q in qs), f'{100 * sum(q > c for q in qs) / len(qs):.2f}%')
                             for c in (0.01, 0.05, 0.1, 0.2)]))
        doc.append('')

    doc += ['## EM assignment', '']
    keys = ('Alignments:', 'Reads:', 'References:', 'Converged in', 'Unique mappers', 'Multi-mappers',
            'EM keep filter', 'Reassigned vs bitscore-best', 'Mean EM margin', 'Mean EM effective refs')
    seen = set()
    em_lines = []
    for l in log:
        s = l.strip()
        for k in keys:
            if s.startswith(k) and k not in seen:
                seen.add(k)
                em_lines.append(s)
    doc += ['```'] + em_lines + ['```', '']
    doc.append(counter_table(ann, 'em_keep', 'em_keep'))
    doc.append('')
    doc.append(counter_table(ann, 'em_ambiguous', 'em_ambiguous'))
    doc.append('')
    rr = collections.Counter(r['em_reject_reason'] for r in ann)
    rows = []
    for v, n in sorted(rr.items(), key=lambda t: -t[1]):
        bits = [name for b, name in REJECT_BITS.items() if v.isdigit() and int(v) & b]
        rows.append((v, ', '.join(bits) if bits else 'kept', n))
    doc.append(md_table(['em_reject_reason', 'meaning', 'reads'], rows))
    doc.append('')

    doc += ['## Per-read damage', '']
    post = [fnum(r['posterior']) for r in ann]
    post = [p for p in post if p is not None]
    hist = collections.Counter(min(int(p * 10), 9) for p in post)
    doc.append(md_table(['posterior', 'reads'], [(f'{b / 10:.1f}-{(b + 1) / 10:.1f}', hist.get(b, 0)) for b in range(10)]))
    doc.append('')
    doc.append(f'posterior >= 0.5: {sum(1 for p in post if p >= 0.5)} of {len(post)} reads with a posterior')
    doc.append('')
    doc.append(counter_table(ann, 'classification', 'classification'))
    doc.append('')
    doc.append(counter_table(ann, 'is_damaged', 'is_damaged'))
    doc.append('')
    doc.append(counter_table(ann, 'tier', 'evidence tier'))
    doc.append('')
    for k in ('Reads with mismatches:', 'Reads with damage-consistent sites:', 'Total damage-consistent sites:',
              'Estimated q_modern:', 'Damage detectability:', 'Library type from damage index:'):
        for l in log:
            if l.strip().startswith(k):
                doc.append('- ' + l.strip())
                break
    doc.append('')

    doc += ['## Proteins', '']
    state = collections.Counter(r['library_damage_state'] for r in prot)
    doc.append('library_damage_state: ' + ', '.join(f'{k} ({n} proteins)' for k, n in state.items()))
    doc.append('')
    doc.append(md_table(['protein_auth_tier', 'proteins'],
                        sorted(collections.Counter(r['protein_auth_tier'] for r in prot).items(), key=lambda t: -t[1])))
    doc.append('')
    doc.append(md_table(['auth_binom_call', 'proteins'],
                        sorted(collections.Counter(r['auth_binom_call'] for r in prot).items(), key=lambda t: -t[1])))
    doc.append('')
    auth = [r for r in prot if r['protein_auth_tier'] in ('detected', 'advisory') and r['auth_binom_call'] == 'ancient']
    cols = ['protein_id', 'n_reads', 'n_damaged', 'damage_sites', 'p_protein_damaged', 'protein_auth_tier', 'auth_binom_call']
    if auth:
        auth.sort(key=lambda r: (-float(r['auth_binom_p_ancient']), -int(r['n_damaged'])))
        doc.append(f'Top proteins called ancient by the authentication test ({len(auth)} in total), by `auth_binom_p_ancient` then `n_damaged`:')
        doc.append('')
        doc.append(md_table(cols, [[r[c] for c in cols] for r in auth[:10]]))
    else:
        doc.append('No protein is called ancient by the authentication test: no protein has `protein_auth_tier` '
                   'detected or advisory with `auth_binom_call = ancient`. The ranking below is ungated, by '
                   '`p_protein_damaged` then `n_damaged`, and is not an authentication.')
        doc.append('')
        ranked = sorted(prot, key=lambda r: (-(fnum(r['p_protein_damaged']) or 0.0), -int(r['n_damaged'])))
        doc.append(md_table(cols, [[r[c] for c in cols] for r in ranked[:10]]))
    doc.append('')

    doc += ['## Genes', '']
    anc = sum(1 for g in genes if (fnum(g['ancient_frac']) or 0.0) >= 0.5)
    doc.append(f'{len(genes)} genes pass the gene-summary read thresholds; {anc} have `ancient_frac >= 0.5`.')
    doc.append('')

    with open(args.out, 'w', encoding='utf-8') as f:
        f.write('\n'.join(doc) + '\n')


if __name__ == '__main__':
    main()
