#!/usr/bin/env python3
"""Compare two damage-annotate per-read tables row by row, keyed by query_id.

usage: compare_annotated.py GOT EXPECTED   (either may be gzipped)
Prints row counts, header equality, reads missing from either side and per-column mismatch
counts. Exits 0 only when the tables are identical per read.
"""
import csv
import gzip
import sys


def load(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        rows = csv.reader(f, delimiter='\t')
        header = next(rows)
        table = {}
        for row in rows:
            if row[0] in table:
                raise SystemExit(f'{path}: duplicate query_id {row[0]}')
            table[row[0]] = row
    return header, table


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    got_h, got = load(sys.argv[1])
    exp_h, exp = load(sys.argv[2])
    bad = 0
    print(f'rows: got={len(got)} expected={len(exp)}')
    print(f'header identical: {got_h == exp_h}')
    if got_h != exp_h:
        bad += 1
        print('  got:      ' + '\t'.join(got_h))
        print('  expected: ' + '\t'.join(exp_h))
    missing = sorted(exp.keys() - got.keys())
    extra = sorted(got.keys() - exp.keys())
    print(f'reads only in expected: {len(missing)}')
    print(f'reads only in got: {len(extra)}')
    for q in missing[:5]:
        print('  missing: ' + q)
    for q in extra[:5]:
        print('  extra: ' + q)
    bad += len(missing) + len(extra)
    per_col = {}
    for q in got.keys() & exp.keys():
        for name, a, b in zip(exp_h, got[q], exp[q]):
            if a != b:
                per_col[name] = per_col.get(name, 0) + 1
    print(f'columns with differences: {len(per_col)}')
    for name, n in sorted(per_col.items(), key=lambda t: -t[1]):
        print(f'  {name}\t{n}')
    bad += sum(per_col.values())
    print('RESULT: ' + ('IDENTICAL' if bad == 0 else f'{bad} DIFFERENCES'))
    return 0 if bad == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
