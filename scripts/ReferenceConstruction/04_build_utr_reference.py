#!/usr/bin/env python3
"""
build_utr_reference.py
======================
Reconstructs a yeast transcriptomic (mRNA) FASTA reference from:

  (1) a UTR-length table  (per-gene 5' and 3' UTR lengths in nt), and
  (2) a curated gene-model backbone GFF3 (gene + CDS + intron features), and
  (3) the S288C genomic FASTA.

It reproduces, byte-for-byte, the original custom reference
`S288C_reference_sequence_R64-4-1_20230823_w.fa` when driven with the
Final UTR table (final_utr.tsv).  The same machinery is reused to
build alternative references (e.g. the Nakagalakshmi UTR definitions) by
swapping the UTR table.

PIPELINE
--------
  step 1  rename genome chromosomes  (NCBI 'ref|NC_00xxxx|' -> I..XVI, mt)
  step 2  overlay UTR features onto the gene-model backbone:
            - UTRs are written as `CDS` features sharing Parent=GENE_mRNA,
              so gffread treats them as exonic and includes them in the
              spliced transcript while excluding `intron` features.
  step 3  gffread <gff> -g <genome> -w <out.fa> -M -A
            -> spliced transcript (mRNA) FASTA, one record per gene model.

UTR COORDINATE RULE (reverse-engineered & validated against the original
build; see METHODS.md and validate_against_original.py)
-----------------------------------------------------------------------
For a gene on [s, e] (1-based, inclusive) with strand:
  * The LEFT (lower-coordinate) UTR is the 5' UTR on '+' genes and the
    3' UTR on '-' genes.  Span = [s - L, s - 1]  (exactly L nt).
  * The RIGHT (higher-coordinate) UTR is the 3' UTR on '+' genes and the
    5' UTR on '-' genes.  Span = [e + 1, e + R + 1]  (R + 1 nt; the extra
    base is a quirk of the original construction and is preserved here).
  * Low coordinate is clamped to 1; high coordinate is clamped to the
    chromosome length.  A UTR whose final genomic span is < 2 nt, or whose
    table length is 0/NA, is emitted as `NA NA` (gffread skips it).
"""
import argparse
import os
import re
import subprocess
import sys

# ------------------------------------------------------------------ helpers
def attr(field9, key):
    m = re.search(key + r'=([^;]+)', field9)
    return m.group(1) if m else None


def rename_genome(in_fa, out_fa):
    """Rename NCBI 'ref|NC_...| [chromosome=X]' headers to bare 'X'
    (and the mitochondrion to 'mt'); sequence lines are copied verbatim.
    If the genome already uses bare chromosome names it is copied as-is."""
    n = 0
    with open(in_fa) as fi, open(out_fa, 'w') as fo:
        for line in fi:
            if line.startswith('>'):
                h = line[1:].strip()
                if 'location=mitochondrion' in h or 'chromosome=mt' in h.lower():
                    name = 'mt'
                else:
                    m = re.search(r'\[chromosome=([^\]]+)\]', h)
                    if m:
                        name = m.group(1)
                    else:
                        # already a bare/simple id -> keep first token
                        name = h.split()[0]
                fo.write('>' + name + '\n')
                n += 1
            else:
                fo.write(line)
    return n


def load_utr_table(path):
    """Read a 3-column UTR-length table (Gene, five_prime, three_prime).
    Accepts TSV or CSV; 'NA'/''/non-numeric -> 0 (no UTR)."""
    utr = {}
    with open(path) as fh:
        first = fh.readline()
        delim = '\t' if '\t' in first else ','
        # Re-read in case there is no header (detect by trying to parse line1)
        def parse_len(x):
            x = x.strip()
            if x == '' or x.upper() == 'NA':
                return 0
            try:
                return max(0, int(float(x)))
            except ValueError:
                return 0
        # decide whether first line is a header
        parts = first.rstrip('\n').split(delim)
        header_like = not parts[1].strip().lstrip('-').isdigit() if len(parts) > 2 else True
        if not header_like and len(parts) >= 3:
            utr[parts[0].strip()] = (parse_len(parts[1]), parse_len(parts[2]))
        for line in fh:
            p = line.rstrip('\n').split(delim)
            if len(p) < 3:
                continue
            utr[p[0].strip()] = (parse_len(p[1]), parse_len(p[2]))
    return utr


def load_chrlen(genome_fa):
    """Chromosome lengths from a FASTA (no .fai required)."""
    chrlen = {}
    name = None
    ln = 0
    with open(genome_fa) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    chrlen[name] = ln
                name = line[1:].strip().split()[0]
                ln = 0
            else:
                ln += len(line.strip())
    if name is not None:
        chrlen[name] = ln
    return chrlen


def utr_span(strand, s, e, five_len, three_len, chrlen):
    """Return {'five':(a,b) or None, 'three':(a,b) or None} per the rule."""
    out = {'five': None, 'three': None}
    if strand == '+':
        left_key, left_len, right_key, right_len = 'five', five_len, 'three', three_len
    else:
        left_key, left_len, right_key, right_len = 'three', three_len, 'five', five_len
    # LEFT UTR: [s-L, s-1], clamp low to 1
    if left_len > 0:
        a, b = max(1, s - left_len), s - 1
        out[left_key] = (a, b) if (b - a + 1) >= 2 else None
    # RIGHT UTR: [e+1, e+R+1], clamp high to chromosome length
    if right_len > 0:
        a, b = e + 1, min(chrlen, e + right_len + 1)
        out[right_key] = (a, b) if (b - a + 1) >= 2 else None
    return out


# ------------------------------------------------------------- gff overlay
def build_gff(template_gff, utr, chrlen, out_gff):
    """Rewrite the UTR (CDS *_UTR) coordinate fields of the template GFF
    using `utr` lengths; gene/CDS/intron lines are copied verbatim so the
    curated gene-model backbone is preserved exactly. Coordinate fields are
    formatted as the original (' %d' with a leading space, 'NA' for none)."""
    # First pass: collect gene coords/strand/chrom from the template.
    gene = {}
    with open(template_gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'gene':
                gid = attr(f[8], 'ID')
                gene[gid] = (int(f[3].strip()), int(f[4].strip()), f[6], f[0])

    n_utr = n_na = 0
    with open(template_gff) as fh, open(out_gff, 'w') as fo:
        for line in fh:
            if line.startswith('#') or not line.strip():
                fo.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                fo.write(line)
                continue
            aid = attr(f[8], 'ID') or ''
            kind = None
            if aid.endswith('_five_prime_UTR'):
                kind = 'five'; g = aid[:-len('_five_prime_UTR')]
            elif aid.endswith('_three_prime_UTR'):
                kind = 'three'; g = aid[:-len('_three_prime_UTR')]
            if kind is None:
                fo.write(line)               # gene / CDS / intron: verbatim
                continue
            # recompute this UTR's coordinates from the table
            s, e, strand, chrom = gene[g]
            five_len, three_len = utr.get(g, (0, 0))
            span = utr_span(strand, s, e, five_len, three_len, chrlen.get(chrom, 10**12))
            v = span[kind]
            if v is None:
                f[3], f[4] = 'NA', 'NA'
                n_na += 1
            else:
                f[3], f[4] = ' %d' % v[0], ' %d' % v[1]
                n_utr += 1
            fo.write('\t'.join(f) + '\n')
    return n_utr, n_na


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--template-gff', required=True,
                    help='Curated gene-model backbone GFF3 with UTR slots '
                         '(Alice_Utrs_Final_forfasta_no_chr.gff3).')
    ap.add_argument('--utr-table', required=True,
                    help='UTR-length table (TSV/CSV: Gene, 5p_len, 3p_len).')
    ap.add_argument('--genome', required=True,
                    help='Genomic FASTA (SGD .fsa with ref| headers, or pre-renamed).')
    ap.add_argument('--out-prefix', required=True,
                    help='Output prefix, e.g. build/final.')
    ap.add_argument('--gffread', default='gffread', help='gffread executable.')
    args = ap.parse_args()

    out_dir = os.path.dirname(args.out_prefix) or '.'
    os.makedirs(out_dir, exist_ok=True)

    genome_chr = args.out_prefix + '_genome_chr_id.fa'
    gff_out    = args.out_prefix + '_UTRs.gff3'
    fa_out     = args.out_prefix + '_w.fa'

    print('[1/3] renaming genome chromosomes -> %s' % genome_chr, file=sys.stderr)
    nseq = rename_genome(args.genome, genome_chr)
    chrlen = load_chrlen(genome_chr)
    print('      %d sequences; lengths: %s' % (nseq, chrlen), file=sys.stderr)

    print('[2/3] overlaying UTRs from %s -> %s' % (args.utr_table, gff_out), file=sys.stderr)
    utr = load_utr_table(args.utr_table)
    n_utr, n_na = build_gff(args.template_gff, utr, chrlen, gff_out)
    print('      table genes: %d; UTR features written: %d; NA (skipped): %d'
          % (len(utr), n_utr, n_na), file=sys.stderr)

    print('[3/3] gffread -M -A -w -> %s' % fa_out, file=sys.stderr)
    cmd = [args.gffread, gff_out, '-g', genome_chr, '-w', fa_out, '-M', '-A']
    print('      ' + ' '.join(cmd), file=sys.stderr)
    r = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    sys.stderr.write(r.stderr)
    if r.returncode != 0:
        sys.exit('gffread failed (%d)' % r.returncode)
    print('Done: %s' % fa_out, file=sys.stderr)


if __name__ == '__main__':
    main()
