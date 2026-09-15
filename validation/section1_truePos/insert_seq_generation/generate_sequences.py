#!/usr/bin/env python3
"""
generate_insertion_sequences.py

Generates random, non-repetitive insertion sequences for VISOR HACk BED entries.
For each required insertion size, produces a sequence with balanced GC content
and verifies absence of simple repeats that would challenge de Bruijn assembly.

Short sequences (< 500 bp) use pure rejection sampling.
Long sequences (>= 500 bp) use constrained Markov chain generation, which
enforces GC balance and dinucleotide diversity by construction before applying
homopolymer and k-mer uniqueness checks as a final pass. This avoids the
geometric inefficiency of rejection sampling at large lengths.

Sequence generation is fully seeded for reproducibility.

Output:
    - insertion_sequences.fasta  : all generated sequences, one per entry
    - insertion_sequences.tsv    : lookup table (locus_id, size_bp, sequence, ...)
    - insertion_sequences.log    : per-sequence QC metrics and SHA256 checksums

Usage:
    python generate_insertion_sequences.py \\
        --sizes 50 200 500 1000 2000 5000 \\
        --replicates 5 \\
        --seed 42 \\
        --outdir insertion_sequences
"""

import argparse
import hashlib
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
import random


# ── QC constants ──────────────────────────────────────────────────────────────

GC_TARGET           = 0.50
GC_TOLERANCE        = 0.05       # accept 45–55% GC
MAX_HOMOPOLYMER     = 5          # reject if any homopolymer run >= this length
MAX_DINUC_FRACTION  = 0.40       # reject if any dinucleotide >= this fraction of all dinucs
KMER_K              = 15
MIN_KMER_UNIQUE     = 0.85       # fraction of k-mers that must be unique
MAX_ATTEMPTS        = 10_000     # per sequence (rejection sampling path)
LONG_SEQ_THRESHOLD  = 500        # bp; at or above this, use Markov generation


# ── Data class ────────────────────────────────────────────────────────────────

@dataclass
class InsertionSequence:
    locus_id: str
    size_bp: int
    sequence: str
    gc_fraction: float
    max_homopolymer_run: int
    max_dinuc_fraction: float
    kmer_unique_fraction: float
    sha256: str


# ── QC metrics ────────────────────────────────────────────────────────────────

def gc_fraction(seq: str) -> float:
    return (seq.count('G') + seq.count('C')) / len(seq)


def max_homopolymer(seq: str) -> int:
    return max((len(m.group()) for m in re.finditer(r'(.)\1+', seq)), default=1)


def max_dinucleotide_fraction(seq: str) -> float:
    if len(seq) < 4:
        return 0.0
    dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
    counts: dict[str, int] = {}
    for d in dinucs:
        counts[d] = counts.get(d, 0) + 1
    return max(counts.values()) / len(dinucs)


def kmer_unique_fraction(seq: str, k: int) -> float:
    if len(seq) < k:
        return 1.0
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    return len(set(kmers)) / len(kmers)


def qc_metrics(seq: str) -> dict:
    return {
        'gc':                  gc_fraction(seq),
        'max_homopolymer':     max_homopolymer(seq),
        'max_dinuc_fraction':  max_dinucleotide_fraction(seq),
        'kmer_unique_fraction': kmer_unique_fraction(seq, KMER_K),
    }


def passes_qc(metrics: dict) -> bool:
    return (
        abs(metrics['gc'] - GC_TARGET) <= GC_TOLERANCE
        and metrics['max_homopolymer'] < MAX_HOMOPOLYMER
        and metrics['max_dinuc_fraction'] < MAX_DINUC_FRACTION
        and metrics['kmer_unique_fraction'] >= MIN_KMER_UNIQUE
    )


# ── Sequence generation ───────────────────────────────────────────────────────

BASES = ['A', 'T', 'G', 'C']

def _rejection_sample(rng: random.Random, size: int) -> str:
    """Pure rejection sampling. Reliable for short sequences."""
    weights = [0.25, 0.25, 0.25, 0.25]
    return ''.join(rng.choices(BASES, weights=weights, k=size))


def _markov_generate(rng: random.Random, size: int) -> str:
    """
    Constrained Markov chain sequence generation.

    Builds the sequence one base at a time. At each step, the transition
    weights are adjusted to:
      1. Steer toward the GC target based on current GC deficit/surplus.
      2. Suppress the last-seen base to limit homopolymer runs.
      3. Suppress the last-seen dinucleotide to limit dinucleotide dominance.

    These soft constraints make the QC filters trivially satisfiable at
    any length without rejection sampling.
    """
    seq: list[str] = []
    gc_count = 0

    for i in range(size):
        remaining = size - i
        current_gc = gc_count / i if i > 0 else GC_TARGET
        gc_deficit = GC_TARGET - current_gc

        # Base weights: start uniform, adjust for GC balance
        w = {'A': 1.0, 'T': 1.0, 'G': 1.0, 'C': 1.0}

        # Steer GC toward target — strength scales with how far off we are
        gc_pull = min(3.0, abs(gc_deficit) * 10)
        if gc_deficit > 0:
            w['G'] += gc_pull
            w['C'] += gc_pull
        elif gc_deficit < 0:
            w['A'] += gc_pull
            w['T'] += gc_pull

        # Suppress current base to limit homopolymer runs
        if seq:
            last = seq[-1]
            run = 1
            j = len(seq) - 2
            while j >= 0 and seq[j] == last:
                run += 1
                j -= 1
            # Exponentially suppress as run length grows
            w[last] *= max(0.01, 1.0 / (2 ** (run - 1)))

        # Suppress last dinucleotide to limit dinucleotide dominance
        if len(seq) >= 1:
            prev = seq[-1]
            for b in BASES:
                dinuc = prev + b
                # Count occurrences of this dinucleotide so far
                count = sum(1 for k in range(len(seq) - 1) if seq[k] == prev and seq[k+1] == b)
                total_dinucs = max(1, len(seq) - 1)
                frac = count / total_dinucs
                if frac > 0.25:
                    w[b] *= max(0.05, 1.0 - (frac - 0.25) * 4)

        weights = [w[b] for b in BASES]
        chosen = rng.choices(BASES, weights=weights, k=1)[0]
        seq.append(chosen)
        if chosen in ('G', 'C'):
            gc_count += 1

    return ''.join(seq)


def make_insertion_sequence(rng: random.Random, size: int, locus_id: str,
                             log: logging.Logger) -> InsertionSequence:
    use_markov = size >= LONG_SEQ_THRESHOLD

    if use_markov:
        # Markov generation: try up to MAX_ATTEMPTS; in practice almost always
        # accepts on the first attempt since GC and dinucleotide constraints
        # are enforced by construction. Homopolymer and k-mer checks remain.
        for attempt in range(1, MAX_ATTEMPTS + 1):
            seq = _markov_generate(rng, size)
            m = qc_metrics(seq)
            if passes_qc(m):
                log.debug(f"{locus_id}: Markov accepted on attempt {attempt} — "
                          f"GC={m['gc']:.3f} homopoly={m['max_homopolymer']} "
                          f"dinuc={m['max_dinuc_fraction']:.3f} "
                          f"kmer_uniq={m['kmer_unique_fraction']:.3f}")
                break
        else:
            raise RuntimeError(
                f"Markov generation failed for {locus_id} (size={size}) "
                f"after {MAX_ATTEMPTS} attempts. This should not happen — "
                "check QC thresholds."
            )
    else:
        for attempt in range(1, MAX_ATTEMPTS + 1):
            seq = _rejection_sample(rng, size)
            m = qc_metrics(seq)
            if passes_qc(m):
                log.debug(f"{locus_id}: rejection sampling accepted on attempt {attempt} — "
                          f"GC={m['gc']:.3f} homopoly={m['max_homopolymer']} "
                          f"dinuc={m['max_dinuc_fraction']:.3f} "
                          f"kmer_uniq={m['kmer_unique_fraction']:.3f}")
                break
        else:
            raise RuntimeError(
                f"Rejection sampling failed for {locus_id} (size={size}) "
                f"after {MAX_ATTEMPTS} attempts. Try relaxing QC thresholds."
            )

    sha = hashlib.sha256(seq.encode()).hexdigest()
    return InsertionSequence(
        locus_id=locus_id,
        size_bp=size,
        sequence=seq,
        gc_fraction=m['gc'],
        max_homopolymer_run=m['max_homopolymer'],
        max_dinuc_fraction=m['max_dinuc_fraction'],
        kmer_unique_fraction=m['kmer_unique_fraction'],
        sha256=sha,
    )


# ── Output writers ────────────────────────────────────────────────────────────

def write_fasta(seqs: list[InsertionSequence], path: Path) -> None:
    with open(path, 'w') as fh:
        for s in seqs:
            fh.write(f">{s.locus_id} size={s.size_bp} gc={s.gc_fraction:.3f} "
                     f"sha256={s.sha256}\n")
            for i in range(0, len(s.sequence), 80):
                fh.write(s.sequence[i:i+80] + '\n')


def write_tsv(seqs: list[InsertionSequence], path: Path) -> None:
    header = '\t'.join(['locus_id', 'size_bp', 'gc_fraction', 'max_homopolymer_run',
                        'max_dinuc_fraction', 'kmer_unique_fraction', 'sha256', 'sequence'])
    with open(path, 'w') as fh:
        fh.write(header + '\n')
        for s in seqs:
            fh.write('\t'.join([
                s.locus_id, str(s.size_bp), f"{s.gc_fraction:.4f}",
                str(s.max_homopolymer_run), f"{s.max_dinuc_fraction:.4f}",
                f"{s.kmer_unique_fraction:.4f}", s.sha256, s.sequence,
            ]) + '\n')


def write_log_summary(seqs: list[InsertionSequence], path: Path, seed: int,
                      sizes: list[int], replicates: int) -> None:
    with open(path, 'w') as fh:
        fh.write("# generate_insertion_sequences.py — QC summary\n")
        fh.write(f"# seed={seed}  sizes={sizes}  replicates={replicates}\n")
        fh.write(f"# total_sequences={len(seqs)}\n")
        fh.write(f"# GC_TARGET={GC_TARGET}  GC_TOLERANCE=±{GC_TOLERANCE}\n")
        fh.write(f"# MAX_HOMOPOLYMER={MAX_HOMOPOLYMER}  "
                 f"MAX_DINUC_FRACTION={MAX_DINUC_FRACTION}\n")
        fh.write(f"# KMER_K={KMER_K}  MIN_KMER_UNIQUE={MIN_KMER_UNIQUE}\n")
        fh.write(f"# LONG_SEQ_THRESHOLD={LONG_SEQ_THRESHOLD} bp "
                 f"(Markov generation at or above this size)\n\n")
        fh.write('\t'.join(['locus_id', 'size_bp', 'gc', 'max_homopoly',
                            'max_dinuc', 'kmer_uniq_frac', 'sha256']) + '\n')
        for s in seqs:
            fh.write('\t'.join([
                s.locus_id, str(s.size_bp), f"{s.gc_fraction:.4f}",
                str(s.max_homopolymer_run), f"{s.max_dinuc_fraction:.4f}",
                f"{s.kmer_unique_fraction:.4f}", s.sha256,
            ]) + '\n')


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sizes', nargs='+', type=int,
        default=[50, 200, 500, 1000, 2000, 5000],
        help='Insertion sizes in bp (default: 50 200 500 1000 2000 5000)')
    p.add_argument('--replicates', type=int, default=5,
        help='Number of independent sequences per size (default: 5)')
    p.add_argument('--seed', type=int, default=42,
        help='Random seed for reproducibility (default: 42)')
    p.add_argument('--outdir', type=str, default='insertion_sequences',
        help='Output directory (default: insertion_sequences)')
    p.add_argument('--verbose', action='store_true',
        help='Log per-attempt debug messages')
    return p.parse_args()


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger(__name__)

    log.info(f"Seed: {args.seed}")
    log.info(f"Sizes: {args.sizes}")
    log.info(f"Replicates per size: {args.replicates}")
    log.info(f"Output directory: {outdir}")
    log.info(f"Generation method: rejection sampling below {LONG_SEQ_THRESHOLD} bp, "
             f"Markov chain at or above")

    rng = random.Random(args.seed)
    seqs: list[InsertionSequence] = []

    for size in sorted(args.sizes):
        for rep in range(1, args.replicates + 1):
            locus_id = f"ins_{size}bp_rep{rep:02d}"
            log.info(f"Generating {locus_id} ...")
            s = make_insertion_sequence(rng, size, locus_id, log)
            seqs.append(s)

    write_fasta(seqs, outdir / 'insertion_sequences.fasta')
    write_tsv(seqs, outdir / 'insertion_sequences.tsv')
    write_log_summary(seqs, outdir / 'insertion_sequences.log',
                      args.seed, args.sizes, args.replicates)

    log.info(f"Done. {len(seqs)} sequences written to {outdir}/")
    log.info("Files: insertion_sequences.fasta, insertion_sequences.tsv, "
             "insertion_sequences.log")


if __name__ == '__main__':
    main()
