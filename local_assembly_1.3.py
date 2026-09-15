#!/usr/bin/env python3
# Author: Chase McFarland

# Version 1.3
# Usage: python3 local_assembly_<version>.py --vcf <clove_variants.vcf> --bam <alignments.bam> \
#                                   --outdir <outdir> --window <1500> --threads <threads>

# I had Claude Opus 4.6 write most of this script because it is basically a wrapper for SPAdes.
# All code was thoroughly evaluated for correctness. - Chase

import os
import sys
import subprocess
import argparse



# Complex SV types not amenable to gap closure
# Variants with these SVTYPE values involve breakpoints on different chromosomes
# or spanning large intra-chromosomal intervals that cannot be bridged by local
# de novo assembly. They are written to effective_windows.tsv for completeness
# but are skipped for read extraction and SPAdes. evaluate_closure.py force-calls
# them OPEN and sets the CX flag.
COMPLEX_SVTYPES = {'ITX', 'CTX', 'INV', 'CVRD', 'CIN'}



# ========= BAM utilities ==========

def get_chrom_lengths(bam_path):
    """
    Parse chromosome lengths from BAM header using samtools view -H.
    Returns dict {chrom: length}.
    """
    result = subprocess.run(
        f"samtools view -H {bam_path}",
        shell=True, capture_output=True, text=True, check=True
    )
    lengths = {}
    for line in result.stdout.splitlines():
        if line.startswith('@SQ'):
            parts = dict(field.split(':', 1) for field in line.split('\t')[1:])
            if 'SN' in parts and 'LN' in parts:
                lengths[parts['SN']] = int(parts['LN'])
    return lengths



# ========== VCF parser ==========

def parse_vcf(vcf_path, pass_only=True):
    """
    Parse VCF and return list of dicts with keys:
        chrom, pos, end, sv_id, svtype, chr2, sup

    END resolution priority:
      1. END= from INFO field (always used when present).
      2. For SVTYPE=DEL only: POS + |SVLEN| as fallback when END is absent.
         For all other types SVLEN is not used to derive END because it
         encodes inserted/inverted sequence length, not reference span.

    pass_only: if True, skip records whose FILTER field is not 'PASS'.
    """
    breakpoints = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            filt = fields[6] if len(fields) > 6 else '.'
            if pass_only and filt != 'PASS':
                continue

            chrom  = fields[0]
            pos    = int(fields[1])
            sv_id  = fields[2] if fields[2] != '.' else f"{chrom}_{pos}"
            info   = fields[7]

            # Parse INFO fields
            info_dict = {}
            for token in info.split(';'):
                if '=' in token:
                    k, v = token.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[token] = None

            svtype = info_dict.get('SVTYPE', 'UNKNOWN').strip()
            chr2   = info_dict.get('CHR2', chrom).strip()
            sup    = info_dict.get('SUP', '.').strip()

            # END resolution
            if 'END' in info_dict:
                end = int(info_dict['END'])
            elif svtype == 'DEL' and 'SVLEN' in info_dict:
                end = pos + abs(int(info_dict['SVLEN']))
            else:
                end = pos  # fallback — window will be centred on POS

            breakpoints.append({
                'chrom':  chrom,
                'pos':    pos,
                'end':    end,
                'sv_id':  sv_id,
                'svtype': svtype,
                'chr2':   chr2,
                'sup':    sup,
            })
    return breakpoints



# ========== Window boundary check ==========

def effective_window(chrom, sv_pos, sv_end, window, chrom_lengths):
    """
    Compute the actual assembly window after clamping to scaffold boundaries.

    Returns:
        region_start  : 0-based start passed to samtools (clamped to 0)
        region_end    : end passed to samtools (clamped to chrom length)
        left_window   : actual bp available to the left  of sv_pos
        right_window  : actual bp available to the right of sv_end
        truncated     : True if either side is shorter than --window
        truncated_side: 'left', 'right', 'both', or None
    """
    chrom_len     = chrom_lengths.get(chrom)
    desired_start = sv_pos - window
    desired_end   = sv_end + window

    region_start  = max(0, desired_start)
    region_end    = min(desired_end, chrom_len) if chrom_len else desired_end

    left_window   = sv_pos  - region_start
    right_window  = region_end - sv_end

    left_trunc  = desired_start < 0
    right_trunc = (chrom_len is not None) and (desired_end > chrom_len)
    truncated   = left_trunc or right_trunc

    if left_trunc and right_trunc:
        truncated_side = 'both'
    elif left_trunc:
        truncated_side = 'left'
    elif right_trunc:
        truncated_side = 'right'
    else:
        truncated_side = None

    return region_start, region_end, left_window, right_window, truncated, truncated_side



# ========= Read depth and mapping quality analysis ========

def mean_read_depth(bam, chrom, region_start, region_end):
    """
    Compute mean read depth over [region_start, region_end] from the
    original WGS BAM using samtools depth -a.
    Returns a float.
    """
    region = f"{chrom}:{region_start}-{region_end}"
    result = subprocess.run(
        f"samtools depth -a -r {region} {bam}",
        shell=True, capture_output=True, text=True, check=True
    )
    depths = [int(line.split('\t')[2])
              for line in result.stdout.splitlines() if line.strip()]
    if not depths:
        return 0.0
    return sum(depths) / len(depths)


def window_mapping_stats(bam_path, chrom, region_start, region_end):
    """
    Compute mean MAPQ and discordant read fraction over the assembly window
    from the original WGS BAM.

    A read is discordant if any of:
      1. Mate maps to a different chromosome (RNEXT != chrom)
      2. NH tag > 1 (aligner-reported multiple mapping positions)
      3. MAPQ == 0 (SAM convention for non-uniquely placed reads)

    Returns (mean_mapq, discordant_frac, n_reads).
    """
    import pysam
    bam        = pysam.AlignmentFile(bam_path, "rb")
    total      = 0
    discordant = 0
    mapq_sum   = 0

    for read in bam.fetch(chrom, max(0, region_start), region_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total    += 1
        mapq_sum += read.mapping_quality
        is_discordant = False

        if (read.is_paired and
                read.next_reference_id >= 0 and
                read.next_reference_name != chrom):
            is_discordant = True

        if not is_discordant:
            try:
                if read.get_tag('NH') > 1:
                    is_discordant = True
            except KeyError:
                pass

        if not is_discordant and read.mapping_quality == 0:
            is_discordant = True

        if is_discordant:
            discordant += 1

    bam.close()

    if total == 0:
        return 0.0, 0.0, 0

    return mapq_sum / total, discordant / total, total



# ======= Read extraction ========


def extract_reads(bam, chrom, region_start, region_end, out_prefix):
    """
    Extract reads and their mates from a pre-computed region, output as
    FASTQ pair. Also computes mean WGS read depth and mapping stats.

    Returns (fq1, fq2, fq_unp, mean_depth, mean_mapq, discordant_frac)
    or False if no reads are found.
    """
    region = f"{chrom}:{region_start}-{region_end}"

    readnames_file = f"{out_prefix}.readnames.txt"
    cmd_names = (
        f"samtools view {bam} {region} | cut -f1 | sort -u > {readnames_file}"
    )
    subprocess.run(cmd_names, shell=True, check=True)

    if os.path.getsize(readnames_file) == 0:
        print(f"  [WARN] No reads found in region {region}, skipping.")
        return False

    depth                         = mean_read_depth(bam, chrom, region_start, region_end)
    mean_mapq, discordant_frac, _ = window_mapping_stats(bam, chrom, region_start, region_end)

    bam_subset = f"{out_prefix}.subset.bam"
    cmd_extract = f"samtools view -b -N {readnames_file} {bam} > {bam_subset}"
    subprocess.run(cmd_extract, shell=True, check=True)

    fq1    = f"{out_prefix}_R1.fastq"
    fq2    = f"{out_prefix}_R2.fastq"
    fq_unp = f"{out_prefix}_unpaired.fastq"

    cmd_fq = (
        f"samtools sort -n {bam_subset} | "
        f"samtools fastq -1 {fq1} -2 {fq2} -0 {fq_unp} -s /dev/null -F 2048 -"
    )
    subprocess.run(cmd_fq, shell=True, check=True)

    return fq1, fq2, fq_unp, depth, mean_mapq, discordant_frac



# ========= SPAdes runner =========

def run_spades(fq1, fq2, fq_unp, out_dir, threads=4):
    """Run SPAdes in careful mode on extracted reads."""
    os.makedirs(out_dir, exist_ok=True)

    paired_ok = (os.path.exists(fq1) and os.path.getsize(fq1) > 0 and
                 os.path.exists(fq2) and os.path.getsize(fq2) > 0)
    unp_ok    = os.path.exists(fq_unp) and os.path.getsize(fq_unp) > 0

    if not paired_ok and not unp_ok:
        print(f"  [WARN] All FASTQ files empty after read extraction — "
              f"skipping SPAdes.")
        return None

    cmd = ["spades.py", "--careful", "-o", out_dir,
           "--threads", str(threads), "--memory", "8"]

    if paired_ok:
        cmd += ["-1", fq1, "-2", fq2]
    if unp_ok:
        cmd += ["-s", fq_unp]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  [ERROR] SPAdes failed:\n{result.stderr[-500:]}")
        return None

    contigs = os.path.join(out_dir, "contigs.fasta")
    if os.path.exists(contigs) and os.path.getsize(contigs) > 0:
        return contigs
    else:
        print(f"  [WARN] SPAdes produced no contigs in {out_dir}")
        return None



# ======= Contig collector =======

def collect_contigs(contig_file, sv_id, output_handle):
    """Append contigs to master output FASTA with SV ID in headers."""
    with open(contig_file) as f:
        for line in f:
            if line.startswith('>'):
                output_handle.write(f">{sv_id}__{line[1:]}")
            else:
                output_handle.write(line)



# ======= Contig alignment =======

def align_contigs(reference, master_fasta, outdir):
    """
    Align the master contig FASTA to the reference with minimap2, sort with
    samtools, and index the resulting BAM.

    Produces {outdir}/contigs_vs_ref.bam and {outdir}/contigs_vs_ref.bam.bai.
    """
    bam_out = os.path.join(outdir, "contigs_vs_ref.bam")
    cmd = (
        f"minimap2 -ax asm5 {reference} {master_fasta} | "
        f"samtools sort -o {bam_out} && "
        f"samtools index {bam_out}"
    )
    subprocess.run(cmd, shell=True, check=True)
    return bam_out



# ========= Main =========

def main():
    parser = argparse.ArgumentParser(
        description="Local assembly at SV breakpoints using samtools + SPAdes"
    )
    parser.add_argument("--vcf",       required=True, help="Input SV VCF file (CLOVE output)")
    parser.add_argument("--bam",       required=True, help="Indexed alignment BAM")
    parser.add_argument("--reference", required=True, help="Reference FASTA for contig alignment")
    parser.add_argument("--outdir",    required=True, help="Output directory")
    parser.add_argument("--window",    type=int, default=1500,
                        help="Bp window around breakpoints (default: 1500)")
    parser.add_argument("--threads",   type=int, default=4,
                        help="Threads per SPAdes run (default: 4)")
    parser.add_argument("--pass-only", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Only process PASS-filtered VCF records (default: True)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    tmp_dir    = os.path.join(args.outdir, "tmp")
    spades_dir = os.path.join(args.outdir, "spades_runs")
    os.makedirs(tmp_dir,    exist_ok=True)
    os.makedirs(spades_dir, exist_ok=True)

    master_fasta   = os.path.join(args.outdir, "all_assembled_contigs.fasta")
    eff_window_tsv = os.path.join(args.outdir, "effective_windows.tsv")

    print("Loading chromosome lengths from BAM header ...")
    chrom_lengths = get_chrom_lengths(args.bam)
    print(f"  {len(chrom_lengths)} chromosomes/scaffolds found.")

    pass_only = args.pass_only
    breakpoints = parse_vcf(args.vcf, pass_only=pass_only)
    print(f"Parsed {len(breakpoints)} breakpoints from VCF "
          f"({'PASS only' if pass_only else 'all records'}).")

    with open(master_fasta, 'w') as master_out, \
         open(eff_window_tsv, 'w') as ew_out:

        ew_out.write(f"# local_assembly.py version 1.2\n")
        ew_out.write(f"# invocation: {' '.join(sys.argv)}\n")
        ew_out.write(
            "SV_ID\tCHROM\tPOS\tEND\t"
            "LEFT_WINDOW\tRIGHT_WINDOW\tTRUNCATED\tTRUNCATED_SIDE\t"
            "MEAN_DEPTH\tMEAN_MAPQ\tDISCORDANT_FRAC\t"
            "SVTYPE\tCLOVE_SUP\tCHR2\n"
        )

        for i, bp in enumerate(breakpoints):
            chrom  = bp['chrom']
            pos    = bp['pos']
            end    = bp['end']
            sv_id  = bp['sv_id']
            svtype = bp['svtype']
            chr2   = bp['chr2']
            sup    = bp['sup']
            safe_id = sv_id.replace('/', '_').replace(':', '_')

            print(f"[{i+1}/{len(breakpoints)}] Processing {sv_id} "
                  f"({chrom}:{pos}-{end}) SVTYPE={svtype}")

            # Complex SV types — write sidecar placeholder and skip assembly
            if svtype in COMPLEX_SVTYPES:
                print(f"  [CX] {sv_id}: SVTYPE={svtype} is a complex SV type — "
                      f"skipping assembly. evaluate_closure.py will force-call OPEN.")
                ew_out.write(
                    f"{sv_id}\t{chrom}\t{pos}\t{end}\t"
                    f"0\t0\tFalse\tnone\t"
                    f"0.00\t0.00\t0.000\t"
                    f"{svtype}\t{sup}\t{chr2}\n"
                )
                continue

            region_start, region_end, left_win, right_win, truncated, trunc_side = \
                effective_window(chrom, pos, end, args.window, chrom_lengths)

            if truncated:
                print(f"  [NOTE] Assembly window truncated at scaffold boundary "
                      f"({trunc_side}): effective window "
                      f"left={left_win} bp, right={right_win} bp "
                      f"(requested {args.window} bp each side)")

            out_prefix = os.path.join(tmp_dir, safe_id)
            spades_out = os.path.join(spades_dir, safe_id)

            result = extract_reads(
                args.bam, chrom, region_start, region_end, out_prefix
            )
            if not result:
                ew_out.write(
                    f"{sv_id}\t{chrom}\t{pos}\t{end}\t"
                    f"{left_win}\t{right_win}\t{truncated}\t"
                    f"{trunc_side or 'none'}\t0.00\t0.00\t0.000\t"
                    f"{svtype}\t{sup}\t{chr2}\n"
                )
                continue

            fq1, fq2, fq_unp, mean_depth, mean_mapq, discordant_frac = result

            ew_out.write(
                f"{sv_id}\t{chrom}\t{pos}\t{end}\t"
                f"{left_win}\t{right_win}\t{truncated}\t"
                f"{trunc_side or 'none'}\t{mean_depth:.2f}\t"
                f"{mean_mapq:.2f}\t{discordant_frac:.3f}\t"
                f"{svtype}\t{sup}\t{chr2}\n"
            )

            contigs = run_spades(fq1, fq2, fq_unp, spades_out,
                                 threads=args.threads)

            if contigs:
                collect_contigs(contigs, sv_id, master_out)
                print(f"  [OK] Contigs written for {sv_id} "
                      f"(mean depth {mean_depth:.1f}x, "
                      f"mean MAPQ {mean_mapq:.1f}, "
                      f"discordant frac {discordant_frac:.3f})")

    print("Aligning contigs to reference ...")
    contigs_bam = align_contigs(args.reference, master_fasta, args.outdir)

    print(f"\nDone.")
    print(f"  Contigs          : {master_fasta}")
    print(f"  Contigs BAM      : {contigs_bam}")
    print(f"  Effective windows: {eff_window_tsv}")


if __name__ == "__main__":
    main()
