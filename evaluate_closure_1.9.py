#!/usr/bin/env python3
# evaluate_closure_1.9.py
#Author: Chase McFarland

# Usage: python3 evaluate_closure_<version>.py \
#           --bam <contigs_vs_ref.bam> \
#           --vcf <variants.vcf> \
#			--eff-windows <effective_windows.tsv> \
#           --outdir <outdir> \
#           --report <output_closure_report.tsv>
#
# Optional parameters:
#   --threshold        <80.0>  Minimum fraction of breakpoint window covered by contig depth. (Float or Int)
#   --flank            <200>   A spanning contig must extend this far beyond each breakpoint end. (Int)
#                              Automatically relaxed at scaffold boundaries via effective_windows.tsv.
#   --max-gap          <50>    Maximum CIGAR deletion tolerated within the breakpoint window. (Int)
#   --window           <0>     Optional extension of depth evaluation region beyond SV boundaries. (Int)
#   --multimap-threshold <0.2> Discordant read fraction above which MM flag is set. (Float or Int)
#   -v / --verbose             Enable some extra diagnostic output about CIGAR deletion details.
#
# FLAGS column codes (comma-separated, empty string if none apply):
#   LT: assembly window truncated on the LEFT  by a scaffold boundary
#   RT: assembly window truncated on the RIGHT by a scaffold boundary
#   MM: discordant read fraction >= --multimap-threshold (see effective_windows.tsv)
#         MM also downgrades CLOSED to PARTIAL
#   SC: split contig: a contig spanning this locus has a supplementary alignment
#         on a different chromosome, indicating a possible inter-chromosomal
#         translocation junction

import argparse
import os
import subprocess
import sys
import pysam

# Complex SV types reported by CLOVE that are not amenable to gap closure.
# These involve breakpoints on different chromosomes or large
# intra-chromosomal events that can't be bridged by local de novo assembly.
# They are force-called OPEN and get the CX flag.
COMPLEX_SVTYPES = {'ITX', 'CTX', 'INV', 'CVRD', 'CIN'}

# Tandem duplication SV types are capped at PARTIAL regardless of spanning
# contig result. These get PARTIAL because SPAdes has a tendency to collapse
# repeat elements into a single copy, so definitive closure calls cannot
# be made and require user eval. IGV is a great tool for this.
TANDEM_CAP_SVTYPES = {'TAN', 'DUP'}



# ============ Helpers =============
 

# Easy runner wrapper for subprocess
def run(cmd, desc=""):
    print(f"  [RUN] {desc or cmd[:80]}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] Command failed:\n{result.stderr[-800:]}", file=sys.stderr)
        sys.exit(1)
    return result.stdout

# VCF parser, going with a list return here instead of dict for iteration purposes
def parse_vcf_loci(vcf):
    loci = []
    with open(vcf) as f:
        for line in f:
            if line.startswith('#'): #skip comment lines
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0]
            pos   = int(fields[1])
            sv_id = fields[2] if fields[2] != '.' else f"{chrom}_{pos}"
            end   = pos
            for item in fields[7].split(';'):
                if item.startswith('END='):
                    end = int(item.split('=')[1])
                elif item.startswith('SVLEN='):
                    svlen = abs(int(item.split('=')[1]))
                    if svlen > 0:
                        end = pos + svlen
            loci.append((chrom, pos, end, sv_id))
    return loci


#Eff. Windows parser
def load_effective_windows(tsv_path):

    eff = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return eff
    with open(tsv_path) as f:
        for line in f:
            if line.startswith('SV_ID'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            sv_id, chrom, pos, end, left_w, right_w, truncated, side = parts[:8]
            mean_depth      = float(parts[8])  if len(parts) >= 9  else None
            mean_mapq       = float(parts[9])  if len(parts) >= 10 else None
            discordant_frac = float(parts[10]) if len(parts) >= 11 else None
            svtype    = parts[11] if len(parts) >= 12 else 'UNKNOWN'
            clove_sup = parts[12] if len(parts) >= 13 else '.'
            chr2      = parts[13] if len(parts) >= 14 else '.'
            eff[sv_id] = {
                'left':            int(left_w),
                'right':           int(right_w),
                'truncated':       truncated.lower() == 'true',
                'side':            side,
                'mean_depth':      mean_depth,
                'mean_mapq':       mean_mapq,
                'discordant_frac': discordant_frac,
                'svtype':          svtype,
                'clove_sup':       clove_sup,
                'chr2':            chr2,
            }
    return eff

# Load the contig depth dictionary after it's written to file, this is efficient
# because most of the depth dictionary is 0 where no contig exists so memory is saved.
# Depth is calculated across the entire genome to get exact start/end coords of the
# contig itself to see if flank boundaries are satisfied rather than relying on
# start/end coordinates reported by the VCF
def load_depth(depth_file):
    depth = {}
    with open(depth_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            c, p, d = parts[0], int(parts[1]), int(parts[2])
            depth[(c, p)] = d
    return depth

# Compute depth across contig using depth dictionary
def compute_depth(depth_map, chrom, start, end):
    span  = max(end - start + 1, 1)
    total, count = 0, 0
    for pos in range(start, end + 1):
        d = depth_map.get((chrom, pos), 0)
        if d > 0:
            count += 1
        total += d
    mean_depth  = total / span
    pct_covered = count / span * 100
    return mean_depth, pct_covered


# Check for CIGAR deletion greater than max-gap
def cigar_has_large_gap(cigar_tuples, ref_start, breakpoint_pos,
                        breakpoint_end, max_gap=50):
    ref_pos = ref_start
    for op, length in cigar_tuples:
        if op in (0, 7, 8):
            ref_pos += length
        elif op in (2, 3):
            gap_start = ref_pos
            gap_end   = ref_pos + length
            if (gap_start <= breakpoint_end and gap_end >= breakpoint_pos
                    and length > max_gap):
                return True
            ref_pos += length
        elif op == 1:
            pass
        elif op in (4, 5):
            pass
    return False

# The next few functions were written by Claude Opus 4.6, the complicated ones are commented by me,
# making sure it did what I told it to
def check_multimapped_reads(sv_id, eff_windows, multimap_threshold=0.2):
    """
    Flag loci where the discordant read fraction recorded in
    effective_windows.tsv exceeds the multimap threshold.

    The discordant read fraction is computed by local_assembly.py over the
    original WGS BAM, using three complementary signals: RNEXT chromosome
    mismatch, NH tag > 1, and MAPQ == 0. Sourcing from the WGS BAM is
    critical — the contig BAM will show high MAPQ even at multimapped loci
    because SPAdes assembles flanking unique reads alongside multimapped
    repeat reads into a confidently-mapping spanning contig.

    MM flag also downgrades CLOSED → PARTIAL in the status assignment.

    Returns
    -------
    flagged        : bool  — True if discordant_frac >= threshold
    multimap_frac  : float — discordant read fraction from effective_windows.tsv
    """
    ew            = eff_windows.get(sv_id, {})
    multimap_frac = ew.get('discordant_frac') or 0.0
    flagged       = multimap_frac >= multimap_threshold
    return flagged, multimap_frac


def check_split_contig(bam_path, chrom, sv_id):
    """
    Detect whether any contig assembled for this locus has a supplementary
    alignment on a different chromosome, indicating the contig spans a
    junction between two chromosomal loci — the hallmark of an
    inter-chromosomal translocation breakpoint.

    minimap2 with -ax asm5 writes chimeric contig alignments as a primary
    record at one location and a supplementary record (FLAG 2048) at the
    other, with the SA tag on the primary encoding the supplementary
    coordinates. This function checks both the supplementary FLAG and the
    SA tag, so it is robust to BAMs produced with or without --secondary=no.

    Contig read names are expected to begin with the SV ID prefix as written
    by local_assembly.py (format: <sv_id>__<spades_contig_name>). Only
    contigs tagged with this locus's SV ID are examined, preventing
    cross-locus contamination in the contig BAM.

    Parameters
    ----------
    bam_path : path to contigs_vs_ref.bam (must be indexed)
    chrom    : chromosome of the variant locus
    sv_id    : variant ID used to identify relevant contigs by name prefix

    Returns
    -------
    split : bool — True if any contig for this locus has a supplementary
                   alignment on a different chromosome
    """
    bam   = pysam.AlignmentFile(bam_path, "rb")
    split = False

    # Collect all contig names belonging to this variant from the alignment BAM.
    # Retrieve using the chrom ID. Contig names are filtered by sv_id prefix below
    try:
        for read in bam.fetch(chrom):
            if not (read.query_name or '').startswith(sv_id + '__'):
                continue

            # Check 1 for whether this record has a supplementary alignment on a
            # different chromosome, implies the primary is elsewhere and this
            # end spans an inter-chromosomal breakend.
            if read.is_supplementary and read.reference_name != chrom:
                split = True
                break

            # Check 2 for SA tag indicating secondary alignments to the 
            # evaluated primary. Each entry is chr,pos,strand,CIGAR,mapQ,NM 
            # separated by commas. Multiple entries are semicolon delim.
            if not read.is_secondary and not read.is_supplementary:
                try:
                    sa = read.get_tag('SA')
                    for entry in sa.rstrip(';').split(';'):
                        sa_chrom = entry.split(',')[0]
                        if sa_chrom != chrom:
                            split = True
                            break
                except KeyError:
                    pass  # No SA tag = not a chimeric alignment

            if split:
                break

    except ValueError:
        pass  # chrom not in BAM index = no contigs mapped here

    bam.close()
    return split


def check_contig_deletion_at_breakpoint(bam_path, chrom, sv_pos, sv_end,
                                          left_flank=200, right_flank=200,
                                          max_gap=50, verbose=False):
    """
    Detect whether a spanning contig contains a deletion CIGAR operation
    overlapping the breakpoint interval [sv_pos, sv_end], and whether the
    true variant deletion length is within the --max-gap tolerance.

    The CIGAR string is used only to confirm that a deletion is present at
    the variant locus — not to measure deletion size. The true deletion
    length is derived from the VCF coordinates as (sv_end - sv_pos), which
    is the ground truth regardless of how minimap2 reports the CIGAR. This
    avoids the left-alignment offset problem: minimap2 shifts deletions as
    far left as possible, producing CIGAR lengths that are systematically
    longer than the true deletion size by a locus-specific offset equal to
    the number of bases by which the deletion was shifted. Using sv_end -
    sv_pos directly bypasses this artefact entirely.

    This check is applied after check_spanning() and check_split_contig().
    When a deletion CIGAR is found overlapping [sv_pos, sv_end] on a
    qualifying spanning alignment, the depth criterion is bypassed and
    closure is determined solely by whether (sv_end - sv_pos) <= max_gap.

    The check is scoped to deletion CIGAR operations that overlap
    [sv_pos, sv_end] only — incidental deletions in flanking sequence do
    not trigger the bypass.

    Parameters
    ----------
    bam_path               : path to contigs_vs_ref.bam
    chrom, sv_pos, sv_end  : variant locus; (sv_end - sv_pos) is the
                             authoritative deletion length for gap comparison
    left_flank, right_flank: flank requirements (same as check_spanning)
    max_gap                : maximum tolerated deletion size in bp

    Returns
    -------
    has_deletion  : bool  — True if any spanning contig has a deletion
                            CIGAR overlapping [sv_pos, sv_end]
    within_gap    : bool  — True if (sv_end - sv_pos) <= max_gap
                            (evaluated only when has_deletion is True)
    true_del_size : int   — sv_end - sv_pos (0 if has_deletion is False)
    cigar_del_size: int   — largest CIGAR deletion size found at the
                            breakpoint, reported for diagnostic purposes
    """
    bam            = pysam.AlignmentFile(bam_path, "rb")
    has_deletion   = False
    cigar_del_size = 0

    # Deletion length from VCF coordinates, takes precedence over CIGAR string
    true_del_size = sv_end - sv_pos
    within_gap    = true_del_size <= max_gap

    fetch_start = max(0, sv_pos - left_flank)
    fetch_end   = sv_end + right_flank

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start is None or read.reference_end is None:
            continue

        # Only examine reads that qualify as spanning alignments
        if not (read.reference_start <= sv_pos - left_flank and
                read.reference_end   >= sv_end  + right_flank):
            continue

        # Walk CIGAR string looking for any deletion overlapping start/end coords.
        # CIGAR length is recorded for diagnostics only 
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 7, 8):       # M / = / X
                ref_pos += length
            elif op in (2, 3):         # D / N
                gap_start = ref_pos
                gap_end   = ref_pos + length
                if gap_start <= sv_end and gap_end >= sv_pos:
                    has_deletion = True
                    if length > cigar_del_size:
                        cigar_del_size = length
                    if verbose:
                        print(f"      [CIGAR-DIAG] {chrom}:{sv_pos}-{sv_end} "
                              f"deletion op: {length}D at ref {gap_start}-{gap_end} "
                              f"overlaps breakpoint interval "
                              f"[true_del={true_del_size} bp, "
                              f"within_gap={within_gap}, max_gap={max_gap}, "
                              f"cigar_len={length} (diagnostic only)]")
                elif verbose:
                    print(f"      [CIGAR-DIAG] {chrom}:{sv_pos}-{sv_end} "
                          f"deletion op: {length}D at ref {gap_start}-{gap_end} "
                          f"does NOT overlap [{sv_pos},{sv_end}] — ignored")
                ref_pos += length
            elif op == 1:              # I
                pass
            elif op in (4, 5):         # S / H
                pass

        if has_deletion:
            break  # one qualifying spanning contig is sufficient

    bam.close()
    if verbose:
        print(f"      [CIGAR-DIAG] {chrom}:{sv_pos}-{sv_end} "
              f"summary: has_deletion={has_deletion} "
              f"true_del_size={true_del_size} within_gap={within_gap} "
              f"cigar_del_size={cigar_del_size} max_gap={max_gap}")
    return has_deletion, within_gap, true_del_size, cigar_del_size


def contig_extent(bam_path, chrom, sv_pos, sv_end, search_flank=2000):
    """
    Return (leftmost_start, rightmost_end) across all primary contig
    alignments overlapping the breakpoint region.
    """
    bam       = pysam.AlignmentFile(bam_path, "rb")
    leftmost  = None
    rightmost = None

    fetch_start = max(0, sv_pos - search_flank)
    fetch_end   = sv_end + search_flank

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start is None or read.reference_end is None:
            continue
        if leftmost  is None or read.reference_start < leftmost:
            leftmost  = read.reference_start
        if rightmost is None or read.reference_end   > rightmost:
            rightmost = read.reference_end

    bam.close()
    return leftmost, rightmost


def check_spanning(bam_path, chrom, sv_pos, sv_end,
                   left_flank=200, right_flank=200, max_gap=50):
    """
    Check whether any primary contig alignment spans the breakpoint with
    sufficient flank on each side and no large internal CIGAR gap.
    """
    bam      = pysam.AlignmentFile(bam_path, "rb")
    spanning = False

    fetch_start = max(0, sv_pos - left_flank)
    fetch_end   = sv_end + right_flank

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start is None or read.reference_end is None:
            continue

        if (read.reference_start <= sv_pos - left_flank and
                read.reference_end >= sv_end + right_flank):

            if not cigar_has_large_gap(
                    read.cigartuples, read.reference_start,
                    sv_pos, sv_end, max_gap):
                spanning = True
                break

    bam.close()
    return spanning


# Prep the depth dictionary for writing to file
def stage_depth(bam, outdir):
    depth_file = os.path.join(outdir, "contig_depth.txt")
    run(f"samtools depth -a {bam} > {depth_file}", "Computing contig depth")
    return depth_file



# =========== Main ===========


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate SV gap closure from contig BAM + original VCF."
    )
    parser.add_argument("--bam",         required=True,
                        help="contigs_vs_ref.bam (must be indexed)")
    parser.add_argument("--vcf",         required=True,
                        help="Original SV breakpoints VCF")
    parser.add_argument("--outdir",      default="closure_eval",
                        help="Working/output directory (default: closure_eval)")
    parser.add_argument("--eff-windows", required=True,
                        help="effective_windows.tsv from local_assembly.py")
    parser.add_argument("--window",      type=int, default=0,
                        help="Optional bp to extend depth evaluation region (default: 0)")
    parser.add_argument("--flank",       type=int, default=200,
                        help="Bp a spanning contig must extend beyond each "
                             "breakpoint end (default: 200)")
    parser.add_argument("--max-gap",     type=int, default=50,
                        help="Max tolerated CIGAR deletion within breakpoint "
                             "(default: 50)")
    parser.add_argument("--threshold",   type=float, default=80.0,
                        help="Min %% of breakpoint window covered by contigs "
                             "(default: 80)")
    parser.add_argument("--report",      default="closure_report.tsv",
                        help="Output TSV path (default: closure_report.tsv)")
    parser.add_argument("--multimap-threshold", type=float, default=0.2,
                        help="Discordant read fraction above which MM flag is "
                             "set and CLOSED is downgraded to PARTIAL (default: 0.2)")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Enable verbose diagnostic output. Activates per-locus "
                             "CIGAR deletion reporting to show "
                             "deletion sizes, breakpoint overlap, and gap threshold "
                             "comparisons. Useful for comparing different --max-gap values.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

	#check for bam index
    if not os.path.exists(args.bam + ".bai"):
        print(f"  [ERROR] BAM index not found: {args.bam}.bai", file=sys.stderr)
        sys.exit(1)

    print("\n=== Stage 1: Computing contig depth (for PCT_COVERED criterion) ===")
    depth_file = stage_depth(args.bam, args.outdir)

    print("\n=== Stage 2: Loading contig depth into memory ===")
    depth_map = load_depth(depth_file)
    print(f"  Loaded {len(depth_map):,} positions.")

    print("\n=== Stage 3: Loading effective window data ===")
    eff_windows = load_effective_windows(args.eff_windows)
    if eff_windows:
        n_trunc = sum(1 for v in eff_windows.values() if v['truncated'])
        print(f"  {len(eff_windows)} loci loaded; {n_trunc} truncated by scaffold boundary.")
    else:
        print("  No effective_windows.tsv provided, parameterized flank used for all loci.")

    print("\n=== Stage 4: Evaluating breakpoints ===")
    loci   = parse_vcf_loci(args.vcf)
    counts = {"CLOSED": 0, "PARTIAL": 0, "OPEN": 0}

    with open(args.report, 'w') as out:
        out.write(f"# evaluate_closure_1.9.py version 1.9\n")
        out.write(f"# invocation: {' '.join(sys.argv)}\n")
        out.write("SV_ID\tCHROM\tPOS\tEND\t"
                  "SVTYPE\tCLOVE_SUP\t"
                  "MEAN_DEPTH\tPCT_COVERED\tSPANNING_CONTIG\t"
                  "LEFT_FLANK_USED\tRIGHT_FLANK_USED\t"
                  "MULTIMAP_FRAC\tFLAGS\tSTATUS\n")

        for chrom, pos, end, sv_id in loci:
            ew        = eff_windows.get(sv_id, {})
            ew_depth  = ew.get('mean_depth')
            mean_d    = ew_depth if ew_depth is not None else 0.0
            svtype    = ew.get('svtype', 'UNKNOWN')
            clove_sup = ew.get('clove_sup', '.')

            w_start = max(0, pos - args.window)
            w_end   = end + args.window
            _, pct_cov = compute_depth(depth_map, chrom, w_start, w_end)
            depth_supported = pct_cov >= args.threshold

            #  Flank adjustment for scaffold-boundary loci
            trunc_side = ew.get('side') if ew.get('truncated') else None

            if ew.get('truncated'):
                c_left, c_right = contig_extent(
                    args.bam, chrom, pos, end, search_flank=args.flank + 2000
                )
                actual_left  = (pos   - c_left)  if c_left  is not None else 0
                actual_right = (c_right - end)    if c_right is not None else 0
                left_flank   = max(1, min(args.flank, actual_left))
                right_flank  = max(1, min(args.flank, actual_right))
                print(f"  Spanning check (boundary-adjusted): {sv_id} "
                      f"({chrom}:{pos}-{end}) "
                      f"contig reach L={actual_left}/R={actual_right} bp → "
                      f"flanks L={left_flank}/R={right_flank} bp")
            else:
                left_flank  = args.flank
                right_flank = args.flank
                print(f"  Spanning check: {sv_id} ({chrom}:{pos}-{end})")

            spanning = check_spanning(
                args.bam, chrom, pos, end,
                left_flank=left_flank, right_flank=right_flank,
                max_gap=args.max_gap
            )

            # Contig deletion at breakpoint check: 
            # Must run after SC check. If a spanning contig contains a
            # deletion CIGAR overlapping start/end coords, the deleted bases
            # will show zero depth, causing the depth criterion to fail
            # even with a valid closure result. In this case, bypass the depth
            # criterion and defer the call to the max-gap value.
            # Only applies to the breakpoint interval only, deletions
            # in flanking sequence don't trigger this method.
            has_bp_del, del_within_gap, true_del_size, cigar_del_size = \
                check_contig_deletion_at_breakpoint(
                    args.bam, chrom, pos, end,
                    left_flank=left_flank, right_flank=right_flank,
                    max_gap=args.max_gap,
                    verbose=args.verbose
                )
            if has_bp_del:
                print(f"  [DEL-CIGAR] {sv_id}: deletion at breakpoint confirmed "
                      f"(true size={true_del_size} bp from VCF, "
                      f"CIGAR reports {cigar_del_size} bp after left-alignment). "
                      f"Bypassing depth criterion, "
                      f"true size {'within' if del_within_gap else 'exceeds'} "
                      f"max-gap threshold ({args.max_gap} bp)")

            # Primary closure call:
            if has_bp_del:
                # When the depth check is bypassed, status is determined by max-gap check.
                # When del_within_gap is True, the CIGAR check has confirmed
                # a contig crosses the breakpoint with a deletion within max-gap,
                # which is sufficient for CLOSED regardless of whether the
                # spanning alignment check clears. The spanning check
                # requires the contig to extend the whole flank value beyond both
                # breakpoint ends, which might not make sense near scaffold boundaries
                # or when the contig ends just inside the flank, but
                # the contig already has breakend resolution.
                # When del_within_gap is False, the deletion exceeds max-gap
                # and the variant is not closed.
                if del_within_gap:
                    status = "CLOSED"
                else:
                    status = "OPEN"
            else:
                if depth_supported and spanning:
                    status = "CLOSED"
                elif depth_supported or spanning:
                    status = "PARTIAL"
                else:
                    status = "OPEN"

            # Init flag list
            flags = []

            # CX: complex SV type not amenable to gap closure and forced OPEN
            if svtype in COMPLEX_SVTYPES:
                flags.append('CX')
                status = "OPEN"
                print(f"  [CX] {sv_id}: SVTYPE={svtype}. Force-calling OPEN")

            if 'CX' not in flags:
                # TD: tandem duplication, cap status at PARTIAL. See above for rationale.
                if svtype in TANDEM_CAP_SVTYPES:
                    flags.append('TD')
                    if status == "CLOSED":
                        status = "PARTIAL"
                    print(f"  [TD] {sv_id}: SVTYPE={svtype}. Capping at PARTIAL")

                # LT/RT: scaffold boundary truncation on either side
                if trunc_side in ('left', 'both'):
                    flags.append('LT')
                if trunc_side in ('right', 'both'):
                    flags.append('RT')

                # MM: multimapped reads, downgrades CLOSED to PARTIAL for user review
                mm_flagged, mm_frac = check_multimapped_reads(
                    sv_id, eff_windows,
                    multimap_threshold=args.multimap_threshold
                )
                if mm_flagged:
                    flags.append('MM')
                    if status == "CLOSED":
                        status = "PARTIAL"
                    print(f"  [MM] {sv_id}: discordant frac {mm_frac:.1%} "
                          f">= threshold {args.multimap_threshold:.0%}")

                # SC: split contig across chromosomes, happens if the contig itself maps
				# to multiple loci due to being over a repeat region, opposite case of MM
				# where reads over repeat regions are used to generate the contig initially
                sc_flagged = check_split_contig(args.bam, chrom, sv_id)
                if sc_flagged:
                    flags.append('SC')
                    print(f"  [SC] {sv_id}: contig has supplementary alignment "
                          f"on a different chromosome. Possible translocation junction")
            else:
                mm_frac = 0.0  # no multimap check performed for CX variants

            flags_str = ','.join(flags) if flags else ''

            counts[status] += 1
            out.write(
                f"{sv_id}\t{chrom}\t{pos}\t{end}\t"
                f"{svtype}\t{clove_sup}\t"
                f"{mean_d:.2f}\t{pct_cov:.1f}\t{spanning}\t"
                f"{left_flank}\t{right_flank}\t"
                f"{mm_frac:.3f}\t{flags_str}\t{status}\n"
            )

    print(f"\nReport written to: {args.report}")
    print(f"  CLOSED:  {counts['CLOSED']}")
    print(f"  PARTIAL: {counts['PARTIAL']}")
    print(f"  OPEN:    {counts['OPEN']}")


if __name__ == "__main__":
    main()
