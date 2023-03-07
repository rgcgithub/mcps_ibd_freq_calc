import argparse
import logging
import numpy as np
import os
import pysam
import sys
import time

from ibd_freq_estimator.ibd_file import IBDFile, InvalidIBDFileError
from ibd_freq_estimator.ancestry_file import AncestryFileRefAlt
from ibd_freq_estimator.fam_file import FamFile

import ibd_freq_estimator.logging_config
logger = logging.getLogger(__name__)


def estimate_frequencies(
    vcf_file,
    ibd_file,
    out_file,
    chrom,
    ancestry_file=None,
    fam_file=None
):
    """Estimate relatedness corrected allele frequencies using IBD segments.

    Parameters
    ----------
        vcf_file : pysam.VariantFile
            pysam file handle for indexed VCF/BCF file
        ibd_file : IBDFile
            file handle to iterate through IBD segments
        out_file : file handle
            file handle to store output results
        chrom : str
            chromosome corresponding to the current VCF file
        ancestry_file : AncestryFileRefAlt
            file handle encoding ancestry estimates using the REF/ALT fields in
            the VCF
        fam_file : FamFile
            file handle that encodes sex information, used for chromosome X
    """
    contigs = list(vcf_file.header.contigs)
    if not chrom in contigs and chrom in ["23", "X"]:
            chrom = "23" if chrom =="X" else "X"

    if not chrom in contigs:
        logger.error(f"chrom {chrom} is not a contig in the VCF")
        sys.exit(1)

    if len(contigs) > 1:
        logger.warning("found multiple contigs in VCF header (expected single chromosome)")
        contigs = [c for c in contigs if c == chrom]
    write_header(out_file, contigs[0], ancestry_file)

    if ancestry_file is not None:
        sample_intersect = intersect_samples(vcf_file, ancestry_file)
    else:
        sample_intersect = None

    variants_processed = 0
    for ibd_graph in ibd_file:
        start, end = ibd_graph.get_span()

        # all connected components greater than size 1
        conn_comp = ibd_graph.get_connected_components()

        if np.isfinite(end):
            region = f"{chrom}:{start}-{end}"
        # this can happen if the IBD file is empty
        elif np.isinf(start) and np.isinf(end):
            region = f"{chrom}:{1}-{1000000000}"
            start = 0
        else:
            region = f"{chrom}:{start}-{1000000000}"

        try:
            iterator = vcf_file.fetch(region=region)
        except IndexError:
            return

        for i,rec in enumerate(iterator):
            # skip indels where the start or end position falls
            # outside the requested region
            if rec.pos < start or rec.pos > end: continue

            if ancestry_file is not None:
                variant_ancestry = ancestry_file.fetch(rec.chrom, rec.pos)

                if variant_ancestry is None:
                    logger.warning(f"missing ancestry estimation at {rec.chrom}:{rec.pos}")

            else:
                variant_ancestry = None

            # (alt) allele count, allele number
            try:
                AC, AN, raw_AC, raw_AN, pop_AC_AN, raw_pop_AC_AN = compute_AC_AN(
                    conn_comp,
                    ibd_graph,
                    rec,
                    ibd_file.hap1_symbol,
                    ibd_file.hap2_symbol,
                    variant_ancestry,
                    sample_intersect,
                    fam_file
                )
            except Exception as e:
                logger.error(f"error at vcf pos {rec.pos}")
                raise e

            write_entry(out_file, rec, AC, AN, pop_AC_AN, raw_AC, raw_AN, raw_pop_AC_AN)
            variants_processed += 1

            if variants_processed % 500 == 0:
                logger.info(f"processed {variants_processed} variants")


def compute_AC_AN(
    conn_comp,
    ibd_graph,
    vcf_rec,
    hap1_symbol,
    hap2_symbol,
    variant_ancestry=None,
    sample_intersect=None,
    fam_file=None
):
    """Compute IBD adjusted alternate allele count and allele number
    """
    AC = 0
    AN = 0
    raw_AC = 0 # uncorrected AC
    raw_AN = 0 # uncorrected AN

    # Uncorrected estimate of allele count to find singletons.
    # The calculation for ancestry specific frequencies among
    # singltons is special. In this case we take the average
    # among the two haplotypes
    singleton_pop_AC = None # valid only when raw_AC = 1

    if variant_ancestry is not None:
        npops = len(variant_ancestry.pops)
        pop_AC = np.zeros(npops)
        pop_AN = np.zeros(npops)
        raw_pop_AC = np.zeros(npops)
        raw_pop_AN = np.zeros(npops)
        use_ancestry = True
    else:
        npops = 0
        pop_AC = np.zeros(0)
        pop_AN = np.zeros(0)
        raw_pop_AC  = np.zeros(0)
        raw_pop_AN = np.zeros(0)
        use_ancestry = False

    # for chromosome X, assume males are pseudo-haploid (i.e. alleles are
    # encoded by 0 or 2). we will use the first haplotype of each
    if fam_file is not None:
        is_chrX = True
    else:
        is_chrX = False

    for comp in conn_comp:
        comp_ac = 0
        comp_an = 0

        # Population effective allele count for ref alleles (0) or alt alleles (1).
        # To compute population specific allele counts / numbers, we take the average
        # ancestry proportion within a connected component. If there is no error in the
        # IBD assignments then only one of these sums will be nonzero. If there is an
        # error and the alleles in the connected component are mixed, then we split
        # the component by haplotypes with the 0 allele and haplotypes with the 1 allele.
        sum_anc0 = np.zeros(npops)
        sum_anc1 = np.zeros(npops)
        n_anc0 = 0
        n_anc1 = 0

        for graph_sid in comp:
            sample_id, hap_id = ibd_graph.get_sample_hap(graph_sid)

            if is_chrX and fam_file.is_missing(sample_id):
                continue

            # skip haplotype 2 for males on chromosome X
            if is_chrX and fam_file.is_male(sample_id) and hap_id == hap2_symbol:
                continue

            if hap_id == hap1_symbol:
                hap_id = 0
            elif hap_id == hap2_symbol:
                hap_id = 1
            else:
                raise ValueError(f"invalid haplotype symbol {hap_id}")

            allele = vcf_rec.samples[sample_id]["GT"][hap_id]
            if allele == 1:
                comp_ac += 1
                comp_an += 1
                raw_AC += 1
                raw_AN += 1
            elif allele == 0:
                comp_an += 1
                raw_AN += 1
            else:
                raise ValueError(
                    f"invalid allele {allele} in VCF: alleles must be 0 or 1"
                )

            sample_has_ancestry = use_ancestry and sample_id in sample_intersect
            if sample_has_ancestry and allele == 1:
                sum_anc1 += variant_ancestry.samples[sample_id][hap_id]
                n_anc1 += 1
                raw_pop_AC += variant_ancestry.samples[sample_id][hap_id]
                raw_pop_AN += variant_ancestry.samples[sample_id][hap_id]
            elif sample_has_ancestry and allele == 0:
                sum_anc0 += variant_ancestry.samples[sample_id][hap_id]
                n_anc0 += 1
                raw_pop_AN += variant_ancestry.samples[sample_id][hap_id]

            if sample_has_ancestry and raw_AC == 1 and allele == 1:
                singleton_pop_AC = np.mean(variant_ancestry.samples[sample_id], axis=0)

        # We need to decide what happens if the genotypes among
        # a connected component are not consistent. If we believe
        # the genotyping over the IBD calling, then we should
        # split the connected components in two by reference and
        # alt alleles. In this case, the component contributes
        # 2 to the allele number instead of 1.
        if comp_ac > 0 and comp_ac == comp_an:
            AC += 1
            AN += 1
        elif comp_ac > 0 and comp_ac != comp_an:
            AC += 1
            AN += 2
        else:
            AN += 1

        if n_anc1 > 0:
            pop_AC += sum_anc1 / n_anc1
            pop_AN += sum_anc1 / n_anc1
        if n_anc0 > 0:
            pop_AN += sum_anc0 / n_anc0

    for sample_id in vcf_rec.samples:

        if is_chrX and fam_file.is_missing(sample_id):
            continue

        alleles = vcf_rec.samples[sample_id]["GT"]
        assert len(alleles) == 2, "samples must be phased and biallelic"

        # Returns -1 if the sample + haplotype not in a connected component.
        hap1_idx = ibd_graph.get_sample_hap_idx(sample_id, hap1_symbol)
        hap2_idx = ibd_graph.get_sample_hap_idx(sample_id, hap2_symbol)

        sample_has_ancestry = use_ancestry and sample_id in sample_intersect
        if hap1_idx == -1 and alleles[0] == 1:
            AC += 1
            AN += 1
            raw_AC += 1
            raw_AN += 1

            if sample_has_ancestry:
                anc = variant_ancestry.samples[sample_id][0]
                pop_AC += anc
                pop_AN += anc
                raw_pop_AC += anc
                raw_pop_AN += anc

            if sample_has_ancestry and raw_AC == 1:
                singleton_pop_AC = np.mean(variant_ancestry.samples[sample_id], axis=0)

        elif hap1_idx == -1 and alleles[0] == 0:
            AN += 1
            raw_AN += 1

            if sample_has_ancestry:
                anc = variant_ancestry.samples[sample_id][0]
                pop_AN += anc
                raw_pop_AN += anc

        # skip haplotype 2 for males on chromosome X
        if is_chrX and fam_file.is_male(sample_id):
            continue

        if hap2_idx == -1 and alleles[1] == 1:
            AC += 1
            AN += 1
            raw_AC += 1
            raw_AN += 1


            if sample_has_ancestry:
                anc = variant_ancestry.samples[sample_id][1]
                pop_AC += anc
                pop_AN += anc
                raw_pop_AC += anc
                raw_pop_AN += anc

            if sample_has_ancestry and raw_AC == 1:
                singleton_pop_AC = np.mean(variant_ancestry.samples[sample_id], axis=0)

        elif hap2_idx == -1 and alleles[1] == 0:
            AN += 1
            raw_AN += 1

            if sample_has_ancestry:
                anc = variant_ancestry.samples[sample_id][1]
                pop_AN += anc
                raw_pop_AN += anc

    if use_ancestry:
        pop_AC = singleton_pop_AC if raw_AC == 1 else pop_AC
        pop_AC_AN = {}
        pops = variant_ancestry.pops
        for pop, pop_AC, pop_AN in zip(pops, pop_AC, pop_AN):
            pop_AC_AN[pop] = [pop_AC, pop_AN]

        raw_pop_AC = singleton_pop_AC if raw_AC == 1 else raw_pop_AC
        raw_pop_AC_AN = {}
        pops = variant_ancestry.pops
        for pop, pop_AC, pop_AN in zip(pops, raw_pop_AC, raw_pop_AN):
            raw_pop_AC_AN[pop] = [pop_AC, pop_AN]
    else:
        pop_AC_AN = {}
        raw_pop_AC_AN = {}

    return AC, AN, raw_AC, raw_AN, pop_AC_AN, raw_pop_AC_AN


def intersect_samples(vcf_file, ancestry_file):
    vcf_samples = set(vcf_file.header.samples)
    ancestry_samples = set(ancestry_file.get_samples())
    intersection = vcf_samples.intersection(ancestry_samples)

    nvcf = len(vcf_samples)
    nancestry = len(ancestry_samples)
    nintersect = len(intersection)
    if nintersect == 0:
        logger.error("the ancestry file and vcf file have no samples in common")
        sys.exit(1)
    if nintersect != nvcf or \
        nintersect != nancestry:
        logger.warning(
            f"the vcf file ({nvcf} samples) and ancestry file " + \
            f"({nancestry} samples) only have {nintersect} samples in common"
        )
    return intersection


def write_header(out_file, contig, ancestry_file=None):
    header = \
f"""##fileformat=VCFv4.3
##INFO=<ID=AC,Number=1,Type=Integer,Description="IBD corrected alternate allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="IBD corrected allele number">
##INFO=<ID=AF,Number=1,Type=Float,Description="IBD corrected allele frequency">
##INFO=<ID=RAW_AC,Number=1,Type=Integer,Description="uncorrected alternate allele count">
##INFO=<ID=RAW_AN,Number=1,Type=Integer,Description="uncorrected allele number">
##INFO=<ID=RAW_AF,Number=1,Type=Float,Description="uncorrected allele frequency">
"""
    header += f"##contig=<ID={contig}>\n"
    out_file.write(header)

    if ancestry_file is not None:
        pops = ancestry_file.get_pops()
        for pop in pops:
            out_file.write(
                f'##INFO=<ID={pop}_AC,Number=1,Type=Float,Description="IBD corrected alternate allele count in population {pop}">\n'
            )
            out_file.write(
                f'##INFO=<ID={pop}_AN,Number=1,Type=Float,Description="IBD corrected allele number in population {pop}">\n'
            )
            out_file.write(
                f'##INFO=<ID={pop}_AF,Number=1,Type=Float,Description="IBD corrected allele frequency in population {pop}">\n'
            )
            out_file.write(
                f'##INFO=<ID=RAW_{pop}_AC,Number=1,Type=Float,Description="uncorrected alternate allele count in population {pop}">\n'
            )
            out_file.write(
                f'##INFO=<ID=RAW_{pop}_AN,Number=1,Type=Float,Description="uncorrected allele number in population {pop}">\n'
            )
            out_file.write(
                f'##INFO=<ID=RAW_{pop}_AF,Number=1,Type=Float,Description="uncorrected allele frequency in population {pop}">\n'
            )

    cpra_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    out_file.write(cpra_line)


def write_entry(
    out_file,
    vcf_rec,
    AC,
    AN,
    pop_AC_AN,
    raw_AC,
    raw_AN,
    raw_pop_AC_AN
):
    info = f"AC={AC};AN={AN};AF={AC/AN:.6g}"
    info += f";RAW_AC={raw_AC};RAW_AN={raw_AN};RAW_AF={raw_AC/raw_AN:.6g}"
    for pop in pop_AC_AN:
        pop_ac = pop_AC_AN[pop][0]
        pop_an = pop_AC_AN[pop][1]
        pop_freq = pop_ac / pop_an if pop_an > 0 else "."
        info += ";" + f"{pop}_AC={pop_ac:.4f};{pop}_AN={pop_an:.4f};{pop}_AF={pop_freq:.6g}"

        raw_pop_ac = raw_pop_AC_AN[pop][0]
        raw_pop_an = raw_pop_AC_AN[pop][1]
        raw_pop_freq = raw_pop_ac / raw_pop_an if raw_pop_an > 0 else "."
        info += ";" + f"RAW_{pop}_AC={raw_pop_ac:.4f};RAW_{pop}_AN={raw_pop_an:.4f};RAW_{pop}_AF={raw_pop_freq:.6g}"

    entry = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        vcf_rec.chrom,
        vcf_rec.pos,
        vcf_rec.id,
        vcf_rec.ref,
        ",".join(vcf_rec.alts),
        ".",
        ".",
        info
    )
    out_file.write(entry)


def main():
    parser = argparse.ArgumentParser(
        description="compute IBD aware allele frequencies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "vcf_file",
        help="path to VCF/BCF file for a single chromosome"
    )
    parser.add_argument(
        "ibd_file",
        help="path to IBD segments for a single chromosome"
    )
    parser.add_argument(
        "chrom",
        help="chromosome string corresponding to VCF and IBD files"
    )
    parser.add_argument(
        "out_prefix",
        help="prefix of output file"
    )
    parser.add_argument(
        "--ancestry-file",
        help="path to VCF/BCF file with local ancestry estimates"
    )
    parser.add_argument(
        "--fam-file",
        help="path to fam file with sample sex (used for chromosome X)"
    )
    parser.add_argument(
        "--sid1-col",
        help="column index in IBD file for sample id 1",
        type=int,
        default=0
    )
    parser.add_argument(
        "--sid2-col",
        help="column index in IBD file for sample id 2",
        type=int,
        default=2
    )
    parser.add_argument(
        "--hid1-col",
        help="column index in IBD file for haplotype id 1",
        type=int,
        default=1
    )
    parser.add_argument(
        "--hid2-col",
        help="column index in IBD file for haplotype id 2",
        type=int,
        default=3
    )
    parser.add_argument(
        "--chrom-col",
        help="column index in IBD file for chromosome",
        type=int,
        default=4
    )
    parser.add_argument(
        "--start-pos-col",
        help="column index in IBD file for segment start position",
        type=int,
        default=5
    )
    parser.add_argument(
        "--end-pos-col",
        help="column index in IBD file for segmend end position",
        type=int,
        default=6
    )
    parser.add_argument(
        "--skip-ibd-header",
        help="if set assumes IBD file has a header line",
        action="store_true"
    )
    parser.add_argument(
        "--hap1-symbol",
        default=None,
        help="Symbol for haplotype 1 in IBD file (e.g. 0, 1)"
    )
    parser.add_argument(
        "--hap2-symbol",
        default=None,
        help="Symbol for haplotype 2 in in IBD file (e.g. 1, 2)"
    )
    args = parser.parse_args()

    for f in [args.vcf_file, args.ibd_file, args.ancestry_file, args.fam_file]:
        if f is not None and not os.path.exists(f):
            logger.error(f"input file {f} does not exist")
            sys.exit(1)


    logger.info("loading IBD file")
    try:
        ibd_file = IBDFile(
            args.ibd_file,
            args.sid1_col,
            args.sid2_col,
            args.hid1_col,
            args.hid2_col,
            args.chrom_col,
            args.start_pos_col,
            args.end_pos_col,
            args.skip_ibd_header,
            args.hap1_symbol,
            args.hap2_symbol
         )
    except InvalidIBDFileError as e:
        logger.error(e)
        sys.exit(1)

    logger.info("loading VCF file")
    vcf_file = pysam.VariantFile(args.vcf_file)

    if args.ancestry_file is not None:
        logger.info("loading ancestry file")
        ancestry_file = AncestryFileRefAlt(args.ancestry_file)
    else:
        ancestry_file = None


    if args.fam_file and args.chrom in ["X", "23"]:
        fam_file = FamFile(args.fam_file)
    elif args.fam_file:
        logger.error("fam file can only be used when chromosome is X or 23")
        sys.exit(1)
        fam_file = None
    else:
        fam_file = None

    logger.info("computing allele frequencies")
    with open(args.out_prefix + ".vcf", "w") as f:
        estimate_frequencies(
            vcf_file,
            ibd_file,
            f,
            args.chrom,
            ancestry_file,
            fam_file
        )


if __name__ == "__main__":
    main()
