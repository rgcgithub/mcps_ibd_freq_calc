import numpy as np
import pysam
import subprocess as sub
import unittest
import xopen


from src.ibd_freq_estimator.estimate_frequencies import estimate_frequencies
from src.ibd_freq_estimator.ancestry_file import AncestryFileRefAlt
from src.ibd_freq_estimator.ibd_file import IBDFile


class TestEstimateFrequencies(unittest.TestCase):


    def test_vanilla_estimate_frequencies(self):
        vcf_file = pysam.VariantFile("./test/data/test.vcf.gz")
        ibd_file = IBDFile(
            "./test/data/ibd_segs.txt",
            sid1_col=0,
            sid2_col=2,
            hid1_col=1,
            hid2_col=3,
            chr_col=4,
            start_pos_col=5,
            end_pos_col=6,
            has_header=True
        )
        out_file = open("./test/data/test_estimates.vcf", "w")
        chrom = "1"
        estimate_frequencies(vcf_file, ibd_file, out_file, chrom)

        ibd_file.close()
        vcf_file.close()
        out_file.close()

        #true_ac_an = get_ac_an("./test/data/test.vcf.gz")
        #est_ac_an = get_ac_an("./test/data/test_estimates.vcf")
        true_ac_an = get_chrom_pos_info("./test/data/test.vcf.gz")
        est_ac_an = get_chrom_pos_info("./test/data/test_estimates.vcf")

        for chrom_pos in true_ac_an:
            self.assertTrue(
                chrom_pos in est_ac_an,
                msg=f"missing {chrom_pos} from estimate file"
            )

            for info_field in true_ac_an[chrom_pos]:
                self.assertTrue(
                    info_field in est_ac_an[chrom_pos],
                    msg=f"missing {info_field} from estimates"
                )
                t = true_ac_an[chrom_pos][info_field]
                e = est_ac_an[chrom_pos][info_field]
                self.assertTrue(
                    t == e,
                    msg=f"true {info_field} {t} != estimated {info_field} {e} at {chrom_pos}"
                )


    def test_ancestry_estimate_frequencies(self):
        vcf_file = pysam.VariantFile("./test/data/test_short.vcf.gz")
        ibd_file = IBDFile(
            "./test/data/ibd_segs.txt",
            sid1_col=0,
            sid2_col=2,
            hid1_col=1,
            hid2_col=3,
            chr_col=4,
            start_pos_col=5,
            end_pos_col=6,
            has_header=True
        )
        ancestry_file = AncestryFileRefAlt("./test/data/test_ancestry.vcf.gz")
        out_file = open("./test/data/test_ancestry_estimates.vcf", "w")
        chrom = "1"

        estimate_frequencies(vcf_file, ibd_file, out_file, chrom, ancestry_file)

        vcf_file.close()
        ibd_file.close()
        out_file.close()
        ancestry_file.close()

        true_chrom_pos_info = get_chrom_pos_info("./test/data/test_ancestry.vcf.gz")
        est_chrom_pos_info = get_chrom_pos_info("./test/data/test_ancestry_estimates.vcf")

        for chrom_pos in true_chrom_pos_info:
            self.assertTrue(
                chrom_pos in est_chrom_pos_info,
                msg="missing {chrom_pos} from ancestry estimates"
            )

            for field in true_chrom_pos_info[chrom_pos]:
                self.assertTrue(
                    true_chrom_pos_info[chrom_pos][field] == est_chrom_pos_info[chrom_pos][field],
                    msg="variant {chrom_pos}: true {field} {true_value} != est {field} {est_value}".format(
                        chrom_pos=chrom_pos,
                        field=field,
                        true_value=true_chrom_pos_info[chrom_pos][field],
                        est_value=est_chrom_pos_info[chrom_pos][field]
                    )
                )


    def test_relatedness_simulations(self):
        vcf_out = open("./test/data/test_rel.vcf", "w")
        ibd_out = open("./test/data/test_rel.ibd", "w")

        ACs, ANs = simulate_related_genotypes(vcf_out, ibd_out, N=100)

        vcf_out.close()
        ibd_out.close()

        sub.check_call(["bgzip", "-f", "./test/data/test_rel.vcf"])
        sub.check_call(["tabix", "./test/data/test_rel.vcf.gz"])

        vcf_file = pysam.VariantFile("./test/data/test_rel.vcf.gz")
        ibd_file = IBDFile(
            "./test/data/test_rel.ibd",
            sid1_col=0,
            sid2_col=2,
            hid1_col=1,
            hid2_col=3,
            chr_col=4,
            start_pos_col=5,
            end_pos_col=6,
            has_header=False
        )

        out_file = open("./test/data/test_rel_est", "w")
        chrom = "1"
        estimate_frequencies(vcf_file, ibd_file, out_file, chrom)
        out_file.close()

        est_AC_ANs = get_ac_an("./test/data/test_rel_est")

        est_ACs = []
        est_ANs = []
        for key in sorted(est_AC_ANs, key=lambda k: int(k.split(":")[1])):
            est_ACs.append(est_AC_ANs[key][0])
            est_ANs.append(est_AC_ANs[key][1])

        for i,(est_ac, ac) in enumerate(zip(est_ACs, ACs)):
            self.assertTrue(
                est_ac == ac,
                msg=f"AC mismatch for variant {i}"
            )

        for i,(est_an, an) in enumerate(zip(est_ANs, ANs)):
            self.assertTrue(
                est_an == an,
                msg=f"AN mismatch for variant {i}"
            )


def get_ac_an(filepath):
    # (chrom,pos) -> (AC, AN)
    ac_an = {}
    with xopen.xopen(filepath) as f:
        for line in f:
            if line[0] == "#": continue
            line = line.strip("\n").split()
            chrom_pos = line[0] + ":" + line[1]
            info = line[7].split(";")
            AC = info[0].split("=")[1]
            AN = info[1].split("=")[1]

            AC = int(AC)
            AN = int(AN)
            ac_an[chrom_pos] = [AC, AN]
    return ac_an


def get_chrom_pos_info(filepath):
    # (chrom_pos) -> {field: value}
    chrom_pos_info = {}
    with xopen.xopen(filepath) as f:
        for line in f:
            if line[0] == "#": continue
            line = line.strip("\n").split()
            chrom_pos = line[0] + ":" + line[1]

            info = line[7].split(";")
            chrom_pos_info[chrom_pos] = {}
            for field in info:
                field = field.split("=")
                chrom_pos_info[chrom_pos][field[0]] = float(field[1])
    return chrom_pos_info


def vcf_header(n_samples):
    header = \
"""##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AC,Number=1,Type=Integer,Description="Relatedness corrected alternate allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Relatedness corrected total number of alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
"""
    chrom_line = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    chrom_line += [str(i) for i in range(n_samples)]
    header += "\t".join(chrom_line) + "\n"
    return header


def simulate_ibd1_genotypes(freq):
    allele1 = np.random.binomial(1,freq)
    allele2 = np.random.binomial(1,freq)
    shared_allele = np.random.binomial(1,freq)

    ac = allele1 + allele2 + shared_allele
    an = 3

    hap_idx1 = np.random.randint(0,2)
    hap_idx2 = np.random.randint(0,2)

    gt1 = [0,0]
    gt2 = [0,0]

    gt1[hap_idx1] = shared_allele
    gt1[1 - hap_idx1] = allele1
    gt2[hap_idx2] = shared_allele
    gt2[1 - hap_idx2] = allele2
    return gt1, gt2, ac, an, hap_idx1, hap_idx2


def simulate_ibd0_genotypes(freq):
    gt1 = [
        np.random.binomial(1, freq),
        np.random.binomial(1, freq)
    ]
    gt2 = [
        np.random.binomial(1, freq),
        np.random.binomial(1, freq)
    ]
    ac = int(np.sum(gt1) + np.sum(gt2))
    an = 4
    return gt1, gt2, ac, an


def simulate_related_genotypes(out_vcf, out_ibd, N):
    """Simulate 2*N related samples.
    """
    np.random.seed(30047)
    out_vcf.write(vcf_header(2*N))

    chrom = 1
    ACs = []
    ANs = []
    positions = []
    genotypes1 = []
    genotypes2 = []
    ibd_segs = []

    freqs = np.linspace(0.1, 0.9, 50)
    for i,f in enumerate(freqs):
        pos = 100*(i+1)
        positions.append(pos)

        gts1 = []
        gts2 = []
        AC = 0
        AN = 0
        for n in range(N):
            start_seg_pos=np.random.randint(pos-50, pos)
            end_seg_pos=np.random.randint(pos, pos+50)

            idx1 = n
            idx2 = N + n

            rel_deg = np.random.randint(0, 4)

            ibd1 = False
            if rel_deg == 0:
                gt1, gt2, ac, an = simulate_ibd0_genotypes(f)
            elif rel_deg == 1:
                ibd1 = True
                gt1, gt2, ac, an, hap_idx1, hap_idx2 = simulate_ibd1_genotypes(f)
            elif rel_deg == 2:
                ibd1_prob = 0.25
                if np.random.uniform() < ibd1_prob:
                    ibd1 = True
                    gt1, gt2, ac, an, hap_idx1, hap_idx2 = simulate_ibd1_genotypes(f)
                else:
                    gt1, gt2, ac, an = simulate_ibd0_genotypes(f)
            else:
                ibd1_prob = 0.125
                if np.random.uniform() < ibd1_prob:
                    ibd1 = True
                    gt1, gt2, ac, an, hap_idx1, hap_idx2 = simulate_ibd1_genotypes(f)
                else:
                    gt1, gt2, ac, an = simulate_ibd0_genotypes(f)

            AC += ac
            AN += an

            gts1.append("{}|{}".format(gt1[0], gt1[1]))
            gts2.append("{}|{}".format(gt2[0], gt2[1]))

            if ibd1:
                ibd_segs.append(
                    [idx1, hap_idx1, idx2, hap_idx2, chrom, start_seg_pos, end_seg_pos]
                )

        genotypes1.append(gts1)
        genotypes2.append(gts2)
        ACs.append(AC)
        ANs.append(AN)

    ibd_segs = sorted(ibd_segs, key=lambda row: row[5])

    for ibd_seg in ibd_segs:
        ibd_seg = [str(s) for s in ibd_seg]
        out_ibd.write("\t".join(ibd_seg) + "\n")

    for i, (gts1, gts2) in enumerate(zip(genotypes1, genotypes2)):
        out_line = [
            "1",
            str(positions[i]),
            str(i),
            "A",
            "G",
            ".",
            ".",
            ".",
            "GT"
        ] + gts1 + gts2
        out_vcf.write("\t".join(out_line) + "\n")

    return ACs, ANs


if __name__ == "__main__":
    unittest.main()
