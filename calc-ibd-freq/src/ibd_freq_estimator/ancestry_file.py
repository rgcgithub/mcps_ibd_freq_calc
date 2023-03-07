import numpy as np
import pysam


class AncestryEstimates:
    """Data structure to store ancestry estimates.
    """

    def __init__(self, chrom, pos, pops):
        self.chrom = chrom
        self.pos = pos
        self.pops = pops
        # sample_id -> [np.array, np.array]
        # where each np.array gives the ancestry
        # estimates of each haplotype
        self.samples = {}


class MissingFormatFieldError(Exception):
    pass


class AncestryFile:
    """Abstract class parse ancestry files.
    """

    def fetch(self, chrom, pos):
        """Get ancestry estimates for haplotypes at chrom, pos

        Returns
        -------
            ancestry_estimates : AncestryEstimates
                Ancestry estimates for each haplotype at chrom-pos
        """
        raise NotImplementedError


    def get_pops(self):
        """Return the list of populations in the ancestry file.
        """
        raise NotImplementedError


    def close(self):
        """Close file handle
        """
        raise NotImplementedError


class AncestryFileRefAlt(AncestryFile):
    """Ancestry parser for estimates stored in a VCF file encoded as genotypes.
    VCF file is required to be indexed.

    Reference allele column enodes one population label. Alternate allele column
    encodes comma-separate list of remaining populations. Genotypes encode
    haplotype assignment probabilities. Haplotype probabilities are separated by
    the | character, lists of probabilities are comma-separated. See
    test/data/test_ancestry.vcf.gz for an example.

    Parameters
    ----------
        ancestry_file : str
            Path to ancestry file
    """

    def __init__(self, ancestry_file):
        self.filepath = ancestry_file
        self.anc_file = pysam.VariantFile(ancestry_file)
        self.pops = None
        self.samples = None


    def fetch(self, chrom, pos):
        # try alternate encodings for chromosome X
        if not chrom in self.anc_file.header.contigs and chrom in ["23", "X"]:
            chrom = "23" if chrom =="X" else "X"
        region = f"{chrom}:{pos}-{pos}"

        variant_ancestry = None
        for rec in self.anc_file.fetch(region=region):
            # fetch returns variants that overlap the target position like indels
            # since ancestry labels are encoded REF/ALT, pysam thinks upstream variants
            # overlap with the requested region. these need to be skipped
            if rec.pos != pos: continue
            variant_ancestry = AncestryEstimates(
                rec.chrom,
                rec.pos,
                [rec.ref] + list(rec.alts)
            )

            for sample in rec.samples:
                try:
                    ap = rec.samples[sample]["AP"]
                except KeyError:
                    raise MissingFormatFieldError(
                        f"missing AP FORMAT field from  {self.filepath}"
                    )
                ap = ",".join(ap)
                hap1_probs, hap2_probs = ap.split("|")
                hap1_probs = hap1_probs.split(",")
                hap2_probs = hap2_probs.split(",")
                hap1_probs = np.array(hap1_probs).astype(float)
                hap2_probs = np.array(hap2_probs).astype(float)

                variant_ancestry.samples[sample] = [hap1_probs, hap2_probs]

        return variant_ancestry


    def get_pops(self):
        if self.pops is not None:
            return self.pops
        else:
            rec = next(self.anc_file.fetch())
            self.pops = [rec.ref] + list(rec.alts)
            return self.pops


    def get_samples(self):
        if self.samples is not None:
            return self.samples
        else:
            self.samples = self.anc_file.header.samples
            return self.samples


    def close(self):
        self.anc_file.close()
