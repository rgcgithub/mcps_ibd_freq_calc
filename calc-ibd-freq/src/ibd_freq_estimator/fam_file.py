class FamFile:
    """Parse fam file to store the sex of each sample.

    Parameters
    ----------
        fam_file : str
            Path to plink fam file
    """

    def __init__(self, fam_file):
        self.sample_sex = self.load_sex(fam_file)

    def load_sex(self, fam_file):
        sample_sex = {}
        with open(fam_file) as f:
            for line in f:
                fid, iid, _, _, sex, _ = line.strip("\n").split()
                if sex == "1":
                    sample_sex[iid] = "male"
                elif sex == "2":
                    sample_sex[iid] = "female"
                elif sex == "0":
                    sample_sex[iid] = "unknown"
                else:
                    raise ValueError("bad fam file (sex should be 0, 1, or 2)")
        return sample_sex

    def is_male(self, sample_id):
        return self.sample_sex[sample_id] == "male"

    def is_female(self, sample_id):
        return self.sample_sex[sample_id] == "female"

    def is_missing(self, sample_id):
        return self.sample_sex[sample_id] == "unknown"
