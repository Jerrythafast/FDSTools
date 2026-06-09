from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class SampleStats(ComponentTestCase):

    def setUp(self) -> None:
        super().setUp()

        self.indir = self.data_dir_root / "_references"
        self.infiles = str(self.data_dir_root / "_references" / "F010?.txt")


    def test_default_param(self):
        """
        -> Single infile and outfile from standard streams
        -> Output is used as input for component tests of seqconvert and findnewalleles.
        -> Output is used (in adapted form) as input for component tests of test_advanced4_param.
        """
        infile = self.data_dir_root / "_references" / "F0107.txt"
        expected_outfile =  self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt"

        fdstools_args = [
        ]

        with open(infile) as infile,\
            patch("sys.stdin", infile), \
            patch("sys.stdout", StringIO()) as generated_outstream:
            test_data = [OutputFile(expected_outfile, generated_outstream)]
            self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)


    def test_advanced1_param(self):
        """
        Single infile and outfile from positional argument
        """
        infile = self.data_dir_root / "_references" / "F0107.txt"
        outfile =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt")

        fdstools_args = [
         str(infile),
         str(outfile.generated_output)
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)


    def test_advanced2_param(self):
        """
        -> Single infiles in batch
        -> Output is used as input for component tests of seqconvert and findnewalleles.
        """
        infile_1 = self.data_dir_root / "_references" / "F0107.txt"
        infile_2 = self.data_dir_root / "_references" / "F0108.txt"
        infile_3 = self.data_dir_root / "_references" / "F0109.txt"

        outfile_F0107 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt")
        outfile_F0108 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0108.txt")
        outfile_F0109 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0109.txt")

        fdstools_args = [
        "--input", str(infile_1), str(infile_2), str(infile_3),
        "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
        ]

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)

    def test_advanced3_param(self):
        """
        Single infiles in batch with regex for input files
        """
        infiles = self.data_dir_root / "_references" / "F010?.txt"
        outfile_F0107 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt")
        outfile_F0108 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0108.txt")
        outfile_F0109 =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0109.txt")

        fdstools_args = [
        "--input", str(infiles),
        "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
        "--tag-expr", "(F010.)"
        ]

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)

    def test_advanced4_param(self):
        """
        -> Changed interpretation settings to differ from default values.
        -> --settings
        -> samplestats_default_F0107_adapated_for_input_adv4-8.txt is an adapted version of
        samplestats_default_F0107, and samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt
        is the bgcorrected version of this (with bgmerge_default_old_20260601.txt as profile).
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS/ "samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt"
        outfile = OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_adv4_F0107.txt")
        settings = self.data_dir_root / Tools.SAMPLESTATS / "settings.txt"

        fdstools_args = [
            "--input", str(infile),
            "--output", str(outfile.generated_output),
            "--settings", str(settings),
            "--min-reads", "10",
            "--min-per-strand", "1",
            "--min-pct-of-max", "1",
            "--min-pct-of-sum", "1",
            "--min-correction", "0",
            "--min-recovery", "0",
            "--min-allele-reads", "1",
            "--max-nonallele-pct", "90",
            "--max-alleles", "3"
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)


    def test_advanced5_param(self):
        """
        --uncall-alleles turned on
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS/ "samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt"
        outfile = OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_adv5_F0107.txt")

        fdstools_args = [
            "--input", str(infile),
            "--output", str(outfile.generated_output),
            "--uncall-alleles",
            "--min-reads", "10",
            "--min-per-strand", "1",
            "--min-pct-of-max", "1",
            "--min-pct-of-sum", "1",
            "--min-correction", "0",
            "--min-recovery", "0",
            "--min-allele-reads", "1",
            "--max-nonallele-pct", "95",
            "--max-alleles", "3"
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)

    def test_advanced6_param(self):
        """
        --uncall-alleles = only
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS/ "samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt"
        outfile =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_adv6_F0107.txt")

        fdstools_args = [
            "--input", str(infile),
            "--output", str(outfile.generated_output),
            "--uncall-alleles", "only",
            "--min-reads", "10",
            "--min-per-strand", "1",
            "--min-pct-of-max", "1",
            "--min-pct-of-sum", "1",
            "--min-correction", "0",
            "--min-recovery", "0",
            "--min-allele-reads", "1",
            "--max-nonallele-pct", "95",
            "--max-alleles", "3"
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)

    def test_advanced7_param(self):
        """
        -> --filter-action = delete
        -> --filter-absolute turned on
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS/ "samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt"
        outfile =  OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_adv7_F0107.txt")

        fdstools_args = [
            "--input", str(infile),
            "--output", str(outfile.generated_output),
            "--filter-action", "delete",
            "--min-reads-filt", "6",
            "--min-per-strand-filt", "1",
            "--min-pct-of-max-filt", "1",
            "--min-pct-of-sum-filt", "1",
            "--min-correction-filt", "0",
            "--min-recovery-filt", "0"
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)

    def test_advanced8_param(self):
        """
        --filter-action = combine
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS/ "samplestats_default_F0107_bgcorrect_for_input_adv4-8.txt"
        outfile = OutputFile(self.data_dir_root / Tools.SAMPLESTATS / "samplestats_adv8_F0107.txt")

        fdstools_args = [
            "--input", str(infile),
            "--output", str(outfile.generated_output),
            "--filter-action", "combine",
            "--min-reads-filt", "6",
            "--min-per-strand-filt", "1",
            "--min-pct-of-max-filt", "1",
            "--min-pct-of-sum-filt", "1",
            "--min-correction-filt", "0",
            "--min-recovery-filt", "0"
        ]

        test_data = [outfile]
        self.componentTest(Tools.SAMPLESTATS, fdstools_args, test_data)
