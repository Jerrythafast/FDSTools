from io import StringIO
from unittest.mock import patch
from os import listdir
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

class TestBGEstimate(ComponentTestCase):

    def setUp(self) -> None:
        super().setUp() 

        self.indir = self.data_dir_root / "_references"
        self.infiles = str(self.indir / "F????.txt")
        self.infile_list = [str(self.indir / infile) for infile in listdir(self.indir)]

        self.library = "ForenSeqA"
        self.allelelist = str(self.data_dir_root / Tools.ALLELEFINDER / "allelefinder_advanced1.txt")
        self.tag_expr = "(F....)"

    def test_default_param(self):
        """
        Single infiles in batch
        """
        expected_output = (self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_default.txt")

        fdstools_args = self.infile_list + [
          "--allelelist", self.allelelist,
          "--library", self.library,
          "--tag-expr", self.tag_expr]

        with patch("sys.stdout", StringIO()) as generated_output:
            test_data = [OutputFile(expected_output, generated_output)]
            self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

    def test_advanced1_param(self):
        """
        Single infiles in batch with regex for input files
        """
        output = OutputFile(self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_default.txt")

        fdstools_args = [
          self.infiles,
          "--output", str(output.generated_output),
          "--allelelist", self.allelelist,
          "--library", self.library,
          "--tag-expr", self.tag_expr
        ]

        test_data = [output]
        self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

    def test_advanced2_param(self):
        """
        -> Arguments differ from default values.
        -> --min-genotypes, "4",
        -> --combine-strand turned on
        -> --min-output-pct, "0"
        -> Output is used as input for component tests of seqconvert and bgmerge.
        """
        output = OutputFile(self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_adv2.txt")

        fdstools_args = self.infile_list + [
            "--output", str(output.generated_output),
            "--allelelist", self.allelelist,
            "--library", self.library,
            "--tag-expr", self.tag_expr,
            "--min-pct", "0",
            "--min-abs", "1",
            "--min-samples", "1",
            "--min-sample-pct", "50",
            "--min-genotypes", "4",
            "--min-output-pct", "0",
            "--combine-strands"
        ]

        test_data = [output]
        self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

    def test_advanced3_param(self):
        """
        --profiles
        """
        start_profile = self.data_dir_root / Tools.BGESTIMATE / "start_profile.txt"
        output = OutputFile(self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_default.txt")

        fdstools_args = [
          self.infiles,
          "--output", str(output.generated_output),
          "--allelelist", self.allelelist,
          "--library", self.library,
          "--tag-expr", self.tag_expr,
          "--profiles", str(start_profile)
        ]

        test_data = [output]
        self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

    def test_advanced4_param(self):
        """
        -> --marker
        -> --min-genotypes, "2",
        """
        output = OutputFile(self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_adv4.txt")

        fdstools_args = [
            self.infiles,
            "--output", str(output.generated_output),
            "--allelelist", self.allelelist,
            "--library", self.library,
            "--tag-expr", self.tag_expr,
            "--min-pct", "0",
            "--min-abs", "1",
            "--min-samples", "1",
            "--min-sample-pct", "50",
            "--min-genotypes", "2",
            "--min-output-pct", "0",
            "--marker", 'D18S51',
            "--combine-strands"
        ]

        test_data = [output]
        self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

    def test_advanced5_param(self):
        """
        --homozygotes turned on
        -> --min-genotypes, "1"
        -> --min-output-pct", "0.6"
        """
        output = OutputFile(self.data_dir_root / Tools.BGESTIMATE /
                           "bgestimate_adv5.txt")

        fdstools_args = self.infile_list + [
            "--output", str(output.generated_output),
            "--allelelist", self.allelelist,
            "--library", self.library,
            "--tag-expr", self.tag_expr,
            "--min-pct", "0",
            "--min-abs", "1",
            "--min-samples", "1",
            "--min-sample-pct", "50",
            "--min-genotypes", "1",
            "--min-output-pct", "0.6",
            "--combine-strands",
            "--homozygotes"
        ]

        test_data = [output]
        self.componentTest(Tools.BGESTIMATE, fdstools_args, test_data)

