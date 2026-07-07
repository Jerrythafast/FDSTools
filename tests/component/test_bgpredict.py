from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestBGPredict(ComponentTestCase):

    def setUp(self) -> None:
        super().setUp()

        self.infile = str(self.data_dir_root / "_references" / "F0107.txt")
        self.library = "ForenSeqA"

    def test_default_param(self):
        """
        Single positional infile and outfile to stdout
        """

        infile_stuttermodel = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2.txt"
        expected_output = (self.data_dir_root / Tools.BGPREDICT /
                           "bgpredict-default.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            "--library", self.library
        ]

        with patch("sys.stdout", StringIO()) as generated_output:
            test_data = [OutputFile(expected_output, generated_output)]
            self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)

    def test_advanced1_param(self):
        """
        Single infile and outfile from positional arguent
        """
        infile_stuttermodel = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2.txt"
        output = OutputFile(self.data_dir_root / Tools.BGPREDICT /
                            "bgpredict-default.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            str(output.generated_output),
            "--library", self.library
        ]

        test_data = [output]
        self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)

    def test_advanced2_param(self):
        """
        -> --combine-strand turned on
        -> Output is used (in adapted form) as input for component tests of seqconvert and bgmerge.
        """
        infile_stuttermodel = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1.txt"
        output = OutputFile(self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv2.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            str(output.generated_output),
            "--library", self.library,
            "--combine-strands"
        ]

        test_data = [output]
        self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)

    def test_advanced3_param(self):
        """
        -> --min-pct and --min-r2 differ from default values.
        -> Output is used (in adapted form) as input for component tests of seqconvert and bgmerge.
        """
        infile_stuttermodel = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1.txt"
        output = OutputFile(self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv3.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            str(output.generated_output),
            "--combine-strands",
            "--library", self.library,
            "--min-pct", "0.1",
            "--min-r2", "0.7"
        ]

        test_data = [output]
        self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)

    def test_advanced4_param(self):
        """
        -> --use-all-data turned on
        -> Output is used (in adapted form) as input for component tests of seqconvert and bgmerge.
        """
        infile_stuttermodel = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1.txt"
        output = OutputFile(self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv4.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            str(output.generated_output),
            "--combine-strands",
            "--use-all-data",
            "--library", self.library,
            "--min-pct", "0.1",
            "--min-r2", "0.7"
        ]

        test_data = [output]
        self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)

    def test_advanced5_param(self):
        """
        -> Only in this test the absense of --no-check-lbounds matters!!
        """
        infile_stuttermodel = self.data_dir_root / Tools.BGPREDICT / "stuttermodel_adapted_adv5.txt" # adapted from an old version of stuttermodel-adv1
        output = OutputFile(self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv5.txt")

        fdstools_args = [
            str(infile_stuttermodel),
            self.infile,
            str(output.generated_output),
            "--combine-strands",
            "--use-all-data",
            "--library", self.library,
            "--min-pct", "0.9",
        ]

        test_data = [output]
        self.componentTest(Tools.BGPREDICT, fdstools_args, test_data)
