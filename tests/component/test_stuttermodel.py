from io import StringIO
from unittest.mock import patch
from os import listdir

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestStuttermodel(ComponentTestCase):
    def setUp(self) -> None:
        super().setUp()

        self.indir = self.data_dir_root / "_references"
        self.infiles = str(self.indir / "F????.txt")
        self.infile_list = [str(self.indir / infile) for infile in listdir(self.indir)]

        self.allelelist = str(self.data_dir_root / "allelefinder" / "allelefinder_advanced1.txt")
        self.base_args = ["--tag-expr", "(F....)", "--allelelist", self.allelelist,
                          "--library", "ForenSeqA"]

    def test_default_param(self):
        """
        Also tests multi input naming specific files.
        Also tests output using stdout for --output.
        """
        fdstools_args = self.infile_list + ["--combine-strands"] + self.base_args
        expected_output = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel-default.txt"

        with patch("sys.stdout", StringIO()) as generated_outstream:
            test_data = [OutputFile(expected_output, generated_outstream)]
            self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_raw_param(self):
        """
        Also tests multi input naming specific files.
        Also tests output using stdout for --raw-outfile.
        """
        fdstools_args = self.infile_list + ["--raw-outfile", "-", "--combine-strands"] + self.base_args
        expected_output = self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel-raw.txt"

        with patch("sys.stdout", StringIO()) as generated_output:
            test_data = [OutputFile(expected_output, generated_output)]
            self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced1_param(self):
        """
        Output is used as input for component tests of bgpredict.
        """
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1-raw.txt")

        fdstools_args = [
            self.infiles,
            "--output", str(output.generated_output),
            "--raw-outfile", str(raw.generated_output),
            "--combine-strands",
            "--min-lengths", "3",
            "--nonnegative"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced2_param(self):
        """
        Output is used as input for component tests of bgpredict.
        """
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2-raw.txt")

        fdstools_args = ([self.infiles, "--output", str(output.generated_output),
                         "--raw-outfile", str(raw.generated_output), "--orphans",
                         "--degree", "1", "--same-shape", "--min-pct",
                         "0.2", "--min-abs", "20", "--min-lengths", "6", "--min-r2", "0.1"] +
                         self.base_args)

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced3_param(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced3.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced3-raw.txt")

        fdstools_args = ([self.infiles, "--output", str(output.generated_output), "--raw-outfile",
                         str(raw.generated_output), "--min-lengths", "3", "--ignore-zeros"] +
                         self.base_args)

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced4_param(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced4.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced4-raw.txt")

        fdstools_args = [
            self.infiles,
            "--output", str(output.generated_output),
            "--raw-outfile", str(raw.generated_output),
            "--combine-strands",
            "--min-lengths", "3",
            "--max-unit-length", "3",
            "--min-r2", "0.1",
            "--min-stutter-pct", "0"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced5_param(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced5.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced5-raw.txt")

        fdstools_args = [self.infiles, "--output", str(output.generated_output), "--raw-outfile",
                         str(raw.generated_output), "--combine-strands", "--min-lengths", "3",
                         "--min-samples", "2"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_marker(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_marker.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_marker-raw.txt")

        fdstools_args = [
            self.infiles,
            "--output", str(output.generated_output),
            "--raw-outfile", str(raw.generated_output),
            "--combine-strands",
            "--min-lengths", "3",
            "--min-stutter-pct", "0",
            "--marker", "D13S317"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)
