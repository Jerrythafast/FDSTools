from io import StringIO
from typing import Tuple
from unittest.mock import patch
from os import listdir

from tests.lib.ComponentTestCase import ComponentTestCase, Tools

from tests.lib.ComponentTestData import OutputFile


def expected_mismatches() -> Tuple[str, str, str, str]:
    # For AGAT D6S1043, lbound values swap between +1 and -1 stutter.
    msg_a = r"AGAT\tD6S1043\t\+1\t%i\t6\t46\ttotal\t0.992\t1.481e-02\t-4.585e-03\t3.527e-04\n"
    msg1 = msg_a % 7
    msg2 = msg_a % 8
    msg_b = r"AGAT\tD6S1043\t-1\t%i\t6\t46\ttotal\t0.992\t1.742e-01\t-5.393e-02\t4.149e-03\n"
    msg3 = msg_b % 7
    msg4 = msg_b % 8

    return msg1, msg2, msg3, msg4


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
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1.txt")
        output.accepted_mismatch_tuple = expected_mismatches()
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced1-raw.txt")

        fdstools_args = [self.infiles, "--output", str(output.generated_output),
                         "--raw-outfile", str(raw.generated_output), "--combine-strands",
                         "--min-lengths", "3"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced2_param(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced2-raw.txt")

        fdstools_args = ([self.infiles, "--output", str(output.generated_output),
                         "--raw-outfile", str(raw.generated_output), "--orphans",
                         "--num-threads", "2", "--degree", "1", "--same-shape", "--min-pct",
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

        fdstools_args = [self.infiles, "--output", str(output.generated_output),
                         "--raw-outfile", str(raw.generated_output), "--combine-strands",
                         "--min-lengths", "3", "--max-unit-length", "3"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_advanced5_param(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced5.txt")
        output.accepted_mismatch_tuple = expected_mismatches()[:4]
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_advanced5-raw.txt")

        fdstools_args = [self.infiles, "--output", str(output.generated_output), "--raw-outfile",
                         str(raw.generated_output), "--combine-strands", "--min-lengths", "3",
                         "--min-samples", "2"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)

    def test_marker(self):
        output = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_marker.txt")
        raw = OutputFile(self.data_dir_root / Tools.STUTTERMODEL / "stuttermodel_marker-raw.txt")

        fdstools_args = [self.infiles, "--output", str(output.generated_output), "--raw-outfile",
                         str(raw.generated_output), "--combine-strands", "--min-lengths", "3",
                         "--marker", "D13S317"] + self.base_args

        test_data = [output, raw]
        self.componentTest(Tools.STUTTERMODEL, fdstools_args, test_data)
