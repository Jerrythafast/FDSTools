import os
from io import StringIO
from shutil import copytree
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestStuttermark(ComponentTestCase):
    def test_default_param(self):
        """Also tests multi in-/output using regex.
           Output is used as input for integration tests of other tools."""
        # repo>tempdir, because output is written to input location
        copytree(self.data_dir_root / "_references", os.getcwd(), dirs_exist_ok=True)

        outdir_expected = self.data_dir_root / Tools.STUTTERMARK / "default"
        test_data = [OutputFile(outdir_expected / e) for e in os.listdir(outdir_expected)]

        fdstools_args = ["--input", "F*.txt"]
        self.componentTest(Tools.STUTTERMARK, fdstools_args, test_data)

    def test_advanced1_param(self):
        """Also tests multi in-/output naming specific files."""
        infiles = [str(self.data_dir_root / "_references" / "F0107.txt"),
                   str(self.data_dir_root / "_references" / "F0809.txt")]

        outdir_expected = self.data_dir_root / Tools.STUTTERMARK / "advanced1"
        test_data = [OutputFile(outdir_expected / e) for e in os.listdir(outdir_expected)]

        fdstools_args = (["--input"] + infiles + ["--output", str(test_data[0].generated_output),
                          str(test_data[1].generated_output), "--stutter=-1:10,+1:8",
                          "--min-reads", "1", "--min-repeats", "2", "--min-report", "2"])
        self.componentTest(Tools.STUTTERMARK, fdstools_args, test_data)

    def test_library_required(self):
        """Also tests single in-/output using stdin and stdout.
           library is required because the input contains raw sequences and stuttermark requires tssv format."""
        infile = self.data_dir_root / Tools.TSSV / "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt"

        expected_output = self.data_dir_root / "stuttermark" / "stuttermark_library_required.txt"

        library = "ID-OmniSTR"
        fdstools_args = ["--library", library]

        with open(infile) as infile,\
                patch("sys.stdin", infile),\
                patch("sys.stdout", StringIO()) as generated_output:
            test_data = [OutputFile(expected_output, generated_output)]
            self.componentTest(Tools.STUTTERMARK, fdstools_args, test_data)

    def test_library_required_error(self):
        """library is required because the input contains raw sequences and stuttermark requires tssv format."""
        infile = self.data_dir_root / Tools.TSSV / "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt"

        with open(infile) as infile,\
                patch("sys.stdin", infile):
            self.componentTestException(
                Tools.STUTTERMARK, [], SystemExit,2,
                test_msg="ValueError: Sequence needs to be converted from raw to tssv, this conversion "
                    "requires a library file")
