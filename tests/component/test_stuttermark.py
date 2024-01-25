import os
import sys
from io import StringIO
from shutil import copytree
from unittest import skipIf
from unittest.mock import patch

from tests.lib.FDSToolsComponentTestCase import FDSToolsComponentTestCase


class Test(FDSToolsComponentTestCase):
    @skipIf(sys.version_info < (3, 8),
            "Skip because shutil.copytree does not support dirs_exist_ok=True in Python < 3.8")
    def test_default_param(self):
        """Also tests multi in-/output using regex.
           Output is used as input for integration tests of other tools."""

        indir = str(self.data_dir / "_references")
        copytree(indir, os.getcwd(), dirs_exist_ok=True)  # input>tempdir, because output is written to input location

        outdir_expected = self.data_dir / "stuttermark" / "default"
        outfile_names = ["F0107-stuttermark.out",
                         "F0108-stuttermark.out",
                         "F0109-stuttermark.out",
                         "F0110-stuttermark.out",
                         "F0207-stuttermark.out",
                         "F0208-stuttermark.out",
                         "F0209-stuttermark.out",
                         "F0210-stuttermark.out",
                         "F0307-stuttermark.out",
                         "F0308-stuttermark.out",
                         "F0309-stuttermark.out",
                         "F0310-stuttermark.out",
                         "F0407-stuttermark.out",
                         "F0408-stuttermark.out",
                         "F0409-stuttermark.out",
                         "F0410-stuttermark.out",
                         "F0507-stuttermark.out",
                         "F0508-stuttermark.out",
                         "F0509-stuttermark.out",
                         "F0510-stuttermark.out",
                         "F0607-stuttermark.out",
                         "F0608-stuttermark.out",
                         "F0609-stuttermark.out",
                         "F0610-stuttermark.out",
                         "F0707-stuttermark.out",
                         "F0708-stuttermark.out",
                         "F0709-stuttermark.out",
                         "F0710-stuttermark.out",
                         "F0807-stuttermark.out",
                         "F0808-stuttermark.out",
                         "F0809-stuttermark.out"]

        fdstools_args = ["--input", "F*.txt"]
        self.subTestToolWorksMultiOutput("stuttermark", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_default_param

    @skipIf(sys.version_info.minor < 8,
            "Skip because shutil.copytree does not support dirs_exist_ok=True in Python < 3.8")
    def test_advanced1_param(self):
        """Also tests multi in-/output naming specific files."""

        indir = str(self.data_dir / "_references")
        copytree(indir, os.getcwd(), dirs_exist_ok=True)  # input>tempdir, because output is written to input location

        outdir_expected = self.data_dir / "stuttermark" / "advanced1"
        outfile_names = ["F0107-stuttermark.out",
                         "F0809-stuttermark.out"]

        fdstools_args = ["--input", "F0107.txt", "F0809.txt", "--output", outfile_names[0], outfile_names[1],
                         "--stutter=-1:10,+1:8", "--min-reads", "1", "--min-repeats", "2", "--min-report", "2"]
        self.subTestToolWorksMultiOutput("stuttermark", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_advanced1_param

    def test_library_required(self):
        """Also tests single in-/output using stdin and stdout.
           library is required because the input contains raw sequences and stuttermark requires tssv format."""
        library = "ID-OmniSTR"
        infile = self.data_dir / "tssv" / "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt"

        outfile_expected = self.data_dir / "stuttermark" / "stuttermark_library_required.txt"

        with open(infile) as infile,\
                patch("sys.stdin", infile),\
                patch("sys.stdout", StringIO()) as outfile_generated:
            fdstools_args = ["--library", library]
            self.assertToolWorks("stuttermark", fdstools_args, outfile_expected,
                                 outfile_generated)
    # test_library_required

    def test_library_required_error(self):
        """library is required because the input contains raw sequences and stuttermark requires tssv format."""
        infile = self.data_dir / "tssv" / "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt"

        with open(infile) as infile,\
                patch("sys.stdin", infile):
            self.subTestToolRaisesException("stuttermark", [], SystemExit, 2, msg="ValueError: Sequence needs to be "
                                            "converted from raw to tssv, this conversion requires a library file")
    # test_library_required_error
