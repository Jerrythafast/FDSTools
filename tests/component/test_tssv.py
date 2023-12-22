from io import StringIO
from pathlib import Path
from unittest import skipUnless
from unittest.mock import patch

from tests.lib.FDSToolsComponentTestCase import FDSToolsComponentTestCase
from tests.lib.FDSToolsTestCase import run_long_tests


def get_outfile_names(outfile, report, outdir):
    dir_file1 = outdir + "/statistics.csv"
    dir_file2 = outdir + "/sequences.csv"
    dir_file3 = outdir + "/AmelogeninX/sequences.csv"
    dir_file4_expected = outdir + "/AmelogeninX/noend.fa"
    dir_file4_generated = outdir + "/AmelogeninX/noend.fa.gz"

    return {outfile: outfile,
            report: report,
            dir_file1: dir_file1,
            dir_file2: dir_file2,
            dir_file3: dir_file3,
            dir_file4_expected: dir_file4_generated}


class Test(FDSToolsComponentTestCase):
    def test_default_param(self):
        library = "ForenSeqA"
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "ForenseqAMixture_R1_tssv_default.txt"
        report = "ForenseqAMixture_R1_tssv_default-report.txt"

        with open(infile) as infile:
            with patch("sys.stdin", infile):
                with patch("sys.stdout", StringIO()) as outfile_generated:
                    with patch("sys.stderr", StringIO()) as report_generated:
                        outfile_names = {outfile: outfile_generated,
                                         report: report_generated}
                        fdstools_args = [library]
                        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected, outfile_names, True)
    # test_default_param

    def test_advanced1_param(self):
        """Output is used as input for component tests of other tools."""
        library = "ForenSeqA"
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "ForenseqAMixture_R1_tssv_advanced1.txt"
        report = "ForenseqAMixture_R1_tssv_advanced1-report.txt"
        outdir = "dir_advanced1"
        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "-R", report, "-D", outdir, "-T", "2",
                         "-a", "10"]
        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)

        self.subTestExpectedSystemFileObjectCount(Path(outdir), 154, 619)
    # test_advanced1_param

    def test_advanced2_param(self):
        library = str(self.data_dir / "_libraries" / "ForenSeqA_with_flanks.ini")
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "ForenseqAMixture_R1_tssv_advanced2.txt"
        report = "ForenseqAMixture_R1_tssv_advanced2-report.txt"
        outdir = "dir_advanced2"
        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "-D", outdir, "-R", report, "-L", "10", "-T", "2",
                         "-F", "allelename", "-n", "3", "-a", "3", "-M", "exclude"]

        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_advanced2_param

    @skipUnless(run_long_tests(), "Test takes longer than 30 seconds.")
    def test_advanced3_param(self):
        library = "ForenSeqA"
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "ForenseqAMixture_R1_tssv_advanced3.txt"
        report = "ForenseqAMixture_R1_tssv_advanced3-report.txt"
        outdir = "dir_advanced3"

        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "-D", outdir, "-R", report, "-T", "2", "-X",
                         "-m", "4", "-n", "1", "-a", "5", "-B"]

        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_advanced3_param

    def test_halt_error(self):
        library = "ForenSeqA"
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        fdstools_args = [library, str(infile), "-M", "halt", "-a", "10", "-B"]

        self.subTestToolRaisesException("tssv", fdstools_args, SystemExit,
                                        2, msg="Marker DXS10103 was not detected!")
    # test_halt

    def test_halt_no_error(self):
        library = "ForenSeqA"
        infile = self.data_dir / "fasta" / "ForenseqAMixture_R1.fasta.gz"

        fdstools_args = [library, str(infile), "-M", "halt"]

        self.subTestToolRaisesException("tssv", fdstools_args, None,
                                        2, msg="No halt exception")
    # test_halt
# Test
