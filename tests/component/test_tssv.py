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
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "OmniSTR_Mixture_R1_4lt2_tssv_default.txt"
        report = "OmniSTR_Mixture_R1_4lt2_tssv_default-report.txt"

        with open(infile) as infile,\
                patch("sys.stdin", infile), \
                patch("sys.stdout", StringIO()) as outfile_generated, \
                patch("sys.stderr", StringIO()) as report_generated:
            outfile_names = {outfile: outfile_generated,
                             report: report_generated}
            fdstools_args = [library]
            self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                             outfile_names, True)
    # test_default_param

    def test_advanced1_param(self):
        """Output is used as input for component tests of other tools."""
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt"
        report = "OmniSTR_Mixture_R1_4lt2_tssv_advanced1-report.txt"
        outdir = "dir_advanced1"
        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "--report", report, "--dir", outdir,
                         "--num-threads", "2", "--minimum", "10"]
        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)

        self.subTestExpectedSystemFileObjectCount(Path(outdir), 31, 127)
    # test_advanced1_param

    def test_advanced2_param(self):
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "OmniSTR_Mixture_R1_4lt2_tssv_advanced2.txt"
        report = "OmniSTR_Mixture_R1_4lt2_tssv_advanced2-report.txt"
        outdir = "dir_advanced2"
        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "--dir", outdir, "--report", report,
                         "--flank-length", "10", "--num-threads", "2",
                         "--sequence-format", "allelename", "--indel-score", "3", "--minimum", "3",
                         "--missing-marker-action", "exclude"]

        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_advanced2_param

    @skipUnless(run_long_tests(), "Test takes longer than 30 seconds.")
    def test_advanced3_param(self):
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        outdir_expected = self.data_dir / "tssv"
        outfile = "OmniSTR_Mixture_R1_4lt2_tssv_advanced3.txt"
        report = "OmniSTR_Mixture_R1_4lt2_tssv_advanced3-report.txt"
        outdir = "dir_advanced3"

        outfile_names = get_outfile_names(outfile, report, outdir)

        fdstools_args = [library, str(infile), outfile, "--dir", outdir, "--report", report,
                         "--num-threads", "2", "--no-deduplicate", "--mismatches", "4",
                         "--indel-score", "1", "--minimum", "5", "--no-aggregate-filtered"]

        self.subTestToolWorksMultiOutput("tssv", fdstools_args, outdir_expected,
                                         outfile_names)
    # test_advanced3_param

    def test_halt_error(self):
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        fdstools_args = [library, str(infile), "--missing-marker-action", "halt", "--minimum", "85",
                         "--no-aggregate-filtered"]

        self.subTestToolRaisesException("tssv", fdstools_args, SystemExit,
                                        2, msg="Marker D19S433 was not detected!")
    # test_halt

    def test_halt_no_error(self):
        library = "ID-OmniSTR"
        infile = self.data_dir / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        fdstools_args = [library, str(infile), "--missing-marker-action", "halt"]

        self.subTestToolRaisesException("tssv", fdstools_args, None,
                                        2, msg="No halt exception")
    # test_halt
# Test
