from io import StringIO
from pathlib import Path
from typing import List
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

class TestTssv(ComponentTestCase):
    library = "ID-OmniSTR"

    def get_test_data(self, outfile: OutputFile, report: OutputFile, outdir: str) \
            -> List[OutputFile]:
        dir_file1 = outdir + "/statistics.csv"
        dir_file2 = outdir + "/sequences.csv"
        dir_file3 = outdir + "/AmelogeninX/sequences.csv"
        dir_file4_expected = outdir + "/AmelogeninX/noend.fa"
        dir_file4_generated = outdir + "/AmelogeninX/noend.fa.gz"

        return [outfile, report,
                OutputFile(self.data_dir_root / Tools.TSSV / dir_file1, Path(dir_file1)),
                OutputFile(self.data_dir_root / Tools.TSSV / dir_file2, Path(dir_file2)),
                OutputFile(self.data_dir_root / Tools.TSSV / dir_file3, Path(dir_file3)),
                OutputFile(self.data_dir_root / Tools.TSSV / dir_file4_expected,
                           Path(dir_file4_generated))]

    def test_default_param(self):
        """Output is used as input for component tests of seqconvert."""
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        expected_output = (self.data_dir_root / Tools.TSSV /
                           "OmniSTR_Mixture_R1_4lt2_tssv_default.txt")
        expected_report = (self.data_dir_root / Tools.TSSV /
                           "OmniSTR_Mixture_R1_4lt2_tssv_default-report.txt")

        fdstools_args = [self.library]

        with open(infile) as infile,\
                patch("sys.stdin", infile), \
                patch("sys.stdout", StringIO()) as generated_output, \
                patch("sys.stderr", StringIO()) as generated_report:

            test_data = [OutputFile(expected_output, generated_output),
                         OutputFile(expected_report, generated_report)]

            self.componentTest(Tools.TSSV, fdstools_args, test_data)

    def test_advanced1_param(self):
        """Output is used as input for component tests of stuttermark."""
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        output = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced1.txt")
        report = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced1-report.txt")
        outdir = "dir_advanced1"

        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--report", str(report.generated_output),
            "--dir", outdir,
            "--num-threads", "2",
            "--minimum", "10"
        ]

        test_data = self.get_test_data(output, report, outdir)
        self.componentTest(Tools.TSSV, fdstools_args, test_data)
        self.componentTestSystemFileCount(
            Path(outdir),
            expected_dir_count=31,
            expected_file_count=127
        )

    def test_advanced2_param(self):
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        output = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced2.txt")
        report = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced2-report.txt")
        outdir = "dir_advanced2"

        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--dir", outdir,
            "--report", str(report.generated_output),
            "--flank-length", "10",
            "--num-threads", "2",
            "--sequence-format", "allelename",
            "--indel-score", "3",
            "--minimum", "3",
            "--missing-marker-action", "exclude"
        ]

        test_data = self.get_test_data(output, report, outdir)
        self.componentTest(Tools.TSSV, fdstools_args, test_data)

    def test_advanced3_param(self):
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        output = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced3.txt")
        report = OutputFile(self.data_dir_root / Tools.TSSV /
                            "OmniSTR_Mixture_R1_4lt2_tssv_advanced3-report.txt")
        outdir = "dir_advanced3"

        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--dir", outdir,
            "--report", str(report.generated_output),
            "--no-deduplicate",
            "--flank-length", "20",
            "--num-threads", "2",
            "--mismatches", "4",
            "--indel-score", "1",
            "--minimum", "5",
            "--no-aggregate-filtered"
        ]

        test_data = self.get_test_data(output, report, outdir)
        self.componentTest(Tools.TSSV, fdstools_args, test_data)

    def test_halt_error(self):
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        fdstools_args = [
            self.library,
            str(infile),
            "--missing-marker-action", "halt",
            "--minimum", "85",
            "--no-aggregate-filtered"
        ]

        self.componentTestException(
            Tools.TSSV, fdstools_args,
            SystemExit,2,
            test_msg="Marker D19S433 was not detected!"
        )

    def test_halt_no_error(self):
        infile = self.data_dir_root / "fasta" / "OmniSTR_Mixture_R1_4lt2.fasta.gz"

        fdstools_args = [
            self.library,
            str(infile),
            "--missing-marker-action", "halt"
        ]

        self.componentTest(Tools.TSSV, fdstools_args, [])
