from io import BytesIO, StringIO, TextIOWrapper
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

class TestMergeReads(ComponentTestCase):

    def test_default_param(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_interleaved_4lt2.fastq.gz"

        expected_output = (self.data_dir_root / Tools.MERGEREADS /
                           "OmniSTR_Mixture_4lt2_mergereads_default.fastq")
        expected_report = (self.data_dir_root / Tools.MERGEREADS /
                           "OmniSTR_Mixture_4lt2_mergereads_default-report.txt")

        fdstools_args = []

        with open(infile) as infile,\
                patch("sys.stdin", infile), \
                patch("sys.stdout", TextIOWrapper(BytesIO())) as generated_output:

            test_data = [OutputFile(expected_output, generated_output),
                         OutputFile(expected_report, self.log)]

            self.componentTest(Tools.MERGEREADS, fdstools_args, test_data)

    def test_advanced_param(self):
        infile1 = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R1_4lt2.fastq.gz"
        infile2 = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R2_4lt2.fastq.gz"

        output = OutputFile(self.data_dir_root / Tools.MERGEREADS /
                            "OmniSTR_Mixture_4lt2_mergereads_advanced.fastq.gz")
        report = OutputFile(self.data_dir_root / Tools.MERGEREADS /
                            "OmniSTR_Mixture_4lt2_mergereads_advanced-report.txt", self.log)
        out_r1 = OutputFile(self.data_dir_root / Tools.MERGEREADS /
                            "OmniSTR_Mixture_4lt2_mergereads_advanced-R1.fastq.gz")
        out_r2 = OutputFile(self.data_dir_root / Tools.MERGEREADS /
                            "OmniSTR_Mixture_4lt2_mergereads_advanced-R2.fastq.gz")

        fdstools_args = [
            str(infile1),
            str(infile2),
            "--min-overlap", "30",
            "--lowercase-overhang",
            "--ambiguous",
            "--qual-eq-threshold", "10",
            "--zip-level", "1",
            "--output", str(output.generated_output),
            "--uncombined-r1", str(out_r1.generated_output),
            "--uncombined-r2", str(out_r2.generated_output),
            "--max-mismatch-density", "0.33",
            "--allow-outies"
        ]

        self.componentTest(Tools.MERGEREADS, fdstools_args, [output, report, out_r1, out_r2])
