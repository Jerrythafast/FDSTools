from io import StringIO
from pathlib import Path
from typing import List
from unittest import skipIf
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

try:
    from pysam import AlignmentFile
    skip_tests = False
except ImportError:
    try:
        from bamnostic import AlignmentFile
        skip_tests = False
    except ImportError:
        skip_tests = True


@skipIf(skip_tests, "Skip because neither pysam (unix systems) nor bamnostic is installed")
class TestBam2Fds(ComponentTestCase):
    library = "ID-OmniSTR"

    def test_default_param(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_Paired.bam"
        expected_output = (self.data_dir_root / Tools.BAM2FDS /
                           "OmniSTR_Mixture_Paired_bam2fds_default.txt")
        fdstools_args = [self.library, str(infile)]
        with patch("sys.stdout", StringIO()) as generated_output:
            test_data = [OutputFile(expected_output, generated_output)]
            self.componentTest(Tools.BAM2FDS, fdstools_args, test_data)

    def test_advanced1_param(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_Paired.bam"
        output = OutputFile(self.data_dir_root / Tools.BAM2FDS /
                            "OmniSTR_Mixture_Paired_bam2fds_advanced1.txt")
        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--num-threads", "2",
            "--minimum", "10"
        ]
        test_data = [output]
        self.componentTest(Tools.BAM2FDS, fdstools_args, test_data)

    def test_advanced2_param(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R1.bam"
        output = OutputFile(self.data_dir_root / Tools.BAM2FDS /
                            "OmniSTR_Mixture_R1_4lt2_bam2fds_advanced2.txt")
        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--num-threads", "2",
            "--sequence-format", "allelename",
            "--minimum", "85",
            "--missing-marker-action", "exclude",
            "--no-aggregate-filtered"
        ]
        test_data = [output]
        self.componentTest(Tools.BAM2FDS, fdstools_args, test_data)

    def test_advanced3_param(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R1.bam"
        output = OutputFile(self.data_dir_root / Tools.BAM2FDS /
                            "OmniSTR_Mixture_R1_4lt2_bam2fds_advanced3.txt")
        fdstools_args = [
            self.library,
            str(infile),
            str(output.generated_output),
            "--num-threads", "2",
            "--minimum", "5",
            "--no-supplementary"
        ]
        test_data = [output]
        self.componentTest(Tools.BAM2FDS, fdstools_args, test_data)

    def test_halt_error(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R1.bam"
        fdstools_args = [
            self.library,
            str(infile),
            "--missing-marker-action", "halt",
            "--minimum", "85",
            "--no-aggregate-filtered"
        ]
        self.componentTestException(
            Tools.BAM2FDS, fdstools_args,
            SystemExit, 2,
            test_msg="Marker D19S433 was not detected!"
        )

    def test_halt_no_error(self):
        infile = self.data_dir_root / "_reads" / "OmniSTR_Mixture_R1.bam"
        fdstools_args = [
            self.library,
            str(infile),
            "--missing-marker-action", "halt"
        ]
        self.componentTest(Tools.BAM2FDS, fdstools_args, [])
