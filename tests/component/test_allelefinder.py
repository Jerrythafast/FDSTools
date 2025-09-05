from io import StringIO
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestAllelefinder(ComponentTestCase):
    def run_default(self, infile: str, name: str):
        expected_output = (self.data_dir_root / Tools.ALLELEFINDER /
                           f"{Tools.ALLELEFINDER}_{name}.txt")
        expected_report = (self.data_dir_root / Tools.ALLELEFINDER /
                           f"{Tools.ALLELEFINDER}_{name}-report.txt")

        with open(infile) as infile,\
                patch("sys.stdin", infile),\
                patch("sys.stdout", StringIO()) as generated_output,\
                patch("sys.stderr", StringIO()) as generated_report:
            test_data = [OutputFile(expected_output, generated_output),
                         OutputFile(expected_report, generated_report)]
            self.componentTest(Tools.ALLELEFINDER, [], test_data)

    def test_default_param(self):
        """Also tests single in-/output using stdin, stdout and stderr."""
        self.run_default(
            str(self.data_dir_root / Tools.STUTTERMARK / "default" / "F0310-stuttermark.out"),
            "default_accepted")
        self.run_default(
            str(self.data_dir_root / Tools.STUTTERMARK / "default" / "F0107-stuttermark.out"),
            "default_rejected")

    def test_advanced1_param(self):
        """Also tests multi input using regex.
           Output is used as input for integration tests of other tools."""
        infiles = str(self.data_dir_root / Tools.STUTTERMARK / "default" / "F?????stuttermark.out")

        output = OutputFile(self.data_dir_root / Tools.ALLELEFINDER / "allelefinder_advanced1.txt")
        report = OutputFile(self.data_dir_root / Tools.ALLELEFINDER /
                           "allelefinder_advanced1-report.txt")
        test_data = [output, report]

        fdstools_args = [
            infiles,
            "--output", str(output.generated_output),
            "--report", str(report.generated_output),
            "--tag-expr", r"(F\d\d\d\d)",
            "--max-noisy", "1"
        ]

        self.componentTest(Tools.ALLELEFINDER, fdstools_args, test_data)

    def test_advanced2_param(self):
        """Also tests multi input naming specific files."""
        infiles = [str(self.data_dir_root / Tools.STUTTERMARK / "default" / "F0107-stuttermark.out"),
                   str(self.data_dir_root / Tools.STUTTERMARK / "default" / "F0108-stuttermark.out")]
        library = "ForenSeqA"

        output = OutputFile(self.data_dir_root / Tools.ALLELEFINDER / "allelefinder_advanced2.txt")
        report = OutputFile(self.data_dir_root / Tools.ALLELEFINDER /
                            "allelefinder_advanced2-report.txt")
        test_data = [output, report]

        fdstools_args = (
                infiles +
                         [
                             "--output", str(output.generated_output),
                             "--report", str(report.generated_output),
                             "--tag-expr", r"(F\d\d\d\d)",
                             "--library", library,
                             "--min-allele-pct", "50",
                             "--max-noise-pct", "5",
                             "--min-reads", "100",
                             "--min-reads-lowest", "20",
                             "--max-alleles", "1",
                             "--max-noisy", "5",
                             "--sequence-format", "raw"
                         ]
        )

        self.componentTest(Tools.ALLELEFINDER, fdstools_args, test_data)
