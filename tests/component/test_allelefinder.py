from io import StringIO
from unittest.mock import patch

from tests.lib.FDSToolsComponentTestCase import FDSToolsComponentTestCase


class Test(FDSToolsComponentTestCase):
    def run_default(self, infile, name):
        outdir_expected = self.data_dir / "allelefinder"
        outfile = "allelefinder_"+name+".txt"
        report = "allelefinder_"+name+"-report.txt"

        with open(infile) as infile,\
                patch("sys.stdin", infile),\
                patch("sys.stdout", StringIO()) as outfile_generated,\
                patch("sys.stderr", StringIO()) as report_generated:
            outfile_names = {outfile: outfile_generated,
                             report: report_generated}
            fdstools_args = []
            self.subTestToolWorksMultiOutput("allelefinder", fdstools_args, outdir_expected, outfile_names)
    # run_default

    def test_default_param(self):
        """Also tests single in-/output using stdin, stdout and stderr."""
        self.run_default(str(self.data_dir / "stuttermark" / "default" / "F0310-stuttermark.out"), "default_accepted")
        self.run_default(str(self.data_dir / "stuttermark" / "default" / "F0107-stuttermark.out"), "default_rejected")
    # test_default_param

    def test_advanced1_param(self):
        """Also tests multi input using regex.
           Output is used as input for integration tests of other tools."""
        infiles = str(self.data_dir / "stuttermark" / "default" / "F?????stuttermark.out")

        outdir_expected = self.data_dir / "allelefinder"
        outfile = "allelefinder_advanced1.txt"
        report = "allelefinder_advanced1-report.txt"
        outfile_names = [outfile, report]

        fdstools_args = [infiles, "--output", outfile, "--report", report, "--tag-expr", "(F....)", "--max-noisy", "1"]
        self.subTestToolWorksMultiOutput("allelefinder", fdstools_args, outdir_expected, outfile_names)
    # test_advanced1_param

    def test_advanced2_param(self):
        """Also tests multi input naming specific files."""
        indir = self.data_dir / "stuttermark" / "default"
        infiles = [str(indir / "F0107-stuttermark.out"),
                   str(indir / "F0108-stuttermark.out")]
        library = "ForenSeqA"

        outdir_expected = self.data_dir / "allelefinder"
        outfile = "allelefinder_advanced2.txt"
        report = "allelefinder_advanced2-report.txt"

        outfile_names = [outfile, report]
        fdstools_args = infiles + ["--output", outfile, "--report", report, "--tag-expr", "(F....)",
                                   "--library", library, "--min-allele-pct", "50", "--max-noise-pct", "5",
                                   "--min-reads", "100", "--min-reads-lowest", "20", "--max-alleles", "1",
                                   "--max-noisy", "5", "--sequence-format", "raw"]
        self.subTestToolWorksMultiOutput("allelefinder", fdstools_args, outdir_expected, outfile_names)
    # test_advanced2_param
