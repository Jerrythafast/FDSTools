from io import TextIOWrapper
from pathlib import Path
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

class TestDemultiplex(ComponentTestCase):

    def get_test_data(self, expected_dir, outdir=None, double_zero=False):
        if not outdir:
            outdir = expected_dir
        expected_path = self.data_dir_root / Tools.DEMULTIPLEX / expected_dir
        generated_path = Path(outdir)
        filenames = [
            "Undetermined_S0_R1_001.fastq.gz",
            "Undetermined_S0_R2_001.fastq.gz",
            "E0049_S1_R1_001.fastq.gz",
            "E0049_S1_R2_001.fastq.gz",
            "E0150_S2_R1_001.fastq.gz",
            "E0150_S2_R2_001.fastq.gz"]
        if double_zero:
            filenames += ["Undetermined_S00_R1_001.fastq.gz", "Undetermined_S00_R2_001.fastq.gz"]
        return [OutputFile(expected_path / file, generated_path / file) for file in filenames]

    def test_default_param(self):
        infile = self.data_dir_root / Tools.DEMULTIPLEX / "Multiplexed_interleaved.fastq"
        barcodes = self.data_dir_root / Tools.DEMULTIPLEX / "barcodes.txt"
        outdir_expected = self.data_dir_root / Tools.DEMULTIPLEX / "default"
        test_data = self.get_test_data("default", ".")
        fdstools_args = [str(barcodes)]
        with open(infile, "rb") as infile, \
                patch("sys.stdin", TextIOWrapper(infile)):
            self.componentTest(Tools.DEMULTIPLEX, fdstools_args, test_data)

    def test_advanced1_param(self):
        # Testing three equivalent commands in the same test method, to also test --truncate.
        infile1 = self.data_dir_root / Tools.DEMULTIPLEX / "Multiplexed_R1.fastq"
        infile2 = self.data_dir_root / Tools.DEMULTIPLEX / "Multiplexed_R2.fastq"
        barcodes = self.data_dir_root / Tools.DEMULTIPLEX / "SampleSheet.csv"
        outdir = "advanced1"
        output_files = self.get_test_data(outdir)
        fdstools_args_basis = [
            str(barcodes),
            str(infile1),
            str(infile2),
            "--barcode-pattern", r":N:0:([ACGTN]+)\+([ACGTN]+)",
            "--dir", outdir,
            "--truncate",
            "--zip-level", "1",
            "--output-scores",
            "--filter-idx-phix",
            "--filter-diff", "0.1",
            "--filter-score", "0.7"]
        fdstools_args1 = fdstools_args_basis + [
            "--matrix", "xleap",
            "--reverse-right",
            "--artefacts", "GGGGGGGGGG",
            "--right-artefacts", "AGATCTCGGT"]
        fdstools_args2 = fdstools_args_basis + [
            "--matrix", "xleap",
            "--reverse-right",
            "--left-artefacts", "GGGGGGGGGG",
            "--right-artefacts", "GGGGGGGGGG", "AGATCTCGGT"]
        fdstools_args3 = fdstools_args_basis + [
            "--xleap",
            "--right-artefacts", "AGATCTCGGT"]
        self.componentTest(Tools.DEMULTIPLEX, fdstools_args1, output_files)
        self.componentTest(Tools.DEMULTIPLEX, fdstools_args2, output_files)
        self.componentTest(Tools.DEMULTIPLEX, fdstools_args3, output_files)


    def test_advanced2_param(self):
        # Testing two commands in the same test method, to also test appending.
        infile1 = self.data_dir_root / Tools.DEMULTIPLEX / "Multiplexed_R1.fastq"
        infile2 = self.data_dir_root / Tools.DEMULTIPLEX / "Multiplexed_R2.fastq"
        barcodes = self.data_dir_root / Tools.DEMULTIPLEX / "barcodes-reversed-right.txt"
        matrix = self.data_dir_root / Tools.DEMULTIPLEX / "xleap-channels.txt"
        outdir = "advanced2"
        outfiles_first = self.get_test_data("advanced2a", outdir)
        outfiles_second = self.get_test_data("advanced2b", outdir, double_zero=True)

        fdstools_args_shared = [
            "--dir", outdir,
            "--zip-level", "1",
            "--output-scores",
            "--matrix", str(matrix),
            "--artefacts", "GGGGGGGGGG",
            "--right-artefacts", "AGATCTCGGT",
            "--reverse-right"]

        # In the first command, we filter PhiX reads.
        fdstools_args1 = [
            str(barcodes),
            str(infile1),
            str(infile2),
            "--filter-idx-phix"] + fdstools_args_shared

        # The second command takes Undetermined output from the first,
        # and appends the PhiX reads to the output written by the first.
        fdstools_args2 = [
            str(barcodes),
            outdir + "/Undetermined_S0_R1_001.fastq.gz",
            outdir + "/Undetermined_S0_R2_001.fastq.gz"] + fdstools_args_shared

        # First filter phix, then append it to samples.
        self.componentTest(Tools.DEMULTIPLEX, fdstools_args1, outfiles_first)
        self.componentTest(Tools.DEMULTIPLEX, fdstools_args2, outfiles_second)
