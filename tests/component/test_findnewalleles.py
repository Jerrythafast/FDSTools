from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestFindNewAlleles(ComponentTestCase):

    def setUp(self) -> None:
        super().setUp()

        self.library = "ForenSeqA"
        self.known = str(self.data_dir_root / Tools.FINDNEWALLELES / "known.txt")

    def test_default_param(self):
        """
        Single infile and outfile from standard streams
        """

        infile = self.data_dir_root / Tools.SAMPLESTATS  / "samplestats_default_F0107.txt"
        expected_outfile =  self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_default.txt"
        fdstools_args = [
            self.known,
            "--library", self.library
        ]

        with open(infile) as infile,\
            patch("sys.stdin", infile), \
            patch("sys.stdout", StringIO()) as generated_outstream:
            test_data = [OutputFile(expected_outfile, generated_outstream)]
            self.componentTest(Tools.FINDNEWALLELES, fdstools_args, test_data)


    def test_advanced2_param(self):
        """
        Single infile and outfile from positional argument
        """
        infile = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt"
        outfile = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_default.txt")

        fdstools_args = [
            self.known,
            str(infile),
            str(outfile.generated_output),
            "--library", self.library
        ]

        test_data = [outfile]
        self.componentTest(Tools.FINDNEWALLELES, fdstools_args, test_data)

    def test_advanced3_param(self):
        """
        Single infiles in batch
        """
        infile_1 = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt"
        infile_2 = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0108.txt"
        infile_3 = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0109.txt"

        outfile_F0107 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0107.txt")
        outfile_F0108 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0108.txt")
        outfile_F0109 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0109.txt")

        fdstools_args = [
            self.known,
            "--input", str(infile_1), str(infile_2), str(infile_3),
            "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
            "--library", self.library
        ]

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.FINDNEWALLELES, fdstools_args, test_data)

    def test_advanced4_param(self):
        """
        Single infiles in batch with regex for input files
        """

        infiles = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F010?.txt"
        outfile_F0107 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0107.txt")
        outfile_F0108 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0108.txt")
        outfile_F0109 = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv3_F0109.txt")

        fdstools_args = [
            self.known,
            "--input", str(infiles),
            "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
            "--library", self.library,
            "--tag-expr", "(F010.)"
        ]

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.FINDNEWALLELES, fdstools_args, test_data)

    def test_advanced5_param(self):
        """
        --remove-allele-flags turned on: should remove all allele tags except of the
        first allele of each marker, as those are present in the known.txt (first alleles)
        """
        infile = self.data_dir_root / Tools.SAMPLESTATS / "samplestats_default_F0107.txt"
        outfile = OutputFile(self.data_dir_root / Tools.FINDNEWALLELES / "newalleles_adv5_F0107.txt")

        fdstools_args = [
            self.known,
            str(infile),
            str(outfile.generated_output),
            "--library", self.library,
            "--remove-allele-flags"
        ]

        test_data = [outfile]
        self.componentTest(Tools.FINDNEWALLELES, fdstools_args, test_data)

