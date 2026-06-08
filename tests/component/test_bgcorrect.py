from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestBGCorrect(ComponentTestCase):

    def setUp(self) -> None:
        super().setUp()

        self.infiles = str(self.data_dir_root / "_references" / "F010?.txt")
        self.library = "ForenSeqA"

    def test_default_param(self):
        """
        Single infile and outfile from standard streams
        """
        infile = self.data_dir_root / "_references" / "F0107.txt"
        expected_outfile =  self.data_dir_root / Tools.BGCORRECT / "bgcorrect_default.txt"

        fdstools_args = [
        str(self.data_dir_root / Tools.BGMERGE / "bgmerge_default.txt"),
        '--library', self.library
        ]

        with open(infile) as infile,\
            patch("sys.stdin", infile), \
            patch("sys.stdout", StringIO()) as generated_outstream:
            test_data = [OutputFile(expected_outfile, generated_outstream)]
            self.componentTest(Tools.BGCORRECT, fdstools_args, test_data)


    def test_advanced1_param(self):
        """
        Single infile and outfile
        """
        infile = self.data_dir_root / "_references" / "F0107.txt"
        outfile =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_default.txt")

        fdstools_args = [
        str(self.data_dir_root / Tools.BGMERGE / "bgmerge_default.txt"),
        "--input", str(infile),
        "--output",str(outfile.generated_output),
        '--library', self.library
        ]

        test_data = [outfile]
        self.componentTest(Tools.BGCORRECT, fdstools_args, test_data)


    def test_advanced2_param(self):
        """
        Single infiles in batch with regex for input files
        """
        outfile_F0107 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv2_F0107.txt")
        outfile_F0108 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv2_F0108.txt")
        outfile_F0109 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv2_F0109.txt")

        fdstools_args = [
        str(self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"),
        "--input", self.infiles,
        "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
        str(outfile_F0109.generated_output),
        '--library', self.library,
        "--tag-expr", "(F010.)"
        ]


        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.BGCORRECT, fdstools_args, test_data)


    def test_advanced3_param(self):
        """
        -> --sequence-format = raw
        -> --combine-strands turned on
        """
        outfile_F0107 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv3_F0107.txt")
        outfile_F0108 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv3_F0108.txt")
        outfile_F0109 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv3_F0109.txt")

        fdstools_args = [
            str(self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"),
            "--input", self.infiles,
            "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
            '--library', self.library,
            "--tag-expr", "(F010.)",
            '--sequence-format', 'raw',
            '--combine-strands'
        ]

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.BGCORRECT, fdstools_args, test_data)


    def test_advanced4_param(self):
        """
        --marker
        """
        outfile_F0107 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv4_F0107.txt")
        outfile_F0108 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv4_F0108.txt")
        outfile_F0109 =  OutputFile(self.data_dir_root / Tools.BGCORRECT / "bgcorrect_adv4_F0109.txt")

        fdstools_args = [
            str(self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"),
            "--input", self.infiles,
            "--output", str(outfile_F0107.generated_output), str(outfile_F0108.generated_output),
            str(outfile_F0109.generated_output),
            '--library', self.library,
            "--tag-expr", "(F010.)",
            '--marker', 'D12S391']

        test_data = [outfile_F0107, outfile_F0108, outfile_F0109]
        self.componentTest(Tools.BGCORRECT, fdstools_args, test_data)