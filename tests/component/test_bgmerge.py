from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


class TestBGMerge(ComponentTestCase):

    def test_default_param(self):
        """
        -> Standard streams
        -> Output is used as input for component tests of bgcorrect.
        """
        infile = self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"
        infile_2 = self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv2.txt"
        expected_outfile = self.data_dir_root / Tools.BGMERGE / "bgmerge_default.txt"

        fdstools_args = [str(infile)]

        with open(infile_2) as infile2,\
            patch("sys.stdin", infile2), \
            patch("sys.stdout", StringIO()) as generated_outstream:
            test_data = [OutputFile(expected_outfile, generated_outstream)]
            self.componentTest(Tools.BGMERGE, fdstools_args, test_data)

    def test_advanced_param(self):
        """

        """
        infile = self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"
        infile_2 = self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv2.txt"
        outfile = OutputFile(self.data_dir_root / Tools.BGMERGE / "bgmerge_default.txt")

        fdstools_args = [str(infile),
                         str(infile_2),
                         "--output", str(outfile.generated_output)
                         ]
        test_data = [outfile]
        self.componentTest(Tools.BGMERGE, fdstools_args, test_data)

