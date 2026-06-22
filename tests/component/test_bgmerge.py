from io import StringIO
from unittest.mock import patch
from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile


def expected_mismatches():
    regex_1 = r"FGA\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t13.15[1]*\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_2 = r"D6S1043\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t6.47[32]\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_3 = r"D18S51\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t100.0\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_4 = r"D6S1043\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t0.55\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_5 = r"D13S317\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t3.675\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_6 = r"D13S317\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t100.0\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_7 = r"D18S51\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t12.963\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_8 = r"FGA\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t100.0\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_9 = r"D6S1043\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t100.0\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_10 = r"D18S51\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t0.509\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_11 = r"FGA\t[ACTG]+\t[ACTG]+\t0.0\t0.0\t0.561\t(?:bgestimate|bgpredict),(?:bgestimate|bgpredict)\n"
    regex_tuple = (regex_1, regex_2, regex_3, regex_4, regex_5, regex_6, regex_7, regex_8,
                   regex_9, regex_10, regex_11)
    return regex_tuple

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
            test_data[0].accepted_mismatch_tuple = expected_mismatches()
            self.componentTest(Tools.BGMERGE, fdstools_args, test_data)

    def test_advanced_param(self):
        """

        """
        infile = self.data_dir_root / Tools.BGESTIMATE / "bgestimate_adv2.txt"
        infile_2 = self.data_dir_root / Tools.BGPREDICT / "bgpredict-adv2.txt"
        outfile = OutputFile(self.data_dir_root / Tools.BGMERGE / "bgmerge_default.txt")
        outfile.accepted_mismatch_tuple = expected_mismatches()

        fdstools_args = [str(infile),
                         str(infile_2),
                         "--output", str(outfile.generated_output)
                         ]
        test_data = [outfile]
        self.componentTest(Tools.BGMERGE, fdstools_args, test_data)

