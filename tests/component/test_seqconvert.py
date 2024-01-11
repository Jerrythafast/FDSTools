from io import StringIO
from unittest.mock import patch

from tests.lib.FDSToolsComponentTestCase import FDSToolsComponentTestCase


class Test(FDSToolsComponentTestCase):
    # todo: test update allele name with --library2
    # todo: test reverse complement with --reverse-complement and --library2
    # todo: test error --reverse-complement without --library2

    def setUp(self):
        super().setUp()
        self.library = "ForenSeqA"
        self.raw = self.data_dir / "tssv" / "ForenseqAMixture_R1_tssv_advanced1.txt"
        self.tssv = self.data_dir / "seqconvert" / "ForenseqAMixture_R1_seqconvert_tssv.txt"
        self.allelename = self.data_dir / "seqconvert" / \
                                          "ForenseqAMixture_R1_seqconvert_allelename.txt"
    # setUp

    def test_raw_to_raw(self):
        """Also tests stdin and stdout."""
        fdstools_args = ["raw"]

        with open(self.raw) as infile,\
                patch("sys.stdin", infile), \
                patch("sys.stdout", StringIO()) as outfile_generated:
            self.assertToolWorks("seqconvert", fdstools_args, self.raw, outfile_generated)
    # test_raw_to_raw

    def test_raw_to_tssv(self):
        msg = "ValueError: missing library"
        fdstools_args = ["tssv", str(self.raw), self.tssv.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.tssv)
    # test_raw_to_tssv

    def test_raw_to_allelename(self):
        msg = "ValueError: missing library"
        fdstools_args = ["allelename", str(self.raw), self.allelename.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.allelename)
    # test_raw_to_allelename

    def test_tssv_to_raw(self):
        msg = "ValueError: missing library"
        fdstools_args = ["raw", str(self.tssv), self.raw.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.raw)
    # test_tssv_to_raw

    def test_tssv_to_tssv(self):
        fdstools_args = ["tssv", str(self.tssv), self.tssv.name]
        self.assertToolWorks("seqconvert", fdstools_args, outfile_expected=self.tssv)
    # test_tssv_to_tssv

    def test_tssv_to_allelename(self):
        msg = "ValueError: missing library"
        fdstools_args = ["allelename", str(self.tssv), self.allelename.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.allelename)
    # test_tssv_to_allelename

    def test_allelename_to_raw(self):
        msg = "ValueError: missing library"
        fdstools_args = ["raw", str(self.allelename), self.raw.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.raw)
    # test_allelename_to_raw

    def test_allelename_to_tssv(self):
        msg = "ValueError: missing library"
        fdstools_args = ["tssv", str(self.allelename), self.tssv.name]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        fdstools_args += ["--library", self.library]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=self.tssv)
    # test_allelename_to_tssv

    def test_allelename_to_allelename(self):
        fdstools_args = ["allelename", str(self.allelename), self.allelename.name]
        self.assertToolWorks("seqconvert", fdstools_args, self.allelename)
    # test_allelename_to_allelename

    def test_allele_column(self):
        """Tests -a, --allele-column and -c, --output-column"""
        sc_dir = self.data_dir / "seqconvert"
        a = self.allelename
        a_t = sc_dir / "ForenseqAMixture_R1_seqconvert-allele-column_allelename+tssv.txt"
        a_a = sc_dir / "ForenseqAMixture_R1_seqconvert-allele-column_allelename+allelename.txt"
        r_a = sc_dir / "ForenseqAMixture_R1_seqconvert-allele-column_raw+allelename.txt"

        fdstools_args_a = ["tssv", str(a), a_t.name, "--output-column", "custom"]
        fdstools_args_c = ["allelename", str(a_t), a_a.name, "--allele-column", "custom"]
        fdstools_args_a_c = ["raw", str(a_a), r_a.name, "--allele-column", "custom", "--output-column", "sequence"]

        msg = "ValueError: missing library"
        self.subTestToolRaisesException("seqconvert", fdstools_args_a, SystemExit, 2, msg=msg)
        self.subTestToolRaisesException("seqconvert", fdstools_args_c, SystemExit, 2, msg=msg)
        self.subTestToolRaisesException("seqconvert", fdstools_args_a_c, SystemExit, 2, msg=msg)

        fdstools_args_a += ["--library", self.library]
        fdstools_args_c += ["--library", self.library]
        fdstools_args_a_c += ["--library", self.library]

        self.subTestToolWorks("seqconvert", fdstools_args_a, outfile_expected=a_t)
        self.subTestToolWorks("seqconvert", fdstools_args_c, outfile_expected=a_a)
        self.subTestToolWorks("seqconvert", fdstools_args_a_c, outfile_expected=r_a)
    # test_allele_column

    def test_marker_column(self):
        """Tests -m, --marker-column"""
        file = self.data_dir / "seqconvert" / "ForenseqAMixture_R1_seqconvert-marker-column.txt"

        # Missing library error is thrown before missing marker column error.
        msg = "ValueError: missing library"
        fdstools_args = ["raw", str(file)]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        # No Error thrown when no conversion is required and marker column is not found.
        fdstools_args = ["allelename", str(file), file.name]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=file)

        # Error thrown when including --library command, even when no conversion is required.
        msg = "ValueError: Column not found in input file: marker"
        fdstools_args += ["--library", self.library]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        # Error disappears when marker column is clarified.
        fdstools_args += ["--marker-column", "custom"]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=file)
    # test_marker_column

    def test_marker(self):
        """Tests -M, --marker"""
        sc_dir = self.data_dir / "seqconvert"
        t_str = sc_dir / "ForenseqAMixture_R1_seqconvert-marker_tssv_STR.txt"
        t_nonstr = sc_dir / "ForenseqAMixture_R1_seqconvert-marker_tssv_NonSTR.txt"
        a_str = sc_dir / "ForenseqAMixture_R1_seqconvert-marker_allelename_STR.txt"
        a_nonstr = sc_dir / "ForenseqAMixture_R1_seqconvert-marker_allelename_NonSTR.txt"

        fdstools_args = ["tssv", str(self.raw), t_str.name, "--library", self.library, "-M", "D13S317"]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=t_str)

        fdstools_args = ["tssv", str(self.raw), t_nonstr.name, "--library", self.library,
                         "--marker", "rs10776839"]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=t_nonstr)

        fdstools_args = ["allelename", str(self.tssv), a_str.name, "--library", self.library,
                         "--marker", "D13S317"]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=a_str)

        fdstools_args = ["allelename", str(self.tssv), a_nonstr.name, "--library", self.library,
                         "--marker", "rs10776839"]
        self.subTestToolWorks("seqconvert", fdstools_args, outfile_expected=a_nonstr)

        # Test below would also break on any non-STR allelename
        msg = "error: Invalid allele name 'REF' encountered!"
        fdstools_args = ["raw", str(self.allelename), "--library", self.library, "--marker", "D13S317"]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)

        # Test below would also break on out of range non-STR (as e.g. SNP naming contains location)
        msg = "error: Unrecognised variant 'CE10_TCTA[10]ATCT[3]'"
        fdstools_args = ["raw", str(self.allelename), "--library", self.library, "--marker", "rs10776839"]
        self.subTestToolRaisesException("seqconvert", fdstools_args, SystemExit, 2, msg=msg)
    # test_marker
