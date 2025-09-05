from io import StringIO
from unittest.mock import patch

from tests.lib.ComponentTestCase import ComponentTestCase, Tools
from tests.lib.ComponentTestData import OutputFile

class SequenceFormat:
    RAW = "raw"
    TSSV = "tssv"
    ALLELENAME = "allelename"

class TestSeqconvert(ComponentTestCase):
    """Use test_tssv.TestTssv().test_default_param() output so there are SNPs present."""
    # todo: test update allele name with --library2
    # todo: test reverse complement with --reverse-complement and --library2
    # todo: test error --reverse-complement without --library2
    library = "ID-OmniSTR"

    def setUp(self):
        super().setUp()
        self.raw = self.data_dir_root / Tools.TSSV / 'OmniSTR_Mixture_R1_4lt2_tssv_default.txt'
        self.tssv = (self.data_dir_root / Tools.SEQCONVERT /
                     'OmniSTR_Mixture_R1_4lt2_seqconvert_tssv.txt')
        self.allelename = (self.data_dir_root / Tools.SEQCONVERT /
                           'OmniSTR_Mixture_R1_4lt2_seqconvert_allelename.txt')

    def test_raw_to_raw(self):
        """Also tests stdin and stdout."""
        fdstools_args = [SequenceFormat.RAW]

        with open(self.raw) as infile,\
                patch("sys.stdin", infile), \
                patch("sys.stdout", StringIO()) as outfile_generated:
            test_data = [OutputFile(expected_output=self.raw, generated_output=outfile_generated)]
            self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_raw_to_tssv(self):
        fdstools_args = [SequenceFormat.TSSV, str(self.raw), self.tssv.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.tssv)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_raw_to_allelename(self):
        fdstools_args = [SequenceFormat.ALLELENAME, str(self.raw), self.allelename.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.allelename)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_tssv_to_raw(self):
        fdstools_args = [SequenceFormat.RAW, str(self.tssv), self.raw.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.raw)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_tssv_to_tssv(self):
        fdstools_args = [SequenceFormat.TSSV, str(self.tssv), self.tssv.name]
        test_data = [OutputFile(expected_output=self.tssv)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_tssv_to_allelename(self):
        fdstools_args = [SequenceFormat.ALLELENAME, str(self.tssv), self.allelename.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.allelename)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_allelename_to_raw(self):
        fdstools_args = [SequenceFormat.RAW, str(self.allelename), self.raw.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.raw)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_allelename_to_tssv(self):
        fdstools_args = [SequenceFormat.TSSV, str(self.allelename), self.tssv.name]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args += ["--library", self.library]
        test_data = [OutputFile(expected_output=self.tssv)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_allelename_to_allelename(self):
        fdstools_args = [SequenceFormat.ALLELENAME, str(self.allelename), self.allelename.name]
        test_data = [OutputFile(expected_output=self.allelename)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_allele_column(self):
        """
        Tests -a, --allele-column and -c, --output-column
        Use ref file as input to reduce run time.
        """
        library = "ForenSeqA"
        sc_dir = self.data_dir_root / Tools.SEQCONVERT
        t = self.data_dir_root / "_references" / "F0107.txt"
        t_a = sc_dir / "F0107_seqconvert-allele-column_tssv+allelename.txt"
        t_t = sc_dir / "F0107_seqconvert-allele-column_tssv+tssv.txt"
        r_t = sc_dir / "F0107_seqconvert-allele-column_raw+tssv.txt"

        fdstools_args_a = [SequenceFormat.ALLELENAME, str(t), t_a.name, "--output-column", "custom"]
        fdstools_args_c = [SequenceFormat.TSSV, str(t_a), t_t.name, "--allele-column", "custom"]
        fdstools_args_a_c = [SequenceFormat.RAW, str(t_t), r_t.name, "--allele-column", "custom",
                             "--output-column", "sequence"]

        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args_a,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args_c,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args_a_c,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        fdstools_args_a += ["--library", library]
        fdstools_args_c += ["--library", library]
        fdstools_args_a_c += ["--library", library]

        self.componentTest(Tools.SEQCONVERT, fdstools_args_a, [OutputFile(expected_output=t_a)])
        self.componentTest(Tools.SEQCONVERT, fdstools_args_c, [OutputFile(expected_output=t_t)])
        self.componentTest(Tools.SEQCONVERT, fdstools_args_a_c, [OutputFile(expected_output=r_t)])

    def test_marker_column(self):
        """Tests -m, --marker-column"""
        file = (self.data_dir_root / Tools.SEQCONVERT /
                "OmniSTR_Mixture_R1_4lt2_seqconvert-marker-column.txt")

        # Missing library error is thrown before missing marker column error.
        fdstools_args = [SequenceFormat.RAW, str(file)]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: missing library"
        )

        # No Error thrown when no conversion is required and marker column is not found.
        fdstools_args = [SequenceFormat.ALLELENAME, str(file), file.name]
        test_data = [OutputFile(expected_output=file)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

        # Error thrown when including --library command, even when no conversion is required.
        fdstools_args += ["--library", self.library]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="ValueError: Column not found in input file: marker"
        )

        # Error disappears when marker column is clarified.
        fdstools_args += ["--marker-column", "custom"]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

    def test_marker(self):
        """Tests -M, --marker"""
        sc_dir = self.data_dir_root / Tools.SEQCONVERT
        t_str = sc_dir / "OmniSTR_Mixture_R1_4lt2_seqconvert-marker_tssv_STR.txt"
        t_non_str = sc_dir / "OmniSTR_Mixture_R1_4lt2_seqconvert-marker_tssv_NonSTR.txt"
        a_str = sc_dir / "OmniSTR_Mixture_R1_4lt2_seqconvert-marker_allelename_STR.txt"
        a_non_str = sc_dir / "OmniSTR_Mixture_R1_4lt2_seqconvert-marker_allelename_NonSTR.txt"

        fdstools_args = [
            SequenceFormat.TSSV,
            str(self.raw),
            t_str.name,
            "--library", self.library,
            "-M", "D13S317"
        ]
        test_data = [OutputFile(expected_output=t_str)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

        fdstools_args = [
            SequenceFormat.TSSV,
            str(self.raw),
            t_non_str.name,
            "--library", self.library,
            "--marker", "AmelogeninY"
        ]
        test_data = [OutputFile(expected_output=t_non_str)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

        fdstools_args = [
            SequenceFormat.ALLELENAME,
            str(self.tssv), a_str.name,
            "--library", self.library,
            "--marker", "D13S317"
        ]
        test_data = [OutputFile(expected_output=a_str)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

        fdstools_args = [
            SequenceFormat.ALLELENAME,
            str(self.tssv),
            a_non_str.name,
            "--library", self.library,
            "--marker", "AmelogeninY"
        ]

        test_data = [OutputFile(expected_output=a_non_str)]
        self.componentTest(Tools.SEQCONVERT, fdstools_args, test_data)

        # Test below would also break on any non-STR allelename
        fdstools_args = [
            SequenceFormat.RAW,
            str(self.allelename),
            "--library", self.library,
            "--marker", "D13S317"
        ]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="error: Invalid allele name 'REF' encountered!"
        )

        # Test below would also break on an STR throwing an 'Unrecognised variant' error.
        fdstools_args = [
            SequenceFormat.RAW,
            str(self.allelename),
            "--library", self.library,
            "--marker", "AmelogeninY"
        ]
        self.componentTestException(
            Tools.SEQCONVERT, fdstools_args,
            SystemExit, 2,
            test_msg="error: Position 11296932 is outside sequence range"
        )
