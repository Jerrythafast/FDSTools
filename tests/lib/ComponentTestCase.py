#!/usr/bin/env python3
import difflib
import os
import re
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from typing import Optional, Type, Tuple, List
from unittest import TestCase

from fdstools import fdstools
from pathlib import Path
from unittest.mock import patch

from tests.lib.ComponentTestData import OutputFile


class Tools:
    ALLELEFINDER = 'allelefinder'
    BGANALYSE = 'bganalyse'
    BGCORRECT = 'bgcorrect'
    BGESTIMATE = 'bgestimate'
    BGHOMRAW = 'bghomraw'
    BGHOMSTATS = 'bghomstats'
    BGMERGE = 'bgmerge'
    BGPREDICT = 'bgpredict'
    FINDNEWALLELES = 'findnewalleles'
    LIBCONVERT = 'libconvert'
    LIBRARY = 'library'
    MPS2CE = 'mps2ce'
    PIPELINE = 'pipeline'
    SAMPLESTATS = 'samplestats'
    SEQCONVERT = 'seqconvert'
    STUTTERMARK = 'stuttermark'
    STUTTERMODEL ='stuttermodel'
    TSSV ='tssv'
    VIS = 'vis'

class ComponentTestCase(TestCase):
    data_dir_root = Path(__file__).parent.parent / "data"

    def setUp(self) -> None:
        self.original_working_directory = os.getcwd()
        self.temp_dir = TemporaryDirectory()
        os.chdir(self.temp_dir.name)

    def tearDown(self) -> None:
        os.chdir(self.original_working_directory)
        self.temp_dir.cleanup()

    def componentTest(
            self,
            tool: str,
            fdstools_args: List[str],
            test_data: List[OutputFile]
    ):
        """
        :param tool: The tool name
        :param fdstools_args: a list containing all arguments for the tool, excluding the tool
        name itself.
        :param test_data: List[OutputFile] contains expected and generated output files.
        """

        with self.subTestExceptionWithArgs(error_type=None):
            self._run_tool(tool, fdstools_args)

        for test_output in test_data:
            file_name = test_output.expected_output.name

            with self.subTest(file_name):
                test_output.update_expected_output_file()
                self._assertEqualOutput(test_output)


    def componentTestSystemFileCount(
            self,
            dir_path: Path,
            expected_dir_count: Optional[int],
            expected_file_count: Optional[int]
    ):
        dir_count, file_count = self._count_dir_and_file_objects_in_tree(dir_path)

        if expected_dir_count:
            with self.subTest("Dir count"):
                self.assertEqual(dir_count, expected_dir_count)
        if expected_file_count:
            with self.subTest("File count"):
                self.assertEqual(file_count, expected_file_count)

    def componentTestException(
            self,
            tool: str,
            fdstools_args: list,
            error_type: Type[BaseException] = SystemExit,
            *error_args,
            test_msg: str = ""
    ):
        with self.subTestExceptionWithArgs(error_type, *error_args, test_msg=test_msg):
            self._run_tool(tool, fdstools_args)

    def _assertEqualOutput(
            self,
            output: OutputFile
    ):
        """
        Compares the generated output to the expected output.
        Performs a subTest for each line that is different, missing or added to the output.
        Ends with an assertion that determines if any unexpected differences were observed in the output.
        Uses the 'accepted_mismatch_tuple' from OutputFile.
        """

        expected_line_list = output.load_expected_output()
        generated_line_list = output.load_generated_output()

        incorrect_or_missing_lines = set(expected_line_list).difference(generated_line_list)
        unexpected_lines = set(generated_line_list).difference(expected_line_list)

        correct_output = True

        for line in incorrect_or_missing_lines:
            msg, acceptable = self._assess_mismatch_type(output.accepted_mismatch_tuple, line)

            similar_line = difflib.get_close_matches(line, unexpected_lines, n=1)
            if similar_line:
                # Expected line differs from generated i
                msg += f" expected: {line}"
                msg += f" generated: {similar_line[0]}"

                # so we don't check it again below
                unexpected_lines.remove(similar_line[0])
            else:
                msg += f" Expected line does not exist in generated output: {line}"

            if not acceptable:
                correct_output = False

            with self.subTest(msg=msg):
                self.assertTrue(acceptable, msg)

        if not output.allow_extra_lines:
            for line in unexpected_lines:
                msg, acceptable = self._assess_mismatch_type(output.accepted_mismatch_tuple, line)
                msg += f" Generated line was not expected in output: {line}"

                if not acceptable:
                    correct_output = False

                with self.subTest(msg=msg):
                    self.assertTrue(acceptable, msg)
        self.assertTrue(correct_output)

    @contextmanager
    def subTestExceptionWithArgs(
            self,
            error_type=None,
            *error_args,
            test_msg=""
    ):
        """
        :param error_type: When 'None', you expect no exception.
        :param error_args: expected contents of BaseException.args
        :param test_msg: Message for the subTest description.
        """
        try:
            yield
        except BaseException as e:
            if error_type is None:
                with self.subTest("Incorrect exception"):
                    self.fail(f"Unexpected exception was raised: {str(type(e))}")
            else:
                with self.subTest(f"{test_msg} [Exception type]"):
                    self.assertIsInstance(e, error_type)
                with self.subTest(f"{test_msg} [Exception args]"):
                    self.assertTupleEqual(e.args, error_args)
        else:
            if error_type is not None:
                with self.subTest("Missing exception"):
                    self.fail(f"No exception raised, but expected: {str(error_type)}")

    @classmethod
    def _count_dir_and_file_objects_in_tree(cls, dir_path: Path or str):
        dir_count = 0
        file_count = 0
        for entry in Path.iterdir(dir_path):
            if Path.is_dir(entry):
                dir_count += 1
                child_dir_count, child_file_count = cls._count_dir_and_file_objects_in_tree(entry)
                dir_count += child_dir_count
                file_count += child_file_count
            elif Path.is_file(entry):
                file_count += 1
        return dir_count, file_count

    @staticmethod
    def _run_tool(tool: str, fdstools_args: List[str]):
        fdstools_args = ["fdstools", tool] + fdstools_args
        with patch("sys.argv", fdstools_args):
            fdstools.main()

    @staticmethod
    def _assess_mismatch_type(accepted_mismatch_tuple: Tuple[str, ...], line: str) \
            -> Tuple[str, bool]:
        """use the accepted_mismatch_tuple to check if this mismatch is acceptable or
        not."""
        if any(re.search(mismatch, line) for mismatch in accepted_mismatch_tuple):
            msg = "Acceptable mismatch:\n"
            acceptable = True
        else:
            msg = "Unacceptable mismatch:\n"
            acceptable = False
        return msg, acceptable