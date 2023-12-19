#!/usr/bin/env python3
import gzip
import os
import sys
from contextlib import contextmanager

from io import StringIO
from tempfile import TemporaryDirectory

from fdstools import fdstools
from pathlib import Path
from unittest.mock import patch

from tests.lib.FDSToolsTestCase import FDSToolsTestCase


def _count_system_file_objects_in_tree(dir_path: Path or str):
    dir_count = 0
    file_count = 0
    for entry in Path.iterdir(dir_path):
        if Path.is_dir(entry):
            dir_count += 1
            child_dir_count, child_file_count = _count_system_file_objects_in_tree(entry)
            dir_count += child_dir_count
            file_count += child_file_count
        elif Path.is_file(entry):
            file_count += 1
    return dir_count, file_count


def _run_tool(tool: str, fdstools_args: list):
    fdstools_args = ["fdstools", tool] + fdstools_args
    with patch("sys.argv", fdstools_args):
        fdstools.main()
    # _run_tool


@contextmanager
def _smart_open(file: Path or StringIO, mode="rt", encoding="UTF-8"):
    """
    So asserts and subTests for files can handle StringIO
    """
    close_file = False
    if isinstance(file, (Path, str)):
        if str(file).endswith(".gz"):
            file = gzip.open(file, mode, encoding=encoding)
        else:
            file = open(file, mode, encoding=encoding)
        close_file = True
    elif isinstance(file, StringIO):
        file.seek(0)
    else:
        raise NotImplementedError(str(type(file)) + " was not implemented for smart_open.")

    yield file.readlines()

    if close_file:
        file.close()
# _smart_open


class FDSToolsComponentTestCase(FDSToolsTestCase):
    data_dir = Path(__file__).parent.parent / "data"
    original_working_directory = os.getcwd()

    def setUp(self) -> None:
        self.temp_dir = TemporaryDirectory()
        os.chdir(self.temp_dir.name)
        self.update_test_data = False
    # setUp

    def tearDown(self) -> None:
        os.chdir(self.original_working_directory)
        self.temp_dir.cleanup()
    # tearDown

    def assertFileEqual(self, outfile_expected: Path, outfile_generated: Path, allow_extra_lines_in_output=False):
        """
            Limitation: when allow_extra_lines_in_output, empty outfile_expected matches outfile_generated.
        """
        self.assertFalse(outfile_expected == outfile_generated,
                         "The expected and generated output file are not the same file.")

        outfile_expected_path = outfile_expected

        is_stderr = outfile_generated is sys.stderr
        with _smart_open(outfile_generated) as outfile_generated:
            try:
                with _smart_open(outfile_expected) as outfile_expected:
                    if allow_extra_lines_in_output or is_stderr:
                        self.assertListContains(sorted(outfile_expected), sorted(outfile_generated))
                    else:
                        self.assertListEqual(sorted(outfile_expected), sorted(outfile_generated))
            except (AssertionError, FileNotFoundError) as e:
                raise e
            finally:
                self.updateOutfileExpected(outfile_expected_path, outfile_generated)

    # assertFileEqual

    def assertToolWorks(self, tool: str, fdstools_args: list, outfile_expected: Path,
                        outfile_generated: Path or StringIO = None, allow_extra_lines_in_output=False):
        """
        fdstools_args: a list containing all arguments for the tool, excluding the tool name itself.
        outfile_generated: if None, assign outfile_expected.name
        """
        if not outfile_generated:
            outfile_generated = Path(outfile_expected.name)

        _run_tool(tool, fdstools_args)

        self.assertFileEqual(outfile_expected, outfile_generated, allow_extra_lines_in_output)
    # assertToolWorks

    def subTestExpectedSystemFileObjectCount(self, dir_path: Path, expected_dir_count: int,
                                             expected_file_count: int, msg=""):
        dir_count, file_count = _count_system_file_objects_in_tree(dir_path)

        if expected_dir_count:
            self.subTestEqual(dir_count, expected_dir_count, "Dir count "+msg)
        if expected_file_count:
            self.subTestEqual(file_count, expected_file_count, "File count "+msg)
    # _count_file_system_objects

    def subTestFileEqual(self, outfile_expected: Path, outfile_generated: Path or StringIO,
                         msg: str = None, allow_extra_lines_in_output=False):
        with self.subTest(msg=msg):
            self.assertFileEqual(outfile_expected, outfile_generated, allow_extra_lines_in_output)
    # subTestFileEqual

    def subTestToolRaisesException(self, tool: str, fdstools_args: list, error_type, *args, msg: str = None):
        with self.subTestExceptionWithArgs(error_type, *args, msg=msg):
            _run_tool(tool, fdstools_args)
    # subTestToolRaisesException

    def subTestToolWorks(self, tool: str, fdstools_args: list, outfile_expected: Path,
                         outfile_generated: Path or StringIO = None, allow_extra_lines_in_output=False):
        with self.subTest(msg=outfile_expected.name):
            self.assertToolWorks(tool, fdstools_args, outfile_expected, outfile_generated, allow_extra_lines_in_output)
    # subTestToolWorks

    def subTestToolWorksMultiOutput(self, tool: str, fdstools_args: list, outdir_expected: Path,
                                    outfile_names: list or dict, allow_extra_lines_in_output=False):
        """
        fdstools_args: a list containing all arguments for the tool, excluding the tool name itself.
        outfile_names: dict[outfile_expected: str, outfile_generated: str or StringIO]
                            or
                       list[str] -> outfile_expected and outfile_generated outfile have the same name.
        """
        _run_tool(tool, fdstools_args)

        if not isinstance(outfile_names, dict):
            outfile_names = dict(zip(outfile_names, outfile_names))

        for outfile_expected, outfile_generated in outfile_names.items():
            msg = outfile_expected
            outfile_expected = outdir_expected / outfile_expected

            self.subTestFileEqual(outfile_expected, outfile_generated, msg, allow_extra_lines_in_output)
    # subTestToolWorksMultiOutput

    def updateOutfileExpected(self, outfile_expected_path: Path or str, outfile_generated: list):
        """
            limitation: It will write OS specific lines from stderr to the expected output file.
                        But this can be caught quite well via git during commit (or review).
        """
        if self.update_test_data:
            Path(outfile_expected_path).parent.mkdir(exist_ok=True, parents=True)
            with open(outfile_expected_path, mode='w') as outfile_expected:
                for line in outfile_generated:
                    outfile_expected.write(line)
    # update_outfile_expected
# TestUtils
