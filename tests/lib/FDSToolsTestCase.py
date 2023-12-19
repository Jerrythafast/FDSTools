#!/usr/bin/env python3
import difflib
import os
from contextlib import contextmanager
from typing import Any
from unittest import TestCase


def run_long_tests():
    """Intended usage: @skipUnless(run_long_tests(), "Test takes longer than 30 seconds.")"""
    try:
        return os.environ['RUN_LONG_TESTS'].lower() in ("true", "t", "1", "yes", "y", "ok", "okay", "yes, please")
    except KeyError:
        return False
# run_long_tests


class FDSToolsTestCase(TestCase):
    def assertListContains(self, output_expected: list, output_generated: list):
        """
            Limitation: An empty output_expected list matches with output_generated
        """
        for line in output_expected:
            if line not in output_generated:
                match = difflib.get_close_matches(line, output_generated, n=1)
                if match:
                    self.assertEqual(line, match[0])
                else:
                    self.fail(msg="Expected line does not exist in output: " + line)
    # assertListContains

    def subTestDictEqual(self, output_expected: dict, output_generated: dict, msg: str = None):
        with self.subTest(msg=msg):
            self.assertDictEqual(output_expected, output_generated)
    # subTestDictEqual

    def subTestEqual(self, output_expected: Any, output_generated: Any, msg: str = None):
        with self.subTest(msg=msg):
            self.assertEqual(output_expected, output_generated)
    # subTestEqual

    @contextmanager
    def subTestExceptionWithArgs(self, error_type=None, *args, msg: str = None):
        """
            When error_type=None, you expect no exception.
        """
        try:
            yield
            if error_type is None:
                self.subTestEqual("", "", "As expected, no exception was raised.")
            else:
                self.subTestFail("Expected " + str(error_type) + " was not raised: " + msg)
        except BaseException as e:
            if error_type is None:
                self.subTestFail("Unexpected exception was raised: " + str(type(e)))
            else:
                with self.subTest(msg=msg):
                    self.assertIsInstance(e, error_type)
                    self.assertTupleEqual(e.args, args)
    # subTestExceptionWithArgs

    def subTestFail(self, msg: str = None):
        with self.subTest(msg=msg):
            self.fail(msg)
    # subTestFail

    def subTestListEqual(self, output_expected: list, output_generated: list, msg: str = None):
        with self.subTest(msg=msg):
            self.assertListEqual(output_expected, output_generated)
    # subTestListEqual
# FDSToolsUnitTestCase
