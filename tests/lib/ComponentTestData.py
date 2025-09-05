import gzip
import os
from contextlib import contextmanager
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import List, AnyStr, Optional, Tuple, Union, Generator

data_dir_root = Path(__file__).parent.parent / "data"

@dataclass
class OutputFile:
    """
    :param expected_output:     included in current code repo: root/tests/data (=data_dir_root)
    :param generated_output:    generated during the test in the current working directory (=tempdir)
    :param update_test_data:    update the expected_output with the generated_output
                                (can also be set through environment variables)
    :param accepted_mismatch_tuple: tuple of regexes matching with lines that are
                                    allowed to mismatch between expected and generated output.
                                    (e.g. a line that sometimes mismatches because a number in
                                    it is inconsistent. Or lines that contain local paths.)
    :param allow_extra_lines:   when True, generated output may contain lines that are not
                                present in the expected output. So you can better control the
                                repository size or remove sensitive/system specific information
                                from the test data.
    """
    expected_output: Path
    generated_output: Optional[Union[Path, StringIO]] = None  # TODO python 3.10: `generated_output: Path | StringIO`

    # component test settings:
    update_test_data: bool = (os.environ.get("UPDATE_TEST_DATA", "false").lower()
                              in
                              ("true", "t", "1", "yes", "y", "ok", "okay", "yes, please"))
    accepted_mismatch_tuple: Optional[Tuple[str, ...]] = field(default_factory=tuple)
    allow_extra_lines: bool = False

    def __post_init__(self):
        if not self.generated_output:
            self.generated_output = Path(self.expected_output.name)

        if self.expected_output == self.generated_output:
            raise ValueError(f"Expected and Generated output should not be the same file: "
                             f"{self.expected_output}")

    def load_expected_output(self) -> List[AnyStr]:
        with self._smart_open(file=self.expected_output) as file:
            return file

    def load_generated_output(self) -> List[AnyStr]:
        with self._smart_open(file=self.generated_output) as file:
            return file

    def update_expected_output_file(self):
        """
            limitation: It will write OS specific lines from stderr to the expected output file.
                        But this can be caught quite well via git during commit (or review).
        """
        if self.update_test_data:
            self.expected_output.parent.mkdir(exist_ok=True, parents=True)
            with self._smart_open(self.expected_output, mode='w') as outfile_expected:
                for line in self.load_generated_output():
                    outfile_expected.write(line)

    @staticmethod
    @contextmanager
    def _smart_open(file: Path or StringIO, mode: str = "rt", encoding: str = "UTF-8") \
            -> Generator[list[bytes], None, None]:
        """So asserts and subTests for files can handle StringIO"""
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
            raise NotImplementedError(f"{str(type(file))} was not implemented for smart_open.")

        yield file.readlines()

        if close_file:
            file.close()