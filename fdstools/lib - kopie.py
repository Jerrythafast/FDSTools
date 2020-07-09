#!/usr/bin/env python

#
# Copyright (C) 2020 Jerry Hoogenboom
#
# This file is part of FDSTools, data analysis tools for Next
# Generation Sequencing of forensic DNA markers.
#
# FDSTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FDSTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FDSTools.  If not, see <http://www.gnu.org/licenses/>.
#

import re

from functools import reduce


# Pattern to split a comma-, semicolon-, or space-separated list.
PAT_SPLIT_QUOTED = re.compile(r""""((?:\\"|[^"])*)"|'((?:\\'|[^'])*)'|(\S+)""")


def split_quoted_string(text):
    return reduce(
        lambda x, y: x + ["".join([
            y[0].replace("\\\"", "\""),
            y[1].replace("\\'", "'"),
            y[2]])],
        PAT_SPLIT_QUOTED.findall(text), [])
#split_quoted_string
