#
# Copyright (C) 2023 Jerry Hoogenboom
#
# This file is part of FDSTools, data analysis tools for Massively
# Parallel Sequencing of forensic DNA markers.
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

"""
Data analysis tools for Massively Parallel Sequencing of forensic DNA markers,
including tools for characterisation and filtering of PCR stutter artefacts and
other systemic noise, and for automatic detection of the alleles in a sample.
"""

__version__ = "2.0.4"
usage = __doc__.split("\n\n\n")


def version(name, toolname=None, toolversion=None):
    """Return a version string for the package or a given tool."""
    verformat = "%s %s"
    toolverformat = "%s (part of %s)"
    if toolname is None:
        return verformat % (name, __version__)
    if toolversion is None:
        return toolverformat % (toolname, verformat % (name, __version__))
    return toolverformat % (verformat % (toolname, toolversion),
                            verformat % (name, __version__))
