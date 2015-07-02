"""
Tools for characterisation and filtering of PCR stutter artefacts and other
systemic noise in Next Generation Sequencing data of forensic STR markers.
"""

__version_info__ = ('0', '0', '1')
__version__ = '.'.join(__version_info__)
usage = __doc__.split("\n\n\n")

def version(name, toolname=None, toolversion=None):
    """Return a version string for the package or a given tool."""
    verformat = "%s %s"
    toolverformat = "%s (part of %s)"
    if toolname is None:
        return verformat % (name, __version__)
    if toolversion is None:
        return toolverformat % (toolname, verformat % (name, __version__))
    return toolverformat % (verformat % (toolname,toolversion),
                            verformat % (name, __version__))
