from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

__all__ = ["experiments", "plotting", "networkanalysis", "stats", "jupyter"]
