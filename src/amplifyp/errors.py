"""Custom exceptions for AmplifyP."""


class DuplicateRepliconfError(ValueError):
    """Exception raised when attempting to add a duplicate Repliconf.

    This error occurs when a `Repliconf` object is added to an `AmpliconGenerator`
    that already contains an identical configuration (same primer and template).
    """
