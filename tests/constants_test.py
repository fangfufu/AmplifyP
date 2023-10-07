# -*- coding: utf-8 -*-
"""Tests for constants.py."""
from dataclasses import FrozenInstanceError

from amplifyp.constants import Nucleotide


def test_nucleotide() -> None:
    """Test the Nucleotide constant."""
    try:
        Nucleotide.single = "XKCD"
    except FrozenInstanceError:
        pass
