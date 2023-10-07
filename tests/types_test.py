# -*- coding: utf-8 -*-
"""Tests for types.py."""
from amplifyp.types import RunLengthWeightTbl


def test_run_length_weight_tbl() -> None:
    """Test the RunLengthWeightTbl class."""
    weight_tbl = RunLengthWeightTbl(5, 0.5)
    assert len(weight_tbl) == 5
    assert weight_tbl[0] == 0.5
    assert weight_tbl[1] == 0.5
    assert weight_tbl[2] == 0.5
    assert weight_tbl[3] == 0.5
    assert weight_tbl[4] == 0.5
    assert str(weight_tbl) == "[0.5, 0.5, 0.5, 0.5, 0.5]"
    assert repr(weight_tbl) == "[0.5, 0.5, 0.5, 0.5, 0.5]"
    assert list(weight_tbl) == [0.5, 0.5, 0.5, 0.5, 0.5]
    weight_tbl[0] = 0.6
    assert weight_tbl[0] == 0.6
    weight_tbl[1] = 0.7
    assert weight_tbl[1] == 0.7
    weight_tbl[2] = 0.8
    assert weight_tbl[2] == 0.8
    weight_tbl[3] = 0.9
    assert weight_tbl[3] == 0.9
    weight_tbl[4] = 1.0
    assert weight_tbl[4] == 1.0
    assert str(weight_tbl) == "[0.6, 0.7, 0.8, 0.9, 1.0]"
    assert repr(weight_tbl) == "[0.6, 0.7, 0.8, 0.9, 1.0]"
    assert list(weight_tbl) == [0.6, 0.7, 0.8, 0.9, 1.0]
