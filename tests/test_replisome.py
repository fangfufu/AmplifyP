# -*- coding: utf-8 -*-
from amplifyp.replisome import DNA


def test_DNA_invalid_symbol() -> None:
    DNA("ACGTX")
