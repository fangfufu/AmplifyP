# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Settings sub-controls package."""

from amplifyp.gui.views.settings.score_table import ScoreTable
from amplifyp.gui.views.settings.replication_tile import ReplicationTile
from amplifyp.gui.views.settings.tm_tile import TmTile
from amplifyp.gui.views.settings.primer_dimer_tile import PrimerDimerTile
from amplifyp.gui.views.settings.appearance_tile import AppearanceTile

__all__ = [
    "ScoreTable",
    "ReplicationTile",
    "TmTile",
    "PrimerDimerTile",
    "AppearanceTile",
]
