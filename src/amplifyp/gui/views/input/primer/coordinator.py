# Copyright (C) 2026 AmplifyP Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""State coordinator for DNA primers input."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.utils.data_helpers import clean_sequence

from .row import PrimerRow
from .validation import (
    get_duplicate_primer_indices,
    reconcile_primer_states,
    validate_primers,
)

if TYPE_CHECKING:
    from .input import PrimerInput


class PrimerCoordinator:
    """Handles state synchronization and data extraction.

    Coordinates validation between UI and central state.
    """

    def __init__(self, owner: PrimerInput) -> None:
        """Initialise the PrimerCoordinator.

        Args:
            owner: The parent PrimerInput component that owns this coordinator.
        """
        self.owner = owner

    def extract_primer_data_from_ui(self) -> list[dict[str, Any]]:
        """Extract primer data from UI controls.

        Returns:
            A list of dicts with keys 'name', 'seq', 'active', 'container',
            and 'checkbox' for each primer row in the UI.
        """
        ui_primers = []
        for row in self.owner.primers_list.controls:
            if not isinstance(row, PrimerRow):
                continue

            name_val = str(row.name_field.value or "").strip()
            seq_val = clean_sequence(str(row.seq_field.value or ""))
            is_active = bool(row.checkbox.value)

            ui_primers.append(
                {
                    "name": name_val,
                    "seq": seq_val,
                    "active": is_active,
                    "container": row,
                    "checkbox": row.checkbox,
                }
            )
        return ui_primers

    def get_duplicate_indices(self) -> set[int]:
        """Find indices of primers with duplicate names or sequences.

        Returns:
            Set of indices corresponding to primers with duplicate names
            or sequences in the central state.
        """
        ignore_inactive_name_dup = self.owner.settings.get(
            "ignore_inactive_name_dup_warn", True
        )
        ignore_inactive_seq_dup = self.owner.settings.get(
            "ignore_inactive_seq_dup_warn", True
        )

        return get_duplicate_primer_indices(
            self.owner.input_data.primers,
            ignore_inactive_name_dup,
            ignore_inactive_seq_dup,
        )

    def sync_to_state(
        self, rebuild_if_needed: bool = False, skip_extract: bool = False
    ) -> bool:
        """Sync current UI controls back to the central state.

        Args:
            rebuild_if_needed: If True, triggers a reconciliation / UI rebuild.
            skip_extract: If True, skips UI extraction and uses existing
                state directly.
        """
        if skip_extract:
            ui_primers = self.owner.input_data.primers
        else:
            ui_primers = self.extract_primer_data_from_ui()

        auto_activate_new = self.owner.settings.get(
            "auto_activate_new_valid_primer", False
        )

        primers = reconcile_primer_states(
            ui_primers,
            self.owner.input_data.primers,
            auto_activate_new=auto_activate_new,
        )

        for reconciled_p, ui_p in zip(primers, ui_primers, strict=True):
            checkbox = ui_p.get("checkbox")
            if checkbox:
                checkbox.value = reconciled_p["active"]

        ignore_inactive_name_dup = self.owner.settings.get(
            "ignore_inactive_name_dup_warn", True
        )
        ignore_inactive_seq_dup = self.owner.settings.get(
            "ignore_inactive_seq_dup_warn", True
        )

        dup_indices = get_duplicate_primer_indices(
            ui_primers, ignore_inactive_name_dup, ignore_inactive_seq_dup
        )
        for p in ui_primers:
            container = p.get("container")
            if container is None:
                continue

            c_idx = container.data
            is_dup = c_idx in dup_indices
            new_color = GUIColours.DUPLICATE_BG if is_dup else None
            if container.bgcolor != new_color:
                container.bgcolor = new_color

        # Run background primer construction/validation
        new_validation_errors = validate_primers(
            primers, ignore_inactive_name_dup, ignore_inactive_seq_dup
        )

        self.owner.input_data.primers = primers
        self.owner.validation_errors = new_validation_errors
        if rebuild_if_needed:
            self.owner.update_ui()
        else:
            # Update error status and duplicate highlights in-place on
            # existing controls
            for idx, row in enumerate(self.owner.primers_list.controls):
                if isinstance(row, PrimerRow) and idx < len(
                    new_validation_errors
                ):
                    row.set_error(new_validation_errors[idx])
                    row.update_tm(self.owner.settings)
            self.owner._update_row_highlights()
            self.owner._update_header_checkbox_state()
            self.owner._update_delete_button_disabled_state()

        return rebuild_if_needed
