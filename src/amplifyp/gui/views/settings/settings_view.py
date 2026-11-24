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

"""Settings View for the Flet application."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING, Any, cast

import flet as ft

from amplifyp.gui.colours import GUIColours
from amplifyp.gui.settings import GUISettings

if TYPE_CHECKING:
    from amplifyp.gui.utils.gui_helpers import BorderedCheckbox
from amplifyp.gui.logger import reconfigure_logging
from amplifyp.gui.views.settings.designer_1d_tile import Designer1DTile
from amplifyp.gui.views.settings.designer_2d_tile import Designer2DTile
from amplifyp.gui.views.settings.diagnostics_tile import DiagnosticsTile
from amplifyp.gui.views.settings.dimer_tile import DimerTile
from amplifyp.gui.views.settings.general_tile import GeneralTile
from amplifyp.gui.views.settings.input_tile import InputTile
from amplifyp.gui.views.settings.replication_tile import ReplicationTile
from amplifyp.settings import ReplicationSettings


class SettingsView(ft.ListView):  # type: ignore[misc]
    """Settings view for configuring PCR replication settings."""

    def __init__(
        self,
        page: ft.Page,
        settings: GUISettings | None = None,
        on_change: Callable[[ft.ControlEvent], None] | None = None,
        on_reset: Callable[[ft.ControlEvent | None], None] | None = None,
        on_update_found: Callable[[str], None] | None = None,
    ) -> None:
        """Initialise the SettingsView."""
        super().__init__(
            expand=True, spacing=20, padding=10, scroll=ft.ScrollMode.ALWAYS
        )
        self.app_page = page
        self.settings = settings if settings is not None else GUISettings()
        self.on_change = on_change
        self.on_reset = on_reset
        self.on_update_found = on_update_found

        self.settings_map: dict[str, Any] = {}

        header_size = self.settings.get("font_size_header", 18)
        font_size_default = self.settings.get("font_size_default", 14)
        font_size_micro = self.settings.get("font_size_micro", 10)
        font_size_table_header = self.settings.get("font_size_table_header", 15)

        # Initialise sub-control tiles
        self.general_tile = GeneralTile(
            page=self.app_page,
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
            font_size_default=font_size_default,
            sync_to_state_callback=self.sync_to_state,
            update_ui_callback=self.update_ui,
            on_update_found=self.on_update_found,
            on_reset_callback=self._on_reset_handler,
            on_external_load_callback=self._on_external_settings_loaded,
        )

        self.input_tile = InputTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,
            header_size=header_size,
            font_size_default=font_size_default,
        )

        self.replication_tile = ReplicationTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        self.dimer_tile = DimerTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
            font_size_default=font_size_default,
            font_size_micro=font_size_micro,
            font_size_table_header=font_size_table_header,
        )

        self.designer_1d_tile = Designer1DTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
        )

        self.designer_2d_tile = Designer2DTile(
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
        )

        self.diagnostics_tile = DiagnosticsTile(
            page=self.app_page,
            settings=self.settings,
            settings_map=self.settings_map,
            on_change_handler=self._on_change_handler,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
            header_size=header_size,
        )

        no_border = ft.RoundedRectangleBorder(side=ft.BorderSide.none())
        tiles: list[ft.ExpansionTile] = [
            self.general_tile,
            self.input_tile,
            self.replication_tile,
            self.dimer_tile,
            self.designer_1d_tile,
            self.designer_2d_tile,
            self.diagnostics_tile,
        ]
        self.controls = []
        for index, tile in enumerate(tiles):
            tile.shape = no_border
            tile.collapsed_shape = no_border
            self.controls.append(tile)
            if index < len(tiles) - 1:
                self.controls.append(
                    ft.Divider(height=1, color=GUIColours.OUTLINE)
                )
        self.controls.extend(
            [
                ft.Divider(height=1, color=GUIColours.OUTLINE),
                self._build_reset_button(),
            ]
        )

        # Sync initial UI state
        self.update_ui()

    # Backwards compatibility properties for unit tests
    @property
    def set_primability_cutoff(self) -> ft.TextField:
        """Get the primability cutoff text field."""
        return self.replication_tile.set_primability_cutoff

    @property
    def set_stability_cutoff(self) -> ft.TextField:
        """Get the stability cutoff text field."""
        return self.replication_tile.set_stability_cutoff

    @property
    def set_amp4_compat(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the amplify4 compatibility mode checkbox."""
        return self.replication_tile.set_amp4_compat

    @property
    def set_tm_dna_conc(self) -> ft.TextField:
        """Get the DNA conc text field."""
        return self.input_tile.set_tm_dna_conc

    @property
    def set_tm_dnap_conc(self) -> ft.TextField:
        """Get the DNA Pol conc text field."""
        return self.input_tile.set_tm_dnap_conc

    @property
    def set_tm_mono_salt(self) -> ft.TextField:
        """Get the monovalent salt conc text field."""
        return self.input_tile.set_tm_mono_salt

    @property
    def set_tm_div_salt(self) -> ft.TextField:
        """Get the divalent salt conc text field."""
        return self.input_tile.set_tm_div_salt

    @property
    def set_tm_dNTP_conc(self) -> ft.TextField:
        """Get the dNTP conc text field."""
        return self.input_tile.set_tm_dNTP_conc

    @property
    def set_tm_method(self) -> ft.Dropdown:
        """Get the Tm method calculation dropdown."""
        return self.input_tile.set_tm_method

    @property
    def set_pd_min_overlap(self) -> ft.TextField:
        """Get the dimer min overlap text field."""
        return self.dimer_tile.set_pd_min_overlap

    @property
    def set_pd_threshold(self) -> ft.TextField:
        """Get the dimer threshold text field."""
        return self.dimer_tile.set_pd_threshold

    @property
    def set_font_family(self) -> ft.Dropdown:
        """Get the font family dropdown."""
        return self.general_tile.set_font_family

    @property
    def set_colour_scheme(self) -> ft.Dropdown:
        """Get the colour scheme dropdown."""
        return self.general_tile.set_colour_scheme

    @property
    def set_dimer_sequence_separator(self) -> ft.Dropdown:
        """Get the dimer sequence separator dropdown."""
        return self.general_tile.set_dimer_sequence_separator

    @property
    def set_sequence_separator(self) -> ft.Dropdown:
        """Get the 5'/3' sequence separator dropdown."""
        return self.general_tile.set_sequence_separator

    @property
    def set_colour_deficient(self) -> ft.Checkbox:
        """Get the colour deficient mode checkbox."""
        return self.general_tile.set_colour_deficient

    @property
    def set_improved_visualisation(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the improved visualisation mode checkbox."""
        return self.replication_tile.set_improved_visualisation

    @property
    def set_show_primer_temperature(self) -> ft.Checkbox | BorderedCheckbox:
        """Get the show primer temperature checkbox."""
        return self.input_tile.set_show_primer_temperature

    @property
    def set_tm_colour_scheme(self) -> ft.Dropdown:
        """Get the Tm colour scheme dropdown."""
        return self.input_tile.set_tm_colour_scheme

    @property
    def set_designer_2d_colour_scheme(self) -> ft.Dropdown:
        """Get the 2D Designer colour scheme dropdown."""
        return self.designer_2d_tile.set_designer_2d_colour_scheme

    @property
    def set_designer_2d_show_rev_fwd(self) -> BorderedCheckbox:
        """Get the 2D Designer show reverse-forward cross-dimer checkbox."""
        return self.designer_2d_tile.set_designer_2d_show_rev_fwd

    @property
    def set_version_checking_frequency(self) -> ft.Dropdown:
        """Get the version checking frequency dropdown."""
        return self.general_tile.set_version_checking_frequency

    @property
    def set_auto_reload_on_startup(self) -> BorderedCheckbox:
        """Get the auto reload on startup checkbox."""
        return self.general_tile.set_auto_reload_on_startup

    @property
    def set_primer_info_panel_position(self) -> ft.Dropdown:
        """Get the primer info panel position dropdown."""
        return self.input_tile.set_primer_info_panel_position

    @property
    def set_primer_info_panel_fixed_height(self) -> BorderedCheckbox:
        """Get the primer info panel fixed height checkbox."""
        return self.input_tile.fixed_height_primer_info_checkbox

    @property
    def set_auto_activate_new_valid_primer(self) -> BorderedCheckbox:
        """Get the auto activate new valid primer checkbox."""
        return self.input_tile.auto_activate_new_valid_primer

    def _build_reset_button(self) -> ft.Row:
        """Build the Reset button Row."""
        return ft.Row(
            list[ft.Control](
                [
                    ft.OutlinedButton(
                        "Reset to Default",
                        icon=ft.Icons.RESTORE,
                        on_click=self.general_tile._on_factory_reset_click,  # pyright: ignore[reportArgumentType, reportAttributeAccessIssue]
                    ),
                ]
            ),
            alignment=ft.MainAxisAlignment.END,
            spacing=10,
        )

    def sync_to_state(self) -> None:
        """Sync current settings UI controls back to the central state."""
        for prefix in ("match_weight_", "run_weight_"):
            stale_keys = [
                k
                for k in self.settings.keys()
                if k.startswith(prefix)
                and k != f"{prefix}default"
                and k not in self.settings_map
            ]
            for k in stale_keys:
                del self.settings[k]
        for k, v in self.settings_map.items():
            if k == "log_rotation_max_bytes":
                try:
                    self.settings[k] = int(v.value) * 1024 * 1024
                except (ValueError, TypeError):
                    self.settings[k] = 5242880
            else:
                self.settings[k] = v.value
        self.general_tile.sync_colour_scheme_to_settings()

    def _reconfigure_logging(self) -> None:
        """Reconfigure logging based on current settings."""
        reconfigure_logging(
            log_level_amplifyp=self.settings.get("log_level_amplifyp", "DEBUG"),
            log_level_flet=self.settings.get("log_level_flet", "INFO"),
            log_console_enabled=self.settings.get("log_console_enabled", True),
            log_file_enabled=self.settings.get("log_file_enabled", True),
            log_file_path=self.settings.get("log_file_path", "(Default)"),
            is_web=self.app_page.web,
            log_rotation_enabled=self.settings.get(
                "log_rotation_enabled", True
            ),
            log_rotation_max_bytes=self.settings.get(
                "log_rotation_max_bytes", 5242880
            ),
        )

    def update_ui(self) -> None:
        """Update Flet UI controls to match the central settings."""
        from amplifyp.gui.utils.gui_helpers import BorderedCheckbox

        for k, v in self.settings.items():
            if k in self.settings_map:
                ctrl = self.settings_map[k]
                if isinstance(ctrl, (ft.Checkbox, BorderedCheckbox)):
                    ctrl.value = bool(v)
                    if isinstance(ctrl, BorderedCheckbox):
                        ctrl.border = ft.Border.all(1, GUIColours.OUTLINE)
                else:
                    ctrl.value = str(v)
                    if isinstance(ctrl, ft.Dropdown):
                        ctrl.border_color = GUIColours.OUTLINE
                    elif isinstance(ctrl, ft.TextField):
                        if isinstance(ctrl.border, ft.NoInputBorder):
                            ctrl.text_style = ft.TextStyle(
                                color=GUIColours.DIAGRAM_BLACK,
                                size=self.settings.get("font_size_default", 14),
                            )
                        else:
                            ctrl.border_color = GUIColours.OUTLINE
        for tile_ctrl in self.controls:
            if isinstance(tile_ctrl, ft.Divider):
                tile_ctrl.color = GUIColours.OUTLINE
        self.general_tile.update_ui()
        self.replication_tile.update_ui()
        self.dimer_tile.update_ui()
        self.diagnostics_tile.update_ui()

    def _on_change_handler(self, e: ft.ControlEvent) -> None:
        """Handle change in settings fields."""
        # Sync the triggering control's value first (Dropdown .value may not
        # be updated yet when on_select fires in Flet).
        ctrl = getattr(e, "control", None)
        if ctrl is not None:
            for k, v in self.settings_map.items():
                if v is ctrl:
                    self.settings[k] = ctrl.value
                    break
        self.sync_to_state()
        self._reconfigure_logging()
        if self.on_change:
            self.on_change(e)

    def _on_external_settings_loaded(self) -> None:
        """Apply externally loaded settings without overwriting from widgets."""
        self._reconfigure_logging()
        if self.on_change:
            self.on_change(cast(Any, None))
        else:
            self.settings.save_to_local(self.app_page)
            self.app_page.update()

    def _on_reset_handler(self, *args: Any) -> None:
        """Handle reset to default button click."""
        self.settings.reset_to_defaults()
        self._reconfigure_logging()
        if self.on_reset:
            self.on_reset(args[0] if args else None)
        self.update_ui()
        if not self.on_reset or self.on_reset != self.on_change:
            self.settings.save_to_local(self.app_page)
            self.app_page.update()

    def get_replication_settings(self) -> ReplicationSettings:
        """Get the current settings as a ReplicationSettings object."""
        self.sync_to_state()
        return self.settings.get_replication_settings()

    def get_state(self) -> dict[str, Any]:
        """Get the current settings for serialisation."""
        self.sync_to_state()
        return self.settings.to_dict()

    def set_state(self, state: dict[str, Any]) -> None:
        """Set the current settings from deserialized data."""
        self.settings.from_dict(state)
        self.update_ui()
        self.app_page.update()

    @property
    def set_designer_1d_show_tm(self) -> BorderedCheckbox:
        """Get the show melting temperature checkbox."""
        return self.designer_1d_tile.set_designer_1d_show_tm

    @property
    def set_designer_1d_show_pct_at(self) -> BorderedCheckbox:
        """Get the show % AT checkbox."""
        return self.designer_1d_tile.set_designer_1d_show_pct_at

    @property
    def set_pcr_amplicon_ranking(self) -> ft.Dropdown:
        """Get the amplicon vertical ranking dropdown."""
        return self.replication_tile.set_amplicon_ranking

    @property
    def set_pcr_primer_label_orientation(self) -> ft.Dropdown:
        """Get the origin map primer label orientation dropdown."""
        return self.replication_tile.set_primer_label_orientation

    @property
    def set_pcr_divider_behaviour(self) -> ft.Dropdown:
        """Get the PCR horizontal divider behaviour dropdown."""
        return self.replication_tile.set_pcr_divider_behaviour
