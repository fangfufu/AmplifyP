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

"""Settings dialog for the PySide6 application."""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from amplifyp.gui2.settings import GUISettings

if TYPE_CHECKING:
    pass


class SettingsView(QDialog):
    """Settings dialog with tabs for different configuration sections."""

    def __init__(
        self,
        settings: GUISettings,
        on_change: Callable[[], None],
        on_reset: Callable[[], None],
        on_update_found: Callable[[str], None],
    ) -> None:
        """Initialize the Settings dialog."""
        super().__init__()
        self.settings = settings
        self.on_change = on_change
        self.on_reset = on_reset
        self.on_update_found = on_update_found
        self.setWindowTitle("Settings")
        self.setMinimumSize(600, 500)
        self._build_ui()

    def _build_ui(self) -> None:
        """Build the settings dialog UI."""
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # Tab widget
        self.tabs = QTabWidget()
        self.tabs.addTab(self._create_general_tab(), "General")
        self.tabs.addTab(self._create_tm_tab(), "Melting Temp")
        self.tabs.addTab(self._create_replication_tab(), "Replication")
        self.tabs.addTab(self._create_dimer_tab(), "Primer Dimers")
        self.tabs.addTab(self._create_primer_list_tab(), "Primer List")
        self.tabs.addTab(self._create_designer_2d_tab(), "Designer 2D")
        self.tabs.addTab(self._create_diagnostics_tab(), "Diagnostics")

        layout.addWidget(self.tabs)

        # Sync/Reset buttons
        btn_layout = QHBoxLayout()
        sync_btn = QPushButton("Sync to PCR/Dimers")
        sync_btn.clicked.connect(self.on_change)
        reset_btn = QPushButton("Reset to Defaults")
        reset_btn.clicked.connect(self.on_reset)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)

        btn_layout.addWidget(sync_btn)
        btn_layout.addWidget(reset_btn)
        btn_layout.addStretch()
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

    def _create_general_tab(self) -> QWidget:
        """Create the general settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        # Font family
        self.font_family = QLineEdit(
            self.settings.get("font_family", "Roboto Mono")
        )
        self.font_family.textChanged.connect(self.on_change)
        layout.addRow("Font Family:", self.font_family)

        # Dark mode
        dark_mode_val = self.settings.get("dark_mode", "system")
        self.dark_mode = QCheckBox("Dark Mode")
        if str(dark_mode_val).lower() not in ("false", "0", "no", "system"):
            self.dark_mode.setChecked(True)
        self.dark_mode.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Dark Mode:", self.dark_mode)

        # Colour deficient mode
        self.colour_deficient = QCheckBox("Colour Deficient Mode")
        if self.settings.get("colour_deficient", False):
            self.colour_deficient.setChecked(True)
        self.colour_deficient.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Colour Deficient:", self.colour_deficient)

        # Auto reload
        self.auto_reload = QCheckBox("Auto-reload on startup")
        if self.settings.get("auto_reload_on_startup", True):
            self.auto_reload.setChecked(True)
        self.auto_reload.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Auto-reload:", self.auto_reload)

        # Backup save/load
        backup_btn_layout = QHBoxLayout()
        save_backup_btn = QPushButton("Save Backup")
        save_backup_btn.clicked.connect(self._save_backup)
        load_backup_btn = QPushButton("Load Backup")
        load_backup_btn.clicked.connect(self._load_backup)
        backup_btn_layout.addWidget(save_backup_btn)
        backup_btn_layout.addWidget(load_backup_btn)
        layout.addRow("Backup:", backup_btn_layout)

        return widget

    def _create_tm_tab(self) -> QWidget:
        """Create the melting temperature settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        # TM method
        self.tm_method = QLineEdit(
            self.settings.get(
                "tm_method", "SantaLucia 1998 / Owczarzy 2008 (Default)"
            )
        )
        self.tm_method.textChanged.connect(self.on_change)
        layout.addRow("TM Method:", self.tm_method)

        # DNA concentration
        self.tm_dna_conc = QLineEdit(
            str(self.settings.get("tm_dna_conc", 50.0))
        )
        self.tm_dna_conc.textChanged.connect(self.on_change)
        layout.addRow("DNA Conc (nM):", self.tm_dna_conc)

        # Monovalent salt
        self.tm_mono_salt = QLineEdit(
            str(self.settings.get("tm_mono_salt", 50.0))
        )
        self.tm_mono_salt.textChanged.connect(self.on_change)
        layout.addRow("Monovalent Salt (mM):", self.tm_mono_salt)

        # Divalent salt
        self.tm_div_salt = QLineEdit(str(self.settings.get("tm_div_salt", 1.5)))
        self.tm_div_salt.textChanged.connect(self.on_change)
        layout.addRow("Divalent Salt (mM):", self.tm_div_salt)

        # dNTP
        self.tm_dNTP_conc = QLineEdit(
            str(self.settings.get("tm_dNTP_conc", 0.8))
        )
        self.tm_dNTP_conc.textChanged.connect(self.on_change)
        layout.addRow("dNTP (mM):", self.tm_dNTP_conc)

        return widget

    def _create_replication_tab(self) -> QWidget:
        """Create the replication settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.primability_cutoff = QLineEdit(
            str(self.settings.get("primability_cutoff", 0.7))
        )
        self.primability_cutoff.textChanged.connect(self.on_change)
        layout.addRow("Primability Cutoff:", self.primability_cutoff)

        self.stability_cutoff = QLineEdit(
            str(self.settings.get("stability_cutoff", 0.0))
        )
        self.stability_cutoff.textChanged.connect(self.on_change)
        layout.addRow("Stability Cutoff:", self.stability_cutoff)

        self.amp4_compat = QCheckBox("Amplify4 Compatibility")
        if self.settings.get("amp4_compat", False):
            self.amp4_compat.setChecked(True)
        self.amp4_compat.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Amp4 Compat:", self.amp4_compat)

        self.improved_vis = QCheckBox("Improved Visualisation")
        if self.settings.get("improved_visualisation", True):
            self.improved_vis.setChecked(True)
        self.improved_vis.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Improved Visualisation:", self.improved_vis)

        return widget

    def _create_dimer_tab(self) -> QWidget:
        """Create the primer dimer settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.pd_min_overlap = QLineEdit(
            str(self.settings.get("pd_min_overlap", 8))
        )
        self.pd_min_overlap.textChanged.connect(self.on_change)
        layout.addRow("Min Overlap:", self.pd_min_overlap)

        self.pd_threshold = QLineEdit(
            str(self.settings.get("pd_threshold", 0.0))
        )
        self.pd_threshold.textChanged.connect(self.on_change)
        layout.addRow("Threshold:", self.pd_threshold)

        return widget

    def _create_primer_list_tab(self) -> QWidget:
        """Create the primer list settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.show_tm = QCheckBox("Show Primer Tm")
        if self.settings.get("show_primer_temperature", False):
            self.show_tm.setChecked(True)
        self.show_tm.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Show Tm:", self.show_tm)

        self.tm_colour_scheme = QLineEdit(
            self.settings.get("tm_colour_scheme", "None")
        )
        self.tm_colour_scheme.textChanged.connect(self.on_change)
        layout.addRow("Tm Colour Scheme:", self.tm_colour_scheme)

        self.ignore_inactive_name = QCheckBox("Ignore inactive name duplicates")
        if self.settings.get("ignore_inactive_name_dup_warn", True):
            self.ignore_inactive_name.setChecked(True)
        self.ignore_inactive_name.stateChanged.connect(
            lambda s: self.on_change()
        )
        layout.addRow("Ignore inactive name dups:", self.ignore_inactive_name)

        self.ignore_inactive_seq = QCheckBox("Ignore inactive seq duplicates")
        if self.settings.get("ignore_inactive_seq_dup_warn", True):
            self.ignore_inactive_seq.setChecked(True)
        self.ignore_inactive_seq.stateChanged.connect(
            lambda s: self.on_change()
        )
        layout.addRow("Ignore inactive seq dups:", self.ignore_inactive_seq)

        return widget

    def _create_designer_2d_tab(self) -> QWidget:
        """Create the designer 2D settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.designer_2d_colour_scheme = QLineEdit(
            self.settings.get("designer_2d_colour_scheme", "Blue-Orange")
        )
        self.designer_2d_colour_scheme.textChanged.connect(self.on_change)
        layout.addRow("2D Colour Scheme:", self.designer_2d_colour_scheme)

        return widget

    def _create_diagnostics_tab(self) -> QWidget:
        """Create the diagnostics settings tab."""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.log_level_amplifyp = QLineEdit(
            self.settings.get("log_level_amplifyp", "INFO")
        )
        self.log_level_amplifyp.textChanged.connect(self.on_change)
        layout.addRow("AmplifyP Log Level:", self.log_level_amplifyp)

        self.log_console = QCheckBox("Enable Console Logging")
        if self.settings.get("log_console_enabled", True):
            self.log_console.setChecked(True)
        self.log_console.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Console Logging:", self.log_console)

        self.log_file = QCheckBox("Enable File Logging")
        if self.settings.get("log_file_enabled", True):
            self.log_file.setChecked(True)
        self.log_file.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("File Logging:", self.log_file)

        self.log_rotation = QCheckBox("Enable Log Rotation")
        if self.settings.get("log_rotation_enabled", True):
            self.log_rotation.setChecked(True)
        self.log_rotation.stateChanged.connect(lambda s: self.on_change())
        layout.addRow("Log Rotation:", self.log_rotation)

        self.log_max_bytes = QLineEdit(
            str(self.settings.get("log_rotation_max_bytes", 5242880))
        )
        self.log_max_bytes.textChanged.connect(self.on_change)
        layout.addRow("Max Log Size (bytes):", self.log_max_bytes)

        self.log_file_path = QLineEdit(
            self.settings.get("log_file_path", "(Default)")
        )
        self.log_file_path.textChanged.connect(self.on_change)
        layout.addRow("Log File Path:", self.log_file_path)

        pick_log_btn = QPushButton("Pick Log File")
        pick_log_btn.clicked.connect(self._pick_log_file)
        layout.addRow("", pick_log_btn)

        return widget

    def _save_backup(self) -> None:
        """Save settings to a backup file."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Settings Backup", "", "YAML Files (*.yaml *.yml)"
        )
        if path:
            import yaml

            with open(path, "w", encoding="utf-8") as f:
                yaml.safe_dump(self.settings.to_dict(), f)

    def _load_backup(self) -> None:
        """Load settings from a backup file."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Settings Backup", "", "YAML Files (*.yaml *.yml)"
        )
        if path:
            import yaml

            with open(path, encoding="utf-8") as f:
                data = yaml.safe_load(f)
            if isinstance(data, dict):
                self.settings.from_dict(data)
                self.on_change()

    def _pick_log_file(self) -> None:
        """Pick a log file path."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Select Log File", "", "Log Files (*.log)"
        )
        if path:
            self.log_file_path.setText(path)
            self.on_change()

    def sync_to_state(self) -> None:
        """Sync UI controls to settings."""

    def update_ui(self) -> None:
        """Update UI from settings."""
