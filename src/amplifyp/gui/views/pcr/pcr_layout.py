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

"""Coordinate transformation and layout solver helper for the PCRView."""

import math
from typing import TYPE_CHECKING

from amplifyp.pcr import PCR

if TYPE_CHECKING:
    from amplifyp.amplicon import Amplicon
    from amplifyp.repliconf import DirIdx, Repliconf


class PCRLayoutSolver:
    """Helper class encapsulating coordinate layout calculations.

    Used by PCRDrawingPanel.
    """

    @staticmethod
    def collect_primer_bindings(
        pcr: PCR, amplicons: list["Amplicon"]
    ) -> tuple[
        dict[tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]],
        dict[tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]],
    ]:
        """Collect and group unique forward and reverse primer binding sites.

        Args:
            pcr: The PCR simulation instance containing repliconf configs.
            amplicons: List of amplicon objects with start/end positions.

        Returns:
            A tuple of (forward_bindings, reverse_bindings) dicts mapping
            (index, name) keys to (name, size, conf, var) tuples.
        """
        raw_fwd, raw_rev = PCR.collect_primer_bindings(pcr, amplicons)
        fwd_bindings = {
            key: (
                binding.name,
                6.0 + (max(0.1, min(1.0, binding.origin.quality)) * 10.0),
                binding.repliconf,
                binding.dir_idx,
            )
            for key, binding in raw_fwd.items()
        }
        rev_bindings = {
            key: (
                binding.name,
                6.0 + (max(0.1, min(1.0, binding.origin.quality)) * 10.0),
                binding.repliconf,
                binding.dir_idx,
            )
            for key, binding in raw_rev.items()
        }
        return fwd_bindings, rev_bindings

    @staticmethod
    def calculate_shifted_x(
        bindings: dict[
            tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]
        ],
        target_length: int,
        t_width: float,
        h_margin: float,
        min_dist: float = 24.0,
    ) -> dict[tuple[int, str], float]:
        """Calculate shifted horizontal pixel positions to prevent overlap.

        Detects clusters of bindings that are too close together and
        redistributes them evenly within the available canvas width.

        Args:
            bindings: Dict mapping (index, name) keys to binding data.
            target_length: Total length of the template in base pairs.
            t_width: Total template drawing width in pixels.
            h_margin: Horizontal margin in pixels.
            min_dist: Minimum pixel distance between shifted bindings.

        Returns:
            Dict mapping (index, name) keys to shifted x-coordinate floats.
        """
        if not bindings:
            return {}

        coords = []
        for key in bindings.keys():
            idx = key[0]
            x = (
                h_margin + (idx / target_length * t_width)
                if target_length
                else h_margin
            )
            coords.append((key, x))
        coords.sort(key=lambda item: item[1])

        n = len(coords)
        shifted = dict(coords)

        clusters = []
        current_cluster = [coords[0]]

        for i in range(1, n):
            key, x = coords[i]
            _, prev_x = coords[i - 1]
            if x - prev_x < min_dist:
                current_cluster.append((key, x))
            else:
                clusters.append(current_cluster)
                current_cluster = [(key, x)]
        clusters.append(current_cluster)

        for cluster in clusters:
            if len(cluster) <= 1:
                continue
            avg_x = sum(x for _, x in cluster) / len(cluster)

            k = len(cluster)
            total_width = (k - 1) * min_dist
            start_x = avg_x - total_width / 2.0

            min_canvas_x = h_margin
            max_canvas_x = h_margin + t_width
            if start_x < min_canvas_x:
                start_x = min_canvas_x
            elif start_x + total_width > max_canvas_x:
                start_x = max_canvas_x - total_width

            for idx_in_cluster, (key, _) in enumerate(cluster):
                shifted[key] = start_x + idx_in_cluster * min_dist

        return shifted

    @staticmethod
    def calculate_horizontal_label_layout(
        bindings: dict[
            tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]
        ],
        shifted_x: dict[tuple[int, str], float],
        c_width: float,
        font_size: int = 13,
        is_fwd: bool = True,
    ) -> tuple[dict[tuple[int, str], tuple[float, int, float]], int]:
        """Calculate horizontal x-positions and non-overlapping row indices.

        Places forward primer labels to the right of their triangles and
        reverse primer labels to the left of their triangles, clamped within
        the canvas bounds. Assigns row 0 when there is no overlap with other
        labels or primer triangles, otherwise stacks into higher rows.

        Args:
            bindings: Dict mapping (index, name) keys to binding data.
            shifted_x: Dict mapping (index, name) keys to shifted x-coords.
            c_width: Total canvas width in pixels.
            font_size: Font size in pixels for primer labels.
            is_fwd: True for forward primers (label on right), False for
                reverse primers (label on left).

        Returns:
            A tuple of (layout_map, num_rows) where layout_map maps
            (index, name) keys to (left, row, label_width) tuples and
            num_rows is the number of rows used.
        """
        if not bindings:
            return {}, 0

        char_width = font_size * 0.72
        gap = 8.0
        min_left = 4.0

        placements: list[
            tuple[float, float, float, float, tuple[int, str], float]
        ] = []
        for key, (name, s_val, _, _) in bindings.items():
            x_render = shifted_x.get(key, 0.0)
            label_width = max(20.0, len(name) * char_width + 4.0)
            max_left = max(min_left, c_width - label_width - 4.0)
            if is_fwd:
                desired_left = x_render + s_val / 2.0 + 4.0
                if desired_left > max_left:
                    alt_left = x_render - s_val / 2.0 - 4.0 - label_width
                    if alt_left >= min_left:
                        desired_left = alt_left
            else:
                desired_left = x_render - s_val / 2.0 - 4.0 - label_width
                if desired_left < min_left:
                    alt_left = x_render + s_val / 2.0 + 4.0
                    if alt_left <= max_left:
                        desired_left = alt_left
            left = max(min_left, min(max_left, desired_left))
            unit_left = min(x_render - s_val / 2.0, left)
            unit_right = max(x_render + s_val / 2.0, left + label_width)
            placements.append(
                (x_render, unit_left, unit_right, left, key, label_width)
            )

        # Forward labels extend right: place right-most primers in lower rows
        # so left-hand primers sit higher and lines avoid labels.
        # Reverse labels extend left: place left-most primers in lower rows
        # so right-hand primers sit further out and lines avoid labels.
        placements.sort(key=lambda item: item[0], reverse=is_fwd)

        row_intervals: list[list[tuple[float, float]]] = []
        layout_map: dict[tuple[int, str], tuple[float, int, float]] = {}

        for _, unit_left, unit_right, left, key, label_width in placements:
            assigned_row = -1
            for idx, intervals in enumerate(row_intervals):
                overlaps = any(
                    not (
                        unit_right + gap <= r_left or unit_left >= r_right + gap
                    )
                    for r_left, r_right in intervals
                )
                if not overlaps:
                    assigned_row = idx
                    intervals.append((unit_left, unit_right))
                    break
            if assigned_row == -1:
                assigned_row = len(row_intervals)
                row_intervals.append([(unit_left, unit_right)])
            layout_map[key] = (left, assigned_row, label_width)

        return layout_map, len(row_intervals)

    @staticmethod
    def calculate_canvas_dimensions(
        target_length: int,
        num_amplicons: int,
        fwd_bindings: dict[
            tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]
        ],
        rev_bindings: dict[
            tuple[int, str], tuple[str, float, "Repliconf", "DirIdx"]
        ],
        page_width: float | None,
        label_orientation: str = "Horizontal",
        font_size_primer: int = 13,
    ) -> tuple[float, float, float, float, float]:
        """Calculate drawing coordinates and baseline positions.

        Args:
            target_length: Total length of the template in base pairs.
            num_amplicons: Number of amplicons to render.
            fwd_bindings: Forward primer binding data dict.
            rev_bindings: Reverse primer binding data dict.
            page_width: Current page width in pixels, or None.
            label_orientation: Primer label orientation ('Horizontal', '45°',
                or 'Vertical').
            font_size_primer: Font size in pixels for primer labels.

        Returns:
            A tuple of (v_target, h_margin, c_width, t_width, v_frag_start)
            representing vertical target position, horizontal margin, canvas
            width, template width, and fragment start vertical position.
        """
        h_margin = 20.0
        c_width = max(600.0, page_width - 80.0) if page_width else 800.0

        if label_orientation in ("Vertical", "45°"):
            char_advance = font_size_primer * 0.72
            if label_orientation == "45°" and target_length > 0:
                cos_45 = math.cos(math.pi / 4)
                needed_margin = h_margin
                for (idx, _), (name, _, _, _) in fwd_bindings.items():
                    frac = max(0.0, min(1.0, idx / target_length))
                    if frac > 0.5:
                        proj_w = (
                            max(20.0, len(name) * char_advance + 4.0) * cos_45
                        )
                        req = (proj_w + 4.0 - (1.0 - frac) * c_width) / (
                            2.0 * frac - 1.0
                        )
                        needed_margin = max(needed_margin, req)
                for (idx, _), (name, _, _, _) in rev_bindings.items():
                    frac = max(0.0, min(1.0, idx / target_length))
                    if frac < 0.5:
                        proj_w = (
                            max(20.0, len(name) * char_advance + 4.0) * cos_45
                        )
                        req = (proj_w + 4.0 - frac * c_width) / (
                            1.0 - 2.0 * frac
                        )
                        needed_margin = max(needed_margin, req)
                h_margin = min(c_width * 0.3, needed_margin)

            t_width = c_width - (2.0 * h_margin)
            max_fwd_len = max(
                (len(name) for name, _, _, _ in fwd_bindings.values()),
                default=0,
            )
            max_rev_len = max(
                (len(name) for name, _, _, _ in rev_bindings.values()),
                default=0,
            )

            scale = (
                1.0
                if label_orientation == "Vertical"
                else math.sin(math.pi / 4)
            )
            fwd_px = max_fwd_len * char_advance * scale
            rev_px = max_rev_len * char_advance * scale

            v_target = max(100.0, 66.0 + fwd_px)
            v_frag_start = v_target + 70.0 + rev_px
        else:
            t_width = c_width - (2.0 * h_margin)
            fwd_shifted = PCRLayoutSolver.calculate_shifted_x(
                fwd_bindings, target_length, t_width, h_margin
            )
            rev_shifted = PCRLayoutSolver.calculate_shifted_x(
                rev_bindings, target_length, t_width, h_margin
            )
            _, num_fwd_rows = PCRLayoutSolver.calculate_horizontal_label_layout(
                fwd_bindings,
                fwd_shifted,
                c_width,
                font_size_primer,
                is_fwd=True,
            )
            _, num_rev_rows = PCRLayoutSolver.calculate_horizontal_label_layout(
                rev_bindings,
                rev_shifted,
                c_width,
                font_size_primer,
                is_fwd=False,
            )
            row_height = font_size_primer + 5.0
            v_target = max(100.0, 100.0 + max(0, num_fwd_rows - 1) * row_height)
            v_frag_start = (
                v_target + 52.0 + max(0, num_rev_rows - 1) * row_height
            )

        return v_target, h_margin, c_width, t_width, v_frag_start
