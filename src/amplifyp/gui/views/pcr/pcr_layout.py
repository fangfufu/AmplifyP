# Copyright (C) 2026 Fufu Fang
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

from typing import Any

from amplifyp.pcr import PCR


class PCRLayoutSolver:
    """Helper class encapsulating coordinate layout calculations.

    Used by PCRDrawingPanel.
    """

    @staticmethod
    def collect_primer_bindings(
        pcr: PCR, amplicons: list[Any]
    ) -> tuple[
        dict[tuple[int, str], tuple[str, float, Any, Any]],
        dict[tuple[int, str], tuple[str, float, Any, Any]],
    ]:
        """Collect and group unique forward and reverse primer binding sites.

        Args:
            pcr: The PCR simulation instance containing repliconf configs.
            amplicons: List of amplicon objects with start/end positions.

        Returns:
            A tuple of (forward_bindings, reverse_bindings) dicts mapping
            (index, name) keys to (name, size, conf, var) tuples.
        """
        fwd_bindings = {}
        rev_bindings = {}
        for amp in amplicons:
            fwd_conf = next(
                (
                    c
                    for c in pcr.amplicon_generator.repliconfs
                    if c.primer is amp.fwd_origin
                ),
                None,
            )
            rev_conf = next(
                (
                    c
                    for c in pcr.amplicon_generator.repliconfs
                    if c.primer is amp.rev_origin
                ),
                None,
            )
            if fwd_conf is None or rev_conf is None:
                continue
            fwd_origin_point = fwd_conf.origin(amp.start)
            rev_origin_point = rev_conf.origin(amp.end)
            if fwd_origin_point is None or rev_origin_point is None:
                continue
            fwd_quality = fwd_origin_point.quality
            rev_quality = rev_origin_point.quality

            # Scale triangle size S based on quality score
            fwd_s = 6.0 + (max(0.1, min(1.0, fwd_quality)) * 10.0)
            rev_s = 6.0 + (max(0.1, min(1.0, rev_quality)) * 10.0)

            fwd_bindings[(amp.start.index, amp.fwd_origin.name)] = (
                amp.fwd_origin.name,
                fwd_s,
                fwd_conf,
                amp.start,
            )
            rev_bindings[(amp.end.index, amp.rev_origin.name)] = (
                amp.rev_origin.name,
                rev_s,
                rev_conf,
                amp.end,
            )
        return fwd_bindings, rev_bindings

    @staticmethod
    def calculate_shifted_x(
        bindings: dict[tuple[int, str], tuple[str, float, Any, Any]],
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
    def calculate_canvas_dimensions(
        target_length: int,
        num_amplicons: int,
        fwd_bindings: dict[tuple[int, str], tuple[str, float, Any, Any]],
        rev_bindings: dict[tuple[int, str], tuple[str, float, Any, Any]],
        page_width: float | None,
    ) -> tuple[float, float, float, float, float]:
        """Calculate drawing coordinates and baseline positions.

        Args:
            target_length: Total length of the template in base pairs.
            num_amplicons: Number of amplicons to render.
            fwd_bindings: Forward primer binding data dict.
            rev_bindings: Reverse primer binding data dict.
            page_width: Current page width in pixels, or None.

        Returns:
            A tuple of (v_target, h_margin, c_width, t_width, v_frag_start)
            representing vertical target position, horizontal margin, canvas
            width, template width, and fragment start vertical position.
        """
        max_fwd_len = max(
            (len(name) for name, _, _, _ in fwd_bindings.values()), default=0
        )
        max_rev_len = max(
            (len(name) for name, _, _, _ in rev_bindings.values()), default=0
        )

        fwd_px = max_fwd_len * 8.0
        rev_px = max_rev_len * 8.0

        v_target = max(100.0, 66.0 + fwd_px)
        h_margin = 20.0
        c_width = max(600.0, page_width - 80.0) if page_width else 800.0
        t_width = c_width - (2.0 * h_margin)

        v_frag_start = v_target + 70.0 + rev_px

        return v_target, h_margin, c_width, t_width, v_frag_start
