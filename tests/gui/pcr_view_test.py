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

"""Tests for PCR View and primer binding site context map popups."""

import math
from typing import Any
from unittest.mock import MagicMock, patch

import flet as ft
import pytest

from amplifyp.dna import DNA, DNADirection, Primer
from amplifyp.gui.colours import GUIColours
from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.pcr import PCRView
from amplifyp.gui.views.pcr.primer_drawing import format_context_lines
from amplifyp.repliconf import Repliconf


def test_pcr_view_click_context_map() -> None:
    """Test that clicking a primer binding site shows a context map."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "tTccACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "10290", "seq": "tTccACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "gTgcAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # The diagram_stack controls should contain some ft.GestureDetectors
    # representing click overlays
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
    ]
    assert len(gesture_detectors) >= 2

    # Verify no cards initially
    assert len(view.result_list.controls) == 0

    # Click on the first gesture detector (forward binding site)
    fwd_detector = gesture_detectors[0]
    assert fwd_detector.on_tap is not None

    # Trigger click
    fwd_detector.on_tap(MagicMock())

    # Verify that a card was added below the overview map
    assert len(view.result_list.controls) == 1
    card = view.result_list.controls[0]
    assert isinstance(card, ft.Card)

    # Verify title text inside card
    title_text = card.content.content.controls[0].controls[0].content.value
    assert "Context Map" in title_text
    assert "10290" in title_text

    # Verify Primeability, Stability and Quality in title row
    title_row = card.content.content.controls[0]
    assert isinstance(title_row, ft.Row)
    metrics_row = title_row.controls[1]
    assert isinstance(metrics_row, ft.Row)
    p_text = metrics_row.controls[0].content.value
    s_text = metrics_row.controls[1].content.value
    q_text = metrics_row.controls[2].content.value
    assert "Primeability: 1.000" in p_text
    assert "Stability: 1.000" in s_text
    assert q_text == "Quality: 1.000"

    # Extract diagram_text
    diagram_text = card.content.content.controls[1].content.controls[0]
    assert isinstance(diagram_text, ft.Text)

    # Check spans content
    text_spans = [span.text for span in diagram_text.spans]
    full_text = "".join(text_spans)

    # Context map should have:
    # 1. 1-indexed positions
    # 2. Down arrows (↓)
    # 3. Primer name (10290) and 5'/3' markers
    # 4. Watson-Crick bonds (||||)
    # 5. Template context sequence
    assert "1" in full_text
    assert "20" in full_text
    assert "V" in full_text
    assert "10290" in full_text
    assert "5'" in full_text
    assert "3'" in full_text
    assert "|||" in full_text
    assert "Context" in full_text
    assert (
        "--------------------" in full_text
    )  # Upstream region (gaps because index 0 is at boundary)
    assert (
        "tTccACTGCGAATCATTAAA" in full_text
    )  # Binding region - must preserve mixed case


def test_pcr_view_click_amplicon() -> None:
    """Test that clicking on an amplicon bar inserts the amplicon details

    card below the map.
    """
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "tTccACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "10290", "seq": "tTccACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "gTgcAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # Verify no cards initially
    assert len(view.result_list.controls) == 0

    # Find the gesture detectors (excluding label detectors which wrap ft.Text).
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Container)
    ]
    # We should have 3 gesture detectors total (2 primers + 1 amplicon)
    assert len(gesture_detectors) == 3

    # Click on the third gesture detector (the amplicon)
    amp_detector = gesture_detectors[2]
    assert amp_detector.on_tap is not None

    # Trigger click
    amp_detector.on_tap(MagicMock())

    # Verify that a card was added below the overview map
    assert len(view.result_list.controls) == 1
    card = view.result_list.controls[0]
    assert isinstance(card, ft.Card)

    # Check contents: should have sequence length and labels
    column = card.content.content
    assert isinstance(column, ft.Column)

    length_text = column.controls[0].controls[0].content.value
    assert "Amplicon: 60 bp" in length_text

    subtitle_text = column.controls[1]
    assert isinstance(subtitle_text, ft.Text)
    spans_content = "".join([span.text for span in subtitle_text.spans])
    assert "— 20 bp —" in spans_content


def test_pcr_view_resize_preserves_cards() -> None:
    """Test that resizing the window does not discard open cards."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "10290", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # The diagram_stack controls should contain some ft.GestureDetectors
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
    ]
    assert len(gesture_detectors) >= 2

    # Click on the first gesture detector (forward binding site)
    fwd_detector = gesture_detectors[0]
    fwd_detector.on_tap(MagicMock())

    # Verify card is present
    assert len(view.result_list.controls) == 1
    assert isinstance(view.result_list.controls[0], ft.Card)

    # Trigger resize event
    view._handle_resize(MagicMock(spec=ft.ControlEvent))

    # Verify that the card is still present
    assert len(view.result_list.controls) == 1
    assert isinstance(view.result_list.controls[0], ft.Card)


def test_pcr_view_no_duplicate_cards() -> None:
    """Test that clicking an element twice does not duplicate its card

    but instead brings the existing card to the top of the list.
    """
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "10290", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # Find the gesture detectors: 2 context maps, 1 amplicon (excluding labels)
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Container)
    ]
    assert len(gesture_detectors) == 3

    # Click first context map (forward binding site)
    gesture_detectors[0].on_tap(MagicMock())
    assert len(view.result_list.controls) == 1
    card1 = view.result_list.controls[0]

    # Click second context map (reverse binding site)
    gesture_detectors[1].on_tap(MagicMock())
    assert len(view.result_list.controls) == 2
    assert view.result_list.controls[0] != card1

    # Click first context map again
    gesture_detectors[0].on_tap(MagicMock())
    # Should still be 2 cards total, and card1 should be moved back to the top
    assert len(view.result_list.controls) == 2
    assert view.result_list.controls[0] == card1

    # Click amplicon
    gesture_detectors[2].on_tap(MagicMock())
    assert len(view.result_list.controls) == 3
    amp_card = view.result_list.controls[0]

    # Click amplicon again
    gesture_detectors[2].on_tap(MagicMock())
    assert len(view.result_list.controls) == 3
    assert view.result_list.controls[0] == amp_card


def test_pcr_view_primer_label_collision() -> None:
    """Test that closely spaced primer labels are shifted horizontally."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    # Three forward-binding primers close to each other
    # and one reverse primer to allow amplicon rendering.
    input_data.primers = [
        {"name": "P1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "P2", "seq": "TCCACTGCGAATCATTAAAG", "active": True},
        {"name": "P3", "seq": "CCACTGCGAATCATTAAAGT", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # The stack should have DrawnPrimer text controls, which are now
    # wrapped in GestureDetector.
    texts = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Text)
        and ctrl.content.color == "blue800"
    ]
    # We should have 3 forward primer label texts
    assert len(texts) == 3

    # Extract their "left" values and sort them
    lefts = sorted([txt.left for txt in texts])

    # Separated by >= 24.0 pixels due to de-collision logic
    assert lefts[1] - lefts[0] >= 24.0
    assert lefts[2] - lefts[1] >= 24.0


def test_pcr_view_primers_same_3prime_end_different_5prime() -> None:
    """Test that primers with same 3' ends but different 5' ends are shown."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    # P1 and P2 have same 3' end sequence (ACTGCGAATCATTAAA) but
    # different 5' ends.
    input_data.primers = [
        {"name": "P1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "P2", "seq": "ACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # Verify that both P1 and P2 labels are rendered
    texts = [
        ctrl.content.value
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Text)
        and ctrl.content.color == "blue800"
    ]
    assert "P1" in texts
    assert "P2" in texts


def test_pcr_view_click_context_map_improved_visualisation() -> None:
    """Test clicking a primer binding site shows RC strand in context map."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "tTccACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "10290", "seq": "tTccACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "gTgcAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    # Turn on improved visualisation
    view.settings["improved_visualisation"] = True
    view.run_pcr()

    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
    ]
    assert len(gesture_detectors) >= 2

    # Click on the first gesture detector (forward binding site)
    fwd_detector = gesture_detectors[0]
    fwd_detector.on_tap(MagicMock())

    assert len(view.result_list.controls) == 1
    card = view.result_list.controls[0]
    assert isinstance(card, ft.Card)

    # Extract diagram_text
    diagram_text = card.content.content.controls[1].content.controls[0]
    assert isinstance(diagram_text, ft.Text)

    # Check spans
    # It should have a span for the comp_line with MUTED_GREY colour
    comp_spans = [
        span
        for span in diagram_text.spans
        if getattr(span, "style", None)
        and span.style.color == GUIColours.MUTED_GREY
    ]
    assert len(comp_spans) == 1
    # Check that it contains "3'-" and "-5'" and the translated comp sequence
    assert "3'-" in comp_spans[0].text
    assert "-5'" in comp_spans[0].text


def test_format_context_lines_alignment_long_label() -> None:
    """Test format_context_lines aligns elements for long primer labels."""
    template = DNA("A" * 100)
    primer = Primer("C" * 20)
    conf = Repliconf(template, primer)
    origin = MagicMock()
    # Mock binding strength string
    origin.binding_strength_str = "|" * 20

    # Long name of length 34
    long_name = "2223 - 108 (a.k.a. 3312) (Reverse)"
    _top, mid, bot = format_context_lines(
        primer_name=long_name,
        padded_idx=50,
        conf=conf,
        origin=origin,
        L=20,
        N=100,
        direction=DNADirection.REV,
    )

    # For DNADirection.REV, the primer label is f"{primer_name} (Reverse)"
    # which is "2223 - 108 (a.k.a. 3312) (Reverse) (Reverse)".
    # Let's count length: 34 + 10 = 44 characters.
    # The label is padded to 44. The sequence prefix is 44 + 3 = 47.
    # Verify that primer sequence starts with '3\'-' at offset 44
    # (default dash)
    assert mid.startswith("2223 - 108 (a.k.a. 3312) (Reverse) (Reverse)3'-")

    # The mid line should be:
    # "2223 - 108 (a.k.a. 3312) (Reverse) (Reverse)3'-CCC...C-5'"
    # The bonds line in bottom line is the first line of bottom_line.
    # The bonds line starts with:
    # 12 + 20 + extra_spaces = 32 + (44 - 29) = 47 spaces.
    lines_bot = bot.split("\n")
    assert lines_bot[0].startswith(" " * 47 + "|")

    # The context line in bottom_line starts with:
    # "Context  " (9 chars) + 15 spaces + "5'-" (3 chars) = 27 spaces,
    # then 20 bp upstream = 47 spaces before binding sequence.
    assert lines_bot[1].startswith("Context  " + " " * 15 + "5'-")

    # Test space separator option
    _top_s_opt, mid_s_opt, bot_s_opt = format_context_lines(
        primer_name=long_name,
        padded_idx=50,
        conf=conf,
        origin=origin,
        L=20,
        N=100,
        direction=DNADirection.REV,
        separator="Space (' ')",
    )
    assert mid_s_opt.startswith(
        "2223 - 108 (a.k.a. 3312) (Reverse) (Reverse)3' "
    )
    lines_bot_s_opt = bot_s_opt.split("\n")
    assert lines_bot_s_opt[1].startswith("Context  " + " " * 15 + "5' ")

    # Test reading separator from settings when separator is None
    from amplifyp.gui.settings import GUISettings

    _top_s, mid_s, _bot_s = format_context_lines(
        primer_name=long_name,
        padded_idx=50,
        conf=conf,
        origin=origin,
        L=20,
        N=100,
        direction=DNADirection.REV,
        settings=GUISettings({"sequence_separator": "Space (' ')"}),
    )
    assert mid_s.startswith("2223 - 108 (a.k.a. 3312) (Reverse) (Reverse)3' ")


def test_pcr_view_shows_binding_sites_when_no_amplicons_found() -> None:
    """Test that PCRView shows binding sites when 0 amplicons are found."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    # Configure 1 forward primer (no reverse primer -> 0 amplicons)
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # Verify that "No amplicons found." text is in result_list controls
    assert len(view.result_list.controls) == 1
    no_amp_text = view.result_list.controls[0]
    assert isinstance(no_amp_text, ft.Text)
    assert no_amp_text.value == "No amplicons found."

    # Verify diagram container is visible
    assert view.diagram_panel.diagram_container.visible is True

    # Verify gesture detectors (binding sites) exist on diagram stack
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
    ]
    assert len(gesture_detectors) >= 1

    # Click binding site and verify context map card appears at top
    gesture_detectors[0].on_tap(MagicMock())
    assert len(view.result_list.controls) == 2
    card = view.result_list.controls[0]
    assert isinstance(card, ft.Card)


def test_pcr_view_all_remaining_branches() -> None:
    """Test all remaining branches across PCR view for 100% coverage."""
    from typing import Any, cast

    from amplifyp.amplicon import Amplicon
    from amplifyp.dir_idx import DirIdx
    from amplifyp.dna import DNAType
    from amplifyp.gui.settings import MAX_AMPLICONS_RENDER, GUISettings
    from amplifyp.gui.views.pcr.amplicon_drawing import (
        AmpliconDetailCard,
        DrawnAmplicon,
    )
    from amplifyp.gui.views.pcr.dismissible_detail_card import (
        DismissibleDetailCard,
    )
    from amplifyp.gui.views.pcr.pcr_layout import PCRLayoutSolver
    from amplifyp.gui.views.pcr.primer_drawing import (
        DrawnPrimer,
        ReplicationContextCard,
        get_template_substring,
    )
    from amplifyp.origin import ReplicationOrigin
    from amplifyp.pcr import PCR
    from amplifyp.settings import ReplicationSettings

    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    settings = GUISettings()
    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev1", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data, settings)

    # 1. Properties
    assert view.diagram_container is not None
    assert view.divider is not None
    # 2. _on_pan_update on diagram_panel
    view.diagram_panel._Control__page = mock_page
    view.diagram_panel.update = MagicMock()
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        view.diagram_panel._on_pan_update(
            MagicMock(local_delta=MagicMock(y=50.0))
        )
    assert view.diagram_panel.diagram_container.height >= 200.0

    # With update raising RuntimeError
    with (
        patch.object(
            view.diagram_panel, "update", side_effect=RuntimeError("Update err")
        ),
        patch.object(ft.Control, "page", new=property(lambda self: mock_page)),
    ):
        view.diagram_panel._on_pan_update(
            MagicMock(local_delta=MagicMock(y=-20.0))
        )

    # 3. open_all_cards when _cached_pcr is None vs present
    view.open_all_cards()  # None -> early return

    view.run_pcr()
    view.open_all_cards()
    assert len(view.result_list.controls) > 0

    # 4. _clear_all_cards
    view._clear_all_cards(MagicMock())
    assert len(view.result_list.controls) == 0

    # 5. Tick interval scaling for length > 10000 and length > 5000
    view.diagram_panel._draw_template_baseline(
        v_target=100.0,
        h_margin=40.0,
        c_width=600.0,
        t_width=500.0,
        target_length=15000,
    )
    view.diagram_panel._draw_template_baseline(
        v_target=100.0,
        h_margin=40.0,
        c_width=600.0,
        t_width=500.0,
        target_length=6000,
    )

    # 6. DrawnAmplicon q_score thresholds (<700, <1500, <4000, >=4000)

    p1 = Primer("A" * 20, name="P1")
    p2 = Primer("T" * 20, name="P2")
    canvas = ft.canvas.Canvas()
    stack = ft.Stack()

    for q in [500.0, 1000.0, 3000.0, 5000.0]:
        amp_mock = MagicMock(spec=Amplicon)
        amp_mock.q_score = q
        amp_mock.circular = False
        amp_mock.start = DirIdx(direction=DNADirection.FWD, index=100)
        amp_mock.end = DirIdx(direction=DNADirection.REV, index=300)
        amp_mock.product = DNA("A" * 200)
        da = DrawnAmplicon(
            amp=amp_mock,
            idx=0,
            target_length=1000,
            t_width=500.0,
            h_margin=40.0,
            v_target=100.0,
            c_width=600.0,
            settings=settings,
            on_click=MagicMock(),
        )
        da.draw(canvas, stack)

    # 7. DrawnAmplicon circular: w_right >= w_left and w_right < w_left
    amp_circ1 = MagicMock(spec=Amplicon)
    amp_circ1.q_score = 100.0
    amp_circ1.circular = True
    amp_circ1.start = DirIdx(direction=DNADirection.FWD, index=800)
    amp_circ1.end = DirIdx(direction=DNADirection.REV, index=100)
    amp_circ1.product = DNA("A" * 300)

    da_circ1 = DrawnAmplicon(
        amp=amp_circ1,
        idx=0,
        target_length=1000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        c_width=600.0,
        settings=settings,
        on_click=MagicMock(),
        v_frag_start=None,
    )
    da_circ1.draw(canvas, stack)

    amp_circ2 = MagicMock(spec=Amplicon)
    amp_circ2.q_score = 100.0
    amp_circ2.circular = True
    amp_circ2.start = DirIdx(direction=DNADirection.FWD, index=900)
    amp_circ2.end = DirIdx(direction=DNADirection.REV, index=400)
    amp_circ2.product = DNA("A" * 300)

    da_circ2 = DrawnAmplicon(
        amp=amp_circ2,
        idx=0,
        target_length=1000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        c_width=600.0,
        settings=settings,
        on_click=MagicMock(),
    )
    da_circ2.draw(canvas, stack)

    # 8. AmpliconDetailCard short sequence (len(full_seq) < fwd_len + rev_len)
    amp_short = MagicMock(spec=Amplicon)
    amp_short.fwd_origin = p1
    amp_short.rev_origin = p2
    amp_short.start = DirIdx(direction=DNADirection.FWD, index=0)
    amp_short.end = DirIdx(direction=DNADirection.REV, index=10)
    amp_short.product = DNA("A" * 10)
    amp_short.q_score = 100.0

    card_short = AmpliconDetailCard(
        amp_short, settings, dismiss_callback=MagicMock()
    )
    assert card_short is not None

    # 9. DismissibleDetailCard close button dismiss callback
    dismissed = False

    def on_dismiss(c: ft.Card) -> None:
        nonlocal dismissed
        dismissed = True

    card_dismiss = DismissibleDetailCard(
        card_id="test_card",
        title="Test Card",
        settings=settings,
        dismiss_callback=on_dismiss,
        body_controls=[ft.Container()],
    )
    # Trigger close button
    close_btn = card_dismiss.content.content.controls[0].controls[-1]
    close_btn.on_click(MagicMock())
    assert dismissed is True

    # 10. PrimerDrawing: REV primer with bent leader line
    # and get_template_substring on circular template
    drawn_rev = DrawnPrimer(
        name="rev1",
        index=100,
        conf=MagicMock(direction=DNADirection.REV),
        var=DirIdx(direction=DNADirection.REV, index=100),
        S=10.0,
        target_length=1000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        settings=settings,
        on_click=MagicMock(),
        x_shifted=150.0,  # different from x_pos to trigger bent leader line
    )
    drawn_rev.draw(canvas, stack)

    tmpl_circ = DNA("ATGCATGC", dna_type=DNAType.CIRCULAR)
    region_circ = get_template_substring(tmpl_circ, 6, 6)
    assert len(region_circ) == 6

    # 11. ReplicationContextCard with amplify4_compatibility_mode
    mock_conf = MagicMock(spec=Repliconf)
    mock_conf.primer = Primer("ATGC", name="P1")
    mock_conf.direction = DNADirection.FWD
    mock_conf.template = DNA("ATGCATGCATGCATGC")
    mock_origin = MagicMock(spec=ReplicationOrigin)
    mock_origin.settings = ReplicationSettings(amplify4_compatibility_mode=True)
    mock_origin.primability = 0.95
    mock_origin.stability = 0.90
    mock_origin.quality = 1.0
    mock_conf.origin.return_value = mock_origin

    ctx_card = ReplicationContextCard(
        primer_name="P1",
        padded_idx=0,
        conf=mock_conf,
        var=DirIdx(direction=DNADirection.FWD, index=0),
        settings=settings,
        dismiss_callback=MagicMock(),
    )
    assert ctx_card is not None
    spans_text = "".join(
        span.text
        for span in ctx_card.content.content.controls[1]
        .content.controls[0]
        .spans
    )
    assert "5'-" in spans_text
    assert "-3'" in spans_text

    space_settings = GUISettings()
    space_settings["sequence_separator"] = "Space (' ')"
    ctx_card_space = ReplicationContextCard(
        primer_name="P1",
        padded_idx=0,
        conf=mock_conf,
        var=DirIdx(direction=DNADirection.FWD, index=0),
        settings=space_settings,
        dismiss_callback=MagicMock(),
    )
    spans_text_space = "".join(
        span.text
        for span in ctx_card_space.content.content.controls[1]
        .content.controls[0]
        .spans
    )
    assert "5' " in spans_text_space
    assert " 3'" in spans_text_space

    # 12. PCRLayoutSolver: clusters hitting boundaries and missing confs
    # Left boundary cluster
    b_left = {
        (0, "p1"): (
            "p1",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=0),
        ),
        (1, "p2"): (
            "p2",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=1),
        ),
        (2, "p3"): (
            "p3",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=2),
        ),
    }
    shifted_left = PCRLayoutSolver.calculate_shifted_x(
        cast(Any, b_left), 1000, 500.0, 40.0
    )
    assert len(shifted_left) == 3

    # Right boundary cluster & multi-cluster
    b_multi = {
        (0, "p1"): (
            "p1",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=0),
        ),
        (998, "p2"): (
            "p2",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=998),
        ),
        (999, "p3"): (
            "p3",
            10.0,
            mock_conf,
            DirIdx(direction=DNADirection.FWD, index=999),
        ),
    }
    shifted_multi = PCRLayoutSolver.calculate_shifted_x(
        cast(Any, b_multi), 1000, 500.0, 40.0
    )
    assert len(shifted_multi) == 3

    # PCRLayoutSolver.collect_primer_bindings with mismatched repliconfs
    # and search not done
    unsearched_conf = MagicMock(spec=Repliconf)
    unsearched_conf.searched = False
    unsearched_conf.primer = Primer("ATGC", name="P_un")
    unsearched_conf.origin_db = MagicMock(
        fwd=[DirIdx(direction=DNADirection.FWD, index=0)],
        rev=[DirIdx(direction=DNADirection.REV, index=10)],
    )
    unsearched_conf.origin.return_value = mock_origin

    pcr_mismatch = MagicMock(spec=PCR)

    pcr_mismatch.amplicons = [amp_short]
    pcr_mismatch.amplicon_generator = MagicMock(repliconfs=[])
    fwd_b, rev_b = PCRLayoutSolver.collect_primer_bindings(
        pcr_mismatch, [amp_short]
    )
    assert len(fwd_b) == 0
    assert len(rev_b) == 0

    pcr_no_amp = MagicMock(spec=PCR)

    pcr_no_amp.amplicons = []
    pcr_no_amp.amplicon_generator = MagicMock(repliconfs=[unsearched_conf])
    fwd_b2, rev_b2 = PCRLayoutSolver.collect_primer_bindings(pcr_no_amp, [])
    assert unsearched_conf.search.called
    assert len(fwd_b2) == 1
    assert len(rev_b2) == 1

    # 13. MAX_AMPLICONS_RENDER warning in run_pcr and render_diagram
    many_amps = [amp_short] * (MAX_AMPLICONS_RENDER + 10)
    pcr_many = MagicMock(spec=PCR)
    pcr_many.amplicons = many_amps
    pcr_many.template = DNA("A" * 1000)
    pcr_many.amplicon_generator = MagicMock(repliconfs=[])

    view._cached_pcr = None
    view._cached_state_key = None
    with patch.object(view, "_execute_pcr_simulation", return_value=pcr_many):
        view.run_pcr()
        assert any(
            isinstance(c, ft.Container)
            and "amplicons" in getattr(getattr(c, "content", None), "value", "")
            for c in view.result_list.controls
        )

    # 14. run_pcr exception handling
    view._cached_pcr = None
    view._cached_state_key = None
    with (
        patch.object(
            view, "_execute_pcr_simulation", side_effect=RuntimeError("PCR err")
        ),
        patch("amplifyp.gui.utils.gui_helpers.show_error_dialog") as mock_err,
    ):
        success = view.run_pcr()
        assert success is False
        mock_err.assert_called_once()

    # 15. Duplicate context map & amplicon cards (moves to top) and dismissal
    view._show_context_map(
        "P1", 0, mock_conf, DirIdx(direction=DNADirection.FWD, index=0)
    )
    assert len(view.result_list.controls) > 0
    first_ctrl = view.result_list.controls[0]
    view._show_context_map(
        "P1", 0, mock_conf, DirIdx(direction=DNADirection.FWD, index=0)
    )
    assert view.result_list.controls[0] == first_ctrl

    # Dismiss context card
    if hasattr(first_ctrl, "_card_id"):
        # invoke its dismiss callback
        close_btn = first_ctrl.content.content.controls[0].controls[-1]
        close_btn.on_click(MagicMock())

    view._show_amplicon_dialog(amp_short)
    amp_ctrl = view.result_list.controls[0]
    view._show_amplicon_dialog(amp_short)
    assert view.result_list.controls[0] == amp_ctrl

    # Dismiss amplicon card
    close_btn2 = amp_ctrl.content.content.controls[0].controls[-1]
    close_btn2.on_click(MagicMock())

    # _draw_amplicons with amplicons=None
    view.diagram_panel._draw_amplicons(
        pcr_many,
        target_length=15000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        c_width=600.0,
        amplicons=None,
    )


def test_drawn_primer_direction_arrows() -> None:
    """Test that forward primers point right and reverse primers point left."""
    import flet.canvas as cv

    from amplifyp.dir_idx import DirIdx
    from amplifyp.gui.settings import GUISettings
    from amplifyp.gui.views.pcr.primer_drawing import DrawnPrimer

    settings = GUISettings()

    # Forward primer test
    canvas_fwd = cv.Canvas(shapes=[])
    stack_fwd = ft.Stack()
    drawn_fwd = DrawnPrimer(
        name="fwd1",
        index=100,
        conf=MagicMock(direction=DNADirection.FWD),
        var=DirIdx(direction=DNADirection.FWD, index=100),
        S=10.0,
        target_length=1000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        settings=settings,
        on_click=MagicMock(),
    )
    drawn_fwd.draw(canvas_fwd, stack_fwd)

    assert len(canvas_fwd.shapes) == 2
    fwd_line = canvas_fwd.shapes[0]
    fwd_triangle = canvas_fwd.shapes[1]

    # Leader line terminates at v_target - 25 - S / 2
    assert fwd_line.elements[-1].y == 100.0 - 25 - 5.0

    # Triangle: tip points right (tip.x > base.x)
    tip_elem = fwd_triangle.elements[0]
    base_elem1 = fwd_triangle.elements[1]
    base_elem2 = fwd_triangle.elements[2]
    assert tip_elem.x > base_elem1.x
    assert tip_elem.x == drawn_fwd.x_pos + 5.0
    assert base_elem1.x == drawn_fwd.x_pos - 5.0
    assert base_elem2.x == drawn_fwd.x_pos - 5.0

    # Reverse primer test
    canvas_rev = cv.Canvas(shapes=[])
    stack_rev = ft.Stack()
    drawn_rev = DrawnPrimer(
        name="rev1",
        index=200,
        conf=MagicMock(direction=DNADirection.REV),
        var=DirIdx(direction=DNADirection.REV, index=200),
        S=10.0,
        target_length=1000,
        t_width=500.0,
        h_margin=40.0,
        v_target=100.0,
        settings=settings,
        on_click=MagicMock(),
    )
    drawn_rev.draw(canvas_rev, stack_rev)

    assert len(canvas_rev.shapes) == 2
    rev_line = canvas_rev.shapes[0]
    rev_triangle = canvas_rev.shapes[1]

    # Leader line terminates at v_target + 25 + S / 2
    assert rev_line.elements[-1].y == 100.0 + 25 + 5.0

    # Triangle: tip points left (tip.x < base.x)
    rev_tip_elem = rev_triangle.elements[0]
    rev_base_elem1 = rev_triangle.elements[1]
    rev_base_elem2 = rev_triangle.elements[2]
    assert rev_tip_elem.x < rev_base_elem1.x
    assert rev_tip_elem.x == drawn_rev.x_pos - 5.0
    assert rev_base_elem1.x == drawn_rev.x_pos + 5.0
    assert rev_base_elem2.x == drawn_rev.x_pos + 5.0


def test_sort_amplicons() -> None:
    """Test amplicon sorting under different ranking modes."""
    from amplifyp.amplicon import Amplicon
    from amplifyp.dir_idx import DirIdx
    from amplifyp.dna import DNA, DNADirection, Primer
    from amplifyp.gui.views.pcr.pcr_diagram_panel import sort_amplicons

    f_primer = Primer("AAAA", name="F")
    r_primer = Primer("TTTT", name="R")

    # amp1: start=10, len=90, q=500
    amp1 = Amplicon(
        product=DNA("A" * 90),
        fwd_origin=f_primer,
        rev_origin=r_primer,
        start=DirIdx(direction=DNADirection.FWD, index=10),
        end=DirIdx(direction=DNADirection.REV, index=100),
        q_score=500.0,
        circular=False,
    )
    # amp2: start=10, len=50, q=300
    amp2 = Amplicon(
        product=DNA("A" * 50),
        fwd_origin=f_primer,
        rev_origin=r_primer,
        start=DirIdx(direction=DNADirection.FWD, index=10),
        end=DirIdx(direction=DNADirection.REV, index=60),
        q_score=300.0,
        circular=False,
    )
    # amp3: start=5, len=75, q=800
    amp3 = Amplicon(
        product=DNA("A" * 75),
        fwd_origin=f_primer,
        rev_origin=r_primer,
        start=DirIdx(direction=DNADirection.FWD, index=5),
        end=DirIdx(direction=DNADirection.REV, index=80),
        q_score=800.0,
        circular=False,
    )

    amps = [amp1, amp2, amp3]

    # 1. Position, then length: amp3 (pos 5), then amp2 (pos 10, len 50),
    # then amp1 (pos 10, len 90)
    res_pos_len = sort_amplicons(amps, "Position, then length")
    assert res_pos_len == [amp3, amp2, amp1]

    # 2. Quality score: amp2 (q 300), then amp1 (q 500), then amp3 (q 800)
    res_quality = sort_amplicons(amps, "Quality score")
    assert res_quality == [amp2, amp1, amp3]

    # 3. Position, then quality
    res_pos_q = sort_amplicons(amps, "Position, then quality")
    assert res_pos_q == [amp3, amp2, amp1]


def test_pcr_view_primer_label_orientation() -> None:
    """Test horizontal (compact) and vertical primer label orientations."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    # Long primer names to verify vertical vs horizontal compactness
    input_data.primers = [
        {
            "name": "very_long_forward_primer_name",
            "seq": "TTCCACTGCGAATCATTAAA",
            "active": True,
        },
        {
            "name": "very_long_reverse_primer_name",
            "seq": "GTGCAGCCTTGGTGAAAACT",
            "active": True,
        },
    ]

    view = PCRView(mock_page, input_data)
    assert view.settings["pcr_primer_label_orientation"] == "Horizontal"
    view.run_pcr()

    horiz_height = view.diagram_panel.diagram_canvas.height
    horiz_labels = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Text)
        and ctrl.content.value
        in ("very_long_forward_primer_name", "very_long_reverse_primer_name")
    ]
    assert len(horiz_labels) == 2
    assert all(lbl.rotate is None for lbl in horiz_labels)
    fwd_lbl = next(
        lbl
        for lbl in horiz_labels
        if lbl.content.value == "very_long_forward_primer_name"
    )
    rev_lbl = next(
        lbl
        for lbl in horiz_labels
        if lbl.content.value == "very_long_reverse_primer_name"
    )
    # Forward primer is at start (left ~ 20), so label is to its right
    assert fwd_lbl.left > 20.0
    # Reverse primer is at end (right ~ 780), so label is to its left
    assert rev_lbl.left + rev_lbl.content.width < 780.0
    assert rev_lbl.content.text_align == ft.TextAlign.RIGHT

    # Switch to 45° orientation
    view.settings["pcr_primer_label_orientation"] = "45°"
    view.run_pcr()

    diag_height = view.diagram_panel.diagram_canvas.height
    diag_labels = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Text)
        and ctrl.content.value
        in ("very_long_forward_primer_name", "very_long_reverse_primer_name")
    ]
    assert len(diag_labels) == 2
    diag_fwd = next(
        lbl
        for lbl in diag_labels
        if lbl.content.value == "very_long_forward_primer_name"
    )
    diag_rev = next(
        lbl
        for lbl in diag_labels
        if lbl.content.value == "very_long_reverse_primer_name"
    )
    assert (
        diag_fwd.rotate is not None
        and diag_fwd.rotate.angle == pytest.approx(-math.pi / 4)
        and diag_fwd.rotate.alignment == ft.Alignment(-1, 0)
    )
    assert (
        diag_rev.rotate is not None
        and diag_rev.rotate.angle == pytest.approx(-math.pi / 4)
        and diag_rev.rotate.alignment == ft.Alignment(1, 0)
    )
    assert diag_rev.left + diag_rev.content.width < 780.0
    assert diag_rev.content.text_align == ft.TextAlign.RIGHT

    # Switch to Vertical orientation
    view.settings["pcr_primer_label_orientation"] = "Vertical"
    view.run_pcr()

    vert_height = view.diagram_panel.diagram_canvas.height
    vert_labels = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
        and isinstance(ctrl.content, ft.Text)
        and ctrl.content.value
        in ("very_long_forward_primer_name", "very_long_reverse_primer_name")
    ]
    assert len(vert_labels) == 2
    assert all(lbl.rotate is not None for lbl in vert_labels)
    assert horiz_height < diag_height < vert_height


def test_horizontal_label_layout_row_stacking_and_clamping() -> None:
    """Test horizontal label layout stacks overlapping labels into rows."""
    from amplifyp.gui.views.pcr.pcr_layout import PCRLayoutSolver

    mock_conf = MagicMock()
    mock_var = MagicMock()
    bindings: dict[tuple[int, str], tuple[str, float, Any, Any]] = {
        (0, "long_primer_alpha"): (
            "long_primer_alpha",
            16.0,
            mock_conf,
            mock_var,
        ),
        (5, "long_primer_beta"): (
            "long_primer_beta",
            16.0,
            mock_conf,
            mock_var,
        ),
        (500, "edge_primer"): ("edge_primer", 16.0, mock_conf, mock_var),
    }
    shifted_x = {
        (0, "long_primer_alpha"): 20.0,
        (5, "long_primer_beta"): 44.0,
        (500, "edge_primer"): 780.0,
    }

    layout, num_rows = PCRLayoutSolver.calculate_horizontal_label_layout(
        bindings, shifted_x, c_width=800.0, font_size=13, is_fwd=True
    )

    assert num_rows == 2
    left_a, row_a, _width_a = layout[(0, "long_primer_alpha")]
    _left_b, row_b, _width_b = layout[(5, "long_primer_beta")]
    left_edge, row_edge, width_edge = layout[(500, "edge_primer")]

    assert row_a == 1
    assert row_b == 0
    assert row_edge == 0
    assert left_a >= 20.0 + 8.0 + 4.0
    assert left_edge + width_edge <= 800.0 - 4.0 + 1e-6

    # Reverse primers: left-hand primer goes in row 0 (closer to baseline),
    # right-hand primer goes in row 1 (further from baseline) so connector
    # lines do not cross labels.
    rev_bindings: dict[tuple[int, str], tuple[str, float, Any, Any]] = {
        (200, "rev_alpha"): ("rev_alpha", 16.0, mock_conf, mock_var),
        (224, "rev_beta"): ("rev_beta", 16.0, mock_conf, mock_var),
    }
    rev_shifted = {
        (200, "rev_alpha"): 200.0,
        (224, "rev_beta"): 224.0,
    }
    rev_layout, rev_rows = PCRLayoutSolver.calculate_horizontal_label_layout(
        rev_bindings, rev_shifted, c_width=800.0, font_size=13, is_fwd=False
    )
    assert rev_rows == 2
    _, rev_row_alpha, _ = rev_layout[(200, "rev_alpha")]
    _, rev_row_beta, _ = rev_layout[(224, "rev_beta")]
    assert rev_row_alpha == 0
    assert rev_row_beta == 1


def test_drawn_primer_triangle_and_label_same_line_when_stacked() -> None:
    """Test that stacked horizontal primers keep triangle and label aligned."""
    from amplifyp.gui.settings import GUISettings
    from amplifyp.gui.views.pcr.primer_drawing import DrawnPrimer

    settings = GUISettings()
    settings["pcr_primer_label_orientation"] = "Horizontal"
    mock_conf = MagicMock()
    mock_var = MagicMock()
    mock_var.direction = DNADirection.FWD

    canvas = MagicMock()
    canvas.shapes = []
    stack = MagicMock()
    stack.controls = []

    drawn = DrawnPrimer(
        name="stacked_fwd",
        index=0,
        conf=mock_conf,
        var=mock_var,
        S=16.0,
        target_length=100,
        t_width=760.0,
        h_margin=20.0,
        v_target=118.0,
        settings=settings,
        on_click=lambda: None,
        label_left=32.0,
        label_row=1,
    )
    drawn.draw(canvas, stack)

    # Triangle is the second shape (index 1); its tip y is at tri_mid:
    # v_target (118) - 25 - row_offset (18) - S/2 (8) = 67.0
    tri_path = canvas.shapes[1]
    tri_mid_y = tri_path.elements[0].y
    assert tri_mid_y == 67.0

    # Label is the first stack control; vertically centred on tri_mid_y
    label_ctrl = stack.controls[0]
    assert label_ctrl.top == 67.0 - (13 + 3.0) / 2.0


def test_pcr_divider_autoresize_viewport_fractions() -> None:
    """Test autoresize fits PCR map up to configured viewport fractions."""
    from amplifyp.gui.settings import GUISettings

    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800
    # 1000 - 100 (nav bar) = 900 available viewport height
    mock_page.height = 1000

    settings = GUISettings()
    assert settings["pcr_divider_behaviour"] == "autoresize_3_4"

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev1", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data, settings)
    view.run_pcr()

    # Single amplicon canvas is smaller than 3/4 of 900 (675),
    # clamped to >= 150.0 and <= 300.0.
    assert view.diagram_container.height is not None
    assert 150.0 <= view.diagram_container.height <= 300.0

    # Simulate a tall canvas (many amplicons): capped at 3/4 of 900 (675.0).
    view.diagram_panel.diagram_canvas.height = 1000.0
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 675.0

    # Viewport grows to 1300 (1200 available) -> 3/4 cap becomes 900.0.
    mock_page.height = 1300
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 900.0

    # 1/4 viewport cap on 1200 available -> 300.0.
    settings["pcr_divider_behaviour"] = "autoresize_1_4"
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 300.0

    # 1/3 viewport cap on 1200 available -> 400.0.
    settings["pcr_divider_behaviour"] = "autoresize"
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 400.0

    # 1/2 viewport cap on 1200 available -> 600.0.
    settings["pcr_divider_behaviour"] = "autoresize_1_2"
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 600.0

    # 2/3 viewport cap on 1200 available -> 800.0.
    settings["pcr_divider_behaviour"] = "autoresize_2_3"
    view.diagram_panel._apply_container_height()
    assert view.diagram_container.height == 800.0


def test_pcr_divider_save_proportionality() -> None:
    """Test saving and restoring PCR map height proportionality."""
    from amplifyp.gui.settings import GUISettings

    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800
    # 1100 - 100 (nav bar) = 1000 available viewport height
    mock_page.height = 1100

    settings = GUISettings()
    settings["pcr_divider_behaviour"] = "save_proportionality"

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.template_circular = False
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev1", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = PCRView(mock_page, input_data, settings)
    view.run_pcr()

    # Default ratio 0.375 on 1000px available viewport -> 375.0px
    assert view.diagram_container.height == 375.0

    # Drag divider down by 125px -> 500px -> ratio 0.5
    view.diagram_panel.update = MagicMock()
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        view.diagram_panel._on_pan_update(
            MagicMock(local_delta=MagicMock(y=125.0))
        )
    assert view.diagram_container.height == 500.0
    assert settings["pcr_map_height_ratio"] == 0.5

    with patch.object(settings, "save_to_local") as mock_save:
        view.diagram_panel._on_pan_end(MagicMock(spec=ft.DragEndEvent))
        mock_save.assert_called_once_with(
            mock_page, keys=("pcr_map_height_ratio",)
        )

    # Resize to 900px (800 available) preserves 0.5 proportionality -> 400px
    mock_page.height = 900
    view._handle_resize(MagicMock())
    assert view.diagram_container.height == 400.0

    # New view in a new session with the same persisted settings restores 0.5
    new_view = PCRView(mock_page, input_data, settings)
    new_view.run_pcr()
    assert new_view.diagram_container.height == 400.0

    # Tiny viewport height (<= 100px -> 0 available) does not raise
    mock_page.height = 80
    with patch.object(ft.Control, "page", new=property(lambda self: mock_page)):
        new_view.diagram_panel._on_pan_update(
            MagicMock(local_delta=MagicMock(y=10.0))
        )


def test_pcr_zeroed_weights_and_visual_cache() -> None:
    """Test zeroed weights do not crash PCR and visual settings hit cache."""
    from amplifyp.origin import ReplicationOrigin
    from amplifyp.settings import LengthWiseWeightTbl, ReplicationSettings

    rep = ReplicationSettings(
        run_weights=LengthWiseWeightTbl(default_weight=0.0, overrides={}),
        match_weight=LengthWiseWeightTbl(default_weight=0.0, overrides={}),
    )
    origin = ReplicationOrigin(target="ATGC", primer="GCAT", settings=rep)
    assert origin.primability == 0.0
    assert origin.stability == 0.0

    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800
    mock_page.height = 900
    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev1", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]
    view = PCRView(mock_page, input_data)
    view.run_pcr()

    # Changing a purely visual setting should reuse cached PCR simulation
    view.settings["pcr_map_height_ratio"] = 0.5
    view.settings["pcr_primer_label_orientation"] = "45°"
    with patch.object(
        view,
        "_execute_pcr_simulation",
        side_effect=AssertionError("Cache miss"),
    ):
        assert view.run_pcr(keep_cards=True) is True


def test_pcr_view_keep_cards_does_not_duplicate_banners() -> None:
    """Test run_pcr(keep_cards=True) preserves cards without cloning banners."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800
    mock_page.height = 900

    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
    ]

    view = PCRView(mock_page, input_data)
    view.run_pcr()
    assert len(view.result_list.controls) == 1
    assert view.cards_header.visible is False

    # Resizing / re-running with keep_cards=True does not duplicate banner
    view.run_pcr(keep_cards=True)
    view.run_pcr(keep_cards=True)
    assert len(view.result_list.controls) == 1
    assert isinstance(view.result_list.controls[0], ft.Text)
    assert view.cards_header.visible is False

    # Opening a binding-site card then resizing keeps card and single banner
    detectors = [
        c
        for c in view.diagram_stack.controls
        if isinstance(c, ft.GestureDetector)
    ]
    detectors[0].on_tap(MagicMock())
    assert len(view.result_list.controls) == 2
    view.run_pcr(keep_cards=True)
    cards = [c for c in view.result_list.controls if isinstance(c, ft.Card)]
    texts = [c for c in view.result_list.controls if isinstance(c, ft.Text)]
    assert len(cards) == 1
    assert len(texts) == 1
    assert isinstance(view.result_list.controls[0], ft.Card)
    assert view.cards_header.visible is True


def test_pcr_edge_primer_labels_45_and_horizontal() -> None:
    """Test edge primer labels stay anchored at 45° and flip horizontally."""
    from amplifyp.dir_idx import DirIdx
    from amplifyp.gui.settings import GUISettings
    from amplifyp.gui.views.pcr.pcr_layout import PCRLayoutSolver
    from amplifyp.gui.views.pcr.primer_drawing import DrawnPrimer

    settings = GUISettings()
    page_width = 680.0
    target_len = 1000
    long_fwd_name = "VeryLongForwardEdgePrimerName"
    long_rev_name = "VeryLongReverseEdgePrimerName"

    fwd_conf: Repliconf = MagicMock(spec=Repliconf, direction=DNADirection.FWD)
    rev_conf: Repliconf = MagicMock(spec=Repliconf, direction=DNADirection.REV)
    fwd_var = DirIdx(DNADirection.FWD, 998)
    rev_var = DirIdx(DNADirection.REV, 2)

    fwd_dict: dict[tuple[int, str], tuple[str, float, Repliconf, DirIdx]] = {
        (998, long_fwd_name): (long_fwd_name, 12.0, fwd_conf, fwd_var)
    }
    rev_dict: dict[tuple[int, str], tuple[str, float, Repliconf, DirIdx]] = {
        (2, long_rev_name): (long_rev_name, 12.0, rev_conf, rev_var)
    }

    # 1. 45° mode expands h_margin so diagonal labels stay anchored at x_pos
    v_target, h_margin, c_width, t_width, _ = (
        PCRLayoutSolver.calculate_canvas_dimensions(
            target_length=target_len,
            num_amplicons=0,
            fwd_bindings=fwd_dict,
            rev_bindings=rev_dict,
            page_width=page_width,
            label_orientation="45°",
            font_size_primer=13,
        )
    )
    assert h_margin > 20.0

    settings["pcr_primer_label_orientation"] = "45°"
    canvas = ft.canvas.Canvas()
    stack = ft.Stack()
    dp_fwd = DrawnPrimer(
        name=long_fwd_name,
        index=998,
        conf=fwd_conf,
        var=fwd_var,
        S=12.0,
        target_length=target_len,
        t_width=t_width,
        h_margin=h_margin,
        v_target=v_target,
        settings=settings,
        on_click=MagicMock(),
    )
    dp_fwd.draw(canvas, stack)
    label_ctrl = next(
        c
        for c in stack.controls
        if isinstance(c, ft.GestureDetector) and c.rotate is not None
    )
    assert label_ctrl.left == dp_fwd.x_pos

    # 2. Horizontal mode flips edge-constrained labels to non-overlapping side
    fwd_x = 20.0 + (998 / target_len) * (c_width - 40.0)
    fwd_map, _ = PCRLayoutSolver.calculate_horizontal_label_layout(
        bindings=fwd_dict,
        shifted_x={(998, long_fwd_name): fwd_x},
        c_width=c_width,
        font_size=13,
        is_fwd=True,
    )
    fwd_lbl_left, _, fwd_w = fwd_map[(998, long_fwd_name)]
    assert fwd_lbl_left + fwd_w <= fwd_x - 6.0 + 1e-6

    rev_x = 20.0 + (2 / target_len) * (c_width - 40.0)
    rev_map, _ = PCRLayoutSolver.calculate_horizontal_label_layout(
        bindings=rev_dict,
        shifted_x={(2, long_rev_name): rev_x},
        c_width=c_width,
        font_size=13,
        is_fwd=False,
    )
    rev_lbl_left, _, _ = rev_map[(2, long_rev_name)]
    assert rev_lbl_left >= rev_x + 6.0 - 1e-6

    # 3. Flipped horizontal labels align toward their triangle (RIGHT when on
    # left side of x_render, left-aligned when on right side)
    settings["pcr_primer_label_orientation"] = "Horizontal"
    stack_h = ft.Stack()
    DrawnPrimer(
        name=long_fwd_name,
        index=998,
        conf=fwd_conf,
        var=fwd_var,
        S=12.0,
        target_length=target_len,
        t_width=c_width - 40.0,
        h_margin=20.0,
        v_target=100.0,
        settings=settings,
        on_click=MagicMock(),
        label_left=fwd_lbl_left,
        label_width=fwd_w,
    ).draw(ft.canvas.Canvas(), stack_h)
    fwd_text = next(
        c.content
        for c in stack_h.controls
        if isinstance(c, ft.GestureDetector) and isinstance(c.content, ft.Text)
    )
    assert fwd_text.text_align == ft.TextAlign.RIGHT
    assert fwd_text.width == fwd_w

    stack_rev_h = ft.Stack()
    DrawnPrimer(
        name=long_rev_name,
        index=2,
        conf=rev_conf,
        var=rev_var,
        S=12.0,
        target_length=target_len,
        t_width=c_width - 40.0,
        h_margin=20.0,
        v_target=100.0,
        settings=settings,
        on_click=MagicMock(),
        label_left=rev_lbl_left,
        label_width=rev_map[(2, long_rev_name)][2],
    ).draw(ft.canvas.Canvas(), stack_rev_h)
    rev_text = next(
        c.content
        for c in stack_rev_h.controls
        if isinstance(c, ft.GestureDetector) and isinstance(c.content, ft.Text)
    )
    assert rev_text.text_align == ft.TextAlign.START
    assert rev_text.width is None

    # 4. Zero-denominator quality and infinite amplicon q_score guards
    from amplifyp.amplicon import Amplicon, AmpliconGenerator
    from amplifyp.gui.views.pcr.amplicon_drawing import AmpliconDetailCard
    from amplifyp.origin import ReplicationOrigin
    from amplifyp.settings import ReplicationSettings

    rep_max = ReplicationSettings(primability_cutoff=1.0, stability_cutoff=1.0)
    origin_max = ReplicationOrigin(
        target="ATGC", primer="GCAT", settings=rep_max
    )
    assert origin_max.quality == 0.0

    mock_tpl = DNA("ATGCATGCATGCATGC")
    amp_gen = AmpliconGenerator(mock_tpl)
    mock_conf_q0 = MagicMock(spec=Repliconf)
    mock_conf_q0.origin.return_value = origin_max
    inf_q = amp_gen.get_amplicon_quality_score(
        mock_conf_q0, mock_conf_q0, fwd_var, rev_var, 100
    )
    assert inf_q == float("inf")
    amp_inf = Amplicon(
        product=DNA("ATGCATGCATGC"),
        fwd_origin=Primer("ATGC", name="F"),
        rev_origin=Primer("GCAT", name="R"),
        start=DirIdx(DNADirection.FWD, 0),
        end=DirIdx(DNADirection.REV, 10),
        q_score=inf_q,
        circular=False,
    )
    assert "very weak" in amp_inf.q_score_report_str(verbose=True)
    inf_card = AmpliconDetailCard(
        amp=amp_inf, settings=settings, dismiss_callback=MagicMock()
    )
    assert inf_card is not None

    # 5. run_pcr(keep_cards=True) preserves cards on resize and rebuilds
    # them only when refresh_cards=True (and discards on cache miss)
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800
    mock_page.height = 900
    input_data = GUIInput()
    input_data.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    input_data.primers = [
        {"name": "fwd1", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev1", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]
    prev_dark = GUIColours.dark_mode
    try:
        GUIColours.dark_mode = False
        pcr_view = PCRView(mock_page, input_data)
        pcr_view.run_pcr()
        pcr_view.open_all_cards()
        detectors = [
            c
            for c in pcr_view.diagram_stack.controls
            if isinstance(c, ft.GestureDetector)
        ]
        detectors[0].on_tap(MagicMock())
        old_cards = list(pcr_view.result_list.controls)

        # Resize path (refresh_cards=False) keeps existing card instances
        pcr_view.run_pcr(keep_cards=True)
        resize_cards = list(pcr_view.result_list.controls)
        assert len(resize_cards) == len(old_cards)
        assert resize_cards[0] is old_cards[0]

        # Theme/settings path (refresh_cards=True) rebuilds open cards
        GUIColours.dark_mode = True
        pcr_view.run_pcr(keep_cards=True, refresh_cards=True)
        new_cards = list(pcr_view.result_list.controls)
        assert len(new_cards) == len(old_cards)
        assert new_cards[0] is not old_cards[0]

        # Replication settings change invalidates cache and discards stale cards
        pcr_view.settings["primability_cutoff"] = 0.99
        pcr_view.run_pcr(keep_cards=True, refresh_cards=True)
        assert not any(
            isinstance(c, ft.Card) for c in pcr_view.result_list.controls
        )

        # Invalid pcr_amplicon_ranking falls back to QUALITY_SCORE when drawing
        pcr_view.settings["primability_cutoff"] = 0.8
        pcr_view.settings._settings["pcr_amplicon_ranking"] = "InvalidRanking"
        assert pcr_view.run_pcr() is True
        pcr_view.open_all_cards()
        assert len(pcr_view.result_list.controls) > 0
    finally:
        GUIColours.dark_mode = prev_dark


def test_relative_q_score_display_and_scaling() -> None:
    """Test Relative Q score scaling, 3-decimal formatting, and tiering."""
    from amplifyp.amplicon import Amplicon
    from amplifyp.gui.settings import GUISettings
    from amplifyp.gui.views.pcr.amplicon_drawing import (
        AmpliconDetailCard,
        DrawnAmplicon,
    )

    settings = GUISettings()
    settings["pcr_relative_q_score"] = True

    amp = MagicMock(spec=Amplicon)
    amp.start.index = 0
    amp.end.index = 50
    amp.product = MagicMock()
    amp.product.seq = "ACGT" * 15
    amp.product.__len__.return_value = 60
    amp.fwd_origin = MagicMock()
    amp.fwd_origin.name = "F1"
    amp.fwd_origin.seq = "ACGT"
    amp.rev_origin = MagicMock()
    amp.rev_origin.name = "R1"
    amp.rev_origin.seq = "TGCA"
    amp.circular = False
    # Raw Q = 2839 (weak), min_q = 2838.0 -> scaled Q = 1.000 (good)
    amp.q_score = 2839.0

    drawn_long = DrawnAmplicon(
        amp=amp,
        idx=0,
        target_length=100,
        t_width=500.0,
        h_margin=20.0,
        v_target=100.0,
        c_width=540.0,
        settings=settings,
        on_click=lambda: None,
        min_q_score=2838.0,
    )
    # Unscaled relative Q = 1.000 sits near the top of [1.0, 8.0]
    assert 7.9 < drawn_long.bar_height <= 8.0

    long_card = AmpliconDetailCard(
        amp=amp,
        settings=settings,
        dismiss_callback=lambda _: None,
        min_q_score=2838.0,
    )
    long_info = long_card.content.content.controls[1]
    long_text = "".join(span.text for span in long_info.spans)
    assert "Relative Q = 1.000 (good amplification" in long_text

    # Raw Q = 1200 (moderate), min_q = 600.0 -> scaled Q = 2.000 (good)
    amp.q_score = 1200.0

    drawn_scaled = DrawnAmplicon(
        amp=amp,
        idx=0,
        target_length=100,
        t_width=500.0,
        h_margin=20.0,
        v_target=100.0,
        c_width=540.0,
        settings=settings,
        on_click=lambda: None,
        min_q_score=600.0,
    )
    # Relative Q = 2.000 is slightly thinner than 1.000, still near 8.0
    assert 7.9 < drawn_scaled.bar_height < drawn_long.bar_height

    card = AmpliconDetailCard(
        amp=amp,
        settings=settings,
        dismiss_callback=lambda _: None,
        min_q_score=600.0,
    )
    info_row = card.content.content.controls[1]
    full_text = "".join(span.text for span in info_row.spans)
    assert "Relative Q = 2.000 (good amplification" in full_text

    # Large relative Q (400.0) scales smoothly down between 1.0 and 8.0
    drawn_okay = DrawnAmplicon(
        amp=amp,
        idx=0,
        target_length=100,
        t_width=500.0,
        h_margin=20.0,
        v_target=100.0,
        c_width=540.0,
        settings=settings,
        on_click=lambda: None,
        min_q_score=3.0,
    )
    assert 1.0 < drawn_okay.bar_height < drawn_scaled.bar_height

    # Toggling setting off and rebuilding restores raw Q score & tier
    settings["pcr_relative_q_score"] = False
    rebuilt = card.rebuild(settings)
    rebuilt_info = rebuilt.content.content.controls[1]
    rebuilt_text = "".join(span.text for span in rebuilt_info.spans)
    assert "Q = 1200 (moderate amplification" in rebuilt_text

    # Infinite Q score remains ∞ even when relative scaling is enabled
    settings["pcr_relative_q_score"] = True
    amp.q_score = float("inf")
    inf_card = AmpliconDetailCard(
        amp=amp,
        settings=settings,
        dismiss_callback=lambda _: None,
        min_q_score=2.0,
    )
    inf_info = inf_card.content.content.controls[1]
    inf_text = "".join(span.text for span in inf_info.spans)
    assert "Relative Q = ∞ (very weak amplification" in inf_text


def test_pcr_view_relative_q_score_end_to_end() -> None:
    """Test PCRView wires min_q_score through diagram and detail cards."""
    from pathlib import Path

    import yaml

    from amplifyp.gui.settings import GUISettings

    fixture_path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "save_states"
        / "lots_of_amplicon.yaml"
    )
    with fixture_path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)

    input_data = GUIInput()
    input_data.from_dict(data["input"])
    settings = GUISettings()
    settings["pcr_relative_q_score"] = True

    mock_page = MagicMock(spec=ft.Page)
    mock_page.width = 800
    mock_page.height = 900

    view = PCRView(mock_page, input_data, settings)
    assert view.run_pcr() is True
    assert view.diagram_panel._min_q_score == pytest.approx(146.0)

    assert view._cached_pcr is not None
    best_amp = min(view._cached_pcr.amplicons, key=lambda a: a.q_score)
    view._show_amplicon_dialog(best_amp)

    card = view.result_list.controls[0]
    info_row = card.content.content.controls[1]
    full_text = "".join(span.text for span in info_row.spans)
    assert "Relative Q = 1.000 (good amplification" in full_text

    settings["pcr_relative_q_score"] = False
    assert view.run_pcr(keep_cards=True, refresh_cards=True) is True
    refreshed_card = view.result_list.controls[0]
    refreshed_info = refreshed_card.content.content.controls[1]
    refreshed_text = "".join(span.text for span in refreshed_info.spans)
    assert "Q = 146 (good amplification" in refreshed_text
