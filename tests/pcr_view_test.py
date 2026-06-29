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

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.user_data import GUIInput
from amplifyp.gui.views.pcr import PCRView


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
    title_text = card.content.content.controls[0].controls[0].value
    assert "Context Map" in title_text
    assert "10290" in title_text

    # Verify Primeability, Stability and Quality display
    stats_row = card.content.content.controls[1]
    assert isinstance(stats_row, ft.Row)
    p_text = stats_row.controls[0].content.value
    s_text = stats_row.controls[1].content.value
    q_text = stats_row.controls[2].content.value
    assert "Primeability: 1.000" in p_text
    assert "Stability: 1.000" in s_text
    assert "Quality: 1.0000" in q_text

    # Extract diagram_text
    diagram_text = card.content.content.controls[2].content.controls[0]
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

    length_text = column.controls[0].controls[0].value
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
    diagram_text = card.content.content.controls[2].content.controls[0]
    assert isinstance(diagram_text, ft.Text)

    # Check spans
    # It should have a span for the comp_line with MUTED_GREY colour
    from amplifyp.gui.colours import GUIColours

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
    from amplifyp.dna import DNA, DNADirection, Primer
    from amplifyp.gui.views.pcr.primer_drawing import format_context_lines
    from amplifyp.repliconf import Repliconf

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
    # So the primer sequence starts at index 47.
    # Verify that primer sequence starts with '3\'-' at offset 44
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
    # Let's check:
    assert lines_bot[1].startswith("Context  " + " " * 15 + "5'-")
