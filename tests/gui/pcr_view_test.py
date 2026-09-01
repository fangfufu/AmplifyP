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
    assert q_text == "Quality: 1"

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
    from unittest.mock import patch

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
