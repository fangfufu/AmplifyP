# Copyright (C) 2026 Fufu Fang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

"""Tests for Result View and primer binding site context map popups."""

from unittest.mock import MagicMock

import flet as ft

from amplifyp.gui.state import GUIState
from amplifyp.gui.views.result_view import ResultView


def test_result_view_click_context_map() -> None:
    """Test that clicking a primer binding site shows a context map."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    state = GUIState()
    state.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    state.template_circular = False
    state.primers = [
        {"name": "10290", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = ResultView(mock_page, state)
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

    # Extract diagram_text
    diagram_text = card.content.content.controls[1].controls[0]
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
    assert "↓" in full_text
    assert "10290" in full_text
    assert "5'" in full_text
    assert "3'" in full_text
    assert "|||" in full_text
    assert "Context" in full_text
    assert (
        "--------------------" in full_text
    )  # Upstream region (gaps because index 0 is at boundary)
    assert "TTCCACTGCGAATCATTAAA" in full_text  # Binding region


def test_result_view_click_amplicon() -> None:
    """Test that clicking on an amplicon bar inserts the amplicon details

    card below the map.
    """
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    state = GUIState()
    state.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    state.template_circular = False
    state.primers = [
        {"name": "10290", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = ResultView(mock_page, state)
    view.run_pcr()

    # Verify no cards initially
    assert len(view.result_list.controls) == 0

    # Find the gesture detectors.
    gesture_detectors = [
        ctrl
        for ctrl in view.diagram_stack.controls
        if isinstance(ctrl, ft.GestureDetector)
    ]
    # We should have 3 gesture detectors total
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


def test_result_view_resize_preserves_cards() -> None:
    """Test that resizing the window does not discard open cards."""
    mock_page = MagicMock(spec=ft.Page)
    mock_page.dialog = None
    mock_page.width = 800

    state = GUIState()
    state.template = (
        "TTCCACTGCGAATCATTAAAGTGGGTATCACAAATTTGGGAGTTTTCACCAAGGCTGCAC"
    )
    state.template_circular = False
    state.primers = [
        {"name": "10290", "seq": "TTCCACTGCGAATCATTAAA", "active": True},
        {"name": "rev_primer", "seq": "GTGCAGCCTTGGTGAAAACT", "active": True},
    ]

    view = ResultView(mock_page, state)
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
    view.handle_resize(MagicMock(spec=ft.ControlEvent))

    # Verify that the card is still present
    assert len(view.result_list.controls) == 1
    assert isinstance(view.result_list.controls[0], ft.Card)
