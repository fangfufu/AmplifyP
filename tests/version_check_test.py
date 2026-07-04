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

"""Tests for version checking utilities."""

from unittest.mock import MagicMock, patch

from amplifyp.gui.utils.version_check import (
    fetch_latest_release_version,
    is_newer_version,
    should_check_for_updates,
)


def test_is_newer_version() -> None:
    """Test is_newer_version comparison correctness."""
    assert is_newer_version("v1.4.0", "1.3.2") is True
    assert is_newer_version("1.4.0", "1.3.2") is True
    assert is_newer_version("v1.3.2", "1.3.2") is False
    assert is_newer_version("1.3.1", "1.3.2") is False
    assert is_newer_version("v2.0.0-beta", "1.3.2") is True
    assert is_newer_version("invalid", "1.3.2") is False
    assert is_newer_version("1.3.2", "1.3.2-beta") is True
    assert is_newer_version("1.3.2-beta2", "1.3.2-beta1") is True
    assert is_newer_version("1.3.2-beta", "1.3.2-alpha") is True
    assert is_newer_version("1.3.2-rc1", "1.3.2-beta2") is True
    assert is_newer_version("1.4.1", "1.4") is True
    assert is_newer_version("1.4.1", "1.4.0") is True
    assert is_newer_version("1.4", "1.4.0") is False
    assert is_newer_version("1.4.0", "1.4") is False
    assert is_newer_version("1.3.2-post1", "1.3.2") is True
    assert is_newer_version("1.3.2-post1", "1.3.2-beta") is True


def test_should_check_for_updates() -> None:
    """Test should_check_for_updates rules."""
    # At Startup always checks
    assert should_check_for_updates("At Startup", 1000.0, 1001.0) is True

    # Disabled never checks
    assert should_check_for_updates("Disabled", 1000.0, 5000000.0) is False

    # Once per Day (86400s)
    assert (
        should_check_for_updates("Once per Day", 1000.0, 1000.0 + 86399)
        is False
    )
    assert (
        should_check_for_updates("Once per Day", 1000.0, 1000.0 + 86400) is True
    )

    # Once per Week (604800s)
    assert (
        should_check_for_updates("Once per Week", 1000.0, 1000.0 + 604799)
        is False
    )
    assert (
        should_check_for_updates("Once per Week", 1000.0, 1000.0 + 604800)
        is True
    )

    # Once per Month (2592000s)
    assert (
        should_check_for_updates("Once per Month", 1000.0, 1000.0 + 2591999)
        is False
    )
    assert (
        should_check_for_updates("Once per Month", 1000.0, 1000.0 + 2592000)
        is True
    )


@patch("urllib.request.urlopen")
def test_fetch_latest_release_version_success(
    mock_urlopen: MagicMock,
) -> None:
    """Test fetch_latest_release_version on successful response."""
    mock_response = MagicMock()
    mock_response.read.return_value = b'{"tag_name": "v1.4.2"}'
    mock_response.__enter__.return_value = mock_response
    mock_urlopen.return_value = mock_response

    assert fetch_latest_release_version() == "v1.4.2"


@patch("urllib.request.urlopen")
def test_fetch_latest_release_version_failure(
    mock_urlopen: MagicMock,
) -> None:
    """Test fetch_latest_release_version on network failure."""
    mock_urlopen.side_effect = Exception("Connection Refused")
    assert fetch_latest_release_version() is None
