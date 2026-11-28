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

"""Unit tests for scripts/bundle_web_deps.py."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import stat
import sys
import tarfile
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import pytest

_PYODIDE_VERSION = "314.0.7"
_MSGPACK_WHEEL = "msgpack-1.1.2-cp314-cp314-pyemscripten_2026_0_wasm32.whl"


class _FakeDistribution:
    """Minimal Distribution stub for testing bundler validation paths."""

    def __init__(
        self,
        version: str,
        files: list[Path],
        base_dir: Path,
        requires: list[str] | None = None,
    ) -> None:
        self.version = version
        self.files = files
        self._base_dir = base_dir
        self.requires = requires or []

    def locate_file(self, path: Path) -> Path:
        return self._base_dir / path


def _load_bundle_module() -> ModuleType:
    script_path = (
        Path(__file__).resolve().parent.parent
        / "scripts"
        / "bundle_web_deps.py"
    )
    spec = importlib.util.spec_from_file_location(
        "bundle_web_deps", script_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["bundle_web_deps"] = module
    spec.loader.exec_module(module)
    return module


def _create_app_archive(path: Path, requirements: list[str]) -> None:
    with tarfile.open(path, "w:gz") as tar:
        main_bytes = b"print('hello')\n"
        main_info = tarfile.TarInfo(name="main.py")
        main_info.size = len(main_bytes)
        tar.addfile(main_info, io.BytesIO(main_bytes))

        req_bytes = "".join(f"{r}\n" for r in requirements).encode("utf-8")
        req_info = tarfile.TarInfo(name="requirements.txt")
        req_info.size = len(req_bytes)
        tar.addfile(req_info, io.BytesIO(req_bytes))


def _msgpack_lock_entry(wheel_bytes: bytes) -> dict[str, object]:
    return {
        "name": "msgpack",
        "version": "1.1.2",
        "file_name": _MSGPACK_WHEEL,
        "sha256": hashlib.sha256(wheel_bytes).hexdigest(),
        "depends": [],
    }


def _setup_mock_dist(
    tmp_path: Path,
    requirements: list[str],
    *,
    lock_packages: dict[str, dict[str, object]] | None = None,
    cached_wheels: dict[str, bytes] | None = None,
) -> tuple[Path, Path, Path]:
    """Create a mock published dist/ directory and optional wheel cache."""
    dist_dir = tmp_path / "dist"
    pyodide_dir = dist_dir / "pyodide"
    cache_root = tmp_path / "cache"
    pyodide_dir.mkdir(parents=True)

    for filename, payload in (cached_wheels or {}).items():
        cached_file = cache_root / _PYODIDE_VERSION / filename
        cached_file.parent.mkdir(parents=True, exist_ok=True)
        cached_file.write_bytes(payload)

    (pyodide_dir / "package.json").write_text(
        json.dumps({"version": _PYODIDE_VERSION}), encoding="utf-8"
    )
    (pyodide_dir / "pyodide-lock.json").write_text(
        json.dumps(
            {
                "info": {"python": "3.14.2"},
                "packages": lock_packages or {},
            }
        ),
        encoding="utf-8",
    )

    archive_path = dist_dir / "app.tar.gz"
    _create_app_archive(archive_path, requirements)
    return dist_dir, pyodide_dir, cache_root


def test_bundle_web_dependencies_real_project_deps(tmp_path: Path) -> None:
    """Verify bundling resolves flet/pyyaml/repath/six and cached msgpack."""
    mod = _load_bundle_module()
    wheel_bytes = b"fake-msgpack-wasm-wheel-bytes"
    dist_dir, pyodide_dir, cache_root = _setup_mock_dist(
        tmp_path,
        ["flet==1.1.3", "pyyaml==6.0.4"],
        lock_packages={
            "msgpack": _msgpack_lock_entry(wheel_bytes),
            "pyyaml": {
                "name": "pyyaml",
                "version": "6.0.3",
                "file_name": "pyyaml-6.0.3-cp314-wasm32.whl",
                "sha256": "deadbeef",
                "depends": [],
            },
        },
        cached_wheels={_MSGPACK_WHEEL: wheel_bytes},
    )

    archive_path = dist_dir / "app.tar.gz"
    index_path = dist_dir / "index.html"
    index_path.write_text(
        '<script>flet.pyodideUrl="pyodide/pyodide.mjs";</script>\n',
        encoding="utf-8",
    )

    mod.bundle_web_dependencies(dist_dir, cache_root=cache_root)

    # Local pyodideUrl should be patched to a valid relative ES module specifier
    assert 'flet.pyodideUrl="./pyodide/pyodide.mjs";' in index_path.read_text(
        "utf-8"
    )

    # msgpack wheel should be copied into dist/pyodide/
    bundled_msgpack = pyodide_dir / _MSGPACK_WHEEL
    assert bundled_msgpack.is_file()
    assert bundled_msgpack.read_bytes() == wheel_bytes

    # Repacked archive should be world-readable (0o644 on POSIX; 0o666 on NTFS)
    mode = stat.S_IMODE(archive_path.stat().st_mode)
    if sys.platform != "win32":
        assert mode == 0o644
    else:
        assert mode & 0o444 == 0o444

    with tarfile.open(archive_path, "r:gz") as tar:
        names = set(tar.getnames())
        assert "main.py" in names
        assert "__pypackages__/flet/__init__.py" in names
        assert "__pypackages__/repath.py" in names
        assert "__pypackages__/six.py" in names
        assert "__pypackages__/yaml/__init__.py" in names
        # Native binaries and bytecode must never enter __pypackages__/
        assert not any(
            n.endswith((".so", ".pyd", ".pyc", ".pyo")) for n in names
        )
        # No runtime requirements remain since msgpack is auto-loaded
        req_file = tar.extractfile("requirements.txt")
        assert req_file is not None
        assert req_file.read().decode("utf-8").strip() == ""

    # Second invocation must reject already-bundled archive without deleting it
    with pytest.raises(RuntimeError, match="already bundled"):
        mod.bundle_web_dependencies(dist_dir, cache_root=cache_root)

    with tarfile.open(archive_path, "r:gz") as tar:
        names_after = set(tar.getnames())
        assert "__pypackages__/flet/__init__.py" in names_after
        assert "__pypackages__/yaml/__init__.py" in names_after


def test_bundle_web_dependencies_rejects_unverified_native_wrapper_pkg(
    tmp_path: Path,
) -> None:
    """Reject unverified native packages even when Python wrappers exist."""
    mod = _load_bundle_module()
    dist_dir, _pyodide_dir, cache_root = _setup_mock_dist(
        tmp_path, ["fake-pillow==13.0.0"]
    )

    wrapper_file = tmp_path / "PIL" / "Image.py"
    wrapper_file.parent.mkdir(parents=True)
    wrapper_file.write_text("# wrapper requiring _imaging\n", encoding="utf-8")

    fake_dist = _FakeDistribution(
        version="13.0.0",
        files=[Path("PIL/Image.py"), Path("PIL/_imaging.so")],
        base_dir=tmp_path,
    )
    with (
        patch("importlib.metadata.distribution", return_value=fake_dist),
        pytest.raises(RuntimeError, match="no verified pure-Python fallback"),
    ):
        mod.bundle_web_dependencies(dist_dir, cache_root=cache_root)


def test_bundle_web_dependencies_rejects_host_version_mismatch(
    tmp_path: Path,
) -> None:
    """Reject copying host files when installed version violates archive pin."""
    mod = _load_bundle_module()
    dist_dir, _pyodide_dir, cache_root = _setup_mock_dist(
        tmp_path, ["pyyaml==9.9.9"]
    )

    with pytest.raises(RuntimeError, match="does not satisfy required"):
        mod.bundle_web_dependencies(dist_dir, cache_root=cache_root)


def test_bundle_web_dependencies_pyodide_wheel_independent_of_host_version(
    tmp_path: Path,
) -> None:
    """Select compatible Pyodide wheel regardless of host package version."""
    mod = _load_bundle_module()
    wheel_bytes = b"fake-msgpack-wasm-wheel"
    dist_dir, pyodide_dir, cache_root = _setup_mock_dist(
        tmp_path,
        ["msgpack==1.1.2"],
        lock_packages={"msgpack": _msgpack_lock_entry(wheel_bytes)},
        cached_wheels={_MSGPACK_WHEEL: wheel_bytes},
    )

    fake_dist = _FakeDistribution(
        version="1.0.0",
        files=[Path("msgpack/_cmsgpack.so")],
        base_dir=tmp_path,
    )
    with patch("importlib.metadata.distribution", return_value=fake_dist):
        mod.bundle_web_dependencies(dist_dir, cache_root=cache_root)

    assert (pyodide_dir / _MSGPACK_WHEEL).read_bytes() == wheel_bytes


def test_ensure_pyodide_wheel_rejects_sha256_mismatch(tmp_path: Path) -> None:
    """Verify corrupted downloads fail SHA-256 verification."""
    mod = _load_bundle_module()
    pyodide_dir = tmp_path / "pyodide"
    cache_dir = tmp_path / "cache"
    pyodide_dir.mkdir()

    pkg_entry = {
        "name": "msgpack",
        "version": "1.1.2",
        "file_name": "msgpack-1.1.2-cp314-wasm32.whl",
        "sha256": "0" * 64,
    }

    with (
        patch(
            "urllib.request.urlopen",
            return_value=io.BytesIO(b"corrupted-bytes"),
        ),
        pytest.raises(RuntimeError, match="SHA-256 mismatch"),
    ):
        mod._ensure_pyodide_wheel(
            pkg_entry, _PYODIDE_VERSION, pyodide_dir, cache_dir
        )
