#!/usr/bin/env python3
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

"""Bundle Python dependencies into a Flet --no-cdn web build for offline use.

Flet's ``flet publish --no-cdn`` bundles CanvasKit, fonts, and the Pyodide core
runtime, but still leaves Python dependencies from ``requirements.txt`` to be
fetched from PyPI / jsDelivr at runtime. This script closes that gap so the
resulting ``dist/`` directory runs completely air-gapped:

1. Resolves the Emscripten dependency graph from ``requirements.txt`` inside
   ``dist/app.tar.gz`` using ``importlib.metadata``.
2. For C-extension packages whose version in ``dist/pyodide/pyodide-lock.json``
   satisfies the required specifiers (such as ``msgpack``), downloads and
   SHA-256-verifies the pre-compiled Pyodide wheel into ``dist/pyodide/``.
3. For pure-Python packages (such as ``flet`` and ``repath``) and packages with
   pure-Python fallbacks whose pinned version is newer than Pyodide's lockfile
   (such as ``pyyaml``), copies their Python files from the active environment
   into ``__pypackages__/`` inside ``dist/app.tar.gz`` (omitting host binaries).
4. Rewrites ``requirements.txt`` inside ``dist/app.tar.gz`` so the browser
   never contacts PyPI or external CDNs at startup.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from collections import deque
from pathlib import Path
from typing import Any

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name
from packaging.version import Version

_CDN_WHEEL_URL = "https://cdn.jsdelivr.net/pyodide/v{version}/full/{filename}"
_NATIVE_SUFFIXES = frozenset({".so", ".pyd", ".dll", ".dylib"})
_SKIPPED_SUFFIXES = frozenset({".pyc", ".pyo", ".pyi"})
# Flet's python-worker.js loads msgpack directly via pyodide_js.loadPackage().
_AUTO_LOADED_PYODIDE_PACKAGES = frozenset({"msgpack"})


def _default_cache_root() -> Path:
    """Return the local Pyodide wheel cache directory."""
    return Path.home() / ".flet" / "cache" / "pyodide"


def _read_pyodide_metadata(
    pyodide_dir: Path,
) -> tuple[str, str, dict[str, dict[str, Any]]]:
    """Read Pyodide version, Python version, and lockfile packages."""
    lock_path = pyodide_dir / "pyodide-lock.json"
    if not lock_path.is_file():
        raise FileNotFoundError(
            f"Expected Pyodide lockfile at {lock_path}. "
            "Run 'flet publish --no-cdn' first."
        )

    lock_data: dict[str, Any] = json.loads(lock_path.read_text("utf-8"))
    info: dict[str, Any] = lock_data.get("info", {})
    python_version = str(info.get("python", "3.14.0"))

    pyodide_version: str | None = None
    pkg_json_path = pyodide_dir / "package.json"
    if pkg_json_path.is_file():
        pkg_data: dict[str, Any] = json.loads(pkg_json_path.read_text("utf-8"))
        raw_ver = pkg_data.get("version")
        if isinstance(raw_ver, str) and raw_ver:
            pyodide_version = raw_ver
    if not pyodide_version:
        raw_info_ver = info.get("version")
        if isinstance(raw_info_ver, str) and raw_info_ver:
            pyodide_version = raw_info_ver
    if not pyodide_version:
        raise RuntimeError(
            f"Could not determine Pyodide version in {pyodide_dir}"
        )

    packages: dict[str, dict[str, Any]] = {}
    for entry in lock_data.get("packages", {}).values():
        if isinstance(entry, dict) and isinstance(entry.get("name"), str):
            packages[canonicalize_name(entry["name"])] = entry

    return pyodide_version, python_version, packages


def _make_emscripten_env(python_version: str) -> dict[str, str]:
    """Build a PEP 508 marker environment for Pyodide / Emscripten."""
    parts = python_version.split(".")
    short_ver = ".".join(parts[:2]) if len(parts) >= 2 else python_version
    return {
        "os_name": "posix",
        "sys_platform": "emscripten",
        "platform_machine": "wasm32",
        "platform_python_implementation": "CPython",
        "platform_system": "Emscripten",
        "python_version": short_ver,
        "python_full_version": python_version,
        "implementation_name": "cpython",
        "implementation_version": python_version,
        "extra": "",
    }


def _read_requirements_from_archive(archive_path: Path) -> list[str]:
    """Extract non-empty requirement lines from requirements.txt in archive."""
    with tarfile.open(archive_path, "r:gz") as tar:
        for member in tar.getmembers():
            if member.name.lstrip("./") == "requirements.txt":
                extracted = tar.extractfile(member)
                if extracted is None:
                    return []
                text = extracted.read().decode("utf-8")
                return [
                    stripped
                    for line in text.splitlines()
                    if (stripped := line.strip())
                    and not stripped.startswith("#")
                ]
    return []


def _resolve_emscripten_graph(
    root_requirements: list[str],
    emscripten_env: dict[str, str],
) -> dict[str, SpecifierSet]:
    """Resolve all Emscripten runtime packages and combined specifiers."""
    specifiers: dict[str, SpecifierSet] = {}
    queue: deque[tuple[Requirement, str]] = deque()
    visited_extras: set[tuple[str, str]] = set()

    for line in root_requirements:
        req = Requirement(line)
        if req.marker is not None and not req.marker.evaluate(emscripten_env):
            continue
        queue.append((req, ""))

    while queue:
        req, _parent_extra = queue.popleft()
        canon_name = canonicalize_name(req.name)
        specifiers[canon_name] = specifiers.get(canon_name, SpecifierSet()) & (
            req.specifier
        )

        extras_to_expand = ["", *sorted(req.extras)]
        dist = importlib.metadata.distribution(req.name)
        for extra in extras_to_expand:
            key = (canon_name, extra)
            if key in visited_extras:
                continue
            visited_extras.add(key)

            env = {**emscripten_env, "extra": extra}
            for dep_str in dist.requires or []:
                dep_req = Requirement(dep_str)
                if dep_req.marker is not None:
                    if not dep_req.marker.evaluate(env):
                        continue
                elif extra:
                    # Unconditional deps were already expanded for extra == "".
                    continue
                queue.append((dep_req, extra))

    return specifiers


def _iter_dist_files(
    dist: importlib.metadata.Distribution,
) -> tuple[list[tuple[Path, str]], bool]:
    """Collect bundleable files and detect whether native binaries exist."""
    files = dist.files or []
    bundleable: list[tuple[Path, str]] = []
    has_native = False

    for pkg_path in files:
        rel_str = str(pkg_path)
        rel_path = Path(rel_str)
        if rel_path.is_absolute() or ".." in rel_path.parts:
            continue
        if "__pycache__" in rel_path.parts:
            continue
        suffix = rel_path.suffix.lower()
        if suffix in _NATIVE_SUFFIXES:
            has_native = True
            continue
        if suffix in _SKIPPED_SUFFIXES:
            continue
        abs_path = Path(str(dist.locate_file(pkg_path)))
        if abs_path.is_file():
            bundleable.append((abs_path, rel_path.as_posix()))

    return bundleable, has_native


def _sha256_file(path: Path) -> str:
    """Compute hex SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _ensure_pyodide_wheel(
    pkg_entry: dict[str, Any],
    pyodide_version: str,
    pyodide_dir: Path,
    cache_dir: Path,
) -> Path:
    """Ensure a Pyodide lockfile wheel exists in pyodide_dir with valid hash."""
    filename = pkg_entry.get("file_name")
    if not isinstance(filename, str) or not filename.endswith(".whl"):
        raise RuntimeError(f"Invalid wheel entry in lockfile: {pkg_entry}")

    expected_sha256 = pkg_entry.get("sha256")
    dest_path = pyodide_dir / filename
    if dest_path.is_file() and (
        not isinstance(expected_sha256, str)
        or _sha256_file(dest_path) == expected_sha256
    ):
        return dest_path

    cache_dir.mkdir(parents=True, exist_ok=True)
    cached_path = cache_dir / filename
    if not cached_path.is_file() or (
        isinstance(expected_sha256, str)
        and _sha256_file(cached_path) != expected_sha256
    ):
        url = _CDN_WHEEL_URL.format(version=pyodide_version, filename=filename)
        if not url.startswith("https://"):
            raise ValueError(f"Refusing non-HTTPS wheel URL: {url}")
        print(f"    Downloading Pyodide wheel {filename}...")
        tmp_download = cached_path.with_suffix(".whl.tmp")
        with (
            urllib.request.urlopen(url, timeout=60) as resp,  # noqa: S310
            tmp_download.open("wb") as out,
        ):
            shutil.copyfileobj(resp, out)
        if (
            isinstance(expected_sha256, str)
            and _sha256_file(tmp_download) != expected_sha256
        ):
            tmp_download.unlink(missing_ok=True)
            raise RuntimeError(
                f"SHA-256 mismatch for downloaded Pyodide wheel {filename}"
            )
        os.replace(tmp_download, cached_path)

    shutil.copy2(cached_path, dest_path)
    return dest_path


def _ensure_pyodide_package_tree(
    package_name: str,
    lock_packages: dict[str, dict[str, Any]],
    pyodide_version: str,
    pyodide_dir: Path,
    cache_dir: Path,
    downloaded: set[str],
) -> None:
    """Download a Pyodide package and all of its lockfile dependencies."""
    queue: deque[str] = deque([canonicalize_name(package_name)])
    while queue:
        name = queue.popleft()
        if name in downloaded:
            continue
        entry = lock_packages.get(name)
        if entry is None:
            raise RuntimeError(
                f"Pyodide lockfile is missing dependency '{name}'"
            )
        _ensure_pyodide_wheel(entry, pyodide_version, pyodide_dir, cache_dir)
        downloaded.add(name)
        for dep in entry.get("depends", []):
            if isinstance(dep, str):
                queue.append(canonicalize_name(dep))


def _repack_app_archive(
    archive_path: Path,
    pypackages_files: dict[str, Path],
    remaining_requirements: list[str],
) -> None:
    """Repack dist/app.tar.gz with __pypackages__/ and updated requirements."""
    fd, tmp_name = tempfile.mkstemp(
        prefix="app-", suffix=".tar.gz", dir=archive_path.parent
    )
    os.close(fd)
    tmp_path = Path(tmp_name)

    try:
        with (
            tarfile.open(archive_path, "r:gz") as src_tar,
            tarfile.open(tmp_path, "w:gz", format=tarfile.GNU_FORMAT) as dst,
        ):
            for member in src_tar.getmembers():
                normalised = member.name.lstrip("./")
                if normalised == "requirements.txt" or normalised.startswith(
                    "__pypackages__/"
                ):
                    continue
                if member.isfile():
                    fileobj = src_tar.extractfile(member)
                    dst.addfile(member, fileobj)
                else:
                    dst.addfile(member)

            for arc_rel, abs_path in sorted(pypackages_files.items()):
                dst.add(abs_path, arcname=f"__pypackages__/{arc_rel}")

            req_bytes = "".join(
                f"{line}\n" for line in remaining_requirements
            ).encode("utf-8")
            req_info = tarfile.TarInfo(name="requirements.txt")
            req_info.size = len(req_bytes)
            req_info.mode = 0o644
            dst.addfile(req_info, io.BytesIO(req_bytes))

        os.chmod(tmp_path, 0o644)
        os.replace(tmp_path, archive_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _patch_local_pyodide_url(index_path: Path) -> None:
    """Ensure local pyodideUrl uses a valid relative ES module specifier.

    Flet 1.1.3 writes ``flet.pyodideUrl="pyodide/pyodide.mjs";`` in
    ``--no-cdn`` mode, which browsers reject in ``import(self.pyodideUrl)``
    inside a module worker because bare specifiers without a leading ``./``
    or ``/`` cannot be resolved.
    """
    if not index_path.is_file():
        return
    content = index_path.read_text("utf-8")
    broken = 'flet.pyodideUrl="pyodide/pyodide.mjs";'
    fixed = 'flet.pyodideUrl="./pyodide/pyodide.mjs";'
    if broken in content:
        index_path.write_text(content.replace(broken, fixed), encoding="utf-8")


def bundle_web_dependencies(
    dist_dir: Path,
    cache_root: Path | None = None,
) -> None:
    """Bundle Python dependencies for an air-gapped Flet web distribution."""
    archive_path = dist_dir / "app.tar.gz"
    pyodide_dir = dist_dir / "pyodide"
    if not archive_path.is_file():
        raise FileNotFoundError(f"Expected web archive at {archive_path}")

    pyodide_version, python_version, lock_packages = _read_pyodide_metadata(
        pyodide_dir
    )
    cache_dir = (cache_root or _default_cache_root()) / pyodide_version
    emscripten_env = _make_emscripten_env(python_version)

    root_requirements = _read_requirements_from_archive(archive_path)
    resolved = _resolve_emscripten_graph(root_requirements, emscripten_env)

    pypackages_files: dict[str, Path] = {}
    remaining_requirements: list[str] = []
    downloaded_pyodide_pkgs: set[str] = set()

    for pkg_name, spec in sorted(resolved.items()):
        dist = importlib.metadata.distribution(pkg_name)
        bundleable, has_native = _iter_dist_files(dist)
        lock_entry = lock_packages.get(pkg_name)

        use_pyodide_wheel = False
        if has_native and lock_entry is not None:
            lock_ver_str = lock_entry.get("version")
            if isinstance(lock_ver_str, str) and spec.contains(
                Version(lock_ver_str), prereleases=True
            ):
                use_pyodide_wheel = True

        if use_pyodide_wheel and lock_entry is not None:
            lock_ver = str(lock_entry["version"])
            print(
                f"  Bundling {pkg_name}=={lock_ver} from Pyodide "
                f"{pyodide_version} into dist/pyodide/"
            )
            _ensure_pyodide_package_tree(
                pkg_name,
                lock_packages,
                pyodide_version,
                pyodide_dir,
                cache_dir,
                downloaded_pyodide_pkgs,
            )
            if pkg_name not in _AUTO_LOADED_PYODIDE_PACKAGES:
                remaining_requirements.append(f"{pkg_name}=={lock_ver}")
            continue

        has_python_code = any(rel.endswith(".py") for _abs, rel in bundleable)
        if not has_python_code:
            raise RuntimeError(
                f"Cannot bundle '{pkg_name}' ({spec}) for air-gapped web "
                f"build: package has native binaries with no pure-Python "
                f"fallback and no compatible wheel in Pyodide "
                f"{pyodide_version}."
            )

        mode = "pure-Python fallback" if has_native else "pure-Python"
        print(
            f"  Bundling {pkg_name}=={dist.version} ({mode}) into "
            f"__pypackages__/"
        )
        for abs_path, rel_posix in bundleable:
            pypackages_files[rel_posix] = abs_path

    # Always ensure msgpack is available locally for python-worker.js.
    for auto_pkg in sorted(_AUTO_LOADED_PYODIDE_PACKAGES):
        if (
            auto_pkg in lock_packages
            and auto_pkg not in downloaded_pyodide_pkgs
        ):
            _ensure_pyodide_package_tree(
                auto_pkg,
                lock_packages,
                pyodide_version,
                pyodide_dir,
                cache_dir,
                downloaded_pyodide_pkgs,
            )

    _repack_app_archive(archive_path, pypackages_files, remaining_requirements)
    _patch_local_pyodide_url(dist_dir / "index.html")
    print(
        f"==> Air-gapped bundle complete: {len(pypackages_files)} files in "
        f"__pypackages__/, {len(downloaded_pyodide_pkgs)} Pyodide wheel(s), "
        f"{len(remaining_requirements)} runtime requirement(s)."
    )


def main() -> None:
    """CLI entry point for bundling web dependencies."""
    parser = argparse.ArgumentParser(
        description=(
            "Bundle Python dependencies into a Flet --no-cdn web build."
        )
    )
    parser.add_argument(
        "dist_dir",
        type=Path,
        nargs="?",
        default=Path(__file__).resolve().parent.parent / "dist",
        help="Path to the Flet published dist/ directory.",
    )
    args = parser.parse_args()
    bundle_web_dependencies(args.dist_dir.resolve())


if __name__ == "__main__":
    main()
