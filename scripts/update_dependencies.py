#!/usr/bin/env python3
"""Script to automatically update dependencies.

Updates both pyproject.toml and requirements.txt to their latest versions.
"""

import json
import os
import re
import ssl
import urllib.request

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PYPROJECT_PATH = os.path.join(ROOT_DIR, "pyproject.toml")
REQUIREMENTS_PATH = os.path.join(ROOT_DIR, "src", "requirements.txt")


def fetch_latest_version(package_name: str) -> str | None:
    """Fetch the latest version of a package from PyPI."""
    # Ignore SSL verification errors if running in restricted environments
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    url = f"https://pypi.org/pypi/{package_name}/json"
    req = urllib.request.Request(  # noqa: S310
        url, headers={"User-Agent": "AmplifyP-Autoupdate"}
    )
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=10) as response:  # noqa: S310
            data = json.loads(response.read().decode())
            latest_version = data["info"]["version"]
            return str(latest_version)
    except Exception as e:
        print(f"Error fetching version for {package_name}: {e}")
        return None


def update_pyproject_toml() -> bool:
    """Analyse and update pyproject.toml dependencies."""
    if not os.path.exists(PYPROJECT_PATH):
        print(f"File not found: {PYPROJECT_PATH}")
        return False

    with open(PYPROJECT_PATH, encoding="utf-8") as f:
        content = f.read()

    # Pattern to match dependency line like:
    # "package>=version" or "package==version"
    pattern = re.compile(
        r'(?P<indent>\s*)"(?P<pkg>[a-zA-Z0-9_-]+)(?P<op>>=|==)(?P<ver>[0-9.]+)"(?P<comma>,?)'
    )

    modified = False

    def replace_dep(match: re.Match[str]) -> str:
        nonlocal modified
        indent = match.group("indent")
        pkg = match.group("pkg")
        op = match.group("op")
        old_ver = match.group("ver")
        comma = match.group("comma")

        latest = fetch_latest_version(pkg)
        if latest and latest != old_ver:
            print(f"Updating pyproject.toml: {pkg} {old_ver} -> {latest}")
            modified = True
            return f'{indent}"{pkg}{op}{latest}"{comma}'
        return str(match.group(0))

    new_content = pattern.sub(replace_dep, content)

    if modified:
        with open(PYPROJECT_PATH, "w", encoding="utf-8") as f:
            f.write(new_content)
        print("Updated pyproject.toml successfully.")
    else:
        print("pyproject.toml is already up to date.")

    return modified


def update_requirements_txt() -> bool:
    """Analyse and update requirements.txt dependencies."""
    if not os.path.exists(REQUIREMENTS_PATH):
        print(f"File not found: {REQUIREMENTS_PATH}")
        return False

    with open(REQUIREMENTS_PATH, encoding="utf-8") as f:
        content = f.read()

    pattern = re.compile(
        r"^(?P<pkg>[a-zA-Z0-9_-]+)(?P<op>>=|==)(?P<ver>[0-9.]+)$", re.MULTILINE
    )

    modified = False

    def replace_dep(match: re.Match[str]) -> str:
        nonlocal modified
        pkg = match.group("pkg")
        op = match.group("op")
        old_ver = match.group("ver")

        latest = fetch_latest_version(pkg)
        if latest and latest != old_ver:
            print(f"Updating requirements.txt: {pkg} {old_ver} -> {latest}")
            modified = True
            return f"{pkg}{op}{latest}"
        return str(match.group(0))

    new_content = pattern.sub(replace_dep, content)

    if modified:
        with open(REQUIREMENTS_PATH, "w", encoding="utf-8") as f:
            f.write(new_content)
        print("Updated requirements.txt successfully.")
    else:
        print("requirements.txt is already up to date.")

    return modified


def main() -> None:
    """Main entry point to execute the update process."""
    print("Starting dependency update analysis...")
    pyproject_updated = update_pyproject_toml()
    requirements_updated = update_requirements_txt()
    if pyproject_updated or requirements_updated:
        print("Dependencies updated.")
    else:
        print("All dependencies are up to date.")


if __name__ == "__main__":
    main()
