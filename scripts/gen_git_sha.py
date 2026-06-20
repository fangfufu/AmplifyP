#!/usr/bin/env python
"""Generate git_sha.py with commit SHA for frozen binary builds."""

import os
import subprocess

SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
SHA_PATH = os.path.join(SRC_DIR, "amplifyp", "gui", "git_sha.py")


def get_sha(cmd: list[str]) -> str:
    """Run a git command and return stripped output or 'unknown'."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=5)  # noqa: S603
        return r.stdout.strip() or "unknown"
    except OSError:
        return "unknown"


def main() -> None:
    """Generate git_sha.py with commit SHA for frozen binary builds."""
    short_sha = get_sha(["git", "rev-parse", "--short", "HEAD"])
    full_sha = get_sha(["git", "rev-parse", "HEAD"])
    os.makedirs(os.path.dirname(SHA_PATH), exist_ok=True)
    with open(SHA_PATH, "w") as f:
        f.write(f'GIT_SHA = "{short_sha}"\n')
        f.write(f'GIT_FULL_SHA = "{full_sha}"\n')
    print(f"Generated {SHA_PATH}: GIT_SHA={short_sha}")


if __name__ == "__main__":
    main()
