#!/usr/bin/env python
"""Helper script to load state and update screenshots in docs/gui/images."""

import os
import subprocess
import sys


def main() -> None:
    """Run AmplifyP to capture screenshots of simple.yaml.

    Captures screenshots at 1600x900 into docs/gui/images.
    """
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    state_file = os.path.join(
        repo_root, "tests", "examples", "save_states", "simple.yaml"
    )
    docs_images_dir = os.path.join(repo_root, "docs", "gui", "images")
    main_script = os.path.join(repo_root, "src", "main.py")

    cmd = [
        sys.executable,
        main_script,
        "--state",
        state_file,
        "--window-width",
        "1600",
        "--window-height",
        "900",
        "--screenshots",
        "--screenshots-dir",
        docs_images_dir,
        "--auto-close",
    ]

    print(f"Executing screenshot update: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=repo_root, check=False)  # noqa: S603
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
