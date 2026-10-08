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
"""Pre-commit check to ensure YAML lines do not exceed the length limit."""

import argparse
import sys
from pathlib import Path


def check_file(path: Path, max_length: int = 80) -> list[str]:
    """Check a file for lines exceeding the maximum length.

    Args:
        path: Path to the YAML file to check.
        max_length: Maximum allowed line length in characters.

    Returns:
        List of formatted error messages for lines exceeding max_length.
    """
    errors: list[str] = []
    try:
        with open(path, encoding="utf-8") as file_handle:
            for line_number, raw_line in enumerate(file_handle, start=1):
                clean_line = raw_line.rstrip("\r\n")
                line_len = len(clean_line)
                if line_len > max_length:
                    errors.append(
                        f"{path}:{line_number}: line exceeds {max_length} "
                        f"characters ({line_len} chars)"
                    )
    except (UnicodeDecodeError, OSError) as err:
        errors.append(f"{path}: failed to read file: {err}")
    return errors


def main(argv: list[str] | None = None) -> int:
    """Entry point for checking YAML file line lengths.

    Args:
        argv: Optional argument list. Defaults to sys.argv[1:].

    Returns:
        Exit code (0 on success, 1 on line length violations).
    """
    parser = argparse.ArgumentParser(
        description="Verify YAML files do not exceed line length limit."
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=80,
        help="Maximum permissible line length in characters (default: 80).",
    )
    parser.add_argument(
        "filenames",
        nargs="*",
        help="Filenames to inspect.",
    )
    parsed_args = parser.parse_args(argv)

    all_errors: list[str] = []
    for filename in parsed_args.filenames:
        target_path = Path(filename)
        if target_path.is_file():
            all_errors.extend(
                check_file(target_path, max_length=parsed_args.max_length)
            )

    if all_errors:
        for error_message in all_errors:
            print(error_message, file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
