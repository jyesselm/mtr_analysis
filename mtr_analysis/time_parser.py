"""
Parse time values from directory names.

This module provides utilities for extracting time durations from
directory naming conventions used in experimental data organization.
"""

import os
import re


def get_minutes_from_dir_name(dir_name: str) -> int:
    """
    Extract the number of minutes from a directory name.

    Handles both 'min' and 'hr' suffixes, converting hours to minutes.
    Also handles bare numbers (assumed to be minutes).

    Args:
        dir_name: Directory name or path containing time information.
            Expected format: '*_t{number}{unit}' where unit is 'min', 'hr', or empty.

    Returns:
        The time value in minutes.

    Raises:
        ValueError: If the time cannot be parsed from the directory name.

    Examples:
        >>> get_minutes_from_dir_name('data/mtr1_mut_lib_t15min')
        15
        >>> get_minutes_from_dir_name('data/mtr1_mut_lib_t4hr')
        240
        >>> get_minutes_from_dir_name('data/mtr1_mut_lib_t0')
        0
    """
    base_name = _extract_base_name(dir_name)
    minutes = _parse_time_with_unit(base_name)
    if minutes is not None:
        return minutes
    minutes = _parse_bare_time(base_name)
    if minutes is not None:
        return minutes
    raise ValueError(f"Could not parse minutes from directory name: {dir_name}")


def _extract_base_name(dir_name: str) -> str:
    """
    Extract the base name from a directory path.

    Args:
        dir_name: Full path or directory name.

    Returns:
        The base name without parent directories.
    """
    if "/" in dir_name:
        return os.path.basename(dir_name)
    return dir_name


def _parse_time_with_unit(base_name: str) -> int | None:
    """
    Parse time value with explicit unit (min or hr).

    Args:
        base_name: Directory base name to parse.

    Returns:
        Time in minutes, or None if pattern not found.
    """
    match = re.search(r"_t(\d+)(min|hr)", base_name)
    if not match:
        return None
    value = int(match.group(1))
    unit = match.group(2)
    if unit == "hr":
        return value * 60
    return value


def _parse_bare_time(base_name: str) -> int | None:
    """
    Parse time value without unit suffix (assumed to be minutes).

    Args:
        base_name: Directory base name to parse.

    Returns:
        Time in minutes, or None if pattern not found.
    """
    match = re.search(r"_t(\d+)$", base_name)
    if not match:
        return None
    return int(match.group(1))
