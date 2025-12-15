"""Tests for time_parser module."""

import pytest

from mtr_analysis.time_parser import (
    _extract_base_name,
    _parse_bare_time,
    _parse_time_with_unit,
    get_minutes_from_dir_name,
)


class TestGetMinutesFromDirName:
    """Tests for get_minutes_from_dir_name function."""

    def test_minutes_suffix(self) -> None:
        """Test parsing directory name with 'min' suffix."""
        result = get_minutes_from_dir_name("data/mtr1_mut_lib_t15min")
        assert result == 15

    def test_hours_suffix(self) -> None:
        """Test parsing directory name with 'hr' suffix."""
        result = get_minutes_from_dir_name("data/mtr1_mut_lib_t4hr")
        assert result == 240

    def test_bare_number(self) -> None:
        """Test parsing directory name with bare number (assumed minutes)."""
        result = get_minutes_from_dir_name("data/mtr1_mut_lib_t0")
        assert result == 0

    def test_no_path(self) -> None:
        """Test parsing directory name without path prefix."""
        result = get_minutes_from_dir_name("mtr1_mut_lib_t30min")
        assert result == 30

    def test_large_time_value(self) -> None:
        """Test parsing large time values."""
        result = get_minutes_from_dir_name("experiment_t1440min")
        assert result == 1440

    def test_hours_conversion(self) -> None:
        """Test that hours are correctly converted to minutes."""
        result = get_minutes_from_dir_name("sample_t2hr")
        assert result == 120

    def test_invalid_format_raises(self) -> None:
        """Test that invalid format raises ValueError."""
        with pytest.raises(ValueError, match="Could not parse minutes"):
            get_minutes_from_dir_name("invalid_directory_name")

    def test_missing_time_raises(self) -> None:
        """Test that missing time pattern raises ValueError."""
        with pytest.raises(ValueError, match="Could not parse minutes"):
            get_minutes_from_dir_name("data/no_time_here")


class TestExtractBaseName:
    """Tests for _extract_base_name helper function."""

    def test_with_path(self) -> None:
        """Test extraction from full path."""
        result = _extract_base_name("path/to/directory")
        assert result == "directory"

    def test_without_path(self) -> None:
        """Test extraction from bare name."""
        result = _extract_base_name("directory")
        assert result == "directory"

    def test_nested_path(self) -> None:
        """Test extraction from deeply nested path."""
        result = _extract_base_name("a/b/c/d/final")
        assert result == "final"


class TestParseTimeWithUnit:
    """Tests for _parse_time_with_unit helper function."""

    def test_minutes(self) -> None:
        """Test parsing minutes suffix."""
        result = _parse_time_with_unit("sample_t45min")
        assert result == 45

    def test_hours(self) -> None:
        """Test parsing hours suffix."""
        result = _parse_time_with_unit("sample_t3hr")
        assert result == 180

    def test_no_match(self) -> None:
        """Test that non-matching strings return None."""
        result = _parse_time_with_unit("no_time_here")
        assert result is None

    def test_zero_value(self) -> None:
        """Test parsing zero time value."""
        result = _parse_time_with_unit("sample_t0min")
        assert result == 0


class TestParseBareTime:
    """Tests for _parse_bare_time helper function."""

    def test_bare_number(self) -> None:
        """Test parsing bare number at end."""
        result = _parse_bare_time("sample_t123")
        assert result == 123

    def test_zero(self) -> None:
        """Test parsing zero."""
        result = _parse_bare_time("sample_t0")
        assert result == 0

    def test_no_match(self) -> None:
        """Test that non-matching strings return None."""
        result = _parse_bare_time("no_time_here")
        assert result is None

    def test_number_not_at_end(self) -> None:
        """Test that numbers not at end don't match."""
        result = _parse_bare_time("sample_t5_extra")
        assert result is None
