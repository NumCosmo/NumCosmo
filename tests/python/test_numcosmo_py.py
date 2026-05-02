#!/usr/bin/env python
#
# test_py_numcosmo_py.py
#
# Sun Apr 27 14:57:44 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_numcosmo_py.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Test NumCosmo numcosmo_py tools."""

import pytest
from numcosmo_py import parse_options_strict, to_camel_case


@pytest.mark.parametrize(
    "input_list, expected",
    [
        (["a=b"], {"a": "b"}),
        (["key=value"], {"key": "value"}),
        (["k=v", "a=b"], {"k": "v", "a": "b"}),
        (["spaced=value with spaces"], {"spaced": "value with spaces"}),
        (["key", "=value"], {"key": "value"}),
        (["k", "=v"], {"k": "v"}),
        (["key", "=", "value"], {"key": "value"}),
        (["k=", "v"], {"k": "v"}),
        (["a=b", "c", "=d", "e=", "f"], {"a": "b", "c": "d", "e": "f"}),
        (["special-key=value-with!@#$%^&*()"], {"special-key": "value-with!@#$%^&*()"}),
    ],
)
def test_valid_cases(input_list, expected):
    """Test valid cases."""
    assert parse_options_strict(input_list) == expected


# Error cases with precise message matching
@pytest.mark.parametrize(
    "input_list, pattern",
    [
        # Empty items
        ([""], r"^Empty item$"),
        (["a", ""], r"^Empty item$"),
        # Malformed items
        (["=value"], r"^Malformed item: =value\.$"),
        (["key=value=extra"], r"^Malformed item: key=value=extra\.$"),
        (["a==b"], r"^Malformed item: a==b\.$"),
        (["=a=b"], r"^Malformed item: =a=b\.$"),
        (["a=b=c"], r"^Malformed item: a=b=c\.$"),
        # Invalid sequences
        (
            ["key", "=value", "extra"],
            r"^Unexpected end of input, expected more items\.$",
        ),
        (["key=", "=value"], r"^Malformed item: =value\.$"),
        (["key", "="], r"^Unexpected end of input, expected more items\.$"),
        (["key", "=", "=value"], r"^Malformed item: =value\.$"),
        (["key="], r"^Unexpected end of input, expected more items\.$"),
        (["key", "value"], r"^Malformed item: value\.$"),
        (["a", "b", "c"], r"^Malformed item: b\.$"),
        # Invalid value formats
        (["key", "=val=ue"], r"^Malformed item: =val=ue\.$"),
        (["key=", "val=ue"], r"^Malformed item: val=ue\.$"),
        (["key", "=", "=value"], r"^Malformed item: =value\.$"),
    ],
)
def test_error_cases(input_list, pattern):
    """Test error cases."""
    with pytest.raises(ValueError, match=pattern):
        parse_options_strict(input_list)


# Special cases with exact matching
@pytest.mark.parametrize(
    "input_list",
    [
        ["a=b", "a=b"],
        ["a=b", "a", "=", "b"],
        ["a=b", "a", "=", "b", "c", "=", "d"],
        ["a=b", "a", "=", "b", "c", "=", "d", "e", "=", "f"],
        ["a", "=b", "c", "=d", "e", "=", "f", "a=b", "c=d", "e=f"],
        ["a=", "b", "c=", "d", "e=", "f", "a=b", "c=d", "e=f"],
    ],
)
def test_duplicate_keys(input_list):
    """Test duplicate keys."""
    with pytest.raises(ValueError, match=r"^Duplicate key: a$"):
        parse_options_strict(input_list)


def test_single_equal_only():
    """Test single equal only."""
    with pytest.raises(ValueError, match=r"^Malformed item: =\.$"):
        parse_options_strict(["="])


def test_equal_only_sequence():
    """Test equal only sequence."""
    with pytest.raises(ValueError, match=r"^Malformed item: =\.$"):
        parse_options_strict(["a", "=", "="])


def test_complex_malformed():
    """Test complex malformed sequence."""
    with pytest.raises(ValueError, match=r"^Malformed item: =b=c\.$"):
        parse_options_strict(["a", "=b=c"])


def test_unexpected_end():
    """Test unexpected end of input."""
    with pytest.raises(
        ValueError, match=r"^Unexpected end of input, expected more items\.$"
    ):
        parse_options_strict(["key"])


def test_equal_right_with_value_containing_equal():
    """Test equal right with value containing equal."""
    with pytest.raises(ValueError, match=r"^Malformed item: val=ue\.$"):
        parse_options_strict(["key=", "val=ue"])


@pytest.mark.parametrize(
    "input_str, expected",
    [
        # Basic cases
        ("hello_world", "HelloWorld"),
        ("single", "Single"),
        ("", ""),
        # Whitespace handling
        ("  leading_space", "LeadingSpace"),
        ("trailing_space  ", "TrailingSpace"),
        ("  both  ", "Both"),
        # Case variations
        ("ALREADY_CAPS", "AlreadyCaps"),
        ("MixED_CaSe", "MixedCase"),
        # Special characters
        ("hello2_world", "Hello2World"),
        ("special_char!", "SpecialChar!"),
        ("num_123_abc", "Num123Abc"),
        # Edge cases with underscores
        ("__leading_underscore", "LeadingUnderscore"),
        ("trailing_underscore__", "TrailingUnderscore"),
        ("multiple__underscores", "MultipleUnderscores"),
        ("_", ""),
        ("___", ""),
        # Complex cases
        ("_partial_CamelCase_", "PartialCamelcase"),
        ("  __multi_word_example__  ", "MultiWordExample"),
    ],
)
def test_to_camel_case(input_str, expected):
    """Test to_camel_case function."""
    assert to_camel_case(input_str) == expected
