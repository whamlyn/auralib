"""
auralib module containing helper functions for reading data from Petrel exports.

Author:   Wes Hamlyn
Created:  26-Apr-2021
Last Mod: 26-Apr-2021

"""

from __future__ import annotations

import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d, PchipInterpolator
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


@dataclass(frozen=True)
class PetrelField:
    """Description of one field in a Petrel ASCII file."""

    name: str
    data_type: str | None
    original_definition: str


@dataclass
class PetrelAsciiFile:
    """Parsed contents of a Petrel ASCII file."""

    comments: list[str]
    version: str | None
    fields: list[PetrelField]
    data: pd.DataFrame

    @property
    def metadata(self) -> dict[str, str]:
        """
        Convert comments of the form:

            # Description: value

        into a metadata dictionary.
        """
        result: dict[str, str] = {}

        for comment in self.comments:
            text = comment.removeprefix("#").strip()

            if ":" in text:
                key, value = text.split(":", maxsplit=1)
                result[key.strip()] = value.strip()

        return result


# A quoted value is recognized only when the token begins with a quote.
# This allows unquoted DMS coordinates such as:
#
#     55°27'53.3161"N
#
# to be treated as one token even though they contain a double quote.
_TOKEN_PATTERN = re.compile(r'"([^"]*)"|(\S+)')

# Petrel types that can appear before a comma in the header.
_KNOWN_TYPES = {
    "FLOAT",
    "DOUBLE",
    "REAL",
    "INT",
    "INTEGER",
    "LONG",
    "STRING",
    "DATE",
    "BOOL",
    "BOOLEAN",
}


def _split_data_line(line: str) -> list[str]:
    """
    Split a Petrel data row while preserving spaces inside quoted strings.

    Examples
    --------
    "Well name with spaces" -> one value
    ""                      -> empty string
    502813.00               -> one unquoted value
    55°27'53.3161"N         -> one unquoted value
    """
    values: list[str] = []

    for match in _TOKEN_PATTERN.finditer(line):
        quoted_value = match.group(1)
        unquoted_value = match.group(2)

        if quoted_value is not None:
            values.append(quoted_value)
        else:
            values.append(unquoted_value)

    return values


def _parse_field_definition(definition: str) -> PetrelField:
    """
    Parse either:

        FLOAT,Surface elevation

    or an untyped field:

        Surface X
    """
    if "," in definition:
        possible_type, field_name = definition.split(",", maxsplit=1)
        possible_type = possible_type.strip().upper()

        if possible_type in _KNOWN_TYPES:
            return PetrelField(
                name=field_name.strip(),
                data_type=possible_type,
                original_definition=definition,
            )

    return PetrelField(
        name=definition.strip(),
        data_type=None,
        original_definition=definition,
    )


def _parse_date(value: str) -> pd.Timestamp | pd.NaT:
    """
    Parse the day/month/year dates used by the supplied Petrel export.
    """
    if not value:
        return pd.NaT

    try:
        return pd.Timestamp(datetime.strptime(value, "%d/%m/%Y"))
    except ValueError:
        return pd.NaT


def _convert_series(
    series: pd.Series,
    data_type: str | None,
    *,
    parse_untyped_values: bool,
) -> pd.Series:
    """Convert a DataFrame column according to its Petrel type."""

    if data_type in {"FLOAT", "DOUBLE", "REAL"}:
        return pd.to_numeric(series, errors="coerce")

    if data_type in {"INT", "INTEGER", "LONG"}:
        return pd.to_numeric(series, errors="coerce").astype("Int64")

    if data_type in {"BOOL", "BOOLEAN"}:
        boolean_values = {
            "true": True,
            "false": False,
            "yes": True,
            "no": False,
            "1": True,
            "0": False,
        }

        return series.map(
            lambda value: (
                boolean_values.get(str(value).strip().lower(), pd.NA)
                if pd.notna(value)
                else pd.NA
            )
        ).astype("boolean")

    if data_type == "DATE":
        return series.map(
            lambda value: _parse_date(str(value))
            if pd.notna(value)
            else pd.NaT
        )

    if data_type == "STRING":
        return series.astype("string")

    if not parse_untyped_values:
        return series

    # For untyped fields, cautiously convert only when every non-null value
    # can be interpreted as numeric. Coordinates written as DMS remain text.
    non_null = series.dropna()

    if non_null.empty:
        return series

    numeric = pd.to_numeric(non_null, errors="coerce")

    if numeric.notna().all():
        converted = pd.to_numeric(series, errors="coerce")

        # Preserve integer columns as nullable integers.
        if converted.dropna().mod(1).eq(0).all():
            return converted.astype("Int64")

        return converted.astype("Float64")

    return series.astype("string")


def read_petrel_headers(
    filename: str | Path,
    *,
    encoding: str = "utf-8-sig",
    null_values: Iterable[str] = ("NULL", "-999"),
    empty_string_is_null: bool = False,
    convert_types: bool = True,
    parse_untyped_values: bool = True,
    strict: bool = True,
    ) -> PetrelAsciiFile:
    """
    Read a Petrel well header export containing comments, a header, and data.

    Parameters
    ----------
    filename:
        Path to the Petrel ASCII file.

    encoding:
        Text encoding used to open the file. ``utf-8-sig`` also handles
        files that begin with a UTF-8 byte-order mark.

    null_values:
        Values to interpret as missing. Remove ``"-999"`` from this
        collection if it should be retained as an actual numeric value.

    empty_string_is_null:
        Interpret quoted empty strings (``""``) as missing values.

    convert_types:
        Convert explicitly typed Petrel fields to pandas-compatible types.

    parse_untyped_values:
        Attempt numeric conversion for fields without explicit Petrel types.

    strict:
        Raise an error when the number of values in a data row differs from
        the number of header fields. If False, short rows are padded and long
        rows are truncated.

    Returns
    -------
    PetrelAsciiFile
        Parsed comments, version, field definitions, and DataFrame.

    Raises
    ------
    ValueError
        If required header markers are missing or a row has the wrong number
        of values in strict mode.
    """
    path = Path(filename)

    with path.open("r", encoding=encoding, errors="replace") as file:
        lines = [line.rstrip("\r\n") for line in file]

    comments: list[str] = []
    version: str | None = None
    header_definitions: list[str] = []
    raw_rows: list[list[str]] = []

    state = "preamble"
    header_found = False
    header_ended = False

    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()

        if not stripped:
            continue

        if state == "preamble":
            if stripped.startswith("#"):
                comments.append(stripped)
                continue

            if stripped.upper().startswith("VERSION"):
                parts = stripped.split(maxsplit=1)
                version = parts[1].strip() if len(parts) == 2 else None
                continue

            if stripped.upper() == "BEGIN HEADER":
                state = "header"
                header_found = True
                continue

            # Permit comments or other non-data preamble records, but fail
            # if a data-looking line appears before BEGIN HEADER.
            if strict:
                raise ValueError(
                    f"Unexpected content before BEGIN HEADER at "
                    f"line {line_number}: {line!r}"
                )

        elif state == "header":
            if stripped.upper() == "END HEADER":
                state = "data"
                header_ended = True
                continue

            header_definitions.append(stripped)

        elif state == "data":
            if stripped.startswith("#"):
                # Preserve comments occurring later in the file as well.
                comments.append(stripped)
                continue

            values = _split_data_line(stripped)
            expected = len(header_definitions)
            actual = len(values)

            if actual != expected:
                message = (
                    f"Line {line_number} contains {actual} values, "
                    f"but the header defines {expected} fields."
                )

                if strict:
                    raise ValueError(message)

                if actual < expected:
                    values.extend([None] * (expected - actual))
                else:
                    values = values[:expected]

            raw_rows.append(values)

    if not header_found:
        raise ValueError("The file does not contain a BEGIN HEADER marker.")

    if not header_ended:
        raise ValueError("The file does not contain an END HEADER marker.")

    if not header_definitions:
        raise ValueError("No fields were found between the header markers.")

    fields = [
        _parse_field_definition(definition)
        for definition in header_definitions
    ]

    column_names = [field.name for field in fields]

    duplicates = pd.Index(column_names)[pd.Index(column_names).duplicated()]
    if len(duplicates) > 0:
        duplicate_names = sorted(set(duplicates))
        raise ValueError(
            "Duplicate field names were found in the header: "
            f"{duplicate_names}"
        )

    frame = pd.DataFrame(raw_rows, columns=column_names)

    null_set = {str(value) for value in null_values}

    def replace_null(value: Any) -> Any:
        if value is None:
            return pd.NA

        if value in null_set:
            return pd.NA

        if empty_string_is_null and value == "":
            return pd.NA

        return value

    frame = frame.map(replace_null)

    if convert_types:
        for field in fields:
            frame[field.name] = _convert_series(
                frame[field.name],
                field.data_type,
                parse_untyped_values=parse_untyped_values,
            )

    return PetrelAsciiFile(
        comments=comments,
        version=version,
        fields=fields,
        data=frame,
    )


def read_well_headers_old(infile):
    """
    Function to read well headers exported from Petrel.
    Note: This is very early code and definitely does not support all
          perturbations of Petrel well header exports.
    """
        
    with open(infile, 'r') as fd:
        buf = fd.readlines()
        
    data = buf[52:]
    
    
    wellname = []
    uwi = []
    utmx = []
    utmy = []
    kb = []
    se = []
    for line in data:
        cur_wellname = line.strip().split('"')[1]
        cur_uwi = line.strip().split('"')[3]
        cur_utmx = float(line.strip().split('"')[4].split(' ')[2])
        cur_utmy = float(line.strip().split('"')[4].split(' ')[3])
        cur_kb = float(line.strip().split('"')[10].strip())
        cur_se = float(line.strip().split('"')[-3].split()[-1])
        
        wellname.append(cur_wellname)
        uwi.append(cur_uwi)
        utmx.append(cur_utmx)
        utmy.append(cur_utmy)
        kb.append(cur_kb)
        se.append(cur_se)
    
    head = pd.DataFrame()
    head['uwi'] = uwi
    head['wellname'] = wellname
    head['kb'] = kb
    head['se'] = se
    head['utmx'] = utmx
    head['utmy'] = utmy

    return head


def read_well_markers(infile):
    """
    Function to read well markers exported from Petrel.
    Note: This is very early code and definitely does not support all
          perturbations of Petrel well marker exports.
    """
        
    with open(infile, 'r') as fd:
        buf = fd.readlines()
    
    data_start_row = np.nonzero(np.array(buf)=='END HEADER\n')[0][0] + 1
    data = buf[data_start_row:]
    
    utmx = []
    utmy = []
    elev = []
    md = []
    name = []
    well = []
    
    for line in data:
        
        cur_utmx = float(line.strip().split('"')[0].split()[0])
        cur_utmy = float(line.strip().split('"')[0].split()[1])
        cur_elev = float(line.strip().split('"')[0].split()[2])
        cur_md = float(line.strip().split('"')[0].split()[6])
        cur_name = line.strip().split('"')[1]
        cur_well = line.strip().split('"')[3]
        
        utmx.append(cur_utmx)
        utmy.append(cur_utmy)
        elev.append(cur_elev)
        md.append(cur_md)
        name.append(cur_name)
        well.append(cur_well)
    
    tops = pd.DataFrame({'name': name,
                         'well': well,
                         'md': md,
                         'elev': elev,
                         'utmx': utmx,
                         'utmy': utmy})
    
    return tops


def parse_deviation(infile):
    """
    Read a deviation survey exported from Petrel.

    Parameters
    ----------
    infile : string
        Full path to input deviation survey file exported from Petrel.

    Returns
    -------
    dev : dictionary
        Dictionary containing uwi, md, tvd, dx, dy, x, y values read from the 
        input directional survey file.
    
    Written by: Wes Hamlyn
    Last Modified: Oct 3, 2022
    
    """

    dev = {}
    
    dev_uwi = []
    dev_md = []
    dev_tvd = []
    dev_dx = []
    dev_dy = []
    dev_x = []
    dev_y = []
    dev_azimtn = []
    dev_azimgn = []
    dev_incl = []
    
    with open(infile, 'r') as fd:
        buf = fd.readlines()
        
        start_of_data = -1
        count = 0
        
        for line in buf:
        
            if line[0] == '#':
                count += 1
                
                if 'WELL NAME:' in line:
                    cur_uwi = line.split(':')[1].strip()
                    dev['uwi'] = cur_uwi
                
                if 'WELL HEAD X-COORDINATE:' in line:
                    cur_x = line.strip().split(':')[1].split()[0]
                    dev['top_x'] = float(cur_x)
                
                if 'WELL HEAD Y-COORDINATE:' in line:
                    cur_y = line.strip().split(':')[1].split()[0]
                    dev['top_y'] = float(cur_y)
                
                if 'WELL DATUM (' in line:
                    cur_kb = line.split(':')[1].split()[0]
                    dev['kb'] = float(cur_kb)
                    
                    # no ground elevations so we'll approximate it
                    dev['grd'] = float(cur_kb)-5
                    
            else:
                start_of_data = count+2
                break
            
        dev_cols = buf[start_of_data-2].strip().split()
        idx_md = dev_cols.index('MD')
        idx_tvd = dev_cols.index('TVD')
        idx_dx = dev_cols.index('DX')
        idx_dy = dev_cols.index('DY')
        idx_x = dev_cols.index('X')
        idx_y = dev_cols.index('Y')
        idx_azimtn = dev_cols.index('AZIM_TN')
        idx_azimgn = dev_cols.index('AZIM_GN')
        idx_incl = dev_cols.index('INCL')
        
        for line in buf[start_of_data:]:
            
            line = line.strip().split()
            
            dev_uwi.append(cur_uwi)
            dev_md.append(float(line[idx_md]))
            dev_tvd.append(float(line[idx_tvd]))
            dev_dx.append(float(line[idx_dx]))
            dev_dy.append(float(line[idx_dy]))
            dev_x.append(float(line[idx_x]))
            dev_y.append(float(line[idx_y]))
            dev_azimtn.append(float(line[idx_azimtn]))
            dev_azimgn.append(float(line[idx_azimgn]))
            dev_incl.append(float(line[idx_incl]))
        
        #dev['dev_uwi'] = dev_uwi
        dev['md'] = np.array(dev_md)
        dev['tvd'] = np.array(dev_tvd)
        dev['dx'] = np.array(dev_dx)
        dev['dy'] = np.array(dev_dy)
        dev['x'] = np.array(dev_x)
        dev['y'] = np.array(dev_y)
        dev['azimtn'] = np.array(dev_azimtn)
        dev['azimgn'] = np.array(dev_azimgn)
        dev['incl'] = np.array(dev_incl)
        
    return dev


class PetrelDev():
    """
    Class to simplify working with directional surveys exported from Petrel
    """
    
    def __init__(self, infile):
        """
        Constructor method for class
        """
        self.infile = infile
        
        try:
            dev = parse_deviation(infile)
            
            self.uwi = dev['uwi']
            self.top_x = dev['top_x']
            self.top_y = dev['top_y']
            self.kb = dev['kb']
            self.grd = dev['grd']
            self.md = dev['md']
            self.tvd = dev['tvd']
            self.dx = dev['dx']
            self.dy = dev['dy']
            self.x = dev['x']
            self.y = dev['y']
            
            if np.min(self.md) < 0.0:
                """
                Remove and deviation points that occur above MD=0.0
                """
                idx = np.nonzero(self.md<0.0)
                self.md = self.md[idx]
                self.dx = self.dx[idx]
                self.dy = self.dy[idx]
                self.x = self.x[idx]
                self.y = self.t[idx]
                self.tvd = self.tvd[idx]
                
            if self.md[0] > 0.0:
                """
                Add a point at MD=0 if the first deviation survey point is at
                some deeper value than MD=0 to minimize the chance of weird
                interpolated values.
                """
                self.md = np.hstack((0.0, self.md))
                self.dx = np.hstack((0.0, self.dx))
                self.dy = np.hstack((0.0, self.dy))
                self.x = np.hstack((self.top_x, self.x))
                self.y = np.hstack((self.top_y, self.y))
                self.tvd = np.hstack((0.0, self.tvd))
                
            
            self.fdx = interp1d(self.md, self.dx, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
            self.fdy = interp1d(self.md, self.dy, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
            self.fx = interp1d(self.md, self.x, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
            self.fy = interp1d(self.md, self.y, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
            self.ftvd = interp1d(self.md, self.tvd, kind='linear',
                                bounds_error=False, fill_value='extrapolate')
        
        except Exception as e:
            print('Failed to initialized Petrel Deviation Object (%s)' % e)
    
    
    def md2xy(self, md):
        """
        Convert md to x/y
        """
        
        x = self.fx(md)
        y = self.fy(md)
        if type(md) in (int, float):
            x = float(x)
            y = float(y)
            
        return x, y
    
    def md2offset(self, md):
        """
        Convert md to dx/dy
        """
        dx = self.fdx(md)
        dy = self.fdy(md)
        if type(md) in (int, float):
            dx = float(dx)
            dy = float(dy)
            
        return dx, dy
       
    def md2tvd(self, md):
        """
        Convert md to tvd
        """
        tvd = self.ftvd(md)
        if type(md) in (int, float):
            tvd = float(tvd)
        
        return tvd
    
    def md2tvdss(self, md):
        """
        Convert md to tvdss
        """
        tvd = self.md2tvd(md)
        tvdss = tvd - self.kb
        return tvdss


def write_petrel_pointset(outfile, x, y, z):
    """
    Write a simple x, y, z pointset in Petrel Pointset format. No attributes
    are supported by this function.
    """
    
    with open(outfile, 'w') as fd:
        # write header
        fd.write('# Petrel Points with attributes\n')
        fd.write('# Unit in X and Y direction: m\n')
        fd.write('# Unit in depth: m\n')
        fd.write('VERSION 1\n')
        fd.write('BEGIN HEADER\n')
        fd.write('X\n')
        fd.write('Y\n')
        fd.write('Z\n')
        fd.write('END HEADER\n')
        
        for i in range(len(x)):
            line = '%.2f %.2f %.2f\n' % (x[i], y[i], z[i])
            fd.write(line)


def read_grid_earthvision(infile):
    """
    Convenience function to read grid points exported from Petrel in
    Earthvision format
    """
        
    with open(infile, 'r') as fd:
        buf = fd.readlines()
    
    x = []
    y = []
    z = []
    c = []
    r = []
    for line in buf:
        if line[0]=='#':
            next
        
        else:
            line = line.strip().split()
            x.append(float(line[0]))
            y.append(float(line[1]))
            z.append(float(line[2]))
            c.append(int(line[3]))
            r.append(int(line[4]))
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    c = np.array(c)
    r = np.array(r)
    
    rmin = r.min()
    rmax = r.max()
    cmin = c.min()
    cmax = c.max()
    
    nr = r.max() - r.min() + 1
    nc = c.max() - c.min() + 1
    
    ridx = r - rmin
    cidx = c - cmin

    data_array = np.full((nr, nc), np.nan)
    data_array[ridx, cidx] = z[:]
    
    return x, y, z, r, c, data_array