"""
segb.py
=======

SEG-B Format-B field-record reader based on:

    Northwood, E. J., Weisinger, R. C. & Bradley, J. J. (1967),
    "Recommended standards for digital tape formats."

This version targets SEG-B Format B and can infer scan geometry from the data stream when legacy/vendor headers do not match the generic 1967 header allocation.

Scope
-----
Implemented:
    * 16-byte common/standard header
    * packed-BCD decoding
    * inferred Format-B scan geometry
    * optional header sizing for Format B
    * 4-byte sync group parsing
    * 10-byte auxiliary group parsing
    * 10-byte seismic groups (4 channels/group)
    * four 4-bit gain codes per seismic group
    * one's- or two's-complement 15-bit seismic mantissa decoding
    * demultiplexing to NumPy arrays of shape (n_channels, n_samples)
    * optional gain scaling according to Q * 2**(-G)
    * access to raw scan bytes and auxiliary words

Important
---------
The 1967 standard describes physical magnetic-tape records and inter-record
gaps. A modern disk image, tape-image container, or vendor conversion may
preserve those physical records in a container-specific way. This module
parses a contiguous SEG-B *record payload*. If your file contains tape-image
framing, record-length words, tape marks, etc., that outer container must be
handled separately.

This module deliberately keeps low-level raw values accessible. That is useful
when diagnosing legacy field data whose implementation differs from the
published recommendation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, Union, Any

import numpy as np


PathLike = Union[str, Path]
NegativeEncoding = Literal["ones", "twos"]


# -----------------------------------------------------------------------------
# BCD helpers
# -----------------------------------------------------------------------------

def _bcd_nibble(n: int) -> int:
    if not 0 <= n <= 9:
        raise ValueError(f"Invalid BCD digit: 0x{n:X}")
    return n


def decode_bcd_byte(value: int) -> int:
    """Decode one packed-BCD byte, e.g. 0x27 -> 27."""
    hi = _bcd_nibble((value >> 4) & 0x0F)
    lo = _bcd_nibble(value & 0x0F)
    return 10 * hi + lo


def decode_bcd(data: bytes) -> int:
    """Decode a sequence of packed-BCD bytes as one integer."""
    value = 0
    for b in data:
        value = value * 100 + decode_bcd_byte(b)
    return value


def decode_bcd_digits(data: bytes) -> str:
    """Return packed BCD as a digit string, preserving leading zeros."""
    out = []
    for b in data:
        hi = _bcd_nibble((b >> 4) & 0x0F)
        lo = _bcd_nibble(b & 0x0F)
        out.extend((str(hi), str(lo)))
    return "".join(out)


# -----------------------------------------------------------------------------
# Signed sample decoding
# -----------------------------------------------------------------------------

def decode_ones_complement_u16(word: int) -> int:
    """
    Decode a 16-bit one's-complement integer.

    Both +0 (0x0000) and -0 (0xFFFF) are returned as zero.
    """
    word &= 0xFFFF
    if word == 0xFFFF:
        return 0
    if word & 0x8000:
        return -((~word) & 0xFFFF)
    return word


def decode_twos_complement_u16(word: int) -> int:
    """Decode a 16-bit two's-complement integer."""
    word &= 0xFFFF
    return word - 0x10000 if word & 0x8000 else word


def decode_signed_nbit(value: int, nbits: int, encoding: NegativeEncoding) -> int:
    """Decode an n-bit one's- or two's-complement signed integer."""
    mask = (1 << nbits) - 1
    sign = 1 << (nbits - 1)
    value &= mask
    if encoding == "ones":
        if value == mask:       # negative zero
            return 0
        return -((~value) & mask) if value & sign else value
    if encoding == "twos":
        return value - (1 << nbits) if value & sign else value
    raise ValueError("encoding must be 'ones' or 'twos'")


def decode_sample_word(word: int, encoding: NegativeEncoding) -> int:
    """Decode a SEG-B seismic sample word.

    SEG-B places Qs,Q1..Q14 in the most-significant 15 bits of the
    two-byte word.  The conventional least-significant bit (document bit 7
    of byte 2) is therefore not part of the seismic mantissa.  Shift it away
    and decode the resulting 15-bit signed quantity.
    """
    mantissa15 = (word & 0xFFFF) >> 1
    return decode_signed_nbit(mantissa15, 15, encoding)


def sample_word_to_fraction(word: int, encoding: NegativeEncoding) -> float:
    """Return Qs.Q1..Q14 as a normalized signed fraction (K excluded)."""
    return decode_sample_word(word, encoding) / float(1 << 14)


# -----------------------------------------------------------------------------
# Format-B auxiliary decoding
# -----------------------------------------------------------------------------

def decode_format_b_time_counter(b0: int, b1: int) -> int:
    """Decode the 14-bit, 1-kHz Format-B auxiliary time counter.

    The SEG document numbers bits 0..7 from MSB to LSB and reserves the
    conventional LSB of the first auxiliary byte as an anti-sync bit.  The
    observed Format-B packing therefore stores the high counter bits in
    ``b0 >> 1`` and the low eight bits in ``b1``.  The returned value is in
    milliseconds.
    """
    return ((int(b0) >> 1) << 8) | int(b1)


def decode_format_b_time_counter_word(word: int) -> int:
    """Decode a raw two-byte Format-B time-counter word to milliseconds."""
    word &= 0xFFFF
    return decode_format_b_time_counter((word >> 8) & 0xFF, word & 0xFF)

# -----------------------------------------------------------------------------
# Header
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class SegBHeader:
    file_number: int
    format_code: int
    data_id: str
    bytes_per_scan: int
    sample_interval_ms: int
    manufacturer_code: int
    equipment_serial: int
    raw_standard: bytes

    @property
    def dt_seconds(self) -> float:
        return self.sample_interval_ms / 1000.0


@dataclass(frozen=True)
class SegBOptionalHeader:
    record_length_seconds: Optional[int]
    gain_control_mode: Optional[int]
    record_type: Optional[int]
    low_cut: Optional[int]
    low_cut_slope_code: Optional[int]
    high_cut: Optional[int]
    high_cut_slope_code: Optional[int]
    special_filter: Optional[int]
    alias_filter: Optional[int]
    channel_gain_bytes: bytes
    raw: bytes


def parse_standard_header(buf: bytes) -> SegBHeader:
    """
    Parse the mandatory 16-byte common SEG-A/SEG-B header.

    Byte numbering below follows the 1967 specification (1-based).
    """
    if len(buf) < 16:
        raise ValueError(
            f"SEG-B standard header requires 16 bytes; received {len(buf)}"
        )

    h = buf[:16]

    file_number = decode_bcd(h[0:2])       # bytes 1-2
    format_code = decode_bcd(h[2:4])       # bytes 3-4
    data_id = decode_bcd_digits(h[4:10])   # bytes 5-10

    # Byte 11 contains B1/B2.  High nibble of byte 12 contains B3.
    # Low nibble of byte 12 contains sample interval I.
    b11 = h[10]
    b12 = h[11]

    b1 = _bcd_nibble((b11 >> 4) & 0x0F)
    b2 = _bcd_nibble(b11 & 0x0F)
    b3 = _bcd_nibble((b12 >> 4) & 0x0F)
    sample_interval_ms = _bcd_nibble(b12 & 0x0F)

    bytes_per_scan = 100 * b1 + 10 * b2 + b3

    manufacturer_code = decode_bcd_byte(h[12])  # byte 13
    equipment_serial = decode_bcd(h[13:16])     # bytes 14-16

    return SegBHeader(
        file_number=file_number,
        format_code=format_code,
        data_id=data_id,
        bytes_per_scan=bytes_per_scan,
        sample_interval_ms=sample_interval_ms,
        manufacturer_code=manufacturer_code,
        equipment_serial=equipment_serial,
        raw_standard=h,
    )


# -----------------------------------------------------------------------------
# Format-B geometry
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class FormatBGeometry:
    bytes_per_scan: int
    n_seismic_groups: int
    n_channels: int
    sync_bytes: int = 4
    auxiliary_bytes: int = 10
    seismic_group_bytes: int = 10

    @classmethod
    def from_bytes_per_scan(cls, bytes_per_scan: int) -> "FormatBGeometry":
        """
        Infer Format-B geometry.

        Format B:
            4-byte sync
          + 10-byte auxiliary group
          + N * 10-byte seismic groups

        Each seismic group contains four channels.
        """
        remainder = bytes_per_scan - 14

        if bytes_per_scan < 24 or remainder % 10 != 0:
            raise ValueError(
                f"{bytes_per_scan} bytes/scan is inconsistent with the "
                "Format-B structure 4 + 10 + N*10."
            )

        n_groups = remainder // 10

        return cls(
            bytes_per_scan=bytes_per_scan,
            n_seismic_groups=n_groups,
            n_channels=4 * n_groups,
        )

    @property
    def binary_gain_and_sync_bytes(self) -> int:
        # 4 sync bytes + one 2-byte gain word per seismic group.
        return 4 + 2 * self.n_seismic_groups

    @property
    def optional_J(self) -> int:
        # Format B definition from the standard:
        # J = bytes/scan - bytes(binary gain and sync) + 2
        return (
            self.bytes_per_scan
            - self.binary_gain_and_sync_bytes
            + 2
        )

    @property
    def extended_header_length(self) -> int:
        # If the optional portion is used, total header length is J + 24.
        return self.optional_J + 24


def parse_optional_header(
    buf: bytes,
    geometry: FormatBGeometry,
) -> SegBOptionalHeader:
    """
    Parse the defined portion of the optional header.

    `buf` must begin at header byte 17, i.e. immediately after the
    mandatory 16-byte standard header.

    The channel gain/identifier area begins at header byte 25. This
    first version preserves that area as raw bytes because interpretation
    depends on channel use and field-system gain configuration.
    """
    total_optional = geometry.extended_header_length - 16

    if len(buf) < total_optional:
        raise ValueError(
            f"Expected {total_optional} optional-header bytes, "
            f"received {len(buf)}"
        )

    h = buf[:total_optional]

    # Header byte 17
    record_length_seconds = decode_bcd_byte(h[0])

    # Header byte 18: high and low nibbles carry J and Z respectively.
    gain_control_mode = (h[1] >> 4) & 0x0F
    record_type = h[1] & 0x0F

    low_cut = decode_bcd_byte(h[2])  # byte 19

    # Byte 20: high nibble is slope code, low nibble reserved zero.
    low_cut_slope_code = (h[3] >> 4) & 0x0F

    # High-cut is three BCD digits across byte 21 and high nibble byte 22.
    hc12 = h[4]
    hc3_hs = h[5]
    hc1 = _bcd_nibble((hc12 >> 4) & 0x0F)
    hc2 = _bcd_nibble(hc12 & 0x0F)
    hc3 = _bcd_nibble((hc3_hs >> 4) & 0x0F)
    high_cut = 100 * hc1 + 10 * hc2 + hc3
    high_cut_slope_code = hc3_hs & 0x0F

    special_filter = decode_bcd_byte(h[6])  # byte 23

    # Byte 24: high nibble alias setting, low nibble reserved zero.
    alias_filter = (h[7] >> 4) & 0x0F

    # Header byte 25 onward.
    channel_gain_bytes = h[8:]

    return SegBOptionalHeader(
        record_length_seconds=record_length_seconds,
        gain_control_mode=gain_control_mode,
        record_type=record_type,
        low_cut=low_cut,
        low_cut_slope_code=low_cut_slope_code,
        high_cut=high_cut,
        high_cut_slope_code=high_cut_slope_code,
        special_filter=special_filter,
        alias_filter=alias_filter,
        channel_gain_bytes=channel_gain_bytes,
        raw=h,
    )


# -----------------------------------------------------------------------------
# Scan structures
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class FormatBSync:
    raw: bytes

    @property
    def sync_lsb_pattern(self) -> tuple[int, int, int, int]:
        """
        Return the least-significant bit of each sync byte.

        In the published diagram/text the sync discriminator is bit 7
        in each byte. With conventional Python byte numbering
        (MSB=bit 7, LSB=bit 0), the SEG document's bit 7 corresponds
        to the low-order bit because the document numbers bits 0..7
        from MSB to LSB.
        """
        return tuple(b & 0x01 for b in self.raw)

    @property
    def looks_valid(self) -> bool:
        return len(self.raw) == 4 and all(self.sync_lsb_pattern)


@dataclass(frozen=True)
class FormatBAuxiliary:
    raw: bytes
    words_u16: tuple[int, int, int, int, int]

    @property
    def time_counter_word(self) -> int:
        """Raw two-byte first auxiliary word."""
        return self.words_u16[0]

    @property
    def time_counter_ms(self) -> int:
        """Decoded 14-bit, 1-kHz Format-B time counter in milliseconds."""
        return decode_format_b_time_counter(self.raw[0], self.raw[1])

    @property
    def uphole_or_water_depth_word(self) -> int:
        return self.words_u16[1]

    @property
    def time_break_word(self) -> int:
        return self.words_u16[2]

    @property
    def timing_word(self) -> int:
        return self.words_u16[3]

    @property
    def special_word(self) -> int:
        return self.words_u16[4]


@dataclass(frozen=True)
class FormatBSeismicGroup:
    gain_word: int
    gain_codes: tuple[int, int, int, int]
    raw_words: tuple[int, int, int, int]
    samples: tuple[int, int, int, int]
    raw: bytes


@dataclass(frozen=True)
class FormatBScan:
    index: int
    sync: FormatBSync
    auxiliary: FormatBAuxiliary
    seismic_groups: tuple[FormatBSeismicGroup, ...]
    raw: bytes

    @property
    def samples(self) -> np.ndarray:
        return np.asarray(
            [x for group in self.seismic_groups for x in group.samples],
            dtype=np.int32,
        )

    @property
    def raw_words(self) -> np.ndarray:
        return np.asarray(
            [x for group in self.seismic_groups for x in group.raw_words],
            dtype=np.uint16,
        )

    @property
    def gain_codes(self) -> np.ndarray:
        return np.asarray(
            [x for group in self.seismic_groups for x in group.gain_codes],
            dtype=np.uint8,
        )


def decode_gain_word(word: int) -> tuple[int, int, int, int]:
    """Split a 16-bit Format-B gain word into four 4-bit gain codes."""
    return (
        (word >> 12) & 0x0F,
        (word >> 8) & 0x0F,
        (word >> 4) & 0x0F,
        word & 0x0F,
    )


def parse_format_b_scan(
    buf: bytes,
    geometry: FormatBGeometry,
    *,
    index: int = 0,
    negative_encoding: NegativeEncoding = "ones",
    validate_sync: bool = True,
) -> FormatBScan:
    """Parse one complete multiplexed Format-B scan."""
    n = geometry.bytes_per_scan

    if len(buf) < n:
        raise ValueError(
            f"Scan {index}: expected {n} bytes, received {len(buf)}"
        )

    raw = buf[:n]

    # 1. Sync group
    sync = FormatBSync(raw=raw[0:4])

    if validate_sync and not sync.looks_valid:
        raise ValueError(
            f"Scan {index}: invalid Format-B sync pattern: "
            f"{sync.raw.hex(' ')}"
        )

    # 2. Auxiliary group
    aux_raw = raw[4:14]
    aux_words = tuple(
        int.from_bytes(aux_raw[i:i + 2], "big")
        for i in range(0, 10, 2)
    )
    auxiliary = FormatBAuxiliary(
        raw=aux_raw,
        words_u16=aux_words,
    )

    # 3. Seismic groups
    groups = []

    for igroup in range(geometry.n_seismic_groups):
        i0 = 14 + 10 * igroup
        graw = raw[i0:i0 + 10]

        gain_word = int.from_bytes(graw[0:2], "big")
        gain_codes = decode_gain_word(gain_word)

        raw_words = tuple(
            int.from_bytes(graw[i:i + 2], "big")
            for i in (2, 4, 6, 8)
        )

        samples = tuple(
            decode_sample_word(word, negative_encoding)
            for word in raw_words
        )

        groups.append(
            FormatBSeismicGroup(
                gain_word=gain_word,
                gain_codes=gain_codes,
                raw_words=raw_words,
                samples=samples,
                raw=graw,
            )
        )

    return FormatBScan(
        index=index,
        sync=sync,
        auxiliary=auxiliary,
        seismic_groups=tuple(groups),
        raw=raw,
    )


# -----------------------------------------------------------------------------
# Format-B stream detection
# -----------------------------------------------------------------------------



@dataclass(frozen=True)
class FormatBStreamAudit:
    """QC summary for a contiguous run of candidate Format-B scans."""
    n_scans_checked: int
    n_valid_prefixes: int
    first_invalid_scan: Optional[int]
    invalid_scan_indices: tuple[int, ...]
    counter_discontinuity_indices: tuple[int, ...]
    counter_start: Optional[int]
    counter_end: Optional[int]
    expected_counter_step_ms: Optional[int]

    @property
    def all_prefixes_valid(self) -> bool:
        return self.n_valid_prefixes == self.n_scans_checked

    @property
    def counter_is_continuous(self) -> bool:
        return len(self.counter_discontinuity_indices) == 0


@dataclass(frozen=True)
class FormatBDetection:
    """Geometry inferred directly from a run of Format-B scans."""
    data_offset: int
    geometry: FormatBGeometry
    n_confirmed_scans: int
    sample_interval_ms: Optional[int]


def _aux_discriminator_valid(aux: bytes) -> bool:
    """
    Check the Format-B auxiliary anti-sync discriminator.

    The 1967 document numbers bits 0..7 from MSB to LSB.  It requires
    document bit 7 (the conventional byte LSB) to be zero in auxiliary
    byte positions 1, 2, 4, 6, 8 and 10.
    """
    if len(aux) < 10:
        return False
    return all((aux[i] & 0x01) == 0 for i in (0, 1, 3, 5, 7, 9))


def _scan_prefix_valid(buf: bytes, offset: int) -> bool:
    """Validate the 4-byte sync plus the 10-byte auxiliary discriminator."""
    if offset < 0 or offset + 14 > len(buf):
        return False
    return (
        FormatBSync(buf[offset:offset + 4]).looks_valid
        and _aux_discriminator_valid(buf[offset + 4:offset + 14])
    )


def detect_format_b_stream(
    buf: bytes,
    *,
    search_start: int = 16,
    search_stop: Optional[int] = None,
    min_scans: int = 4,
    max_bytes_per_scan: int = 4096,
) -> FormatBDetection:
    """
    Infer first scan offset and bytes/scan from repeated Format-B framing.

    Detection uses both constraints defined by the standard:
      * sync bit set in all four sync bytes; and
      * the auxiliary anti-sync bit clear in byte positions 1,2,4,6,8,10.

    A candidate scan spacing must also satisfy 4 + 10 + N*10 bytes.
    This is deliberately independent of the BCD bytes/scan header field,
    because surviving legacy files can contain implementation-specific
    headers.
    """
    if min_scans < 2:
        raise ValueError("min_scans must be at least 2")

    n = len(buf)
    if search_stop is None:
        search_stop = min(n, max(search_start + 65536, 65536))
    search_stop = min(search_stop, n - 14)

    candidates = [
        i for i in range(max(0, search_start), max(0, search_stop) + 1)
        if _scan_prefix_valid(buf, i)
    ]

    best = None
    for ia, start in enumerate(candidates):
        for nxt in candidates[ia + 1:]:
            spacing = nxt - start
            if spacing > max_bytes_per_scan:
                break
            try:
                geometry = FormatBGeometry.from_bytes_per_scan(spacing)
            except ValueError:
                continue

            count = 1
            pos = start + spacing
            while pos + 14 <= n and _scan_prefix_valid(buf, pos):
                count += 1
                pos += spacing

            if count < min_scans:
                continue

            score = (count, -start)
            if best is None or score > best[0]:
                best = (score, start, geometry, count)

    if best is None:
        raise ValueError(
            "Could not detect a repeating SEG-B Format-B scan structure. "
            "Specify header_length and/or bytes_per_scan explicitly."
        )

    _, start, geometry, count = best

    # The first auxiliary word is normally a 1-kHz time counter.  Estimate
    # dt only when several consecutive increments agree and are plausible.
    counters = []
    for k in range(min(count, 16)):
        o = start + k * geometry.bytes_per_scan + 4
        counters.append(decode_format_b_time_counter(buf[o], buf[o + 1]))

    dt_ms = None
    if len(counters) >= 3:
        diffs = [((b - a) & 0x3FFF) for a, b in zip(counters, counters[1:])]
        if diffs and len(set(diffs)) == 1 and 0 < diffs[0] <= 1000:
            dt_ms = diffs[0]

    return FormatBDetection(
        data_offset=start,
        geometry=geometry,
        n_confirmed_scans=count,
        sample_interval_ms=dt_ms,
    )


# -----------------------------------------------------------------------------
# Record / file reader
# -----------------------------------------------------------------------------


@dataclass
class SegBRecord:
    """One logical SEG-B seismic record inside a disk file."""
    parent: "SegBFile"
    index: int
    header_offset: int
    data_offset: int
    n_samples: int
    geometry: FormatBGeometry
    header: SegBHeader
    sample_interval_ms: int

    @property
    def n_channels(self) -> int:
        return self.geometry.n_channels

    @property
    def record_length_seconds(self) -> float:
        return self.n_samples * self.dt_seconds

    @property
    def raw_preamble(self) -> bytes:
        """Complete bytes preceding this record's first data scan."""
        return self.parent._bytes[self.header_offset:self.data_offset]

    @property
    def raw_header(self) -> bytes:
        """Alias for raw_preamble; retained as a convenient user-facing name."""
        return self.raw_preamble

    @property
    def data_id(self) -> str:
        return self.header.data_id

    @property
    def vendor_extension(self) -> bytes:
        """Recorder-specific portion of the preamble not interpreted by this module.

        The 1967 SEG allocation is not sufficient to identify the observed 420-byte
        preamble unambiguously, so bytes beyond the 16-byte standard header are
        intentionally preserved rather than assigned speculative meanings.
        """
        return self.raw_preamble[16:]

    @property
    def metadata(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "data_id": self.data_id,
            "n_channels": self.n_channels,
            "n_samples": self.n_samples,
            "sample_interval_ms": self.sample_interval_ms,
            "record_length_seconds": self.record_length_seconds,
            "bytes_per_scan": self.bytes_per_scan,
            "header_offset": self.header_offset,
            "data_offset": self.data_offset,
            "preamble_bytes": self.data_offset - self.header_offset,
        }

    @property
    def bytes_per_scan(self) -> int:
        return self.geometry.bytes_per_scan

    @property
    def dt_seconds(self) -> float:
        return self.sample_interval_ms / 1000.0

    @property
    def time_seconds(self) -> np.ndarray:
        return np.arange(self.n_samples, dtype=np.float64) * self.dt_seconds

    @property
    def end_offset(self) -> int:
        return self.data_offset + self.n_samples * self.bytes_per_scan

    def read_scan(self, index: int) -> FormatBScan:
        if index < 0:
            index += self.n_samples
        if not 0 <= index < self.n_samples:
            raise IndexError(f"scan index {index} outside [0, {self.n_samples})")
        n = self.bytes_per_scan
        i0 = self.data_offset + index * n
        return parse_format_b_scan(
            self.parent._bytes[i0:i0+n], self.geometry, index=index,
            negative_encoding=self.parent.negative_encoding,
            validate_sync=self.parent.validate_sync,
        )

    def read(self, *, gain_correct: bool = False, normalized: bool = False,
             scale_constant: float = 1.0, dtype=np.float64) -> np.ndarray:
        """Read seismic channels as ``(n_channels, n_samples)``.

        With ``gain_correct=False`` the returned values are the decoded
        15-bit signed mantissa integers.  ``normalized=True`` divides the
        mantissa by 2**14, corresponding to the SEG-B Qs.Q1..Q14 binary
        fraction.  ``gain_correct=True`` additionally applies 2**(-GGGG).
        ``scale_constant`` supplies K in the published SEG-B equation; its
        physical value is recorder-dependent and is therefore 1.0 by default.
        """
        floating = gain_correct or normalized or scale_constant != 1.0
        out = np.empty((self.n_channels, self.n_samples),
                       dtype=np.float64 if floating else np.int32)
        for i in range(self.n_samples):
            scan = self.read_scan(i)
            v = scan.samples.astype(np.float64) if floating else scan.samples
            if normalized or gain_correct:
                v = v / float(1 << 14)
            if gain_correct:
                v = v * np.exp2(-scan.gain_codes.astype(np.float64))
            if scale_constant != 1.0:
                v = v * float(scale_constant)
            out[:, i] = v
        return out.astype(dtype, copy=False) if floating else out

    def read_gain_codes(self) -> np.ndarray:
        out = np.empty((self.n_channels, self.n_samples), dtype=np.uint8)
        for i in range(self.n_samples):
            out[:, i] = self.read_scan(i).gain_codes
        return out

    def read_auxiliary_words(self) -> np.ndarray:
        out = np.empty((5, self.n_samples), dtype=np.uint16)
        for i in range(self.n_samples):
            out[:, i] = self.read_scan(i).auxiliary.words_u16
        return out

    def read_time_counter(self) -> np.ndarray:
        """Return the decoded 14-bit auxiliary time counter in milliseconds."""
        out = np.empty(self.n_samples, dtype=np.uint16)
        for i in range(self.n_samples):
            out[i] = self.read_scan(i).auxiliary.time_counter_ms
        return out

    def audit_time_counter(self) -> tuple[int, ...]:
        """Return scan indices whose counter increment differs from dt.

        Counter differences are evaluated modulo 2**14 so a normal 14-bit
        wrap is not flagged.  Index ``i`` identifies the second sample in the
        discontinuous pair (i-1, i).
        """
        if self.sample_interval_ms <= 0:
            raise ValueError("Sample interval is unknown")
        t = self.read_time_counter().astype(np.int64)
        if t.size < 2:
            return ()
        d = (np.diff(t) % (1 << 14))
        return tuple((np.flatnonzero(d != self.sample_interval_ms) + 1).tolist())

    def validate(self) -> dict[str, Any]:
        """Validate scan framing and auxiliary timing for this record.

        Returns a dictionary suitable for QC/logging. ``valid`` is True only
        when every scan has a valid Format-B prefix and the decoded 14-bit
        auxiliary time counter advances by the record sample interval.
        """
        bad_sync = []
        for i in range(self.n_samples):
            off = self.data_offset + i * self.bytes_per_scan
            if not _scan_prefix_valid(self.parent._bytes, off):
                bad_sync.append(i)
        bad_time = self.audit_time_counter()
        return {
            "record": self.index,
            "n_samples": self.n_samples,
            "bytes_per_scan": self.bytes_per_scan,
            "invalid_scan_indices": tuple(bad_sync),
            "counter_discontinuity_indices": bad_time,
            "valid": not bad_sync and not bad_time,
        }

    def get_trace(self, channel: int, *, gain_correct: bool = False, one_based: bool = True) -> np.ndarray:
        ich = channel - 1 if one_based else channel
        if not 0 <= ich < self.n_channels:
            raise IndexError(f"channel {channel} outside available channel range")
        return self.read(gain_correct=gain_correct)[ich]

    def __repr__(self) -> str:
        return (f"SegBRecord(index={self.index}, channels={self.n_channels}, "
                f"samples={self.n_samples}, dt={self.sample_interval_ms} ms, "
                f"header_offset={self.header_offset}, data_offset={self.data_offset})")

class SegBFile:
    """Reader for a SEG-B Format-B file or contiguous record payload.

    In ``header_mode='auto'`` the reader first parses the 16-byte header for
    metadata, but it does not blindly trust the header bytes/scan field.  It
    detects repeated Format-B scan framing in the data stream and uses that
    geometry when the two disagree.
    """

    def __init__(
        self,
        filename: PathLike,
        *,
        negative_encoding: NegativeEncoding = "ones",
        header_length: Optional[int] = None,
        header_mode: Literal["auto", "standard", "extended"] = "auto",
        validate_sync: bool = True,
        bytes_per_scan: Optional[int] = None,
        sample_interval_ms: Optional[int] = None,
    ):
        self.filename = Path(filename)
        self.negative_encoding = negative_encoding
        self.validate_sync = validate_sync
        self._bytes = self.filename.read_bytes()

        if negative_encoding not in ("ones", "twos"):
            raise ValueError("negative_encoding must be 'ones' or 'twos'")
        if header_mode not in ("auto", "standard", "extended"):
            raise ValueError("header_mode must be 'auto', 'standard', or 'extended'")
        if len(self._bytes) < 16:
            raise ValueError("File is too short to contain a SEG-B header")

        # Header parsing is retained because its identification fields are
        # useful, but geometry can be recovered independently from the scans.
        self.header = parse_standard_header(self._bytes[:16])
        self.detection: Optional[FormatBDetection] = None
        self.geometry_source = "header"

        header_geometry = None
        try:
            header_geometry = FormatBGeometry.from_bytes_per_scan(
                self.header.bytes_per_scan
            )
        except ValueError:
            pass

        if bytes_per_scan is not None:
            self.geometry = FormatBGeometry.from_bytes_per_scan(int(bytes_per_scan))
            self.geometry_source = "explicit"
        elif header_mode == "auto":
            try:
                self.detection = detect_format_b_stream(self._bytes)
                self.geometry = self.detection.geometry
                self.geometry_source = "detected"
            except ValueError:
                if header_geometry is None:
                    raise ValueError(
                        f"Header reports {self.header.bytes_per_scan} bytes/scan, "
                        "which is not valid Format-B geometry, and scan geometry "
                        "could not be detected from the data stream."
                    )
                self.geometry = header_geometry
        else:
            if header_geometry is None and bytes_per_scan is None:
                raise ValueError(
                    f"Header reports {self.header.bytes_per_scan} bytes/scan, "
                    "which is inconsistent with Format B. Use header_mode='auto' "
                    "or specify bytes_per_scan explicitly."
                )
            self.geometry = header_geometry

        if header_length is not None:
            self.header_length = int(header_length)
        elif header_mode == "auto" and self.detection is not None:
            self.header_length = self.detection.data_offset
        elif header_mode == "standard":
            self.header_length = 16
        elif header_mode == "extended":
            self.header_length = self.geometry.extended_header_length
        else:
            self.header_length = self._infer_header_length()

        if not 0 <= self.header_length < len(self._bytes):
            raise ValueError("header_length/data offset is outside the file")

        self.optional_header: Optional[SegBOptionalHeader] = None
        expected = self.geometry.extended_header_length
        if self.header_length == expected and self.header_length > 16:
            try:
                self.optional_header = parse_optional_header(
                    self._bytes[16:self.header_length], self.geometry
                )
            except ValueError:
                # Preserve access to the data even when a legacy optional
                # header does not follow the generic allocation exactly.
                self.optional_header = None

        if sample_interval_ms is not None:
            self.sample_interval_ms = int(sample_interval_ms)
            self.sample_interval_source = "explicit"
        elif self.detection is not None and self.detection.sample_interval_ms is not None:
            self.sample_interval_ms = self.detection.sample_interval_ms
            self.sample_interval_source = "auxiliary_time_counter"
        elif self.header.sample_interval_ms > 0:
            self.sample_interval_ms = self.header.sample_interval_ms
            self.sample_interval_source = "header"
        else:
            self.sample_interval_ms = 0
            self.sample_interval_source = "unknown"

        data_bytes = len(self._bytes) - self.header_length
        self.n_complete_scans = data_bytes // self.geometry.bytes_per_scan
        self.trailing_bytes = data_bytes % self.geometry.bytes_per_scan

        # Discover logical records.  A disk image may contain many SEG-B
        # records concatenated with their per-record headers and no tape gaps.
        self.records = self._discover_records()

    def _count_scan_run(self, data_offset: int, geometry: FormatBGeometry) -> int:
        nscan = 0
        step = geometry.bytes_per_scan
        pos = data_offset
        while pos + step <= len(self._bytes) and _scan_prefix_valid(self._bytes, pos):
            nscan += 1
            pos += step
        return nscan

    def _discover_records(self) -> list[SegBRecord]:
        records: list[SegBRecord] = []
        geometry = self.geometry
        data_offset = self.header_length
        header_offset = 0
        idx = 0

        while data_offset + geometry.bytes_per_scan <= len(self._bytes):
            nscan = self._count_scan_run(data_offset, geometry)
            if nscan == 0:
                break

            try:
                hdr = parse_standard_header(self._bytes[header_offset:header_offset+16])
            except Exception:
                hdr = self.header if idx == 0 else self.header

            records.append(SegBRecord(
                parent=self, index=idx, header_offset=header_offset,
                data_offset=data_offset, n_samples=nscan, geometry=geometry,
                header=hdr, sample_interval_ms=self.sample_interval_ms,
            ))

            boundary = data_offset + nscan * geometry.bytes_per_scan
            if boundary + 16 >= len(self._bytes):
                break

            # In the observed disk representation, the next header starts
            # exactly at the previous data end.  Do not assume its length;
            # find the next repeating Format-B stream nearby.
            try:
                det = detect_format_b_stream(
                    self._bytes, search_start=boundary + 16,
                    search_stop=min(len(self._bytes)-14, boundary + 8192),
                    min_scans=3, max_bytes_per_scan=4096,
                )
            except ValueError:
                break

            if det.geometry.bytes_per_scan != geometry.bytes_per_scan:
                geometry = det.geometry
            header_offset = boundary
            data_offset = det.data_offset
            idx += 1

        return records

    @property
    def n_records(self) -> int:
        return len(self.records)

    def record(self, index: int) -> SegBRecord:
        return self.records[index]

    def read_record(self, index: int, **kwargs) -> np.ndarray:
        """Convenience wrapper for ``record(index).read(**kwargs)``."""
        return self.record(index).read(**kwargs)

    def read_trace(self, record: int, channel: int, *, one_based: bool = True,
                   **kwargs) -> np.ndarray:
        """Read one channel from one logical record."""
        rec = self.record(record)
        ich = channel - 1 if one_based else channel
        if not 0 <= ich < rec.n_channels:
            raise IndexError(f"channel {channel} outside available channel range")
        return rec.read(**kwargs)[ich]

    def summary(self, *, max_records: Optional[int] = None, print_output: bool = True) -> str:
        """Return a compact table describing the logical records in the file."""
        rows = self.records if max_records is None else self.records[:int(max_records)]
        lines = [
            "SEG-B file",
            "-" * 78,
            f"File: {self.filename}",
            f"Records: {self.n_records}",
            f"Seismic channels: {self.n_channels}",
            f"Sample interval: {self.sample_interval_ms} ms",
            f"Bytes/scan: {self.geometry.bytes_per_scan}",
            "",
            "Index  Data ID       Samples  Length(s)  Header offset  Data offset",
            "-----  ------------  -------  ---------  -------------  -----------",
        ]
        for r in rows:
            lines.append(
                f"{r.index:5d}  {r.data_id:12s}  {r.n_samples:7d}  "
                f"{r.record_length_seconds:9.3f}  {r.header_offset:13d}  {r.data_offset:11d}"
            )
        if max_records is not None and self.n_records > len(rows):
            lines.append(f"... {self.n_records - len(rows)} more records")
        text = "\n".join(lines)
        if print_output:
            print(text)
        return text

    def validate(self) -> dict[str, Any]:
        """Validate every discovered logical record in the file.

        Checks Format-B scan framing and the decoded 14-bit auxiliary time
        counter independently for each record.
        """
        records = [r.validate() for r in self.records]
        return {
            "n_records": self.n_records,
            "n_valid_records": sum(int(r["valid"]) for r in records),
            "valid": all(r["valid"] for r in records),
            "records": records,
        }

    def _sync_at(self, offset: int) -> bool:
        return _scan_prefix_valid(self._bytes, offset)

    def _infer_header_length(self) -> int:
        candidates = [16, self.geometry.extended_header_length]
        good = []
        for offset in candidates:
            if not self._sync_at(offset):
                continue
            nxt = offset + self.geometry.bytes_per_scan
            if nxt + 14 <= len(self._bytes) and not self._sync_at(nxt):
                continue
            good.append(offset)
        if good:
            return max(good)
        raise ValueError(
            "Could not infer the first Format-B scan offset. Specify "
            "header_length explicitly or use header_mode='auto'."
        )

    def __repr__(self) -> str:
        return (
            f"SegBFile({str(self.filename)!r}, channels={self.n_channels}, "
            f"dt={self.sample_interval_ms} ms, records={self.n_records}, "
            f"bytes_per_scan={self.geometry.bytes_per_scan}, "
            f"data_offset={self.header_length}, geometry={self.geometry_source!r})"
        )

    @property
    def n_channels(self) -> int:
        return self.geometry.n_channels

    @property
    def n_samples(self) -> int:
        return self.n_complete_scans

    @property
    def dt_seconds(self) -> float:
        return self.sample_interval_ms / 1000.0

    @property
    def time_seconds(self) -> np.ndarray:
        if self.sample_interval_ms <= 0:
            raise ValueError(
                "Sample interval is unknown. Supply sample_interval_ms explicitly."
            )
        return np.arange(self.n_samples, dtype=np.float64) * self.dt_seconds


    def audit_stream(self, *, max_scans: Optional[int] = None,
                     max_reported_errors: int = 32) -> FormatBStreamAudit:
        """
        Check scan framing and the 14-bit auxiliary time counter.

        This is useful before treating every complete 314-byte block (or other
        detected geometry) as one continuous seismic record.  Counter
        differences are evaluated modulo 2**14, so normal 14-bit wraparound
        is not reported as a discontinuity.
        """
        ncheck = self.n_complete_scans
        if max_scans is not None:
            ncheck = min(ncheck, int(max_scans))

        invalid = []
        nvalid = 0
        discontinuities = []
        counters = []
        step = self.sample_interval_ms if self.sample_interval_ms > 0 else None
        n = self.geometry.bytes_per_scan

        for i in range(ncheck):
            off = self.header_length + i * n
            valid = _scan_prefix_valid(self._bytes, off)
            if valid:
                nvalid += 1
            elif len(invalid) < max_reported_errors:
                invalid.append(i)

            if off + 6 <= len(self._bytes):
                counter = decode_format_b_time_counter(self._bytes[off + 4], self._bytes[off + 5])
                counters.append(counter)
                if step is not None and len(counters) >= 2:
                    diff = (counters[-1] - counters[-2]) & 0x3FFF
                    if diff != step and len(discontinuities) < max_reported_errors:
                        discontinuities.append(i)

        return FormatBStreamAudit(
            n_scans_checked=ncheck,
            n_valid_prefixes=nvalid,
            first_invalid_scan=invalid[0] if invalid else None,
            invalid_scan_indices=tuple(invalid),
            counter_discontinuity_indices=tuple(discontinuities),
            counter_start=counters[0] if counters else None,
            counter_end=counters[-1] if counters else None,
            expected_counter_step_ms=step,
        )

    def read_scan(self, index: int) -> FormatBScan:
        if index < 0:
            index += self.n_complete_scans
        if not 0 <= index < self.n_complete_scans:
            raise IndexError(f"scan index {index} outside [0, {self.n_complete_scans})")
        n = self.geometry.bytes_per_scan
        i0 = self.header_length + index * n
        return parse_format_b_scan(
            self._bytes[i0:i0 + n], self.geometry, index=index,
            negative_encoding=self.negative_encoding,
            validate_sync=self.validate_sync,
        )

    def read(self, *, gain_correct: bool = False, dtype=np.float64) -> np.ndarray:
        data = np.empty(
            (self.n_channels, self.n_samples),
            dtype=np.float64 if gain_correct else np.int32,
        )
        for iscan in range(self.n_samples):
            scan = self.read_scan(iscan)
            values = scan.samples
            if gain_correct:
                values = values.astype(np.float64)
                values *= np.exp2(-scan.gain_codes.astype(np.float64))
            data[:, iscan] = values
        return data.astype(dtype, copy=False) if gain_correct else data

    def read_raw_words(self) -> np.ndarray:
        data = np.empty((self.n_channels, self.n_samples), dtype=np.uint16)
        for iscan in range(self.n_samples):
            data[:, iscan] = self.read_scan(iscan).raw_words
        return data

    def read_gain_codes(self) -> np.ndarray:
        gain = np.empty((self.n_channels, self.n_samples), dtype=np.uint8)
        for iscan in range(self.n_samples):
            gain[:, iscan] = self.read_scan(iscan).gain_codes
        return gain

    def read_auxiliary_words(self) -> np.ndarray:
        aux = np.empty((5, self.n_samples), dtype=np.uint16)
        for iscan in range(self.n_samples):
            aux[:, iscan] = self.read_scan(iscan).auxiliary.words_u16
        return aux

    def read_time_counter(self) -> np.ndarray:
        """Return decoded time-counter milliseconds for the legacy first stream."""
        out = np.empty(self.n_samples, dtype=np.uint16)
        for iscan in range(self.n_samples):
            out[iscan] = self.read_scan(iscan).auxiliary.time_counter_ms
        return out

    def get_trace(self, channel: int, *, gain_correct: bool = False,
                  one_based: bool = True) -> np.ndarray:
        ich = channel - 1 if one_based else channel
        if not 0 <= ich < self.n_channels:
            raise IndexError(f"channel {channel} outside available channel range")
        return self.read(gain_correct=gain_correct)[ich]


# -----------------------------------------------------------------------------
# Convenience function
# -----------------------------------------------------------------------------

def read_segb(
    filename: PathLike,
    *,
    negative_encoding: NegativeEncoding = "ones",
    gain_correct: bool = False,
    header_length: Optional[int] = None,
    header_mode: Literal["auto", "standard", "extended"] = "auto",
    validate_sync: bool = True,
    bytes_per_scan: Optional[int] = None,
    sample_interval_ms: Optional[int] = None,
) -> tuple[np.ndarray, SegBFile]:
    """
    Convenience reader.

    Returns
    -------
    data, segb
        `data` has shape (n_channels, n_samples).
        `segb` provides header, geometry, gains and auxiliary information.
    """
    segb = SegBFile(
        filename,
        negative_encoding=negative_encoding,
        header_length=header_length,
        header_mode=header_mode,
        validate_sync=validate_sync,
        bytes_per_scan=bytes_per_scan,
        sample_interval_ms=sample_interval_ms,
    )
    return segb.read(gain_correct=gain_correct), segb


__all__ = [
    "SegBFile",
    "SegBRecord",
    "SegBHeader",
    "SegBOptionalHeader",
    "FormatBGeometry",
    "FormatBSync",
    "FormatBAuxiliary",
    "FormatBSeismicGroup",
    "FormatBScan",
    "FormatBDetection",
    "FormatBStreamAudit",
    "detect_format_b_stream",
    "decode_bcd_byte",
    "decode_bcd",
    "decode_bcd_digits",
    "decode_gain_word",
    "decode_ones_complement_u16",
    "decode_twos_complement_u16",
    "decode_signed_nbit",
    "sample_word_to_fraction",
    "decode_sample_word",
    "parse_standard_header",
    "parse_optional_header",
    "parse_format_b_scan",
    "read_segb",
]
