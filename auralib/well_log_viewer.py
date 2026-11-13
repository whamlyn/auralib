"""Matplotlib-based interactive well-log viewer.

Dependencies: numpy, matplotlib
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.widgets import RangeSlider

ArrayLike = Union[np.ndarray, Sequence[float]]


@dataclass
class LogStyle:
    color: str = "black"
    linestyle: str = "-"
    linewidth: float = 1.0
    alpha: float = 1.0
    label: Optional[str] = None
    zorder: float = 3.0


@dataclass
class TrackSpec:
    title: str = ""
    xscale: str = "linear"          # "linear" or "log"
    xlim: Optional[Tuple[float, float]] = None
    width_ratio: float = 1.0
    grid: bool = True


@dataclass
class LogCurve:
    name: str
    values: np.ndarray
    depth: np.ndarray
    track: int
    style: LogStyle = field(default_factory=LogStyle)


@dataclass
class FillSpec:
    kind: str                        # solid_left, solid_right, solid_between, cmap_left, cmap_right, cmap_between
    track: int
    curve1: str
    curve2: Optional[str] = None
    baseline: Optional[float] = None
    color: str = "0.8"
    alpha: float = 0.65
    cmap: str = "viridis"
    value_curve: Optional[str] = None
    vmin: Optional[float] = None
    vmax: Optional[float] = None
    log_norm: bool = False
    zorder: float = 1.0


class WellLogViewer:
    """Interactive multi-track well-log display.

    Depth is shared across all tracks and increases downward. A vertical
    RangeSlider selects the visible depth interval. Mouse wheel zooms around
    the cursor; right-button drag pans vertically.
    """

    def __init__(
        self,
        n_tracks: int,
        track_specs: Optional[Sequence[TrackSpec]] = None,
        figsize: Tuple[float, float] = (12.0, 8.0),
        depth_label: str = "Depth",
    ) -> None:
        if n_tracks < 1:
            raise ValueError("n_tracks must be >= 1")
        self.n_tracks = int(n_tracks)
        self.track_specs = list(track_specs or [TrackSpec() for _ in range(n_tracks)])
        if len(self.track_specs) != self.n_tracks:
            raise ValueError("track_specs length must equal n_tracks")
        self.figsize = figsize
        self.depth_label = depth_label
        self.curves: Dict[str, LogCurve] = {}
        self.fills: List[FillSpec] = []
        self.fig = None
        self.axes = []
        self.depth_slider = None
        self._depth_min = None
        self._depth_max = None
        self._drag_y = None
        self._updating_slider = False

    @staticmethod
    def _array(a: ArrayLike, name: str) -> np.ndarray:
        out = np.asarray(a, dtype=float).reshape(-1)
        if out.size == 0:
            raise ValueError(f"{name} cannot be empty")
        return out

    def add_log(
        self,
        name: str,
        values: ArrayLike,
        depth: ArrayLike,
        track: int,
        *,
        color: str = "black",
        linestyle: str = "-",
        linewidth: float = 1.0,
        alpha: float = 1.0,
        label: Optional[str] = None,
        zorder: float = 3.0,
    ) -> "WellLogViewer":
        if not 0 <= track < self.n_tracks:
            raise IndexError("track is zero-based and must be within n_tracks")
        v = self._array(values, "values")
        d = self._array(depth, "depth")
        if v.size != d.size:
            raise ValueError("values and depth must have the same length")
        order = np.argsort(d)
        self.curves[name] = LogCurve(
            name=name, values=v[order], depth=d[order], track=track,
            style=LogStyle(color, linestyle, linewidth, alpha, label or name, zorder),
        )
        return self

    def add_solid_fill(
        self,
        track: int,
        curve1: str,
        *,
        side: str = "left",
        curve2: Optional[str] = None,
        baseline: Optional[float] = None,
        color: str = "0.8",
        alpha: float = 0.65,
        zorder: float = 1.0,
    ) -> "WellLogViewer":
        kind = "solid_between" if curve2 else f"solid_{side}"
        if kind not in {"solid_left", "solid_right", "solid_between"}:
            raise ValueError("side must be 'left' or 'right'")
        self.fills.append(FillSpec(kind, track, curve1, curve2, baseline, color, alpha, zorder=zorder))
        return self

    def add_colormapped_fill(
        self,
        track: int,
        curve1: str,
        *,
        side: str = "left",
        curve2: Optional[str] = None,
        baseline: Optional[float] = None,
        value_curve: Optional[str] = None,
        cmap: str = "viridis",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        log_norm: bool = False,
        alpha: float = 0.85,
        zorder: float = 1.0,
    ) -> "WellLogViewer":
        kind = "cmap_between" if curve2 else f"cmap_{side}"
        if kind not in {"cmap_left", "cmap_right", "cmap_between"}:
            raise ValueError("side must be 'left' or 'right'")
        self.fills.append(FillSpec(
            kind=kind, track=track, curve1=curve1, curve2=curve2,
            baseline=baseline, alpha=alpha, cmap=cmap,
            value_curve=value_curve or curve1, vmin=vmin, vmax=vmax,
            log_norm=log_norm, zorder=zorder,
        ))
        return self

    def _validate(self) -> None:
        if not self.curves:
            raise RuntimeError("Add at least one log before show()")
        for f in self.fills:
            for name in (f.curve1, f.curve2, f.value_curve):
                if name is not None and name not in self.curves:
                    raise KeyError(f"Unknown curve in fill: {name}")
            if self.curves[f.curve1].track != f.track:
                raise ValueError(f"{f.curve1} is not in fill track {f.track}")
            if f.curve2 and self.curves[f.curve2].track != f.track:
                raise ValueError(f"{f.curve2} is not in fill track {f.track}")

    @staticmethod
    def _interp_to(depth_target: np.ndarray, curve: LogCurve) -> np.ndarray:
        good = np.isfinite(curve.depth) & np.isfinite(curve.values)
        if good.sum() < 2:
            return np.full(depth_target.shape, np.nan)
        return np.interp(depth_target, curve.depth[good], curve.values[good], left=np.nan, right=np.nan)

    def _track_edge(self, track: int, side: str) -> float:
        spec = self.track_specs[track]
        if spec.xlim is not None:
            return spec.xlim[0 if side == "left" else 1]
        vals = np.concatenate([c.values[np.isfinite(c.values)] for c in self.curves.values() if c.track == track])
        return float(np.nanmin(vals) if side == "left" else np.nanmax(vals))

    def _draw_cmap_fill(self, ax, f: FillSpec) -> None:
        c1 = self.curves[f.curve1]
        d = c1.depth
        x1 = c1.values
        if f.curve2:
            x2 = self._interp_to(d, self.curves[f.curve2])
        else:
            side = "left" if f.kind.endswith("left") else "right"
            edge = f.baseline if f.baseline is not None else self._track_edge(f.track, side)
            x2 = np.full_like(x1, edge)
        z = self._interp_to(d, self.curves[f.value_curve])
        if d.size < 2:
            return
        de = np.empty(d.size + 1)
        de[1:-1] = 0.5 * (d[:-1] + d[1:])
        de[0] = d[0] - 0.5 * (d[1] - d[0])
        de[-1] = d[-1] + 0.5 * (d[-1] - d[-2])
        X = np.column_stack((x1, x2))
        Xe = np.empty((d.size + 1, 2))
        Xe[1:-1] = 0.5 * (X[:-1] + X[1:])
        Xe[0], Xe[-1] = X[0], X[-1]
        Y = np.column_stack((de, de))
        C = z[:, None]
        finite = z[np.isfinite(z)]
        if finite.size == 0:
            return
        vmin = f.vmin if f.vmin is not None else float(np.nanmin(finite))
        vmax = f.vmax if f.vmax is not None else float(np.nanmax(finite))
        norm = LogNorm(vmin=max(vmin, np.finfo(float).tiny), vmax=vmax) if f.log_norm else Normalize(vmin, vmax)
        ax.pcolormesh(Xe, Y, C, cmap=f.cmap, norm=norm, shading="flat", alpha=f.alpha, zorder=f.zorder)

    def _draw(self) -> None:
        self._validate()
        all_depth = np.concatenate([c.depth[np.isfinite(c.depth)] for c in self.curves.values()])
        self._depth_min, self._depth_max = float(all_depth.min()), float(all_depth.max())
        width_ratios = [s.width_ratio for s in self.track_specs]
        self.fig, axs = plt.subplots(
            1, self.n_tracks, sharey=True, figsize=self.figsize,
            gridspec_kw={"width_ratios": width_ratios}, squeeze=False,
        )
        self.axes = list(axs[0])
        self.fig.subplots_adjust(left=0.08, right=0.91, bottom=0.08, top=0.90, wspace=0.08)

        for i, (ax, spec) in enumerate(zip(self.axes, self.track_specs)):
            ax.set_xscale(spec.xscale)
            if spec.xlim is not None:
                ax.set_xlim(*spec.xlim)
            ax.set_title(spec.title or f"Track {i + 1}")
            ax.xaxis.set_label_position("top")
            ax.xaxis.tick_top()
            ax.grid(spec.grid, which="both", color="0.85", linewidth=0.6)
            if i == 0:
                ax.set_ylabel(self.depth_label)
            else:
                ax.tick_params(labelleft=False)

        for f in self.fills:
            ax = self.axes[f.track]
            c1 = self.curves[f.curve1]
            if f.kind.startswith("solid"):
                if f.curve2:
                    x2 = self._interp_to(c1.depth, self.curves[f.curve2])
                else:
                    side = "left" if f.kind.endswith("left") else "right"
                    edge = f.baseline if f.baseline is not None else self._track_edge(f.track, side)
                    x2 = np.full_like(c1.values, edge)
                mask = np.isfinite(c1.depth) & np.isfinite(c1.values) & np.isfinite(x2)
                ax.fill_betweenx(c1.depth, c1.values, x2, where=mask, color=f.color,
                                 alpha=f.alpha, interpolate=True, zorder=f.zorder)
            else:
                self._draw_cmap_fill(ax, f)

        for curve in self.curves.values():
            s = curve.style
            self.axes[curve.track].plot(curve.values, curve.depth, color=s.color,
                                        linestyle=s.linestyle, linewidth=s.linewidth,
                                        alpha=s.alpha, label=s.label, zorder=s.zorder)
        for ax in self.axes:
            if ax.lines:
                ax.legend(loc="upper center", fontsize=8, framealpha=0.85)
            ax.set_ylim(self._depth_max, self._depth_min)

        slider_ax = self.fig.add_axes([0.935, 0.12, 0.025, 0.72])
        self.depth_slider = RangeSlider(
            slider_ax, "Depth", self._depth_max, self._depth_min,
            valinit=(self._depth_max, self._depth_min), orientation="vertical",
        )
        self.depth_slider.on_changed(self._on_slider)
        self.fig.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.fig.canvas.mpl_connect("button_press_event", self._on_press)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.fig.canvas.mpl_connect("button_release_event", self._on_release)

    def _apply_depth_window(self, top: float, base: float, update_slider: bool = True) -> None:
        top, base = sorted((max(self._depth_min, top), min(self._depth_max, base)))
        if base <= top:
            return
        for ax in self.axes:
            ax.set_ylim(base, top)
        if update_slider and self.depth_slider is not None and not self._updating_slider:
            self._updating_slider = True
            self.depth_slider.set_val((top, base))
            self._updating_slider = False
        self.fig.canvas.draw_idle()

    def _on_slider(self, val) -> None:
        if not self._updating_slider:
            self._apply_depth_window(float(val[0]), float(val[1]), update_slider=False)

    def _on_scroll(self, event) -> None:
        if event.inaxes not in self.axes or event.ydata is None:
            return
        top, base = sorted(self.depth_slider.val)
        factor = 0.80 if event.button == "up" else 1.25
        new_span = min((base - top) * factor, self._depth_max - self._depth_min)
        frac = (event.ydata - top) / max(base - top, np.finfo(float).eps)
        new_top = event.ydata - frac * new_span
        new_base = new_top + new_span
        if new_top < self._depth_min:
            new_base += self._depth_min - new_top; new_top = self._depth_min
        if new_base > self._depth_max:
            new_top -= new_base - self._depth_max; new_base = self._depth_max
        self._apply_depth_window(new_top, new_base)

    def _on_press(self, event) -> None:
        if event.button == 3 and event.inaxes in self.axes and event.ydata is not None:
            self._drag_y = event.ydata

    def _on_motion(self, event) -> None:
        if self._drag_y is None or event.inaxes not in self.axes or event.ydata is None:
            return
        top, base = sorted(self.depth_slider.val)
        shift = self._drag_y - event.ydata
        span = base - top
        new_top = np.clip(top + shift, self._depth_min, self._depth_max - span)
        self._apply_depth_window(new_top, new_top + span)
        self._drag_y = event.ydata

    def _on_release(self, event) -> None:
        self._drag_y = None

    def show(self, block: bool = True):
        if self.fig is None:
            self._draw()
        plt.show(block=block)
        return self.fig, self.axes


if __name__ == "__main__":
    depth = np.linspace(1000.0, 2000.0, 2500)
    rng = np.random.default_rng(7)
    gr = np.clip(75 + 35*np.sin(depth/37) + rng.normal(0, 8, depth.size), 5, 180)
    rt = np.clip(np.exp(2.2 + 0.8*np.sin(depth/55) + rng.normal(0, .15, depth.size)), 0.2, 200)
    phi = np.clip(0.18 + 0.07*np.sin(depth/43) + rng.normal(0, .015, depth.size), 0.01, 0.35)
    sw = np.clip(0.55 + 0.35*np.sin(depth/90), 0, 1)

    viewer = WellLogViewer(3, [
        TrackSpec("Gamma Ray", "linear", (0, 180)),
        TrackSpec("Resistivity", "log", (0.2, 200)),
        TrackSpec("Porosity / Sw", "linear", (0, 1)),
    ])
    viewer.add_log("GR", gr, depth, 0, color="forestgreen", linewidth=0.9)
    viewer.add_log("RT", rt, depth, 1, color="firebrick", linewidth=0.9)
    viewer.add_log("PHI", phi, depth, 2, color="royalblue")
    viewer.add_log("SW", sw, depth, 2, color="black", linestyle="--")
    viewer.add_solid_fill(0, "GR", side="left", baseline=0, color="gold", alpha=0.35)
    viewer.add_colormapped_fill(1, "RT", side="right", baseline=200,
                                value_curve="GR", cmap="turbo", vmin=0, vmax=180)
    viewer.add_colormapped_fill(2, "PHI", curve2="SW", value_curve="SW",
                                cmap="Blues", vmin=0, vmax=1, alpha=0.8)
    viewer.show()
