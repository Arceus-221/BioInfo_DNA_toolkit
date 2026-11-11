"""Matplotlib plots saved as PNG files."""
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def safe_name(text: str) -> str:
    return re.sub(r"[^\w.-]", "_", text)


METRICS = [
    ("gc_content", "GC content", "GC content (%)", None),
    ("gc_skew", "GC skew", "GC skew (G-C)/(G+C), unitless", 0),
    ("at_skew", "AT skew", "AT skew (A-T)/(A+T), unitless", 0),
    ("entropy", "Shannon entropy", "Entropy (bits)", None),
]


def plot_window_metrics(rid: str, windows: list, window: int, step: int, outdir: str) -> list:
    """Four PNGs per sequence: GC content, GC skew, AT skew, entropy vs position."""
    os.makedirs(outdir, exist_ok=True)
    paths = []
    xs = [(w["start"] + w["end"]) / 2 for w in windows]
    for key, title, ylabel, baseline in METRICS:
        ys = [float("nan") if w[key] is None else w[key] for w in windows]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(xs, ys, marker="o" if len(xs) < 40 else None, linewidth=1.6,
                label=f"{title} (window {window}, step {step})")
        if baseline is not None:
            ax.axhline(baseline, color="gray", linestyle="--", linewidth=1, label="Zero line")
        ax.set_title(f"{title} along {rid}")
        ax.set_xlabel("Genomic position (bp, window midpoint, 1-based)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(loc="best")
        fig.tight_layout()
        path = os.path.join(outdir, f"{safe_name(rid)}_{key}.png")
        fig.savefig(path, dpi=130)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_benchmark(rows: list, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    for scenario, style in (("random", "-"), ("repeat", "--")):
        sub = [r for r in rows if r["scenario"] == scenario]
        if not sub:
            continue
        sizes = [r["size"] for r in sub]
        ax.plot(sizes, [r["naive_ms"] for r in sub], style, marker="o", color="tab:red",
                label=f"Naive ({scenario})")
        ax.plot(sizes, [r["kmp_ms"] for r in sub], style, marker="s", color="tab:blue",
                label=f"KMP ({scenario})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Naive vs KMP pattern matching: median runtime")
    ax.set_xlabel("Sequence length (bases)")
    ax.set_ylabel("Median time (ms)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    path = os.path.join(outdir, "benchmark_naive_vs_kmp.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
