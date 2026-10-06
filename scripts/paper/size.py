from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {0.0: "#3B6EA8", 1.0: "#7B3294"}
BANDS = {0.0: (0.15, 0.30), 1.0: (0.15, 0.30)}


def sem(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(values.std(ddof=1) / np.sqrt(values.size))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_completion(run_dir: Path) -> None:
    completion = json.loads((run_dir / "completion.json").read_text())
    for name, record in completion["files"].items():
        path = run_dir / name
        if not path.is_file():
            raise RuntimeError(f"missing run file: {path}")
        if path.stat().st_size != int(record["bytes"]):
            raise RuntimeError(f"size mismatch: {path}")
        if sha256(path) != record["sha256"]:
            raise RuntimeError(f"checksum mismatch: {path}")


def load_campaign(campaign: Path) -> tuple[dict, pd.DataFrame]:
    manifest = json.loads((campaign / "manifest.json").read_text())
    summary = json.loads((campaign / "analysis" / "summary.json").read_text())
    if summary["n_planned_runs"] != 30 or summary["n_complete_runs"] != 30:
        raise RuntimeError("the v2 campaign is not complete")
    if summary["missing_runs"]:
        raise RuntimeError(f"missing v2 runs: {summary['missing_runs']}")
    conditions = pd.DataFrame(summary["conditions"])
    return manifest, conditions


def analyze_runs(campaign: Path, manifest: dict, conditions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = {
        (int(row.n_chains), float(row.kappa)): float(row.primary_k_star)
        for row in conditions.itertuples()
    }
    rows = []
    health = []
    for spec in manifest["runs"]:
        run_dir = campaign / "runs" / spec["run_id"]
        verify_completion(run_dir)
        snapshots = pd.read_csv(run_dir / "snapshots.csv")
        numeric = snapshots.select_dtypes(include=[np.number]).to_numpy(dtype=float)
        if numeric.size == 0 or not np.isfinite(numeric).all():
            raise RuntimeError(f"nonfinite snapshots: {spec['run_id']}")
        meta = json.loads((run_dir / "meta.json").read_text())
        sequence = np.asarray(meta["sequence"], dtype=int)
        overlap = meta["initial_overlap_relaxation"]
        if sequence.size != int(spec["n_beads"]) or not np.isclose(sequence.mean(), 0.5):
            raise RuntimeError(f"composition mismatch: {spec['run_id']}")
        if not overlap["reached_target"]:
            raise RuntimeError(f"overlap relaxation failed: {spec['run_id']}")

        with np.load(run_dir / "direct_structure_factor.npz", allow_pickle=False) as archive:
            window = archive["window"].astype("U5")
            final = window == "final"
            early = window == "early"
            q = np.asarray(archive["q"], dtype=float)
            shell_q = np.asarray(archive["shell_q"], dtype=float)
            mode = np.asarray(archive["S_psi_psi_over_2"], dtype=float)
            shell = np.asarray(archive["S_psi_psi_over_2_shell"], dtype=float)
            expected = 0.5 * (
                np.asarray(archive["S_AA"], dtype=float)
                + np.asarray(archive["S_BB"], dtype=float)
                - 2.0 * np.asarray(archive["S_AB"], dtype=float)
            )
            np.testing.assert_allclose(mode, expected, rtol=2e-13, atol=2e-13)
        key = (int(spec["n_chains"]), float(spec["kappa"]))
        q_star = selected[key]
        shell_index = int(np.argmin(np.abs(shell_q - q_star)))
        if not np.isclose(shell_q[shell_index], q_star, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"selected shell missing: {spec['run_id']}")
        low, high = BANDS[float(spec["kappa"])]
        band_mask = (q >= low) & (q < high)
        if not band_mask.any():
            raise RuntimeError(f"empty physical-q band: {spec['run_id']}")
        peak_final = float(shell[final, shell_index].mean())
        peak_early = float(shell[early, shell_index].mean())
        band_final = float(mode[final][:, band_mask].mean())
        band_early = float(mode[early][:, band_mask].mean())
        rows.append(
            {
                "run_id": spec["run_id"],
                "n_chains": int(spec["n_chains"]),
                "n_beads": int(spec["n_beads"]),
                "box_size_sigma": float(spec["box_size"]),
                "bead_density_sigma_minus_3": float(spec["bead_density"]),
                "kappa": float(spec["kappa"]),
                "pi": float(spec["pi"]),
                "seed": int(spec["seed"]),
                "q_min_sigma_minus_1": float(2.0 * np.pi / float(spec["box_size"])),
                "q_star_sigma_minus_1": q_star,
                "peak_final": peak_final,
                "peak_early": peak_early,
                "peak_early_to_final_change": peak_final - peak_early,
                "band_q_low_sigma_minus_1": low,
                "band_q_high_sigma_minus_1": high,
                "band_final": band_final,
                "band_early": band_early,
                "band_early_to_final_change": band_final - band_early,
            }
        )
        health.append(
            {
                "run_id": spec["run_id"],
                "n_chains": int(spec["n_chains"]),
                "kappa": float(spec["kappa"]),
                "seed": int(spec["seed"]),
                "n_snapshots": int(len(snapshots)),
                "last_step": int(snapshots["step"].iloc[-1]),
                "temperature_min": float(snapshots["temperature_inst"].min()),
                "temperature_max": float(snapshots["temperature_inst"].max()),
                "energy_min": float(snapshots["total_energy"].min()),
                "energy_max": float(snapshots["total_energy"].max()),
                "radius_gyration_min": float(snapshots["mean_Rg"].min()),
                "radius_gyration_max": float(snapshots["mean_Rg"].max()),
                "overlap_min_before_sigma": float(overlap["min_separation_before"]),
                "overlap_min_after_sigma": float(overlap["min_separation_after"]),
                "overlap_target_reached": bool(overlap["reached_target"]),
                "completion_checksums_valid": True,
            }
        )
    return pd.DataFrame(rows).sort_values(["n_chains", "kappa", "seed"]), pd.DataFrame(health)


def aggregate(by_seed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (n_chains, kappa), group in by_seed.groupby(["n_chains", "kappa"], sort=True):
        rows.append(
            {
                "n_chains": int(n_chains),
                "kappa": float(kappa),
                "n_seeds": int(len(group)),
                "box_size_sigma": float(group["box_size_sigma"].iloc[0]),
                "q_min_sigma_minus_1": float(group["q_min_sigma_minus_1"].iloc[0]),
                "q_star_sigma_minus_1": float(group["q_star_sigma_minus_1"].iloc[0]),
                "peak_mean": float(group["peak_final"].mean()),
                "peak_sem": sem(group["peak_final"].to_numpy()),
                "band_q_low_sigma_minus_1": float(group["band_q_low_sigma_minus_1"].iloc[0]),
                "band_q_high_sigma_minus_1": float(group["band_q_high_sigma_minus_1"].iloc[0]),
                "band_mean": float(group["band_final"].mean()),
                "band_sem": sem(group["band_final"].to_numpy()),
                "band_early_mean": float(group["band_early"].mean()),
                "band_early_to_final_change_mean": float(group["band_early_to_final_change"].mean()),
                "band_early_to_final_change_sem": sem(group["band_early_to_final_change"].to_numpy()),
                "band_relative_change_percent": float(
                    100.0 * group["band_early_to_final_change"].mean() / group["band_final"].mean()
                ),
            }
        )
    result = pd.DataFrame(rows)
    ratios = []
    for n_chains, group in result.groupby("n_chains"):
        k0 = group[np.isclose(group["kappa"], 0.0)].iloc[0]
        k1 = group[np.isclose(group["kappa"], 1.0)].iloc[0]
        ratios.append(
            {
                "n_chains": int(n_chains),
                "exact_shell_mean_ratio_kappa1_over_kappa0": float(k1.peak_mean / k0.peak_mean),
                "fixed_band_mean_ratio_kappa1_over_kappa0": float(k1.band_mean / k0.band_mean),
            }
        )
    return result.merge(pd.DataFrame(ratios), on="n_chains", how="left")


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.annotate(
        label,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(-9, 7),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        clip_on=False,
    )


def plot(by_seed: pd.DataFrame, conditions: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 8.2,
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.0,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), constrained_layout=True)
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")
    for kappa in (0.0, 1.0):
        color = COLORS[kappa]
        condition = conditions[np.isclose(conditions["kappa"], kappa)].sort_values("n_chains")
        seed_rows = by_seed[np.isclose(by_seed["kappa"], kappa)]
        for seed_row in seed_rows.itertuples():
            offset = np.exp((seed_row.seed - 3) * 0.010 + (-0.018 if kappa == 0.0 else 0.018))
            axes[0].plot(
                seed_row.n_chains * offset,
                seed_row.peak_final,
                "o",
                ms=2.8,
                color=color,
                alpha=0.38,
                markeredgewidth=0,
                zorder=1,
            )
            axes[1].plot(
                seed_row.n_chains * offset,
                seed_row.band_final,
                "o",
                ms=2.8,
                color=color,
                alpha=0.38,
                markeredgewidth=0,
                zorder=1,
            )
        axes[0].errorbar(
            condition["n_chains"],
            condition["peak_mean"],
            yerr=condition["peak_sem"],
            marker="o",
            ms=5.2,
            lw=1.5,
            elinewidth=1.5,
            capsize=2.5,
            capthick=1.0,
            markeredgewidth=1.0,
            color=color,
            label=fr"$\kappa={kappa:g}$",
            zorder=3,
        )
        axes[1].errorbar(
            condition["n_chains"],
            condition["band_mean"],
            yerr=condition["band_sem"],
            marker="o",
            ms=5.2,
            lw=1.5,
            elinewidth=1.5,
            capsize=2.5,
            capthick=1.0,
            markeredgewidth=1.0,
            color=color,
            zorder=3,
        )
        axes[2].plot(
            condition["n_chains"],
            condition["q_star_sigma_minus_1"],
            marker="o",
            ms=5.2,
            lw=1.5,
            markeredgewidth=1.0,
            color=color,
            zorder=3,
        )
    size_rows = conditions[np.isclose(conditions["kappa"], 0.0)].sort_values("n_chains")
    axes[2].plot(
        size_rows["n_chains"],
        size_rows["q_min_sigma_minus_1"],
        "--",
        marker="x",
        ms=4.2,
        color="0.45",
        lw=1.2,
        markeredgewidth=1.0,
        label=r"$q_{\min}=2\pi/L$",
        zorder=4,
    )
    axes[0].set_ylabel(r"$S_{\psi\psi}^{(N)}(q_\star)/2$")
    axes[1].set_ylabel(r"fixed-band $S_{\psi\psi}^{(N)}/2$")
    axes[2].set_ylabel(r"$q_\star\sigma$")
    axes[2].set_ylim(0.0, 1.25)
    for axis in axes:
        axis.set_xscale("log")
        axis.set_xticks([144, 288, 576])
        axis.set_xticklabels(["144", "288", "576"])
        axis.minorticks_off()
        axis.set_xlabel(r"chains $M$ at fixed $\rho$")
        axis.grid(alpha=0.18, linewidth=1.0)
    for axis, label in zip(axes, "abc"):
        panel_label(axis, label)
    handles_a, labels_a = axes[0].get_legend_handles_labels()
    handles_c, labels_c = axes[2].get_legend_handles_labels()
    fig.legend(
        handles_a + handles_c,
        labels_a + labels_c,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.0),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output.with_suffix(".pdf"),
        bbox_inches="tight",
        pad_inches=0.025,
        facecolor="white",
        transparent=False,
    )
    fig.savefig(
        output.with_suffix(".png"),
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.025,
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)


def _verify_conditions(by_seed: pd.DataFrame, conditions: pd.DataFrame) -> None:
    expected = aggregate(by_seed).set_index(["n_chains", "kappa"]).sort_index()
    archived = conditions.set_index(["n_chains", "kappa"]).sort_index()
    pd.testing.assert_index_equal(expected.index, archived.index)
    columns = expected.columns.intersection(archived.columns)
    np.testing.assert_allclose(expected[columns], archived[columns], rtol=0, atol=1e-11)


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-data", type=Path, default=repo / "analysis" / "paper" / "source_data" / "fixed_density")
    parser.add_argument("--out", type=Path, default=repo / "output" / "paper_figures" / "submission" / "pdf" / "FigS8_finite_size")
    parser.add_argument("--campaign", type=Path)
    parser.add_argument("--tables", type=Path)
    args = parser.parse_args()
    if args.campaign is None:
        by_seed = pd.read_csv(args.source_data / "finite_size_by_seed.csv")
        conditions = pd.read_csv(args.source_data / "finite_size_conditions.csv")
        _verify_conditions(by_seed, conditions)
        plot(by_seed, conditions, args.out)
        return
    manifest, archived_conditions = load_campaign(args.campaign)
    by_seed, health = analyze_runs(args.campaign, manifest, archived_conditions)
    conditions = aggregate(by_seed)
    archived = archived_conditions[
        [
            "n_chains",
            "kappa",
            "primary_k_star",
            "primary_amplitude_seed_mean_at_selected_shell",
            "primary_amplitude_seed_sem_at_selected_shell",
        ]
    ].copy()
    merged = conditions.merge(archived, on=["n_chains", "kappa"], how="left")
    np.testing.assert_allclose(merged["q_star_sigma_minus_1"], merged["primary_k_star"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        merged["peak_mean"], merged["primary_amplitude_seed_mean_at_selected_shell"], rtol=0, atol=1e-11
    )
    np.testing.assert_allclose(
        merged["peak_sem"], merged["primary_amplitude_seed_sem_at_selected_shell"], rtol=0, atol=1e-11
    )
    tables = args.tables or args.out.parent / "finite_size_tables"
    if tables.resolve() == args.source_data.resolve():
        raise ValueError("reanalysis tables must use a different directory from source data")
    tables.mkdir(parents=True, exist_ok=True)
    by_seed.to_csv(tables / "finite_size_by_seed.csv", index=False)
    conditions.to_csv(tables / "finite_size_conditions.csv", index=False)
    health.to_csv(tables / "finite_size_health.csv", index=False)
    summary = {
        "campaign": manifest["campaign_id"],
        "n_runs": int(len(by_seed)),
        "all_runs_included": bool(len(by_seed) == 30),
        "final_frames": 5,
        "final_steps": [242000, 244000, 246000, 248000, 250000],
        "seed_replicates_per_condition": 5,
        "peak_definition": "maximum of the across-seed mean exact-shell S_psi_psi^(N)(q)/2 spectrum",
        "fixed_bands_sigma_minus_1": {
            "kappa_0": list(BANDS[0.0]),
            "kappa_1": list(BANDS[1.0]),
        },
        "conditions": conditions.to_dict(orient="records"),
    }
    (tables / "finite_size_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot(by_seed, conditions, args.out)


if __name__ == "__main__":
    main()
