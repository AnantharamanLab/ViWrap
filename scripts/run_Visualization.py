#!/usr/bin/env python3

"""
Cody Martin
University of Wisconsin-Madison
Department of Bacteriology
Anantharaman lab
"""
import argparse
from functools import partial
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

FILES = {
    "ko_rel_abun": Path("KO_ID_relative_abundance.txt"),
    "ko_met_abun": Path("KO_metabolism_relative_abundance.txt"),
    "vir_fam_abun": Path("virus_family_relative_abundance.txt"),
    "vir_stats": Path("virus_statistics.txt"),
}

PLOT_NAMES = {
    name: [file.with_suffix(".png"), file.with_suffix(".pdf")]
    for name, file in FILES.items()
}

FIGSIZES = {
    "ko_rel_abun": (12, 6),
    "ko_met_abun": (6, 6),
    "vir_fam_abun": (6, 6),
    "vir_stats": (6, 6),
}


def read_files(inputdir: Path) -> Dict[str, pd.DataFrame]:
    """Read in the 4 visualization summary files to generate plots later.

    Args:
        inputdir (Path): directory to plotting inputs

    Returns:
        Dict[str, pd.DataFrame]: maps dataframe name to dataframe
    """
    files = {name: inputdir.joinpath(file) for name, file in FILES.items()}

    ko_rel_abun = (
        pd.read_table(files["ko_rel_abun"])
        .assign(KO_num=lambda df: df["KO"].str.lstrip("K"))
        .astype({"KO_num": int})
        .sort_values(by="KO_num")
        .reset_index(drop=True)
        .drop("KO_num", axis=1)
    )
    ko_met_abun = pd.read_table(files["ko_met_abun"])

    vir_fam_abun = (
        pd.read_table(files["vir_fam_abun"])
        .assign(family=lambda df: df["family"].str.split(";").str[0])
        .rename(lambda c: c.replace(" ", "_"), axis=1)
    )

    vir_stats = (
        pd.read_table(files["vir_stats"])
        .melt(var_name="metric", value_name="Count")
        .drop_duplicates()
    )

    data = {
        name: df
        for name, df in zip(
            files.keys(), (ko_rel_abun, ko_met_abun, vir_fam_abun, vir_stats,)
        )
    }

    return data


def _plot_KO_abundance(df: pd.DataFrame, ax: plt.Axes):
    sns.barplot(data=df, x="KO", y="relative abundance", ax=ax, edgecolor="black")
    ax.xaxis.set_tick_params(rotation=45)


def _plot_KO_metabolism_abundance(
    df: pd.DataFrame, ax: plt.Axes, wedgeprops: Dict[str, str]
):
    data = df.query("`relative abundance` > 0")
    ax.pie(
        data["relative abundance"],
        labels=data["KO metabolism"],
        autopct="%.2f%%",
        wedgeprops=wedgeprops,
    )


def _plot_virus_families(df: pd.DataFrame, ax: plt.Axes, wedgeprops: Dict[str, str]):
    labels = [
        f"{row.family} ({row.relative_abundance * 100:.2f}%)" for row in df.itertuples()
    ]
    ax.pie(
        df["relative_abundance"],
        labels=labels,
        rotatelabels=True,
        wedgeprops=wedgeprops,
    )


def _plot_virus_stats(df: pd.DataFrame, ax: plt.Axes):
    sns.barplot(data=df, x="metric", y="Count", ax=ax, edgecolor="black")
    ax.xaxis.set_tick_params(rotation=45)
    ax.set(xlabel="", ylim=(0, ax.get_ylim()[1] * 1.1))
    ax.bar_label(ax.containers[0])


def plot(data: Dict[str, pd.DataFrame], outdir: Path):
    wedgeprops = {"edgecolor": "black"}
    PLTFUNCS = {
        "ko_rel_abun": _plot_KO_abundance,
        "ko_met_abun": partial(_plot_KO_metabolism_abundance, wedgeprops=wedgeprops),
        "vir_fam_abun": partial(_plot_virus_families, wedgeprops=wedgeprops),
        "vir_stats": _plot_virus_stats,
    }

    for name, df in data.items():
        figsize = FIGSIZES[name]
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)
        func = PLTFUNCS[name]
        func(df=df, ax=ax)
        fig.tight_layout()
        for outfile in PLOT_NAMES[name]:
            output = outdir.joinpath(outfile)
            fig.savefig(output.as_posix())


def main(inputdir: Path, outdir: Path):
    outdir.mkdir(exist_ok=True)
    data = read_files(inputdir)
    plot(data, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ViWrap summary plotting")
    parser.add_argument(
        "-i",
        "--input-dir",
        default="./09_Virus_statistics_visualization/Result_visualization_inputs",
        help="path to input files (default: %(default)s)",
    )
    parser.add_argument(
        "-r",
        "--viwrap-resultsdir",
        required=True,
        help="path to entire results directory from ViWrap output",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="./09_Virus_statistics_visualization/Result_visualization_outputs",
        help="path to output plot files (default: %(default)s)",
    )

    args = parser.parse_args()
    inputdir = Path(args.input_dir)
    outdir = Path(args.viwrap_resultsdir).joinpath(args.outdir)
    main(inputdir, outdir)
