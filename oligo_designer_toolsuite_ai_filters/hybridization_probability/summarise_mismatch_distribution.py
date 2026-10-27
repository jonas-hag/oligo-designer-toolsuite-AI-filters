"""Summarise the mismatch distribution over all regions.

The input are the "*_mismatches_per_alignment.csv" files written by
determine_mismatch_distribution.py, i.e. one file per region. The mismatch counts of all
regions are pooled and the same metrics as in determine_mismatch_distribution.py are
calculated per oligo length bin, separately for

    - intergenic regions (region names starting with "InterReg")
    - genes (all other region names)
    - all regions

The region name is extracted from the file name.
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd

# imported from the script in the same directory, so that the summary uses exactly the
# same metrics as determine_mismatch_distribution.py (importing the package instead
# would pull in torch, which is not needed here)
from determine_mismatch_distribution import (
    BIN_LABELS,
    BIN_LABEL_OUT_OF_RANGE,
    distribution_by_bin,
    summarize_by_bin,
)

FILE_SUFFIX = "_mismatches_per_alignment"
PATTERN_INTERGENIC = re.compile(r"^InterReg")
# the sequences are not needed for the summary and are skipped to save memory
COLUMNS_USED = [
    "query",
    "oligo_length",
    "alignment_length",
    "n_matches",
    "n_mutations",
    "n_insertions",
    "n_deletions",
    "identical",
    "length_bin",
]


def determine_region(mismatch_file: Path) -> str:
    """Extract the region name from the name of the mismatch file."""
    region = mismatch_file.stem
    if region.endswith(FILE_SUFFIX):
        region = region[: -len(FILE_SUFFIX)]
    # the file names of the alignments start with the region name,
    # e.g. ABCA13_filtered_sampled_blast_results_filtered_sample
    region_match = re.search(r"^(.*?)(?=_filtered_sampled)", region)
    if region_match:
        region = region_match.group(1)

    return region


def determine_region_type(region: str) -> str:
    """Assign a region to intergenic regions or genes."""
    return "intergenic" if PATTERN_INTERGENIC.match(region) else "gene"


def read_mismatch_files(mismatch_files: list) -> pd.DataFrame:
    """Read the mismatch files of all regions and add the region and the region type."""
    mismatches = []
    n_regions_without_alignments = 0
    for one_mismatch_file in mismatch_files:
        one_mismatches = pd.read_csv(one_mismatch_file, usecols=COLUMNS_USED)
        if one_mismatches.empty:
            # regions without any alignment contribute nothing to the counts, but would
            # turn the count columns into object columns when they are concatenated
            n_regions_without_alignments += 1
            continue
        region = determine_region(one_mismatch_file)
        one_mismatches["region"] = region
        one_mismatches["region_type"] = determine_region_type(region)
        mismatches.append(one_mismatches)

    if n_regions_without_alignments:
        print(f"Skipped {n_regions_without_alignments} regions without any alignment.")
    if not mismatches:
        raise ValueError("None of the files contains an alignment.")

    mismatches = pd.concat(mismatches, axis=0, ignore_index=True)
    # restore the order of the length bins that is lost when writing to csv
    mismatches["length_bin"] = pd.Categorical(
        mismatches["length_bin"],
        categories=BIN_LABELS + [BIN_LABEL_OUT_OF_RANGE],
        ordered=True,
    )

    return mismatches


def summarize_region_type(mismatches: pd.DataFrame, region_type: str) -> pd.DataFrame:
    """Calculate the summary per length bin for one region type."""
    summary = summarize_by_bin(mismatches)
    n_regions = (
        mismatches.groupby("length_bin", observed=True)["region"]
        .nunique()
        .reset_index(name="n_regions")
    )
    summary = summary.merge(n_regions, on="length_bin", how="left")
    # keep the region information at the beginning of the table
    summary.insert(1, "n_regions", summary.pop("n_regions"))
    summary.insert(0, "region_type", region_type)

    return summary


def distribution_region_type(
    mismatches: pd.DataFrame, region_type: str
) -> pd.DataFrame:
    """Calculate the distribution of the mismatch counts per length bin for one region type."""
    distributions = distribution_by_bin(mismatches)
    distributions.insert(0, "region_type", region_type)

    return distributions


def main():
    #########################
    # read in arguments #
    #########################

    parser = argparse.ArgumentParser(
        prog="Summarise Mismatch Distribution",
        usage="summarise_mismatch_distribution [options]",
        description=__doc__,
    )
    parser.add_argument(
        "-i",
        "--indir",
        help="directory with the *_mismatches_per_alignment.csv files of the single regions",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="directory the results are written to (default: the input directory)",
        default=None,
    )
    parser.add_argument(
        "-p",
        "--pattern",
        help="pattern used to collect the mismatch files",
        default=f"*{FILE_SUFFIX}.csv",
    )
    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir
    os.makedirs(outdir, exist_ok=True)

    mismatch_files = sorted(f for f in indir.glob(args.pattern) if f.is_file())
    if not mismatch_files:
        raise ValueError(f"No files matching {args.pattern} found in {indir}.")

    #########################
    # summarise mismatches #
    #########################

    mismatches = read_mismatch_files(mismatch_files)

    summaries = []
    distributions = []
    for region_type in ["gene", "intergenic"]:
        mismatches_region_type = mismatches[mismatches["region_type"] == region_type]
        if mismatches_region_type.empty:
            print(f"No alignments found for region type {region_type}.")
            continue
        summaries.append(summarize_region_type(mismatches_region_type, region_type))
        distributions.append(
            distribution_region_type(mismatches_region_type, region_type)
        )
    summaries.append(summarize_region_type(mismatches, "all"))
    distributions.append(distribution_region_type(mismatches, "all"))

    summaries = pd.concat(summaries, axis=0, ignore_index=True)
    distributions = pd.concat(distributions, axis=0, ignore_index=True)

    #########################
    # write results #
    #########################

    file_summary = outdir / "mismatch_summary_by_bin_and_region_type.csv"
    file_distributions = outdir / "mismatch_distribution_by_bin_and_region_type.csv"

    summaries.to_csv(file_summary, index=False)
    distributions.to_csv(file_distributions, index=False)

    n_regions = mismatches["region"].nunique()
    n_intergenic = mismatches[mismatches["region_type"] == "intergenic"][
        "region"
    ].nunique()
    print(
        f"Summarised {mismatches.shape[0]} alignments from {len(mismatch_files)} files "
        f"({n_regions} regions: {n_regions - n_intergenic} genes, {n_intergenic} intergenic)"
    )
    print(summaries.to_string(index=False))
    print(f"Written to:\n{file_summary}\n{file_distributions}")


if __name__ == "__main__":
    main()
