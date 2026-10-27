"""Determine the distribution of mismatch types per oligo length bin.

The input is a csv file with the alignments of the sampled oligos (e.g. the output of
generate_real_dataset_select_alignments_updated.py) that contains the columns
"query", " gapped_query", " gapped_reference" and (optionally) "oligo_length".

For every alignment of a query sequence the gapped query and the gapped reference are
compared position by position and the following mismatch types are counted:
    - mutation: both strands have a base at this position, but the bases differ
    - insertion: the reference has a "-" at this position (the query has an extra base)
    - deletion: the query has a "-" at this position (the query is missing a base)

The counts are aggregated per oligo length bin and written to three csv files.

In addition, it is counted which base of the query is aligned to which base of the
reference for all mutations (e.g. A>T, A>G, A>C, ...). These substitutions are counted
over all alignments of the file, i.e. they are not split up by oligo length bin, and are
written to a fourth csv file.
"""

import argparse
import os
from collections import Counter
from pathlib import Path

import pandas as pd

# oligos per region per oligo length
LENGTH_BINS = [
    (15, 20),
    (21, 30),
    (31, 40),
    (41, 50),
    (51, 60),
    (61, 70),
    (71, 80),
    (81, 90),
    (91, 100),
]
BIN_LABELS = [f"{lower}-{upper}" for lower, upper in LENGTH_BINS]
BIN_LABEL_OUT_OF_RANGE = "out_of_range"

MISMATCH_TYPES = {
    "mutation": "n_mutations",
    "insertion": "n_insertions",
    "deletion": "n_deletions",
}
# the counts determined for every alignment by classify_alignment
MISMATCH_COLUMNS = [
    "alignment_length",
    "n_matches",
    "n_mutations",
    "n_insertions",
    "n_deletions",
]

BASES = ["A", "C", "G", "T"]
SUBSTITUTION_COLUMNS = [
    "substitution",
    "query_base",
    "reference_base",
    "n_substitutions",
]


def assign_length_bin(oligo_length: int) -> str:
    """Assign an oligo length to one of the length bins."""
    for lower, upper in LENGTH_BINS:
        if lower <= oligo_length <= upper:
            return f"{lower}-{upper}"
    return BIN_LABEL_OUT_OF_RANGE


def classify_alignment(
    gapped_query: str, gapped_reference: str, substitutions: Counter = None
) -> dict:
    """Count matches, mutations, insertions and deletions of one alignment.

    Sequences are expected to be upper case and of equal length. Ambiguous bases ("N")
    are treated like any other base, i.e. they count as a mutation if the two strands
    differ at this position.

    If a Counter is passed as substitutions, the (query base, reference base) pair of
    every mutation is added to it. The same Counter is used for all alignments of one
    file, so that the substitutions are counted over all alignments.
    """
    n_matches = 0
    n_mutations = 0
    n_insertions = 0
    n_deletions = 0

    for base_query, base_reference in zip(gapped_query, gapped_reference):
        if base_query == "-" and base_reference == "-":
            # gap in both strands, should not occur in the blast output
            continue
        elif base_reference == "-":
            n_insertions += 1
        elif base_query == "-":
            n_deletions += 1
        elif base_query == base_reference:
            n_matches += 1
        else:
            n_mutations += 1
            if substitutions is not None:
                substitutions[(base_query, base_reference)] += 1

    return {
        "alignment_length": len(gapped_query),
        "n_matches": n_matches,
        "n_mutations": n_mutations,
        "n_insertions": n_insertions,
        "n_deletions": n_deletions,
    }


def read_alignments(alignment_file: Path) -> pd.DataFrame:
    """Read the alignment file and clean up column names and sequences.

    The column names of the alignment files have a leading whitespace
    (" gapped_query", " gapped_reference") and so do their values, therefore both are
    stripped. Lower case bases (soft-masked repeat regions) are converted to upper case
    so that they are not counted as mutations.
    """
    alignments = pd.read_csv(alignment_file, index_col=0)
    alignments.columns = alignments.columns.str.strip()

    required_columns = ["query", "gapped_query", "gapped_reference"]
    missing_columns = [c for c in required_columns if c not in alignments.columns]
    if missing_columns:
        raise ValueError(
            f"The following required columns are missing in {alignment_file}: {missing_columns}"
        )

    for column in required_columns:
        alignments[column] = alignments[column].str.strip().str.upper()

    if "oligo_length" not in alignments.columns:
        # the query column contains the un-gapped oligo sequence
        alignments["oligo_length"] = alignments["query"].str.len()
    alignments["oligo_length"] = alignments["oligo_length"].astype(int)

    return alignments


def count_substitutions(substitutions: Counter) -> pd.DataFrame:
    """Turn the counted substitutions into a table.

    All 12 combinations of the four bases are always part of the table, even if they
    did not occur. Substitutions with other bases (e.g. "N" in the reference genome)
    are appended, so that they are not lost.
    """
    combinations = [
        (query, reference)
        for query in BASES
        for reference in BASES
        if query != reference
    ]
    combinations += sorted(set(substitutions) - set(combinations))

    substitution_counts = pd.DataFrame(
        [
            {
                "substitution": f"{query}>{reference}",
                "query_base": query,
                "reference_base": reference,
                "n_substitutions": substitutions[(query, reference)],
            }
            for query, reference in combinations
        ],
        columns=SUBSTITUTION_COLUMNS,
    )

    return substitution_counts


def determine_mismatches(alignments: pd.DataFrame) -> tuple:
    """Count the mismatch types of every alignment and assign the oligo length bin.

    Returns the mismatch counts per alignment and the substitution counts over all
    alignments of the file.
    """
    unequal_length = alignments["gapped_query"].str.len() != alignments[
        "gapped_reference"
    ].str.len()
    if unequal_length.any():
        print(
            f"Skipping {unequal_length.sum()} alignments where gapped_query and "
            "gapped_reference have a different length."
        )
        alignments = alignments[~unequal_length]

    # the substitutions of all alignments of this file are collected in one Counter
    substitutions = Counter()

    # the columns are set explicitly, otherwise alignment files without any alignment
    # would result in a DataFrame without any column
    mismatches = pd.DataFrame(
        [
            classify_alignment(gapped_query, gapped_reference, substitutions)
            for gapped_query, gapped_reference in zip(
                alignments["gapped_query"], alignments["gapped_reference"]
            )
        ],
        index=alignments.index,
        columns=MISMATCH_COLUMNS,
    ).astype(int)

    mismatches = pd.concat(
        [
            alignments[["query", "gapped_query", "gapped_reference", "oligo_length"]],
            mismatches,
        ],
        axis=1,
    )
    mismatches["identical"] = (
        mismatches["gapped_query"] == mismatches["gapped_reference"]
    )
    mismatches["length_bin"] = pd.Categorical(
        mismatches["oligo_length"].apply(assign_length_bin),
        categories=BIN_LABELS + [BIN_LABEL_OUT_OF_RANGE],
        ordered=True,
    )

    return mismatches, count_substitutions(substitutions)


def summarize_by_bin(mismatches: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the mismatch counts per oligo length bin."""
    summary = (
        mismatches.groupby("length_bin", observed=True)
        .agg(
            n_alignments=("query", "size"),
            n_queries=("query", "nunique"),
            n_identical_alignments=("identical", "sum"),
            n_aligned_positions=("alignment_length", "sum"),
            n_mutations=("n_mutations", "sum"),
            n_insertions=("n_insertions", "sum"),
            n_deletions=("n_deletions", "sum"),
            mean_mutations=("n_mutations", "mean"),
            mean_insertions=("n_insertions", "mean"),
            mean_deletions=("n_deletions", "mean"),
            median_mutations=("n_mutations", "median"),
            median_insertions=("n_insertions", "median"),
            median_deletions=("n_deletions", "median"),
            max_mutations=("n_mutations", "max"),
            max_insertions=("n_insertions", "max"),
            max_deletions=("n_deletions", "max"),
        )
        .reset_index()
    )

    # per-base rates, so that bins with different oligo lengths can be compared
    for mismatch_type, column in MISMATCH_TYPES.items():
        summary[f"{mismatch_type}_rate_per_base"] = (
            summary[column] / summary["n_aligned_positions"]
        )

    return summary


def distribution_by_bin(mismatches: pd.DataFrame) -> pd.DataFrame:
    """Count how many alignments have how many mutations/insertions/deletions per bin."""
    n_alignments_per_bin = mismatches.groupby("length_bin", observed=True).size()

    distributions = []
    for mismatch_type, column in MISMATCH_TYPES.items():
        one_distribution = (
            mismatches.groupby(["length_bin", column], observed=True)
            .size()
            .reset_index(name="n_alignments")
            .rename(columns={column: "n_bases"})
        )
        one_distribution.insert(1, "mismatch_type", mismatch_type)
        distributions.append(one_distribution)

    distributions = pd.concat(distributions, axis=0, ignore_index=True)
    distributions["fraction_of_alignments_in_bin"] = distributions[
        "n_alignments"
    ] / distributions["length_bin"].map(n_alignments_per_bin)
    distributions = distributions.sort_values(
        ["length_bin", "mismatch_type", "n_bases"]
    ).reset_index(drop=True)

    return distributions


def main():
    #########################
    # read in arguments #
    #########################

    parser = argparse.ArgumentParser(
        prog="Mismatch Distribution",
        usage="determine_mismatch_distribution [options]",
        description=__doc__,
    )
    parser.add_argument(
        "-i",
        "--input",
        help="path to the csv file with the alignments",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="directory the results are written to (default: directory of the input file)",
        default=None,
    )
    args = parser.parse_args()

    alignment_file = Path(args.input)
    outdir = Path(args.outdir) if args.outdir else alignment_file.parent
    os.makedirs(outdir, exist_ok=True)

    #########################
    # determine mismatches #
    #########################

    alignments = read_alignments(alignment_file)
    mismatches, substitution_counts = determine_mismatches(alignments)
    if mismatches.empty:
        print(f"No alignments found in {alignment_file}, writing empty result files.")
    summary = summarize_by_bin(mismatches)
    distributions = distribution_by_bin(mismatches)

    #########################
    # write results #
    #########################

    file_mismatches = outdir / f"{alignment_file.stem}_mismatches_per_alignment.csv"
    file_summary = outdir / f"{alignment_file.stem}_mismatch_summary_by_bin.csv"
    file_distributions = (
        outdir / f"{alignment_file.stem}_mismatch_distribution_by_bin.csv"
    )
    file_substitutions = outdir / f"{alignment_file.stem}_substitution_counts.csv"

    mismatches.to_csv(file_mismatches, index=False)
    summary.to_csv(file_summary, index=False)
    distributions.to_csv(file_distributions, index=False)
    substitution_counts.to_csv(file_substitutions, index=False)

    print(f"Analysed {mismatches.shape[0]} alignments from {alignment_file}")
    print(summary.to_string(index=False))
    print(
        f"Written to:\n{file_mismatches}\n{file_summary}\n{file_distributions}"
        f"\n{file_substitutions}"
    )


if __name__ == "__main__":
    main()
