import os
from pathlib import Path
import argparse
import yaml
import pandas as pd
import numpy as np


def select_alignments(alignment_file, config, seed):
    alignment_file_long = os.path.join(
        config["alignment_files_long"], alignment_file.name
    )
    alignments_short = pd.read_csv(alignment_file)
    alignments_long = pd.read_csv(alignment_file_long)

    # ignore lower-case matches
    # I messed up the column names, gapped_query and gapped_reference have a leading whitespace
    alignments_short[" gapped_reference"] = alignments_short[
        " gapped_reference"
    ].str.upper()
    alignments_long[" gapped_reference"] = alignments_long[
        " gapped_reference"
    ].str.upper()

    # remove duplicated alignments
    alignments_short.drop_duplicates(
        subset=["query", " gapped_reference"], inplace=True
    )
    alignments_long.drop_duplicates(subset=["query", " gapped_reference"], inplace=True)

    # calculate original oligo length and select alignments
    alignments = pd.concat([alignments_short, alignments_long], axis=0)
    alignments["oligo_length"] = alignments["query"].str.len()
    interval_data = config["interval_config"]

    all_sampled_alignments = []
    for interval_name, interval_data in interval_data.items():
        print(interval_name)
        temp_data = alignments[
            (alignments["oligo_length"] >= interval_data["lower"])
            & (alignments["oligo_length"] <= interval_data["upper"])
        ]
        # filter out alignments to itself
        temp_data = temp_data[
            temp_data[" gapped_query"] != temp_data[" gapped_reference"]
        ]
        data_for_oligo_selection = temp_data.value_counts("query").reset_index(
            name="count"
        )
        # only use oligos with at least 5 alignments as this is what we are looking for
        # keep the rest as backup if not enough alignments are available
        data_for_oligo_selection_backup = data_for_oligo_selection[
            data_for_oligo_selection["count"] < 5
        ]
        data_for_oligo_selection = data_for_oligo_selection[
            data_for_oligo_selection["count"] >= 5
        ]
        oligos = data_for_oligo_selection[["query"]]
        oligos_backup = data_for_oligo_selection_backup[["query"]]
        print(data_for_oligo_selection)
        print(interval_data["n"])
        if oligos.shape[0] >= interval_data["n"]:
            oligos = oligos.sample(n=interval_data["n"], random_state=seed[0])
            temp_data_selected_alignments = temp_data.merge(
                oligos, how="inner", on="query"
            )
            temp_data_selected_alignments = (
                temp_data_selected_alignments.groupby(["query"])
                .sample(n=5, random_state=seed[1])
                .reset_index(drop=True)
            )
        else:
            temp_data_selected_alignments_1 = temp_data.merge(
                oligos, how="inner", on="query"
            )
            temp_data_selected_alignments_1 = (
                temp_data_selected_alignments_1.groupby(["query"])
                .sample(n=5, random_state=seed[2])
                .reset_index(drop=True)
            )
            n_remaining_oligos = interval_data["n"] - oligos.shape[0]
            if n_remaining_oligos > oligos_backup.shape[0]:
                print(
                    f"Not enough oligos with fewer than 5 alignments: {n_remaining_oligos} required, but only {oligos_backup.shape[0]} available."
                )
                n_remaining_oligos = oligos_backup.shape[0]
            oligos_backup = oligos_backup.sample(
                n=n_remaining_oligos, random_state=seed[3]
            )
            temp_data_selected_alignments_2 = temp_data.merge(
                oligos_backup, how="inner", on="query"
            )
            temp_data_selected_alignments_2 = (
                temp_data_selected_alignments_2.groupby(["query"])[
                    ["query", " gapped_query", " gapped_reference", "oligo_length"]
                ]
                .apply(
                    lambda g: g.sample(n=min(len(g), 5), random_state=seed[4]),
                    include_groups=False,
                )
                .reset_index(drop=True)
            )
            temp_data_selected_alignments = pd.concat(
                [temp_data_selected_alignments_1, temp_data_selected_alignments_2],
                axis=0,
                ignore_index=True,
            )
        all_sampled_alignments.append(temp_data_selected_alignments)
    all_sampled_alignments = pd.concat(
        all_sampled_alignments, axis=0, ignore_index=True
    )

    return all_sampled_alignments


def main():
    # 1. filter alignments: no lowercase bases and no duplicates
    # 2. select oligos according to the length schema:
    # 15-20: 20
    # 21-30: 40
    # 31-40: 40
    # 41-50: 40
    # 51-60: 12
    # 61-70: 12
    # 71-80: 12
    # 81-90: 12
    # 91-100: 12
    # 2. select 5 hits per oligo -> what to do when selected oligos don't have enough alignments?
    # 3. calculate NUPACK binding prop -> for 6 temperatures
    # 4. distribute to train/test/validation split

    #########################
    # read in arguments #
    #########################

    parser = argparse.ArgumentParser(
        prog="Real Dataset",
        usage="generate_real_dataset [options]",
        description=main.__doc__,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the configuration file",
        default="config/generate_real_dataset_blastn.yaml",
    )
    args = parser.parse_args()
    with open(args.config, "r") as handle:
        config = yaml.safe_load(handle)

    # as default use results from short oligos
    alignment_files_short = [
        f for f in Path(config["alignment_files_short"]).glob("*.csv") if f.is_file()
    ]

    np.random.seed(config["seed"])
    # in the function, up to 5 seeds are needed per region
    list_of_seeds = [
        np.random.randint(1e8, size=5) for i in range(len(alignment_files_short))
    ]

    selected_alignments = [
        select_alignments(
            alignment_file=one_alignment_file,
            config=config,
            seed=one_seed,
        )
        for one_alignment_file, one_seed in zip(alignment_files_short, list_of_seeds)
    ]

    for one_df, alignment_file in zip(selected_alignments, alignment_files_short):
        one_df.to_csv(
            Path(config["outdir"]) / (alignment_file.stem + "_filtered_sample.csv")
        )


if __name__ == "__main__":
    main()
