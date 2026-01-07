import os
from pathlib import Path
import argparse
import yaml
import pandas as pd

def select_alignments(alignment_file, config):
    alignment_file_long = os.path.join(config["alignment_files_long"], alignment_file.name)
    alignments_short = pd.read_csv(alignment_file)
    alignments_long = pd.read_csv(alignment_file_long)

    # soft-mask the found alignments
    # I messed up the column names, gapped_query and gapped_reference have a leading whitespace
    alignments_short = alignments_short[~alignments_short[' gapped_reference'].str.contains(r'[a-z]', na=False)]
    alignments_long = alignments_long[~alignments_long[' gapped_reference'].str.contains(r'[a-z]', na=False)]
    
    # remove duplicated alignments
    alignments_short.drop_duplicates(subset=["query", " gapped_reference"], inplace=True)
    alignments_long.drop_duplicates(subset=["query", " gapped_reference"], inplace=True)

    # calculate original oligo length and select alignments
    alignments = pd.concat([alignments_short, alignments_long], axis=0)
    alignments["oligo_length"] = alignments["query"].str.len()
    interval_data = config["interval_config"]

    data_intervals = []
    for interval_name, interval_data in interval_data.items():
        temp_data = alignments[(alignments["oligo_length"] >= interval_data["lower"]) &
                                       (alignments["oligo_length"] <= interval_data["upper"])]
        # filter out alignments to itself
        temp_data = temp_data[temp_data[" gapped_query"] != temp_data[" gapped_reference"]]
        # data_for_oligo_selection = temp_data.drop_duplicates(subset=["query"])
        data_for_oligo_selection = temp_data.value_counts("query").reset_index(name="count")
        data_for_oligo_selection = data_for_oligo_selection[data_for_oligo_selection["count"] >= 5]
        oligos = data_for_oligo_selection[["query"]]

        region_name = alignment_file.name.replace("_blast_results.csv", "")
        number_oligos_at_least_5_assignments = oligos.shape[0]
        total_number_alignments = temp_data.shape[0]
        one_dataset = {"region": region_name, "interval": interval_name,
                       "n_oligos_at_least_5_assignments": number_oligos_at_least_5_assignments,
                       "n_total_alignments": total_number_alignments}
        data_intervals.append(one_dataset)

    return data_intervals

        

def main():
    # check how many alignments are available in the following intervals:
        # 15-20: 20
        # 21-30: 40
        # 31-40: 40
        # 41-50: 40
        # 51-60: 12
        # 61-70: 12
        # 71-80: 12
        # 81-90: 12
        # 91-100: 12

    #########################
    # read in arguments #
    #########################

    parser = argparse.ArgumentParser(
        prog="Real Dataset",
        usage="generate_real_dataset [options]",
        description=main.__doc__,
    )
    parser.add_argument("-c", "--config", help="path to the configuration file", default="config/generate_real_dataset_blastn.yaml")
    args = parser.parse_args()
    with open(args.config, "r") as handle:
        config = yaml.safe_load(handle)

    # as default use results from short oligos
    alignment_files_short = [f for f in Path(config["alignment_files_short"]).glob("*.csv") if f.is_file()]

    information_alignments = [select_alignments(
            alignment_file=one_alignment_file,
            config=config
        )
        for one_alignment_file in alignment_files_short
    ]

    information_alignments = [item for sublist in information_alignments for item in sublist]

    df_alignments = pd.DataFrame(information_alignments)
    out_path = Path(config["outfile"]) / "information_alignments.csv"
    df_alignments.to_csv(out_path)

if __name__ == "__main__":
    main()