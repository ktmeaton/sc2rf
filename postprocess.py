#!/usr/bin/env python3

import pandas as pd
import click
import os
import logging
import requests

NO_DATA_CHAR = "NA"
LAPIS_URL_BASE = "https://lapis.cov-spectrum.org/open/v1/sample/aggregated?fields=pangoLineage&nucMutations="
LINEAGE_PROP_THRESHOLD = 0.01

@click.command()
@click.option("--csv", help="CSV output from sc2rf.", required=True)
@click.option("--ansi", help="ANSI output from sc2rf.", required=False)
@click.option("--motifs", help="TSV of breakpoint motifs", required=False)
@click.option("--prefix", help="Prefix for output files.", required=False, default="sc2rf.recombinants")
@click.option("--min-len", help="Minimum region length.", required=False, default=1000)
@click.option(
    "--max-parents", help="Maximum number of parents.", required=False, default=2
)
@click.option("--outdir", help="Output directory", required=False, default=".")
@click.option(
    "--aligned",
    help="Extract recombinants from this alignment (Note: requires seqkit)",
    required=False,
)
@click.option(
    "--nextclade",
    help="The TSV output of Nextclade, to use the substitutions columns to find the parent lineages.",
    required=False,
)
@click.option(
    "--max-breakpoints",
    help="The maximum number of breakpoints",
    required=False,
    default=2,
)
@click.option(
    "--issues", help="Issues TSV metadata from pango-designation (https://github.com/ktmeaton/ncov-recombinant/raw/master/resources/issues.tsv)", required=False
)
@click.option("--log", help="Path to a log file", required=False, default="postprocess.log")
def main(
    csv,
    ansi,
    prefix,
    min_len,
    outdir,
    aligned,
    nextclade,
    max_parents,
    issues,
    max_breakpoints,
    motifs,
    log,
):
    """Detect recombinant seqences from sc2rf. Dependencies: pandas, click"""

    # Check for directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # create logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    fh = logging.FileHandler(log)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    # -----------------------------------------------------------------------------
    # Import Dataframe

    # sc2rf output (required)
    logging.info("Reading input data: {}".format(csv))

    df = pd.read_csv(csv, sep=",", index_col=0)
    df.fillna("", inplace=True)
    df["sc2rf_status"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_details"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_lineage"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_clades_filter"] = [NO_DATA_CHAR] * len(df)   
    df["sc2rf_regions_filter"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_regions_length"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_breakpoints_filter"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_num_breakpoints_filter"] = [NO_DATA_CHAR] * len(df)
    df["sc2rf_breakpoints_motif"] = [NO_DATA_CHAR] * len(df)
    df["cov-spectrum_parents"] = [NO_DATA_CHAR] * len(df)
    df["cov-spectrum_parents_confidence"] = [NO_DATA_CHAR] * len(df)
    df["cov-spectrum_parents_subs"] = [NO_DATA_CHAR] * len(df)

    # if using issues.tsv of pango-designation issues (optional)
    # does lineage assignment by parent+breakpoint matching
    if issues:

        logging.info("Parsing issues: {}".format(issues))

        breakpoint_col = "breakpoints_curated"
        parents_col = "parents_curated"
        breakpoint_df = pd.read_csv(issues, sep="\t")
        breakpoint_df.fillna(NO_DATA_CHAR, inplace=True)
        drop_rows = breakpoint_df[breakpoint_df[breakpoint_col] == NO_DATA_CHAR].index
        breakpoint_df.drop(drop_rows, inplace=True)

        # Convert CSV to lists
        breakpoint_df[breakpoint_col] = [
            bp.split(",") for bp in breakpoint_df[breakpoint_col]
        ]
        breakpoint_df[parents_col] = [p.split(",") for p in breakpoint_df[parents_col]]

        # Consider a breakpoint match if within 50 base pairs
        breakpoint_approx_bp = 50

    # (Optional) motifs dataframe
    if motifs:
        logging.info("Parsing motifs: {}".format(motifs))
        motifs_df = pd.read_csv(motifs, sep="\t")

    # (Optional) nextclade tsv dataframe
    if nextclade:
        logging.info("Parsing nextclade: {}".format(nextclade))
        nextclade_df = pd.read_csv(nextclade, sep="\t", index_col=0)
        nextclade_df.fillna(NO_DATA_CHAR, inplace=True)

    # Initialize a dictionary of false_positive strains
    # key: strain, value: reason
    false_positives = {}

    logging.info("Post-processing table: {}".format(csv))

    for rec in df.iterrows():

        strain = rec[0]

        regions_str = rec[1]["regions"]
        regions_split = regions_str.split(",")

        # Keys are going to be the start coord of the region
        regions_filter = {}
        breakpoints_filter = []

        prev_clade = None
        prev_start_coord = 0
        prev_end_coord = 0

        for region in regions_split:
            coords = region.split("|")[0]
            clade = region.split("|")[1]
            start_coord = int(coords.split(":")[0])
            end_coord = int(coords.split(":")[1])
            region_len = (end_coord - start_coord) + 1
            coord_list = list(regions_filter)
            coord_list.reverse()

            # -----------------------------------------------------------------
            # FIRST PASS

            # Just ignore singletons, no calculation necessary
            if region_len == 1:
                continue

            # Is this the first region?
            if not prev_clade:
                regions_filter[start_coord] = {"clade": clade, "end": end_coord}
                prev_clade = clade
                prev_start_coord = start_coord

            # -----------------------------------------------------------------
            # SECOND PASS: REVERSE ITER COLLAPSE

            for coord in coord_list:
                prev_start_coord = coord
                prev_end_coord = regions_filter[prev_start_coord]["end"]
                prev_region_len = (prev_end_coord - prev_start_coord) + 1
                prev_clade = regions_filter[coord]["clade"]

                # If the previous region was too short AND from a different clade
                # Delete that previous region, it's an intermission
                if prev_region_len < min_len and clade != prev_clade:
                    del regions_filter[prev_start_coord]

                # Collapse the current region into the previous one
                elif clade == prev_clade:
                    regions_filter[prev_start_coord]["end"] = end_coord
                    break

                # Otherwise, clades differ and this is the start of a new region
                else:
                    regions_filter[start_coord] = {"clade": clade, "end": end_coord}
                    break

            # Check if the reveres iter collapse wound up deleting all the regions
            if len(regions_filter) == 0:
                regions_filter[start_coord] = {"clade": clade, "end": end_coord}

            # These get updated regardless of condition
            prev_clade = clade
            prev_end_coord = end_coord

        # Check the last region for length
        if len(regions_filter) > 1:
            start_coord = list(regions_filter)[-1]
            end_coord = regions_filter[start_coord]["end"]
            region_len = end_coord - start_coord
            if region_len < min_len:
                del regions_filter[start_coord]

        # Check if all the regions were collapsed
        if len(regions_filter) < 2:
            false_positives[rec[0]] = "single parent"

        # -----------------------------------------------------------------
        # THIRD PASS: BREAKPOINT DETECTION

        prev_start_coord = None
        for start_coord in regions_filter:

            end_coord = regions_filter[start_coord]["end"]

            # Skip the first record for breakpoints
            if prev_start_coord:
                breakpoint_start = prev_end_coord + 1
                breakpoint_end = start_coord - 1
                breakpoint = "{}:{}".format(breakpoint_start, breakpoint_end)
                breakpoints_filter.append(breakpoint)

            prev_start_coord = start_coord
            prev_end_coord = end_coord

        # check if the number of breakpoints changed
        # the filtered breakpoints should only ever be equal or less
        # 2022-06-17: Why? Under what conditions does the filtered breakpoints increase?
        # Except! If the breakpoints were initially 0
        #num_breakpoints = df["breakpoints"][strain]
        num_breakpoints_filter = len(breakpoints_filter)

        # Check for too many breakpoints
        if num_breakpoints_filter > max_breakpoints:
            false_positives[strain] = "{} breakpoints > {} max breakpoints".format(
                num_breakpoints_filter,
                max_breakpoints,
            )

        # Check if postprocessing increased the number of breakpoints
        # Why is this bad again? Should be fine as long as we're under the max?
        # if (num_breakpoints > 0) and (num_breakpoints_filter > num_breakpoints):
        #     false_positives[
        #         strain
        #     ] = "{} filtered breakpoints > {} raw breakpoints".format(
        #         num_breakpoints_filter, num_breakpoints
        #     )

        # Identify the new filtered clades
        clades_filter = [regions_filter[s]["clade"] for s in regions_filter]
        #clades_filter_csv = ",".join(clades_filter)
        num_parents = len(set(clades_filter))
        if num_parents > max_parents:
            false_positives[strain] = "{} parents > {}".format(num_parents, max_parents)

        # Extract the lengths of each region
        regions_length = [str(regions_filter[s]["end"] - s) for s in regions_filter]

        # Construct the new filtered regions
        regions_filter = [
            "{}:{}|{}".format(s, regions_filter[s]["end"], regions_filter[s]["clade"])
            for s in regions_filter
        ]

        # Identify lineage based on breakpoint and parents!
        # But only if we've suppled the issues.tsv for pango-designation
        if issues:
            sc2rf_lineage = ""
            sc2rf_lineages = {bp_s: [] for bp_s in breakpoints_filter}

            for bp_s in breakpoints_filter:
                start_s = int(bp_s.split(":")[0])
                end_s = int(bp_s.split(":")[1])

                match_found = False

                for bp_rec in breakpoint_df.iterrows():

                    # Skip over this potential lineage if parents are wrong
                    bp_parents = bp_rec[1][parents_col]
                    if bp_parents != clades_filter:
                        continue

                    for bp_i in bp_rec[1][breakpoint_col]:

                        start_i = int(bp_i.split(":")[0])
                        end_i = int(bp_i.split(":")[1])
                        start_diff = abs(start_s - start_i)
                        end_diff = abs(end_s - end_i)

                        if (
                            start_diff <= breakpoint_approx_bp
                            and end_diff <= breakpoint_approx_bp
                        ):

                            sc2rf_lineages[bp_s].append(bp_rec[1]["lineage"])
                            match_found = True

                if not match_found:
                    sc2rf_lineages[bp_s].append(NO_DATA_CHAR)

            # if len(sc2rf_lineages) == num_breakpoints_filter:
            collapse_lineages = []
            for bp in sc2rf_lineages.values():
                for lineage in bp:
                    collapse_lineages.append(lineage)

            collapse_lineages = list(set(collapse_lineages))

            # When there are multiple breakpoint, a match must be the same for all!
            collapse_lineages_filter = []
            for lin in collapse_lineages:

                if lin == NO_DATA_CHAR:
                    continue
                # By default, assume they all match
                matches_all_bp = True
                for bp_s in sc2rf_lineages:
                    # If the lineage is missing, it's not in all bp
                    if lin not in sc2rf_lineages[bp_s]:
                        matches_all_bp = False
                        break

                # Check if we should drop it
                if matches_all_bp:
                    collapse_lineages_filter.append(lin)

            if len(collapse_lineages_filter) == 0:
                collapse_lineages_filter = [NO_DATA_CHAR]

            sc2rf_lineage = ",".join(collapse_lineages_filter)
            df.at[strain, "sc2rf_lineage"] = sc2rf_lineage

        # check for breakpoint motifs, to override lineage call
        # all breakpoints must include a motif!
        # ---------------------------------------------------------------------
        if motifs:
            breakpoints_motifs = []
            for bp in breakpoints_filter:
                bp_motif = False
                bp_start = int(bp.split(":")[0])
                bp_end = int(bp.split(":")[1])

                # Add buffers
                bp_start = bp_start - breakpoint_approx_bp
                bp_end = bp_end + breakpoint_approx_bp

                for motif_rec in motifs_df.iterrows():
                    motif_start = motif_rec[1]["start"]
                    motif_end = motif_rec[1]["end"]

                    # Is motif contained within the breakpoint
                    # Allow fuzzy matching
                    if motif_start >= bp_start and motif_end <= bp_end:
                        #print("\t\t", motif_start, motif_end)
                        bp_motif = True

                breakpoints_motifs.append(bp_motif)

            # If there's a lone "False" value, it gets re-coded to an empty string on export
            # To prevent that, force it to be the NO_DATA_CHAR (ex. 'NA')
            if len(breakpoints_motifs) == 0:
                breakpoints_motifs_str = [NO_DATA_CHAR]
            else:
                breakpoints_motifs_str = [str(m) for m in breakpoints_motifs]

            df.at[strain,"sc2rf_breakpoints_motif"] = ",".join(breakpoints_motifs_str)

            # Override the linaege call if one breakpoint had not motif
            if False in breakpoints_motifs:
                df.at[strain, "sc2rf_lineage"] = "false_positive"
                df.at[strain, "sc2rf_status"]  = "false_positive"

        # ---------------------------------------------------------------------
        # Identify parent lineages by querying cov-spectrum mutations

        if nextclade:

            parent_lineages = []
            parent_lineages_confidence = []
            parent_lineages_subs = []

            substitutions = nextclade_df["substitutions"][strain].split(",")
            unlabeled_privates = nextclade_df["privateNucMutations.unlabeledSubstitutions"][strain].split(",")

            # Remove NA char
            if NO_DATA_CHAR in substitutions:
                substitutions.remove(NO_DATA_CHAR)
            if NO_DATA_CHAR in unlabeled_privates:
                unlabeled_privates.remove(NO_DATA_CHAR)

            # Exclude privates from mutations to query
            for private in unlabeled_privates:
                # Might not be in there if it's an indel
                if private in substitutions:
                    substitutions.remove(private) 

            # Split mutations by region
            for region in regions_filter:
                region_coords = region.split("|")[0]
                region_start = int(region_coords.split(":")[0])
                region_end = int(region_coords.split(":")[1])
                region_subs = []

                for sub in substitutions:
                    sub_coord = int(sub[1:-1])
                    if sub_coord >= region_start and sub_coord <= region_end:
                        region_subs.append(sub)

                region_subs_csv = ",".join(region_subs)
                # Query cov-spectrum
                url = LAPIS_URL_BASE + region_subs_csv
                r = requests.get(url)
                result = r.json()
                lineage_data = result["data"]

                # Have keys be counts
                lineage_dict = {}
                for rec in lineage_data:
                    lineage = rec["pangoLineage"]
                    count = rec["count"]
                    lineage_dict[count] = lineage

                # Sort in order
                lineage_dict = {k:lineage_dict[k] for k in sorted(lineage_dict, reverse=True)}

                if len(lineage_dict) == 0:
                    max_lineage = NO_DATA_CHAR
                    max_prop = NO_DATA_CHAR
                else:
                    # Temporarily set to fake data
                    total_lineages = sum(lineage_dict.keys())
                    max_count = max(lineage_dict.keys())
                    max_prop = max_count / total_lineages                
                    max_lineage = lineage_dict[max_count]

                    # Don't want to report recombinants as parents yet
                    while max_lineage.startswith("X") or max_lineage == "Unassigned":
                        lineage_dict = {
                            count:lineage for count,lineage in lineage_dict.items() 
                            if lineage != max_lineage
                        }
                        # If there are no other options, set to NA
                        if len(lineage_dict) == 0:
                            max_lineage = NO_DATA_CHAR
                            break
                        # Otherwise try again!
                        else:
                            # For now, deliberately don't update total_lineages
                            max_count = max(lineage_dict.keys())
                            max_prop = max_count / total_lineages                
                            max_lineage = lineage_dict[max_count]               
                
                parent_lineages_sub_str = "{}|{}".format(region_subs_csv, max_lineage)

                parent_lineages.append(max_lineage)
                parent_lineages_confidence.append(max_prop)
                parent_lineages_subs.append(parent_lineages_sub_str)

            # Update the dataframe columns
            df.at[strain, "cov-spectrum_parents"] = ",".join(parent_lineages)
            df.at[strain, "cov-spectrum_parents_confidence"] = ",".join(
                str(round(c,3)) if type(c) == float else NO_DATA_CHAR
                for c in parent_lineages_confidence
            )
            df.at[strain, "cov-spectrum_parents_subs"] = ";".join(parent_lineages_subs)

        df.at[strain, "sc2rf_clades_filter"] = ",".join(clades_filter)
        df.at[strain, "sc2rf_regions_filter"] = ",".join(regions_filter)
        df.at[strain, "sc2rf_regions_length"] = ",".join(regions_length)
        df.at[strain, "sc2rf_breakpoints_filter"] = ",".join(breakpoints_filter)
        df.at[strain, "sc2rf_num_breakpoints_filter"] = num_breakpoints_filter
        if strain in false_positives:
            df.at[strain, "sc2rf_status"]  = "false_positive"
            df.at[strain, "sc2rf_details"] = false_positives[strain]
            df.at[strain, "sc2rf_breakpoints_filter"]  = NO_DATA_CHAR
        else:
            df.at[strain, "sc2rf_status"]  = "positive"

    # write exclude strains
    outpath_exclude = os.path.join(outdir, prefix + ".exclude.tsv")
    if len(false_positives) > 0:
        with open(outpath_exclude, "w") as outfile:
            for strain, reason in false_positives.items():
                outfile.write(strain + "\t" + reason + "\n")
    else:
        cmd = "touch {outpath}".format(outpath=outpath_exclude)
        os.system(cmd)
        
    # drop strains
    #false_positives = set(false_positives.keys())
    #df.drop(false_positives, inplace=True)



    # -------------------------------------------------------------------------
    # Add in the Negatives (if alignment was specified)
    # Avoiding the Bio module, I just need names not sequence

    if aligned:
        logging.info("Reporting non-recombinants in the alignment: {}".format(aligned))        
        with open(aligned) as infile:
            aligned_content = infile.read()
            for line in aligned_content.split("\n"):
                if line.startswith(">"):
                    strain = line.replace(">","")
                    # Ignore this strain if it's already in dataframe (it's a recombinant)
                    if strain in df.index: continue
                    # Otherwise add it, with no data as default           
                    df.loc[strain] = NO_DATA_CHAR
                    df.at[strain, "sc2rf_status"]  = "negative"
                    df.at[strain, "sc2rf_details"]  = "no recombination detected"

    # -------------------------------------------------------------------------
    # write output table

    # Drop old columns
    df.drop(
        ["examples", "intermissions", "breakpoints", "regions"],
        axis="columns",
        inplace=True,
    )
    df.insert(loc=0, column="strain", value=df.index)
    df.rename(
        {
            "sc2rf_clades_filter": "sc2rf_parents",
            "sc2rf_regions_filter": "sc2rf_regions",
            "sc2rf_breakpoints_filter": "sc2rf_breakpoints",
            "sc2rf_num_breakpoints_filter": "sc2rf_num_breakpoints",
        },
        axis="columns",
        inplace=True,
    )
    outpath_rec = os.path.join(outdir, prefix + ".tsv")

    logging.info("Writing the output table: {}".format(outpath_rec))    
    df.to_csv(outpath_rec, sep="\t", index=False)

    # -------------------------------------------------------------------------
    # write output strains
    outpath_strains = os.path.join(outdir, prefix + ".txt")
    strains_df = df[
        (df["sc2rf_status"] != "negative") 
        & (df["sc2rf_status"] != "false_positive")
    ]
    strains = list(strains_df.index)
    strains_txt = "\n".join(strains)
    with open(outpath_strains, "w") as outfile:
        outfile.write(strains_txt)

    # -------------------------------------------------------------------------
    # filter the ansi output
    if ansi:
        
        outpath_ansi = os.path.join(outdir, prefix + ".ansi.txt")
        logging.info("Writing filtered ansi input: {}".format(outpath_ansi))        
        if len(false_positives) > 0:
            cmd = "cut -f 1 {exclude} | grep -v -f - {inpath} > {outpath}".format(
                exclude=outpath_exclude,
                inpath=ansi,
                outpath=outpath_ansi,
            )
        else:
            cmd = "cp -f {inpath} {outpath}".format(
                inpath=ansi,
                outpath=outpath_ansi,
            )
        os.system(cmd)

    # -------------------------------------------------------------------------
    # write alignment
    if aligned:
        outpath_fasta = os.path.join(outdir, prefix + ".fasta")
        logging.info("Writing filtered alignment: {}".format(outpath_fasta)) 

        cmd = "seqkit grep -f {outpath_strains} {aligned} > {outpath_fasta};".format(
            outpath_strains=outpath_strains,
            aligned=aligned,
            outpath_fasta=outpath_fasta,
        )
        os.system(cmd)


if __name__ == "__main__":
    main()