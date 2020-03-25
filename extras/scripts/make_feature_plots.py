from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import sys, os, json

def ret_color_selector():

    color_selector = { "cas1": "lightblue",
                        "cas2": "seagreen",
                        "cas3": "gold",
                        "cas4": "springgreen",
                        "cas5": "darkred",
                        "cas6": "thistle",
                        "cas7": "coral",
                        "cas8": "red",
                        "cas9": "palegreen",
                        "cas10": "gray",
                        "cas11": "tan",
                        "cas12": "orange",
                        "cas13": "saddlebrown",
                        "tnsA": "navy",
                        "tnsB": "blue",
                        "tnsC": "purple",
                        "tnsD": "white",
                        "tniQ": "teal"
                        }
    
    return color_selector

def get_neighborhood_count(results_all):

    count = 0
    for result in results_all.values():
        count += len(result)

    return count

def merge_results(file_paths):

    results_all = {}
    with open(file_paths) as f:
        for path in f:
            path = path.strip()
            if path != "":
                sample_name = path.split("/")[-2]
                contig_name = path.split("/")[-1]
                result_id = "{}_{}".format(sample_name, contig_name)

                with open(path) as s:
                    results_all[result_id] = json.load(s)
    
    return results_all

def get_graphic_record(neighborhood):

    seq_len = neighborhood["n_stop"] - neighborhood["n_start"]
    color_selector = ret_color_selector()
    feature_list = []
    for hit in neighborhood["hits"].values():
        
        if "ORFID" in hit.keys():
            if hit["ORFID"].split("|")[-1] == "1":
                strand = +1
                start = int(hit["Start"])
                end = int(hit["Stop"])
            else:
                strand = -1
                start = int(hit["Stop"])
                end = int(hit["Start"])
            color = color_selector[hit["Name"]]
            feature = GraphicFeature(start=start, end=end, strand=strand, label=hit["Name"], color=color)
            feature_list.append(feature)
        
        else:
            label = "CRISPR_{}".format(hit["Copies"])
            end_pos = int(hit["Position"]) + int(hit["Length"])
            feature = GraphicFeature(start=int(hit["Position"]), end=end_pos, label=label)
            feature_list.append(feature)

    return GraphicRecord(sequence_length=seq_len, features=feature_list, first_index=int(neighborhood["n_start"]))

def make_feature_plots(file_paths):
    """Plot results from crisposon jobs using DNA features viewer.

    Args:
        file_paths (str): Text file containing paths to crisposon results.
            Expects one path per line.

    Figures contain a maximum of 20 subplots (due to matplotlib figsize constraints)
    where each subplot represents one gene neighborhood.
    """
    results_all = merge_results(file_paths)
    #total = get_neighborhood_count(results_all)
    #print(total)
   
    fig, axs = plt.subplots(20, 1, figsize=(15, 70))

    ax_index = 0
    fig_index = 0
    for contig in results_all.keys():

        count = len(results_all[contig].values())
        
        if ax_index + count > 20:
            plt.savefig("batch_{}".format(fig_index), bbox_inches="tight")
            fig, axs = plt.subplots(20, 1, figsize=(15, 70))
            fig_index += 1
            ax_index = 0

        for i, neighborhood in enumerate(results_all[contig].values()):

            if i == 0:
                axs[ax_index].set_title(contig)
            
            graphic_record = get_graphic_record(neighborhood)
            graphic_record.plot(ax=axs[ax_index], figure_width=5)
            ax_index += 1

    plt.savefig("batch_{}".format(fig_index), bbox_inches="tight")

if __name__ == "__main__":

    make_feature_plots(sys.argv[1])
