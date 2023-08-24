from Methods import *
from matplotlib import pyplot
import os


def Plot_Skew_Diagram(genome):
    skew_values = Skew_Diagram(genome)

    x_vals = [i for i in range(len(skew_values))]
    y_vals = skew_values

    pyplot.plot(x_vals, y_vals)
    pyplot.ylabel("Skew")
    pyplot.xlabel("Position")

    pyplot.show();

def Find_Ori(genome, window_size, d):
    minimum_skew = Minimun_Skew(genome)[0]
    ref_text = genome[minimum_skew:minimum_skew+window_size]
    dnaA_boxes = Frequent_Words_With_Mismathces_And_Reverse_Complements(ref_text, 9, 1)

    oris = []
    for dnaA_box in dnaA_boxes:
        if ref_text.find(dnaA_box) > 0:
            oris.append(dnaA_box)

    return oris


        



