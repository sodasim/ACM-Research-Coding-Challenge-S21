# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 16:18:02 2021

@author: sodasim
"""
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SeqIO.read("Genome.gb", "genbank")

gd_first_set = GenomeDiagram.FeatureSet()
gd_second_set = GenomeDiagram.FeatureSet()
gd_total_set = GenomeDiagram.FeatureSet()


for feature in record.features:
    firstTrack = True
    if feature.type != "gene":
        continue
    for other_feature in gd_first_set.get_features():
        if other_feature.type != "gene":
            continue
        # Determines if feature should be placed on the second track
        if feature.location.start < other_feature.location.end and feature.location.start > other_feature.location.start:
            if len(gd_second_set) % 2 == 0:
                color = colors.blue
            else:
                color = colors.lightblue
            gd_second_set.add_feature(feature, color=color, label=True, sigil="OCTO", label_size=25)    
            if len(gd_total_set) % 2 == 0:
                color = colors.red
            else:
                color = colors.pink
            gd_total_set.add_feature(feature, color=color, label=True, sigil="OCTO", label_size=25)    
            firstTrack = False
    if firstTrack:
        if len(gd_first_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        gd_first_set.add_feature(feature, color=color, label=True, sigil="OCTO", label_size=25)
        if len(gd_total_set) % 2 == 0:
            color = colors.red
        else:
            color = colors.pink
        gd_total_set.add_feature(feature, color=color, label=True, sigil="OCTO", label_size=25)
    
first_track_for_features = GenomeDiagram.Track()
second_track_for_features = GenomeDiagram.Track()
total_track = GenomeDiagram.Track()
gd_diagram = GenomeDiagram.Diagram("Tomato curly stunt virus")
gd_diagram_overlap = GenomeDiagram.Diagram("Tomato curly stunt virus Overlap")

total_track.add_set(gd_total_set)
gd_diagram_overlap.add_track(total_track, 1)

first_track_for_features.add_set(gd_first_set)
gd_diagram.add_track(first_track_for_features, 2)
second_track_for_features.add_set(gd_second_set)
gd_diagram.add_track(second_track_for_features, 1)

gd_diagram.draw(
    format="circular",
    circular=True,
    pagesize=(40 * cm, 40 * cm),
    start=0,
    end=len(record),
    circle_core=0.7,
)

gd_diagram_overlap.draw(
    format="circular",
    circular=True,
    pagesize=(40 * cm, 40 * cm),
    start=0,
    end=len(record),
    circle_core=0.7,
)

gd_diagram.write("tomato_circular.png", "PNG")
gd_diagram_overlap.write("tomato_circular_overlap.png", "PNG")