# Copyright (c) Streamlit Inc. (2018-2022) Snowflake Inc. (2022)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###############################
# Importing necessary Modules #
###############################

import streamlit as st
from streamlit.logger import get_logger
import os, io, random
import string
import numpy as np

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

#import panel as pn
import pandas as pd

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, LinearColorMapper
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

LOGGER = get_logger(__name__)

def run():
    st.write('Input: fastq alignment')

    # Need two functionality
    # -  Select fastq sequence / browse
    # - Select reference sequence

    mutSig = extractMutSig('example/test.fastq')
    p = view_alignment(mutSig)
    st.bokeh_chart(p)

        

    ##############################################
    # Definign necessary Functions and variables #
    ##############################################




    st.set_page_config(
        page_title="Hello",
        page_icon="👋",
    )

    st.write("# Welcome to Streamlit! 👋")

    st.sidebar.success("Select a demo above.")

color_mapper = LinearColorMapper(
    palette='Magma256',
    low=0,
    high=40)

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    quality = [i for rec in seqs for i in rec.letter_annotations['phred_quality'] ]
    return quality


def view_alignment(aln, fontsize="9pt"):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(aln)
    N = len(seqs[0])
    S = len(seqs)
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, #plot_width= plot_width, plot_height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color={'field':'colors','transform':color_mapper},
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, #plot_width=plot_width, plot_height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color={'field':'colors','transform':color_mapper},
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p],[p1]], toolbar_location='below')
    return p

# Extract mutation signature and quality-score for all reads in an alignment

def extractMutSig(fastq):
    extractedFastq = []
    for seq_record in SeqIO.parse(fastq, 'fastq'):
        seq = ''
        qual = []
        for i in range(len(seq)):
            base = seq_record.seq[i-1]
            quality = seq_record.letter_annotations['phred_quality'][i - 1]
            seq += base
            qual.append(quality)
        signatureSeq = SeqRecord(
            Seq(seq),
            id = seq_record.id)
        signatureSeq.letter_annotations["phred_quality"] = qual
        extractedFastq.append(signatureSeq)
    return extractedFastq


if __name__ == "__main__":
    run()
