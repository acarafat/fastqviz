###############################
# Importing necessary Modules #
###############################

import streamlit as st
from streamlit.logger import get_logger
import numpy as np
from io import StringIO

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d, LinearColorMapper
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

LOGGER = get_logger(__name__)

def run():
    st.set_page_config(
        page_title="fastqViz",
        page_icon="🧬",
    )

    st.markdown("""
                # FastqViz
                ### A *fastq* alignment quality visualization tool
                
                Please upload a fastq file to start:
                """)

    uploaded_file = st.file_uploader("Choose a fast file:")

    if uploaded_file is None:
        st.write("""
                 Here showing an example fastq file output.
                 """)
        p = view_alignment('example/test.fastq')
    else:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        p = view_alignment(stringio)
        
    st.bokeh_chart(p,  use_container_width=True)


color_mapper = LinearColorMapper(
    palette='Magma256',
    low=40,
    high=0)

def get_colors(seqs):
    """make colors for bases in sequence"""
    quality = [i for rec in seqs for i in rec.letter_annotations['phred_quality'] ]
    return quality


def view_alignment(aln, fontsize="8pt", plot_width=800):
    """Bokeh sequence alignment view"""

    # make sequence and id lists from the aln object
    recs = [rec for rec in SeqIO.parse(aln, 'fastq')]
    seqs = [rec.seq for rec in recs]
    ids = [rec.id for rec in recs]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(recs)
    N = len(seqs[0])
    S = len(seqs)

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
        viewlen=80
    else:
        viewlen=N

    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width= plot_width, plot_height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color={'field':'colors','transform':color_mapper},
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
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




if __name__ == "__main__":
    run()
