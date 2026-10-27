from collections import defaultdict
from matplotlib.lines import Line2D
import matplotlib as plt
from pygenomeviz import GenomeViz
from gffpandas import gffpandas as gffpd
from matplotlib import cm
import matplotlib.colors as mcolors
import regex as re
import pandas as pd

##### Parameters
annotation_path = r'data/annotations/SARS-CoV2.gff3'
refseq_path = r'data/references/SARS-CoV2.fasta'
greac_region_path = r'data/greac_regions/SARS-CoV2.csv'
# annotation_path = r'data/annotations/MPOX.gff3'
# refseq_path = r'data/references/MPOX.fasta'
# greac_region_path = r'data/greac_regions/MPOX.csv'

HIGHLIGTH_COLOR = "#2ca02c" # verde

##### FUNCTIONS
def read_gff3(annotation_path:str):
    df_annotation = gffpd.read_gff3(annotation_path).attributes_to_columns()
    df_annotation = df_annotation[['type', 'start', 'end', 'Name', 'product']]
    return df_annotation

def get_seq_infos(refseq_path:str):
    seq_len = 0
    seq_name = refseq_path.split('.fasta')[0].split('/')[-1]
    with open(refseq_path, 'r') as seq:
        for line in seq:
            if line.startswith('>'):
                continue
            seq_len += len(line.strip())
    return seq_name, seq_len
        
def get_greac_regions(greac_region_path: str) -> list[tuple[int,int]]:
    greac_regions = []
    
    with open(greac_region_path, 'r') as greac_file:
        for line in greac_file:
            next(greac_file)
            for line in greac_file:
                _, start, end = line.split(',')
                greac_regions.append((int(start), int(end)))
    return greac_regions

def get_regions_and_colors(df):
    regions = defaultdict(list)

    # Limpar e filtrar nomes válidos
    def is_valid_name(name):
        if pd.isna(name):
            return False
        name = str(name).strip()
        if not name:
            return False
        # if name.upper() == "R":
        #     return False
        if re.match(r"^[NY]P_\d+\.\d+$", name):  # NP_057849.4, YP_...
            return False
        return True

    # Filtra somente os nomes válidos
    valid_df = df[df['Name'].apply(is_valid_name)]

    # Obtém os nomes únicos
    unique_names = valid_df['Name'].unique()

    # Gera paleta de cores
    cmap = cm.get_cmap('tab20', len(unique_names))

    def to_hex_safe(rgb):
        hex_color = mcolors.to_hex(rgb)
        return "#000000" if hex_color.lower() == "#2ca02c" else hex_color

    colors_by_name = {
        name: to_hex_safe(cmap(i))
        for i, name in enumerate(unique_names)
    }

    # Preenche regiões apenas para nomes válidos
    for _, row in valid_df.iterrows():
        regions[row['Name']].append((int(row['start']), int(row['end'])))

    return regions, colors_by_name

def make_plot(seq_name:str, seq_len:int,
              greac_regions:list[tuple[int,int]],
              regions:defaultdict[list[tuple[int,int]]],
              colors_by_name):
    handles = []
    
    gv = GenomeViz(fig_track_height=0.5, feature_track_ratio=2)
    # gv.set_scale_bar(ymargin=1.5)
    gv.set_scale_xticks()
    
    track = gv.add_feature_track(seq_name, seq_len)
    
    track.add_sublabel()
    
    for greac_region in greac_regions:
        start, end = greac_region[0], greac_region[1]
        end = seq_len if end > seq_len else end
        track.add_feature(start,
                          end,
                          -1,
                          plotstyle="box",
                          fc=HIGHLIGTH_COLOR)
                          # label="MicroSats",
                          # )#,
                          #text_kws=dict(rotation=-45, vpos="bottom"))
    handles.append(Line2D([], [],
                          marker="s",
                          color=HIGHLIGTH_COLOR,
                          label="GREAC Region",
                          ms=12,
                          ls="none"))
    
    for k, v in regions.items():
        if k in colors_by_name.keys():
            start, end = v[0][0], v[0][1]
            end = seq_len if end > seq_len else end
            track.add_feature(start,
                              end,
                              1,
                              plotstyle="box",
                              fc=colors_by_name.get(k))
            handles.append(Line2D([], [],
                                  marker="s",
                                  color=colors_by_name.get(k),
                                  label=k,
                                  ms=12,
                                  ls="none"))
        
    fig = gv.plotfig()
    
    # Plot legend for groups
    _ = fig.legend(
        handles=handles,
        fontsize=12,
        title="Regions",
        title_fontsize=12,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        handlelength=1.0,
    )

    plot_name = 'graphs/' + seq_name
    fig.savefig(plot_name, dpi=fig.dpi)
    fig.clf()

#####
seq_name, seq_len = get_seq_infos(refseq_path)
greac_regions = get_greac_regions(greac_region_path)
df = read_gff3(annotation_path)
regions, colors_by_name = get_regions_and_colors(df)
make_plot(seq_name, seq_len, greac_regions, regions, colors_by_name)