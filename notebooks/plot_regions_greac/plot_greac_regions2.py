from collections import defaultdict
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from pygenomeviz import GenomeViz
from gffpandas import gffpandas as gffpd
import matplotlib.colors as mcolors
import regex as re
import pandas as pd
import numpy as np
import csv
import os # Para criar a pasta de gráficos
import sys


##### Parameters
# (anotação, referência, regiões por variante exportadas pelo wdwclassifier.ipynb, eventos ou None)
ORGANISMS = {
    'SARS-CoV2': ('data/annotations/SARS-CoV2.gff3', 'data/references/SARS-CoV2.fasta',
                  'data/greac_regions/SARS-CoV2.csv', 'data/annotations/SARS-CoV2-events.tsv'),
    'MPOX': ('data/annotations/MPOX2.gff3', 'data/references/MPOX.fasta',
             'data/greac_regions/MPOX.csv', 'data/annotations/MPOX-events.tsv'),
    'DENV': ('data/annotations/DENV2.gff3', 'data/references/DENV.fasta',
             'data/greac_regions/DENV.csv', 'data/annotations/DENV-events.tsv'),
    'HBV': ('data/annotations/HBV2.gff3', 'data/references/HBV.fasta',
            'data/greac_regions/HBV.csv', 'data/annotations/HBV-events.tsv'),
    'HIV': ('data/annotations/HIV.gff3', 'data/references/HIV.fasta',
            'data/greac_regions/HIV.csv', None),
}
# MPOX tem ~180 genes: mostra só os que têm eventos (lista 'teste' comentada no fim do arquivo)
GENES_FILTER = {
    'MPOX': ['OPG174', 'OPG047', 'OPG125', 'OPG210', 'OPG057', 'OPG083', 'OPG002_dup', 'OPG015', 'OPG025',
             'OPG122', 'OPG185', 'OPG198', 'OPG074', 'OPG021', 'OPG133', 'OPG001', 'OPG001_dup', 'OPG188', 'OPG190'],
}

HIGHLIGTH_COLOR = "#2ca02c" # verde

# Bloco novo e corrigido
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.title_fontsize': 11,
    'legend.fontsize': 10
})

# Cria o diretório para os gráficos se não existir
if not os.path.exists('graphs'):
    os.makedirs('graphs')

##### FUNCTIONS
def read_events(events_file: str, seq_len: int) -> list[int]:
    events = [0] * seq_len
    with open(events_file, 'r') as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            try:
                base_str, val_str = line.strip().split('\t')
                base = int(base_str)
                val = int(val_str)
                if 1 <= base <= seq_len:
                    events[base - 1] = val
            except ValueError as e:
                print(f"Aviso: Ignorando linha mal formatada: '{line.strip()}' - Erro: {e}")
    return events


def read_gff3(annotation_path:str):
    df_annotation = gffpd.read_gff3(annotation_path).df
    # attributes_to_columns() quebra com ';' no fim (DENV2.gff3), que também só tem 'gene' em vez de 'Name'
    attrs = df_annotation['attributes'].map(lambda a: dict(kv.split('=', 1) for kv in a.split(';') if '=' in kv))
    df_annotation['Name'] = attrs.map(lambda d: d.get('Name', d.get('gene')))
    return df_annotation[['type', 'start', 'end', 'Name']]

def get_seq_infos(refseq_path:str):
    seq_len = 0
    seq_name = os.path.splitext(os.path.basename(refseq_path))[0]
    with open(refseq_path, 'r') as seq:
        for line in seq:
            if line.startswith('>'):
                continue
            seq_len += len(line.strip())
    return seq_name, seq_len

def get_greac_regions(greac_region_path: str) -> dict[str, list[tuple[int, int]]]:
    """Regiões de cada variante (colunas variant,region_init,region_end)."""
    greac_regions = defaultdict(list)
    with open(greac_region_path, 'r') as greac_file:
        for row in csv.DictReader(greac_file):
            greac_regions[row['variant']].append((int(row['region_init']), int(row['region_end'])))
    return greac_regions

def get_regions_and_colors(df):
    regions = defaultdict(list)

    def is_valid_name(name):
        if pd.isna(name): return False
        name = str(name).strip()
        if not name: return False
        if re.match(r"^[NY]P_\d+\.\d+$", name): return False
        return True

    valid_df = df[df['Name'].apply(is_valid_name)].copy()
    unique_names = valid_df['Name'].unique()

    cmap = plt.get_cmap('tab20', len(unique_names))

    def to_hex_safe(rgb):
        hex_color = mcolors.to_hex(rgb)
        return "#000000" if hex_color.lower() == HIGHLIGTH_COLOR else hex_color

    colors_by_name = {
        name: to_hex_safe(cmap(i))
        for i, name in enumerate(unique_names)
    }

    for _, row in valid_df.iterrows():
        regions[row['Name']].append((int(row['start']), int(row['end'])))

    return regions, colors_by_name

def make_plot(seq_name:str, seq_len:int,
              greac_regions:dict[str, list[tuple[int,int]]],
              regions:defaultdict[list[tuple[int,int]]],
              colors_by_name: dict,
              events: list[int] | None):
    handles = []

    gv = GenomeViz(fig_track_height=0.2, feature_track_ratio=0.8)

    # Create tracks - pygenomeviz doesn't use track_ratio parameter
    events_track = gv.add_feature_track("", seq_len) if events else None
    feature_track = gv.add_feature_track(" ", seq_len)
    feature_track.add_sublabel()

    for k, v in regions.items():
        if k in colors_by_name:
            start, end = v[0][0], v[0][1]
            ## included alpha value
            feature_track.add_feature(start, min(end, seq_len), 1, plotstyle="box", fc=colors_by_name.get(k))#, alpha=0.5)
            handles.append(Line2D([], [], marker="s", color=colors_by_name.get(k), label=k, ms=12, ls="none"))

    # Uma linha por variante, com as regiões da máscara de cada uma
    variant_tracks = []
    for variant, variant_regions in greac_regions.items():
        track = gv.add_feature_track(variant, seq_len, labelsize=11)
        for start, end in variant_regions:
            if start <= seq_len:
                track.add_feature(start, min(end, seq_len), 1, plotstyle="box", fc=HIGHLIGTH_COLOR)
        variant_tracks.append(track)
    handles.insert(0, Line2D([], [], marker="s", color=HIGHLIGTH_COLOR, label="GREAC Region", ms=12, ls="none"))

    fig = gv.plotfig()

    fig.subplots_adjust(right=0.85, top=0.95)

    if events_track is not None:
        # Trilha de eventos 3x mais alta, no topo; as demais descem juntas logo abaixo
        events_ax = events_track.ax
        events_pos = events_ax.get_position()
        new_events_height = events_pos.height * 3.0
        new_events_y = 0.95 - new_events_height
        events_ax.set_position([events_pos.x0, new_events_y, events_pos.width, new_events_height])

        first_pos = feature_track.ax.get_position()
        shift = (new_events_y - 0.02) - first_pos.y1
        for track in [feature_track] + variant_tracks:
            pos = track.ax.get_position()
            track.ax.set_position([pos.x0, pos.y0 + shift, pos.width, pos.height])

        x_coords = np.arange(1, seq_len + 1)
        event_indices = np.where(np.array(events) > 0)[0]

        if len(event_indices) > 0:
            event_values = np.array(events)[event_indices]
            x_coords_mut = x_coords[event_indices]

            # Mantendo sua escolha de paleta de cores
            cmap_mutations = plt.get_cmap('turbo')
            norm = mcolors.Normalize(vmin=min(event_values), vmax=max(event_values))

            events_ax.vlines(
                x=x_coords_mut, ymin=0, ymax=event_values,
                colors=cmap_mutations(norm(event_values)),
                alpha=0.9, linewidth=2.0
            )

            sm = plt.cm.ScalarMappable(cmap=cmap_mutations, norm=norm)
            sm.set_array([])

            # Colorbar à esquerda da trilha de eventos
            cax = fig.add_axes([events_pos.x0 - 0.06, new_events_y, 0.01, new_events_height])
            cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
            cbar.ax.tick_params(labelsize=8)

        # Estilização (como na sua versão)
        events_ax.set_ylabel("")
        events_ax.set_yticklabels([])
        events_ax.tick_params(axis='y', length=0)
        events_ax.grid(True, axis='y', linestyle=':', alpha=0.6, color='grey')

        events_ax.set_xlim(0, seq_len)
        if max(events) > 0:
            events_ax.set_ylim(0, max(events) * 1.1)
        else:
            events_ax.set_ylim(0, 1)
        events_ax.spines[['top', 'right', 'left']].set_visible(False)
        events_ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    feature_track.ax.tick_params(axis='x', labelsize=9)

    # Legenda das regiões
    _ = fig.legend(
        handles=handles, fontsize=10, title="Regions", title_fontsize=10,
        loc="center left", bbox_to_anchor=(1.01, 0.5), handlelength=1.0,
        ncol=1 if len(handles) <= 13 else 2,
    )

    plot_name = f'graphs/{seq_name}_variants.png'
    fig.savefig(plot_name, dpi=300, bbox_inches='tight', pad_inches=0.02)
    print(f"Gráfico salvo em: {plot_name}")
    plt.close(fig)

#####
# python plot_greac_regions2.py [ORGANISMO ...]  (sem argumentos: todos)
for organism in sys.argv[1:] or ORGANISMS:
    annotation_path, refseq_path, greac_region_path, events_file = ORGANISMS[organism]
    seq_name, seq_len = get_seq_infos(refseq_path)
    greac_regions = get_greac_regions(greac_region_path)
    df = read_gff3(annotation_path)
    if organism in GENES_FILTER:
        df = df[df['Name'].isin(GENES_FILTER[organism])]
    events = read_events(events_file, seq_len) if events_file else None
    regions, colors_by_name = get_regions_and_colors(df)

    make_plot(seq_name, seq_len, greac_regions, regions, colors_by_name, events)

# genes = [
#     {"Name": "C",    "start": 102,  "end": 440},
#     {"Name": "pr",   "start": 441,  "end": 713},
#     {"Name": "M",    "start": 441,  "end": 938},
#     {"Name": "E",    "start": 939,  "end": 2423},
#     {"Name": "NS1",  "start": 2424, "end": 3479},
#     {"Name": "NS2A", "start": 3480, "end": 4133},
#     {"Name": "NS2B", "start": 4134, "end": 4523},
#     {"Name": "NS3",  "start": 4524, "end": 6377},
#     {"Name": "NS4A", "start": 6378, "end": 6758},
#     {"Name": "2K",   "start": 6759, "end": 6827},
#     {"Name": "NS4B", "start": 6828, "end": 7562},
#     {"Name": "NS5",  "start": 7563, "end": 10262},
# ]
# df = pd.DataFrame(genes, columns=["Name", "start", "end"])
#HBV
# genes = [{'Name': 'X', 'end': 1840, 'start': 1376},
#  {'Name': 'C', 'end': 4807, 'start': 1816},
#  {'Name': 'pre-capsid', 'end': 2454, 'start': 1816},
#  {'Name': 'capsid', 'end': 2454, 'start': 1903},
#  {'Name': 'pol', 'end': 4807, 'start': 2309},
#  {'Name': 'S', 'end': 4019, 'start': 2850},
#  {'Name': 'envL', 'end': 4019, 'start': 2850},
#  {'Name': 'envM', 'end': 4019, 'start': 3174},
#  {'Name': 'envS', 'end': 837, 'start': 157}]

# df = pd.DataFrame(genes, columns=["Name", "start", "end"])

# teste = [
#     {
#         'locus_tag': 'NBT03_gp150',
#         'Name': 'OPG174',
#         'product': 'Hydroxysteroid dehydrogenase',
#         'start': 148355,
#         'end': 149395,
#         'max_event_frequency': 2,
#         'event_position': 148412
#     },
#     {
#         'locus_tag': 'NBT03_gp027',
#         'Name': 'OPG047',
#         'product': 'Kelch-like protein (2)',
#         'start': 29746,
#         'end': 31194,
#         'max_event_frequency': 18,
#         'event_position': 29763
#     },
#     {
#         'locus_tag': 'NBT03_gp027',
#         'Name': 'OPG047',
#         'product': 'Kelch-like protein (2)',
#         'start': 29746,
#         'end': 31194,
#         'max_event_frequency': 9,
#         'event_position': 29758
#     },
#     {
#         'locus_tag': 'NBT03_gp104',
#         'Name': 'OPG125',
#         'product': 'Rifampicin resistance protein',
#         'start': 105999,
#         'end': 107654,
#         'max_event_frequency': 11,
#         'event_position': 107670
#     },
#     {
#         'locus_tag': 'NBT03_gp173',
#         'Name': 'OPG210',
#         'product': 'B22R family protein',
#         'start': 181356,
#         'end': 186998,
#         'max_event_frequency': 11,
#         'event_position': 172087
#     },
#     {
#         'locus_tag': 'NBT03_gp037',
#         'Name': 'OPG057',
#         'product': 'Palmytilated EEV membrane protein',
#         'start': 39077,
#         'end': 40195,
#         'max_event_frequency': 12,
#         'event_position': 39072
#     },
#     {
#         'locus_tag': 'NBT03_gp037',
#         'Name': 'OPG057',
#         'product': 'Palmytilated EEV membrane protein',
#         'start': 39077,
#         'end': 40195,
#         'max_event_frequency': 8,
#         'event_position': 39070
#     },
#     {
#         'locus_tag': 'NBT03_gp062',
#         'Name': 'OPG083',
#         'product': 'Viral core cysteine proteinase',
#         'start': 63445,
#         'end': 64716,
#         'max_event_frequency': 6,
#         'event_position': 63811
#     },
#     {
#         'locus_tag': 'NBT03_gp178',
#         'Name': 'OPG002_dup',
#         'product': 'Crm-B secreted TNF-alpha-receptor-like protein',
#         'start': 194459,
#         'end': 195508,
#         'max_event_frequency': 6,
#         'event_position': 18107
#     },
#     {
#         'locus_tag': 'NBT03_gp004',
#         'Name': 'OPG015',
#         'product': 'Ankyrin repeat protein (39)',
#         'start': 4843,
#         'end': 6156,
#         'max_event_frequency': 5,
#         'event_position': 3987
#     },
#     {
#         'locus_tag': 'NBT03_gp010',
#         'Name': 'OPG025',
#         'product': 'Ankyrin repeat protein (14)',
#         'start': 13375,
#         'end': 15267,
#         'max_event_frequency': 5,
#         'event_position': 13686
#     },
#     {
#         'locus_tag': 'NBT03_gp101',
#         'Name': 'OPG122',
#         'product': 'MutT motif protein',
#         'start': 102428,
#         'end': 103174,
#         'max_event_frequency': 6,
#         'event_position': 102694
#     },
#     {
#         'locus_tag': 'NBT03_gp156',
#         'Name': 'OPG185',
#         'product': 'Hemagglutinin',
#         'start': 158917,
#         'end': 159858,
#         'max_event_frequency': 5,
#         'event_position': 159415
#     },
#     {
#         'locus_tag': 'NBT03_gp166',
#         'Name': 'OPG198',
#         'product': 'Ser/thr kinase',
#         'start': 170038,
#         'end': 170886,
#         'max_event_frequency': 4,
#         'event_position': 170573
#     },
#     {
#         'locus_tag': 'NBT03_gp053',
#         'Name': 'OPG074',
#         'product': 'Iev morphogenesis protein',
#         'start': 55128,
#         'end': 57125,
#         'max_event_frequency': 4,
#         'event_position': 51559
#     },
#     {
#         'locus_tag': 'NBT03_gp006',
#         'Name': 'OPG021',
#         'product': 'Zinc finger-like protein (2)',
#         'start': 9587,
#         'end': 10315,
#         'max_event_frequency': 4,
#         'event_position': 8599
#     },
#     {
#         'locus_tag': 'NBT03_gp112',
#         'Name': 'OPG133',
#         'product': 'Early transcription factor 82 kDa subunit',
#         'start': 113564,
#         'end': 115696,
#         'max_event_frequency': 4,
#         'event_position': 113201
#     },
#     {
#         'locus_tag': 'NBT03_gp001',
#         'Name': 'OPG001',
#         'product': 'Chemokine binding protein',
#         'start': 835,
#         'end': 1575,
#         'max_event_frequency': 4,
#         'event_position': 3818
#     },
#     {
#         'locus_tag': 'NBT03_gp179',
#         'Name': 'OPG001_dup',
#         'product': 'Chemokine binding protein',
#         'start': 195635,
#         'end': 196375,
#         'max_event_frequency': 4,
#         'event_position': 189258
#     },
#     {
#         'locus_tag': 'NBT03_gp158',
#         'Name': 'OPG188',
#         'product': 'Schlafen (1)',
#         'start': 161596,
#         'end': 163107,
#         'max_event_frequency': 4,
#         'event_position': 161575
#     },
#     {
#         'locus_tag': 'NBT03_gp160',
#         'Name': 'OPG190',
#         'product': 'EEV type-I membrane glycoprotein',
#         'start': 165122,
#         'end': 166075,
#         'max_event_frequency': 4,
#         'event_position': 165304
#     }
# ]

# df = pd.DataFrame(teste)
