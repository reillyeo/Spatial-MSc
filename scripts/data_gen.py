# This script reads baysor results table and converts to single-cell csv that can be read into seurat

import argparse
import pandas as pd

# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("data", help="path to files containing baysor results")
#parser.add_argument("section_prefix", help="character string of file prefix for section to be analysed")
parser.add_argument("output_file", help="Output CSV file path")
args = parser.parse_args()

df = pd.read_csv(args.data)
stats = pd.read_csv(args.data+'_cell_stats.csv')

df = df[df.cell != 0]
df = df[df.assignment_confidence > 0.8]
stats = stats[['cell', 'density', 'area', 'n_transcripts', 'x', 'y']]
df = (df.pivot_table(index=['cell'], columns="gene", values="molecule_id", aggfunc='count', fill_value=0))
df = pd.merge(df, stats, on = 'cell', how = 'inner')
for column in ['n_transcripts', 'area', 'density', 'cell', 'y', 'x']:
    df.insert(0, column, df.pop(column))
df.to_csv(args.output_file, index=False)





