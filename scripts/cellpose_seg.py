import argparse
import pandas as pd
from cellpose import models, utils, io
import skimage.io as skio
import pandas as pd
import numpy as np
import stardist


# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("dapi", help="path to dapi image file")
parser.add_argument("reads", help="path to csv file containing transcript reads and coordinates")
parser.add_argument("output_file", help="output csv file path")
args = parser.parse_args()
prefix = args.output_file[:-4]


dapi = skio.imread(args.dapi)
reads = pd.read_csv(args.reads)


cp_model = models.Cellpose(gpu=True, model_type='nuclei')
channels = [[0,0]]
masks,flows,styles,diams = cp_model.eval(dapi, channels=channels, flow_threshold=3, cellprob_threshold=-6)
io.masks_flows_to_seg(dapi, masks, flows, diams, prefix, channels)
cp_masks = np.load(prefix+"_seg.npy", allow_pickle=True).item()['masks']


reads['int_coords'] = reads.apply(lambda x: (int(np.rint(x['Y'])), int(np.rint(x['X']))), axis=1)
reads['cp_label'] = reads.apply(lambda x: cp_masks[x['int_coords']], axis = 1)


reads.to_csv(args.output_file)