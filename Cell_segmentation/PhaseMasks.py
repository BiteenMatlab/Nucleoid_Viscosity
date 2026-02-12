# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 20:17:00 2021

@author: azaldegc

Python script to get phase masks from phase contrast images for single molecule
analysis using SMALL LABS
Activate SingleMolecule environment to run
conda activate SingleMolecule
python PhaseMasks.py pathtofiles\*.tif
Make sure that the only tif files are the phase constrast captures
"""

import sys
import numpy as np
import glob

import tifffile as tif
from PIL import Image
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
from cellpose import models, io, plot
from skimage import (
    color, feature, filters, measure, morphology, segmentation, util
)


# function to look at multiple files
def filepull(directory):
    '''Finds images files in directory and returns them'''
    # create list of image files in directory
    filenames = [img for img in glob.glob(directory)]   
    
    return filenames


def segment_cells(model, filename, img, user_diam, PlotFig=False):
    
    blur_img = filters.gaussian(img, sigma=0.1)# blur image before cellpose
  
    #channels = [0,0]
    channels = [0,0]
    masks, flows, styles, = model.eval(blur_img, diameter=user_diam, 
                                              channels=channels)
    
    #io.masks_flows_to_seg(blur_img, masks, flows, filename, channels)
    
    masks = morphology.remove_small_objects(masks, 200) 
    masks = segmentation.clear_border(masks)
    # results for verification, if running through data set as false
    if PlotFig == True:
        
        fig = plt.figure(figsize=(12,5))
        plot.show_segmentation(fig, img, masks, flows[0], channels=channels)
        plt.tight_layout()
        plt.show()

    return masks, blur_img

# foldername that contains all of the fits.mat files
directory = sys.argv[1]
# pull in all of the fits.mat files
tif_files = filepull(directory)

print("{} .tif files found".format(len(tif_files)))

#---------------------------------------
plot_fig = False
# name = 'Cwx2565_18h_230622'

model_nuc = models.Cellpose(gpu=False, model_type='nuclei')
model_ec = models.CellposeModel(gpu=False, pretrained_model=r'C:\Users\xiaofend\Documents\Xiaofeng\CODE\Cellpose\cellpose_residual_on_style_on_concatenation_off_train_2021_07_22_00_30_46.822532')
model_yh = models.CellposeModel(gpu=False, pretrained_model=r'C:\Users\xiaofend\Documents\Xiaofeng\CODE\Cellpose\cellpose_residual_on_style_on_concatenation_off_train_2021_07_22_00_30_46.822532')
model_try = models.CellposeModel(gpu=False, pretrained_model=r'C:\Users\xiaofend\Documents\Xiaofeng\CODE\Cellpose\cellpose_residual_on_style_on_concatenation_off_ecolitraining_2021_08_08_14_40_46.082526')

model = model_ec
#-----------------------------------------

f = 1
for ii,file in enumerate(tif_files):
    
    
    img = Image.open(file)
    
    img_inv = np.invert(img)
    img_inv = np.invert(img_inv)
  
    segmented_cells, img_process = segment_cells(model, file, img_inv,
                                     user_diam = None)
    
    tif.imsave(file[:-4] + '_PhaseMask.tif'.format(num=f),
              segmented_cells)
    
    if plot_fig == True:
        fig, ax = plt.subplots(ncols=3, figsize=(6, 3), 
                               sharey=True, sharex=True)
        ax[0].imshow(img, cmap='gray')
        ax[0].contour(segmented_cells, z=0, levels=1, linewidth=0.05, colors='r')
        ax[0].set_title('Original Image')
        ax[0].axis('off')
        ax[1].imshow(img_process, cmap='gray')
        ax[1].contour(segmented_cells, z=0, levels=1, linewidth=0.05, colors='r')
        ax[1].set_title('Pre-processed')
        ax[1].axis('off') 
        ax[2].imshow(segmented_cells, cmap='turbo')
        ax[2].set_title('Segmented cells')
        ax[2].axis('off')
    
        fig.tight_layout()
        plt.show()
    
    print(file[:-4] + '_PhaseMask.tif'.format(num=f))
    f += 1
    
