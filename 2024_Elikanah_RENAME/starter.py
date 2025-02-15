# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 09:38:33 2024

@author: piercetf
"""

from skimage import io as skio

from skimage import filters, morphology, exposure, segmentation, measure

FILENAME = "C:\\Users\\piercetf\\OneDrive - National Institutes of Health\\Pictures\\T5 after injection example.jpg"

image_arr = skio.imread(FILENAME)

#skio.imshow(image_arr)

s = image_arr[:,250:900,:]

skio.imshow(s)

