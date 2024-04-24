#!/usr/bin/env python3

from PIL import Image
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class honours_plot:
    """
        Plots an RM grid as a bubble plot.

        Must enter plt.show after executing this function.
    """
    def plot_RMs(rms, scale=1):
        rms_pos = np.array(list(filter(lambda val: val[2] > 0, rms.transpose()))).transpose()
        rms_neg = np.array(list(filter(lambda val: val[2] < 0, rms.transpose()))).transpose()

        plt.scatter(rms_pos[0] * scale, rms_pos[1] * scale, s=rms_pos[2], c=(1, 1, 1, 0), edgecolors='red', linewidth=2)
        plt.scatter(rms_neg[0] * scale, rms_neg[1] * scale, s=-1 * rms_neg[2], c=(1, 1, 1, 0), edgecolors='blue', linewidth=2)

    def hi():
        print("Hello World!")

    def plot_image_crop(image, image_cropped, pix_up, pix_down):
        plt.imshow(image.data)
        plt.gca().add_patch(Rectangle(
            (pix_up[0], pix_down[1]),
            np.abs(pix_up[0]-pix_down[0]),
            np.abs(pix_up[1]-pix_down[1]),
            linewidth=1,edgecolor='r',facecolor='none'
            ))
        plt.show()

        plt.imshow(image_cropped)
        plt.show()