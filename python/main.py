# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:02:51 2023

@author: burro
"""

import cv2
import BZGame
import numpy as np
import matplotlib.pyplot as plt
    
def main(img_name = 'ultrasound.jpg'):
    
    clean = cv2.imread(img_name)
    clean = cv2.cvtColor(clean, cv2.COLOR_BGR2GRAY)
    clean = clean / 255.0
    sigma = 15 #noise level
    #add noise:
    im = clean + np.random.normal(0, sigma/255, clean.shape)
    im = np.clip(im,0,1)
    
    lambda1 = 15
    mu1 = 1
    sigma1 = 2
    xi1 = 100
    lambda2 = 20
    mu2 = 1
    sigma2 = 2
    xi2 = 1000
    lambda3 = 20
    mu3 = 1
    sigma3 = 2
    xi3 = 1000
    
    
    #Run optimisation
    u,v,z = BZGame.game_threeplayer(im,lambda1,mu1,sigma1,xi1,lambda2,mu2,sigma2,xi2,lambda3,mu3,sigma3,xi3)
    
    th = 0.49
    bin_seg = u > th
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.imshow(im,cmap='gray')
    plt.subplot(2,2,2)
    plt.imshow(u,cmap='gray')
    plt.subplot(2,2,3)
    plt.imshow(v, cmap='gray')
    plt.subplot(2,2,4)
    plt.imshow(z, cmap='gray')
    
    plt.figure()
    plt.imshow(bin_seg)
    
    
    
    
if __name__ == "__main__":
    main()
    
    