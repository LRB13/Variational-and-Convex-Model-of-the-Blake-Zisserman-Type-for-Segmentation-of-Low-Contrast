# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:06:15 2023

@author: burro
"""
import numpy as np

def game_threeplayer(im,lambda1,mu1,sigma1,xi1,lambda2,mu2,sigma2,xi2,lambda3,mu3,sigma3,xi3):
    
    n,m = im.shape
    
    u = im.copy()
    wx = np.zeros(im.shape)
    wy = np.zeros(im.shape)
    bx = np.zeros(im.shape)
    by = np.zeros(im.shape)
    
    v = np.zeros(im.shape)
    w2x = np.zeros(im.shape)
    w2y = np.zeros(im.shape)
    b2x = np.zeros(im.shape)
    b2y = np.zeros(im.shape)
    
    z = np.zeros(im.shape)
    w3x = np.zeros(im.shape)
    w3y = np.zeros(im.shape)
    b3x = np.zeros(im.shape)
    b3y = np.zeros(im.shape)
    
    maxit = 200
    
    #imx, imy = Fx(im), Fy(im)
    imx, imy = np.gradient(im)
    
    #Build kernels for FFT algorithm (5 point stencil)
    uker = np.zeros(im.shape)
    uker[0,0] = 4
    uker[0,1] = -1
    uker[1,0] = -1
    uker[-1,0] = -1
    uker[0,-1] = -1
    
    uker1 = lambda1+(mu1+sigma1 + xi1)*np.fft.fft2(uker)
    uker2 = lambda2 + xi2 +(mu2+sigma2)*np.fft.fft2(uker)
    uker3 = lambda3 + xi3 +(mu3+sigma3)*np.fft.fft2(uker)
    
    lhs_u = lambda1 + 4*(mu1+xi1+sigma1)
    lhs_u = lambda1 + 4*sigma1
    for k in range(maxit):
        
        #u problem: Solve using GaussNewton..
        # rhs = lambda1*im - xi1*Bx(v) - xi1*By(z) - sigma1*Bx(wx-bx) - sigma1*By(wy-by) + (mu1 + xi1 + sigma1)*(Fx(Bx(u)) + Fy(By(u)) + 4*u)
        # for k2 in range(5):
        #     rhs = lambda1*im - sigma1*Bx(wx-bx) - sigma1*By(wy-by) + sigma1*(Bx(Fx(u))+By(Fy(u)) + 4*u)
        #     u = rhs/lhs_u
        oldu = u.copy()
        
        rhs1 = lambda1*im - sigma1*Bx(wx-bx) - sigma1*By(wy-by) - xi1*Bx(v) + xi1*By(z)
        u = np.fft.ifft2(np.fft.fft2(rhs1)/uker1)
        u = np.real(u)
        
        err = np.linalg.norm(u-oldu)/np.linalg.norm(u)
        if err < 1e-4:
            break
        
        
        #wx, wy:
        temp1 = Fx(u) + bx
        temp2 = Fy(u) + by
        wx,wy = shrink2(temp1,temp2,1/sigma1)
        
        bx = bx + Fx(u) - wx
        by = by + Fy(u) - wy
        
        #ux, uy = Fx(u), Fy(u)
        ux, uy = np.gradient(u)
        
        rhs2 = lambda2*imx + xi2*ux - sigma2*Bx(w2x-b2x) - sigma2*By(w2y-b2y)
        v = np.fft.ifft2(np.fft.fft2(rhs2)/uker2)
        v = np.real(v)

        
        rhs3 = lambda3*imy + xi3*uy - sigma3*Bx(w3x-b3x) - sigma3*By(w3y-b3y)
        z = np.fft.ifft2(np.fft.fft2(rhs3)/uker3)
        z = np.real(z)
        
        #w2x, w2y
        temp1=Fx(v)+b2x
        temp2=Fy(v)+b2y
        w2x,w2y =shrink2(temp1,temp2,1/sigma2)
        #bproblem
        b2x = b2x + Fx(v) - w2x
        b2y = b2y + Fy(v) - w2y
        
        #w3x,w3y problem
        
        temp1=Fx(z)+b3x
        temp2=Fy(z)+b3y
        w3x,w3y =shrink2(temp1,temp2,1/sigma3)
        #bproblem
        b3x = b3x + Fx(z) - w3x
        b3y = b3y + Fy(z) - w3y
 
        
    
    return u,v,z


def circshift(x, amount, axis):
    return np.roll(x,amount,axis=axis)

#Finite differences
def Fx(u):
    return circshift(u,-1,1) - u
def Fy(u):
    return circshift(u,-1,0) - u
def Bx(u):
    return - circshift(u,1,1) + u
def By(u):
    return -circshift(u,1,0) + u

def shrink2(sx,sy,sig):
    s = np.sqrt(sx**2 + sy**2 + 1e-9)
    t = np.maximum(s-sig,0)
    wx = t*(sx/s)
    wy = t*(sy/s)
    return wx,wy
