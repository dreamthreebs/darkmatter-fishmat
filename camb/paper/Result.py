from PIL import Image
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator
from labellines import labelLine, labelLines

minors = (10** np.arange(-4.0,6.0)[:,None] * np.arange(0.1,1,0.1)).ravel()
print(minors)

minors_y = (10** np.arange(-4.0,6.0)[:,None] * np.arange(0.1,1,0.1)).ravel()
print(minors_y)

def Fill_Array(xmin, xmax, nx):
    '''
    Get an array evenly spaced in [xmin, xmax] with nx elements
    ---- inputs ----
    xmin: left boundary
    xmax: right boundary
    nx  : number of x elements
    '''
    dx = (xmax - xmin) / (nx - 1)
    x = np.empty(nx)
    x[0] = xmin
    for idx in np.arange(1, nx):
        x[idx] = x[idx - 1] + dx
    '''if(x[nx-1]>=xmax):
        x[nx-1] = xmax'''
    x[nx - 1] = xmax
    return x


def MCR(xmin, xmax, ymin, ymax, nx, Filename, mode, Smooth_nx,trans):
    '''
    Extract x and y data from a getdist plot image
    ---- inputs ----
    xmin : left boundary of input image
    xmax : right boundary of input image
    ymin : lower boundary of input image
    ymax : upper boundary of input image
    nx   : number of x in extracted data
    Filename : name of input image file
    mode : what you want as outputs
           1 : extracted data, spline interpoaltion
           2 : extracted data, linear interpoaltion
           3 : extracted data, no interpolation
           4 : raw pixel data, useful for debugging
    Smooth_nx : number of x in smoothening downsample, no smoothening if Smooth_nx < 2
    '''

    im = Image.open(Filename)
    width = im.size[0]
    height = im.size[1]
    im = im.convert('RGB')

    # Step 1 : Get pixel locations (x1, y1)
    x1 = np.arange(0, width)
    y1 = np.empty(width)
    for xid in range(0, width):
        r, g, b = im.getpixel((xid, 0))
        yid = 1
        PROCEED = 1
        while PROCEED == 1:
            R, G, B = im.getpixel((xid, yid))
            cf = ((R == r) & (G == g) & (B == b))
            if cf != 0:
                yid = yid + 1
            else:
                PROCEED = 0
        y1[xid] = yid

    # Step 2 : Get physical locations (x2, y2)
    x2 = Fill_Array(xmin, xmax, width)
    y2 = np.empty(width)
    for i in range(0, width):
        y2[i] = (height - y1[i]) * (ymax - ymin) / height + ymin

    # Step 3 : Get low-res samples for smoothened results
    if Smooth_nx < 2:
        x3, y3 = x2, y2
    else:
        x3 = Fill_Array(xmin, xmax, Smooth_nx)
        y3 = interp1d(x2, y2)(x3)

    # Step 4 : interpolate
    x = Fill_Array(xmin, xmax, nx)

    if mode == 1:
        return 10**x, (10 ** (trans)) * interp1d(x3, y3, kind='cubic')(x)
    elif mode == 2:
        return 10**x, (10 ** (trans)) * interp1d(x3, y3)(x)
    elif mode == 3:
        return x2, y2
    elif mode == 4:
        return x1, y1
    else:
        raise Exception("Wrong choice of mode")


#Ay
[x,y]=MCR(-5,3.7,0,2.1,500,'/Users/zhangzixuan/Desktop/Work/MCR/11_ttt.png', 1, 50, -27)# Planck
[x1,y1]=MCR(-5,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/170.png', 1, 50, -28) # AWP, HILC, v1
[x2,y2]=MCR(-5,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/171.png',1, 50, -28) # AWP, HILC, v2
[x3,y3]=MCR(-5,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/172.png',1, 50, -28) # AWP, NILC, v1
[x4,y4]=MCR(-5,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/173.png', 1, 50, -28) # AWP, NILC, v2
[x5,y5]=MCR(-5,3.7,0,101,500,'/Users/zhangzixuan/Desktop/Work/MCR/134.png',1, 50, -29) # PICO, HILC, v1
[x6,y6]=MCR(-5,3.7,0,15,500,'/Users/zhangzixuan/Desktop/Work/MCR/154.png',1, 50, -28) # PICO, HILC, v2
[x7,y7]=MCR(-5,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/146.png',1, 50, -28) # PICO, NILC, v1
[x8,y8]=MCR(-5,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/150.png',1, 50, -28) # PICO, NILC, v2
[x9,y9]=MCR(-5,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/138.png',1, 50, -28) # CVL
x10 = x
y10 = 3 * 10**(-26) * np.ones(len(x10)) / x10

xa_ay=10**3.7
xi_ay=10**(-5)
ya_ay=4*10**(-27)
yi_ay=10**(-28)

#Ae
[x11,y11]=MCR(-3.2,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/15_t.png', 1, 50, -28) # Planck
[x12,y12]=MCR(-3.2,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/178.png', 1, 50, -28) # AWP, HILC, v1
[x13,y13]=MCR(-3.2,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/179.png', 1, 50, -28) # AWP, HILC, v2
[x14,y14]=MCR(-3.2,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/180.png', 1, 50, -28) # AWP, NILC, v1
[x15,y15]=MCR(-3.2,3.7,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/181.png', 1, 50, -28) # AWP, NILC, v2
[x16,y16]=MCR(-3.2,3.7,0,100,500,'/Users/zhangzixuan/Desktop/Work/MCR/131.png', 1, 50, -29) # PICO, HILC, v1
[x17,y17]=MCR(-3.2,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/153.png', 1, 50, -28) # PICO, HILC, v2
[x18,y18]=MCR(-3.2,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/145.png', 1, 50, -28) # PICO, NILC, v1
[x19,y19]=MCR(-3.2,3.7,0,10,500,'/Users/zhangzixuan/Desktop/Work/MCR/149.png', 1, 50, -28) # PICO, NILC, v2
[x20,y20]=MCR(-3.2,3.7,0,7.5,500,'/Users/zhangzixuan/Desktop/Work/MCR/139.png', 1, 50, -28) # CVL
x21 = x11
y21 = 3 * 10**(-26) * np.ones(len(x21)) / x21

xa_ae=10**3.7
xi_ae=10**(-3.2)
ya_ae=3*10**(-27)
yi_ae=9*10**(-29)


#Dy
[x22,y22]=MCR(-5,4,0,12,500,'/Users/zhangzixuan/Desktop/Work/MCR/29.png', 1, 50, -25) # Planck
[x23,y23]=MCR(-5,4,0,80,500,'/Users/zhangzixuan/Desktop/Work/MCR/174.png', 1, 50, -26) # AWP, HILC, v1
[x24,y24]=MCR(-5,4,0,80,500,'/Users/zhangzixuan/Desktop/Work/MCR/175.png', 1, 50, -26) # AWP, HILC, v2
[x25,y25]=MCR(-5,4,0,80,500,'/Users/zhangzixuan/Desktop/Work/MCR/176.png', 1, 50, -26) # AWP, NILC, v1
[x26,y26]=MCR(-5,4,0,80,500,'/Users/zhangzixuan/Desktop/Work/MCR/177.png', 1, 50, -26) # AWP, NILC, v2
[x27,y27]=MCR(-5,4,0,30,500,'/Users/zhangzixuan/Desktop/Work/MCR/132.png', 1, 50, -26) # PICO, HILC, v1
[x28,y28]=MCR(-5,4,0,30,500,'/Users/zhangzixuan/Desktop/Work/MCR/156.png', 1, 50, -26) # PICO, HILC, v2
[x29,y29]=MCR(-5,4,0,25.5,500,'/Users/zhangzixuan/Desktop/Work/MCR/132.png', 1, 50, -26) # PICO, NILC, v1
[x30,y30]=MCR(-5,4,0,30,500,'/Users/zhangzixuan/Desktop/Work/MCR/152_c.png', 1, 50, -26) # PICO, NILC, v2
[x31,y31]=MCR(-5,4,0,19,500,'/Users/zhangzixuan/Desktop/Work/MCR/140.png', 1, 50, -26) # CVL

xa_dy=10**4
xi_dy=10**(-5)
ya_dy=2*10**(-24)
yi_dy=10**(-26)

#De
[x32,y32]=MCR(-2.9,4,0,80,500,'/Users/zhangzixuan/Desktop/Work/MCR/28.png', 1, 30, -26) # Planck
[x33,y33]=MCR(-2.9,4,0,60,500,'/Users/zhangzixuan/Desktop/Work/MCR/182.png', 1, 30, -26) # AWP, HILC, v1
[x34,y34]=MCR(-2.9,4,0,60,500,'/Users/zhangzixuan/Desktop/Work/MCR/183.png', 1, 30, -26) # AWP, HILC, v2
[x35,y35]=MCR(-2.9,4,0,60,500,'/Users/zhangzixuan/Desktop/Work/MCR/184.png', 1, 30, -26) # AWP, NILC, v1
[x36,y36]=MCR(-2.9,4,0,60,500,'/Users/zhangzixuan/Desktop/Work/MCR/185.png', 1, 30, -26) # AWP, NILC, v2
[x37,y37]=MCR(-2.9,4,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/133.png', 1, 30, -26) # PICO, HILC, v1
[x38,y38]=MCR(-2.9,4,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/155.png', 1, 30, -26) # PICO, HILC, v2
[x39,y39]=MCR(-2.9,4,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/147.png', 1, 30, -26) # PICO, NILC, v1
[x40,y40]=MCR(-2.9,4,0,20,500,'/Users/zhangzixuan/Desktop/Work/MCR/151.png', 1, 30, -26) # PICO, NILC, v2
[x41,y41]=MCR(-2.9,4,0,15,500,'/Users/zhangzixuan/Desktop/Work/MCR/141_check.png', 1, 50, -26) # CVL

xa_de=10**4
xi_de=10**(-2.9)
ya_de=10**(-24)
yi_de=2 * 10**(-27)

titlesize = 25
fontsize = 25

plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.default'] = 'regular'

plt.figure(figsize=(40,25))
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

ax1 = plt.subplot(2,2,1)
plt.title(r'$\chi\chi \rightarrow \gamma\gamma$',fontsize=titlesize, pad=10)
ax1.plot(x10, y10, linestyle = ':', color='k',label=r'$Thermal\;\langle\sigma\nu\rangle$')
labelLines(plt.gca().get_lines(), zorder=0, yoffsets=1e-27, xvals=24, size=fontsize, font='Times New Roman')
ax1.plot(x,y,label="Planck", color='k')
ax1.plot(x1,y1,label="Ground observation, HILC, V1", color='r')
ax1.plot(x2,y2,label="Ground observation, HILC, V2", color='r', linestyle='--')
ax1.plot(x3,y3,label="Ground observation, NILC, V1", color='b')
ax1.plot(x4,y4,label="Ground observation, NILC, V2", color='b', linestyle='--')
ax1.plot(x5, y5, label='PICO, HILC, V1', color='g')
ax1.plot(x6, y6, label='PICO, HILC, V2', color='g', linestyle='--')
ax1.plot(x7, y7, label='PICO, NILC, V1', color='m')
ax1.plot(x8, y8, label='PICO, NILC, V2', color='m', linestyle='--')
ax1.plot(x9, y9, label='CVL')
plt.xlabel(r'$m_\chi[GeV]$',fontsize=fontsize)
plt.ylabel(r'$\langle\sigma\nu\rangle/m_\chi\quad[cm^3s^{-1}GeV^{-1}]$',fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xticks(10**np.arange(-5.0,5.0), [r'$10^{-5}$','',r'$10^{-3}$','',r'$10^{-1}$','',r'$10^{1}$','',r'$10^{3}$',''])
plt.xlim(xi_ay,xa_ay)
plt.ylim(yi_ay,ya_ay)
ax1.xaxis.set_minor_locator(FixedLocator(minors))
ax1.tick_params(pad=8, top='on',right='on',which='both', labelsize=fontsize)

ax2 = plt.subplot(2,2,2)
plt.title(r'$\chi\chi \rightarrow e^-e^+$',fontsize=titlesize, pad=10)
ax2.plot(x21, y21, linestyle = ':', color='k',label=r'$Thermal\;\langle\sigma\nu\rangle$')
labelLines(plt.gca().get_lines(), zorder=0, yoffsets=4e-28, xvals=30, size=fontsize, font='Times New Roman')
ax2.plot(x11,y11,label="Planck", color='k')
ax2.plot(x12,y12,label="Ground observation, HILC, V1", color='r')
ax2.plot(x13,y13,label="Ground observation, HILC, V2", color='r', linestyle='--')
ax2.plot(x14,y14,label="Ground observation, NILC, V1", color='b')
ax2.plot(x15,y15,label="Ground observation, NILC, V2", color='b', linestyle='--')
ax2.plot(x16, y16, label='PICO, HILC, V1', color='g')
ax2.plot(x17, y17, label='PICO, HILC, V2', color='g', linestyle='--')
ax2.plot(x18, y18, label='PICO, NILC, V1', color='m')
ax2.plot(x19, y19, label='PICO, NILC, V2', color='m', linestyle='--')
ax2.plot(x20, y20, label='CVL')
plt.xlabel(r'$m_\chi[GeV]$',fontsize=fontsize)
plt.ylabel(r'$\langle\sigma\nu\rangle/m_\chi\quad[cm^3s^{-1}GeV^{-1}]$',fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xlim(xi_ae,xa_ae)
plt.ylim(yi_ae,ya_ae)
ax2.xaxis.set_minor_locator(FixedLocator(minors))
ax2.tick_params(pad=8, top='on',right='on',which='both', labelsize=fontsize)

ax3 = plt.subplot(2,2,3)
plt.title(r'$\chi \rightarrow \gamma\gamma$',fontsize=titlesize, pad=10)
ax3.plot(x22,y22,label="Planck", color='k')
ax3.plot(x23,y23,label="Ground observation, HILC, V1", color='r')
ax3.plot(x24,y24,label="Ground observation, HILC, V2", color='r', linestyle='--')
ax3.plot(x25,y25,label="Ground observation, NILC, V1", color='b')
ax3.plot(x26,y26,label="Ground observation, NILC, V2", color='b', linestyle='--')
ax3.plot(x27, y27, label='PICO, HILC, V1', color='g')
ax3.plot(x28, y28, label='PICO, HILC, V2', color='g', linestyle='--')
ax3.plot(x29, y29, label='PICO, NILC, V1', color='m')
ax3.plot(x30, y30, label='PICO, NILC, V2', color='m', linestyle='--')
ax3.plot(x31, y31, label='CVL')
plt.xlabel(r'$m_\chi[GeV]$',fontsize=fontsize)
plt.ylabel(r'$\Gamma_\chi[s^{-1}]$',fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xticks(10**np.arange(-5.0,5.0), [r'$10^{-5}$','',r'$10^{-3}$','',r'$10^{-1}$','',r'$10^{1}$','',r'$10^{3}$',''])
plt.xlim(xi_dy,xa_dy)
plt.ylim(yi_dy,ya_dy)
ax3.xaxis.set_minor_locator(FixedLocator(minors))
ax3.tick_params(pad=8, top='on',right='on',which='both', labelsize=fontsize)

ax4 = plt.subplot(2,2,4)
plt.title(r'$\chi \rightarrow e^-e^+$',fontsize=titlesize, pad=10)
ax4.plot(x32,y32,label="Planck", color='k')
ax4.plot(x33,y33,label="Ground observation, HILC, V1", color='r')
ax4.plot(x34,y34,label="Ground observation, HILC, V2", color='r', linestyle='--')
ax4.plot(x35,y35,label="Ground observation, NILC, V1", color='b')
ax4.plot(x36,y36,label="Ground observation, NILC, V2", color='b', linestyle='--')
ax4.plot(x37, y37, label='PICO, HILC, V1', color='g')
ax4.plot(x38, y38, label='PICO, HILC, V2', color='g', linestyle='--')
ax4.plot(x39, y39, label='PICO, NILC, V1', color='m')
ax4.plot(x40, y40, label='PICO, NILC, V2', color='m', linestyle='--')
ax4.plot(x41, y41, label='CVL')
plt.xlabel(r'$m_\chi[GeV]$',fontsize=fontsize)
plt.ylabel(r'$\Gamma_\chi[s^{-1}]$',fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.xlim(xi_de,xa_de)
plt.ylim(yi_de,ya_de)
ax4.xaxis.set_minor_locator(FixedLocator(minors))
ax4.tick_params(pad=8, top='on',right='on',which='both', labelsize=fontsize)
plt.legend(loc='lower right', fontsize=fontsize)

plt.savefig('/Users/zhangzixuan/Desktop/Results.pdf', bbox_inches='tight')
plt.show()