#!/usr/bin/env python3
import os
#from pathlib import path
import copy
import math
from collections import defaultdict
import argparse
import pandas as pd
from collections import defaultdict
from itertools import combinations
import numpy as np
import subprocess
import matplotlib
# Use the pgf backend (must be set before pyplot imported)
matplotlib.use('pgf')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoLocator
from matplotlib.collections import PathCollection

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from datetime import datetime
import glob
import logging
cm = 1/2.54  # centimeters in inches
ratio_edge_width=0.35/np.sqrt(np.pi)
log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(filename)s:%(lineno)04d - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
log.handlers=[console_handler]

norm_phase = colors.Normalize(vmin=-np.pi, vmax=np.pi)
zorder_scatter=3#either 2 or 3 to be below or above gridlines and labels
colours = [(0, 1, 0), (1, 0, 0)] # first color is green, last is red
cmapgaugeerrors= LinearSegmentedColormap.from_list("Custom", colours, N=20)
angleforradiallabel={**{"goldentriangle":112,"p31mtriangle":30},**{i:180 for i in ['p31mkites','p31mkitesturn','p6mhexa','p6mkites']}}
angularticks={**{"goldentriangle":[0,0.9046933/2,0.9046933,np.pi/2,(np.pi/2+np.pi-0.666131)/2,np.pi-0.666131,np.pi,np.pi+0.9046933/2,np.pi+0.9046933,3/2*np.pi,(3/2*np.pi+2*np.pi-0.666131)/2,2*np.pi-0.666131],"penrose-centered":np.linspace(0,2*np.pi*9/10,10),"penrose":np.linspace(0,2*np.pi*9/10,10),"ammann":np.linspace(0,2*np.pi*7/8,8),"shield":np.linspace(0,2*np.pi/12*11,12),"dodecagonalqp":np.linspace(0,2*np.pi/12*11,12)},**{i:np.linspace(0,2*np.pi/6*5,6) for i in ['p31mtriangle','p31mkites','p31mkitesturn','p6mhexa','p6mkites','hexagonalqp','question_christelle_1','question_christelle_2']}}
colorsampl = [(0, 'blue'), (0.5, (0,1,0)), (1, 'red')]
cmap_amplitudeerrors = cmapgaugeerrors#LinearSegmentedColormap.from_list('blue_green_red', colorsampl)
offset_to_annotations={("fibonaccisquares",3):np.array((-4,0)),("fibonaccisquares",4):np.array((0,-4))}
def add_pixels(parser,args):
    parser.add_argument('--pixels', type=int,nargs=2,default={"small":[1001,1001],"medium":[3201,3201],"large":[6401,6401],"xl":[6401,6401],"xL":[6401,6401],"xxl":[6401,6401],"xxL":[6401,6401],"xxxl":[10001,10001],"xxxL":[10001,10001],"xxxxl":[10001,10001]}[args.size],help='two numbers describing width and height')
def parse_cmd_line():
    parser = argparse.ArgumentParser(description='plot fourier transform of given pattern')
    parser.add_argument('--file', type=str,help='file containing studied grayscale image')
    parser.add_argument('--titlesize',help='fontsize to set title of plots' ,type=str,default='small')
    parser.add_argument('--intervall',help='how many pixels to include in diffraction diagram: square region from -intervall to +intervall' ,type=int)
    parser.add_argument('--studywaveveccolor',help='color for selected subset of integer multiples of basis vectors in overlay plot' ,type=str,default='gray')
    parser.add_argument('--threshold',help='threshfold of amplitude of a peak in order to be shown in a diffraction diagramm' ,type=float,default=0.01)
    parser.add_argument('--handletextpad',help='spacing inside the legend',type=float,default=0)
    parser.add_argument('--ylabelpad',help='padding for k_2 y-label of diffraction diagram plots' ,type=float,default=-2.)
    parser.add_argument('--nbins',help='number of bins for histogram' ,type=int,default=20)
    parser.add_argument('--thresholdgauge',help='threshfold of amplitude of a peak in order to be considered for gauge calculations' ,type=float,default=0.01)
    parser.add_argument('--legendfontsize',help='fontsize of legend' ,type=str,default='x-small')
    parser.add_argument('--customticks',action='store_true',help='wether to use ticks indicated by plot_fftpattern.py:tickmarkers or not')
    parser.add_argument('--figurewidth',type=float,help='width of figure in cm',default=9)
    parser.add_argument('--figureheight',type=float,help='height of figure in cm',default=9)
    parser.add_argument('--figurewidthhisto',type=float,help='width of figure in cm',default=11)
    parser.add_argument('--figureheighthisto',type=float,help='height of figure in cm',default=7.5)
    parser.add_argument('--histogramvlinerotations',type=float,help='position of red vertical line drawn in histogram of rotations',default=0.114)#source: p4m_0.02_invert/large
    parser.add_argument('--histogramvlinemirrors',type=float,help='position of red vertical line drawn in histogram of mirrors',default=0.03)#source: p4m_0.02_invert/large
    parser.add_argument('--scalesize',type=float,help='scales size of diffraction peaks',default=1/2)
    parser.add_argument('--circles',type=float,help='radii of circles to add on diffraction image',nargs='*',default=[])
    parser.add_argument('--polarticks',type=float,help='ticks in polar direction',nargs='*')
    parser.add_argument('--nolegend',type=str,help='specifies which plots shall not have a legend attached to them',default=[],action='append',choices=['amplitudeerrors','gaugeerrors','gauges','diffractionfft','diffractionexact','combinederrors'])
    parser.add_argument('--noyticks',type=str,help='specifies which plots shall not have a y-ticks (radial)',default=[],action='append',choices=['amplitudeerrors','gaugeerrors','gauges','diffractionfft','diffractionexact','combinederrors'])
    parser.add_argument('--rank',help='rank of lattice up to which we verify gauge linearity' ,type=int,default=4)
    parser.add_argument('--fouriertitle',help='title to place on diffractiondiagrams',type=str,default=r'Fourier transform $\hat{\rho}(\ve{k})$ at considered wavevectors $\latticefourier{}$')
    parser.add_argument('--fontsizelabels',help='fontsize for labels and annotations',type=int,default=10)
    parser.add_argument('--fontsizeticks',help='fontsize for ticks',type=int,default=9)
    parser.add_argument('--annotatediffraction',help='wether to indicate basis vectors in diffraction diagramms',action='store_true')
    parser.add_argument('--noannotategauge',help='wether to indicate basis vectors in gauge plots',action='store_true')
    parser.add_argument('--radiallabels',help='do not plot radial labels',action='store_true')
    parser.add_argument('--annotatestyle',help='choose arrows if you want to indicate basevectors using arrows',type=str,default='circle',choices=['circle','arrow'])
    parser.add_argument('--plottype',help='fontsize for ticks',type=str,default='polar',choices=['polar','cartesian'])
    args = parser.parse_args()
    return args
def plot_figure(fig,basepath,close=True,transparent=False):
    log.info(f'exporting {basepath}')
    fig.savefig(basepath+'.png',bbox_inches="tight",transparent=transparent,dpi=600)
#    fig.savefig(basepath+'.pdf',bbox_inches="tight",transparent=transparent)
    if os.environ.get('Cluster',0)!='1' and False:
        for ax in fig.get_axes():
            ax.set_title(ax.get_title().replace('_',r'\_'))
            ax.set_xlabel(ax.get_xlabel().replace('_',r'\_'))
            ax.set_ylabel(ax.get_ylabel().replace('_',r'\_'))
            if ax.get_legend() is not None:
                for text in ax.get_legend().get_texts():
                    text.set_text(text.get_text().replace('_',r'\_'))
    if close:
        plt.close(fig)

def read_gaugeerror_file(args,ge_file):
    data=np.genfromtxt(ge_file,delimiter=',')
    symmetryname=' '.join(ge_file.split('_')[3:])[:-4]
    basename=ge_file.replace('.csv','')
    symmetrysign=r'\mathbf{M}' if 'mirror' in symmetryname else r'\mathbf{R}'

    wavevectors=data[:,0:2]
    gauges=data[:,2]
    gaugeerrors=np.abs(data[:,3])
    amplitudes=data[:,4]
    amplitudeerrors=np.abs(data[:,5])
    return symmetryname,basename,symmetrysign,wavevectors,gauges,gaugeerrors,amplitudes,amplitudeerrors
def round_first_sig(x):
    if x == 0:
        return 0
    n = math.floor(math.log10(abs(x)))
    factor = 10**n
    return round(x / factor) * factor

def calculate_histogram_spacing(largest_value,nbins_approx=20):
    binsize=round_first_sig(largest_value/nbins_approx)
    nr_bins=math.ceil(largest_value/binsize)
    return binsize,nr_bins
def calculate_ticks_histogram(binsize,nr_bins):
    major_ticks=[i*binsize for i in range(0,nr_bins,4)]
    minor_ticks=[i*binsize for i in range(2,nr_bins,4)]
    return major_ticks,minor_ticks

def plot_gaugeerrors_histogram(args,ge_file):
    symmetryname,basename,symmetrysign,wavevectors,gauges,gaugeerrors,amplitudes,amplitudeerrors=read_gaugeerror_file(args,ge_file)
    basename=basename+'_histogram'
    if 'mirror' in symmetryname:
        #filtering points left invariant by mirror
        vline_pos=args.histogramvlinemirrors
    else:
        vline_pos=args.histogramvlinerotations
    plot_histogramm(basename.replace('gauge','amplitude'),amplitudeerrors,symmetryname,symmetrysign,'amplitudeerrors',args.nbins,(args.figurewidthhisto*cm,args.figureheighthisto*cm))
    plot_histogramm(basename,gaugeerrors,symmetryname,symmetrysign,'gaugeerrors',args.nbins,figsize=(args.figurewidthhisto*cm,args.figureheighthisto*cm),limitx=0.5)

    plot_histogramm(basename.replace('gauge','combined'),1-(1-2*gaugeerrors)*(1-amplitudeerrors),symmetryname,symmetrysign,'combinederrors',args.nbins,figsize=(args.figurewidthhisto*cm,args.figureheighthisto*cm),histogramvline=vline_pos)
    return 
def plot_histogramm(basename,errors,symmetryname,symmetrysign,typ,nbins_approx,figsize,limitx=1,histogramvline=None):
    if (np.sum(errors)==0):
        return
    print(figsize)
    fig,ax=plt.subplots(figsize=figsize)

    maxval=np.max(errors) if histogramvline is None else max(np.max(errors),histogramvline)
    (binsize,nr_bins)=calculate_histogram_spacing(maxval,nbins_approx)
    (major_ticks,minor_ticks)=calculate_ticks_histogram(binsize,nr_bins)

    ax.hist(errors,[i*binsize for i in range(nr_bins+1)],edgecolor='black')

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks,minor=True)

    if typ=="gaugeerrors":
        ax.set_xlabel(r'$\Delta\Phi_{'+symmetrysign+'}$',fontsize='x-large')
    elif typ=="amplitudeerrors":
        ax.set_xlabel(r'$\Delta \left|\hat{\rho}\right|$',fontsize='x-large')
    elif typ=="combinederrors":
        ax.set_xlabel(r"$\Delta^{\mathrm{sym}}$",fontsize='x-large')

    if histogramvline is not None:
       ax.axvline(histogramvline,color='red') 

       ax.text(histogramvline, 0.95, f'threshold={histogramvline}', color='r', ha='right', va='top',rotation=90,
            transform=ax.get_xaxis_transform(),size='x-large')
    ax.set_title(f'Deviation measures of {len(errors)} analyzed frequencies')
    plot_figure(fig,basename,close=True)

def plot_gaugeerrors(ge_file,args,maxsize,number_bases,diffractiondata):
    symmetryname,basename,symmetrysign,wavevectors,gauges,gaugeerrors,amplitudes,amplitudeerrors=read_gaugeerror_file(args,ge_file)

    plot_gaugefunction_spatial(wavevectors,amplitudes,amplitudeerrors,args,symmetryname,symmetrysign,basename.replace('gauge','amplitude')+'_spatial',maxsize,number_bases,typ='amplitudeerrors',diffractiondata=diffractiondata)
    plot_gaugefunction_spatial(wavevectors,amplitudes,np.mod(gauges,1),args,symmetryname,symmetrysign,basename.replace('error','')+'_spatial',maxsize,number_bases,diffractiondata=diffractiondata)
    plot_gaugefunction_spatial(wavevectors,amplitudes,gaugeerrors,args,symmetryname,symmetrysign,basename+'_spatial',maxsize,number_bases,typ='gaugeerrors',diffractiondata=diffractiondata)
    plot_gaugefunction_spatial(wavevectors,amplitudes,1-(1-2*gaugeerrors)*(1-amplitudeerrors),args,symmetryname,symmetrysign,basename.replace('gauge','combined')+'_spatial',maxsize,number_bases,typ='combinederrors',diffractiondata=diffractiondata)

    if False:
        fig,ax=plt.subplots()
        ax.scatter(amplitudes,gaugeerrors)
        ax.set_xlabel(r'$\left|{\hat{\rho}\right|(\ve{k})$')
        ax.set_xscale('log')
        ax.set_ylabel(r'$\Delta\Phi_{'+symmetrysign+r'}(\ve{k})$')
        #ax.set_title(f'Deviation of phase function for {symmetryname}')
        ax.axhline(y=0, color='gray', linestyle='--')
        plot_figure(fig,basename)
        plt.close(fig)

        fig,ax=plt.subplots()
        ax.axhline(y=0, color='black', linestyle='--')
        ax.scatter(amplitudes,gaugeerrors)
        ax.set_xlabel(r"$\frac{\left|{\left|\hat{\rho}\right|(\ve{k})-\left|{\hat{\rho}}\right|("+symmetrysign+r"\ve{k})}\right|}{\max(\left|{\hat{\rho}}\right|(\ve{k}),\left|{\hat{\rho}}\right|("+symmetrysign+r"\ve{k}))}$")
        ax.set_ylabel(r'$\Delta\Phi_{'+symmetrysign+r'}(\ve{k})$')
        #ax.set_title(f'Deviation of phase function  function for {symmetryname}')
        plot_figure(fig,basename+'amplitudeerrors')
        plt.close(fig)
def sizefun(maxamplitude,args,lim):
    scale_size=lambda s: 3000*args.scalesize*(args.figureheight*args.figurewidth*(s/maxamplitude)/lim**2)**0.7
    scale_size=lambda s: 100*args.scalesize*(s)**0.7
    return scale_size

def plot_fourier_at_idces(args,basename,maxamplitude):
    diffractiondatabase=np.genfromtxt(f"fourier_coefficients_base_{args.file.replace('.png','')}.csv",delimiter=',')
    diffractiondataother=np.genfromtxt(f"fourier_coefficients_other_{args.file.replace('.png','')}.csv",delimiter=',')
    lim=lim_of_diff_plot(args,np.vstack((diffractiondatabase[:,:2],diffractiondataother[:,:2])))
    scale_size=sizefun(maxamplitude,args,lim)
    settings_basevectors= {'facecolors':plt.cm.hsv(norm_phase(diffractiondatabase[:,3])),'edgecolors':'black','linewidths':np.sqrt(scale_size(diffractiondatabase[:,2]))*ratio_edge_width} if args.annotatediffraction else {'facecolors':plt.cm.hsv(norm_phase(diffractiondatabase[:,3])),'linewidths':0}
    log.debug(f'settings_basevectors={settings_basevectors}')

    plot_diffraction_image([
        [diffractiondataother,{'facecolors':plt.cm.hsv(norm_phase(diffractiondataother[:,3])),'edgecolors':plt.cm.hsv(norm_phase(diffractiondataother[:,3])),'linewidths':0}],
        [diffractiondatabase,settings_basevectors]],
        args,
        basename,maxamplitude,circle_radii=args.circles,title=args.fouriertitle,do_legend=("diffractionexact" not in args.nolegend),annotate_bases=args.annotatediffraction)
def lim_of_diff_plot(args,wavevecs):
    if args.plottype=='cartesian':
        return args.intervall #np.max(wavevecs)*1.1 if plotsize shall change depending on actually plottet points
    elif args.plottype=='polar':
        return np.max(np.sqrt(wavevecs[:,0]**2 + wavevecs[:,1]**2)) #np.max(wavevecs)*1.1 if plotsize shall change depending on actually plottet points
def transform_wavevecs(args,wavevecbasex,wavevecbasey):
    if args.plottype=='cartesian':
        return wavevecbasex,wavevecbasey
    elif args.plottype=='polar':
        return  np.arctan2(wavevecbasey, wavevecbasex),np.sqrt(wavevecbasex**2 + wavevecbasey**2)
def plot_considered_vectors(args,wavevecbase,wavevecother,ax,colors=['gray','black'],label1=r'Fundamental frequencies $\mathbf{k}_i$',label2=r'Their integer multiples $\mathcal{F}^{\mathrm{s}}$'):
    scatter_base=ax.scatter(*transform_wavevecs(args,wavevecbase[:,0],wavevecbase[:,1]),s=6*args.scalesize,c=colors[1],label=label1,zorder=zorder_scatter)
    scatter_other=ax.scatter(*transform_wavevecs(args,wavevecother[:,0],wavevecother[:,1]),s=6*args.scalesize,color=colors[0],linewidths=0.5,zorder=zorder_scatter,label=label2)
    return scatter_base,scatter_other
def plot_peaks_and_fft(args,wavevecbase,wavevecother,diffractiondata):
    lim=lim_of_diff_plot(args,wavevecother)
    fig,ax=figax_for_diffraction(args,lim,doylabel=True)
    if args.plottype=='polar':
        ax.set_yticklabels([])
        set_angularticks(ax,args)
    ax.grid(linewidth=0.3)

    scatter_other,scatter_base=plot_considered_vectors(args,wavevecbase,wavevecother,ax,[args.studywaveveccolor,args.studywaveveccolor],label1=None,label2=None)
    plot_figure(fig,f'fourier_module_{args.file.replace(".png","")}',close=False)
    #overlay with diffractiondata
    minamplitude=args.threshold
    maxamplitude=max(diffractiondata[:,2])

    lim=lim_of_diff_plot(args,diffractiondata[:,0:2])
    set_limit_diffraction(ax,lim,args,yticks=True)

    sizes=diffractiondata[:,2]
    scale_size=sizefun(maxamplitude,args,lim)
    sc=ax.scatter(*transform_wavevecs(args,diffractiondata[:,0],diffractiondata[:,1]),s=scale_size(sizes),facecolors='none',linewidths=0.15,zorder=zorder_scatter,edgecolors='k')

    leg=add_legend(args,diffractiondata[:,2],scale_size,ax,legendfacecolor='none')

    leg2 = plt.legend([scatter_other],[r"$\mathcal{M}$"],handletextpad=0,borderpad=0.15)
    #leg2 = plt.legend([scatter_other],[r"$\mathcal{L}$"],handletextpad=0,borderpad=0.15)
    ax.add_artist(leg)
    leg.set_loc('upper left')
    leg.set_bbox_to_anchor((0, 1))
    ax.add_artist(leg2)
    leg2.set_loc('upper right')
    leg2.set_bbox_to_anchor((1, 1))

    plot_figure(fig,f'overlay_considered_vectors_diffraction_{args.file.replace(".png","")}',transparent=False)
def plot_cbar(sc,cmap,args,cax=None,ax=None,title=None,cbar_limits=[-np.pi, np.pi],pad=0): 
    cbar = plt.colorbar(cmap,ax=ax, cax=cax, orientation='vertical',pad=pad)
    if title is  not None:
        #cbar.set_label(title)
        cbar.ax.set_title(title)
    if cbar_limits[0] == -np.pi and cbar_limits[1] == np.pi:
        cbar.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        cbar.set_ticklabels([r'\scalebox{0.6}{$-\pi$}', r'\scalebox{0.6}{$-\frac{\pi}{2}$}', r'\scalebox{0.6}{0}', r'\scalebox{0.6}{$\frac{\pi}{2}$}', r'\scalebox{0.6}{$\pi$}'],fontsize=args.fontsizeticks*1.3)
    else:
        tick_positions = [round(i,1) for i in cbar.get_ticks()]
        cbar.set_ticks(tick_positions)
        cbar.ax.set_yticklabels([f'{i}' for i in tick_positions],fontsize=args.fontsizeticks)
    sc.set_clim(cbar_limits)
def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    for ax in fig.get_axes():
        ax.axis('off')
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)
    #legend.set_bbox_to_anchor(old_anchor)
    for ax in fig.get_axes():
        ax.axis('on')
def add_legend(args,amplitudes,scale_size,ax,legendtitle=r'\hat{\rho}',legendfacecolor='black'):
    legend_sizes = np.geomspace(np.min(amplitudes), np.max(amplitudes), 4)
    legend_scaled_sizes=scale_size(legend_sizes)
    for legend_size, legend_scaled_size in zip(legend_sizes, legend_scaled_sizes):
       ax.scatter([], [], s=legend_scaled_size, alpha=1.0, label=f'{legend_size:.2g}',facecolors=legendfacecolor, edgecolors='k',linewidths=0.15)
    legend=ax.legend(title=r'$\left|{'+legendtitle+r'}\right|$',handletextpad=args.handletextpad,bbox_to_anchor=(0,1),loc='upper center',fontsize=args.legendfontsize)
    return legend

def plot_gaugefunction_spatial(wavevecs,amplitudes,gauges,args,symmetryname,symmetrysign,basename,maxsize,number_bases,typ='gauges',diffractiondata=None):
    #gauges are interpreted either as errors from actual and propagated gauge or as actual gauge values according to is_errors 
    cm = 1/2.54  # centimeters in inches
    dir='.'
    #Fourier space
    minamplitude=np.min(amplitudes)
    maxamplitude=maxsize
    log.debug(f'len(wavevecs)={len(wavevecs)}') 
    lim=lim_of_diff_plot(args,wavevecs)
    fig,ax=figax_for_diffraction(args,lim,doylabel=(typ not in args.noyticks))
    scale_size=sizefun(maxamplitude,args,lim)

    if typ=="gaugeerrors":
        colours = [(0, 1, 0), (1, 0, 0)] # first color is green, last is red
        cmap = cmapgaugeerrors
        cbar_limits=[0,0.5]
        title=r'$\Delta\Phi_{'+symmetrysign+'}$'
    elif typ=="gauges":
        cmap=plt.cm.twilight
        cbar_limits=[0,1.0]
        title=r'$\Phi_{'+symmetrysign+'}$'
    elif typ=="amplitudeerrors":
        cmap=cmap_amplitudeerrors
        cbar_limits=[0.0,1.0]
        title=r"$\Delta \left|{\hat{\rho}}\right|$"
    elif typ=="combinederrors":
        cmap=cmap_amplitudeerrors
        cbar_limits=[0.0,1.0]
        title=r"$\Delta^{\mathrm{sym}}$"

    for rad in args.circles:
        circle = plt.Circle((0, 0), rad, color='gray', fill=False, linewidth=1,transform=ax.transData._b)
        ax.add_artist(circle)

    norm = plt.Normalize(*tuple(cbar_limits))


    opts={"facecolors":cmap(norm(gauges[:number_bases])),"linewidths":np.sqrt(scale_size(amplitudes[:number_bases]))*ratio_edge_width,"edgecolors":(cmap(norm(gauges[:number_bases])) if args.noannotategauge else 'black')}
    #plot non-base vectors
    sc2=ax.scatter(*transform_wavevecs(args,wavevecs[number_bases:,0],wavevecs[number_bases:,1]),s=scale_size(amplitudes[number_bases:]),c=gauges[number_bases:],cmap=cmap,norm=norm,linewidths=0,zorder=zorder_scatter)
    #plot basevectors
    sc1=ax.scatter(*transform_wavevecs(args,wavevecs[0:number_bases,0],wavevecs[0:number_bases,1]),s=scale_size(amplitudes[0:number_bases]),zorder=zorder_scatter,**opts) 

    if diffractiondata is not None:
        sizes=diffractiondata[:,2]
        scale_size=sizefun(maxamplitude,args,lim)
        sc=ax.scatter(*transform_wavevecs(args,diffractiondata[:,0],diffractiondata[:,1]),s=scale_size(sizes),facecolors='none',linewidths=0.15,zorder=zorder_scatter,edgecolors='k')

    if not args.noannotategauge:
        annotate_basevecs(args,ax,wavevecs[:number_bases,0:2])
    sc=PathCollection(sc1.get_paths() + sc2.get_paths(),np.concatenate((sc1.get_sizes(), sc2.get_sizes())))

    if typ not in args.nolegend:
        add_legend(args,amplitudes,scale_size,ax)#we need to add legend before setting the limit, otherwise the polar plots become distorted
        log.info(f'adding legend to {basename}')
    set_limit_diffraction(ax,lim,args,yticks=(None if (typ not in args.noyticks) or len(args.circles)==0 else args.circles))

    cax = get_cax(args,fig,ax)
    plot_cbar(sc,plt.cm.ScalarMappable(norm,cmap = cmap),args,cbar_limits=cbar_limits,cax=cax,ax=ax,title=title) 
    if args.plottype=='polar':
        fig.tight_layout()
    plot_figure(fig,os.path.join(dir,basename),close=True)
    plt.close(fig)
def annotate_basevecs(args,ax,basevecs):
    for (ib,basevec) in enumerate(basevecs):
        offset=np.dot(np.array([[np.cos(np.pi/4),-np.sin(np.pi/4)],[np.sin(np.pi/4),np.cos(np.pi/4)]]),np.array([basevec[0],basevec[1]]))/np.linalg.norm(basevec)*10
        ax.annotate(r'$\ve{k}_'+str(ib+1)+'$',transform_wavevecs(args,basevec[0],basevec[1]),fontsize=args.fontsizelabels,weight='bold',xytext=offset,textcoords='offset points',ha='center',va='center',zorder=5)
        if args.annotatestyle=='arrow':
            ax.arrow(0, 0, *transform_wavevecs(args,basevec[0],basevec[1]))
def get_cax(args,fig,ax):
    if args.plottype=='cartesian':
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
    else:
        cax=fig.add_axes([0.97,0.15,0.03,0.7])
    return cax
def figax_for_diffraction(args,lim,doylabel=True):
    if args.plottype =='cartesian':
        fig,ax= plt.subplots(figsize=(args.figurewidth*cm,args.figureheight*cm))
        if doylabel:
            ax.set_ylabel(r'$k_y$',fontsize=args.fontsizelabels,labelpad=args.ylabelpad)
        else:
            ax.set_yticks([])
        ax.set_xlabel(r'$k_x$',fontsize=args.fontsizelabels)
        ax.tick_params(axis='x', labelsize=0.8*args.figurewidth)
        ax.tick_params(axis='y', labelsize=0.8*args.figureheight)
        ax.set_aspect('equal')
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)
    else:
        fig,ax= plt.subplots(figsize=(args.figurewidth*cm,args.figureheight*cm),subplot_kw={'projection': 'polar'})
    return fig,ax
def set_angularticks(ax,args):
    xticks_radian=args.polarticks if args.polarticks is not None else [2*np.pi/8*i for i in range(8)]
    xticks_deg=[round(i*180/np.pi,1) for i in xticks_radian]
    xticks_radian=[i/180*np.pi for i in xticks_deg]
    ax.set_xticks(xticks_radian)
    ax.set_xticklabels([str((int(i) if i%1==0 else i))+'°' for i in xticks_deg],fontsize='medium')
    ax.tick_params(axis='x',pad=-3)
    for label in ax.get_xticklabels():
        angle = label.get_position()[0]  # Get the angle position of the label
        
        label.set_position((angle, label.get_position()[1] * 0.9))  # Scale radius to bring labels closer
def set_limit_diffraction(ax,lim,args,yticks=True):
    if args.plottype =='cartesian':
        ax.set_axisbelow(True)
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)
        ax.xaxis.set_major_locator(AutoLocator())
        x_ticks = ax.get_xticks()

        # Set x-ticks as usual and apply these positions for grid lines on both axes
        ax.set_yticks(x_ticks)

        ax.set_ylim(-lim,lim)
        ax.grid(zorder=-100, linestyle='--', linewidth=0.5)
        if yticks is None:
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_visible(False)
                tick.tick2line.set_visible(False)
                tick.label1.set_visible(False)
                tick.label2.set_visible(False)

    else:
        set_angularticks(ax,args)
        ax.set_axisbelow(True)
        rmax=lim*1.05
        ax.set_rmax(rmax)
        numticks=5
        #yticks control radial ticks it seems
        delta=round(rmax/(numticks+1))

        yticks=np.around(np.linspace(delta,delta*numticks,numticks)) 
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(int(i)) for i in yticks], verticalalignment='center', horizontalalignment='center',fontsize='x-small',zorder=3000)
        for label in ax.get_yticklabels():
            label.set_zorder(300)  # Set zorder higher than the scatter plot if needed
        ax.set_rlabel_position(-22.5)
        ax.grid(zorder=-100, linestyle='--', linewidth=0.5)
        ax.set_axisbelow
        if not args.radiallabels:
            ax.set_yticklabels([])
        print(f'xticks={ax.get_xticks()}')
def plot_diffraction_image(datas,args,diffraction_name,maxsize,annotate_bases=False,circle_radii=[],title='',legendtitle=r'\hat{\rho}',do_legend=True,legendfacecolor='black',do_yticks=True):
    #datas=list of pairs data and keywords to plot data, keywords should include facecolors,edgecolors,linewidths
    #data contains 4 columns, x,y,absolute value and phase of peaks
    cm = 1/2.54  # centimeters in inches
    dir='.'
    #Fourier space
    alldata=np.concatenate([i[0] for i in datas], axis=0)
    minamplitude=args.threshold
    maxamplitude=maxsize
    lim=lim_of_diff_plot(args,alldata[:,0:2])
    scale_size=sizefun(maxamplitude,args,lim)
    for colour in [True,False]:
        basename=diffraction_name+'_'+args.file.replace('.png','')+f'{"amplitudes" if colour==False else ""}'
        fig,ax=figax_for_diffraction(args,lim,doylabel=do_yticks)
        for rad in circle_radii:
            circle = plt.Circle((0, 0), rad, color='gray', fill=False, linewidth=1,transform=ax.transData._b)
            ax.add_artist(circle)

        #ax.set_title(title,fontsize=args.titlesize)
        #ax.grid(True, which='both', axis='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        for idat,(data,options) in enumerate(reversed(datas)):#we reverse so that basevectors and their annotation are plotted on top
            wavevecs=data[:,0:2]
            sizes=data[:,2]
            phases=data[:,3]
            if colour==False:
                optionstouse=options if ('facecolors' not in options or options['facecolors'] is 'none') else {**options,**{"facecolors":'k'}}
                #if linewidths is already given, it stays the same
                optionstouse={**{"linewidths":np.sqrt(scale_size(sizes))*ratio_edge_width},**optionstouse}
                ax.scatter(*transform_wavevecs(args,wavevecs[:,0],wavevecs[:,1]),s=scale_size(sizes),**{**optionstouse,**{"edgecolors":'k'}},zorder=zorder_scatter)
            else:
                optionstouse={**{"linewidths":np.sqrt(scale_size(sizes))*ratio_edge_width},**options}
                sc=ax.scatter(*transform_wavevecs(args,wavevecs[:,0],wavevecs[:,1]),s=scale_size(sizes),**optionstouse,zorder=zorder_scatter)
                if idat==0 and annotate_bases:
                    annotate_basevecs(args,ax,wavevecs[:,0:2])
        if do_legend:#we need to add legend before setting the limit, otherwise the polar plots become distorted
            log.info('adding legend')
            add_legend(args,alldata[:,2],scale_size,ax,legendfacecolor=legendfacecolor)
        set_limit_diffraction(ax,lim,args,yticks=(None if not do_yticks or len(args.circles)==0 else args.circles))

        if colour==True:
            cax=get_cax(args,fig,ax)
            plot_cbar(sc,plt.cm.ScalarMappable(colors.Normalize(-np.pi, np.pi),cmap = plt.cm.hsv),args,title=r'\scalebox{1}{$\mathrm{Arg}(\hat{\rho})$}',cax=cax)
            if args.plottype=='polar':
                fig.tight_layout()

        plot_figure(fig,os.path.join(dir,basename),close=True)
        plt.close(fig)

def plot_fft_heatmap(args,fft_heatmap):
    lim=np.floor(fft_heatmap.shape[0]/2)
    maxval=np.max(fft_heatmap)
 
    fig,ax=plt.subplots()
    im=ax.imshow(fft_heatmap,'plasma',extent=(-lim,lim,-lim,lim))
    cax=get_cax(args,fig,ax)
    cbar=fig.colorbar(im,orientation='vertical',ax=ax, cax=cax)
    cbar.ax.set_title(r'\scalebox{0.6}{$\left|{\hat{\rho}}\right|$}')
    tick_positions=[round(i,2) for i in [0,maxval/2,maxval]]
    cbar.set_ticks(tick_positions)
    cbar.ax.tick_params(labelsize=12)
    plot_figure(fig,f"FFT_heatmap_{args.file.replace('.png','')}")
def main():
    args=parse_cmd_line()
    log.info(f'args={args}')

    diffractiondata=np.genfromtxt(f"diffractiondiagram_{args.file.replace('.png','')}.csv",delimiter=',')
    maxsize=max(diffractiondata[:,2])
    plot_diffraction_image([[diffractiondata,{'facecolors':'none','edgecolors':plt.cm.hsv(norm_phase(diffractiondata[:,3])),'linewidths':0.2}]],args,"diffraction_threshold",maxsize,title=r'FFT-Values $\hat{\rho}^{c}$ with amplitudes above c='+str(args.threshold),legendtitle=r'\hat{\rho}',do_legend=("diffractionfft" not in args.nolegend),do_yticks=("diffractionfft" not in args.noyticks),legendfacecolor='none')

    if os.path.isfile(f"fourier_coefficients_base_{args.file.replace('.png','')}.csv"):
        wavevecbase=np.genfromtxt(f"fourier_coefficients_base_{args.file.replace('.png','')}.csv",delimiter=',')[:,0:2]
        number_bases=wavevecbase.shape[0]
        wavevecother=np.genfromtxt(f"fourier_coefficients_other_{args.file.replace('.png','')}.csv",delimiter=',')[:,0:2]

        plot_peaks_and_fft(args,wavevecbase,wavevecother,diffractiondata)
        plot_fourier_at_idces(args,"wavevectors",maxsize)

        #threshholded diffractiondata
        oldplot=args.plottype
        args.plottype='cartesian'
        args.plottype=oldplot

        #gauges and associated errors
        gaugeerrorfiles=glob.glob('gaugeerrors*.csv')
        for ge_file in gaugeerrorfiles:
            log.info(f'analysing {ge_file}')
            plot_gaugeerrors_histogram(args,ge_file)
            plot_gaugeerrors(ge_file,args,maxsize,number_bases,diffractiondata)#pass diffractiondata here if you want to overlay gaugeerrors and fft


if __name__=="__main__":
    main()
