import os
import numpy as np
from matplotlib import pyplot as plt
from tools.antenna import antenna_info
from matplotlib.colors import LogNorm

Antenna = antenna_info()[0]

color16=['firebrick','saddlebrown', 'deeppink','tomato', 'lightsalmon','greenyellow','forestgreen'
         ,'olive','cyan','steelblue','dodgerblue','navy' ,'purple','lightslategray','gray','black']

def plot_16_log_theta(xlabel,ylabel,title
                ,x_data,y_data
                ,x_data1,y_data1
                ,x_data2,y_data2
                ,x_data3,y_data3
                ,ymin,ymax
                ,d_path,file_name
                ,message):
    
    fig = plt.figure(figsize=(24, 18)) # figure size

    ax = fig.add_subplot(111) # want to have a one big XY label for all 16 plots
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel(xlabel, labelpad=30,fontsize=40)
    ax.set_ylabel(ylabel,labelpad=50, fontsize=40)

    plt.title(title, y=1.04,fontsize=35)

    for b in range(16):
   
        ax = fig.add_subplot(4,4,b+1)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=20)
        #ax.set_xlim(-200,800)
        #ax.set_ylim(-0.8,0.8)
        ax.set_xlim(0,1)
        #ax.set_ylim(1e-4,1e2)
        ax.set_ylim(ymin,ymax)
        ax.set_yscale('log')
        ax.grid(linestyle=':')
        ax.set_title(Antenna[b],fontsize=25)
        
        if b == 15:
            ax.plot(x_data,y_data[:,b],'-',lw=2,color='red',alpha=0.7,label='fft')
            ax.plot(x_data1,y_data1[:,b],'-',lw=2,color='navy',alpha=0.7,label='fft w/ 0 deg band')
            ax.plot(x_data2,y_data2[:,b],'-',lw=2,color='blue',alpha=0.7,label='fft w/ -36 deg band')
            ax.plot(x_data3,y_data3[:,b],'-',lw=2,color='cyan',alpha=0.7,label='fft w/ -72 deg band')

            plt.legend(loc='best',numpoints = 1 ,fontsize=15)

        else:
            ax.plot(x_data,y_data[:,b],'-',lw=2,color='red',alpha=0.7)
            ax.plot(x_data1,y_data1[:,b],'-',lw=2,color='navy',alpha=0.7)
            ax.plot(x_data2,y_data2[:,b],'-',lw=2,color='blue',alpha=0.7)
            ax.plot(x_data3,y_data3[:,b],'-',lw=2,color='cyan',alpha=0.7)
        
    plt.tight_layout()

    # saving png into output path
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    print(message)

def plot_16(xlabel,ylabel,title
            ,x_data,y_data,legend
            ,xmin,xmax
            ,d_path,file_name
            ,message):
    
    fig = plt.figure(figsize=(24, 18)) # figure size

    ax = fig.add_subplot(111) # want to have a one big XY label for all 16 plots
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel(xlabel, labelpad=30,fontsize=40)
    ax.set_ylabel(ylabel,labelpad=50, fontsize=40)

    plt.title(title, y=1.04,fontsize=35)

    for b in range(16):
   
        ax = fig.add_subplot(4,4,b+1)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=20)
        #ax.set_xlim(-200,800)
        ax.set_xlim(xmin,xmax)
        #ax.set_ylim(-0.8,0.8)
        ax.grid(linestyle=':')
        ax.set_title(Antenna[b],fontsize=25)

        ax.plot(x_data,y_data[:,b],'-',lw=2,color='red',alpha=0.7,label='Vpeak:'+str(legend[b]))

        plt.legend(loc='best',numpoints = 1 ,fontsize=15)

    plt.tight_layout()

    # saving png into output path
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    print(message)

def plot_16_indi(xlabel,ylabel,title
            ,x_data,y_data
            ,d_path,file_name
            ,message
            ,xmin = None, xmax = None
            ,ymin = None, ymax = None
            ,y_scale = None
            ,vpeak = None
            ,absol = None):

    fig = plt.figure(figsize=(24, 18)) # figure size

    ax = fig.add_subplot(111) # want to have a one big XY label for all 16 plots
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel(xlabel, labelpad=30,fontsize=40)
    ax.set_ylabel(ylabel,labelpad=50, fontsize=40)

    plt.title(title, y=1.04,fontsize=35)

    #for b in range(16):
    for b in range(20):

        #ax = fig.add_subplot(4,4,b+1)
        ax = fig.add_subplot(5,4,b+1)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=20)
        #ax.set_xlim(-200,800)
        if xmin is not None and xmax is not None:
            ax.set_xlim(xmin,xmax)
        if ymin is not None and ymax is not None:
            ax.set_ylim(ymin,ymax)
        if y_scale is not None:
            ax.set_yscale('log')    
        #ax.set_ylim(-0.8,0.8)
        ax.grid(linestyle=':')
        #ax.set_title(Antenna[b],fontsize=25)
        ax.set_title(f'ch{b}',fontsize=25)

        if absol is not None:
            y_dat = np.abs(y_data[b])
        else:
            y_dat = y_data[b]

        if vpeak is not None:
            ax.plot(x_data[b],y_dat,'-',lw=2,color='red',alpha=0.7,label='Vpeak:'+str(np.round(np.nanmax(np.abs(y_data[b])),2)))
            plt.legend(loc='best',numpoints = 1 ,fontsize=15)
        else: 
            ax.plot(x_data[b],y_dat,'-',lw=2,color='red',alpha=0.7)

    plt.tight_layout()

    # saving png into output path
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    print(message)

def plot_16_3(xlabel,ylabel,title
            ,x_data,y_data,legend
            ,x_data1,y_data1,legend1
            ,x_data2,y_data2,legend2
            ,d_path,file_name
            ,message):
    
    fig = plt.figure(figsize=(24, 18)) # figure size

    ax = fig.add_subplot(111) # want to have a one big XY label for all 16 plots
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel(xlabel, labelpad=30,fontsize=40)
    ax.set_ylabel(ylabel,labelpad=50, fontsize=40)

    plt.title(title, y=1.04,fontsize=35)

    for b in range(16):
   
        ax = fig.add_subplot(4,4,b+1)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=20)
        #ax.set_xlim(-200,800)
        #ax.set_ylim(-0.8,0.8)
        ax.grid(linestyle=':')
        ax.set_title(Antenna[b],fontsize=25)

        ax.plot(x_data,y_data[:,b],'-',lw=2,color='red',alpha=0.7,label=legend[b])
        ax.plot(x_data1,y_data1[:,b],'-',lw=2,color='blue',alpha=0.7,label=legend1[b])
        ax.plot(x_data2,y_data2[:,b],'-',lw=2,color='green',alpha=0.7,label=legend2[b])

        plt.legend(loc='best',numpoints = 1 ,fontsize=15)

    plt.tight_layout()

    # saving png into output path
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    print(message)

def plot_16_overlap(xlabel,ylabel,title
            ,x_data,y_data
            ,d_path,file_name
            ,pol
            ,message):
 
    if pol == 'Vpol':
        Ant = Antenna[:8]
    if pol == 'Hpol':
        Ant = Antenna[8:]

    fig = plt.figure(figsize=(10, 7))
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    plt.grid(linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=15)
    #plt.yscale('log')
    #plt.xlim(140,180)
    #plt.ylim(140,180)

    for b in range(8):

        plt.plot(x_data,y_data[:,b],'-',lw=3,color=color16[b],alpha=0.5,label=Ant[b])

    plt.legend(loc='lower center',bbox_to_anchor=(1.15,0), numpoints = 1 ,fontsize=15)
    #plt.legend(loc='best', numpoints = 1 ,fontsize=15)
    
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()
  
    del Ant 

    print(message)

def plot_1(xlabel,ylabel,title
                ,x_data,y_data
                ,d_path,file_name
                ,message):

    fig = plt.figure(figsize=(10, 7))
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    plt.grid(linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=15)
    #plt.yscale('log')
    #plt.xlim(140,180)
    #plt.ylim(0,25)

    plt.plot(x_data,y_data,'-',lw=3,color='red',alpha=0.5,label=str(np.round(np.nanmax(y_data),2)))

    #plt.legend(loc='lower center',bbox_to_anchor=(1.23,0), numpoints = 1 ,fontsize=15)
    plt.legend(loc='best', numpoints = 1 ,fontsize=15)
    
    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close() 

    print(message)

def sky_map(evt,nadir_range, phi_range, angle_width
                , title
                , d_path, file_name
                , message):

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.ylabel(r'Nadir Angle [ $deg.$ ]', fontsize=25)
    plt.xlabel(r'Azimuthal Angle [ $deg.$ ]', fontsize=25)
    plt.grid(linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=20)
    plt.ylim(180,0)
    plt.xlim(0,360)

    cc = plt.imshow(evt, interpolation='none',cmap='viridis',extent=[phi_range.min()-angle_width/2, phi_range.max()+angle_width/2, nadir_range.max()+angle_width/2, nadir_range.min()-angle_width/2], aspect='auto')

    evt_max = np.nanmax(evt)
    evt_loc = np.where(evt == evt_max)
    evt_max_round = np.round(evt_max,2)

    plt.plot(evt_loc[1][0]*angle_width+angle_width/2,evt_loc[0][0]*angle_width+angle_width/2,'*',markersize=7, color='red'
             , label='S:'+str(evt_max_round)+', N:'+str(evt_loc[0][0]*angle_width+angle_width/2)+'deg, P:'+str(evt_loc[1][0]*angle_width+angle_width/2)+'deg')

    cbar1 = plt.colorbar(cc, ax=ax)
    cbar1.ax.tick_params(axis='y', labelsize=15)
    cbar1.ax.set_ylabel(r'Averaged Event-wise SNR [ $V/RMS$ ]', fontsize=15)#, color= 'red')

    plt.legend(loc='lower center',bbox_to_anchor=(0.5,-0.35), numpoints = 1 ,fontsize=15)

    plt.tight_layout()

    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    del evt_max, evt_loc, evt_max_round

    print(message)

def hist_map(xlabel,ylabel,title
                ,trig,evt_snr_v,evt_snr_h
                ,d_path,file_name
                ,message):

    rf_loc = np.where(trig == 0)[0]
    cal_loc = np.where(trig == 1)[0]
    soft_loc = np.where(trig == 2)[0]

    evt_w_snr_v_rf = evt_snr_v[rf_loc]
    evt_w_snr_v_cal = evt_snr_v[cal_loc]
    evt_w_snr_v_soft = evt_snr_v[soft_loc]
    evt_w_snr_h_rf = evt_snr_h[rf_loc]
    evt_w_snr_h_cal = evt_snr_h[cal_loc]
    evt_w_snr_h_soft = evt_snr_h[soft_loc]

    #plot
    fig = plt.figure(figsize=(10, 7))
    plt.xlabel(r'Averaged & event-wise SNRs [ $V/RMS$ ]', fontsize=25)
    plt.ylabel(r'# of events', fontsize=25)
    plt.grid(linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=15)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(0,300)
    #plt.ylim(0.9,500)
    plt.xlim(1,1e2)
    plt.ylim(1e-1,1e4)

    bin_range_step=np.arange(-10,10000,0.05)

    plt.hist(evt_w_snr_v_soft,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='cyan',alpha=0.7,label='A2 Vpol soft (1run)')
    plt.hist(evt_w_snr_h_soft,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='orange',alpha=0.7,label='A2 Hpol soft (1run)')
    plt.hist(evt_w_snr_v_rf,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='navy',alpha=0.7,label='A2 Vpol rf (1run)')
    plt.hist(evt_w_snr_h_rf,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='red',alpha=0.7,label='A2 Hpol rf (1run)')
    plt.hist(evt_w_snr_v_cal,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='dodgerblue',alpha=0.7,label='A2 Vpol cal (1run)')
    plt.hist(evt_w_snr_h_cal,bins=bin_range_step,histtype='step',linewidth=3,linestyle='-',color='orangered',alpha=0.7,label='A2 Hpol cal (1run)')

    plt.legend(loc='best', numpoints = 1 ,fontsize=15)

    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()

    del rf_loc, cal_loc, evt_w_snr_v_rf, evt_w_snr_v_cal, evt_w_snr_h_rf, evt_w_snr_h_cal, bin_range_step

    print(message)

def hist_2d(title
            , xlabel, ylabel
            , d_path
            , file_name
            , cmap_c = 'viridis'
            , cbar_legend = None
            , x_range = None, y_range = None
            , x_dat = None, y_dat = None
            , x_dat1 = None, y_dat1 = None, dat1_legend = None
            , y_dat1_err = None, dat1_err_legend = None
            , x_dat2 = None, y_dat2 = None, dat2_legend = None
            , y_dat2_err = None, dat2_err_legend = None
            , xmin = None, xmax = None
            , ymin = None, ymax = None
            , message = None):

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.ylabel(ylabel, fontsize=25)
    plt.xlabel(xlabel, fontsize=25)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=15)

    if x_dat is not None:
        if x_range is not None:
            cc = plt.hist2d(x_dat, y_dat, bins=(x_range, y_range), cmap=cmap_c,norm=LogNorm())
        else:
            cc = plt.hist2d(x_dat, y_dat, cmap=cmap_c,norm=LogNorm())
    
        cbar1 = plt.colorbar(cc[3], ax=ax)
        cbar1.ax.tick_params(axis='y', labelsize=15)
        if cbar_legend is not None:
            cbar1.ax.set_ylabel(cbar_legend, fontsize=15)#, color= 'red')

    if x_dat1 is not None:
        plt.plot(x_dat1,y_dat1,'-',lw=3,color='red',alpha=0.8,label=dat1_legend)
        if y_dat1_err is not None:
            plt.fill_between(x_dat1, y_dat1-y_dat1_err, y_dat1+y_dat1_err, linestyle = '--', linewidth = 2, edgecolor = 'red', facecolor='none', alpha=0.8,label=dat1_err_legend)

    if x_dat2 is not None:
        plt.plot(x_dat2,y_dat2,'-',lw=3,color='magenta',alpha=0.8,label=dat2_legend)
        if y_dat2_err is not None:
            plt.fill_between(x_dat2, y_dat2-y_dat2_err, y_dat2+y_dat2_err, linestyle = '--', linewidth = 2, edgecolor = 'magenta', facecolor='none', alpha=0.8,label=dat2_err_legend)

    if dat1_legend is not None or dat2_legend is not None or dat1_err_legend is not None or dat2_err_legend is not None:
        plt.legend(loc='best',numpoints = 1 ,fontsize=10)

    plt.grid(linestyle=':')
    if ymin is not None:
        plt.ylim(ymin,ymax)
    if xmin is not None:
        plt.xlim(xmin,xmax)
    plt.tight_layout()

    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()
    del fig#, ax, cc, cbar1
 
    if message is not None:
        print(message)

def imshow_2d(title
            , xlabel, ylabel
            , d_path
            , file_name
            , cmap_c = 'viridis'
            , cbar_legend = None
            , x_range = None, y_range = None
            , x_width = None, y_width = None
            , xy_dat = None
            , cmin = None, cmax = None
            , xmin = None, xmax = None
            , ymin = None, ymax = None
            , message = None):

    fig, ax = plt.subplots(figsize=(12, 6))
    plt.ylabel(ylabel, fontsize=25)
    plt.xlabel(xlabel, fontsize=25)
    plt.grid(linestyle=':')
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    plt.title(title, y=1.02,fontsize=15)

    if xy_dat is not None:
        if x_range is not None:
            cc = plt.imshow(xy_dat, interpolation='none',cmap= cmap_c
                    #,norm=LogNorm()
                    ,extent=[x_range.min()-x_width/2, x_range.max()+x_width/2, y_range.min()-y_width/2, y_range.max()+y_width/2]
                    , aspect='auto')
            if cmin is not None:
                cc = plt.imshow(xy_dat, interpolation='none',cmap= cmap_c
                        #,norm=LogNorm()
                        ,vmin = cmin, vmax = cmax
                        ,extent=[x_range.min()-x_width/2, x_range.max()+x_width/2, y_range.min()-y_width/2, y_range.max()+y_width/2]
                        , aspect='auto')
        else:
            cc = plt.imshow(xy_dat, interpolation='none',cmap= cmap_c
                    #,norm=LogNorm()
                    #,extent=[x_range.min()-x_width/2, x_range.max()+x_width/2, y_range.min()-y_width/2, y_range.max()+y_width/2]
                    , aspect='auto')
            if cmin is not None:
                cc = plt.imshow(xy_dat, interpolation='none',cmap= cmap_c
                        #,norm=LogNorm()
                        ,vmin = cmin, vmax = cmax
                        #,extent=[x_range.min()-x_width/2, x_range.max()+x_width/2, y_range.min()-y_width/2, y_range.max()+y_width/2]
                        , aspect='auto')

        cbar1 = plt.colorbar(cc, ax=ax)
        cbar1.ax.tick_params(axis='y', labelsize=15)
        if cbar_legend is not None:
            cbar1.ax.set_ylabel(cbar_legend, fontsize=15)#, color= 'red')

    plt.grid(linestyle=':')
    if ymin is not None:
        plt.ylim(ymin,ymax)
    if xmin is not None:
        plt.xlim(xmin,xmax)
    plt.tight_layout()

    os.chdir(d_path)
    fig.savefig(file_name,bbox_inches='tight')
    #plt.show()
    plt.close()
    del fig#, ax, cc, cbar1

    if message is not None:
        print(message)















