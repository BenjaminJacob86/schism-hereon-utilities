#!/usr/bin/env python
# Copyright: This document has been placed in the public domain.

""" Collection of formulas and plotting diagrams for validation """


__author__ = "Benjamin Jacob"
__copyright__   = "Copyright 2018 - 03\2021 Helmholtz-Zentrum Geesthacht"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"

__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"






############### Temporal interpolation ###############
def interp_to_data_time(tdata,vdata,tmod,vmod,min_maxdate=None):
	""" interpolates inputs to common time of tdata
	tinterp,vdata_interp,vmod_interp=interp_to_data_time(tdata,vdata,tmod,vmod,min_maxdate=None)  """
	# convert inputs to dateime 64 to avoid compatability issues
	from scipy.interpolate import interp1d	
	tdata=np.asarray(tdata,np.datetime64)
	tmod=np.asarray(tmod,np.datetime64)

	if min_maxdate !=None:
		date0,date1=np.asarray(min_maxdate,np.datetime64)
	else:
		date0=[]
		date1=[]
	date0=np.max([tdata[0],tmod[0],date0])
	date1=np.min([tdata[-1],tmod[-1],date1])

	idata=[(tdata >= date0) & (tdata<= date1)]
	tdata_intp=tdata[idata]
	vdata_intp=vdata[idata]
	
	tin=(tmod-tmod[0])/(tmod[1]-tmod[0])
	tout=(tdata_intp-tmod[0])/(tmod[1]-tmod[0])

	fintp=interp1d(tin,vmod)
	vmod_intp=fintp(tout)

	return tdata_intp,vdata_intp,vmod_intp



############### QQ PLot ##################################
def QQplot(dataObs,dataMod,stationname='',obsname='bouy',modname='SCHISM WWM',parameter='Hs',unit=' [m]'):	
	from scipy.stats import gaussian_kde
	import matplotlib as mpl
	xy = np.vstack([dataObs,dataMod])
	density = gaussian_kde(xy)(xy)#(xy)
	
	# stats
	entries=len(dataObs)
	meanR=dataObs.mean()
	meanM=dataMod.mean()
	SDR=np.std(dataObs)
	SDM=np.std(dataMod)
	RMSE=np.sqrt(np.mean((dataObs-dataMod)**2))
	E=dataMod-dataObs
	SE=np.sqrt(1/(entries-1)*np.sum((E-E.mean())**2))
	SI=SE/dataObs.mean()
	bias=np.mean(dataObs-dataMod)
	cor=np.corrcoef(dataObs,dataMod)[0,1]

	txt=' Entries = {:d} \n Mean R ={:.3} {:s} \n Mean M ={:.3} {:s} \n SD R ={:.3} {:s} \n SD M ={:.3} {:s} \n RMSE ={:.3} {:s} \n SI ={:.3} \n bias ={:.3} {:s} \n CORR ={:.3}'.format(entries,meanR,unit,meanM,unit,SDR,unit,SDM,unit,RMSE,unit,SI,bias,unit,cor)


	coef = np.polyfit(dataObs,dataMod,1)
	poly1d_fn = np.poly1d(coef) ## poly1d_fn is now a function which takes in x and returns an estimate for y 
	lbl=str('{:.2f}x + {:.2f}'.format(coef[0],coef[1]))

	# scatter
	plt.cla()
	plt.scatter(dataObs,dataMod,s=0.4,c=density)
	norm = mpl.colors.Normalize(vmin = np.min(density), vmax = np.max(density))
	#cbar = fig.colorbar(mpl.cm.ScalarMappable(norm = norm), ax=plt.gca())
	cbar=plt.colorbar(mpl.cm.ScalarMappable(norm = norm),pad=-0.1)
	cbar.ax.set_ylabel('Density')
	plt.axis('square')
	plt.xlabel(obsname+ ' ' +parameter + unit)
	plt.ylabel(modname+ ' ' +parameter + unit)
	#plt.title(stationname + ' ' + str(date0)[:10] + '- '+ str(date1)[:11])
	ylim=plt.ylim()
	xlim=plt.xlim()
	plt.text(0,ylim[1]/1.8,txt)
	plt.plot([0,xlim[1]],[0,ylim[1]],'k')
	ph=plt.plot([0,xlim[1]], poly1d_fn([0,ylim[1]]), 'r')	
	cbar.set_label('Density')
	plt.legend(ph,[lbl],loc='lower right',frameon=False)	



##########################################################


##############  Tayloer diagramm ########################
"""
Taylor diagram (Taylor, 2001) implementation.

Note: If you have found these software useful for your research, I would
appreciate an acknowledgment.
"""

__version__ = "Time-stamp: <2018-12-06 11:43:41 ycopin>"
__author__ = "Yannick Copin <yannick.copin@laposte.net>"

import numpy as np
import matplotlib.pyplot as plt


class TaylorDiagram(object):
    """
    Taylor diagram.

    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd,
                 fig=None, rect=111, label='_', srange=(0, 1.5), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.

        Parameters:

        * refstd: reference standard deviation to be compared to
        * fig: i'salt'ut Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        from matplotlib.projections import PolarAxes
        import mpl_toolkits.axisartist.floating_axes as FA
        import mpl_toolkits.axisartist.grid_finder as GF

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1])
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = GF.FixedLocator(tlocs)    # Positions
        tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = FA.GridHelperCurveLinear(
            tr,
            extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)

        if fig is None:
            fig = plt.figure()

        ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("Correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("Standard deviation")

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)          # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, 'k*',
                          ls='', ms=10, label=label)
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        l, = self.ax.plot(np.arccos(corrcoef), stddev,
                          *args, **kwargs)  # (theta, radius)
        self.samplePoints.append(l)

        return l

    def add_grid(self, *args, **kwargs):
        """Add a grid."""

        self._ax.grid(*args, **kwargs)

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(np.linspace(self.smin, self.smax),
                             np.linspace(0, self.tmax))
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rs**2 - 2*self.refstd*rs*np.cos(ts))

        contours = self.ax.contour(ts, rs, rms, levels, **kwargs)

        return contours


def plotTaylor(samples,stdref=1,extend=False):
    """
    generate talyor plot calling taylor diagram class.
	input stderef: referecnce standard deviation
	samples list of list each item containing standard deviation (normalized), corelation coefficient
	and a name to be show within legend
    """

    # Reference std
    stdref =1 

    # Samples std,rho,name
	# std corcoef name
	# surface data
	#samples=[ [compare[tag]['temp']['std2'][0]/compare[tag]['temp']['std1'][0],compare[tag]['temp']['cor'][0] ,tag+' '+str(Zsalt[tag][0])+'m'] for tag in names ]
    #samples = [[25.939, 0.385, "Model A"],
    #           [29.593, 0.509, "Model B"],
    #           [33.125, 0.585, "Model C"],
    #           [29.593, 0.509, "Model D"],
    #           [71.215, 0.473, "Model E"],
    #           [27.062, 0.360, "Model F"],
    #           [38.449, 0.342, "Model G"],
    #           [35.807, 0.609, "Model H"],
    #           [17.831, 0.360, "Model I"]]

    #fig = plt.figure()
    fig=plt.figure(figsize=[10,8]) # size in inch
    dia = TaylorDiagram(stdref, fig=fig, label='Reference', extend=extend)
    dia.samplePoints[0].set_color('r')  # Mark reference point as a red star

    # Add models to Taylor diagram
    for i, (stddev, corrcoef, name) in enumerate(samples):
        dia.add_sample(stddev, corrcoef,
                       marker='$%d$' % (i+1), ms=10, ls='',
                       mfc='k', mec='k',
                       label=name)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=7, colors='0.5')  # 5 levels in grey
    plt.clabel(contours, inline=1, fontsize=10, fmt='%.2f')

    dia.add_grid()                                  # Add grid
    dia._ax.axis[:].major_ticks.set_tick_out(True)  # Put ticks outward

    # Add a figure legend and title
    ncol=int(np.ceil(len(samples)/10))
    fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ],
               numpoints=1, prop=dict(size='small'), loc='upper right',ncol=ncol)
    fig.suptitle("Taylor diagram", size='x-large')  # Figure title
    plt.tight_layout()
    return dia






