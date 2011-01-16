#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

version="2.35"

#Classes: fig->data->line, my_function

#import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy
import pylab
import sys
import zipfile

class config:
	'''
	configuration class for the script
	'''


#	figsize=(8,6) # size of the figure
#	adjust_bottom=0.11
	figsize=(8,5) # size of the figure
	figsize2d=(7,5) # size of the figure
	fonts=20
	dpi=200 #None
	axis_width=2 # axis line width
	markeredgewidth=2 # axis markers
	markersize=5 # axis marker length
	linewidth=2 # plot line with
	MarkerSize=linewidth*5
	adjust_left=0.11 # - parameters of the plot
	adjust_bottom=0.13
	adjust_right=0.96
	adjust_top=0.97
	x_label="X" # defauld x label
	y_label="Y" # defauld y label
#	colorformat='$10^{%d}$' # format of the color scale
	colorformat='%d' # format of the color scale
	cmap=pylab.cm.jet # default colormap of the 2d plot
	interpolation2d='nearest' # 2d interpolation
	origin2d='lower'
	whongdata="Wrong data in line: " # info message for wrong input data


#~ linestyle 	description
#~ '-' 	solid
#~ ''' 	dashed
#~ '-.' 	dash_dot
#~ ':' 	dotted

#~ marker 	description
#~ '.' 	point
#~ ',' 	pixel
#~ 'o' 	circle
#~ 'v' 	triangle_down
#~ '^' 	triangle_up
#~ '<' 	triangle_left
#~ '>' 	triangle_right
#~ '1' 	tri_down
#~ '2' 	tri_up
#~ '3' 	tri_left
#~ '4' 	tri_right
#~ 's' 	square
#~ 'p' 	pentagon
#~ '*' 	star
#~ 'h' 	hexagon1
#~ 'H' 	hexagon2
#~ '+' 	plus
#~ 'x' 	x
#~ 'D' 	diamond
#~ 'd' 	thin_diamond
#~ '|' 	vline
#~ '_' 	hline
#~ TICKLEFT 	tickleft
#~ TICKRIGHT 	tickright
#~ TICKUP 	tickup
#~ TICKDOWN 	tickdown
#~ CARETLEFT 	caretleft
#~ CARETRIGHT 	caretright
#~ CARETUP 	caretup
#~ CARETDOWN 	caretdown
#~ 'None' 	nothing
#~ ' ' 	nothing
#~ '' 	nothing
#~ ACCEPTS: [ '+' | '*' | ',' | '.' | '1' | '2' | '3' | '4'| '<' | '>' | 'D' | 'H' | '^' | '_' | 'd' | 'h' | 'o' | 'p' | 's' | 'v' | 'x' | '|' | TICKUP | TICKDOWN | TICKLEFT | TICKRIGHT | 'None' | ' ' | '' ]

def say(what):
	#print what
	pass

class zip: # load zip (can be used directly or from figure by loadzip)
	'''
	work with zip archives
	'''
	def __init__(self, zipname, mode="r"):
		self.zipname=zipname
		self.ffile = zipfile.ZipFile( zipname, mode, zipfile.ZIP_DEFLATED)
		self.fname=self.zipname+".txt"

	def read( self, filename ):
		''' read zip file with filename and return it's content '''
		self.filename=filename
		d=self.ffile.read(filename)
		return d

	def writestr(self, str, fname=None):
		''' write string to a file fname in the zip archive and close the file.'''
		if fname==None:
			fname=self.fname
		else:
			self.fname=fname
		self.ffile.writestr(fname, str)

	def close(self):
		''' close the file. '''
		self.ffile.close()

	def __del__(self):
		''' destructor '''
		self.close()

class fit: 
	''' class of functions used for fitting. Needed as a prototype. '''
	
	class Line: #line
		''' Line function used for fitting. Needed as a prototype. '''
		
		def peval(self, x, p): # evaluate the function (must be present)
			''' evaluate the function (must be present and should be modified according to needs)'''
			return (p[0]+(p[1]*x))

		def residuals(self, p, x, y):
			''' find deviation (must be present) '''
			err = y-self.peval(x,p)
			return err

	class Pseudo_Voigt(Line): 
		''' Pseudo-Voigt function used for fitting. '''

		def peval(self, x, p): 
			return (p[0]+p[1]*(abs(p[2])*self.lagr(x,p)+(1-p[2])*self.gaussian(x,p)))

		def lagr(self, x, p):
			return (1/(1+(((x-p[3])/p[4])**2)))
			
		def gaussian(self, x, p):
			return (numpy.exp(-numpy.log(2)*((x-p[3])/p[4])**2))

	class Gaussian(Line): 
		''' Gaussian function used for fitting. '''

		def peval(self, x, p): 
			return p[0]+p[1]*numpy.exp(-0.5*((x-p[2])/p[3])**2)

class data: 
	''' class holding the data (is used if you need to plot multiple columns of the sama data file) '''

	def __init__(self, ax, lw):
		self.ax=ax
		self.lw=lw
		self.data=None
		self.z=None

	def load(self, filename): 
		''' Loads data from the file into a new array, plots columns x_col versus y_col, with label scaled by scale '''
		self.data=self.read(filename)
		return self
		#return self.draw(self.data[:,x_col], self.data[:,y_col]*scale, self.lw, label)

	def summ(self, a, b):
		return numpy.vstack((a,b)) 

	def load_slow(self, filename, scan_nr=None): 
		''' Loads data from the file into a new array, plots columns x_col versus y_col, with label scaled by scale '''
		if (scan_nr==None) or (self.data==None):
			self.data=self.read_slow(filename, scan_nr)
		else:
			self.data=numpy.vstack((self.data, self.read_slow(filename, scan_nr)))
		return self

	def loadzip(self, zipname, filename, scan_nr=None): 
		''' Loads data from the file in a zip file returns new data object '''
		if ((self.z==None) or (self.z.zipname!=zipname)):
			del self.z
			self.z=zip(zipname)
		dd=self.z.read(filename)
		data=self.work_read(dd.split("\n"), scan_nr)
		if (scan_nr==None) or (self.data==None):
			self.add(data)
		else:
			self.data=numpy.vstack((self.data, data))
		return self

	def add(self, data): 
		''' add data to the class '''
		self.data=numpy.array(data, dtype='float')
		return self

	def add1d(self, data, xrange=[0,0]): 
		''' add data to the class '''
		data=numpy.array(data, dtype='float')
		rng=numpy.arange(0,data.shape[0])
		if xrange!=[0,0]:
			rng=xrange[0]+(rng*(xrange[1]-xrange[0])/float(data.shape[0]-1))
		self.data=numpy.column_stack((rng, data))
		return self

	def fftsmoothme(self, lowpass=10): 
		''' FFT frequency cut makes the data smooth (modyfies the data) '''
		gauss=fit.Gaussian().peval
		fft=numpy.fft.fft(self.data[:,1])
		fftl=len(fft)
		for i in range(1, fftl):
			fft[i]=fft[i]*(gauss(i, [0, 1, 0, lowpass])+gauss(i, [0, 1, fftl, lowpass]))
		gauss=None
		self.data[:,1]=numpy.fft.ifft(fft).real
		return self

	def smoothme(self, every): 
		''' running average of the data in this class (modyfies the data) '''
		size=numpy.shape(self.data)[0]
		for n in range(1, every):
			s=0
			for i in range(0,n):
				s+=self.data[i]
			self.data[n]=s/n
		for n in range(0, size-every):
			s=0
			for i in range(0,every):
				s+=self.data[n+i]
			self.data[n]=s/every
		for n in range(size-every, size):
			s=0
			for i in range(n,size):
				s+=self.data[i]
			self.data[n]=s/(size-n)
		return self

	def averageme(self, every): 
		''' average simplify data on the data of the class (reduces ammount of points) '''
		self.data=self.average(self.data, every)
		return self

	def average(self, data, every): 
		''' average simplify data (reduces ammount of points) '''
		line=0
		f1=[]
		s=[]
		for n in data:
			for m in range(n.size):
				try:
					s[m]+=n[m]
				except:
					s.append(0)
					s[m]+=n[m]
			if line==every:
				ff=[]
				for m in range(len(s)):
					ff.append(s[m]/(every+1))
				f1.append(ff)
				s=[]
				line=0
			else:
				line+=1
		return numpy.array(f1, dtype='float')

	def fit(self, range=[0,0], x_y_col=[0,1], p0=None, maximumfittingcycles=20000, function=None): 
		try:
			from scipy.optimize import leastsq
		except:
			print("no scipy")
		''' fit the data with function (returns the set of parameters) '''
		if range==[0,0]:
			data=self.data
		else:
			data=self.data[(self.data[:,x_y_col[0]]>range[0]) & (self.data[:,x_y_col[0]]<range[1]),:]
		if range==[0,0]:
			range=[data[0,x_y_col[0]],data[data[:,x_y_col[0]].size-1,x_y_col[0]]]
		if p0==None:
			self.p0_orig=[1.18643310e+02, 3.96555414e+02, 4.77081488e-06, 1.96415331e+01, 8.80491880e-02]
			#print(numpy.argmax(data[:,y_col]), numpy.size(data))
			n=numpy.argmax(data[:,x_y_col[1]])
			self.p0_orig[0]=data[0,x_y_col[1]] # background			
			self.p0_orig[1]=data[n,x_y_col[1]] # intensity
			self.p0_orig[2]=0.5 # ratio			
			self.p0_orig[3]=data[n,x_y_col[0]] # position
			self.p0_orig[4]=(range[1]-range[0])/3. # FWHM			
		else:
			self.p0_orig=p0
		if function==None:
			self.fitfunc=Pseudo_Voigt()
		else:
			self.fitfunc=function
		self.p0=leastsq(self.fitfunc.residuals, self.p0_orig, args=(data[:,x_y_col[0]], data[:,x_y_col[1]]), maxfev=maximumfittingcycles)[0]
		return self.p0

	def plotfit(self, name="fit", x_y_col=[0, 1], **kwargs):
		''' plot the result of fitting '''
		data1=self.data.copy()
		data1[:,x_y_col[1]]=self.fitfunc.peval(data1[:,x_y_col[0]],self.p0)
		#y=self.fitfunc.peval(self.data[:,x],self.p0)
		self.plot(x_y_col[0], x_y_col[1], name,"", data=data1, **kwargs)

	def plot(self, x_col=0, y_col=1, label=None, marker='', scale=1, data=None, log=0, **kwargs):
		''' plot the data '''
		if data==None:
			data=self.data
		if scale==0:
			scale=1/self.data[:,y_col].max()
		if log==1:
			plotf=self.ax.semilogy
		elif log==2:
			plotf=self.ax.semilogx
		elif log==3:
			plotf=self.ax.loglog
		else:
			plotf=self.ax.plot
		return self.draw(data[:,x_col], data[:,y_col]*scale, marker, self.lw, label, scale, plotf=plotf, **kwargs)

	def draw(self, x, y, marker='', lw=config.linewidth, l=None, scale=1, plotf=None, **kwargs):
		''' draw the data (pass parameters to the plot directrly)'''
		if plotf==None:
			plotf=self.ax.plot
		l, = plotf(x, y, marker, lw=lw, label=l, markersize=config.MarkerSize, **kwargs)
		return l

	def read(self, filename):
		''' Loads data from the file into a new array '''
		return numpy.loadtxt(filename)

	def read_slow(self, filename, scan_nr=None):
		''' Loads data from the file into a new array - slow but error prune '''
		self.filename=filename
		fr=open(filename,'r')
		data=self.work_read(fr, scan_nr)
		fr.close()
		return numpy.array(data, dtype='float')

	def work_read(self, fr, scan_nr=None):
		''' helper function for read_slow'''
		data=[]
		n=1
		first=True
		for line in fr:
			if line != "" and line != "\n" and line[0]!="#":
				line=line.replace(",",".")
				dat=[]
				good=True
				l=line.split()
				if len(l)<1:
					good=False
				for d in l:
					try:
						fd=float(d)
						dat.append(float(d))
					except:
						good=False
				if good:
					if first:
						datlen=len(dat)
						first=False
					else:
						if datlen!=len(dat):
							good=False
				if good:
					if scan_nr!=None:
						dat.append(float(scan_nr))
					#print dat
					data.append(dat)
				else:
					say(config.whongdata+str(n))
			n+=1
		return data

class fig:
	''' 2D figure and operations on it. '''

	def __init__(self, xt=config.x_label, yt=config.y_label, xlimit=None, ylimit=None, lw=config.linewidth, fonts=config.fonts): 
		''' Initialize new canvas '''
		self.xt=xt # x title
		self.yt=yt # y title
		self.lx=xlimit # x range eg. [0,1]
		self.ly=ylimit # y range eg. [0,1]
		self.lw=lw # line width
		self.fonts=fonts # font size
		self.plotsetup()

	def plotsetup(self):
		''' finalize the plot setup '''
		self.fig = plt.figure(figsize=config.figsize)
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlabel(self.xt, fontsize=self.fonts)
		self.ax.set_ylabel(self.yt, fontsize=self.fonts)
		self.fig.subplots_adjust(left=config.adjust_left, bottom=config.adjust_bottom, right=config.adjust_right, top=config.adjust_top, wspace=None, hspace=None)

		for spine in self.ax.spines.itervalues():
		      spine.set_linewidth(config.axis_width)
		for label in self.ax.get_xticklabels() + self.ax.get_yticklabels():
			label.set_fontsize(self.fonts) 
		for line in self.ax.xaxis.get_ticklines() + self.ax.yaxis.get_ticklines():
#			line.set_color('green')
			line.set_markersize(config.markersize)
			line.set_markeredgewidth(config.markeredgewidth)

	def data(self):
		''' returns a data object '''
		return data(self.ax, self.lw)

	def load(self, filename):
		''' Loads data from the file into a new data object '''
		d=data(self.ax, self.lw)
		return d.load(filename)

	def loadzip(self, zipname, filename, scan_nr=None):
		''' Loads data from the file in a zip file returns new data object '''
		do=data(self.ax, self.lw)
		return do.loadzip(zipname, filename, scan_nr)

	def legend(self, *args, **kwargs):
		''' Adds a legend. Use loc. 1 to 10 to change location '''
		leg=plt.legend(*args, **kwargs)
		for t in leg.get_texts():
			t.set_fontsize(self.fonts) # the legend text fontsize

	def axis(self):
		''' set axis limits '''
		if self.lx!=None:
			self.ax.set_xlim(self.lx[0],self.lx[1])
		if self.ly!=None:
			self.ax.set_ylim(self.ly[0],self.ly[1])

	def show(self):
		''' show the plot '''
		self.axis()
		plt.show()

	def save(self, filename):
		''' save the plot to a file '''
		self.axis()
		plt.savefig(filename, dpi = (config.dpi))

	def label(self, x, y, text, dir=0, **kwargs):
		''' add a label to the plot '''
		self.ax.text(x,y, text, rotation=dir, size=self.fonts, **kwargs) #, color='red'

	def plot(self, *args, **kwargs):
		''' plot the data '''
		data(self.ax, self.lw).draw(*args, **kwargs)

	def close(self):
		plt.close()

class fig2d:
	''' 2D+color figure and operations on it. '''
	def __init__(self, xt="X", yt="Y", cbt="Z", xlimit=None, ylimit=None, zlimit=None, linewidth=2, fonts=config.fonts, colorformat=config.colorformat, aspect=1, extent=None): 
		''' Initialize new canvas '''
		self.xt=xt # x title
		self.yt=yt # y title
		self.cbt=cbt # y title
		self.lx=xlimit # x range eg. [0,1]
		self.ly=ylimit # y range eg. [0,1]
		self.lz=zlimit # y range eg. [0,1]
		self.lw=linewidth # line width
		self.fonts=fonts # font size
		self.colorformat=colorformat
		self.aspect=aspect
		self.extent=extent
	
	def plotsetup(self):
		''' finalize the plot setup '''
		self.ax.set_xlabel(self.xt, fontsize=self.fonts)
		self.ax.set_ylabel(self.yt, fontsize=self.fonts)
		self.cb = plt.colorbar(format=pylab.FormatStrFormatter(self.colorformat)) # draw colorbar
		self.cb.set_label(self.cbt, rotation=-90, fontsize=self.fonts)
		self.fig.subplots_adjust(left=config.adjust_left, bottom=config.adjust_bottom, right=config.adjust_right-0.2, top=config.adjust_top, wspace=None, hspace=None)

		for label in self.ax.get_xticklabels() + self.ax.get_yticklabels():
			label.set_fontsize(self.fonts) 
		for line in self.ax.xaxis.get_ticklines() + self.ax.yaxis.get_ticklines():
			line.set_markeredgewidth(3)
		for t in self.cb.ax.get_yticklabels():
			t.set_fontsize(self.fonts)

	def add(self, data):
		''' add data to the plot (must be a three dimensional numpy array) '''
		self.data=data
		return self

	def addgrid(self, dqx, dqy, di, xpix=1500, ypix=1000):
		dqx=numpy.array(dqx)
		dqy=numpy.array(dqy)
		di=numpy.array(di)
		xi = numpy.linspace(dqx.min(),dqx.max(),xpix)
		yi = numpy.linspace(dqy.min(),dqy.max(),ypix)
		self.data=griddata(dqx,dqy,di,xi,yi)
		return self

	def plotdata(self, data):
		''' plot the data '''
		self.add(data)
		self.plot()

	def plot(self):
		''' plot the figure '''
		self.fig = plt.figure(figsize=config.figsize2d)
		self.ax = self.fig.add_subplot(111)
		self.ax.hold()
#		self.ax = plt.axes()

		if self.lz!=None:
			l = plt.imshow(self.data,vmin=self.lz[0],vmax=self.lz[1],cmap=config.cmap, interpolation=config.interpolation2d,origin=config.origin2d, aspect=self.aspect, extent=self.extent)
		else:
			l = plt.imshow(self.data,cmap=config.cmap, interpolation=config.interpolation2d,origin=config.origin2d, aspect=self.aspect, extent=self.extent)
		self.plotsetup()
#		self.ax.set_xticks(xtricks)
		return self

	def show(self):
		''' show the figure '''
		plt.show()

	def save(self, *args, **kwargs):
		''' save the figure '''
		plt.savefig( *args, **kwargs)

	def label(self, x, y, text, dir=0):
		''' add a label '''
		self.ax.text(x,y, text, rotation=dir, size=self.fonts) #, color='red'

	def close(self):
		plt.close()

if __name__ == "__main__":
	
	# do not use psyco - (returns an error!)

	if len(sys.argv)==2:
		i=fig()
		i.load(sys.argv[1]).plot(0, 1)
		i.show()
	elif len(sys.argv)==4:
		i=fig()
		i.load(sys.argv[1]).plot(sys.argv[2], sys.argv[3])
		i.show()
	elif len(sys.argv)==6:
		i=fig(sys.argv[4], sys.argv[5])
		i.load(sys.argv[1]).plot(sys.argv[2], sys.argv[3])
		i.show()
	else:
		print("Usage: plot_me filename x_collumn y_collumn x_label y_label") # eg. plot_me test.txt 0 1 time intensity


#--- examples:

#	i=fig(r"$2\theta$ (degrees)", "Intensity [arb. units]", xlimit=[27,80], ylimit=[0.165,1])
#	i.load("9104682_grad_02_0001a.sfrm.Q-int.dat").plot(0, 1, "first","",0)

#	d=i.data()
#	dw=d.load("9104682_grad_02_0001a.sfrm.Q-int.dat")

#	line=d.add(d.smooth(dw,2)).plot(0, 1, "second","",0)

#	d.add(dw).plot(0, 1, "second","",line.scale)

#	d=i.loadzip("CIGS-gradient1_00003.TIF_th10-62.zip", name).plot(0, 1, "a","",0)
#	i.save(name+".eps")
#	i.show()

	#~ import numpy
	#~ from plot_me import fig2d
	#~ from numpy.random import randn
	#~ Z = numpy.clip(randn(250, 250), -1, 1)
	#~ i=fig2d("x [pixel]","y[pixel]" ,"Intensity [arb. units]", zlimit=(0,10))
	#~ i.add(Z).plot().show()
