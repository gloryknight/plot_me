import plot_me # works with version 2.36 (print plot_me.version)
#import pylab
#import matplotlib as mpl

import numpy as np
from numpy.random import randn

# random 2D
Z = np.clip(randn(250, 250), -1, 1)
i=plot_me.fig2d("x (pixel)","y(pixel)" ,"Intensity (arb. units)", zlimit=(0,10), colorformat="%f", extent=[0, 100, 0, 100], fixcbsize=6, cbpad=0.08)
#i.ncmap='hot'
i.ncolors=10
#i.colorformat='$10^{%d}$'
i.add(Z).plot()
i.ax.text(0.5,1, "Voltage=1V", fontsize=20,bbox=dict(facecolor='white', alpha=0.5))
i.add_line(np.array([[10,40], [90, 90], [80, 50]]), lw=1, color="red")

i.fig.subplots_adjust(left=0.14, right=0.69, bottom=0.13, top=0.94)

plot_me.config.bbox_inches=None
#i.save("test1.png")
plot_me.config.bbox_inches='tight'

plot_me.config.ncmap='Paired'
# 2D from arbitrary data (convert to grid first)
i=plot_me.fig2d("x (pixel)","y(pixel)" ,"Intensity (arb. units)", colorformat="%g")
i.hcorient=1
i.addgrid(np.array([10,20,30,40]),np.array([10,20,25,10]),np.array([15,11,12,10]), xpix=150, ypix=100).plot()
#i.show()

plot_me.config.ncmap='Dark2'
plot_me.config.ncolors=11 # this trick also works for 2D

# Variour 1D plots with loading the data form txt and zip
i=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])
d=i.data().load_slow("plot_me_test.txt")
i.load("plot_me_test.txt").plot(0, 1, "load","", 1, log=0)
#i.ly=[90,140]
i.data().load_slow("plot_me_test.txt").plot(0, 1, "load_slow","", 1.1, log=0)
i.data().load_slow("plot_me_test.txt", 1).plot(0, 1, "load_slow_nr","", 1.11, log=0) # use -1 for the scan number
d=i.loadzip("plot_me_test.txt.zip","plot_me_test.txt")
d.plot(0, 1, "load_zip","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").fftsmoothme(10).plot(0, 1, "fftsmooth10","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").smoothme(10).plot(0, 1, "smooth10","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").averageme(10).plot(0, 1, "average10","", 1.2, log=0)
#i.save("1.png")

#Simple multidimentional plot of data with independant fit
td=np.array([[80,91], [87,90], [89,93], [84,91], [89,80]]) # create y data in two columns
nd=np.column_stack(([[5, 8, 9, 6, 5],td])) # add x column in front
ttd=plot_me.data(1,2)

ttd.data=nd
fpar=ttd.fit(p0=[1, 1], range=[0,0], x_y_col=[0,1], function=plot_me.fit.Line())
print "The fit par:", fpar
nd=np.column_stack(([nd,plot_me.fit.Line().peval(nd[:,0],fpar[0])])) # add x column in front

ni=plot_me.fig("y (arb. units)", "y (arb. units)")
ni.data().add(nd).plot(0,[1,2,3],"multy_dim") # y is a tupple of the colluns to be ploted
ni.legend()

# Fitting 1D data
class line(plot_me.fit.Line): #line
	def peval(self, x, p): # evaluate the function (must be present)
		return (p[0]+(p[1]*x))

d=i.loadzip("plot_me_test.txt.zip","plot_me_test.txt")
try:
	print "parameters leastsq:", (d.fit(p0=[1, 1], range=[0,0], x_y_col=[0,1], function=line())) # print the parameters
	d.plotfit("fit") # plot the fit
	err1=d.getfiterr()
	print "error:", err1
	print "parameters fmin_slsqp:", (d.fit(p0=[91, 1.2], range=[0,0], x_y_col=[0,1], function=plot_me.fit.Line(), leastsq=False, maxerr=1e-70, debug=0, boundaries=[(90, 92),(0.5, 1.5)])) # print the parameters
	d.plotfit("fit") # plot the fit
	err2=d.getfiterr()
	print "error:", err2, "\n"
	print "error:", (err1-err2)
except:
	pass

i.legend()

# Drawind warious data
d=i.data()
d.add1d([81, 83, 87], xrange=[10, 15]).plot(0,1, "add1d", log=0)
d.draw([5, 8, 9, 6, 5], [80, 81, 85, 84, 80], l="draw")

i4=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])
d4=i4.data()
plot_me.mpl.rcParams['lines.linewidth']=4 # configuration see the "site-packages\matplotlib\mpl-data\matplotlibrc" file
i4.ax.plot([5, 8, 9, 6, 5], td, label="multy_draw")
d4.draw([5, 8, 9, 6, 5], [81, 85, 85, 80, 83], l="draw1")
#plot_me.config.ltrans=1 #set transparancy of the legend to zero.
i.legend() # add the legend


#i.close()
#i=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])

#i.save(name1+".png")
#i.save(name1+".pdf")
#i.save(name1+".eps")
#i.save(name1+".svg")

i.show()