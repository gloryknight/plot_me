import plot_me # works with version 2.22 (print plot_me.version)

import numpy
from numpy.random import randn
# random 2D
Z = numpy.clip(randn(250, 250), -1, 1)
i=plot_me.fig2d("x (pixel)","y(pixel)" ,"Intensity (arb. units)", zlimit=(0,10), colorformat="%d", extent=[0, 100, 0, 100])
i.add(Z).plot()

# 2D from arbitrary data (convert to grid first)
i=plot_me.fig2d("x (pixel)","y(pixel)" ,"Intensity (arb. units)", colorformat="%d")
i.addgrid([10,20,30,40],[10,20,25,10],[15,11,12,10], xpix=150, ypix=100).plot()
#i.show()

# Variour 1D plots with loading the data form txt and zip
i=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])
d=i.data().load_slow("plot_me_test.txt")
i.load("plot_me_test.txt").plot(0, 1, "load","", 1, log=0)
i.data().load_slow("plot_me_test.txt").plot(0, 1, "load_slow","", 1.1, log=0)
i.data().load_slow("plot_me_test.txt", 1).plot(0, 1, "load_slow_nr","", 1.11, log=0) # use -1 for the scan number
d=i.loadzip("plot_me_test.txt.zip","plot_me_test.txt")
d.plot(0, 1, "load_zip","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").fftsmoothme(10).plot(0, 1, "fftsmooth10","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").smoothme(10).plot(0, 1, "smooth10","", 1.2, log=0)
i.loadzip("plot_me_test.txt.zip","plot_me_test.txt").averageme(10).plot(0, 1, "average10","", 1.2, log=0)
#i.legend()
#i.save("1.png")

# Fitting 1D data
class line(plot_me.fit.Line): #line
	def peval(self, x, p): # evaluate the function (must be present)
		return (p[0]+(p[1]*x))

d=i.loadzip("plot_me_test.txt.zip","plot_me_test.txt")
try:
	print (d.fit(p0=[1, 1], range=[0,0], x_y_col=[0,1], function=line())) # print the parameters
	d.plotfit("fit") # plot the fit
except:
	pass

# Drawind warious data
d=i.data()
d.add1d([81, 83, 85]).plot(0,1, "add1d", log=0)
d.draw([5, 7, 2], [80, 84, 85], l="draw")

i.legend() # add the legend

#i.close()
#i=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])

#i.save(name1+".png")
#i.save(name1+".pdf")
#i.save(name1+".eps")
#i.save(name1+".svg")

i.show()