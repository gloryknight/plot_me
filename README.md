plot_me: Production quality plotting library.
=======

## What is it
It is a frontend for producing plots suitable for publication in scientific journals.

## Main Features
- 1D plots (fig)
- 2D plots (fig2d)
- Add data from arrays (add, add1d)
- Read data from txt and inside of zip files (read, loadzip)
- Running averaging, FFT frequency cut and windowing (smoothme, fftsmoothme, smooth)
- 1D fitting (fit, plotfit, getfiterr)
- Polynomial 2D fitting (polyfit2d)


## Dependencies
- [NumPy](http://www.numpy.org): 1.6.1 or higher
- [matplotlib](http://matplotlib.sourceforge.net/): the backend
- [SciPy](http://www.scipy.org): fitting
- standard libraries

## Installation from sources
Just copy "plot_me.py" to you **site-packages** library directory.

## Usage
```python
import plot_me
i=plot_me.fig("Time (seconds)", "Power (arb. units)") #, ylimit=[30,180])
i.load("plot_me_test.txt").plot(0, 1, "series1","", 1, log=0)
i.save("Power.png")
```
See "plot_me_test.txt.py" for further examples.

## License
BSD

## Background
Work on ``plot_me`` started as a private collection of help functions in 2008 and
has been under active development since then.
