# 2D NLL Minimiser

A large dataset for the decay times of the D0 Meson from the LHC was processed to obtain an estimate for the mean lifetime of the particle. This was initially done using a 1D approach, and then extended to a 2D minimiser after taking into account the errors arising from false background readings.

In order to get the results as shown below, simply run the `main.py` module. This module imports all of the other modules and runs a script to produce the results. Make sure that all the files are kept in one folder, including the `lifetime.txt` file.

Remember to change the directory to the folder containing all the files so that the `lifetime.txt` file can be identified and read.

NOTE: AFTER RUNNING THE MODULE, IT CAN TAKE UP TO 35 SECONDS FOR ALL RESULTS AND PLOTS TO APPEAR.

The following information should be output by default:

```
Area under initial PDF (without background): 1.00007564437
Area under refined PDF (with background): 1.00008090162
NLL at Minimum 1D: 6220.44689279
Tau at Minimum 1D: 0.404545876571 picoseconds
Positive Tau Error 1D: 0.00473 picoseconds
Negative Tau Error 1D: 0.00467 picoseconds
Mean Tau Error 1D: 0.0047 picoseconds
Last Parabolic Estimate Error of Tau: 0.00470713686503 picoseconds
Intersect of extrapolation gives approximately 226589 readings required for accuracy of 0.001ps
NLL at Minimum 2D: 6218.39444549
Tau at Minimum 2D: 0.409705504195 picoseconds
Mean Tau Error from contour at NLLMin + 0.5: 0.00347057135406 picoseconds
a at Minimum 2D: 0.98361897218
Mean a Error from contour at NLLMin + 0.5: 0.00488295445041
Number of Background Readings out of 10,000: 164
```

Additionally, 5 figures will be produced:
1) Histogram with two PDFs overlaid
2) NLL fit for 1D, along with minimum indicated
3) Log(Error in tau) vs Log(Number of Readings)
4) 3D surface of 2D NLL (Negative Log-Likelihood) fit
5) Contour plot of 2D NLL near minimum

You can also play around with some of the methods in the `functions.py` and `minimiser.py` files to change the results.

## Example Graphs

![A histrogram showing the raw data from the lifetime.txt file, as well as two PDFs.](/images/hist.png?raw=true)
![3D visualitation of a 2D NLL fit.](/images/3d.png?raw=true)
![Contour plot of 2D NLL near minimum.](/images/contour.png?raw=true)
