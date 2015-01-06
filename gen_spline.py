#!/usr/bin/python

import matplotlib.pyplot as mpl;
import numpy as np;
from scipy import interpolate;
from bisect import bisect_left;
from scipy import signal;

degree = 5;
arg = np.linspace(0.0, 2.0 * np.pi, 100);
val = np.sin(arg);
s = interpolate.splrep(arg, val, k=degree, s=0, per=1);

degree = s[2]
knots = s[0];
cp = s[1];

f = open('spline', 'w');
f.write(str(degree) + '\n');
f.write(str(knots.tolist()) + '\n');
f.write(str(cp.tolist()) + '\n');
f.close();
