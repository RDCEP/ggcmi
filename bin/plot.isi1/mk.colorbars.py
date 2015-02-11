#!/usr/bin/env python

import matplotlib
from numpy import arange
import matplotlib.pyplot as plt

cdict_beta = {'red': ((0.0, 0.0, 0.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 0.0, 0.0)),
              'green': ((0.0, 0.0, 0.0),
                        (0.5, 1.0, 1.0),
                        (1.0, 0.5, 0.0)),
              'blue': ((0.0, 0.0, 0.5),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0))
             }

cdict_lambda = {'red': ((0.0, 0.5, 0.5),
                        (0.5, 1.0, 1.0),
                        (1.0, 1.0, 1.0)),
                'green': ((0.0, 0.0, 0.0),
                          (0.5, 0.0, 0.0),
                          (1.0, 1.0, 1.0)),
                'blue': ((0.0, 0.0, 0.0),
                         (0.5, 0.0, 0.0),
                         (1.0, 1.0, 1.0))
               }

cdict_count = {'red': ((0.0, 0.0, 0.0),
                       (0.5, 0.375, 0.375),
                       (1.0, 0.75, 0.75)),
               'green': ((0.0, 0.75, 0.75),
                         (0.5, 0.375, 0.375),
                         (1.0, 0.0, 0.0)),
               'blue': ((0.0, 0.0, 0.0),
                        (0.5, 0.0, 0.0),
                        (1.0, 0.0, 0.0))
              }

# beta
f = plt.figure(1, figsize = (1.6, 5))
ax = f.add_axes([0.2, 0.05, 0.2, 0.9])

cmap = matplotlib.colors.LinearSegmentedColormap('BlueGreen', cdict_beta)
norm = matplotlib.colors.Normalize(vmin = -100, vmax = 100)

cb = matplotlib.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
cb.set_ticks(arange(-100, 125, 25))
cb.set_label(r'%')

plt.show()
plt.savefig('cb_beta.png')
plt.close()

# lambda
f = plt.figure(2, figsize = (1.6, 5))
ax = f.add_axes([0.2, 0.05, 0.2, 0.9])

cmap = matplotlib.colors.LinearSegmentedColormap('RedScale', cdict_lambda)
norm = matplotlib.colors.Normalize(vmin = -100, vmax = 0)

cb = matplotlib.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
cb.set_ticks(arange(-100, 20, 20))
cb.set_label(r'%')

plt.show()
plt.savefig('cb_lambda.png')
plt.close()

# delta yield
f = plt.figure(3, figsize = (1.3, 5))
ax = f.add_axes([0.2, 0.05, 0.25, 0.9])

cmap = matplotlib.cm.seismic_r
norm = matplotlib.colors.Normalize(vmin = -8, vmax = 8)

cb = matplotlib.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
cb.set_ticks(arange(-8, 10, 2))
cb.set_label(r'Gcal/ha')

plt.show()
plt.savefig('cb_delta_yield.png')
plt.close()

# count
f = plt.figure(4, figsize = (1.3, 5))
ax = f.add_axes([0.2, 0.05, 0.25, 0.9])

cmap = matplotlib.colors.LinearSegmentedColormap('RedScale2', cdict_count)
norm = matplotlib.colors.Normalize(vmin = 0, vmax = 100)

cb = matplotlib.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
cb.set_ticks(arange(0, 120, 20))
cb.set_label(r'%')

plt.show()
plt.savefig('cb_count.png')
plt.close()