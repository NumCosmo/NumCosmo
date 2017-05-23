#!/usr/bin/python2 

from astropy.io import fits
from plot2DCorner import plotCorner
import corner
import matplotlib.pyplot as pl
import numpy as np
import sys
import os.path

sys.argv.pop (0)

files = []
col_name_att = {}

while len (sys.argv) > 2:

  fitsfile      = sys.argv.pop (0)
  burning_in    = int (sys.argv.pop (0))
  high_prec     = int (sys.argv.pop (0))
  (prefix, ext) = os.path.splitext (fitsfile)
  pindex = []
  if (len (sys.argv) > 0):
    pindex = map (int, sys.argv)
    sys.argv = []

  print ("# Opening: %s." % fitsfile)
  hdulist = fits.open (fitsfile)
  cols = hdulist[1].columns
  ncols = len (cols)
  col_names = []
  col_symbols = []
  col_data = []
  padding_cols = hdulist[1].header["NADDVAL"]
  print ("# Number of columns %d padding columns %d burning in %d." % (ncols, padding_cols, burning_in))
  if (len (pindex) == 0):
    pindex = range (ncols)

  for i in pindex:
    col = hdulist[1].data.field (i)
    col = col[burning_in:]
    col_data.append (col)
    col_name = cols.names[i]
    col_names.append (col_name)
    col_symbol = None
    
    if i < padding_cols:
      col_symbol = r'$%s$' % hdulist[1].header["ASYMB%d" % (i + 1)]
    else:
      col_symbol = r'$%s$' % hdulist[1].header["FSYMB%d" % (i - padding_cols + 1)]
    col_symbols.append (col_symbol)

    epsilon = 1.e-7
    col_min = np.percentile (col, epsilon) # range is defined computing the 6\sigma confidence interval
    col_max = np.percentile (col, 100.0 - epsilon) 

    #col_min = min (col)
    #col_max = max (col)
    #print col_min, col_max
    interval = col_max - col_min
    col_min = col_min - interval * 0.03
    col_max = col_max + interval * 0.03
    
    if col_name in col_name_att:
      cur_min = col_name_att[col_name]['min']
      cur_max = col_name_att[col_name]['max']      
      col_name_att[col_name] = {'min' : min (col_min, cur_min), 'max' : max (col_max, cur_max)}
    else:
      col_name_att[col_name] = {'min' : col_min, 'max' : col_max}

  col_data = np.asarray (col_data)
  col_data = np.matrix.transpose (np.asarray (col_data))
  col_data = col_data[~np.all (np.abs (col_data) < 1e-6, axis=1)]
  print col_data.shape
  
  files.append ({'fitsfile' : fitsfile, 'pad' : padding_cols, 'prefix' : prefix, 'ncols' : len (pindex), 'cols' : col_data, 'col_names' : col_names, 'col_symbols' : col_symbols, 'high_prec' : high_prec})
  hdulist.close()
  del hdulist
  del fitsfile
  del padding_cols
  del prefix
  del ext

for file in files:
  scatterstyle = {'color' : 'r', 'alpha' : 0.5}
  styleargs    = {'color' : 'k', 'scatterstyle' : scatterstyle}
  padding_cols = file['pad']
  ncols        = file['ncols']
  cols         = np.array (file['cols'])[:,padding_cols:]
  col_symbols  = file['col_symbols'][padding_cols:]

  params_range = []
  for i in range (padding_cols, ncols):
    param_lim = [col_name_att[file['col_names'][i]]['min'], col_name_att[file['col_names'][i]]['max']]
    params_range.append (param_lim)
      
  #print params_range
  #print ncols
  #print len (file['cols'][0,])
  #print len (cols)
  #print len (col_symbols)
  #print len (params_range)
  #print file['cols'].shape
  #print params_range

  fig = None
  if file['high_prec'] == 0:  
    fig = corner.corner (cols, labels = col_symbols, range = params_range, show_titles=True)
  else:
    fig = plotCorner (cols, bins = 50, labels = col_symbols, label_size = 20, range = params_range)

  figfile = "%s_corner.pdf" % (file['prefix'])
  print ("#  Saving figure: %s." % figfile)
  pl.savefig (figfile)

  pl.close ()
