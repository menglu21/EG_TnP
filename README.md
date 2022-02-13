# EG_TnP

step1: use makehist.py to make histograms for fit
step2: use eleSF.py for fit

N.B.: 
when fitting DY sample, the shape is extracted directly from DY histogram, and smearing with a GAUSS function. and no bkg is needed.

when fitting Data, the signal shape is the one from DY, and the bkg shape is defined in RooCMSShape.cc (which is used in the official TnP method)
