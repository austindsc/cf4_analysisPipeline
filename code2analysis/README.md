### Code 2 organization

The new goal is as follows. 

NCA - get NCA from estimated Vover (already done)

### Alpha Data Analysis ###
1. Get list of all good alpha data. This includes:
- ?different heights?
- ?multiple voltages?
- different filters (bare, lp)

### Gamma Data Analysis ###
1. get list of all good gamma data. This includes:
- gamma sources
- background runs
- different heights
- ?multiple voltages?
- different filters (bare, lp)

2. Somehow organize this data.
3. write some usable code to analyze all the data. Output bkg corrected spectra. save them to files (remain organized)

repeat analysis with many voltage data later

## Analysis Steps ##
__LP data__
1. background + spectra building script
   1. save BKG histogram at TH1
2. run script in notebook. load .root file
3. then commence to fitting MC data

---
Variables to constrain:

alpha: Vis yield, UV yield
gamma: Vis yield, UV yield. possible UV absorption

### Alpha Simulation Comparison (comp2 directory)
Simulate different vis fraction (1D) w/ LY = 0.1 ph/keV,  
matching data for each run

1. load in all desired alpha data --> histograms (no bkg or noise subtraction)
2. loop through MC data (vis Fraction) - 1D scan
3. fit MC histograms with gaussians. scale means by NCA
4. minimize chi2 for variables (2x vertical scaling of MC to match data, 1x horizontal LY scaling)
  1. minimize chi2
  2. two vertical scalings because rates make no sense (saturated)

### Gamma Simulation Comparison

Longpass first:
Simulate data with vis frac = 100, many LY values (1D)
1. load in pre-built LP data histograms (bkg subtracted).
2. loop through MC data, scale by NCA
  1. minimize chi2 for vertical 'rate scaling' x1 per source so x2.
3. find best vis LY for LP data

Then bare:
repeat simulations with LY_vis = from above, scanning UV frac (1D)
1. load in pre-built bare data histograms (bkg subtracted).
2. loop through MC data, scale by NCA
  1. minimize chi2 for vertical 'rate scaling' x1 per source so x2.
3. find best UC LY for bare data

total is answer!

if UV LY gets higher when excluding 5.5cm data, then UV absorption is needed (or vice versa for closer data)

---
To do:
- check that alpha spectra with LP filter dont care about UV light yield (have to scale LY by Vis Frac)