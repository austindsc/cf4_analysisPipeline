# Data-Monte Carlo Comparison Workflow

## Overview
This directory contains the complete data-MC comparison system for CF4 bubble chamber analysis. The workflow processes experimental data and Monte Carlo simulations to extract physics parameters through chi-squared fitting.

## Core Classes Architecture

### Data Processing Pipeline
```
Raw Data → newProcVars.py → Reduced Data
Reduced Data → RunClass  → RealDataSet → rcMCpair → Comparator → Chi2 Analysis
Raw G4       → FullMCRun → MCdataSet   ↗
```

RealDataSet is not necessary?

## Class Descriptions

### `RunClass` (from anaStructs.py)
**Purpose**: Core experimental data analysis container (loads reduced data)
- **Key Attributes**: 
  - `tt`: ROOT TTree with processed/reduced data
  - `mBin_rawCounts`, `mBin_edges`: Histogram data
  - `mHeight`: Source height (cm)
- **Key Methods**: 
  - `GetRate()`: Calculate event rates
  - `GetCorrSpectrum()`: Background-subtracted spectrum
  - `print()`: Display run summary

### `RealDataSet` (helpers/RealDataSet.py)
**Purpose**: Container for all experimental runs grouped by source/filter type
- **Key Attributes**: 
  - `mCsBare[]`: Cesium-137 bare runs
  - `mCsBLP[]`: Cesium-137 LP-filtered runs  
  - `mCoBare[]`: Cobalt-57 bare runs
  - `mCoLP[]`: Cobalt-57 LP-filtered runs
- **Key Methods**: 
  - `findRunByHeight()`: Retrieve run by source height
  - `print()`: Tabulated run summary with rates

### `FullMCRun` (from G4_runClasses.py)
**Purpose**: Monte Carlo simulation run container
- **Key Attributes**: 
  - `mRunPars`: Simulation parameters (height, source, filter)
  - `mRaw_Counts`: Raw histogram counts
  - `Nbeam`: Number of simulated photons
- **Key Methods**: 
  - `doCAcorrection_*()`: Apply charge avalanche corrections

### `MCdataSet` (helpers/MCdataSet.py)
**Purpose**: Container for all MC runs with automatic sorting
- **Key Attributes**: 
  - `mCsBareList[]`: Sorted Cs-137 bare simulations
  - `mCoBareList[]`: Sorted Co-57 bare simulations
  - `gammaMap{}`: Dictionary mapping for easy access
- **Key Methods**: 
  - `buildOutDataArrays()`: Sort runs by source/filter type
  - `findRunByHeight()`: Retrieve MC run by height

### `rcMCpair` (helpers/Comparator.py)
**Purpose**: Pairs experimental run with corresponding MC simulation
- **Key Attributes**: 
  - `rc`: RunClass (experimental data)
  - `mc`: FullMCRun (MC simulation)
  - `height`: Source height
  - `firstBin`, `lastBin`: Analysis range
  - `ndof`: Degrees of freedom
- **Key Methods**: 
  - `SetBinRanges()`: Determine valid bin range (>10 counts)
  - `GuessScaling()`: Analytical scaling estimate
  - `computeChi2_reduced()`: Calculate chi-squared

### `Comparator` (helpers/Comparator.py)
**Purpose**: Manages all data-MC pairs and performs global fitting
- **Key Attributes**: 
  - `mRealData`: RealDataSet object
  - `mMCDataSet`: MCdataSet object
  - `csBareMap[]`: Cs-137 data-MC pairs
  - `coBareMap[]`: Co-57 data-MC pairs
  - `csGuess`, `coGuess`: Scaling parameters
- **Key Methods**: 
  - `do_mapping_fromJson()`: Create pairs from config mapping
  - `SetBinRanges()`: Set analysis ranges for all pairs
  - `GuessOverallScaling()`: Calculate initial scaling estimates
  - `getCsBarechi2()`, `getCoBarechi2()`: Compute chi-squared by source

## Configuration System

### config.json Structure
```json
{
  "global_analysis": {
    "csBareMap": {"mc_height": data_height, ...},
    "coBareMap": {"mc_height": data_height, ...},
    "binLow": 0.5, "binHigh": 100.5, "nBins": 50
  },
  "mc_analysis": {
    "NCA": 0.17  // Charge avalanche parameter
  },
  "data_analysis": {
    "qvar": "qFlat",  // PE calculation method
    "peCut": 0.5      // PE threshold
  }
}
```

## Workflow Steps

1. **Raw Data Processing** (`newProcVars.py`)
   - Process raw experimental data files
   - Apply calibrations and corrections
   - Generate reduced data files

2. **Data Loading** (`makeDataObj_cs.py`)
   - Load reduced experimental ROOT files
   - Create RunClass objects with PE normalization
   - Group into RealDataSet

3. **MC Loading**
   - Load Geant4 ROOT files
   - Create FullMCRun objects
   - Group into MCdataSet with automatic sorting

4. **Pairing** (`Comparator.do_mapping_fromJson()`)
   - Map MC heights to data heights via config
   - Create rcMCpair objects

5. **Analysis Range Setting**
   - Find valid bin ranges (>10 counts in both data/MC)
   - Calculate degrees of freedom

6. **Scaling Estimation**
   - Analytical guess based on integrated counts
   - Account for runtime differences

7. **Chi-Squared Calculation**
   - Per-pair chi-squared with Poisson statistics
   - Global chi-squared for parameter fitting

## Key Files

- **newProcVars.py**: Raw data processing script (creates reduced data)
- **c1_gammasBare.ipynb**: Main analysis notebook
- **makeDataObj_cs.py**: Reduced data loading script  
- **helpers/Comparator.py**: Core comparison classes
- **helpers/RealDataSet.py**: Experimental data container
- **helpers/MCdataSet.py**: MC data container
- **config.json**: Analysis configuration

## Usage Example

```python
# Load data and MC
realData = RealDataSet()
mcData = MCdataSet("/path/to/mc/folder")

# Create comparator and map runs
comp = Comparator(realData)
comp.loadNewG4(mcData)
comp.do_mapping_fromJson()

# Set analysis ranges and compute chi2
comp.SetBinRanges(firstBin=1)
comp.GuessOverallScaling()
chi2_cs = comp.getCsBarechi2()
chi2_co = comp.getCoBarechi2()
```

## Dependencies

- **ROOT/PyROOT**: Histogram handling and file I/O
- **numpy**: Numerical operations and chi-squared calculations
- **json5**: Configuration file parsing
- **scipy.optimize**: Parameter fitting (future use)

## Notes

- All heights in cm, with data heights measured from chamber bottom
- MC scaling accounts for runtime differences and beam intensity
- LP filter blocks <490nm UV while transmitting 95% visible light
- Charge avalanche corrections applied to MC using VUV4-Q-50 SiPM parameters
