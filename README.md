# Export and plot Raman spectra from Gaussian 16 output file

A Python toolkit for extracting, processing, and visualizing Raman spectroscopy data from Gaussian 16 computational chemistry output files. This package supports both resonance Raman (RR) and non-resonance Raman (NRR) calculations with multiple incident light frequencies.

## Features

- **Extract Raman data** from Gaussian 16 `.log` files
- **Support for multiple calculation types**:
  - Resonance Raman (RR) spectroscopy 
  - Non-resonance Raman (NRR) spectroscopy with multiple incident light frequencies
- **Flexible spectral broadening** with Gaussian and Lorentzian line shapes
- **Automatic spectrum generation** and CSV export
- **Interactive plotting** with matplotlib

## Installation

Clone the repository and ensure you have the required dependencies:

```bash
git clone https://github.com/Alchemist-Aloha/raman_from_gaussian.git
cd raman_from_gaussian
```

### Dependencies

- Python 3
- Jupyter Notebook (for examples)
- NumPy
- Matplotlib
- Pathlib (standard library)
- re (standard library)


## Usage

### Quick Start

Follow this gaussian document[https://gaussian.com/freq/?tabid=3#Freq_keyword__ROA_option] to set up your Gaussian input files for pre-resonance Raman or resonance Raman calculations.

```python
from core import generate_rr_spectrum, generate_nrr_spectrum

# For resonance Raman calculations
generate_rr_spectrum(
    min=0, max=4000, nstep=8000, 
    input_filename="your_rr_calculation.log",
    lorentzian_width=10, gaussian_width=10, 
    lorentzian_percentage=1
)

# For non-resonance Raman calculations with multiple incident lights
generate_nrr_spectrum(
    min=0, max=4000, nstep=8000,
    input_filename="your_nrr_calculation.log", 
    lorentzian_width=10, gaussian_width=10,
    lorentzian_percentage=1
)
```

### Functions Overview

#### Data Extraction Functions

**`extract_rr(file_path)`**
- Extracts frequency and Raman activity data from resonance Raman calculations
- Returns: NumPy array with columns [frequency, Raman_activity]

**`extract_nrr(file_path)`** 
- Extracts frequency and Raman activities for multiple incident light frequencies
- Dynamically handles any number of incident light frequencies
- Returns: NumPy array with columns [frequency, RamAct_Fr1, RamAct_Fr2, ...] and incident light frequencies

#### Spectrum Generation Functions

**`generate_rr_spectrum(min, max, nstep, input_filename, lorentzian_width, gaussian_width, lorentzian_percentage)`**
- Complete workflow for resonance Raman: extract data → broaden peaks → plot → save CSV
- Creates Lorentzian/Gaussian broadened spectrum

**`generate_nrr_spectrum(min, max, nstep, input_filename, lorentzian_width, gaussian_width, lorentzian_percentage)`**
- Complete workflow for non-resonance Raman with multiple incident wavelengths (set with rdfreq in Gaussian input)
- Creates Lorentzian/Gaussian broadened spectrum
- Generates separate spectra for each incident light frequency
- Creates individual plots and CSV files for each incident light

#### Plotting Function

**`plot_raman_spectrum(min, max, nstep, freqs, intensities, lorentzian_width, gaussian_width, lorentzian_percentage)`**
- Low-level plotting function for custom workflows
- Plot broadened spectra from Raman activity data
- Shows comparison between Gaussian, Lorentzian, and combined broadening

### Parameters

- **`min, max`**: Frequency range for the spectrum (cm⁻¹)
- **`nstep`**: Number of points in the frequency grid
- **`lorentzian_width`**: FWHM for Lorentzian broadening (cm⁻¹)
- **`gaussian_width`**: FWHM for Gaussian broadening (cm⁻¹) 
- **`lorentzian_percentage`**: Mix ratio (0=pure Gaussian, 1=pure Lorentzian)

### Output Files

The functions automatically generate:
- **Plots**: Interactive matplotlib figures showing spectral comparisons
- **CSV files**: 
  - RR: `filename_rrspec_lorentzian{width}.csv` and `filename_rrspec_gaussian{width}.csv`
  - NRR: `filename_nrrspec_lorentzian{width}.csv` and `filename_nrrspec_gaussian{width}.csv`

## Example Workflow

See `raman_from_gaussian.ipynb` for detailed examples:

## Gaussian 16 Calculation Requirements

### For Resonance Raman (RR)
- Optimize ground state geometry
- Optimize excited state geometry or calculate force with tddft at ground state geometry
- See `Freq=ReadFCHT` in Gaussian documentation for resonance Raman calculations

### For Non-Resonance Raman (NRR)  
- Use `freq=raman` with `cphf=rdfreq` with multiple incident light wavelengths at the end of the input file
- Output should contain sections like:
  ```
  Incident light (cm**-1):       0.00  12500.00  20366.60
  ... (more lines) ...
  Frequencies --  [values]
  ... (more lines) ...
  RamAct Fr= 1--  [values]
  RamAct Fr= 2--  [values] 
  RamAct Fr= 3--  [values]
  ```


## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for:
- Bug fixes
- New features  
- Documentation improvements
- Additional Gaussian output format support

## License

This project is open source with MIT license. Please cite appropriately if used in academic work.
