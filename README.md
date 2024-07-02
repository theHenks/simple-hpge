# Simple HPGe
This packages provides simple scripts to analyze HPGe data. It is based on the [Juleana](https://github.com/legend-exp/Juleana.jl) sofwtare stack and enables the user to perform basic data analysis tasks such as energy calibration, peak fitting, and efficiency calibration including Digital Signal processing steps.
## Installation
To install the package, simply clone the repository and install the required packages. The following commands will install the package and the required dependencies:
``` julia
using Pkg
pkg"registry add General https://github.com/legend-exp/LegendJuliaRegistry.git"
Pkg.activate("./")
Pkg.instantiate()
```
## Requirements
All data needs to be able to be read by the [LegendHDF5IO](https://github.com/legend-exp/LegendHDF5IO.jl) package. The data should be stored in a HDF5 file with the following structure. At the moment only FlashCam and Struc ADCs are supported. However, other ADCs can be added by extending the `LegendHDF5IO` package and parsing DAQ raw data into the required format.
In addition, the scripts require a run structure where the DAQ data is taken as a run consisiting of one or multiple `HDF5` files in a folder. The user can specify the run folder and further `IO` folders in the `data_config.json` file.


## How-to-use
All scripts can be found in the `src/` folder and can be executed individually. However, there is an intrinsic order in which the scripts should be executed. The following list provides an overview of the scripts and their order of execution:
1. `split_peak.jl`: This script splits the data into individual peaks of interest. It analyzes all waveforms and extracts the peaks of interest. The output is saved in peaks file containing the waveforms in a specific key. 
2. `decay_time.jl`: This script performs the decay time analysis of the HPGe detector. It extractts the decay for a subset of waveforms to determine the decay time of the Charge Sensitive Amplifier (CSA). The out is saved in a `pars.json` file.
3. `optimization.jl`: This script performs the optimization of the energy calibration. It uses the peaks file and the `pars.json` file to optimize the energy filters used in the DSP. The user can specificy the filter to be optimized within the script. The output is saved in a `pars.json` file.
4. `dsp.jl`: This script performs the DSP steps on the data. It uses the `pars.json` file to apply the optimized filters to the data. The output is saved in a `*_dsp.lh5` files. The full DSP is optimized to be run on a single server machine in one or multiple `julia` processes. The user can specify the number of processes in the script.

**Attention**: All further scripts are at the moment WIP...

## Configuration
Vairous configuration files are needed to perform the indivual steps mentioned above. Examples for that can be found in the `config/` folder. 
Depending on the waveform input data and the used ADCs, the user needs to specify the configuration files accordingly. The `config/` folder contains the following files:
- `dsp_config.json`: This file contains the configuration for the DSP steps. The user can specify the filters to be used in the DSP and the optimization steps.
- `energy_config.json`: This file contains the configuration for the energy calibration. The user can specify the energy calibration parameters.
- `qc_config.json`: This file contains the configuration for the quality cuts. The user can specify the quality cut parameters.

