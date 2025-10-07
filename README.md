<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/user-attachments/assets/3b0fb24b-c23a-44c6-93da-9739b594da17" width="200">
    <img alt="Logo" src="https://github.com/user-attachments/assets/d35a5529-2c19-4dbb-89b6-ece1fc890b8a" width="200">
  </picture>
</div>

## Overview
**PRISM** (**P**hoton-statistic **R**etrieval v**I**a **S**PAD arrays **M**odelling) is a library for simulating SPAD array photo-detection and retrieving input light statistics. This document outlines the steps to install dependencies, compile, and run the project on a Linux environment.

## Prerequisites

### Dependencies
1. **C++ Compiler**: A C++17 compatible compiler, `g++` or `clang`.
2. **CMake**: Version 3.10 or higher.
3. **YAML-CPP**: A C++ library for parsing YAML files. Install `yaml-cpp` with a package manager.
4. **Eigen**: A C++ template library for linear algebra. Copy the `Eigen` folder to `/usr/local/include`

### Additional Tools
- `git` (optional, for cloning the repository)

## Installation

### Clone the Repository
First, clone the PRISM project repository:
```bash
git clone git@github.com:ricalbr/prism.git
cd prism
```

### Install Dependencies

#### Linux
On Debian/Ubuntu-based distributions:
```bash
sudo apt-get update
sudo apt-get install build-essential cmake libyaml-cpp-dev
```

### Build Instructions
1. Create a build directory and navigate into it:
   ```bash
   mkdir build && cd build
   ```
2. Run CMake to configure the project:
   ```bash
   cmake ..
   ```
3. Compile the project:
   ```bash
   cmake --build .
   ```
4. Run the executable:
   ```bash
   ./bin/prism_simulator
   ```

### Plotting Instructions
Create a Python virtual environment and install the required packages with `pip`
```bash
python3 -m venv $HOME/.venvs/prism
source $HOME/.venvs/prism/bin/activate
pip install numpy matplotlib scipy pyyaml pyqt6
```


## Windows
Setup Linux on WSL following the instructions at the following [link](https://learn.microsoft.com/en-us/windows/wsl/install) and follow the Linux setup procedure.

## File Structure
- **`prism/`**: Contains the library source and headers.
- **`config.yaml`**: Contains the configuration file for running the simulations.
- **`main.cpp`**: The main file to run the simulation.
- **`CMakeLists.txt`**: Configuration files for building the project.

## Notes
* SPAD measuring matrix comprising a finite number of detectors ($D$), SPAD efficiency ($\eta$), and dark-count rate ($k$) is retrieved in a novel recursive approach. Similar works:
  - Miatto et al., "Explicit formulas for photon-number discrimination with on/off detectors" [https://doi.org/10.1364/AO.57.006750](https://doi.org/10.1364/AO.57.006750)
  - Nehra et al., "Photon-number-resolving segmented detectors based on single-photon avalanche-photodiodes" [https://doi.org/10.1364/OE.380416](https://doi.org/10.1364/OE.380416)
* Cross-talk is modeled as in Afek et al., "Quantum state measurements using multipixel photon detectors" [https://doi.org/10.1103/PhysRevA.79.043830](https://doi.org/10.1103/PhysRevA.79.043830).
  - Note: [https://doi.org/10.1364/OE.15.014539](https://doi.org/10.1364/OE.15.014539) adopt a recursive approach for the computation of the cross-talk.
* The EME algorithm is a slightly modified version of the code given in Hlousek et al., "Accurate detection of arbitrary photon statistics" [https://doi.org/10.1103/PhysRevLett.123.153604](https://doi.org/10.1103/PhysRevLett.123.153604).

