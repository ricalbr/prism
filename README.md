<p align="center">
  <img src="https://github.com/user-attachments/assets/d35a5529-2c19-4dbb-89b6-ece1fc890b8a" alt="prism_logo" align="center" width="200"/>
</p>

## Overview
**PRISM** (**P**hoton-statistic **R**etrieval v**I**a **S**PAD arrays **M**odelling) is a library for simulating SPAD array photo-detection and retrieving input light statistics. This document outlines the steps to install dependencies, compile, and run the project on Linux and Windows.

## Prerequisites

### Dependencies
1. **C++ Compiler**: A C++17 compatible compiler is required.
   - On Linux: `g++` or `clang`
   - On Windows: Microsoft Visual Studio (MSVC) with C++ tools
2. **CMake**: Version 3.10 or higher.
3. **YAML-CPP**: A C++ library for parsing YAML files.
   - On Linux: Install via your package manager (e.g., `apt`, `dnf`)
   - On Windows: Install via vcpkg or build from source

### Additional Tools
- `git` (optional, for cloning the repository)

## Installation

### Clone the Repository
First, clone the PRISM project repository:
```bash
git clone <repository_url>
cd <repository_name>
```

### Install Dependencies

#### Linux
On Debian/Ubuntu-based distributions:
```bash
sudo apt-get update
sudo apt-get install build-essential cmake libyaml-cpp-dev
```
On Red Hat/Fedora-based distributions:
```bash
sudo dnf install gcc-c++ cmake yaml-cpp-devel
```

#### Windows
1. Install [Visual Studio](https://visualstudio.microsoft.com/) with the C++ development workload.
2. Install CMake from [cmake.org](https://cmake.org/).
3. Install `yaml-cpp` using [vcpkg](https://vcpkg.io/):
   ```bash
   vcpkg install yaml-cpp
   ```

Ensure `vcpkg` integration is enabled by running:
```bash
vcpkg integrate install
```

### Build Instructions

#### Linux
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
   ./prism_main
   ```

#### Windows
1. Open a `Developer Command Prompt` for Visual Studio or use `PowerShell`.
2. Create a build directory and navigate into it:
   ```bash
   mkdir build && cd build
   ```
3. Run CMake to configure the project:
   ```bash
   cmake .. -G "Visual Studio 16 2019" -A x64
   ```
4. Open the generated `prism.sln` file in Visual Studio and build the solution, or build from the command line:
   ```bash
   cmake --build .
   ```
5. Run the executable:
   ```
   prism_main.exe
   ```

## File Structure
- **`prism/`**: Contains the library source and headers.
- **`main.cpp`**: The main file to run the simulation.
- **`CMakeLists.txt`**: Configuration files for building the project.

## Notes
- Use the `-DCMAKE_BUILD_TYPE=Release` flag with CMake to enable optimizations for production builds.
- For debugging, use `-DCMAKE_BUILD_TYPE=Debug`.

## License
This project is distributed under the [MIT License](LICENSE).

For issues or contributions, feel free to open a pull request or file an issue in the repository. Happy coding!


