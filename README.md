# OIC-293

The goal of this project is to quantify how different high-fat diets change the morphological characteristics of different organelles in a cell. The dataset used here are scans of prepared tissue slides. Each slide has two channels: a marker for the organelle of interest, and DAPI (nucleus).

This is a project with Yujin Leong (Lein Lab).

## Getting started

### Prerequisites

- [Python](https://www.python.org/downloads/) version 3.14.0

### Installation

1. Download or clone the GitHub repository
   ```bash
   git clone git@github.com:vaioic/OIC-293.git
   cd OIC-293
   ```

2. Create a python virtual environment
   ```bash
   python -m venv venv
   ```

python3.13 -m venv .venv
source .venv/bin/activate

   ```powershell
   python3.13 -m venv .venv
   .\.venv\Scripts\Activate.ps1
   ```

3. Activate the virtual environment
   ```bash
   .\venv\Scripts\activate
   ```

4. Install the dependencies using Pip
   ```bash
   python -m pip install --upgrade pip
   python -m pip install -r .\requirements.txt
   ```

### Running the code

1. Start the virtual environment if not already loaded
   ```bash
   .\venv\Scripts\activate
   python -m analysis.20260616_batch1
   ```

TBD

## Issues

If you encounter any issues with running the code or have any questions, please create an [Issue](https://github.com/vaioic/OIC-293/issues) or send an email to opticalimaging@vai.org. If you are reporting a programmatic bug, please include any error messages to aid with troubleshooting.

## Contributors

<a href="https://github.com/vaioic/OIC-293/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=vaioic/OIC-293" />
</a>

## Acknowledgements

### Dependencies

This project relies on the following packages:

**Note:** For full dependency list, see [requirements.txt](requirements.txt).

