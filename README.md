# Nucleic Acid Analysis Application

[![tests](https://github.com/chelseadack/Nucleic-Acid-Analysis-Application/actions/workflows/tests.yml/badge.svg)](https://github.com/chelseadack/Nucleic-Acid-Analysis-Application/actions/workflows/tests.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

This is a Python-based graphical user interface for the NUPACK<sup>1,2</sup> package, designed to aid nucleic acid sequence design and analysis. As demonstrated below, it allows users to input nucleic acid sequences and parameters, run thermodynamic simulations, and visualise the resulting structures via NucDraw<sup>3</sup>.

Developed by Chelsea Dack (chelsea.dack.24@ucl.ac.uk), PhD student in the [Booth group](https://boothlab.uk/).

## How the model works
It predicts how a population of DNA or RNA strands will fold and bind together in a test tube at equilibrium, showing you not just the single most stable shape, but the probability of all shapes and the final concentration of every complex.
1. You give the model nucleic acid sequences, their concentrations, reaction temperature and salt condition. It then predicts what structures individual strands/complexes could form.
2. For each structure, it uses a nearest neighbour thermodynamic model to calculate the free energies of individual structural motifs. The parameters for this are derived experimentally e.g. by UV melting curves to derive entropy and enthalpy, and hence the free energies.
3. It sums these motif free energies to get the free energy of that specific structure.
4. It uses dynamic programming algorithms to efficiently calculate the partition function of these structures. It uses that and the Boltzmann distribution to find the probability that a base pairing exists.
5. It upscales this analysis, calculating the partition coefficient of each strand/complex and uses free energy minimisation to determine the concentrations of these at equilibrium.
6. The interface this script encodes allows non-programmers to access the model. It will return a bar chart figure of the dominant complexes that are most likely to form at equilibrium. You can then choose a complex to visualise its most likely structure.

## In Development
- Integration with Multistrand<sup>4,5</sup> for kinetic simulations
- Update to latest NUPACK<sup>1,2</sup> version to enable DNA/RNA hybrids, OMe RNA, and wider salt conditions

## Installation
1. First, make sure the following requirements are installed on your host system:
   - Python 3.10+
   - [NUPACK 4.0.2+](https://www.nupack.org/download/overview) (see NUPACK installation instructions [here](https://docs.nupack.org/start/#installation-requirements))
2. Install the other Python packages (matplotlib, PySide6, nucdraw, pytest):
    - `pip install -r requirements.txt`

## Use
To launch Pack on MacOS:
   - Open a new Terminal window and type this line, substituting the location of the GUI_run_script.py script (`GUI_logic.py` must be in the same folder):
     - `python3 ~/insert_path_to/GUI_run_script.py`
   - Note: A standalone application (not included here) has also been developed for MacOS. This allows users to launch the GUI without using Terminal. The NUPACK module still requires separate installation however, due to licence restrictions.

To launch the GUI on Windows (using Ubuntu, assuming NUPACK installation as described [here](https://docs.nupack.org/start/#installation-requirements)):
   - Move the GUI_run_script.py and GUI_logic.py scripts into Ubuntu's Linux home directory
   - Open a new Ubuntu window and type this line:
     - `python3 GUI_run_script.py`

To analyse nucleic acid sequences:
1. Upon launching the script, this window will appear, where you can input your parameters and nucleic acid sequence(s):
<img src="example_imgs/user_interface.png" width="250" height="350">
2. Click 'Run Analysis' for calculation of equilibrium concentrations (via NUPACK<sup>1,2</sup>). This window will appear:
<img src="example_imgs/example_equilibrium.png" width="600" height="150">
3. Click 'Plot Structure' for structure visualisation of a chosen complex (using NucDraw<sup>3</sup>). This window will appear, where individual bases / strand backbones are coloured by identity and base pairs are coloured by MFE probability:
<img src="example_imgs/example_structure.png" width="400" height="300">
Figure appearance can be customised, for example by changing base size, structure rotation, and figure zoom.

## Project structure
| File | Purpose |
|---|---|
| `GUI_run_script.py` | The GUI: user input, NUPACK analysis, equilibrium bar chart and NucDraw structure plot |
| `GUI_logic.py` | GUI-independent logic: sequence validation and conversion, complex labels, dot-bracket parsing |
| `test_GUI_logic.py` | Unit tests for `GUI_logic.py` |
| `example_imgs/` | Screenshots used in this README |

Input is checked before it is sent to NUPACK. Empty sequences, invalid bases, repeated strand names and names containing `+ ( ) < >` are reported in a pop-up instead of crashing the app or giving ambiguous complex labels.

## Running the tests
The core logic is covered by unit tests, which run automatically on GitHub Actions on every push. To run them locally:
```
python3 -m pytest
```

## References:
1. M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed. ACS Synth Biol, 9:2665-2678, 2020.
2. R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce. Thermodynamic analysis of interacting nucleic acid strands. SIAM Rev, 49:65-88, 2007.
3. 10.5281/zenodo.15352138
4. Joseph M. Schaeffer (2013) Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. Dissertation (Ph.D.), California Institute of Technology. [https://thesis.library.caltech.edu/7460/](https://thesis.library.caltech.edu/7460/)
5. Joseph M. Schaeffer, Chris Thachuk, Erik Winfree (2015) Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. DNA Computing and Molecular Programming (DNA21), Lecture Notes in Computer Science (LNCS) volume 9211, pp 194-211. [https://doi.org/10.1007/978-3-319-21999-8_13](https://doi.org/10.1007/978-3-319-21999-8_13)
