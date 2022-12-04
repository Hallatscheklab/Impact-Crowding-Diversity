# Impact-Crowding-Diversity

## Introduction
This repository contains codes, data, and CAD blueprints for the manuscript: Carl F. Schreck, Diana Fusco, Yuya Karita, Stephen Martis, Jona Kayser, Marie-Cécilia Duvernoy, Oskar Hallatschek, "Impact of crowding on the diversity of expanding populations", bioRxiv, doi: https://doi.org/10.1101/743534.
Contact ohallats[at]berkeley.edu if you have any questions.

## CAD descriptions
CAD blueprints were drawn with QCAD. The leftmost two holes are for midium inlets, the middle hole is for a cell inlet, and the right hole is for a waste outlet.

## Data descriptions
### Microfluidic experiments
All the data were taken by Olympus IX81 microscope with a 10x objective every 10 minutes. Pixel-to-um conversion is 1 um = 1.55 pixel.
### Plate experiments
Colony collision experiments were recorded by Zeiss Axio Zoom microscope. Luria-Delbrück experiments were recorded by a digital camera.
## Codes desciptions
### Matlab codes
"Drift_Correct.m" is a code for correcting the drift of a microscope stage, used for data pre-processing. "Semi_manual_tracking.m" is a code for tracking the growth and displacement of color-switched cells.
### Agent-based simulaitons
The files named "budding.damped.linear_front.serial_mut.*" correspond to agent-based simulations where cells grow via budding, the colony has periodic boundary conditions in the x-direction, and mutations are introduced at a regular time-interval. The file with extension ".f" is the Fortran 90 code, the file with extension ".o" is the corresponding executable compiled with the Intel compiler ifort, and the file with extension ".sh" is a wrapper for the exectuable to simplify running the code. We provide two example scripts with input parameters: run.production.sh corresponds to a run with 100 neutral mutations for data-collection purposes, and run.movie.sh, corresponds to a run with 10 neutral mutation for the purposes of creating a movie of the simulation (with frequent system configuration outputs). The output of these programs is in the folder outs/

