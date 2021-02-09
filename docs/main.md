# AirQuant Docs
By Ashkan Pakzad https://ashkanpakzad.github.io

See [/readme.md](/readme.md) before reading this document regarding install, about and purpose.


This software primarily revolves around the [library/AirQuant.m] class. Every CT case will exist as a new AirQuant class object. When a case is loaded in, it is first initialised, performing a number of checks and short processes.
Currently there is only one method for taking airway measurements, the FWHM_esl method [ref]. This first relies on the CT scan to be interpolated at right angles to the principle airway axis. For reading of the concept, please see K Quan et al. [ref] on which this software is based.

## Initialisation
Extracting the skeleton, and parsing every airway into a graph tree structure.

### Skeletonisation

### Lobe Classification

## Traversing the airway: Interpolating the CT

## Full Width at Half Maximum - Edge-cued Segmentation Limited Technique

## Visualisation methods

## Metrics
### Long airway tapering
### Airway segment tapering
