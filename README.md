AirwaySkel is a MATLAB class for extracting airway measurements from a fully segmented chest CT of the desired airways.

## Dependencies
* MATLAB
* MATLAB Signal Processing toolbox

## Install
* In Terminal:
```
$ git clone https://github.com/ashkanpakzad/AirQuant.git
```
* Open MATLAB and navigate into AirQuant folder.
* In MATLAB Command Window:
```
% Adding AirQuant to MATLAB path
current_path = pwd;
addpath(genpath(current_path))
```

## Use


## Packaged Dependencies
* skel2graph3d-matlab (https://github.com/phi-max/skel2graph3d-matlab)
* ellipse (https://uk.mathworks.com/matlabcentral/fileexchange/289-ellipse-m)

## Acknowledgments
* Based on Kin Quan "Airway Tapering In CT" https://github.com/quan14/AirwayTaperingInCT [1]
* skel2graph3d-matlab by Philip Kollmannsberger [2]
* ellipse.m by David Long [3]

## References
[1] K. Quan et al., “Tapering analysis of airways with bronchiectasis,” Proc. SPIE 10574, 105742G (2018). https://arxiv.org/abs/1909.06604

[2] Kollmannsberger, Kerschnitzki et al., "The small world of osteocytes: connectomics of the lacuno-canalicular network in bone." New Journal of Physics 19:073019, 2017.

[3]  David Long (2020). ellipse.m (https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m), MATLAB Central File Exchange. Retrieved January 21, 2020.
