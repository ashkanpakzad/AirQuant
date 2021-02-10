AirQuant is a software tool based in MATLAB primarily for extracting airway measurements from fully segmented airways of a chest CT.

## Dependencies
* MATLAB 2020b (>9.9)
* MATLAB Signal Processing Toolbox (>8.3)
* MATLAB Image Processing Toolbox (>11.0)
* MATLAB Curve Fitting Toolbox (>3.5.10)

## Install
* In Terminal:
```
$ git clone https://github.com/ashkanpakzad/AirQuant.git
```
* To use recommended skeletonisation/centreline algorithm, please install Pulmonary Toolkit, by Tom Doel. https://github.com/tomdoel/pulmonarytoolkits.
* Example use: see [scripts/example.m](scripts/example.m)

## Development
See [docs/main.md](docs/main.md).

## Packaged Dependencies
* skel2graph3d-matlab (https://github.com/phi-max/skel2graph3d-matlab)
* ellipse (https://uk.mathworks.com/matlabcentral/fileexchange/289-ellipse-m)
* linspecer (https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap)

## Acknowledgments
* Based on Kin Quan "Airway Tapering In CT" https://github.com/quan14/AirwayTaperingInCT [1]
* skel2graph3d-matlab by Philip Kollmannsberger [2]
* ellipse.m by David Long [3]
* linspecer.m by Jonathan C. Lansey [4]
* Tom Doel's Pulmonary Toolkit [5]

## Sample Data Acknowledgments
The CT images provided from QIN LUNG CT (https://wiki.cancerimagingarchive.net/display/Public/QIN+LUNG+CT#b8d88cce4fd14620bef4e5e35ec3d589) under Creative Commons Attribution 3.0 Unported License (https://creativecommons.org/licenses/by/3.0/). The citations are:

* Goldgof, Dmitry, Hall, Lawrence, Hawkins, Samuel, Schabath, Matthew, Stringfield, Olya, Garcia, Alberto, … Gillies, Robert. (2015). Data From QIN_LUNG_CT. The Cancer Imaging Archive. http://doi.org/10.7937/K9/TCIA.2015.NPGZYZBZ
* Jayashree Kalpathy-Cramer, Sandy Napel, Dmitry Goldgof, Binsheng Zhao. QIN multi-site collection of Lung CT data with Nodule Segmentations.  http://dx.doi.org/10.7937/K9/TCIA.2015.1BUVFJR7
* Clark K, Vendt B, Smith K, Freymann J, Kirby J, Koppel P, Moore S, Phillips S, Maffitt D, Pringle M, Tarbox L, Prior F. The Cancer Imaging Archive (TCIA): Maintaining and Operating a Public Information Repository, Journal of Digital Imaging, Volume 26, Number 6, December, 2013, pp 1045-1057. (https://doi.org/10.1007/s10278-013-9622-7)

## References
[1] K. Quan et al., “Tapering analysis of airways with bronchiectasis,” Proc. SPIE 10574, 105742G (2018). https://arxiv.org/abs/1909.06604

[2] Kollmannsberger, Kerschnitzki et al., "The small world of osteocytes: connectomics of the lacuno-canalicular network in bone." New Journal of Physics 19:073019, 2017.

[3]  David Long (2020). ellipse.m (https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m), MATLAB Central File Exchange. Retrieved January 21, 2020.

[4]  Jonathan C. Lansey (2020). Beautiful and distinguishable line colors + colormap (https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap), MATLAB Central File Exchange. Retrieved July 24, 2020.

[5] Pulmonary Toolkit, Tom Doel. https://github.com/tomdoel/pulmonarytoolkits
