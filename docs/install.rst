Installation
============
See `Latest Releases <https://github.com/ashkanpakzad/AirQuant/releases>`__:

.. code-block:: RST

  $ git clone https://github.com/ashkanpakzad/AirQuant.git

**For skeletonisation** dependency, clone this specific version of `PTK <https://github.com/ashkanpakzad/pulmonarytoolkit/releases/tag/ForAirQuant1.0>`__:

.. code-block:: RST

  $ git clone https://github.com/ashkanpakzad/pulmonarytoolkit.git
  $ cd pulmonarytoolkit
  $ git checkout nifti-gz-load-import

**MATLAB path** configuring your ```startup.m``` file to find AirQuant and PTK. Run in matlab:

.. code-block:: RST

  >>> edit ~/Documents/MATLAB/startup.m

Add lines to run the addpath functions of each package, edit paths to the correct locations.

.. code-block:: RST

  >>> run ~/AirQuant/AirQuantAddPath.m
  >>> run ~/pulmonarytoolkit/PTKAddPaths.m

Save and run ```startup.m``` or restart MATLAB to take effect.

Required Dependencies
---------------------
* MATLAB 2022a (>9.12)
* MATLAB Signal Processing Toolbox (>8.3)
* MATLAB Image Processing Toolbox (>11.0)
* MATLAB Curve Fitting Toolbox (>3.5.10)
* MATLAB Statistics and Machine Learning Toolbox (>12.0)

.. todo:
  * Quick check tool for dependencies.

Optional Dependencies
---------------------

* `Pulmonary Toolkit <https://github.com/ashkanpakzad/pulmonarytoolkit/releases/tag/ForAirQuant1.0>`__ - For in-built skeletonisation algorithm.


Packaged Dependencies
---------------------

Please note that each package has its own licenses.

* `skel2graph3d-matlab`_ - convert binary skeleton images into graphs.

* `ellipse`_ - robustly fit ellipses to points.

* `linspecer`_ - beautiful contrasting colours.

* `MatlabProgressBar`_ - tqdm progress style feedback.

* `skeleton3d-matlab`_ - A 3D thinning skeletonisation algorithm.

* `ellipsePerimeter`_ - ellipse perimeter computation.

.. _skel2graph3d-matlab: https://github.com/phi-max/skel2graph3d-matlab/releases/tag/v1.2
.. _ellipse: https://www.mathworks.com/matlabcentral/fileexchange/289-ellipse-m
.. _linspecer: https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
.. _MatlabProgressBar: https://www.mathworks.com/help/matlab/ref/waitbar.html
.. _skeleton3d-matlab: https://github.com/phi-max/skeleton3d-matlab/releases/tag/v1.1
.. _ellipsePerimeter: https://www.mathworks.com/matlabcentral/fileexchange/66647-ellipse-perimeter
