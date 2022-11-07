Tutorials
=========

**Basic**

.. code-block:: matlab

    AirQuantAddPath();
    AQdownload_data('chestct');
    skel = logical(niftiread('chestct_airway_PTKskel.nii.gz'));
    AQnet = ClinicalAirways(skel);


Tutorials
---------

.. grid-item-card:: Quickstart: Airways

  .. button-link:: _static/CA_quickstart.html
      :color: dark
      :align: left
      :expand:

      Open tutorial

  .. button-link:: https://raw.githubusercontent.com/ashkanpakzad/AirQuant/master/tutorials/CA_quickstart.m
      :color: secondary
      :align: right
      :expand:

      Open raw


.. grid-item-card:: Measure airways using FWHMesl

  .. button-link:: _static/CA_FWHMesl.html
      :color: dark
      :align: left
      :expand:

      Open tutorial

  .. button-link:: https://raw.githubusercontent.com/ashkanpakzad/AirQuant/master/tutorials/CA_FWHMesl.m
      :color: secondary
      :align: right
      :expand:

      Open raw
