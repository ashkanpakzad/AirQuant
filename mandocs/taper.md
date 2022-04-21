# Tapering Metrics

One of the ultimate goals of this intense data processing is to measure the gradient of tapering of airways. There are two 'schools of thought' built into AirQuant so far. One we will call **long tapering**, this is where you can look at the gradient of tapering from the carina to the distal point of the outermost airways, the other we call **segmental tapering**, where we consider each airway graph edge a segment as an individual unit to measure.

## Long Tapering

With this analysis we are looking at the tapering gradient from the carina to the most distal point of every airway end terminal. So we will have the same number of measurements as terminal airway end nodes (i.e nodes with one edge minus the top trachea node).

### AllTaperResults = ComputeTaperAll(obj)

We have one high level function which will compute the tapering gradient of every long airway path for all three sets of airway measurements and store it in a table. This is done relatively fast and can be saved in a new variable which the user can export to csv, xlsx etc. if they wish. Furthermore, the table is also stored within the AQ property structure `obj.Specs`.

*Notes*

This method calls on `ListTerminalNodes` to get a list of terminal nodes and process each of them through  `ConstructTaperPath`.

This is the measure that K Quan et al. considers and analysed extensively [Ref 1](https://doi.org/10.1117/12.2292306), [Ref 2](https://doi.org/10.1117/1.jmi.6.3.034003) and [Ref 3](http://arxiv.org/abs/1906.12225). Though their method is not exactly the same to process airway measurements.

*Example*

```
% reloading processed AQ object.
savename = 'results/github_demo/github_demo_AQ.m'
AQ = AirQuant(savename);
% call function
AllTaperResults = ComputeTaperAll(AQ)
```


## Segmental Tapering

This analysis is looking at the tapering of each individual airway branch as a single unit. Thus we will have the same number of measurements as there are branches minus the trachea branch. Specifically, there are two measures considered here **Intratapering** and **Intertapering**. The former considers percentage change in airway diameter along the segment relative to the first airway diameter measurement, the latter considers percentage change in average airway diameter relative to the previous segment.

### SegmentTaperResults = SegmentTaperAll(obj, prunelength)

We have one high level function which will compute the tapering gradient of every branch for all three sets of airway measurements and store it in a table. This is done relatively fast and can be saved in a new variable which the user can export to csv, xlsx etc. if they wish. Furthermore, the table is also stored within the AQ property structure `obj.Specs`.

The second argument is optional and sets the length in mm to prune either end of the airway branches, it is a two element array. This is sometimes decided upon in order to avoid natural diameter changes due to bifurcations. For no pruning set the `prunelength` to `[0 0]`.

*Notes*

This method calls on `ComputeIntraTaper`, `ComputeIntraTaperAll` and `ComputeInterTaper`.

This is a measure that [Kuo et al.](doi.org/10.1007/s00330-019-06606-w) considers extensively.

*Example*

```
% reloading processed AQ object.
savename = 'results/github_demo/github_demo_AQ.m'
AQ = AirQuant(savename);
% call function
SegmentTaperResults = SegmentTaperAll(AQ, [0 0]])
```
