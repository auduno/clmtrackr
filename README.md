clmtrackr
======

clmtrackr is a javascript library for precise tracking of facial features in videos or images. It currently is an implementation of *constrained local models* fitted by *regularized landmark mean-shift*, as described in [Jason M. Saragih's paper](http://dl.acm.org/citation.cfm?id=1938021). Due to the heavy calculations involved, this library requires WebGL with support for floating-point textures.

The model currently used in this library is based on the model provided by [Xiaoguang Yan](http://ccxgyan.wordpress.com/) in his [implementation of Constrained Local Models](https://sites.google.com/site/xgyanhome/home/projects/clm-implementation), which he kindly allowed us to use. The aim is to provide a more generic model at some point.

The library requires ccv.js (for initial face detection) and parts of google closure library (for matrix math).

For some more information about Constrained Local Models, take a look at Xiaoguang Yan's excellent tutorial, which was of great help in implementing this library.

License
=======

clmtrackr is distributed under the [MIT License](http://www.opensource.org/licenses/MIT)