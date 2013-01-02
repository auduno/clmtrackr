clmtrackr
======

![tracked face](http://auduno.github.com/clmtrackr/media/clmtrackr01.jpg)

clmtrackr is a work-in-progress javascript library for precise tracking of facial features in videos or images. It currently is an implementation of *constrained local models* fitted by *regularized landmark mean-shift*, as described in [Jason M. Saragih's paper](http://dl.acm.org/citation.cfm?id=1938021). 

The library provides a generic face model that was trained on [the MUCT database](http://www.milbo.org/muct/). [Xiaoguang Yan](https://sites.google.com/site/xgyanhome/my-information) also kindly allowed us to use the trained model from [his implementation of Constrained Local Models](https://sites.google.com/site/xgyanhome/home/projects/clm-implementation). The aim is to provide a model that is trained on a larger set of facial expressions, as well as a model builder for building your own facial models.

The library requires [jsfeat](https://github.com/inspirit/jsfeat) (for initial face detection) and [numeric.js](http://numericjs.com) (for matrix math).

For some more information about Constrained Local Models, take a look at Xiaoguang Yan's [excellent tutorial](https://sites.google.com/site/xgyanhome/home/projects/clm-implementation/ConstrainedLocalModel-tutorial%2Cv0.7.pdf?attredirects=0), which was of great help in implementing this library.

The current version of the library implements two different types of filters for the landmarks, either SVM or MOSSE filters. Due to the heavy calculations involved in the SVM filter version, these require WebGL with support for floating-point textures.

### Examples ###

* [Tracking via webcam (support for getUserMedia needed)](http://auduno.github.com/clmtrackr/clm_video_gum.html)
* [Tracking in image](http://auduno.github.com/clmtrackr/clm_small.html)
* [Tracking in video](http://auduno.github.com/clmtrackr/clm_video.html)

### License ###

clmtrackr is distributed under the [MIT License](http://www.opensource.org/licenses/MIT)