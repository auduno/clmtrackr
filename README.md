clmtrackr
======

![tracked face](http://auduno.github.com/clmtrackr/media/clmtrackr01.jpg)

clmtrackr is a work-in-progress javascript library for precise tracking of facial features in videos or images. It currently is an implementation of *constrained local models* fitted by *regularized landmark mean-shift*, as described in [Jason M. Saragih's paper](http://dl.acm.org/citation.cfm?id=1938021). Due to the heavy calculations involved, the current version of the library requires WebGL with support for floating-point textures.

[Xiaoguang Yan](https://sites.google.com/site/xgyanhome/my-information) kindly allowed us to use the trained model from [his implementation of Constrained Local Models](https://sites.google.com/site/xgyanhome/home/projects/clm-implementation). The aim is to provide a more generic model that will fit general faces.

The library requires [ccv](https://github.com/liuliu/ccv) (for initial face detection) and parts of [google closure library](https://developers.google.com/closure/library/) (for matrix math).

For some more information about Constrained Local Models, take a look at Xiaoguang Yan's [excellent tutorial](https://sites.google.com/site/xgyanhome/home/projects/clm-implementation/ConstrainedLocalModel-tutorial%2Cv0.7.pdf?attredirects=0), which was of great help in implementing this library.

### Examples ###

* [Tracking in image](http://auduno.github.com/clmtrackr/clm_small.html)
* [Tracking in video](http://auduno.github.com/clmtrackr/clm_video.html)

Note that currently the code does not contain any mechanism for detecting and reinitializing when the tracking fails, so the tracking will fail (badly) at around 0:25 in the video.

### License ###

clmtrackr is distributed under the [MIT License](http://www.opensource.org/licenses/MIT)