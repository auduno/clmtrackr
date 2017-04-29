var jsfeat_face = function(video, maxWorkSize) {
  var videoWidth = video.width;
  var videoHeight = video.height;

  var classifier = jsfeat.haar.frontalface;

  // scale down canvas we do detection on (to reduce noisy detections)
  var scale = Math.min(maxWorkSize/videoWidth, maxWorkSize/videoHeight);
  var w = (videoWidth*scale)|0;
  var h = (videoHeight*scale)|0;

  var img_u8 = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
  var edg = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
  var work_canvas = document.createElement('canvas');
  work_canvas.width = w;
  work_canvas.height = h;
  var work_ctx = work_canvas.getContext('2d');
  var ii_sum = new Int32Array((w+1)*(h+1));
  var ii_sqsum = new Int32Array((w+1)*(h+1));
  var ii_tilted = new Int32Array((w+1)*(h+1));
  var ii_canny = new Int32Array((w+1)*(h+1));

  this.findFace = function (params) {
    work_ctx.drawImage(video, 0, 0, work_canvas.width, work_canvas.height);
    var imageData = work_ctx.getImageData(0, 0, work_canvas.width, work_canvas.height);

    jsfeat.imgproc.grayscale(imageData.data, work_canvas.width, work_canvas.height, img_u8);

    // possible params
    if(params.equalizeHistogram) {
      jsfeat.imgproc.equalize_histogram(img_u8, img_u8);
    }
    //jsfeat.imgproc.gaussian_blur(img_u8, img_u8, 3);

    jsfeat.imgproc.compute_integral_image(img_u8, ii_sum, ii_sqsum, classifier.tilted ? ii_tilted : null);

    if(params.useCanny) {
      jsfeat.imgproc.canny(img_u8, edg, 10, 50);
      jsfeat.imgproc.compute_integral_image(edg, ii_canny, null, null);
    }

    jsfeat.haar.edgesDensity = params.edgesDensity;
    var rects = jsfeat.haar.detect_multi_scale(ii_sum, ii_sqsum, ii_tilted, params.useCanny? ii_canny : null, img_u8.cols, img_u8.rows, classifier, params.scaleFactor, params.minScale);
    rects = jsfeat.haar.group_rectangles(rects, 1);

    var rl = rects.length;

    if(rl == 0) {
      return false;
    }

    var best = rects[0];
    for (var i = 1; i < rl; i++) {
      if (rects[i].neighbors > best.neighbors) {
        best = rects[i]
      } else if (rects[i].neighbors == best.neighbors) {
        // if (rects[i].width > best.width) best = rects[i]; // use biggest rect
        if (rects[i].confidence > best.confidence) best = rects[i]; // use most confident rect
      }
    }

    var sc = videoWidth / img_u8.cols;
    best.x = (best.x*sc)|0;
    best.y = (best.y*sc)|0;
    best.width = (best.width*sc)|0;
    best.height = (best.height*sc)|0;

    return best;
  }
}