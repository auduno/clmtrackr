// simple wrapper for jsfeat face detector
// requires jsfeat

var jsfeat_face = function(image) {
  
  var img_u8,work_canvas,work_ctx,ii_sum,ii_sqsum,ii_tilted,edg;
  
  var w = image.width;
  var h = image.height;
  
  if (image.tagName == 'VIDEO' || image.tagName == 'IMG') {
    work_canvas = document.createElement('canvas');
    work_canvas.height = h;
    work_canvas.width = w;
    work_ctx = work_canvas.getContext('2d');
  } else if (image.tagName == 'CANVAS') {
    work_ctx = image.getContext('2d');
  }
  
  img_u8 = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
  ii_sum = new Int32Array((w+1)*(h+1));
  ii_sqsum = new Int32Array((w+1)*(h+1));
  ii_tilted = new Int32Array((w+1)*(h+1));
  
  var classifier = jsfeat.haar.frontalface;
  
  this.findFace = function () {
    if (image.tagName == 'VIDEO' || image.tagName == 'IMG') {
      work_ctx.drawImage(image, 0, 0);
    } 
    var imageData = work_ctx.getImageData(0, 0, w, h);
                  
    jsfeat.imgproc.grayscale(imageData.data, img_u8.data);
    
    jsfeat.imgproc.equalize_histogram(img_u8, img_u8);
    
    jsfeat.imgproc.compute_integral_image(img_u8, ii_sum, ii_sqsum, null);

    var rects = jsfeat.haar.detect_multi_scale(ii_sum, ii_sqsum, ii_tilted, null, img_u8.cols, img_u8.rows, classifier, 1.15, 2);
    
    rects = jsfeat.haar.group_rectangles(rects, 1);
    
    var rl = rects.length;
    
    if (rl > 0) {
      var best = rects[0];
      for (var i = 1;i < rl;i++) {
          if (rects[i].neighbors > best.neighbors) {
              best = rects[i]
          } else if (rects[i].neighbors == best.neighbors) {
              if (rects[i].confidence > best.confidence) best = rects[i];
          }
      }
      return [best];
    } else {
      return false;
    }
  }
  
}