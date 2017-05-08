//Curtousy of stackoverflow this function
Worker.createURL = function(func_or_string){
    var str = (typeof func_or_string === 'function')?func_or_string.toString():func_or_string;
    var blob = new Blob(['\'use strict\';\nself.onmessage ='+str], { type: 'text/javascript' });
    return window.URL.createObjectURL(blob);
};

Worker.create = function(func_or_string){
  return new Worker(Worker.createURL(func_or_string));
};

// simple wrapper for jsfeat face detector that runs as a webworker
// requires jsfeat
var jsfeat_face = function(video, maxWorkSize) {
  var videoWidth = video.width;
  var videoHeight = video.height;

  // scale down canvas we do detection on (to reduce noisy detections)
  var scale = Math.min(maxWorkSize/videoWidth, maxWorkSize/videoHeight);
  var w = (videoWidth*scale)|0;
  var h = (videoHeight*scale)|0;

  var work_canvas = document.createElement('canvas');
  work_canvas.height = h;
  work_canvas.width = w;
  var work_ctx = work_canvas.getContext('2d');
  // img_u8 = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
  // ii_sum = new Int32Array((w+1)*(h+1));
  // ii_sqsum = new Int32Array((w+1)*(h+1));
  // ii_tilted = new Int32Array((w+1)*(h+1));

  // var classifier = jsfeat.haar.frontalface;

  this.findFace = function (params, callback) {
    work_ctx.drawImage(video, 0, 0, work_canvas.width, work_canvas.height);
    var imageData = work_ctx.getImageData(0, 0, work_canvas.width, work_canvas.height);

    var worker = Worker.create(findFaceWorker);

    worker.addEventListener('message', function (e) {
      this.faceDetected(e, callback);
    }.bind(this), false);

    worker.postMessage({
      w: work_canvas.width,
      h: work_canvas.height,
      videoWidth: videoWidth,
      imageData:imageData,
      params: params
    });

  };

};
