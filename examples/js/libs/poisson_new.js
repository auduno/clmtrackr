/**
 * Poisson Image Editing module
 * http://rest-term.com
 */
(function() {
  var Poisson;       // top-level namaspace
  var _root = this;  // reference to 'window' or 'global'

  if(typeof exports !== 'undefined') {
    Poisson = exports;   // for CommonJS
  } else {
    Poisson = _root.Poisson = {};
  }

  // core operations
  var EPS = 1.0E-08;
  var ctx = null,
      images = [new Image(), new Image(), new Image()],
      files = [],
      data = [],
      loadedCount = 0;

  var sumf = new Float64Array(3),
      sumfstar = new Float64Array(3),
      sumvq = new Float64Array(3),
      fp = new Float64Array(3),
      fq = new Float64Array(3),
      gp = new Float64Array(3),
      gq = new Float64Array(3),
      subf = new Float64Array(3),
      naddr = new Int32Array(4),
      subg = new Float64Array(3);
  function pt1(y, offsetY, x, offsetX, w, h, sumf, blendData, m, sumvq, srcData, l, naddr) {
      if(y + offsetY >= 0 && x + offsetX >= 0 &&
         y + offsetY < h && x + offsetX < w) {
        for(n=0; n<4; n++) {
          for(var c=0; c<3; c++) {
            sumf[c] += blendData[naddr[n] + m + c];
            sumvq[c] += srcData[l + c] - srcData[naddr[n] + c];
          }
        }
      }
  }
  function pt2(y, offsetY, x, offsetX, w, h, fp, gp, dstData, fq, gq, sumfstar, subf, subg, sumf, blendData, m, sumvq, srcData, l, naddr) {
    if(y + offsetY >= 0 && x + offsetX >= 0 &&
       y + offsetY < h && x + offsetX < w) {
      fp[0] = dstData[l + m];
      fp[1] = dstData[l + m + 1];
      fp[2] = dstData[l + m + 2];
      gp[0] = srcData[l];
      gp[1] = srcData[l + 1];
      gp[2] = srcData[l + 2];
      for(n=0; n<4; n++) {
        for(c=0; c<3; c++) {
          fq[c] = dstData[naddr[n] + m + c];
          // modification : we ignore pixels outside face mask, since these cause artifacts
          //gq[c] = srcData[naddr[n] + c];
          gq[c] = srcData[l+c];
          sumfstar[c] += fq[c];
          subf[c] = fp[c] - fq[c];
          subf[c] = subf[c] > 0 ? subf[c] : -subf[c];
          subg[c] = gp[c] - gq[c];
          subg[c] = subg[c] > 0 ? subg[c] : -subg[c];
          if(subf[c] > subg[c]) {
            sumvq[c] += subf[c];
          } else {
            sumvq[c] += subg[c];
          }
        }
      }
    }
  }
  function pt3(fp, c, sumf, sumfstar, sumvq, terminate, EPS, blendData, l, m, srcData, paste) {
    fp[c] = (sumf[c] + sumfstar[c] + sumvq[c])*0.25; // division 4
    error = Math.floor(fp[c] - blendData[l + m + c]);
    error = error > 0 ? error : -error;
    if(terminate[c] && error >
       EPS*(1 + (fp[c] > 0 ? fp[c] : -fp[c]))) {
      terminate[c] = false;
      //log[l+m+c] = {x:x,y:y,c:c,error:error};
    }
    blendData[l + m + c] = fp[c];
    if(paste) blendData[l + m + c] = srcData[l + c];
  }
  
  var core = {
    // loads images
    load: function(srcCanvas, dstCanvas, maskCanvas, onComplete, onError) {
      ctx = document.createElement('canvas').getContext('2d');
      // load complete handler
      
      var imgData = srcCanvas.getContext('2d').getImageData(0, 0, srcCanvas.width, srcCanvas.height);
      data[0] = imgData;
      imgData = dstCanvas.getContext('2d').getImageData(0, 0, dstCanvas.width, dstCanvas.height);
      data[1] = imgData;
      imgData = maskCanvas.getContext('2d').getImageData(0, 0, maskCanvas.width, maskCanvas.height);
      data[2] = imgData;
      ctx.drawImage(dstCanvas, 0, 0);
      data[3] = dstCanvas.getContext('2d').getImageData(0, 0, dstCanvas.width, dstCanvas.height);
      if(typeof onComplete === 'function') {
        onComplete(data);
      }

    },
    reset: function() {
      ctx.drawImage(dstCanvas, 0, 0);
      data[3] = ctx.getImageData(0, 0, dstCanvas.width, dstCanvas.height);
      return data[3];
    },
    // applies poisson image editing
    blend: function(iteration, offsetX, offsetY, paste) {
      var w = data[0].width,
          h = data[0].height,
          len = w*h*4,
          srcData = data[0].data,
          dstData = data[1].data,
          maskData = data[2].data,
          blendData = data[3].data,
          edge = false,
          error = 0.0,
          threashold = 128,
          terminate = [],
          step, l, m;
      var log = {};
      // validation
      if(!(parseInt(iteration) && typeof(offsetX) == "number" && typeof(offsetY) == "number")) {
        throw TypeError('invalid parameter type');
      }
      // core operation
      for(var i=0; i<iteration; i++) {
        terminate = [true, true, true];
        for(var y=1; y<h-1; y++) {
          step = y*w << 2;
          for(var x=1; x<w-1; x++) {
            l = step + (x << 2);
            m = offsetY*w + offsetX << 2;
            naddr[0] = l - (w << 2);
            naddr[1] = l - 4;
            naddr[2] = l + 4;
            naddr[3] = l + (w << 2);
            if(maskData[l] > threashold) { // on the mask
              sumf[0] = 0;
              sumf[1] = 0;
              sumf[2] = 0;
              sumfstar[0] = 0;
              sumfstar[1] = 0;
              sumfstar[2] = 0;
              sumvq[0] = 0;
              sumvq[1] = 0;
              sumvq[2] = 0;
              edge = false;
              for(var n=0; n<4; n++) {
                if(maskData[naddr[n]] <= threashold) {
                  edge = true;
                  break;
                }
              }
              if(!edge) {
                pt1(y, offsetY, x, offsetX, w, h, sumf, blendData, m, sumvq, srcData, l, naddr);
              } else {
                pt2(y, offsetY, x, offsetX, w, h, fp, gp, dstData, fq, gq, sumfstar, subf, subg, sumf, blendData, m, sumvq, srcData, l, naddr);
              }
              for(c=0; c<3; c++) {
                pt3(fp, c, sumf, sumfstar, sumvq, terminate, EPS, blendData, l, m, srcData, paste)
              }
            } // end mask
          } // end x loop
        } // end y loop
        if(terminate[0] && terminate[1] && terminate[2]) break;
      } // end iteration
      //core.debugPlot(blendData, log);
      return data[3];
    },
    debugPlot: function(img, data) {
      for(var i in data) {
        console.log(data[i]);
        img[i] = 255;
      }
    }
  };
  // aliases (public APIs)
  Poisson.load = core.load;
  Poisson.reset = core.reset;
  Poisson.blend = core.blend;
}).call(this);