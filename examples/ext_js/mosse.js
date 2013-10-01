"use strict";
/*
 * MOSSE correlation filter
 *
 * Optional parameters to constructor:
 *   drawResponse {canvasElement} : draws the correlation filter output on the given canvas element (default is none)
 *   psrThreshold {number} : peak-to-sidelobe-ratio threshold to use when updating filter while tracking (default is 10)
 *   eta {number} : adjusts how much new input affects the mosse filter, when updating filter while tracking
 *     number should be between 0 and 1 (default is 0.1)
 *   convertToGrayscale {boolean} : whether to convert canvas output to grayscale (default is true)
 *     if this is set to false, we assume all channels are equal and only grab values from red channel
 *
 * @author auduno / github.com/auduno
 */ 

function mosseFilter(params) {
    
    var _filter, _top, _bottom;
    var _fft;
    var _w,_h;
    var _im_part;
    var _arrlen;
    var _cc;
    var _image_array;
    
    this.psr_prev = undefined;
    this.peak_prev = undefined;
    var peak = 0.0;
    var updateable = false;
    
    if (!params) params = {};
    // setup of canvas for drawing responses, if given
    if (params.drawResponse === undefined) {
        params.drawResponse = false;
    } else {
        if (params.drawResponse.tagName != 'CANVAS') {
            params.drawResponse = false;
        } else {
            var responseContext = params.drawResponse.getContext('2d');
        }
    }
    if (params.psrThreshold === undefined) params.psrThreshold = 10;
    if (params.eta === undefined) params.eta = 0.10;
    if (params.convertToGrayscale === undefined) params.convertToGrayscale = true;
    
    this.load = function(filter) {
        // initialize filter width and height
        _w = filter.width;
        _h = filter.height;
        _arrlen = _w*_h;
        _filter = [filter.real, filter.imag];
        // handling top and bottom when they're not present
        if (filter.top && filter.bottom) {
          updateable = true;
          _top = [filter.top.real, filter.top.imag];
          _bottom = [filter.bottom.real, filter.bottom.imag];
        }
        
        // initialize fft to given width
        _fft = new FFT();
        _fft.init(filter.width);
        
        // set up temporary variables
        if(typeof Float64Array !== 'undefined') {
            _im_part = new Float64Array(_arrlen);
            _image_array = new Float64Array(_arrlen);
        } else {
            _im_part = new Array(_arrlen);
            _image_array = new Array(_arrlen);
        }
        var canvas = document.createElement("canvas");
        canvas.setAttribute('width', _w);
        canvas.setAttribute('height', _h);
        _cc = canvas.getContext('2d');
    }
    
    this.init = function(w,h) {
        // initialize filter width and height for a blank filter
        _w = w;
        _h = h;
        _arrlen = _w*_h;
        
        _filter = [[],[]];
        _top = [[],[]];
        _bottom = [[],[]];
        for (var i = 0;i < _arrlen;i++) {
            _filter[0][i] = 0;
            _filter[1][i] = 0;
            _top[0][i] = 0;
            _top[1][i] = 0;
            _bottom[0][i] = 0;
            _bottom[1][i] = 0;
        }
        updateable = true;
        
        // initialize fft to given width
        _fft = new FFT();
        _fft.init(w);
        
        // set up temporary variables
        if(typeof Float64Array !== 'undefined') {
            _im_part = new Float64Array(_arrlen);
        } else {
            _im_part = new Array(_arrlen);
        }
        var canvas = document.createElement("canvas");
        canvas.setAttribute('width', _w);
        canvas.setAttribute('height', _h);
        _cc = canvas.getContext('2d');
    }
    
    // fft function
    this.fft = function(array) {
        // not in-place
        
        var cn = new Array(_arrlen);
        for (var i = 0;i < _arrlen;i++) {
          cn[i] = 0.0;
        }
        
        _fft.fft2d(array,cn)
        return [array, cn];
    }
    
    // fft function
    this.fft_inplace = function(array) {
        // in-place
        
        for (var i = 0;i < _arrlen;i++) {
          _im_part[i] = 0.0;
        }
        
        _fft.fft2d(array,_im_part)
        return [array, _im_part];
    }
    
    this.ifft = function(rn, cn) {
        // in-place
        _fft.ifft2d(rn, cn);
        return rn;
    }

    // peak to sidelobe ratio function (optional)
    this.psr = function(array) {
        // proper
        var sum = 0;
        var max = 0;
        var maxpos = [];
        var sdo = 0;
        var val;
        for (var x = 0;x < _w;x++) {
            for (var y = 0;y < _h;y++) {
                val = array[(y*_w)+x];
                sum += val;
                sdo += (val*val);
                if (max < val) {
                    max = val;
                    maxpos = [x,y];
                }
            }
        }
        
        // subtract values around peak
        for (var x = -5;x < 6;x++) {
            for (var y = -5;y < 6;y++) {
                if (Math.sqrt(x*x+y*y) < 5) {
                    val = array[((maxpos[1]+y)*_w)+(maxpos[0]+x)]
                    sdo -= (val*val);
                    sum -= val;
                }
            }
        }
        
        var mean = sum/array.length;
        var sd = Math.sqrt((sdo/array.length)-(mean*mean));
        
        // get mean/variance of output around peak
        var psr = (max-mean)/sd;
        return psr;
    }
    
    this.getResponse = function(imageData) {
        // in-place
        
        // preprocess
        var prepImage = preprocess(imageData);
        prepImage = cosine_window(prepImage);
        
        // filter
        var res = this.fft_inplace(prepImage);
        
        // elementwise multiplication with filter
        complex_mult_inplace(res, _filter);
        
        // do inverse 2d fft
        var filtered = this.ifft(res[0],res[1]);
        return filtered;
    }
    
    this.track = function(input, left, top, width, height, updateFilter, gaussianPrior, calcPSR) {
        // finds position of filter in input image
        
        if (!_filter) {
            console.log("Mosse-filter needs to be initialized or trained before starting tracking.");
            return false;
        }
        
        if (input.tagName == "VIDEO" || input.tagName == "IMG") {
            // scale selection according to original source image
            var videoLeft = Math.round((left/input.width)*input.videoWidth);
            var videoTop = Math.round((top/input.height)*input.videoHeight);
            var videoWidth = Math.round((width/input.width)*input.videoWidth);
            var videoHeight = Math.round((height/input.height)*input.videoHeight);
            _cc.drawImage(input, videoLeft, videoTop, videoWidth, videoHeight, 0, 0, _w, _h);
        } else if (input.tagName == "CANVAS") {
            _cc.drawImage(input, left, top, width, height, 0, 0, _w, _h);
        }
        
        var image = _cc.getImageData(0,0,_w,_h);
        var id = image.data;
        
        if (params.convertToGrayscale) {
            // convert to grayscale
            for (var i = 0;i < _arrlen;i++) {
                _image_array[i] = id[(4*i)]*0.3;
                _image_array[i] += id[(4*i)+1]*0.59;
                _image_array[i] += id[(4*i)+2]*0.11;
            } 
        } else {
            // use only one channel
            for (var i = 0;i < _arrlen;i++) {
                _image_array[i] = id[(4*i)];
            } 
        }
        
        // preprocess
        var prepImage = preprocess(_image_array);
        prepImage = cosine_window(prepImage);
        
        // filter
        var res = this.fft_inplace(prepImage);
        // elementwise multiplication with filter
        var nures = complex_mult(res, _filter);
        // do inverse 2d fft
        var filtered = this.ifft(nures[0],nures[1]);
        
        // find max and min
        var max = 0;
        var min = 0;
        var maxpos = [];
        
        //method using centered gaussian prior
        if (gaussianPrior) {
            var prior, dx, dy;
            var variance = 128;
            for (var x = 0;x < _w;x++) {
                for (var y = 0;y < _h;y++) {
                    dx = x - _w/2;
                    dy = y - _h/2;
                    prior = Math.exp(-0.5*((dx*dx)+(dy*dy))/variance)
                    if ((filtered[(y*_w)+x]*prior) > max) {
                        max = filtered[(y*_w)+x]*prior;
                        maxpos = [x,y];
                    }
                    if (filtered[(y*_w)+x] < min) {
                        min = filtered[(y*_w)+x];
                    }
                }
            }
        } else {
            for (var x = 0;x < _w;x++) {
                for (var y = 0;y < _h;y++) {
                    if (filtered[(y*_w)+x] > max) {
                        max = filtered[(y*_w)+x];
                        maxpos = [x,y];
                    }
                    if (filtered[(y*_w)+x] < min) {
                        min = filtered[(y*_w)+x];
                    }
                }
            }
        }
        this.peak_prev = max;
        
        if (params.drawResponse) {
            // draw response
            var diff = max-min;
            var dc = document.createElement('canvas');
            dc.setAttribute('width', 32);
            dc.setAttribute('height', 32);
            var dcc = dc.getContext('2d');
            var psci = dcc.createImageData(32, 32);
            var pscidata = psci.data;
            for (var j = 0;j < 32*32;j++) {
                //draw with priors
                //var val = filtered[j]*Math.exp(-0.5*(((j%_w - _w/2)*(j%_w -_w/2))+((Math.floor(j/_h)-(_h/2))*(Math.floor(j/_h)-(_h/2))))/128);
                var val = filtered[j];
                val = Math.round((val+Math.abs(min))*(255/diff));
                pscidata[j*4] = val;
                pscidata[(j*4)+1] = val;
                pscidata[(j*4)+2] = val;
                pscidata[(j*4)+3] = 255;
            }
            dcc.putImageData(psci, 0, 0);
            responseContext.drawImage(dc, left, top, width, width);
        }
        
        if (calcPSR) {
          this.psr_prev = this.psr(filtered);
        }
        
        if (updateFilter) {
            if (!updateable) {
                console.log("The loaded filter does not support updating. Ignoring parameter 'updateFilter'.");
            } else {
                if (calcPSR) {
                  var psr = this.psr_prev;
                } else {
                  var psr = this.psr(filtered);
                }
                
                if (psr > params.psrThreshold) {
                    // create target
                    var target = [];
                    var nux = maxpos[0];
                    var nuy = maxpos[1];
                    for (var x = 0;x < _w;x++) {
                        for (var y = 0;y < _h;y++) {
                            target[(y*_w)+x] = Math.exp(-(((x-nux)*(x-nux))+((y-nuy)*(y-nuy)))/(2*2));
                        }
                    }
                    
                    //fft target
                    target = this.fft(target);
                    
                    // create filter
                    var res_conj = complex_conj(res);
                    var fuTop = complex_mult(target,res_conj);
                    var fuBottom = complex_mult(res,res_conj);
                    
                    // add up
                    var eta = params.eta;
                    for (var i = 0;i < _arrlen;i++) {
                        _top[0][i] = eta*fuTop[0][i] + (1-eta)*_top[0][i];
                        _top[1][i] = eta*fuTop[1][i] + (1-eta)*_top[1][i];
                        _bottom[0][i] = eta*fuBottom[0][i] + (1-eta)*_bottom[0][i];
                        _bottom[1][i] = eta*fuBottom[1][i] + (1-eta)*_bottom[1][i];
                    }
                    
                    _filter = complex_div(_top,_bottom);
                }
            }
        }
        
        /*if (psr < 5) {
          maxpos = [_w/2,_h/2]; 
        }*/
        
        maxpos[0] = maxpos[0]*(width/_w);
        maxpos[1] = maxpos[1]*(width/_h);
        
        // check if output is strong enough
        // if not, return false?
        if (max < 0) {
          return false;
        } else {
          return maxpos;
        }
    }
    
    this.train = function(input, left, top, width, height) {
        
        if (!updateable) {
          console.log("The loaded filter does not support updating. Unable to do training.");
          return false;
        }
        
        if (input.tagName == "VIDEO" || input.tagName == "IMG") {
            // scale selection according to original source image
            var videoLeft = Math.round((left/input.width)*input.videoWidth);
            var videoTop = Math.round((top/input.height)*input.videoHeight);
            var videoWidth = Math.round((width/input.width)*input.videoWidth);
            var videoHeight = Math.round((height/input.height)*input.videoHeight);
            _cc.drawImage(input, videoLeft, videoTop, videoWidth, videoHeight, 0, 0, _w, _h);
        } else if (input.tagName == "CANVAS") {
            _cc.drawImage(input, left, top, width, height, 0, 0, _w, _h);
        }
        
        var image = _cc.getImageData(0,0,_w,_h);
        var id = image.data;
         
        // convert to grayscale
        for (var i = 0;i < _arrlen;i++) {
            _image_array[i] = id[(4*i)]*0.3;
            _image_array[i] += id[(4*i)+1]*0.59;
            _image_array[i] += id[(4*i)+2]*0.11;
        }
        
        // preprocess
        var prepImage = preprocess(_image_array);
        prepImage = cosine_window(prepImage);
        
        // create target
        var target = [];
        var nux = _w/2;
        var nuy = _h/2;
        for (var x = 0;x < _w;x++) {
            for (var y = 0;y < _h;y++) {
                target[(y*_w)+x] = Math.exp(-(((x-nux)*(x-nux))+((y-nuy)*(y-nuy)))/(2*2));
            }
        }
        
        //fft target
        target = this.fft(target);
        
        // filter
        var res = this.fft(prepImage);
        // create filter
        var res_conj = complex_conj(res);
        var fuTop = complex_mult(target,res_conj);
        var fuBottom = complex_mult(res,res_conj);
        
        // add up
        var eta = params.eta;
        for (var i = 0;i < _arrlen;i++) {
            _top[0][i] = eta*fuTop[0][i] + (1-eta)*_top[0][i];
            _top[1][i] = eta*fuTop[1][i] + (1-eta)*_top[1][i];
            _bottom[0][i] = eta*fuBottom[0][i] + (1-eta)*_bottom[0][i];
            _bottom[1][i] = eta*fuBottom[1][i] + (1-eta)*_bottom[1][i];
        }
        
        _filter = complex_div(_top,_bottom);
        
        return true;
    }
    
    var preprocess = function(array) {
        // in-place
        
        // log adjusting
        for (var i = 0;i < _arrlen;i++) {
          array[i] = Math.log(array[i]+1);
        }
        
        // normalize to mean 0 and norm 1
        var mean = 0;
        for (var i = 0;i < _arrlen;i++) {
          mean += array[i];
        }
        mean /= _arrlen;
        
        for (var i = 0;i < _arrlen;i++) {
          array[i] -= mean;
        }
        var norm = 0.0;
        for (var i = 0;i < _arrlen;i++) {
          norm += (array[i]*array[i]);
        }
        norm = Math.sqrt(norm);
        for (var i = 0;i < _arrlen;i++) {
          array[i] /= norm;
        }
        
        return array;
    }
    
    var cosine_window = function(array) {
        // calculate rect cosine window (in-place)
        var pos = 0;
        for (var i = 0;i < _w;i++) {
            for (var j = 0;j < _h;j++) {
                //pos = (i%_w)+(j*_w);
                var cww = Math.sin((Math.PI*i)/(_w-1))
                var cwh = Math.sin((Math.PI*j)/(_h-1))
                array[pos] = Math.min(cww,cwh)*array[pos];
                pos++;
            }
        }
        
        return array;
    }
    
    var complex_mult = function(cn1, cn2) {
        // not in-place
        var re_part = new Array(_w);
        var im_part = new Array(_w);
        var nucn = [re_part, im_part];
        for (var r = 0;r < _arrlen;r++) {
            nucn[0][r] = (cn1[0][r]*cn2[0][r]) - (cn1[1][r]*cn2[1][r]);
            nucn[1][r] = (cn1[0][r]*cn2[1][r]) + (cn1[1][r]*cn2[0][r]);
        }
        return nucn;
    }
    
    var complex_mult_inplace = function(cn1, cn2) {
        // in-place
        var temp1, temp2;
        for (var r = 0;r < _arrlen;r++) {
            temp1 = (cn1[0][r]*cn2[0][r]) - (cn1[1][r]*cn2[1][r]);
            temp2 = (cn1[0][r]*cn2[1][r]) + (cn1[1][r]*cn2[0][r]);
            cn1[0][r] = temp1;
            cn1[1][r] = temp2;
        }
    }
    
    var complex_conj = function(cn) {
        // not in-place (TODO)
        var nucn = [[],[]];
        for (var i = 0;i < _arrlen;i++) {
            nucn[0][i] = cn[0][i]
            nucn[1][i] = -cn[1][i];
        }
        return nucn;
    }
    
    var complex_div = function(cn1, cn2) {
        // not in-place (TODO)
        var nucn = [[],[]];
        for (var r = 0;r < _arrlen;r++) {
            nucn[0][r] = ((cn1[0][r]*cn2[0][r])+(cn1[1][r]*cn2[1][r])) / ((cn2[0][r]*cn2[0][r]) + (cn2[1][r]*cn2[1][r]));
            nucn[1][r] = ((cn1[1][r]*cn2[0][r])-(cn1[0][r]*cn2[1][r])) / ((cn2[0][r]*cn2[0][r]) + (cn2[1][r]*cn2[1][r]));
        }
        return nucn;
    }
}

/**
 * Fast Fourier Transform
 * 1D-FFT/IFFT, 2D-FFT/IFFT (radix-2)
 * 
 * @author ryo / github.com/wellflat
 * Based on https://github.com/wellflat/jslib with some tiny optimizations
 */

function FFT() {
  
  var _n = 0,          // order
      _bitrev = null,  // bit reversal table
      _cstb = null;    // sin/cos table
  var _tre, _tim;
  
  this.init = function (n) {
    if(n !== 0 && (n & (n - 1)) === 0) {
      _n = n;
      _setVariables();
      _makeBitReversal();
      _makeCosSinTable();
    } else {
      throw new Error("init: radix-2 required");
    }
  }
    
  // 1D-FFT
  this.fft1d = function (re, im) {
    fft(re, im, 1);
  }
    
  // 1D-IFFT
  this.ifft1d = function (re, im) {
    var n = 1/_n;
    fft(re, im, -1);
    for(var i=0; i<_n; i++) {
      re[i] *= n;
      im[i] *= n;
    }
  }
  
  // 2D-FFT
  this.fft2d = function (re, im) {
    var i = 0;
    // x-axis
    for(var y=0; y<_n; y++) {
      i = y*_n;
      for(var x1=0; x1<_n; x1++) {
        _tre[x1] = re[x1 + i];
        _tim[x1] = im[x1 + i];
      }
      this.fft1d(_tre, _tim);
      for(var x2=0; x2<_n; x2++) {
        re[x2 + i] = _tre[x2];
        im[x2 + i] = _tim[x2];
      }
    }
    // y-axis
    for(var x=0; x<_n; x++) {
      for(var y1=0; y1<_n; y1++) {
        i = x + y1*_n;
        _tre[y1] = re[i];
        _tim[y1] = im[i];
      }
      this.fft1d(_tre, _tim);
      for(var y2=0; y2<_n; y2++) {
        i = x + y2*_n;
        re[i] = _tre[y2];
        im[i] = _tim[y2];
      }
    }
  }
  
  // 2D-IFFT
  this.ifft2d = function (re, im) {
    var i = 0;
    // x-axis
    for(var y=0; y<_n; y++) {
      i = y*_n;
      for(var x1=0; x1<_n; x1++) {
        _tre[x1] = re[x1 + i];
        _tim[x1] = im[x1 + i];
      }
      this.ifft1d(_tre, _tim);
      for(var x2=0; x2<_n; x2++) {
        re[x2 + i] = _tre[x2];
        im[x2 + i] = _tim[x2];
      }
    }
    // y-axis
    for(var x=0; x<_n; x++) {
      for(var y1=0; y1<_n; y1++) {
        i = x + y1*_n;
        _tre[y1] = re[i];
        _tim[y1] = im[i];
      }
      this.ifft1d(_tre, _tim);
      for(var y2=0; y2<_n; y2++) {
        i = x + y2*_n;
        re[i] = _tre[y2];
        im[i] = _tim[y2];
      }
    }
  }
  
  // core operation of FFT
  function fft(re, im, inv) {
    var d, h, ik, m, tmp, wr, wi, xr, xi,
        n4 = _n >> 2;
    // bit reversal
    for(var l=0; l<_n; l++) {
      m = _bitrev[l];
      if(l < m) {
        tmp = re[l];
        re[l] = re[m];
        re[m] = tmp;
        tmp = im[l];
        im[l] = im[m];
        im[m] = tmp;
      }
    }
    // butterfly operation
    for(var k=1; k<_n; k<<=1) {
      h = 0;
      d = _n/(k << 1);
      for(var j=0; j<k; j++) {
        wr = _cstb[h + n4];
        wi = inv*_cstb[h];
        for(var i=j; i<_n; i+=(k<<1)) {
          ik = i + k;
          xr = wr*re[ik] + wi*im[ik];
          xi = wr*im[ik] - wi*re[ik];
          re[ik] = re[i] - xr;
          re[i] += xr;
          im[ik] = im[i] - xi;
          im[i] += xi;
        }
        h += d;
      }
    }
  }
  
  // set variables
  function _setVariables() {
    if(typeof Uint8Array !== 'undefined') {
      _bitrev = new Uint8Array(_n);
    } else {
      _bitrev = new Array(_n);
    }
    if(typeof Float64Array !== 'undefined') {
      _cstb = new Float64Array(_n*1.25);
      _tre = new Float64Array(_n*_n);
      _tim = new Float64Array(_n*_n);
    } else {
      _cstb = new Array(_n*1.25);
      _tre = new Array(_n*_n);
      _tim = new Array(_n*_n);
    }
  }
  
  // make bit reversal table
  function _makeBitReversal() {
    var i = 0,
        j = 0,
        k = 0;
    _bitrev[0] = 0;
    while(++i < _n) {
      k = _n >> 1;
      while(k <= j) {
        j -= k;
        k >>= 1;
      }
      j += k;
      _bitrev[i] = j;
    }
  }
  
  // make trigonometric function table
  function _makeCosSinTable() {
    var n2 = _n >> 1,
        n4 = _n >> 2,
        n8 = _n >> 3,
        n2p4 = n2 + n4,
        t = Math.sin(Math.PI/_n),
        dc = 2*t*t,
        ds = Math.sqrt(dc*(2 - dc)),
        c = _cstb[n4] = 1,
        s = _cstb[0] = 0;
    t = 2*dc;
    for(var i=1; i<n8; i++) {
      c -= dc;
      dc += t*c;
      s += ds;
      ds -= t*s;
      _cstb[i] = s;
      _cstb[n4 - i] = c;
    }
    if(n8 !== 0) {
      _cstb[n8] = Math.sqrt(0.5);
    }
    for(var j=0; j<n4; j++) {
      _cstb[n2 - j]  = _cstb[j];
    }
    for(var k=0; k<n2p4; k++) {
      _cstb[k + n2] = -_cstb[k];
    }
  }
}
