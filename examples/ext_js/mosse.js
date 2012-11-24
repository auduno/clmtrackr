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
    
    var preprocess = function(array) {
        // log adjusting
        var processed = array.map(function(a) {return Math.log(a+1);});
        // normalize to mean 0 and norm 1
        var mean = 0;
        for (var i = 0;i < processed.length;i++) {
          mean += processed[i];
        }
        mean /= processed.length;
        processed = processed.map(function(x) {return (x-mean);});
        norm = processed.reduce(function(a,b) {return a+(b*b)});
        norm = Math.sqrt(norm);
        processed = processed.map(function(x) {return x/norm;});
        
        return processed;
    }
    
    var cosine_window = function(array) {
        // calculate rect cosine window
        var r,x,y,pos;
        var newArray = [];
        for (var i = 0;i < _w;i++) {
            for (var j = 0;j < _h;j++) {
                x = i-(_w/2);
                y = j-(_h/2);
                pos = (i%_w)+(j*_w);
                var cww = Math.sin((Math.PI*i)/(_w-1))
                var cwh = Math.sin((Math.PI*j)/(_h-1))
                newArray[pos] = Math.min(cww,cwh)*array[pos];
            }
        }
        
        return newArray;
    }
    
    var complex_mult = function(cn1, cn2) {
        var nucn = [[],[]];
        var cnlen = cn1[0].length;
        for (var r = 0;r < cnlen;r++) {
            nucn[0][r] = (cn1[0][r]*cn2[0][r]) - (cn1[1][r]*cn2[1][r]);
            nucn[1][r] = (cn1[0][r]*cn2[1][r]) + (cn1[1][r]*cn2[0][r]);
        }
        return nucn;
    }
    
    var complex_conj = function(cn) {
        var nucn = [[],[]];
        cnlen = cn[1].length;
        for (var i = 0;i < cnlen;i++) {
            nucn[0][i] = cn[0][i]
            nucn[1][i] = -cn[1][i];
        }
        return nucn;
    }
    
    var complex_div = function(cn1, cn2) {
        var nucn = [[],[]];
        var cnlen = cn1[0].length;
        for (var r = 0;r < cnlen;r++) {
            nucn[0][r] = ((cn1[0][r]*cn2[0][r])+(cn1[1][r]*cn2[1][r])) / ((cn2[0][r]*cn2[0][r]) + (cn2[1][r]*cn2[1][r]));
            nucn[1][r] = ((cn1[1][r]*cn2[0][r])-(cn1[0][r]*cn2[1][r])) / ((cn2[0][r]*cn2[0][r]) + (cn2[1][r]*cn2[1][r]));
        }
        return nucn;
    }
    
    this.load = function(filter) {
        // initialize filter width and height
        _w = filter.width;
        _h = filter.height;
        _filter = [filter.real, filter.imag];
        _top = [filter.top.real, filter.top.imag];
        _bottom = [filter.bottom.real, filter.bottom.imag];
        
        // initialize fft to given width
        _fft = FFT;
        _fft.init(filter.width);
    }
    
    this.init = function(w,h) {
        // initialize filter width and height
        _w = w;
        _h = h;
        
        _filter = [[],[]];
        _top = [[],[]];
        _bottom = [[],[]];
        for (var i = 0;i < _w*_h;i++) {
            _filter[0][i] = 0;
            _filter[1][i] = 0;
            _top[0][i] = 0;
            _top[1][i] = 0;
            _bottom[0][i] = 0;
            _bottom[1][i] = 0;
        }
        
        // initialize fft to given width
        _fft = FFT;
        _fft.init(w);
    }
    
    // fft function
    this.fft = function(array) {
        var cn = [];
        var rn = [];
        
        for (var i = 0;i < array.length;i++) {
          cn[i] = 0.0;
          rn[i] = array[i];
        }
        _fft.fft2d(rn,cn)
        return [rn, cn];
    }
    
    this.ifft = function(rn, cn) {
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
    
    this.track = function(input, left, top, width, height, updateFilter, gaussianPrior) {
        // finds position of filter in input image
        
        if (!_filter) {
            console.log("Mosse-filter needs to be initialized or trained before starting tracking.");
            return false;
        }
        
        var canvas = document.createElement("canvas");
        canvas.setAttribute('width', _w);
        canvas.setAttribute('height', _h);
        var cc = canvas.getContext('2d');
        
        if (input.tagName == "VIDEO" || input.tagName == "IMG") {
            // scale selection according to original source image
            var videoLeft = Math.round((left/input.width)*input.videoWidth);
            var videoTop = Math.round((top/input.height)*input.videoHeight);
            var videoWidth = Math.round((width/input.width)*input.videoWidth);
            var videoHeight = Math.round((height/input.height)*input.videoHeight);
            cc.drawImage(input, videoLeft, videoTop, videoWidth, videoHeight, 0, 0, _w, _h);
        } else if (input.tagName == "CANVAS") {
            cc.drawImage(input, left, top, width, height, 0, 0, _w, _h);
        }
        
        var image = cc.getImageData(0,0,_w,_h);
        var id = image.data;
        
        var imagesize = _w*_h;
        var newImage = [];  
        if (params.convertToGrayscale) {
            // convert to grayscale
            for (var i = 0;i < imagesize;i++) {
                newImage[i] = id[(4*i)]*0.3;
                newImage[i] += id[(4*i)+1]*0.59;
                newImage[i] += id[(4*i)+2]*0.11;
            } 
        } else {
            // use only one channel
            for (var i = 0;i < imagesize;i++) {
                newImage[i] = id[(4*i)];
            } 
        }
        
        // preprocess
        var prepImage = preprocess(newImage);
        prepImage = cosine_window(prepImage);
        
        // filter
        var res = this.fft(prepImage);
        // elementwise multiplication with filter
        var nures = complex_mult(res, _filter);
        // do inverse 2d fft
        var filtered = this.ifft(nures[0],nures[1]);
        
        // find max and min
        var max = 0;
        var min = 0;
        maxpos = [];
        
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
        
        
        if (updateFilter) {
            var psr = this.psr(filtered);
            
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
                var fulen = fuTop[0].length;
                for (var i = 0;i < fulen;i++) {
                    _top[0][i] = eta*fuTop[0][i] + (1-eta)*_top[0][i];
                    _top[1][i] = eta*fuTop[1][i] + (1-eta)*_top[1][i];
                    _bottom[0][i] = eta*fuBottom[0][i] + (1-eta)*_bottom[0][i];
                    _bottom[1][i] = eta*fuBottom[1][i] + (1-eta)*_bottom[1][i];
                }
                
                _filter = complex_div(_top,_bottom);
            }
        }
        
        /*if (psr < 5) {
          maxpos = [_w/2,_h/2]; 
        }*/
        
        //maxpos[0] = Math.round(maxpos[0]*(width/_w));
        //maxpos[1] = Math.round(maxpos[1]*(width/_h));
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
        
        var canvas = document.createElement("canvas");
        canvas.setAttribute('width', _w);
        canvas.setAttribute('height', _h);
        var cc = canvas.getContext('2d');
        
        if (input.tagName == "VIDEO" || input.tagName == "IMG") {
            // scale selection according to original source image
            var videoLeft = Math.round((left/input.width)*input.videoWidth);
            var videoTop = Math.round((top/input.height)*input.videoHeight);
            var videoWidth = Math.round((width/input.width)*input.videoWidth);
            var videoHeight = Math.round((height/input.height)*input.videoHeight);
            cc.drawImage(input, videoLeft, videoTop, videoWidth, videoHeight, 0, 0, _w, _h);
        } else if (input.tagName == "CANVAS") {
            cc.drawImage(input, left, top, width, height, 0, 0, _w, _h);
        }
        
        var image = cc.getImageData(0,0,_w,_h);
        var id = image.data;
        
        // convert to grayscale
        var imagesize = _w*_h;
        var newImage = [];  
        for (var i = 0;i < imagesize;i++) {
            newImage[i] = id[(4*i)]*0.3;
            newImage[i] += id[(4*i)+1]*0.59;
            newImage[i] += id[(4*i)+2]*0.11;
        }
        
        // preprocess
        var prepImage = preprocess(newImage);
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
        var fulen = fuTop[0].length;
        for (var i = 0;i < fulen;i++) {
            _top[0][i] = eta*fuTop[0][i] + (1-eta)*_top[0][i];
            _top[1][i] = eta*fuTop[1][i] + (1-eta)*_top[1][i];
            _bottom[0][i] = eta*fuBottom[0][i] + (1-eta)*_bottom[0][i];
            _bottom[1][i] = eta*fuBottom[1][i] + (1-eta)*_bottom[1][i];
        }
        
        _filter = complex_div(_top,_bottom);
        
        return true;
    }
}