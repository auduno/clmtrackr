/**
 * Fast Fourier Transform module
 * 1D-FFT/IFFT, 2D-FFT/IFFT (radix-2)
 * https://github.com/wellflat/jslib
 */
(function() {
  var FFT;           // top-level namespace
  var _root = this;  // reference to 'window' or 'global'

  if(typeof exports !== 'undefined') {
    FFT = exports;   // for CommonJS
  } else {
    FFT = _root.FFT = {};
  }

  var version = {
    release: '0.2.0',
    date: '2012-01'
  };
  FFT.toString = function() {
    return "version " + version.release + ", released " + version.date;
  };

  // core operations
  var _n = 0,          // order
      _bitrev = null,  // bit reversal table
      _cstb = null;    // sin/cos table
  var core = {
    init : function(n) {
      if(n !== 0 && (n & (n - 1)) === 0) {
        _n = n;
        core._setVariables();
        core._makeBitReversal();
        core._makeCosSinTable();
      } else {
        throw new Error("init: radix-2 required");
      }
    },
    // 1D-FFT
    fft1d : function(re, im) {
      core.fft(re, im, 1);
    },
    // 1D-IFFT
    ifft1d : function(re, im) {
      var n = 1/_n;
      core.fft(re, im, -1);
      for(var i=0; i<_n; i++) {
        re[i] *= n;
        im[i] *= n;
      }
    },
    // 2D-FFT
    fft2d : function(re, im) {
      var tre = [],
          tim = [],
          i = 0;
      // x-axis
      for(var y=0; y<_n; y++) {
        i = y*_n;
        for(var x1=0; x1<_n; x1++) {
          tre[x1] = re[x1 + i];
          tim[x1] = im[x1 + i];
        }
        core.fft1d(tre, tim);
        for(var x2=0; x2<_n; x2++) {
          re[x2 + i] = tre[x2];
          im[x2 + i] = tim[x2];
        }
      }
      // y-axis
      for(var x=0; x<_n; x++) {
        for(var y1=0; y1<_n; y1++) {
          i = x + y1*_n;
          tre[y1] = re[i];
          tim[y1] = im[i];
        }
        core.fft1d(tre, tim);
        for(var y2=0; y2<_n; y2++) {
          i = x + y2*_n;
          re[i] = tre[y2];
          im[i] = tim[y2];
        }
      }
    },
    // 2D-IFFT
    ifft2d : function(re, im) {
      var tre = [],
          tim = [],
          i = 0;
      // x-axis
      for(var y=0; y<_n; y++) {
        i = y*_n;
        for(var x1=0; x1<_n; x1++) {
          tre[x1] = re[x1 + i];
          tim[x1] = im[x1 + i];
        }
        core.ifft1d(tre, tim);
        for(var x2=0; x2<_n; x2++) {
          re[x2 + i] = tre[x2];
          im[x2 + i] = tim[x2];
        }
      }
      // y-axis
      for(var x=0; x<_n; x++) {
        for(var y1=0; y1<_n; y1++) {
          i = x + y1*_n;
          tre[y1] = re[i];
          tim[y1] = im[i];
        }
        core.ifft1d(tre, tim);
        for(var y2=0; y2<_n; y2++) {
          i = x + y2*_n;
          re[i] = tre[y2];
          im[i] = tim[y2];
        }
      }
    },
    // core operation of FFT
    fft : function(re, im, inv) {
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
    },
    // set variables
    _setVariables : function() {
      _bitrev = [];
      _cstb = [];
      // if(typeof Uint8Array !== 'undefined') {
      //   _bitrev = new Uint8Array(_n);
      // } else {
      //   _bitrev = new Array(_n);
      // }
      // if(typeof Float64Array !== 'undefined') {
      //   _cstb = new Float64Array(_n*1.25);
      // } else {
      //   _cstb = new Array(_n*1.25);
      // }
    },
    // make bit reversal table
    _makeBitReversal : function() {
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
    },
    // make trigonometiric function table
    _makeCosSinTable : function() {
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
  };
  // aliases (public APIs)
  var apis = ['init', 'fft1d', 'ifft1d', 'fft2d', 'ifft2d'];
  for(var i=0; i<apis.length; i++) {
    FFT[apis[i]] = core[apis[i]];
  }
  FFT.fft = core.fft1d;
  FFT.ifft = core.ifft1d;
}).call(this);