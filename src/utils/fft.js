/**
 * Fast Fourier Transform
 * 1D-FFT/IFFT, 2D-FFT/IFFT (radix-2)
 *
 * @author ryo / github.com/wellflat
 * Based on https://github.com/wellflat/javascript-labs with some tiny optimizations
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
			throw new Error('init: radix-2 required');
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

	// 2D-IFFT, real-valued
	// only outputs the real valued part
	this.real_ifft2d = function (re, im) {
		var i2;
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
		var halfn = _n/2;
		var rowIdx = 0;
		for(var x=0; x<_n; x+=2) {
			//untangle
			i = x;
			i2 = x+1;
			_tre[0] = re[0 + i];
			_tim[0] = re[0 + i2];
			_tre[_n/2] = re[(halfn*_n) + i];
			_tim[_n/2] = re[(halfn*_n) + i2];
			for (var x2=1;x2<halfn;x2++) {
				rowIdx = x2*_n;
				_tre[x2] = re[rowIdx+i] - im[rowIdx + i2];
				_tre[_n - x2] = re[rowIdx+i] + im[rowIdx + i2];
				_tim[x2] = im[rowIdx+i] + re[rowIdx+i2];
				_tim[_n - x2] = re[rowIdx+i2] - im[rowIdx+i];
			}
			this.ifft1d(_tre, _tim);
			for(var y2=0; y2<_n; y2++) {
				i = x + y2*_n;
				i2 = (x + 1) + y2*_n;
				re[i] = _tre[y2];
				re[i2] = _tim[y2];
			}
		}
	}

	// 2D-FFT, real-valued only
	// ignores the imaginary input
	//   see:
	// http://www.inf.fu-berlin.de/lehre/SS12/SP-Par/download/fft1.pdf
	// http://cnx.org/content/m12021/latest/
	// http://images.apple.com/acg/pdf/g4fft.pdf
	// http://www.ti.com/lit/an/spra291/spra291.pdf
	this.real_fft2d = function (re, im) {
		var i = 0, i2 = 0;
		// x-axis
		for(var y=0; y<_n; y += 2) {
			i = y*_n;
			i2 = (y+1)*_n;
			// tangle
			for(var x1=0; x1<_n; x1++) {
				_tre[x1] = re[x1 + i];
				_tim[x1] = re[x1 + i2];
			}
			this.fft1d(_tre, _tim);
			// untangle
			re[0 + i] = _tre[0];
			re[0 + i2] = _tim[0];
			im[0 + i] = 0;
			im[0 + i2] = 0;
			re[_n/2 + i] = _tre[_n/2];
			re[_n/2 + i2] = _tim[_n/2];
			im[_n/2 + i] = 0;
			im[_n/2 + i2] = 0;
			for(var x2=1;x2<(_n/2);x2++) {
				re[x2 + i] = 0.5 * (_tre[x2] + _tre[_n - x2]);
				im[x2 + i] = 0.5 * (_tim[x2] - _tim[_n - x2]);
				re[x2 + i2] = 0.5 * (_tim[x2] + _tim[_n - x2]);
				im[x2 + i2] = -0.5 * (_tre[x2] - _tre[_n - x2]);
				re[(_n-x2) + i] = re[x2 + i];
				im[(_n-x2) + i] = -im[x2 + i];
				re[(_n-x2) + i2] = re[x2 + i2];
				im[(_n-x2) + i2] = -im[x2 + i2];
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
		//butfly(re,im,inv,n4);
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

	function butfly(re, im, inv, n4) {
		var h,d,wr,wi,ik,xr,xi;
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
			_tre = new Float64Array(_n);
			_tim = new Float64Array(_n);
		} else {
			_cstb = new Array(_n*1.25);
			_tre = new Array(_n);
			_tim = new Array(_n);
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

export default FFT;