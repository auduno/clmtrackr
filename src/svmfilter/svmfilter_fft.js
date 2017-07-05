import FFT from '../utils/fft.js';

var svmFilter = function() {

	var _fft, fft_filters, responses, biases;
	var fft_size, filterLength, filter_width, search_width, num_patches;
	var temp_imag_part, temp_real_part;

	// fft function
	this.fft_inplace = function(array, _im_part) {
		// in-place

		if (typeof _im_part == 'undefined') {
			_im_part = temp_imag_part;
		}

		for (var i = 0;i < filterLength;i++) {
			_im_part[i] = 0.0;
		}

		_fft.real_fft2d(array,_im_part);

		return [array, _im_part];
	}

	this.ifft = function(rn, cn) {
		// in-place
		_fft.real_ifft2d(rn, cn);
		return rn;
	}

	var complex_mult_inplace = function(cn1, cn2) {
		// in-place, cn1 is the one modified
		var temp1, temp2;
		for (var r = 0;r < filterLength;r++) {
			temp1 = (cn1[0][r]*cn2[0][r]) - (cn1[1][r]*cn2[1][r]);
			temp2 = (cn1[0][r]*cn2[1][r]) + (cn1[1][r]*cn2[0][r]);
			cn1[0][r] = temp1;
			cn1[1][r] = temp2;
		}
	}

	this.init = function(filter_input, bias_input, numPatches, filterWidth, searchWidth) {
		// calculate needed size of fft (has to be power of two)
		fft_size = upperPowerOfTwo(filterWidth-1+searchWidth);
		filterLength = fft_size*fft_size;
		_fft = new FFT();
		_fft.init(fft_size);
		fft_filters = Array(numPatches);
		var fft_filter;
		var edge = (filterWidth-1)/2;

		for (var i = 0;i < numPatches;i++) {
			var flar_fi0 = new Float64Array(filterLength);
			var flar_fi1 = new Float64Array(filterLength);

			// load filter
			var xOffset, yOffset;
			for (var j = 0;j < filterWidth;j++) {
				for (var k = 0;k < filterWidth;k++) {
					// TODO : rotate filter

					xOffset = k < edge ? (fft_size-edge) : (-edge);
					yOffset = j < edge ? (fft_size-edge) : (-edge);
					flar_fi0[k+xOffset+((j+yOffset)*fft_size)] = filter_input[i][(filterWidth-1-j)+((filterWidth-1-k)*filterWidth)];

					/*xOffset = k < edge ? (fft_size-edge) : (-edge);
					yOffset = j < edge ? (fft_size-edge) : (-edge);
					flar_fi0[k+xOffset+((j+yOffset)*fft_size)] = filter_input[i][k+(j*filterWidth)];*/

					//console.log(k + ','+ j+':' + (k+xOffset+((j+yOffset)*fft_size)))
				}
			}

			// fft it and store
			fft_filter = this.fft_inplace(flar_fi0, flar_fi1);
			fft_filters[i] = fft_filter;

		}

		// set up biases
		biases = new Float64Array(numPatches);
		for (var i = 0;i < numPatches;i++) {
			biases[i] = bias_input[i];
		}

		responses = Array(numPatches);
		temp_imag_part = Array(numPatches);
		for (var i = 0;i < numPatches;i++) {
			responses[i] = new Float64Array(searchWidth*searchWidth);
			temp_imag_part[i] = new Float64Array(searchWidth*searchWidth);
		}
		temp_real_part = new Float64Array(filterLength);

		num_patches = numPatches;
		filter_width = filterWidth;
		search_width = searchWidth;
	}

	this.getResponses = function(patches) {
		var response, edge;
		var patch_width = filter_width-1+search_width;
		for (var i = 0;i < num_patches;i++) {
			// reset zeroes in temp_real_part
			for (var j = 0;j < fft_size*fft_size;j++) {
				temp_real_part[j] = 0.0;
			}

			// normalize patches to 0-1
			patches[i] = normalizePatches(patches[i]);

			// patch must be padded (with zeroes) to match fft size
			for (var j = 0;j < patch_width;j++) {
				for (var k = 0;k < patch_width;k++) {
					temp_real_part[j + (fft_size*k)] = patches[i][k + (patch_width*j)];
				}
			}

			//drawData(document.getElementById('sketch').getContext('2d'), temp_real_part, 32, 32, false, 0, 0);

			// fft it
			response = this.fft_inplace(temp_real_part);

			// multiply pointwise with filter
			complex_mult_inplace(response, fft_filters[i]);

			// inverse fft it
			response = this.ifft(response[0], response[1]);

			// crop out edges
			edge = (filter_width-1)/2;
			for (var j = 0;j < search_width;j++) {
				for (var k = 0;k < search_width;k++) {
					responses[i][j + (k*search_width)] = response[edge + k + ((j+edge)*(fft_size))];
				}
			}

			// add bias
			for (var j = 0;j < search_width*search_width;j++) {
				responses[i][j] += biases[i];
			}

			// logistic transformation
			responses[i] = logisticResponse(responses[i]);

			/*responses[i] = new Float64Array(32*32)
			for (var j = 0;j < 32;j++) {
				for (var k = 0;k < 32;k++) {
					responses[i][k + (j*(32))] = response[k + (j*(32))]
				}
			}*/

			// normalization?
			inplaceNormalizeFilterMatrix(responses[i]);
		}

		return responses;
	}

	var normalizePatches = function(patch) {
		var patch_width = filter_width-1+search_width;
		var max = 0;
		var min = 1000;
		var value;
		for (var j = 0;j < patch_width;j++) {
			for (var k = 0;k < patch_width;k++) {
				value = patch[k + (patch_width*j)];
				if (value < min) {
					min = value;
				}
				if (value > max) {
					max = value;
				}
			}
		}
		var scale = max-min;
		for (var j = 0;j < patch_width;j++) {
			for (var k = 0;k < patch_width;k++) {
				patch[k + (patch_width*j)] = (patch[k + (patch_width*j)]-min)/scale;
			}
		}
		return patch;
	}

	var logisticResponse = function(response) {
		// create probability by doing logistic transformation
		for (var j = 0;j < search_width;j++) {
			for (var k = 0;k < search_width;k++) {
				response[j + (k*search_width)] = 1.0/(1.0 + Math.exp(- (response[j + (k*search_width)] - 1.0 )));
			}
		}
		return response;
	}

	var upperPowerOfTwo = function(x) {
		x--;
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;
		x++;
		return x;
	}

	var inplaceNormalizeFilterMatrix = function(response) {
		// normalize responses to lie within [0,1]
		var msize = response.length;
		var max = 0;
		var min = 1;

		for (var i = 0;i < msize;i++) {
			max = response[i] > max ? response[i] : max;
			min = response[i] < min ? response[i] : min;
		}
		var dist = max-min;

		if (dist == 0) {
			//console.log('a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.');
		} else {
			for (var i = 0;i < msize;i++) {
				response[i] = (response[i]-min)/dist;
			}
		}
	}
}

export default svmFilter;
