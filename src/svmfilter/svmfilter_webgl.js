
import { setupWebGL, loadShader, createProgram} from '../../examples/js/libs/webgl-utils.js';

var webglFilter = function() {

	/*
	 * Textures:
	 * 0 : raw filter
	 * 1 : patches
	 * 2 : finished response
	 * 3 : grad/lbp treated patches
	 * 4 : sobel filter
	 * 5 : lbp filter
	 *
	 * Routing:
	 *         (              )  0/4/5 --\
	 *         (              )          _\|
	 * 1 ----> ( ---------->3 ) ----------> 2
	 *         lbpResponse/      patchResponse
	 *         gradientResponse
	 */

	var gl, canvas;
	var filterWidth, filterHeight, patchWidth, patchHeight, numPatches, canvasWidth, canvasHeight;
	var patchResponseProgram, patchDrawProgram;
	var fbo, numBlocks, patchTex;
	var drawRectBuffer, drawLayerBuffer, drawImageBuffer, rttTexture;
	var texCoordBuffer, texCoordLocation, apositionBuffer;
	var newCanvasWidth, newCanvasBlockHeight, newCanvasHeight;
	var drawOutRectangles, drawOutImages, drawOutLayer;
	var patchCells, textureWidth, textureHeight, patchSize, patchArray;
	var biases;

	var lbpResponseProgram;
	var lbpTexCoordLocation, lbpTexCoordBuffer, lbpPositionLocation, lbpAPositionBuffer;

	var gradientResponseProgram;
	var gbo, gradTexCoordLocation, gradTexCoordBuffer, gradPositionLocation, gradAPositionBuffer;

	var lbpInit = false;
	var sobelInit = false;
	var rawInit = false;

	var lbpResponseVS = [
		'attribute vec2 a_texCoord;',
		'attribute vec2 a_position;',
		'',
		'varying vec2 v_texCoord;',
		'',
		'void main() {',
		'   // transform coordinates to regular coordinates',
		'   gl_Position = vec4(a_position,0.0,1.0);',
		' ',
		'   // pass the texCoord to the fragment shader',
		'   v_texCoord = a_texCoord;',
		'}'
	].join('\n');
	var lbpResponseFS;

	var gradientResponseVS = [
		'attribute vec2 a_texCoord;',
		'attribute vec2 a_position;',
		'',
		'varying vec2 v_texCoord;',
		'',
		'void main() {',
		'   // transform coordinates to regular coordinates',
		'   gl_Position = vec4(a_position,0.0,1.0);',
		' ',
		'   // pass the texCoord to the fragment shader',
		'   v_texCoord = a_texCoord;',
		'}'
	].join('\n');
	var gradientResponseFS;

	var patchResponseVS;
	var patchResponseFS;

	var drawResponsesVS = [
		'attribute vec2 a_texCoord_draw;',
		'attribute vec2 a_position_draw;',
		'attribute float a_patchChoice_draw;',
		'',
		'uniform vec2 u_resolutiondraw;',
		'',
		'varying vec2 v_texCoord;',
		'varying float v_select;',
		'',
		'void main() {',
		'   // convert the rectangle from pixels to 0.0 to 1.0',
		'   vec2 zeroToOne = a_position_draw / u_resolutiondraw;',
		'',
		'   // convert from 0->1 to 0->2',
		'   vec2 zeroToTwo = zeroToOne * 2.0;',
		'',
		'   // convert from 0->2 to -1->+1 (clipspace)',
		'   vec2 clipSpace = zeroToTwo - 1.0;',
		'   ',
		'   // transform coordinates to regular coordinates',
		'   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);',
		'',
		'   // pass the texCoord to the fragment shader',
		'   v_texCoord = a_texCoord_draw;',
		'   ',
		'   v_select = a_patchChoice_draw;',
		'}'
	].join('\n');

	var drawResponsesFS = [
		'precision mediump float;',
		'',
		'// our responses',
		'uniform sampler2D u_responses;',
		'',
		'// the texCoords passed in from the vertex shader.',
		'varying vec2 v_texCoord;',
		'varying float v_select;',
		'',
		'const vec4 bit_shift = vec4(256.0*256.0*256.0, 256.0*256.0, 256.0, 1.0);',
		'const vec4 bit_mask  = vec4(0.0, 1.0/256.0, 1.0/256.0, 1.0/256.0);',
		'',
		'// packing code from here http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work',
		'void main() {',
		'  vec4 colorSum = texture2D(u_responses, v_texCoord);',
		'  float value = 0.0;',
		'  if (v_select < 0.1) {',
		'    value = colorSum[0];',
		'  } else if (v_select > 0.9 && v_select < 1.1) {',
		'    value = colorSum[1];',
		'  } else if (v_select > 1.9 && v_select < 2.1) {',
		'    value = colorSum[2];',
		'  } else if (v_select > 2.9 && v_select < 3.1) {',
		'    value = colorSum[3];',
		'  } else {',
		'    value = 1.0;',
		'  }',
		'  ',
		'  vec4 res = fract(value * bit_shift);',
		'  res -= res.xxyz * bit_mask;',
		'  ',
		'  //gl_FragColor = vec4(value, value, value, value);',
		'  //gl_FragColor = vec4(1.0, value, 1.0, 1.0);',
		'  gl_FragColor = res;',
		'}'
	].join('\n');

	this.init = function(filters, bias, nP, pW, pH, fW, fH) {
		// we assume filterVector goes from left to right, rowwise, i.e. row-major order

		if (fW != fH) {
			alert('filter width and height must be same size!');
			return;
		}

		// if filter width is not odd, alert
		if (fW % 2 == 0 || fH % 2 == 0) {
			alert('filters used in svm must be of odd dimensions!');
			return;
		}

		// setup variables
		biases = bias;
		filterWidth = fW;
		filterHeight = fH;
		patchWidth = pW;
		patchHeight = pH;
		numPatches = nP;
		numBlocks = Math.floor(numPatches / 4) + Math.ceil((numPatches % 4)/4);
		canvasWidth = patchWidth;
		canvasHeight = patchHeight*numBlocks;
		newCanvasWidth = patchWidth-filterWidth+1;
		newCanvasBlockHeight = patchHeight-filterWidth+1;
		newCanvasHeight = newCanvasBlockHeight*numPatches;
		patchCells = (Math.floor(numPatches / 4) + Math.ceil((numPatches % 4)/4));
		textureWidth = patchWidth;
		textureHeight = patchHeight*patchCells;
		patchSize = patchWidth*patchHeight;
		patchArray = new Float32Array(patchSize*patchCells*4);
		var opp = [1/patchWidth, 1/(patchHeight*numBlocks)];

		// write out shaders
		patchResponseFS = [
			'precision mediump float;',
			'',
			'const vec2 u_onePixelPatches = vec2('+(1/patchWidth).toFixed(10)+','+(1/(patchHeight*numBlocks)).toFixed(10)+');',
			'const vec2 u_onePixelFilters = vec2('+(1/filterWidth).toFixed(10)+','+(1/(filterHeight*numBlocks)).toFixed(10)+');',
			'const float u_halffilterwidth = '+((filterWidth-1.0)/2).toFixed(1)+';',
			'const float u_halffilterheight = '+((filterHeight-1.0)/2).toFixed(1)+';',
			'',
			'// our patches',
			'uniform sampler2D u_patches;',
			'// our filters',
			'uniform sampler2D u_filters;',
			'',
			'// the texCoords passed in from the vertex shader.',
			'varying vec2 v_texCoord;',
			'varying vec2 v_texCoordFilters; // this should give us correct filter',
			'',
			'void main() {',
			'  vec4 colorSum = vec4(0.0, 0.0, 0.0, 0.0);',
			'  vec4 maxn = vec4(0.0, 0.0, 0.0, 0.0);',
			'  vec4 minn = vec4(256.0, 256.0, 256.0, 256.0);',
			'  vec4 scale = vec4(0.0, 0.0, 0.0, 0.0);',
			'  vec4 patchValue = vec4(0.0, 0.0, 0.0, 0.0);',
			'  vec4 filterValue = vec4(0.0, 0.0, 0.0, 0.0);',
			'  vec4 filterTemp = vec4(0.0, 0.0, 0.0, 0.0);',
			'  for (int w = 0;w < '+filterWidth+';w++) {',
			'    for (int h = 0;h < '+filterHeight+';h++) {',
			'      patchValue = texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(float(w)-u_halffilterwidth, float(h)-u_halffilterheight));',
			'      filterValue = texture2D(u_filters, v_texCoordFilters + u_onePixelFilters * vec2(float(w)-u_halffilterwidth, float(h)-u_halffilterheight));',
			'      maxn = max(patchValue, maxn);',
			'      minn = min(patchValue, minn);',
			'      colorSum += patchValue*filterValue;',
			'      filterTemp += filterValue;',
			'    } ',
			'  }',
			'  scale = maxn-minn;',
			'  colorSum = (colorSum-(minn*filterTemp))/scale;',
			'  // logistic transformation',
			'  colorSum = 1.0/(1.0 + exp(- (colorSum) ));',
			'  gl_FragColor = colorSum;',
			'}'
		].join('\n');

		patchResponseVS = [
			'attribute vec2 a_texCoord;',
			'attribute vec2 a_position;',
			'',
			'const vec2 u_resolution = vec2('+canvasWidth.toFixed(1)+','+canvasHeight.toFixed(1)+');',
			'const float u_patchHeight = '+(1/numBlocks).toFixed(10)+';',
			'const float u_filterHeight = '+(1/numBlocks).toFixed(10)+';',
			'const vec2 u_midpoint = vec2(0.5 ,'+(1/(numBlocks*2)).toFixed(10)+');',
			'',
			'varying vec2 v_texCoord;',
			'varying vec2 v_texCoordFilters;',
			'',
			'void main() {',
			'   // convert the rectangle from pixels to 0.0 to 1.0',
			'   vec2 zeroToOne = a_position / u_resolution;',
			'',
			'   // convert from 0->1 to 0->2',
			'   vec2 zeroToTwo = zeroToOne * 2.0;',
			'',
			'   // convert from 0->2 to -1->+1 (clipspace)',
			'   vec2 clipSpace = zeroToTwo - 1.0;',
			'   ',
			'   // transform coordinates to regular coordinates',
			'   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);',
			' ',
			'   // pass the texCoord to the fragment shader',
			'   v_texCoord = a_texCoord;',
			'   ',
			'   // set the filtertexture coordinate based on number filter to use',
			'   v_texCoordFilters = u_midpoint + vec2(0.0, u_filterHeight * floor(a_texCoord[1]/u_patchHeight));',
			'}'
		].join('\n');

		if ('lbp' in filters) {
			// lbpResponseFragment
			lbpResponseFS = [
				'precision mediump float;',
				'',
				'uniform vec2 u_onePixelPatches;',
				'',
				'// our patches',
				'uniform sampler2D u_patches;',
				'',
				'// the texCoords passed in from the vertex shader.',
				'varying vec2 v_texCoord;',
				'',
				'void main() {',
				'  vec4 topLeft = texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', -'+opp[1].toFixed(5)+'));',
				'  vec4 topMid = texture2D(u_patches, v_texCoord + vec2(0.0, -'+opp[1].toFixed(5)+'));',
				'  vec4 topRight = texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', -'+opp[1].toFixed(5)+'));',
				'  vec4 midLeft = texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', 0.0));',
				'  vec4 midMid = texture2D(u_patches, v_texCoord);',
				'  vec4 midRight = texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', 0.0));',
				'  vec4 bottomLeft = texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', '+opp[1].toFixed(5)+'));',
				'  vec4 bottomMid = texture2D(u_patches, v_texCoord + vec2(0.0, '+opp[1].toFixed(5)+'));',
				'  vec4 bottomRight = texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', '+opp[1].toFixed(5)+'));',
				'  vec4 lbp = step(midMid, midRight)*1.0 + step(midMid, topRight)*2.0 + step(midMid, topMid)*4.0;',
				'  lbp = lbp + step(midMid, topLeft)*8.0 + step(midMid, midLeft)*16.0 + step(midMid, bottomLeft)*32.0;',
				'  lbp = lbp + step(midMid, bottomMid)*64.0 + step(midMid, bottomRight)*128.0;',
				'  gl_FragColor = lbp;',
				'}'
			].join('\n');
		}

		if ('sobel' in filters) {
			// gradResponseFragment
			gradientResponseFS = [
				'precision mediump float;',
				'',
				'uniform vec2 u_onePixelPatches;',
				'',
				'// our patches',
				'uniform sampler2D u_patches;',
				'',
				'// the texCoords passed in from the vertex shader.',
				'varying vec2 v_texCoord;',
				'',
				'void main() {',
				'  vec4 bottomLeft = texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', '+opp[1].toFixed(5)+'));',
				'  vec4 bottomRight = texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', '+opp[1].toFixed(5)+'));',
				'  vec4 topLeft = texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', -'+opp[1].toFixed(5)+'));',
				'  vec4 topRight = texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', -'+opp[1].toFixed(5)+'));',
				'  vec4 dx = (',
				'    bottomLeft +',
				'    (texture2D(u_patches, v_texCoord + vec2(-'+opp[0].toFixed(5)+', 0.0))*vec4(2.0,2.0,2.0,2.0)) +',
				'    topLeft -',
				'    bottomRight -',
				'    (texture2D(u_patches, v_texCoord + vec2('+opp[0].toFixed(5)+', 0.0))*vec4(2.0,2.0,2.0,2.0)) -',
				'    topRight)/4.0;',
				'  vec4 dy = (',
				'    bottomLeft +',
				'    (texture2D(u_patches, v_texCoord + vec2(0.0, '+opp[1].toFixed(5)+'))*vec4(2.0,2.0,2.0,2.0)) +',
				'    bottomRight -',
				'    topLeft -',
				'    (texture2D(u_patches, v_texCoord + vec2(0.0, -'+opp[1].toFixed(5)+'))*vec4(2.0,2.0,2.0,2.0)) -',
				'    topRight)/4.0;',
				'  vec4 gradient = sqrt((dx*dx) + (dy*dy));',
				'  gl_FragColor = gradient;',
				'}'
			].join('\n');
		}

		//create webglcanvas
		canvas = document.createElement('canvas');
		canvas.setAttribute('width', (patchWidth-filterWidth+1)+'px');
		canvas.setAttribute('height', ((patchHeight-filterHeight+1)*numPatches)+'px');
		canvas.setAttribute('id', 'renderCanvas');
		canvas.setAttribute('style', 'display:none;');
		//document.body.appendChild(canvas);
		gl = setupWebGL(canvas, {
			premultipliedAlpha: false,
			preserveDrawingBuffer : true,
			antialias : false
		});


		// check for float textures support and fail if not
		if (!gl.getExtension('OES_texture_float')) {
			alert('Your graphics card does not support floating point textures! :(');
			return;
		}

		/** insert filters into textures **/
		if ('raw' in filters) {
			insertFilter(filters['raw'], gl.TEXTURE0);
			rawInit = true;
		}
		if ('sobel' in filters) {
			insertFilter(filters['sobel'], gl.TEXTURE4);
			sobelInit = true;
		}
		if ('lbp' in filters) {
			insertFilter(filters['lbp'], gl.TEXTURE5);
			lbpInit = true;
		}

		/** calculate vertices for calculating responses **/

		// vertex rectangles to draw out
		var rectangles = [];
		var halfFilter = (filterWidth-1)/2;
		var yOffset;
		for (var i = 0;i < numBlocks;i++) {
			yOffset = i*patchHeight;
			//first triangle
			rectangles = rectangles.concat([
				halfFilter, yOffset+halfFilter,
				patchWidth-halfFilter, yOffset+halfFilter,
				halfFilter, yOffset+patchHeight-halfFilter
			]);
			//second triangle
			rectangles = rectangles.concat([
				halfFilter, yOffset+patchHeight-halfFilter,
				patchWidth-halfFilter, yOffset+halfFilter,
				patchWidth-halfFilter, yOffset+patchHeight-halfFilter
			]);
		}
		rectangles = new Float32Array(rectangles);

		// image rectangles to draw out
		var irectangles = [];
		for (var i = 0;i < rectangles.length;i++) {
			if (i % 2 == 0) {
				irectangles[i] = rectangles[i]/canvasWidth;
			} else {
				irectangles[i] = rectangles[i]/canvasHeight;
			}
		}
		irectangles = new Float32Array(irectangles);

		if ('lbp' in filters || 'sobel' in filters) {
			var topCoord = 1.0 - 2/(patchHeight*numBlocks);
			var bottomCoord = 1.0 - 2/numBlocks + 2/(patchHeight*numBlocks);
			var yOffset;
			// calculate position of vertex rectangles for gradient/lbp program
			var gradRectangles = [];
			for (var i = 0;i < numBlocks;i++) {
				yOffset = i * (2/numBlocks);
				//first triangle
				gradRectangles = gradRectangles.concat([
					-1.0, topCoord - yOffset,
					1.0, topCoord - yOffset,
					-1.0, bottomCoord - yOffset
				]);
				//second triangle
				gradRectangles = gradRectangles.concat([
					-1.0, bottomCoord - yOffset,
					1.0, topCoord - yOffset,
					1.0, bottomCoord - yOffset
				]);
			}
			gradRectangles = new Float32Array(gradRectangles);

			topCoord = 1.0 - 1/(patchHeight*numBlocks);
			bottomCoord = 1.0 - 1/numBlocks + 1/(patchHeight*numBlocks);
			// calculate position of image rectangles to draw out
			var gradIRectangles = [];
			for (var i = 0;i < numBlocks;i++) {
				yOffset = i * (1/numBlocks);
				//first triangle
				gradIRectangles = gradIRectangles.concat([
					0.0, topCoord - yOffset,
					1.0, topCoord - yOffset,
					0.0, bottomCoord - yOffset
				]);
				//second triangle
				gradIRectangles = gradIRectangles.concat([
					0.0, bottomCoord - yOffset,
					1.0, topCoord - yOffset,
					1.0, bottomCoord - yOffset
				]);
			}
			gradIRectangles = new Float32Array(gradIRectangles);
		}

		// vertices for drawing out responses

		// drawOutRectangles
		drawOutRectangles = new Float32Array(12*numPatches);
		var yOffset, indexOffset;
		for (var i = 0;i < numPatches;i++) {
			yOffset = i*newCanvasBlockHeight;
			indexOffset = i*12;

			//first triangle
			drawOutRectangles[indexOffset] = 0.0;
			drawOutRectangles[indexOffset+1] = yOffset;
			drawOutRectangles[indexOffset+2] = newCanvasWidth;
			drawOutRectangles[indexOffset+3] = yOffset;
			drawOutRectangles[indexOffset+4] = 0.0;
			drawOutRectangles[indexOffset+5] = yOffset+newCanvasBlockHeight;

			//second triangle
			drawOutRectangles[indexOffset+6] = 0.0;
			drawOutRectangles[indexOffset+7] = yOffset+newCanvasBlockHeight;
			drawOutRectangles[indexOffset+8] = newCanvasWidth;
			drawOutRectangles[indexOffset+9] = yOffset;
			drawOutRectangles[indexOffset+10] = newCanvasWidth;
			drawOutRectangles[indexOffset+11] = yOffset+newCanvasBlockHeight;
		}

		// images
		drawOutImages = new Float32Array(numPatches*12);
		var halfFilterWidth = ((filterWidth-1)/2)/patchWidth;
		var halfFilterHeight = ((filterWidth-1)/2)/(patchHeight*patchCells);
		var patchHeightT = patchHeight / (patchHeight*patchCells);
		for (var i = 0;i < numPatches;i++) {
			yOffset = Math.floor(i / 4)*patchHeightT;
			indexOffset = i*12;

			//first triangle
			drawOutImages[indexOffset] = halfFilterWidth;
			drawOutImages[indexOffset+1] = yOffset+halfFilterHeight;
			drawOutImages[indexOffset+2] = 1.0-halfFilterWidth;
			drawOutImages[indexOffset+3] = yOffset+halfFilterHeight;
			drawOutImages[indexOffset+4] = halfFilterWidth;
			drawOutImages[indexOffset+5] = yOffset+patchHeightT-halfFilterHeight;

			//second triangle
			drawOutImages[indexOffset+6] = halfFilterWidth;
			drawOutImages[indexOffset+7] = yOffset+patchHeightT-halfFilterHeight;
			drawOutImages[indexOffset+8] = 1.0-halfFilterWidth;
			drawOutImages[indexOffset+9] = yOffset+halfFilterHeight;
			drawOutImages[indexOffset+10] = 1.0-halfFilterWidth;
			drawOutImages[indexOffset+11] = yOffset+patchHeightT-halfFilterHeight;
		}

		// layer
		drawOutLayer = new Float32Array(numPatches*6);
		var layernum;
		for (var i = 0;i < numPatches;i++) {
			layernum = i % 4;
			indexOffset = i*6;
			drawOutLayer[indexOffset] = layernum;
			drawOutLayer[indexOffset+1] = layernum;
			drawOutLayer[indexOffset+2] = layernum;
			drawOutLayer[indexOffset+3] = layernum;
			drawOutLayer[indexOffset+4] = layernum;
			drawOutLayer[indexOffset+5] = layernum;
		}

		/** set up programs and load attributes etc **/

		if ('sobel' in filters) {
			var grVertexShader = loadShader(gl, gradientResponseVS, gl.VERTEX_SHADER);
			var grFragmentShader = loadShader(gl, gradientResponseFS, gl.FRAGMENT_SHADER);
			gradientResponseProgram = createProgram(gl, [grVertexShader, grFragmentShader]);
			gl.useProgram(gradientResponseProgram);

			// set up vertices with rectangles
			gradPositionLocation = gl.getAttribLocation(gradientResponseProgram, 'a_position');
			gradAPositionBuffer = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, gradAPositionBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, gradRectangles, gl.STATIC_DRAW);
			gl.enableVertexAttribArray(gradPositionLocation);
			gl.vertexAttribPointer(gradPositionLocation, 2, gl.FLOAT, false, 0, 0);

			// set up texture positions
			gradTexCoordLocation = gl.getAttribLocation(gradientResponseProgram, 'a_texCoord');
			gradTexCoordBuffer = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, gradTexCoordBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, gradIRectangles, gl.STATIC_DRAW);
			gl.enableVertexAttribArray(gradTexCoordLocation);
			gl.vertexAttribPointer(gradTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

			// set up patches texture in gradientResponseProgram
			gl.uniform1i(gl.getUniformLocation(gradientResponseProgram, 'u_patches'), 1);
		}
		if ('lbp' in filters) {
			var lbpVertexShader = loadShader(gl, lbpResponseVS, gl.VERTEX_SHADER);
			var lbpFragmentShader = loadShader(gl, lbpResponseFS, gl.FRAGMENT_SHADER);
			lbpResponseProgram = createProgram(gl, [lbpVertexShader, lbpFragmentShader]);
			gl.useProgram(lbpResponseProgram);

			// set up vertices with rectangles
			lbpPositionLocation = gl.getAttribLocation(lbpResponseProgram, 'a_position');
			lbpAPositionBuffer = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, lbpAPositionBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, gradRectangles, gl.STATIC_DRAW);
			gl.enableVertexAttribArray(lbpPositionLocation);
			gl.vertexAttribPointer(lbpPositionLocation, 2, gl.FLOAT, false, 0, 0);

			// set up texture positions
			gradTexCoordLocation = gl.getAttribLocation(lbpResponseProgram, 'a_texCoord');
			lbpTexCoordBuffer = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, lbpTexCoordBuffer);
			gl.bufferData(gl.ARRAY_BUFFER, gradIRectangles, gl.STATIC_DRAW);
			gl.enableVertexAttribArray(lbpTexCoordLocation);
			gl.vertexAttribPointer(lbpTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

			// set up patches texture in lbpResponseProgram
			gl.uniform1i(gl.getUniformLocation(lbpResponseProgram, 'u_patches'), 1);
		}

		// setup patchdraw program
		var drVertexShader = loadShader(gl, drawResponsesVS, gl.VERTEX_SHADER);
		var drFragmentShader = loadShader(gl, drawResponsesFS, gl.FRAGMENT_SHADER);
		patchDrawProgram = createProgram(gl, [drVertexShader, drFragmentShader]);
		gl.useProgram(patchDrawProgram);

		// set the resolution/dimension of the canvas
		var resolutionLocation = gl.getUniformLocation(patchDrawProgram, 'u_resolutiondraw');
		gl.uniform2f(resolutionLocation, newCanvasWidth, newCanvasHeight);

		// set u_responses
		var responsesLocation = gl.getUniformLocation(patchDrawProgram, 'u_responses');
		gl.uniform1i(responsesLocation, 2);

		// setup patchresponse program
		var prVertexShader = loadShader(gl, patchResponseVS, gl.VERTEX_SHADER);
		var prFragmentShader = loadShader(gl, patchResponseFS, gl.FRAGMENT_SHADER);
		patchResponseProgram = createProgram(gl, [prVertexShader, prFragmentShader]);
		gl.useProgram(patchResponseProgram);

		// set up vertices with rectangles
		var positionLocation = gl.getAttribLocation(patchResponseProgram, 'a_position');
		apositionBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, rectangles, gl.STATIC_DRAW);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		texCoordLocation = gl.getAttribLocation(patchResponseProgram, 'a_texCoord');
		texCoordBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, irectangles, gl.STATIC_DRAW);
		gl.enableVertexAttribArray(texCoordLocation);
		gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

		if ('lbp' in filters || 'sobel' in filters) {
			// set up gradient/lbp buffer (also used for lbp)
			gl.activeTexture(gl.TEXTURE3);
			var gradients = gl.createTexture();
			gl.bindTexture(gl.TEXTURE_2D, gradients);
			gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, patchWidth, patchHeight*numBlocks, 0, gl.RGBA, gl.FLOAT, null);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

			// set up gradient/lbp framebuffer
			gbo = gl.createFramebuffer();
			gl.bindFramebuffer(gl.FRAMEBUFFER, gbo);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, gradients, 0);
		}

		// set up buffer to draw to
		gl.activeTexture(gl.TEXTURE2);
		rttTexture = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, rttTexture);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, patchWidth, patchHeight*numBlocks, 0, gl.RGBA, gl.FLOAT, null);

		// set up response framebuffer
		fbo = gl.createFramebuffer();
		gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
		gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, rttTexture, 0);

		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		/* initialize some textures and buffers used later on */

		patchTex = gl.createTexture();
		drawRectBuffer = gl.createBuffer();
		drawImageBuffer = gl.createBuffer();
		drawLayerBuffer = gl.createBuffer();
	}

	this.getRawResponses = function(patches) {
		// TODO: check patches correct length/dimension

		insertPatches(patches);

		// switch to correct program
		gl.useProgram(patchResponseProgram);

		// set u_patches to point to texture 1
		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_patches'), 1);

		// set u_filters to point to correct filter
		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_filters'), 0);

		// set up vertices with rectangles
		var positionLocation = gl.getAttribLocation(patchResponseProgram, 'a_position');
		gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		var texCoordLocation = gl.getAttribLocation(patchResponseProgram, 'a_texCoord');
		gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
		gl.enableVertexAttribArray(texCoordLocation);
		gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

		// set framebuffer to the original one if not already using it
		gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);

		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		// draw to framebuffer
		gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);

		//gl.finish();

		var responses = drawOut('raw');

		return responses;
	}

	this.getSobelResponses = function(patches) {
		// check that it is initialized
		if (!sobelInit) return;

		insertPatches(patches);

		/* do sobel filter on patches */

		// switch to correct program
		gl.useProgram(gradientResponseProgram);

		// set up vertices with rectangles
		var gradPositionLocation = gl.getAttribLocation(gradientResponseProgram, 'a_position');
		gl.bindBuffer(gl.ARRAY_BUFFER, gradAPositionBuffer);
		gl.enableVertexAttribArray(gradPositionLocation);
		gl.vertexAttribPointer(gradPositionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		var gradTexCoordLocation = gl.getAttribLocation(gradientResponseProgram, 'a_texCoord');
		gl.bindBuffer(gl.ARRAY_BUFFER, gradTexCoordBuffer);
		gl.enableVertexAttribArray(gradTexCoordLocation);
		gl.vertexAttribPointer(gradTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

		// set framebuffer to the original one if not already using it
		gl.bindFramebuffer(gl.FRAMEBUFFER, gbo);

		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		// draw to framebuffer
		gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);

		/* calculate responses */

		gl.useProgram(patchResponseProgram);

		// set patches and filters to point to correct textures
		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_filters'), 4);
		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_patches'), 3);

		var positionLocation = gl.getAttribLocation(patchResponseProgram, 'a_position');
		gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		var texCoordLocation = gl.getAttribLocation(patchResponseProgram, 'a_texCoord');
		gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
		gl.enableVertexAttribArray(texCoordLocation);
		gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

		gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		// draw to framebuffer
		gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);

		/* get the responses */

		var responses = drawOut('sobel');

		return responses;
	}

	this.getLBPResponses = function(patches) {
		// check that it is initialized
		if (!lbpInit) return;

		insertPatches(patches);

		/* do sobel filter on patches */

		// switch to correct program
		gl.useProgram(lbpResponseProgram);

		// set up vertices with rectangles
		var lbpPositionLocation = gl.getAttribLocation(lbpResponseProgram, 'a_position');
		gl.bindBuffer(gl.ARRAY_BUFFER, lbpAPositionBuffer);
		gl.enableVertexAttribArray(lbpPositionLocation);
		gl.vertexAttribPointer(lbpPositionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		var lbpTexCoordLocation = gl.getAttribLocation(lbpResponseProgram, 'a_texCoord');
		gl.bindBuffer(gl.ARRAY_BUFFER, lbpTexCoordBuffer);
		gl.enableVertexAttribArray(lbpTexCoordLocation);
		gl.vertexAttribPointer(lbpTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

		// set framebuffer to the original one if not already using it
		gl.bindFramebuffer(gl.FRAMEBUFFER, gbo);

		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		// draw to framebuffer
		gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);

		/* calculate responses */

		gl.useProgram(patchResponseProgram);

		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_filters'), 5);
		gl.uniform1i(gl.getUniformLocation(patchResponseProgram, 'u_patches'), 3);

		var positionLocation = gl.getAttribLocation(patchResponseProgram, 'a_position');
		gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

		// set up texture positions
		var texCoordLocation = gl.getAttribLocation(patchResponseProgram, 'a_texCoord');
		gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
		gl.enableVertexAttribArray(texCoordLocation);
		gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

		gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
		gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		// draw to framebuffer
		gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);

		/* get the responses */

		var responses = drawOut('lbp');

		return responses;
	}

	var insertPatches = function(patches) {
		// pass patches into texture, each patch in either r, g, b or a
		var patchArrayIndex = 0;
		var patchesIndex1 = 0;
		var patchesIndex2 = 0;
		for (var i = 0;i < patchCells;i++) {
			for (var j = 0;j < patchHeight;j++) {
				for (var k = 0;k < patchWidth;k++) {
					patchesIndex1 = i*4;
					patchesIndex2 = (j*patchWidth) + k;
					patchArrayIndex = ((patchSize*i) + patchesIndex2)*4;

					//set r with first patch
					if (patchesIndex1 < numPatches) {
						patchArray[patchArrayIndex] = patches[patchesIndex1][patchesIndex2];
					} else {
						patchArray[patchArrayIndex] = 0;
					}
					//set g with 2nd patch
					if (patchesIndex1+1 < numPatches) {
						patchArray[patchArrayIndex + 1] = patches[patchesIndex1+1][patchesIndex2];
					} else {
						patchArray[patchArrayIndex + 1] = 0;
					}
					//set b with 3rd patch
					if (patchesIndex1+2 < numPatches) {
						patchArray[patchArrayIndex + 2] = patches[patchesIndex1+2][patchesIndex2];
					} else {
						patchArray[patchArrayIndex + 2] = 0;
					}
					//set a with 4th patch
					if (patchesIndex1+3 < numPatches) {
						patchArray[patchArrayIndex + 3] = patches[patchesIndex1+3][patchesIndex2];
					} else {
						patchArray[patchArrayIndex + 3] = 0;
					}
				}
			}
		}

		// pass texture into an uniform
		gl.activeTexture(gl.TEXTURE1);
		gl.bindTexture(gl.TEXTURE_2D, patchTex);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, textureWidth, textureHeight, 0, gl.RGBA, gl.FLOAT, patchArray);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
	}

	var insertFilter = function(filter, textureNum) {
		var filterSize = filterWidth*filterHeight;
		var filterArray = new Float32Array(filterSize*(numBlocks)*4);
		for (var i = 0;i < numBlocks;i++) {
			for (var j = 0;j < filterHeight;j++) {
				for (var k = 0;k < filterWidth;k++) {
					//set r with first filter
					if (i*4 < filter.length) {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = filter[i*4][(j*filterWidth) + k];
					} else {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = 0;
					}
					//set g with 2nd filter
					if ((i*4 + 1) < filter.length) {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = filter[(i*4)+1][(j*filterWidth) + k];
					} else {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = 0;
					}
					//set b with 3rd filter
					if ((i*4 + 2) < filter.length) {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = filter[(i*4)+2][(j*filterWidth) + k];
					} else {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = 0;
					}
					//set a with 4th filter
					if ((i*4 + 3) < filter.length) {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 3] = filter[(i*4)+3][(j*filterWidth) + k];
					} else {
						filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 3] = 0;
					}
				}
			}
		}

		gl.activeTexture(textureNum);
		var filterTexture = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, filterTexture);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, filterWidth, filterHeight*numBlocks, 0, gl.RGBA, gl.FLOAT, filterArray);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
	}

	var drawOut = function(type) {
		// switch programs
		gl.useProgram(patchDrawProgram);

		// bind canvas buffer
		gl.bindFramebuffer(gl.FRAMEBUFFER, null);
		gl.viewport(0, 0, newCanvasWidth, newCanvasHeight);

		gl.clearColor(0.0, 0.0, 0.0, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER);

		gl.bindBuffer(gl.ARRAY_BUFFER, drawRectBuffer);
		gl.bufferData(
			gl.ARRAY_BUFFER,
			drawOutRectangles,
			gl.STATIC_DRAW);
		var positionLocation = gl.getAttribLocation(patchDrawProgram, 'a_position_draw');
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

		gl.bindBuffer(gl.ARRAY_BUFFER, drawImageBuffer);
		gl.bufferData(
			gl.ARRAY_BUFFER,
			drawOutImages,
			gl.STATIC_DRAW);
		var textureLocation = gl.getAttribLocation(patchDrawProgram, 'a_texCoord_draw');
		gl.enableVertexAttribArray(textureLocation);
		gl.vertexAttribPointer(textureLocation, 2, gl.FLOAT, false, 0, 0);

		gl.bindBuffer(gl.ARRAY_BUFFER, drawLayerBuffer);
		gl.bufferData(
			gl.ARRAY_BUFFER,
			drawOutLayer,
			gl.STATIC_DRAW);
		var layerLocation = gl.getAttribLocation(patchDrawProgram, 'a_patchChoice_draw');
		gl.enableVertexAttribArray(layerLocation);
		gl.vertexAttribPointer(layerLocation, 1, gl.FLOAT, false, 0, 0);

		// draw out
		gl.drawArrays(gl.TRIANGLES, 0, numPatches*6);

		var responses = getOutput();
		responses = unpackToFloat(responses);
		responses = splitArray(responses, numPatches);
		responses = addBias(responses, biases[type]);

		// normalize responses to lie within [0,1]
		var rl = responses.length;
		for (var i = 0;i < rl;i++) {
			responses[i] = normalizeFilterMatrix(responses[i]);
		}

		return responses;
	}

	var addBias = function(responses, bias) {
		// do a little trick to add bias in the logit function
		var biasMult;
		for (var i = 0;i < responses.length;i++) {
			biasMult = Math.exp(bias[i]);
			for (var j = 0;j < responses[i].length;j++) {
				responses[i][j] = 1/(1+((1-responses[i][j])/(responses[i][j]*biasMult)));
			}
		}
		return responses;
	}

	var splitArray = function(array, parts) {
		var sp = [];
		var al = array.length;
		var splitlength = al/parts;
		var ta = [];
		for (var i = 0;i < al;i++) {
			if (i % splitlength == 0) {
				if (i != 0) {
					sp.push(ta);
				}
				ta = [];
			}
			ta.push(array[i]);
		}
		sp.push(ta);
		return sp;
	}

	var getOutput = function() {
		// get data
		var pixelValues = new Uint8Array(4*canvas.width*canvas.height);
		gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, pixelValues);
		return pixelValues;
	}

	var unpackToFloat = function(array) {
		// convert packed floats to proper floats : see http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work
		var newArray = [];
		var al = array.length;
		for (var i = 0;i < al;i+=4) {
			newArray[(i / 4) >> 0] = ((array[i]/(256*256*256*256))+(array[i+1]/(256*256*256))+(array[i+2]/(256*256))+(array[i+3]/256));
		}
		return newArray;
	}

	var normalizeFilterMatrix = function(response) {
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
			response = response.map(function() {return 1});
		} else {
			for (var i = 0;i < msize;i++) {
				response[i] = (response[i]-min)/dist;
			}
		}

		return response;
	}
};

export default webglFilter;