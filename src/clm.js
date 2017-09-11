/**
 * clmtrackr library (https://www.github.com/auduno/clmtrackr/)
 *
 * Copyright (c) 2013, Audun Mathias Øygard
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

import numeric from 'numeric';
import raf from 'raf';
import Promise from 'promise-polyfill';

import emitEvent from './utils/events.js'
import faceDetection from './facedetector/faceDetection.js';
import svmFilter from './svmfilter/svmfilter_fft.js';
import webglFilter from './svmfilter/svmfilter_webgl.js';
import mosseFilterResponses from './mossefilter/mosseFilterResponses.js';
import pModel from '../models/model_pca_20_svm.js';
import canRenderToFloatTexture from './utils/webgl_tests.js';

import { version } from '../package.json';

//import { drawPatches } from './utils.debugging.js';

var DEFAULT_MODEL = pModel;

// polyfills
raf.polyfill();
if (!window.Promise) window.Promise = Promise;

var clm = {
	tracker : function(params) {

		if (!params) params = {};
		if (params.constantVelocity === undefined) params.constantVelocity = true;
		if (params.searchWindow === undefined) params.searchWindow = 11;
		if (params.useWebGL === undefined) params.useWebGL = true;
		if (params.scoreThreshold === undefined) params.scoreThreshold = 0.5;
		if (params.stopOnConvergence === undefined) params.stopOnConvergence = false;
		if (params.weightPoints === undefined) params.weightPoints = undefined;
		if (params.sharpenResponse === undefined) params.sharpenResponse = false;
		if (params.faceDetection === undefined) params.faceDetection = {};

		/** @type {Number} Minimum convergence before firing `clmtrackrConverged` event. */
		var convergenceThreshold = 0.5;

		var numPatches, patchSize, numParameters, patchType;
		var gaussianPD;
		var eigenVectors, eigenValues;
		var sketchCC, sketchW, sketchH, sketchCanvas;
		var weights, model, biases;

		var sobelInit = false;
		var lbpInit = false;

		var currentParameters = [];
		var currentPositions = [];
		var previousParameters = [];
		var previousPositions = [];

		var patches = [];
		var responses = [];
		var meanShape = [];

		var responseMode = 'single';
		var responseList = ['raw'];
		var responseIndex = 0;

		/*
		It's possible to experiment with the sequence of variances used for the finding the maximum in the KDE.
		This sequence is pretty arbitrary, but was found to be okay using some manual testing.
		*/
		var varianceSeq = [10,5,1];
		//var varianceSeq = [3,1.5,0.75];
		//var varianceSeq = [6,3,0.75];
		var PDMVariance = 0.7;

		var relaxation = 0.1;

		var first = true;
		var detectingFace = false;

		var convergenceLimit = 0.01;

		var searchWindow;
		var modelWidth, modelHeight;
		var halfSearchWindow, vecProbs, responsePixels;

		if(typeof Float64Array !== 'undefined') {
			var updatePosition = new Float64Array(2);
			var vecpos = new Float64Array(2);
		} else {
			var updatePosition = new Array(2);
			var vecpos = new Array(2);
		}
		var pw, pl, pdataLength;

		var facecheck_count = 0;

		var webglFi, svmFi, mosseCalc;

		var scoringCanvas = document.createElement('canvas');
		var scoringContext = scoringCanvas.getContext('2d');
		var msxmin, msymin, msxmax, msymax;
		var msmodelwidth, msmodelheight;
		var scoringWeights, scoringBias;
		var scoringHistory = [];
		var meanscore = 0;

		var runnerTimeout, runnerElement, runnerBox;

		var pointWeights;

		var halfPI = Math.PI/2;

		var faceDetector;

		/*
		 *	load model data, initialize filters, etc.
		 *
		 *	@param	<Object>	pdm model object
		 */
		this.init = function(pdmmodel) {
			// default model is pca 20 svm model
			if (pdmmodel === undefined) pdmmodel = DEFAULT_MODEL;

			model = pdmmodel;

			// load from model
			patchType = model.patchModel.patchType;
			numPatches = model.patchModel.numPatches;
			patchSize = model.patchModel.patchSize[0];
			if (patchType == 'MOSSE') {
				searchWindow = patchSize;
			} else {
				searchWindow = params.searchWindow;
			}
			numParameters = model.shapeModel.numEvalues;
			modelWidth = model.patchModel.canvasSize[0];
			modelHeight = model.patchModel.canvasSize[1];

			// set up canvas to work on
			sketchCanvas = document.createElement('canvas');
			sketchCC = sketchCanvas.getContext('2d');

			sketchW = sketchCanvas.width = modelWidth + (searchWindow-1) + patchSize-1;
			sketchH = sketchCanvas.height = modelHeight + (searchWindow-1) + patchSize-1;

			// load eigenvectors
			eigenVectors = numeric.rep([numPatches*2,numParameters],0.0);
			for (var i = 0;i < numPatches*2;i++) {
				for (var j = 0;j < numParameters;j++) {
					eigenVectors[i][j] = model.shapeModel.eigenVectors[i][j];
				}
			}

			// load mean shape
			for (var i = 0; i < numPatches;i++) {
				meanShape[i] = [model.shapeModel.meanShape[i][0], model.shapeModel.meanShape[i][1]];
			}

			// get max and mins, width and height of meanshape
			msxmax = msymax = 0;
			msxmin = msymin = 1000000;
			for (var i = 0;i < numPatches;i++) {
				if (meanShape[i][0] < msxmin) msxmin = meanShape[i][0];
				if (meanShape[i][1] < msymin) msymin = meanShape[i][1];
				if (meanShape[i][0] > msxmax) msxmax = meanShape[i][0];
				if (meanShape[i][1] > msymax) msymax = meanShape[i][1];
			}
			msmodelwidth = msxmax-msxmin;
			msmodelheight = msymax-msymin;

			// get scoringweights if they exist
			if (model.scoring) {
				scoringWeights = new Float64Array(model.scoring.coef);
				scoringBias = model.scoring.bias;
				scoringCanvas.width = model.scoring.size[0];
				scoringCanvas.height = model.scoring.size[1];
			}

			// load eigenvalues
			eigenValues = model.shapeModel.eigenValues;

			weights = model.patchModel.weights;
			biases = model.patchModel.bias;

			// precalculate gaussianPriorDiagonal
			gaussianPD = numeric.rep([numParameters+4, numParameters+4],0);
			// set values and append manual inverse
			for (var i = 0;i < numParameters;i++) {
				if (model.shapeModel.nonRegularizedVectors.indexOf(i) >= 0) {
					gaussianPD[i+4][i+4] = 1/10000000;
				} else {
					gaussianPD[i+4][i+4] = 1/eigenValues[i];
				}
			}

			for (var i = 0;i < numParameters+4;i++) {
				currentParameters[i] = 0;
			}

			if (patchType == 'SVM') {
				var webGLContext;
				var webGLTestCanvas = document.createElement('canvas');
				if (window.WebGLRenderingContext) {
					webGLContext = webGLTestCanvas.getContext('webgl') || webGLTestCanvas.getContext('experimental-webgl');
					if (!webGLContext || !webGLContext.getExtension('OES_texture_float')) {
						webGLContext = null;
					} else {
						// test whether it's possible to render to float texture
						if (!canRenderToFloatTexture(webGLContext)) {
							webGLContext = null;
						}
					}
				}

				if (webGLContext && params.useWebGL && (typeof(webglFilter) !== 'undefined')) {
					webglFi = new webglFilter();
					try {
						webglFi.init(weights, biases, numPatches, searchWindow+patchSize-1, searchWindow+patchSize-1, patchSize, patchSize);
						if ('lbp' in weights) lbpInit = true;
						if ('sobel' in weights) sobelInit = true;
					}
					catch(err) {
						console.error(err);
						alert('There was a problem setting up webGL programs, falling back to slightly slower javascript version. :(');
						webglFi = undefined;
						svmFi = new svmFilter();
						svmFi.init(weights['raw'], biases['raw'], numPatches, patchSize, searchWindow);
					}
				} else if (typeof(svmFilter) !== 'undefined') {
					// use fft convolution if no webGL is available
					svmFi = new svmFilter();
					svmFi.init(weights['raw'], biases['raw'], numPatches, patchSize, searchWindow);
				} else {
					throw new Error('Could not initiate filters, please make sure that svmfilter.js or svmfilter_conv_js.js is loaded.');
				}
			} else if (patchType == 'MOSSE') {
				mosseCalc = new mosseFilterResponses();
				mosseCalc.init(weights, numPatches, patchSize, patchSize);
			}

			if (patchType == 'SVM') {
				pw = pl = patchSize+searchWindow-1;
			} else {
				pw = pl = searchWindow;
			}
			pdataLength = pw*pl;
			halfSearchWindow = (searchWindow-1)/2;
			responsePixels = searchWindow*searchWindow;
			if(typeof Float64Array !== 'undefined') {
				vecProbs = new Float64Array(responsePixels);
				for (var i = 0;i < numPatches;i++) {
					patches[i] = new Float64Array(pdataLength);
				}
			} else {
				vecProbs = new Array(responsePixels);
				for (var i = 0;i < numPatches;i++) {
					patches[i] = new Array(pdataLength);
				}
			}

			if (params.weightPoints) {
				// weighting of points
				pointWeights = [];
				for (var i = 0;i < numPatches;i++) {
					if (i in params.weightPoints) {
						pointWeights[(i*2)] = params.weightPoints[i];
						pointWeights[(i*2)+1] = params.weightPoints[i];
					} else {
						pointWeights[(i*2)] = 1;
						pointWeights[(i*2)+1] = 1;
					}
				}
				pointWeights = numeric.diag(pointWeights);
			}

			faceDetector = new faceDetection(model, params.faceDetection);
		}

		/*
		 *	starts the tracker to run on a regular interval
		 */
		this.start = function(element, box) {
			// check if model is initalized, else return false
			if (typeof(model) === 'undefined') {
				console.log('tracker needs to be initalized before starting to track.');
				return false;
			}
			//check if a runnerelement already exists, if not, use passed parameters
			if (typeof(runnerElement) === 'undefined') {
				runnerElement = element;
				runnerBox = box;
			}

			faceDetector.init(element);

			// start named timeout function
			runnerTimeout = requestAnimationFrame(runnerFunction);
		}

		var runnerFunction = function() {
			runnerTimeout = requestAnimationFrame(runnerFunction);
			// schedule as many iterations as we can during each request
			var startTime = (new Date()).getTime();
			while (((new Date()).getTime() - startTime) < 16) {
				var tracking = this.track(runnerElement, runnerBox);
				if (!tracking) break;
			}
		}.bind(this);

		/*
		 *	stop the running tracker
		 */
		this.stop = function() {
			// stop the running tracker if any exists
			cancelAnimationFrame(runnerTimeout);
		}

		/*
		 *  element : canvas or video element
		 *  TODO: should be able to take img element as well
		 */
		this.track = function(element, box) {
			emitEvent('clmtrackrBeforeTrack');

			var scaling, translateX, translateY, rotation;
			var ptch, px, py;

			if (first) {
				if (!detectingFace) {
					detectingFace = true;

					// this returns a Promise
					faceDetector.getInitialPosition(box)
						.then(function (result) {
							scaling = result[0];
							rotation = result[1];
							translateX = result[2];
							translateY = result[3];

							currentParameters[0] = (scaling*Math.cos(rotation))-1;
							currentParameters[1] = (scaling*Math.sin(rotation));
							currentParameters[2] = translateX;
							currentParameters[3] = translateY;

							currentPositions = calculatePositions(currentParameters, true);

							first = false;
							detectingFace = false;
						})
						.catch(function (error) {
							// send an event on no face found
							emitEvent('clmtrackrNotFound');

							detectingFace = false;
						});
				}

				return false;
			} else {
				facecheck_count += 1;

				if (params.constantVelocity) {
					// calculate where to get patches via constant velocity prediction
					if (previousParameters.length >= 2) {
						for (var i = 0;i < currentParameters.length;i++) {
							currentParameters[i] = (relaxation)*previousParameters[1][i] + (1-relaxation)*((2*previousParameters[1][i]) - previousParameters[0][i]);
							//currentParameters[i] = (3*previousParameters[2][i]) - (3*previousParameters[1][i]) + previousParameters[0][i];
						}
					}
				}

				// change translation, rotation and scale parameters
				rotation = halfPI - Math.atan((currentParameters[0]+1)/currentParameters[1]);
				if (rotation > halfPI) {
					rotation -= Math.PI;
				}
				scaling = currentParameters[1] / Math.sin(rotation);
				translateX = currentParameters[2];
				translateY = currentParameters[3];
			}

			// copy canvas to a new dirty canvas
			sketchCC.save();

			// clear canvas
			sketchCC.clearRect(0, 0, sketchW, sketchH);

			sketchCC.scale(1/scaling, 1/scaling);
			sketchCC.rotate(-rotation);
			sketchCC.translate(-translateX, -translateY);

			sketchCC.drawImage(element, 0, 0, element.width, element.height);

			sketchCC.restore();
			//	get cropped images around new points based on model parameters (not scaled and translated)
			var patchPositions = calculatePositions(currentParameters, false);

			// check whether tracking is ok
			if (scoringWeights && (facecheck_count % 10 == 0)) {
				if (!checkTracking()) {
					// reset all parameters
					resetParameters();

					// send event to signal that tracking was lost
					emitEvent('clmtrackrLost');

					return false;
				}
			}


			var pdata, pmatrix, grayscaleColor;
			for (var i = 0; i < numPatches; i++) {
				px = patchPositions[i][0]-(pw/2);
				py = patchPositions[i][1]-(pl/2);
				ptch = sketchCC.getImageData(Math.round(px), Math.round(py), pw, pl);
				pdata = ptch.data;

				// convert to grayscale
				pmatrix = patches[i];
				for (var j = 0;j < pdataLength;j++) {
					grayscaleColor = pdata[j*4]*0.3 + pdata[(j*4)+1]*0.59 + pdata[(j*4)+2]*0.11;
					pmatrix[j] = grayscaleColor;
				}
			}

			// draw weights for debugging
			//drawPatches(sketchCC, weights, patchSize, patchPositions, function(x) {return x*2000+127});

			// draw patches for debugging
			//drawPatches(sketchCC, patches, pw, patchPositions, false, [27,32,44,50]);

			if (patchType == 'SVM') {
				if (typeof(webglFi) !== 'undefined') {
					responses = getWebGLResponses(patches);
				} else if (typeof(svmFi) !== 'undefined') {
					responses = svmFi.getResponses(patches);
				} else {
					throw new Error('SVM-filters do not seem to be initiated properly.');
				}
			} else if (patchType == 'MOSSE') {
				responses = mosseCalc.getResponses(patches);
			}

			// option to increase sharpness of responses
			if (params.sharpenResponse) {
				for (var i = 0;i < numPatches;i++) {
					for (var j = 0;j < responses[i].length;j++) {
						responses[i][j] = Math.pow(responses[i][j], params.sharpenResponse);
					}
				}
			}

			// draw responses for debugging
			//drawPatches(sketchCC, responses, searchWindow, patchPositions, function(x) {return x*255});

			// iterate until convergence or max 10, 20 iterations?:
			var originalPositions = currentPositions;
			var jac;
			var meanshiftVectors = [];

			for (var i = 0; i < varianceSeq.length; i++) {

				// calculate jacobian
				jac = createJacobian(currentParameters, eigenVectors);

				// for debugging
				//var debugMVs = [];

				var opj0, opj1;

				for (var j = 0;j < numPatches;j++) {
					opj0 = originalPositions[j][0]-((searchWindow-1)*scaling/2);
					opj1 = originalPositions[j][1]-((searchWindow-1)*scaling/2);

					// calculate PI x gaussians
					var vpsum = gpopt(searchWindow, currentPositions[j], updatePosition, vecProbs, responses, opj0, opj1, j, varianceSeq[i], scaling);

					// calculate meanshift-vector
					gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1, scaling);
					//var debugMatrixMV = gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1);

					meanshiftVectors[j] = [vecpos[0] - currentPositions[j][0], vecpos[1] - currentPositions[j][1]];

					//debugMVs[j] = debugMatrixMV;
				}

				// draw meanshiftVector for debugging
				//drawPatches(sketchCC, debugMVs, searchWindow, patchPositions, function(x) {return x*255*500});

				var meanShiftVector = numeric.rep([numPatches*2, 1],0.0);
				for (var k = 0;k < numPatches;k++) {
					meanShiftVector[k*2][0] = meanshiftVectors[k][0];
					meanShiftVector[(k*2)+1][0] = meanshiftVectors[k][1];
				}

				// compute pdm parameter update
				//var prior = numeric.mul(gaussianPD, PDMVariance);
				var prior = numeric.mul(gaussianPD, varianceSeq[i]);
				var jtj;
				if (params.weightPoints) {
					jtj = numeric.dot(numeric.transpose(jac), numeric.dot(pointWeights, jac));
				} else {
					jtj = numeric.dot(numeric.transpose(jac), jac);
				}
				var cpMatrix = numeric.rep([numParameters+4, 1],0.0);
				for (var l = 0;l < (numParameters+4);l++) {
					cpMatrix[l][0] = currentParameters[l];
				}
				var priorP = numeric.dot(prior, cpMatrix);
				var jtv;
				if (params.weightPoints) {
					jtv = numeric.dot(numeric.transpose(jac), numeric.dot(pointWeights, meanShiftVector));
				} else {
					jtv = numeric.dot(numeric.transpose(jac), meanShiftVector);
				}
				var paramUpdateLeft = numeric.add(prior, jtj);
				var paramUpdateRight = numeric.sub(priorP, jtv);

				var paramUpdate = numeric.dot(numeric.inv(paramUpdateLeft), paramUpdateRight);
				//var paramUpdate = numeric.solve(paramUpdateLeft, paramUpdateRight, true);

				var oldPositions = currentPositions;

				// update estimated parameters
				for (var k = 0;k < numParameters+4;k++) {
					currentParameters[k] -= paramUpdate[k];
				}

				// clipping of parameters if they're too high
				var clip;
				for (var k = 0;k < numParameters;k++) {
					clip = Math.abs(3*Math.sqrt(eigenValues[k]));
					if (Math.abs(currentParameters[k+4]) > clip) {
						if (currentParameters[k+4] > 0) {
							currentParameters[k+4] = clip;
						} else {
							currentParameters[k+4] = -clip;
						}
					}

				}

				// update current coordinates
				currentPositions = calculatePositions(currentParameters, true);

				// check if converged
				// calculate norm of parameterdifference
				var positionNorm = 0;
				var pnsq_x, pnsq_y;
				for (var k = 0;k < currentPositions.length;k++) {
					pnsq_x = (currentPositions[k][0]-oldPositions[k][0]);
					pnsq_y = (currentPositions[k][1]-oldPositions[k][1]);
					positionNorm += ((pnsq_x*pnsq_x) + (pnsq_y*pnsq_y));
				}

				// if norm < limit, then break
				if (positionNorm < convergenceLimit) {
					break;
				}

			}

			if (params.constantVelocity) {
				// add current parameter to array of previous parameters
				previousParameters.push(currentParameters.slice());
				if (previousParameters.length == 3) {
					previousParameters.shift();
				}
			}

			// store positions, for checking convergence
			if (previousPositions.length == 10) {
				previousPositions.shift();
			}
			previousPositions.push(currentPositions.slice(0));

			// send an event on each iteration
			emitEvent('clmtrackrIteration');

			// we must get a score before we can say we've converged
			if (scoringHistory.length >= 5 && this.getConvergence() < convergenceThreshold) {
				if (params.stopOnConvergence) {
					this.stop();
				}

				emitEvent('clmtrackrConverged');
			}

			// return new points
			return currentPositions;
		}

		function resetParameters() {
			first = true;
			scoringHistory = [];
			previousParameters = [];
			currentPositions = [];
			previousPositions = [];
			for (var i = 0;i < currentParameters.length;i++) {
				currentParameters[i] = 0;
			}
		}

		/*
		 *	reset tracking, so that track() will start a new detection
		 */
		this.reset = function() {
			resetParameters();
			runnerElement = undefined;
			runnerBox = undefined;
		}

		/*
		 *	draw model on given canvas
		 */
		this.draw = function(canvas, pv, path) {
			// if no previous points, just draw in the middle of canvas

			var params;
			if (pv === undefined) {
				params = currentParameters.slice(0);
			} else {
				params = pv.slice(0);
			}

			var cc = canvas.getContext('2d');
			cc.fillStyle = 'rgb(200,200,200)';
			cc.strokeStyle = 'rgb(130,255,50)';
			//cc.lineWidth = 1;

			var paths;
			if (path === undefined) {
				paths = model.path.normal;
			} else {
				paths = model.path[path];
			}

			for (var i = 0;i < paths.length;i++) {
				if (typeof(paths[i]) == 'number') {
					drawPoint(cc, paths[i], params);
				} else {
					drawPath(cc, paths[i], params);
				}
			}
		}

		/*
		 * 	get the score of the current model fit
		 *	(based on svm of face according to current model)
		 */
		this.getScore = function() {
			return meanscore;
		}

		/*
		 *	calculate positions based on parameters
		 */
		this.calculatePositions = function(parameters) {
			return calculatePositions(parameters, true);
		}

		/*
		 *	get coordinates of current model fit
		 */
		this.getCurrentPosition = function() {
			if (first) {
				return false;
			} else {
				return currentPositions;
			}
		}

		/*
		 *	get parameters of current model fit
		 */
		this.getCurrentParameters = function() {
			return currentParameters;
		}

		/*
		 *	Get the average of recent model movements
		 *	Used for checking whether model fit has converged
		 */
		this.getConvergence = function() {
			if (previousPositions.length < 10) return 999999;

			var prevX = 0.0;
			var prevY = 0.0;
			var currX = 0.0;
			var currY = 0.0;

			// average 5 previous positions
			for (var i = 0;i < 5;i++) {
				for (var j = 0;j < numPatches;j++) {
					prevX += previousPositions[i][j][0];
					prevY += previousPositions[i][j][1];
				}
			}
			prevX /= 5;
			prevY /= 5;

			// average 5 positions before that
			for (var i = 5;i < 10;i++) {
				for (var j = 0;j < numPatches;j++) {
					currX += previousPositions[i][j][0];
					currY += previousPositions[i][j][1];
				}
			}
			currX /= 5;
			currY /= 5;

			// calculate difference
			var diffX = currX-prevX;
			var diffY = currY-prevY;
			var msavg = ((diffX*diffX) + (diffY*diffY));
			msavg /= previousPositions.length;
			return msavg;
		}

		/*
		 * Set response mode (only useful if webGL is available)
		 * mode : either "single", "blend" or "cycle"
		 * list : array of values "raw", "sobel", "lbp"
		 */
		this.setResponseMode = function(mode, list) {
			// clmtrackr must be initialized with model first
			if (typeof(model) === 'undefined') {
				console.log('Clmtrackr has not been initialized with a model yet. No changes made.');
				return;
			}
			// must check whether webGL or not
			if (typeof(webglFi) === 'undefined') {
				console.log('Responsemodes are only allowed when using webGL. In pure JS, only "raw" mode is available.');
				return;
			}
			if (['single', 'blend', 'cycle'].indexOf(mode) < 0) {
				console.log('Tried to set an unknown responsemode : "'+mode+'". No changes made.');
				return;
			}
			if (!(list instanceof Array)) {
				console.log('List in setResponseMode must be an array of strings! No changes made.');
				return;
			} else {
				for (var i = 0;i < list.length;i++) {
					if (['raw', 'sobel', 'lbp'].indexOf(list[i]) < 0) {
						console.log('Unknown element in responsemode list : "'+list[i]+'". No changes made.');
					}
					// check whether filters are initialized
					if (list[i] == 'sobel' && sobelInit == false) {
						console.log('The sobel filters have not been initialized! No changes made.');
					}
					if (list[i] == 'lbp' && lbpInit == false) {
						console.log('The LBP filters have not been initialized! No changes made.');
					}
				}
			}
			// reset index
			responseIndex = 0;
			responseMode = mode;
			responseList = list;
		}


		var getWebGLResponsesType = function(type, patches) {
			if (type == 'lbp') {
				return webglFi.getLBPResponses(patches);
			} else if (type == 'raw') {
				return webglFi.getRawResponses(patches);
			} else if (type == 'sobel') {
				return webglFi.getSobelResponses(patches);
			}
		}

		var getWebGLResponses = function(patches) {
			if (responseMode == 'single') {
				return getWebGLResponsesType(responseList[0], patches);
			} else if (responseMode == 'cycle') {
				var response = getWebGLResponsesType(responseList[responseIndex], patches);
				responseIndex++;
				if (responseIndex >= responseList.length) responseIndex = 0;
				return response;
			} else {
				// blend
				var responses = [];
				for (var i = 0;i < responseList.length;i++) {
					responses[i] = getWebGLResponsesType(responseList[i], patches);
				}
				var blendedResponses = [];
				var searchWindowSize = searchWindow * searchWindow;
				for (var i = 0;i < numPatches;i++) {
					var response = Array(searchWindowSize);
					for (var k = 0;k < searchWindowSize;k++) response[k] = 0;
					for (var j = 0;j < responseList.length;j++) {
						for (var k = 0;k < searchWindowSize;k++) {
							response[k] += (responses[j][i][k]/responseList.length);
						}
					}
					blendedResponses[i] = response;
				}
				return blendedResponses;
			}
		}

		// generates the jacobian matrix used for optimization calculations
		var createJacobian = function(parameters, eigenVectors) {

			var jacobian = numeric.rep([2*numPatches, numParameters+4],0.0);
			var j0,j1;
			for (var i = 0;i < numPatches;i ++) {
				// 1
				j0 = meanShape[i][0];
				j1 = meanShape[i][1];
				for (var p = 0;p < numParameters;p++) {
					j0 += parameters[p+4]*eigenVectors[i*2][p];
					j1 += parameters[p+4]*eigenVectors[(i*2)+1][p];
				}
				jacobian[i*2][0] = j0;
				jacobian[(i*2)+1][0] = j1;
				// 2
				j0 = meanShape[i][1];
				j1 = meanShape[i][0];
				for (var p = 0;p < numParameters;p++) {
					j0 += parameters[p+4]*eigenVectors[(i*2)+1][p];
					j1 += parameters[p+4]*eigenVectors[i*2][p];
				}
				jacobian[i*2][1] = -j0;
				jacobian[(i*2)+1][1] = j1;
				// 3
				jacobian[i*2][2] = 1;
				jacobian[(i*2)+1][2] = 0;
				// 4
				jacobian[i*2][3] = 0;
				jacobian[(i*2)+1][3] = 1;
				// the rest
				for (var j = 0;j < numParameters;j++) {
					j0 = parameters[0]*eigenVectors[i*2][j] - parameters[1]*eigenVectors[(i*2)+1][j] + eigenVectors[i*2][j];
					j1 = parameters[0]*eigenVectors[(i*2)+1][j] + parameters[1]*eigenVectors[i*2][j] + eigenVectors[(i*2)+1][j];
					jacobian[i*2][j+4] = j0;
					jacobian[(i*2)+1][j+4] = j1;
				}
			}

			return jacobian;
		}

		// calculate positions from parameters
		var calculatePositions = function(parameters, useTransforms) {
			var x, y, a, b;
			var numParameters = parameters.length;
			var positions = [];
			for (var i = 0;i < numPatches;i++) {
				x = meanShape[i][0];
				y = meanShape[i][1];
				for (var j = 0;j < numParameters-4;j++) {
					x += model.shapeModel.eigenVectors[(i*2)][j]*parameters[j+4];
					y += model.shapeModel.eigenVectors[(i*2)+1][j]*parameters[j+4];
				}
				if (useTransforms) {
					a = parameters[0]*x - parameters[1]*y + parameters[2];
					b = parameters[0]*y + parameters[1]*x + parameters[3];
					x += a;
					y += b;
				}
				positions[i] = [x,y];
			}

			return positions;
		}

		// part one of meanshift calculation
		var gpopt = function(responseWidth, currentPositionsj, updatePosition, vecProbs, responses, opj0, opj1, j, variance, scaling) {
			var pos_idx = 0;
			var vpsum = 0;
			var dx, dy;
			for (var k = 0;k < responseWidth;k++) {
				updatePosition[1] = opj1+(k*scaling);
				for (var l = 0;l < responseWidth;l++) {
					updatePosition[0] = opj0+(l*scaling);

					dx = currentPositionsj[0] - updatePosition[0];
					dy = currentPositionsj[1] - updatePosition[1];
					vecProbs[pos_idx] = responses[j][pos_idx] * Math.exp(-0.5*((dx*dx)+(dy*dy))/(variance*scaling));

					vpsum += vecProbs[pos_idx];
					pos_idx++;
				}
			}

			return vpsum;
		}

		// part two of meanshift calculation
		var gpopt2 = function(responseWidth, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1, scaling) {
			//for debugging
			//var vecmatrix = [];

			var pos_idx = 0;
			var vecsum = 0;
			vecpos[0] = 0;
			vecpos[1] = 0;
			for (var k = 0;k < responseWidth;k++) {
				updatePosition[1] = opj1+(k*scaling);
				for (var l = 0;l < responseWidth;l++) {
					updatePosition[0] = opj0+(l*scaling);
					vecsum = vecProbs[pos_idx]/vpsum;

					//vecmatrix[k*responseWidth + l] = vecsum;

					vecpos[0] += vecsum*updatePosition[0];
					vecpos[1] += vecsum*updatePosition[1];
					pos_idx++;
				}
			}
			//return vecmatrix;
		}

		// calculate score of current fit
		var checkTracking = function() {
			var trackingImgW = 20;
			var trackingImgH = 22;

			scoringContext.drawImage(sketchCanvas, Math.round(msxmin+(msmodelwidth/4.5)), Math.round(msymin-(msmodelheight/12)), Math.round(msmodelwidth-(msmodelwidth*2/4.5)), Math.round(msmodelheight-(msmodelheight/12)), 0, 0, trackingImgW, trackingImgH);
			// getImageData of canvas
			var imgData = scoringContext.getImageData(0,0,trackingImgW,trackingImgH);
			// convert data to grayscale
			var trackingImgSize = trackingImgW * trackingImgH;
			var scoringData = new Array(trackingImgSize);
			var scdata = imgData.data;
			var scmax = 0;
			for (var i = 0;i < trackingImgSize;i++) {
				scoringData[i] = scdata[i*4]*0.3 + scdata[(i*4)+1]*0.59 + scdata[(i*4)+2]*0.11;
				scoringData[i] = Math.log(scoringData[i]+1);
				if (scoringData[i] > scmax) scmax = scoringData[i];
			}

			if (scmax > 0) {
				// normalize & multiply by svmFilter
				var mean = 0;
				for (var i = 0;i < trackingImgSize;i++) {
					mean += scoringData[i];
				}
				mean /= trackingImgSize;
				var sd = 0;
				for (var i = 0;i < trackingImgSize;i++) {
					sd += (scoringData[i]-mean)*(scoringData[i]-mean);
				}
				sd /= trackingImgSize;
				sd = Math.sqrt(sd);

				var score = 0;
				for (var i = 0;i < trackingImgSize;i++) {
					scoringData[i] = (scoringData[i]-mean)/sd;
					score += (scoringData[i])*scoringWeights[i];
				}
				score += scoringBias;
				score = 1/(1+Math.exp(-score));

				if (scoringHistory.length == 5) {
					scoringHistory.shift();
				}
				scoringHistory.push(score);

				if (scoringHistory.length > 4) {
					// get average
					meanscore = 0;
					for (var i = 0;i < 5;i++) {
						meanscore += scoringHistory[i];
					}
					meanscore /= 5;
					// if below threshold, then reset (return false)
					if (meanscore < params.scoreThreshold) return false;
				}
			}
			return true;
		}

		// draw a parametrized line on a canvas
		var drawPath = function(canvasContext, path, dp) {
			canvasContext.beginPath();
			var i, x, y, a, b;
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = meanShape[i/2][0];
				y = meanShape[i/2][1];
				for (var j = 0;j < numParameters;j++) {
					x += model.shapeModel.eigenVectors[i][j]*dp[j+4];
					y += model.shapeModel.eigenVectors[i+1][j]*dp[j+4];
				}
				a = dp[0]*x - dp[1]*y + dp[2];
				b = dp[0]*y + dp[1]*x + dp[3];
				x += a;
				y += b;

				if (i == 0) {
					canvasContext.moveTo(x,y);
				} else {
					canvasContext.lineTo(x,y);
				}
			}
			canvasContext.moveTo(0,0);
			canvasContext.closePath();
			canvasContext.stroke();
		}

		// draw a point on a canvas
		function drawPoint(canvasContext, point, dp) {
			var i, x, y, a, b;
			i = point*2;
			x = meanShape[i/2][0];
			y = meanShape[i/2][1];
			for (var j = 0;j < numParameters;j++) {
				x += model.shapeModel.eigenVectors[i][j]*dp[j+4];
				y += model.shapeModel.eigenVectors[i+1][j]*dp[j+4];
			}
			a = dp[0]*x - dp[1]*y + dp[2];
			b = dp[0]*y + dp[1]*x + dp[3];
			x += a;
			y += b;
			canvasContext.beginPath();
			canvasContext.arc(x, y, 1, 0, Math.PI*2, true);
			canvasContext.closePath();
			canvasContext.fill();
		}

		return true;
	},
	version : version
}

export default clm;
