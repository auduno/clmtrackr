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

"use strict";
//requires: ccv.js, numeric.js

var clm = {
	tracker : function(params) {
		
		if (!params) params = {};
		if (params.constantVelocity === undefined) params.constantVelocity = true;
		if (params.searchWindow === undefined) params.searchWindow = 11;
		if (params.useWebGL === undefined) params.useWebGL = true;
		if (params.scoreThreshold === undefined) params.scoreThreshold = 0.30;
		if (params.stopOnConvergence === undefined) params.stopOnConvergence = false;
		if (params.weightPoints === undefined) params.weightPoints = undefined;
		if (params.sharpenResponse === undefined) params.sharpenResponse = false;
		
		var numPatches, patchSize, numParameters, patchType;
		var gaussianPD;
		var eigenVectors, eigenValues;
		var sketchCC, sketchW, sketchH, sketchCanvas;
		var candidate;
		var weights, model;
		
		var currentParameters = [];
		var currentPositions = [];
		var previousParameters = [];
		var previousPositions = [];
		
		var patches = [];
		var responses = [];
		var meanShape = [];
		
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
		
		var convergenceLimit = 0.01;
		
		var learningRate = [];
		var stepParameter = 1.25;
		var prevCostFunc = []
		
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
		//document.body.appendChild(scoringCanvas);
		var scoringContext = scoringCanvas.getContext('2d');
		var msxmin, msymin, msxmax, msymax;
		var msmodelwidth, msmodelheight;
		var scoringWeights, scoringBias;
		var scoringHistory = [];
		var meanscore = 0;

		var mossef_lefteye, mossef_righteye, mossef_nose;
		var right_eye_position = [0.0,0.0];
		var left_eye_position = [0.0,0.0];
		var nose_position = [0.0,0.0];
		var lep, rep, mep;
		var runnerTimeout, runnerElement, runnerBox;
		
		var pointWeights;

		var halfPI = Math.PI/2;
		
		/*
		 *	load model data, initialize filters, etc.
		 *
		 *	@param	<Object>	pdm model object
		 */
		this.init = function(pdmmodel) {
			
			model = pdmmodel;
			
			// load from model
			patchType = model.patchModel.patchType;
			numPatches = model.patchModel.numPatches;
			patchSize = model.patchModel.patchSize[0];
			if (patchType == "MOSSE") {
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
			
			if (model.hints && mosseFilter && left_eye_filter && right_eye_filter && nose_filter) {
				//var mossef_lefteye = new mosseFilter({drawResponse : document.getElementById('overlay2')});
				mossef_lefteye = new mosseFilter();
				mossef_lefteye.load(left_eye_filter);
				//var mossef_righteye = new mosseFilter({drawResponse : document.getElementById('overlay2')});
				mossef_righteye = new mosseFilter();
				mossef_righteye.load(right_eye_filter);
				//var mossef_nose = new mosseFilter({drawResponse : document.getElementById('overlay2')});
				mossef_nose = new mosseFilter();
				mossef_nose.load(nose_filter);
			} else {
				console.log("MOSSE filters not found, using rough approximation for initialization.");
			}

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
			
			if (patchType == "SVM") {
				var webGLContext;
				var webGLTestCanvas = document.createElement('canvas');
				if (window.WebGLRenderingContext) {
					webGLContext = webGLTestCanvas.getContext('webgl') || webGLTestCanvas.getContext('experimental-webgl');
					if (!webGLContext || !webGLContext.getExtension('OES_texture_float')) {
						webGLContext = null;
					}
				} 
				
				if (webGLContext && params.useWebGL && (typeof(webglFilter) !== "undefined")) {
					webglFi = new webglFilter();
					webglFi.init(weights, numPatches, searchWindow+patchSize-1, searchWindow+patchSize-1, patchSize, patchSize, true);
				} else if (typeof(svmFilter) !== "undefined") {
					// use fft convolution if no webGL is available
					svmFi = new svmFilter();
					svmFi.init(weights, numPatches, patchSize, searchWindow);
				} else {
					throw "Could not initiate filters, please make sure that svmfilter.js or svmfilter_conv_js.js is loaded."
				}
			} else if (patchType == "MOSSE") {
				mosseCalc = new mosseFilterResponses();
				mosseCalc.init(weights, numPatches, patchSize, patchSize);
			}
			
			if (patchType == "SVM") {
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
			
			for (var i = 0;i < numPatches;i++) {
				learningRate[i] = 1.0;
				prevCostFunc[i] = 0.0;
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
		}
		
		/*
		 *	starts the tracker to run on a regular interval
		 */
		this.start = function(element, box) {
			// check if model is initalized, else return false
			if (typeof(model) === "undefined") {
				console.log("tracker needs to be initalized before starting to track.");
				return false;
			}
			//check if a runnerelement already exists, if not, use passed parameters
			if (typeof(runnerElement) === "undefined") {
				runnerElement = element;
				runnerBox = box;
			}
			// start named timeout function
			runnerTimeout = requestAnimFrame(runnerFunction);
		}

		/*
		 *	stop the running tracker
		 */
		this.stop = function() {
			// stop the running tracker if any exists
			cancelRequestAnimFrame(runnerTimeout);
		}

		/*
		 *  element : canvas or video element
		 *  TODO: should be able to take img element as well
		 */
		this.track = function(element, box) {
			
			var scaling, translateX, translateY, rotation;
			var croppedPatches = [];
			var ptch, px, py;
						
			if (first) {
				// do viola-jones on canvas to get initial guess, if we don't have any points
				var gi = getInitialPosition(element, box);
				if (!gi) {
					// send an event on no face found
					var evt = document.createEvent("Event");
					evt.initEvent("clmtrackrNotFound", true, true);
					document.dispatchEvent(evt)
					
					return false;
				}
				scaling = gi[0];
				rotation = gi[1];
				translateX = gi[2];
				translateY = gi[3];
				
				first = false;
			} else {
				facecheck_count += 1;
				
				// TODO : do cross-correlation/correlation-filter or similar to find translation of face
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
					first = true;
					scoringHistory = [];
					for (var i = 0;i < currentParameters.length;i++) {
						currentParameters[i] = 0;
						previousParameters = [];
					}
					
					// send event to signal that tracking was lost
					var evt = document.createEvent("Event");
					evt.initEvent("clmtrackrLost", true, true);
					document.dispatchEvent(evt)
					
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
			
			/*print weights*/
			/*sketchCC.clearRect(0, 0, sketchW, sketchH);
			var nuWeights;
			for (var i = 0;i < numPatches;i++) {
				nuWeights = weights[i].map(function(x) {return x*2000+127;});
				drawData(sketchCC, nuWeights, patchSize, patchSize, false, patchPositions[i][0]-(patchSize/2), patchPositions[i][1]-(patchSize/2));
			}*/
			
			// print patches
			/*sketchCC.clearRect(0, 0, sketchW, sketchH);
			for (var i = 0;i < numPatches;i++) {
				if ([27,32,44,50].indexOf(i) > -1) {
					drawData(sketchCC, patches[i], pw, pl, false, patchPositions[i][0]-(pw/2), patchPositions[i][1]-(pl/2));
				}
			}*/
			
			if (patchType == "SVM") {
				if (typeof(webglFi) !== "undefined") {
					var responses = webglFi.getResponses(patches, true);
				} else if (typeof(svmFi) !== "undefined"){
					var responses = svmFi.getResponses(patches);
				} else {
					throw "SVM-filters do not seem to be initiated properly."
				}
			} else if (patchType == "MOSSE") {
				var responses = mosseCalc.getResponses(patches);
			}

			// option to increase sharpness of responses
			if (params.sharpenResponse) {
				for (var i = 0;i < numPatches;i++) {
					for (var j = 0;j < responses[i].length;j++) {
						responses[i][j] = Math.pow(responses[i][j], params.sharpenResponse);
					}
				}
			}

			// print responses
			/*sketchCC.clearRect(0, 0, sketchW, sketchH);
			var nuWeights;
			for (var i = 0;i < numPatches;i++) {
		
				nuWeights = [];
				for (var j = 0;j < responses[i].length;j++) {
					nuWeights.push(responses[i][j]*255);
				}
				
				//drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[i][0]-((searchWindow-1)/2), patchPositions[i][1]-((searchWindow-1)/2));
				if ([27,32,44,50].indexOf(i) > -1) {
					drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[i][0]-((searchWindow-1)/2), patchPositions[i][1]-((searchWindow-1)/2));
				}
			}*/
			
			// iterate until convergence or max 10, 20 iterations?:
			var originalPositions = currentPositions;
			var jac;
			var meanshiftVectors = [];
			
			for (var i = 0; i < varianceSeq.length; i++) {
				
				// calculate jacobian
				jac = createJacobian(currentParameters, eigenVectors);

				// for debugging
				//var debugMVs = [];
				//
				
				var opj0, opj1;
				
				for (var j = 0;j < numPatches;j++) {
					opj0 = originalPositions[j][0]-((searchWindow-1)*scaling/2);
					opj1 = originalPositions[j][1]-((searchWindow-1)*scaling/2);
					
					// calculate PI x gaussians
					var vpsum = gpopt(searchWindow, currentPositions[j], updatePosition, vecProbs, responses, opj0, opj1, j, varianceSeq[i], scaling);
					
					// calculate meanshift-vector
					gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1, scaling);
					
					// for debugging
					//var debugMatrixMV = gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1);
					
					// evaluate here whether to increase/decrease stepSize
					/*if (vpsum >= prevCostFunc[j]) {
						learningRate[j] *= stepParameter;
					} else {
						learningRate[j] = 1.0;
					}
					prevCostFunc[j] = vpsum;*/
					
					// compute mean shift vectors
					// extrapolate meanshiftvectors
					/*var msv = [];
					msv[0] = learningRate[j]*(vecpos[0] - currentPositions[j][0]);
					msv[1] = learningRate[j]*(vecpos[1] - currentPositions[j][1]);
					meanshiftVectors[j] = msv;*/
					meanshiftVectors[j] = [vecpos[0] - currentPositions[j][0], vecpos[1] - currentPositions[j][1]];
					
					//if (isNaN(msv[0]) || isNaN(msv[1])) debugger;
					
					//for debugging
					//debugMVs[j] = debugMatrixMV;
					//
				}
				
				// draw meanshiftVector
				/*sketchCC.clearRect(0, 0, sketchW, sketchH);
				var nuWeights;
				for (var npidx = 0;npidx < numPatches;npidx++) {
					nuWeights = debugMVs[npidx].map(function(x) {return x*255*500;});
					drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[npidx][0]-((searchWindow-1)/2), patchPositions[npidx][1]-((searchWindow-1)/2));
				}*/
				
				var meanShiftVector = numeric.rep([numPatches*2, 1],0.0);
				for (var k = 0;k < numPatches;k++) {
					meanShiftVector[k*2][0] = meanshiftVectors[k][0];
					meanShiftVector[(k*2)+1][0] = meanshiftVectors[k][1];
				}
				
				// compute pdm parameter update
				//var prior = numeric.mul(gaussianPD, PDMVariance);
				var prior = numeric.mul(gaussianPD, varianceSeq[i]);
				if (params.weightPoints) {
					var jtj = numeric.dot(numeric.transpose(jac), numeric.dot(pointWeights, jac));
				} else {
					var jtj = numeric.dot(numeric.transpose(jac), jac);
				}
				var cpMatrix = numeric.rep([numParameters+4, 1],0.0);
				for (var l = 0;l < (numParameters+4);l++) {
					cpMatrix[l][0] = currentParameters[l];
				}
				var priorP = numeric.dot(prior, cpMatrix);
				if (params.weightPoints) {
					var jtv = numeric.dot(numeric.transpose(jac), numeric.dot(pointWeights, meanShiftVector));
				} else {
					var jtv = numeric.dot(numeric.transpose(jac), meanShiftVector);
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
				//console.log("positionnorm:"+positionNorm);
				
				// if norm < limit, then break
				if (positionNorm < convergenceLimit) {
					break;
				}
			
			}
			
			if (params.constantVelocity) {
				// add current parameter to array of previous parameters
				previousParameters.push(currentParameters.slice());
				previousParameters.splice(0, previousParameters.length == 3 ? 1 : 0);
			}
			
			// store positions, for checking convergence
			previousPositions.splice(0, previousPositions.length == 10 ? 1 : 0);
			previousPositions.push(currentPositions.slice(0));
			
			// send an event on each iteration
			var evt = document.createEvent("Event");
			evt.initEvent("clmtrackrIteration", true, true);
			document.dispatchEvent(evt)
			
			if (this.getConvergence() < 0.5) {
				if (params.stopOnConvergence) {
					this.stop();
				}
				var evt = document.createEvent("Event");
				evt.initEvent("clmtrackrConverged", true, true);
				document.dispatchEvent(evt)
			}
			
			// return new points
			return currentPositions;
		}

		/*
		 *	reset tracking, so that track() will start a new detection
		 */
		this.reset = function() {
			first = true;
			scoringHistory = [];
			for (var i = 0;i < currentParameters.length;i++) {
				currentParameters[i] = 0;
				previousParameters = [];
			}
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
			cc.fillStyle = "rgb(200,200,200)";
			cc.strokeStyle = "rgb(130,255,50)";
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
			msavg /= previousPositions.length
			return msavg;
		}

		var runnerFunction = function() {
			runnerTimeout = requestAnimFrame(runnerFunction);
			// schedule as many iterations as we can during each request
			var startTime = (new Date()).getTime();
			while (((new Date()).getTime() - startTime) < 16) {
				var tracking = this.track(runnerElement, runnerBox);
				if (!tracking) continue;
			}
		}.bind(this);

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
		
		// detect position of face on canvas/video element
		var detectPosition = function(el) {
			var canvas = document.createElement('canvas');
			canvas.width = el.width;
			canvas.height = el.height;
			var cc = canvas.getContext('2d');
			cc.drawImage(el, 0, 0, el.width, el.height);
			
			// do viola-jones on canvas to get initial guess, if we don't have any points
			var comp = ccv.detect_objects(
				ccv.grayscale(canvas), ccv.cascade, 5, 1
			);
			
			//var jf = new jsfeat_face(canvas);
			//var comp = jf.findFace();
			
			if (comp.length > 0) {
				candidate = comp[0];
			} else {
				return false;
			}
			
			for (var i = 1; i < comp.length; i++) {
				if (comp[i].confidence > candidate.confidence) {
					candidate = comp[i];
				}
			}
			
			return candidate;
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
					
					//for debugging
					//vecmatrix[k*responseWidth + l] = vecsum;
					
					vecpos[0] += vecsum*updatePosition[0];
					vecpos[1] += vecsum*updatePosition[1];
					pos_idx++;
				}
			}
			// for debugging
			//return vecmatrix;
		}
		
		// calculate score of current fit
		var checkTracking = function() {			
			scoringContext.drawImage(sketchCanvas, Math.round(msxmin), Math.round(msymin), Math.round(msmodelwidth), Math.round(msmodelheight), 0, 0, 20, 22);
			// getImageData of canvas
			var imgData = scoringContext.getImageData(0,0,20,22);
			// convert data to grayscale
			var scoringData = new Array(20*22);
			var scdata = imgData.data;
			var scmax = 0;
			var scmin = 255;
			for (var i = 0;i < 20*22;i++) {
			  scoringData[i] = scdata[i*4]*0.3 + scdata[(i*4)+1]*0.59 + scdata[(i*4)+2]*0.11;
			  if (scoringData[i] > scmax) scmax = scoringData[i];
			  if (scoringData[i] < scmin) scmin = scoringData[i];
			}
			
			if (scmax > 0) {
				// normalize & multiply by svmFilter
				var score = 0;
				var newim = scoringContext.createImageData(20,22);
				var newimdata = newim.data;
				var scdaata = new Array(20*22);
				var scdamin = 10000000;
				var scdamax = -10000000;
				for (var i = 0;i < 20*22;i++) {
					scoringData[i] = (scoringData[i]-scmin)/(scmax-scmin);
					scdaata[i] = scoringData[i]*scoringWeights[i];
					if (scdaata[scdamin] < i) scdamin = scdaata[i];
					if (scdaata[i] > scdamax) scdamax = scdaata[i];
				}

				for (var i = 0;i < 20*22;i++) {
					var nid = 255*((scdaata[i]-scdamin)/(scdamax-scdamin));
					newimdata[i*4] = nid;
					newimdata[(i*4)+1] = nid;
					newimdata[(i*4)+2] = nid;
					newimdata[(i*4)+3] = 255;
				}
				scoringContext.putImageData(newim,0,0);

				for (var i = 0;i < 20*22;i++) {
					score += (scoringData[i])*-scoringWeights[i];
				}
				score += scoringBias;

				scoringHistory.splice(0, scoringHistory.length == 5 ? 1 : 0);
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
		
		// get initial starting point for model
		var getInitialPosition = function(element, box) {
			var translateX, translateY, scaling, rotation;
			if (box) {
				candidate = {x : box[0], y : box[1], width : box[2], height : box[3]};
			} else {
				var det = detectPosition(element);
				if (!det) {
					// if no face found, stop.
					return false;
				}
			}
			
			if (model.hints && mosseFilter && left_eye_filter && right_eye_filter && nose_filter) {
				var noseFilterWidth = candidate.width * 4.5/10;
				var eyeFilterWidth = candidate.width * 6/10;
				
				// detect position of eyes and nose via mosse filter
				//
				/*element.pause();
				
				var canvasContext = document.getElementById('overlay2').getContext('2d')
				canvasContext.clearRect(0,0,500,375);
				canvasContext.strokeRect(candidate.x, candidate.y, candidate.width, candidate.height);*/
				//
				
				var nose_result = mossef_nose.track(element, Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(5/8)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth, false);
				var right_result = mossef_righteye.track(element, Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth, false);
				var left_result = mossef_lefteye.track(element, Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth, false);
				right_eye_position[0] = Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2))+right_result[0];
				right_eye_position[1] = Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2))+right_result[1];
				left_eye_position[0] = Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2))+left_result[0];
				left_eye_position[1] = Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2))+left_result[1];
				nose_position[0] = Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2))+nose_result[0];
				nose_position[1] = Math.round(candidate.y+candidate.height*(5/8)-(noseFilterWidth/2))+nose_result[1];
				
				//
				/*canvasContext.strokeRect(Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
				canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
				//canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(3/4)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth);
				canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(5/8)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth);
				
				canvasContext.fillStyle = "rgb(0,0,250)";
				canvasContext.beginPath();
				canvasContext.arc(left_eye_position[0], left_eye_position[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();
				
				canvasContext.beginPath();
				canvasContext.arc(right_eye_position[0], right_eye_position[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();
				
				canvasContext.beginPath();
				canvasContext.arc(nose_position[0], nose_position[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();
				
				debugger;
				element.play()
				canvasContext.clearRect(0,0,element.width,element.height);*/
				//
				
				// get eye and nose positions of model
				var lep = model.hints.leftEye;
				var rep = model.hints.rightEye;
				var mep = model.hints.nose;
				
				// get scaling, rotation, etc. via procrustes analysis
				var procrustes_params = procrustes([left_eye_position, right_eye_position, nose_position], [lep, rep, mep]);
				translateX = procrustes_params[0];
				translateY = procrustes_params[1];
				scaling = procrustes_params[2];
				rotation = procrustes_params[3];
				
				//element.play();
				
				//var maxscale = 1.10;
				//if ((scaling*modelHeight)/candidate.height < maxscale*0.7) scaling = (maxscale*0.7*candidate.height)/modelHeight;
				//if ((scaling*modelHeight)/candidate.height > maxscale*1.2) scaling = (maxscale*1.2*candidate.height)/modelHeight;
				
				/*var smean = [0,0];
				smean[0] += lep[0];
				smean[1] += lep[1];
				smean[0] += rep[0];
				smean[1] += rep[1];
				smean[0] += mep[0];
				smean[1] += mep[1];
				smean[0] /= 3;
				smean[1] /= 3;
				
				var nulep = [(lep[0]*scaling*Math.cos(-rotation)+lep[1]*scaling*Math.sin(-rotation))+translateX, (lep[0]*scaling*(-Math.sin(-rotation)) + lep[1]*scaling*Math.cos(-rotation))+translateY];
				var nurep = [(rep[0]*scaling*Math.cos(-rotation)+rep[1]*scaling*Math.sin(-rotation))+translateX, (rep[0]*scaling*(-Math.sin(-rotation)) + rep[1]*scaling*Math.cos(-rotation))+translateY];
				var numep = [(mep[0]*scaling*Math.cos(-rotation)+mep[1]*scaling*Math.sin(-rotation))+translateX, (mep[0]*scaling*(-Math.sin(-rotation)) + mep[1]*scaling*Math.cos(-rotation))+translateY];
				
				canvasContext.fillStyle = "rgb(200,10,100)";
				canvasContext.beginPath();
				canvasContext.arc(nulep[0], nulep[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();
				
				canvasContext.beginPath();
				canvasContext.arc(nurep[0], nurep[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();
				
				canvasContext.beginPath();
				canvasContext.arc(numep[0], numep[1], 3, 0, Math.PI*2, true);
				canvasContext.closePath();
				canvasContext.fill();*/
				
				currentParameters[0] = (scaling*Math.cos(rotation))-1;
				currentParameters[1] = (scaling*Math.sin(rotation));
				currentParameters[2] = translateX;
				currentParameters[3] = translateY;
				
				//this.draw(document.getElementById('overlay'), currentParameters);
				
			} else {
				scaling = candidate.width/modelheight;
				//var ccc = document.getElementById('overlay').getContext('2d');
				//ccc.strokeRect(candidate.x,candidate.y,candidate.width,candidate.height);
				translateX = candidate.x-(xmin*scaling)+0.1*candidate.width;
				translateY = candidate.y-(ymin*scaling)+0.25*candidate.height;
				currentParameters[0] = scaling-1;
				currentParameters[2] = translateX;
				currentParameters[3] = translateY;
			}
		
			currentPositions = calculatePositions(currentParameters, true);
			
			return [scaling, rotation, translateX, translateY];
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
		
		// procrustes analysis
		function procrustes(template, shape) {
			// assume template and shape is a vector of x,y-coordinates
			//i.e. template = [[x1,y1], [x2,y2], [x3,y3]];
			var templateClone = [];
			var shapeClone = [];
			for (var i = 0;i < template.length;i++) {
				templateClone[i] = [template[i][0], template[i][1]];
			}
			for (var i = 0;i < shape.length;i++) {
				shapeClone[i] = [shape[i][0], shape[i][1]];
			}
			shape = shapeClone;
			template = templateClone;
			
			// calculate translation
			var templateMean = [0.0, 0.0];
			for (var i = 0;i < template.length;i++) {
				templateMean[0] += template[i][0];
				templateMean[1] += template[i][1];
			}
			templateMean[0] /= template.length;
			templateMean[1] /= template.length;
			
			var shapeMean = [0.0, 0.0];
			for (var i = 0;i < shape.length;i++) {
				shapeMean[0] += shape[i][0];
				shapeMean[1] += shape[i][1];
			}
			shapeMean[0] /= shape.length;
			shapeMean[1] /= shape.length;
			
			var translationX = templateMean[0] - shapeMean[0];
			var translationY = templateMean[1] - shapeMean[1];
			
			// centralize
			for (var i = 0;i < shape.length;i++) {
				shape[i][0] -= shapeMean[0];
				shape[i][1] -= shapeMean[1];
			}
			for (var i = 0;i < template.length;i++) {
				template[i][0] -= templateMean[0];
				template[i][1] -= templateMean[1];
			}
			
			// scaling
			
			var scaleS = 0.0;
			for (var i = 0;i < shape.length;i++) {
				scaleS += ((shape[i][0])*(shape[i][0]));
				scaleS += ((shape[i][1])*(shape[i][1]));
			}
			scaleS = Math.sqrt(scaleS/shape.length);
			
			var scaleT = 0.0;
			for (var i = 0;i < template.length;i++) {
				scaleT += ((template[i][0])*(template[i][0]));
				scaleT += ((template[i][1])*(template[i][1]));
			}
			scaleT = Math.sqrt(scaleT/template.length);
			
			var scaling = scaleT/scaleS;
			
			for (var i = 0;i < shape.length;i++) {
				shape[i][0] *= scaling;
				shape[i][1] *= scaling;
			}
				
			// rotation
			
			var top = 0.0;
			var bottom = 0.0;
			for (var i = 0;i < shape.length;i++) {
				top += (shape[i][0]*template[i][1] - shape[i][1]*template[i][0]);
				bottom += (shape[i][0]*template[i][0] + shape[i][1]*template[i][1]);
			}
			var rotation = Math.atan(top/bottom);
			
			translationX += (shapeMean[0]-(scaling*Math.cos(-rotation)*shapeMean[0])-(scaling*shapeMean[1]*Math.sin(-rotation)));
			translationY += (shapeMean[1]+(scaling*Math.sin(-rotation)*shapeMean[0])-(scaling*shapeMean[1]*Math.cos(-rotation)));
			
			//returns rotation, scaling, transformx and transformx
			return [translationX, translationY, scaling, rotation];
		}
		
		// function to draw pixeldata on some canvas, only used for debugging
		var drawData = function(canvasContext, data, width, height, transposed, drawX, drawY) {
			var psci = canvasContext.createImageData(width, height);
			var pscidata = psci.data;
			for (var j = 0;j < width*height;j++) {
				if (!transposed) {
					var val = data[(j%width)+((j/width) >> 0)*width];
				} else {
					var val = data[(j%height)*height+((j/height) >> 0)];
				}
				val = val > 255 ? 255 : val;
				val = val < 0 ? 0 : val;
				pscidata[j*4] = val;
				pscidata[(j*4)+1] = val;
				pscidata[(j*4)+2] = val;
				pscidata[(j*4)+3] = 255;
			}
			canvasContext.putImageData(psci, drawX, drawY);
		}
		
		var requestAnimFrame = (function() {
			return window.requestAnimationFrame ||
			window.webkitRequestAnimationFrame ||
			window.mozRequestAnimationFrame ||
			window.oRequestAnimationFrame ||
			window.msRequestAnimationFrame ||
			function(/* function FrameRequestCallback */ callback, /* DOMElement Element */ element) {
				return window.setTimeout(callback, 1000/60);
			};
		})();
		
		var cancelRequestAnimFrame = (function() {
			return window.cancelCancelRequestAnimationFrame ||
				window.webkitCancelRequestAnimationFrame ||
				window.mozCancelRequestAnimationFrame ||
				window.oCancelRequestAnimationFrame ||
				window.msCancelRequestAnimationFrame ||
				window.clearTimeout;
		})();
		
		return true;
	}
}
"use strict";

var webglFilter = function() {
  var gl, canvas;
  var filterWidth, filterHeight, patchWidth, patchHeight, numPatches, canvasWidth, canvasHeight;
  var patchResponseProgram, patchDrawProgram;
  var fbo, numBlocks, patchTex;
  var first = true;
  var drawRectBuffer, drawLayerBuffer, drawImageBuffer, rttTexture, filters;
  var texCoordBuffer, texCoordLocation, apositionBuffer;
  var newCanvasWidth, newCanvasBlockHeight, newCanvasHeight;
  var drawOutRectangles, drawOutImages, drawOutLayer;
  var patchCells, textureWidth, textureHeight, patchSize, patchArray;
  
  this.init = function(filterVector, nP, pW, pH, fW, fH, drawOut) {
    // we assume filterVector goes from left to right, rowwise, i.e. row-major order
    
    if (fW != fH) {
      alert("filter width and height must be same size!");
      return;
    }
    
    // if filter width is not odd, alert
    if (fW % 2 == 0 || fH % 2 == 0) {
      alert("filters used in svm must be of odd dimensions!")
      return;
    }
    
    filterWidth = fW;
    filterHeight = fH;
    patchWidth = pW;
    patchHeight = pH;
    numPatches = nP;

    patchResponseFS = [
      "precision mediump float;",
      "",
      "uniform vec2 u_onePixelFilters;",
      "uniform vec2 u_onePixelPatches;",
      "const float u_filterwidth = "+filterWidth.toFixed(1)+";",
      "const float u_filterheight = "+filterHeight.toFixed(1)+";",
      "const float u_halffilterwidth = "+((filterWidth-1.0)/2).toFixed(1)+";",
      "const float u_halffilterheight = "+((filterHeight-1.0)/2).toFixed(1)+";",
      "const float u_filtersize = "+(filterWidth*filterHeight)+".0;",
      "const float u_divfiltersize = "+(1/(filterWidth*filterHeight)).toFixed(10)+";",
      "",
      "// our patches",
      "uniform sampler2D u_patches;",
      "// our filters",
      "uniform sampler2D u_filters;",
      "",
      "// the texCoords passed in from the vertex shader.",
      "varying vec2 v_texCoord;",
      "varying vec2 v_texCoordFilters; // this should give us correct filter",
      "",
      "void main() {",
      "  vec4 colorSum = vec4(0.0, 0.0, 0.0, 0.0);",
      "  vec4 maxn = vec4(0.0, 0.0, 0.0, 0.0);",
      "  vec4 minn = vec4(256.0, 256.0, 256.0, 256.0);",
      "  vec4 scale = vec4(0.0, 0.0, 0.0, 0.0);",
      "  vec4 value = vec4(0.0, 0.0, 0.0, 0.0);",
      "  for (float w = 0.0;w < u_filterwidth;w++) {",
      "    for (float h = 0.0;h < u_filterheight;h++) {",
      "      value = texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(w-u_halffilterwidth, h-u_halffilterheight));",
      "      maxn = max(value, maxn);",
      "      minn = min(value, minn);",
      "    } ",
      "  }",
      "  scale = maxn-minn;",
      "  for (float w = 0.0;w < u_filterwidth;w++) {",
      "    for (float h = 0.0;h < u_filterheight;h++) {",
      "      colorSum += ",
      "        ((texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(w-u_halffilterwidth, h-u_halffilterheight))-minn)/(scale)) * ",
      "        texture2D(u_filters, v_texCoordFilters + u_onePixelFilters * vec2(w-u_halffilterwidth, h-u_halffilterheight)); ",
      "    } ",
      "  }",
      "  // logistic transformation",
      "  colorSum = 1.0/(1.0 + exp(- (colorSum-1.0) ));",
      "  gl_FragColor = colorSum;",
      "}"
    ].join('\n');
    
    numBlocks = Math.floor(numPatches / 4) + Math.ceil(numPatches % 4);
    canvasWidth = patchWidth;
    canvasHeight = patchHeight*numBlocks;
    newCanvasWidth = patchWidth-filterWidth+1;
    newCanvasBlockHeight = patchHeight-filterWidth+1;
    newCanvasHeight = newCanvasBlockHeight*numPatches;
    
    //create canvas
    canvas = document.createElement('canvas')
    canvas.setAttribute('width', (patchWidth-filterWidth+1)+"px");
    canvas.setAttribute('height', ((patchHeight-filterHeight+1)*numPatches)+"px");
    canvas.setAttribute('id', 'renderCanvas');
    canvas.setAttribute('style', 'display:none;');
    document.body.appendChild(canvas);
    
    // TODO : isolate this library from webgl-util.js
    gl = setupWebGL(canvas, {premultipliedAlpha: false, preserveDrawingBuffer : true});
    
    // check for float textures support and fail if not
    if (!gl.getExtension("OES_texture_float")) {
      alert("Your graphics card does not support floating point textures! :(");
      return;
    }
    
    // TODO : alternatively set up fallback to packing and unpacking floats:
    //   http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work
    
    // calculate position of vertex rectangles to draw out
    var rectangles = [];
    var halfFilter = (filterWidth-1)/2;
    var yOffset;
    
    for (var i = 0;i < numBlocks;i++) {
      yOffset = i*patchHeight;
      //first triangle
      rectangles = rectangles.concat(
        [halfFilter, yOffset+halfFilter, 
        patchWidth-halfFilter, yOffset+halfFilter,
        halfFilter, yOffset+patchHeight-halfFilter]
      );
      //second triangle
      rectangles = rectangles.concat(
        [halfFilter, yOffset+patchHeight-halfFilter, 
        patchWidth-halfFilter, yOffset+halfFilter,
        patchWidth-halfFilter, yOffset+patchHeight-halfFilter]
      );
    }
    rectangles = new Float32Array(rectangles);
    
    // calculate position of image rectangles to draw out
    var irectangles = [];
    for (var i = 0;i < rectangles.length;i++) {
      if (i % 2 == 0) {
        irectangles[i] = rectangles[i]/canvasWidth;
      } else {
        irectangles[i] = rectangles[i]/canvasHeight;
      }
    }
    irectangles = new Float32Array(irectangles);
    
    // setup patchdraw program
    var drVertexShader = loadShader(gl, drawResponsesVS, gl.VERTEX_SHADER);
    var drFragmentShader = loadShader(gl, drawResponsesFS, gl.FRAGMENT_SHADER);
    patchDrawProgram = createProgram(gl, [drVertexShader, drFragmentShader]);
    gl.useProgram(patchDrawProgram);
    
    // set the resolution/dimension of the canvas
    var resolutionLocation = gl.getUniformLocation(patchDrawProgram, "u_resolutiondraw");
    gl.uniform2f(resolutionLocation, newCanvasWidth, newCanvasHeight);
    
    // set u_responses
    var responsesLocation = gl.getUniformLocation(patchDrawProgram, "u_responses");
    gl.uniform1i(responsesLocation, 2);
    
    
    
    
    // setup patchresponse program
    var prVertexShader = loadShader(gl, patchResponseVS, gl.VERTEX_SHADER);
    var prFragmentShader = loadShader(gl, patchResponseFS, gl.FRAGMENT_SHADER);
    patchResponseProgram = createProgram(gl, [prVertexShader, prFragmentShader]);
    gl.useProgram(patchResponseProgram);
    
    // set the resolution/dimension of the canvas
    var resolutionLocation = gl.getUniformLocation(patchResponseProgram, "u_resolution");
    gl.uniform2f(resolutionLocation, canvasWidth, canvasHeight);
    
    // set the patchHeight
    var patchHeightLocation = gl.getUniformLocation(patchResponseProgram, "u_patchHeight");
    gl.uniform1f(patchHeightLocation, 1/numBlocks);
    
    // set the filterHeight
    var filterHeightLocation = gl.getUniformLocation(patchResponseProgram, "u_filterHeight");
    gl.uniform1f(filterHeightLocation, 1/numBlocks);
    
    // set the midpoint
    var midpointLocation = gl.getUniformLocation(patchResponseProgram, "u_midpoint");
    gl.uniform2f(midpointLocation, 0.5, (1/(numBlocks*2.0)) );
    
    // set the onepixel size for patches
    var onePixelPatchLocation = gl.getUniformLocation(patchResponseProgram, "u_onePixelPatches");
    gl.uniform2f(onePixelPatchLocation, 1/patchWidth, 1/(patchHeight*numBlocks));
    
    // set the onepixel size for filters
    var onePixelFilterLocation = gl.getUniformLocation(patchResponseProgram, "u_onePixelFilters");
    gl.uniform2f(onePixelFilterLocation, 1/filterWidth, 1/(filterHeight*numBlocks));
    
    // set up vertices with rectangles
    var positionLocation = gl.getAttribLocation(patchResponseProgram, "a_position");
    apositionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
    gl.bufferData(
      gl.ARRAY_BUFFER, 
      rectangles, 
      gl.STATIC_DRAW);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
    
    // set up texture positions
    texCoordLocation = gl.getAttribLocation(patchResponseProgram, "a_texCoord");
    // provide texture coordinates for the rectangle.
    texCoordBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, irectangles, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(texCoordLocation);
    gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);
    
    // insert filters into float32array and pass to texture via this mechanism:
    // http://stackoverflow.com/questions/7709689/webgl-pass-array-shader
    var filterSize = filterWidth*filterHeight;
    var filterArray = new Float32Array(filterSize*(numBlocks)*4);
    for (var i = 0;i < numBlocks;i++) {
      for (var j = 0;j < filterHeight;j++) {
        for (var k = 0;k < filterWidth;k++) {
          //set r with first filter
          if (i*4 < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = filterVector[i*4][(j*filterWidth) + k];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = 0;
          }
          //set g with 2nd filter
          if ((i*4 + 1) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = filterVector[(i*4)+1][(j*filterWidth) + k];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = 0;
          }
          //set b with 3rd filter
          if ((i*4 + 2) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = filterVector[(i*4)+2][(j*filterWidth) + k];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = 0;
          }
          //set a with 4th filter
          if ((i*4 + 3) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 3] = filterVector[(i*4)+3][(j*filterWidth) + k];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 3] = 0;
          }
        }
      }
    }
    
    gl.activeTexture(gl.TEXTURE0);
    filters = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, filters);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, filterWidth, filterHeight*numBlocks, 0, gl.RGBA, gl.FLOAT, filterArray);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.uniform1i(gl.getUniformLocation(patchResponseProgram, "u_filters"), 0);
    
    // set up buffer to draw to
    gl.activeTexture(gl.TEXTURE2);
    rttTexture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, rttTexture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, patchWidth, patchHeight*numBlocks, 0, gl.RGBA, gl.FLOAT, null);
    
    // set this buffer as framebuffer
    fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, rttTexture, 0);
    
    gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);
    
    // preinitialization of drawOut variables
    
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
    patchCells = (Math.floor(numPatches / 4) + Math.ceil(numPatches % 4));
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
    
    // initialize patchvariables
    textureWidth = patchWidth;
    textureHeight = patchHeight*patchCells;
    patchSize = patchWidth*patchHeight;
    patchArray = new Float32Array(patchSize*patchCells*4);
    
  }

  this.getResponses = function(patches, drawOut) {
    // TODO: check patches correct length/dimension
    
    // switch to response generation program if we're not already using it
    if (!first) {
      gl.useProgram(patchResponseProgram);
      
      // set up vertices with rectangles
      var positionLocation = gl.getAttribLocation(patchResponseProgram, "a_position");
      gl.bindBuffer(gl.ARRAY_BUFFER, apositionBuffer);
      gl.enableVertexAttribArray(positionLocation);
      gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
      
      // set up texture positions
      var texCoordLocation = gl.getAttribLocation(patchResponseProgram, "a_texCoord");
      // provide texture coordinates for the rectangle.
      gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
      gl.enableVertexAttribArray(texCoordLocation);
      gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);
      
      // set framebuffer to the original one if not already using it
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      
      gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);
    }
    
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
    if (first) {
      patchTex = gl.createTexture();
    }
    gl.bindTexture(gl.TEXTURE_2D, patchTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, textureWidth, textureHeight, 0, gl.RGBA, gl.FLOAT, patchArray);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.uniform1i(gl.getUniformLocation(patchResponseProgram, "u_patches"), 1);
    
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER)
    
    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, patchCells*6);
    
    //gl.finish();
    
    if (drawOut) {
      
      // switch programs
      gl.useProgram(patchDrawProgram);
      
      // bind canvas buffer
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, newCanvasWidth, newCanvasHeight);
      
      gl.clearColor(0.0, 0.0, 0.0, 1.0);
      gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER)
      
      if (first) {
        drawRectBuffer = gl.createBuffer();
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, drawRectBuffer);
      gl.bufferData(
        gl.ARRAY_BUFFER, 
        drawOutRectangles, 
        gl.STATIC_DRAW);
      var positionLocation = gl.getAttribLocation(patchDrawProgram, "a_position_draw");
      gl.enableVertexAttribArray(positionLocation);
      gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
      
      if (first) {
        drawImageBuffer = gl.createBuffer();
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, drawImageBuffer);
      gl.bufferData(
        gl.ARRAY_BUFFER, 
        drawOutImages, 
        gl.STATIC_DRAW);
      var textureLocation = gl.getAttribLocation(patchDrawProgram, "a_texCoord_draw");
      gl.enableVertexAttribArray(textureLocation);
      gl.vertexAttribPointer(textureLocation, 2, gl.FLOAT, false, 0, 0);
      
      if (first) {
        drawLayerBuffer = gl.createBuffer();
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, drawLayerBuffer);
      gl.bufferData(
        gl.ARRAY_BUFFER, 
        drawOutLayer, 
        gl.STATIC_DRAW);
      var layerLocation = gl.getAttribLocation(patchDrawProgram, "a_patchChoice_draw");
      gl.enableVertexAttribArray(layerLocation);
      gl.vertexAttribPointer(layerLocation, 1, gl.FLOAT, false, 0, 0);
      
      // draw out
      gl.drawArrays(gl.TRIANGLES, 0, numPatches*6);
    }

    var responses = getOutput();
    responses = unpackToFloat(responses);
    
    // split
    responses = splitArray(responses, numPatches);
    
    // normalize responses to lie within [0,1]
    var rl = responses.length;
    
    for (var i = 0;i < rl;i++) {
      responses[i] = normalizeFilterMatrix(responses[i]);
    }
    
    first = false;
    
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
    var data = gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, pixelValues);
    // return
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
      console.log("a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.")
      response = response.map(function() {return 1});
    } else {
      for (var i = 0;i < msize;i++) {
        response[i] = (response[i]-min)/dist;
      }
    }
    
    return response
  }
  
  var patchResponseVS = [
    "attribute vec2 a_texCoord;",
    "attribute vec2 a_position;",
    "",
    "uniform vec2 u_resolution;",
    "uniform float u_patchHeight;",
    "uniform float u_filterHeight;",
    "uniform vec2 u_midpoint;",
    "",
    "varying vec2 v_texCoord;",
    "varying vec2 v_texCoordFilters;",
    "",
    "void main() {",
    "   // convert the rectangle from pixels to 0.0 to 1.0",
    "   vec2 zeroToOne = a_position / u_resolution;",
    "",
    "   // convert from 0->1 to 0->2",
    "   vec2 zeroToTwo = zeroToOne * 2.0;",
    "",
    "   // convert from 0->2 to -1->+1 (clipspace)",
    "   vec2 clipSpace = zeroToTwo - 1.0;",
    "   ",
    "   // transform coordinates to regular coordinates",
    "   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);",
    " ",
    "   // pass the texCoord to the fragment shader",
    "   v_texCoord = a_texCoord;",
    "   ",
    "   // set the filtertexture coordinate based on number filter to use",
    "   v_texCoordFilters = u_midpoint + vec2(0.0, u_filterHeight * floor(a_texCoord[1]/u_patchHeight));",
    "}"
  ].join('\n');
  
  var patchResponseFS;
  
  var drawResponsesVS = [
    "attribute vec2 a_texCoord_draw;",
    "attribute vec2 a_position_draw;",
    "attribute float a_patchChoice_draw;",
    "",
    "uniform vec2 u_resolutiondraw;",
    "",
    "varying vec2 v_texCoord;",
    "varying float v_select;",
    "",
    "void main() {",
    "   // convert the rectangle from pixels to 0.0 to 1.0",
    "   vec2 zeroToOne = a_position_draw / u_resolutiondraw;",
    "",
    "   // convert from 0->1 to 0->2",
    "   vec2 zeroToTwo = zeroToOne * 2.0;",
    "",
    "   // convert from 0->2 to -1->+1 (clipspace)",
    "   vec2 clipSpace = zeroToTwo - 1.0;",
    "   ",
    "   // transform coordinates to regular coordinates",
    "   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);",
    "",
    "   // pass the texCoord to the fragment shader",
    "   v_texCoord = a_texCoord_draw;",
    "   ",
    "   v_select = a_patchChoice_draw;",
    "}"
  ].join('\n');
  
  var drawResponsesFS = [
    "precision mediump float;",
    "",
    "// our responses",
    "uniform sampler2D u_responses;",
    "",
    "// the texCoords passed in from the vertex shader.",
    "varying vec2 v_texCoord;",
    "varying float v_select;",
    "",
    "const vec4 bit_shift = vec4(256.0*256.0*256.0, 256.0*256.0, 256.0, 1.0);",
    "const vec4 bit_mask  = vec4(0.0, 1.0/256.0, 1.0/256.0, 1.0/256.0);",
    "",
    "// packing code from here http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work",
    "void main() {",
    "  vec4 colorSum = texture2D(u_responses, v_texCoord);",
    "  float value = 0.0;",
    "  if (v_select < 0.1) {",
    "    value = colorSum[0];",
    "  } else if (v_select > 0.9 && v_select < 1.1) {",
    "    value = colorSum[1];",
    "  } else if (v_select > 1.9 && v_select < 2.1) {",
    "    value = colorSum[2];",
    "  } else if (v_select > 2.9 && v_select < 3.1) {",
    "    value = colorSum[3];",
    "  } else {",
    "    value = 1.0;",
    "  }",
    "  ",
    "  vec4 res = fract(value * bit_shift);",
    "  res -= res.xxyz * bit_mask;",
    "  ",
    "  //gl_FragColor = vec4(value, value, value, value);",
    "  //gl_FragColor = vec4(1.0, value, 1.0, 1.0);",
    "  gl_FragColor = res;",
    "}"
  ].join('\n');
};

// The rest of the code is based on webgl-utils.js authored by Gregg Tavares, license below:
/*
 * Copyright (c) 2011, Gregg Tavares
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 *  * Neither the name of greggman.com nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

(function() {
  
  /**
   * Wrapped logging function.
   * @param {string} msg The message to log.
   */
  var log = function(msg) {
    if (window.console && window.console.log) {
      window.console.log(msg);
    }
  };

  /**
   * Wrapped logging function.
   * @param {string} msg The message to log.
   */
  var error = function(msg) {
    if (window.console) {
      if (window.console.error) {
        window.console.error(msg);
      }
      else if (window.console.log) {
        window.console.log(msg);
      }
    }
  };

  /**
   * Turn off all logging.
   */
  var loggingOff = function() {
    log = function() {};
    error = function() {};
  };

  /**
   * Check if the page is embedded.
   * @return {boolean} True of we are in an iframe
   */
  var isInIFrame = function() {
    return window != window.top;
  };

  /**
   * Converts a WebGL enum to a string
   * @param {!WebGLContext} gl The WebGLContext to use.
   * @param {number} value The enum value.
   * @return {string} The enum as a string.
   */
  var glEnumToString = function(gl, value) {
    for (var p in gl) {
      if (gl[p] == value) {
        return p;
      }
    }
    return "0x" + value.toString(16);
  };

  /**
   * Creates the HTLM for a failure message
   * @param {string} canvasContainerId id of container of th
   *        canvas.
   * @return {string} The html.
   */
  var makeFailHTML = function(msg) {
    return '' +
      '<table style="background-color: #8CE; width: 100%; height: 100%;"><tr>' +
      '<td align="center">' +
      '<div style="display: table-cell; vertical-align: middle;">' +
      '<div style="">' + msg + '</div>' +
      '</div>' +
      '</td></tr></table>';
  };

  /**
   * Mesasge for getting a webgl browser
   * @type {string}
   */
  var GET_A_WEBGL_BROWSER = '' +
    'This page requires a browser that supports WebGL.<br/>' +
    '<a href="http://get.webgl.org">Click here to upgrade your browser.</a>';

  /**
   * Mesasge for need better hardware
   * @type {string}
   */
  var OTHER_PROBLEM = '' +
    "It doesn't appear your computer can support WebGL.<br/>" +
    '<a href="http://get.webgl.org/troubleshooting/">Click here for more information.</a>';

  /**
   * Creates a webgl context. If creation fails it will
   * change the contents of the container of the <canvas>
   * tag to an error message with the correct links for WebGL.
   * @param {Element} canvas. The canvas element to create a
   *     context from.
   * @param {WebGLContextCreationAttirbutes} opt_attribs Any
   *     creation attributes you want to pass in.
   * @return {WebGLRenderingContext} The created context.
   */
  var setupWebGL = function(canvas, opt_attribs) {
    function showLink(str) {
      var container = canvas.parentNode;
      if (container) {
        container.innerHTML = makeFailHTML(str);
      }
    };

    if (!window.WebGLRenderingContext) {
      //showLink(GET_A_WEBGL_BROWSER);
      return null;
    }

    var context = create3DContext(canvas, opt_attribs);
    if (!context) {
      //showLink(OTHER_PROBLEM);
      return null;
    }
    return context;
  };

  /**
   * Creates a webgl context.
   * @param {!Canvas} canvas The canvas tag to get context
   *     from. If one is not passed in one will be created.
   * @return {!WebGLContext} The created context.
   */
  var create3DContext = function(canvas, opt_attribs) {
    var names = ["webgl", "experimental-webgl"];
    var context = null;
    for (var ii = 0; ii < names.length; ++ii) {
      try {
        context = canvas.getContext(names[ii], opt_attribs);
      } catch(e) {}
      if (context) {
        break;
      }
    }
    return context;
  }

  var updateCSSIfInIFrame = function() {
    if (isInIFrame()) {
      document.body.className = "iframe";
    }
  };

  /**
   * Gets a WebGL context.
   * makes its backing store the size it is displayed.
   */
  var getWebGLContext = function(canvas) {
    if (isInIFrame()) {
      updateCSSIfInIFrame();

      // make the canvas backing store the size it's displayed.
      canvas.width = canvas.clientWidth;
      canvas.height = canvas.clientHeight;
    }

    var gl = setupWebGL(canvas);
    return gl;
  };

  /**
   * Loads a shader.
   * @param {!WebGLContext} gl The WebGLContext to use.
   * @param {string} shaderSource The shader source.
   * @param {number} shaderType The type of shader.
   * @param {function(string): void) opt_errorCallback callback for errors.
   * @return {!WebGLShader} The created shader.
   */
  var loadShader = function(gl, shaderSource, shaderType, opt_errorCallback) {
    var errFn = opt_errorCallback || error;
    // Create the shader object
    var shader = gl.createShader(shaderType);

    // Load the shader source
    gl.shaderSource(shader, shaderSource);

    // Compile the shader
    gl.compileShader(shader);

    // Check the compile status
    var compiled = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
    if (!compiled) {
      // Something went wrong during compilation; get the error
      lastError = gl.getShaderInfoLog(shader);
      errFn("*** Error compiling shader '" + shader + "':" + lastError);
      gl.deleteShader(shader);
      return null;
    }

    return shader;
  }

  /**
   * Creates a program, attaches shaders, binds attrib locations, links the
   * program and calls useProgram.
   * @param {!Array.<!WebGLShader>} shaders The shaders to attach
   * @param {!Array.<string>} opt_attribs The attribs names.
   * @param {!Array.<number>} opt_locations The locations for the attribs.
   */
  var loadProgram = function(gl, shaders, opt_attribs, opt_locations) {
    var program = gl.createProgram();
    for (var ii = 0; ii < shaders.length; ++ii) {
      gl.attachShader(program, shaders[ii]);
    }
    if (opt_attribs) {
      for (var ii = 0; ii < opt_attribs.length; ++ii) {
        gl.bindAttribLocation(
            program,
            opt_locations ? opt_locations[ii] : ii,
            opt_attribs[ii]);
      }
    }
    gl.linkProgram(program);

    // Check the link status
    var linked = gl.getProgramParameter(program, gl.LINK_STATUS);
    if (!linked) {
        // something went wrong with the link
        lastError = gl.getProgramInfoLog (program);
        error("Error in program linking:" + lastError);

        gl.deleteProgram(program);
        return null;
    }
    return program;
  };

  /**
   * Loads a shader from a script tag.
   * @param {!WebGLContext} gl The WebGLContext to use.
   * @param {string} scriptId The id of the script tag.
   * @param {number} opt_shaderType The type of shader. If not passed in it will
   *     be derived from the type of the script tag.
   * @param {function(string): void) opt_errorCallback callback for errors.
   * @return {!WebGLShader} The created shader.
   */
  var createShaderFromScript = function(
      gl, scriptId, opt_shaderType, opt_errorCallback) {
    var shaderSource = "";
    var shaderType;
    var shaderScript = document.getElementById(scriptId);
    if (!shaderScript) {
      throw("*** Error: unknown script element" + scriptId);
    }
    shaderSource = shaderScript.text;

    if (!opt_shaderType) {
      if (shaderScript.type == "x-shader/x-vertex") {
        shaderType = gl.VERTEX_SHADER;
      } else if (shaderScript.type == "x-shader/x-fragment") {
        shaderType = gl.FRAGMENT_SHADER;
      } else if (shaderType != gl.VERTEX_SHADER && shaderType != gl.FRAGMENT_SHADER) {
        throw("*** Error: unknown shader type");
        return null;
      }
    }

    return loadShader(
        gl, shaderSource, opt_shaderType ? opt_shaderType : shaderType,
        opt_errorCallback);
  };

  /* export functions */
  window.setupWebGL = setupWebGL;
  window.createProgram = loadProgram;
  window.createShaderFromScriptElement = createShaderFromScript;
  window.getWebGLContext = getWebGLContext;
  window.updateCSSIfInIFrame = updateCSSIfInIFrame;
  window.loadShader = loadShader;

}());

"use strict";

var svmFilter = function() {
  
  var _fft, fft_filters, responses;
  var fft_size, filterLength, filter_width, search_width, num_patches;
  var temp_imag_part, temp_real_part;
  
  // fft function
  this.fft_inplace = function(array, _im_part) {
      // in-place
      
      if (typeof _im_part == "undefined") {
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
  
  this.init = function(filter_input, numPatches, filterWidth, searchWidth) {
    
    var temp, fft, offset;
    
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
          
          //console.log(k + ","+ j+":" + (k+xOffset+((j+yOffset)*fft_size)))
        }
      }

      // fft it and store
      fft_filter = this.fft_inplace(flar_fi0, flar_fi1);
      fft_filters[i] = fft_filter;
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
    var response, temp, edge;
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
        value = patch[k + (patch_width*j)]
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
    return response
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
      console.log("a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.")
    } else {
      for (var i = 0;i < msize;i++) {
        response[i] = (response[i]-min)/dist;
      }
    }
  }

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
          rowIdx = x2*_n
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
      var fftlen = (_n*_n)-1;
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
}

// requires mosse.js

var mosseFilterResponses = function() {
  
  var filters = [];
  var responses = [];
  var num_Patches = 0;
  
  this.init = function(filter_input, numPatches, filterWidth, filterHeight) {
    // load filters, make fft ready
    
    for (var i = 0;i < numPatches;i++) {
      var temp = {};
      temp.width = filterWidth;
      temp.height = filterHeight;
      var filterLength = filterWidth*filterHeight
      var flar_fi0 = new Float64Array(filterLength);
      var flar_fi1 = new Float64Array(filterLength);
      for (var j = 0;j < filterLength;j++) {
        flar_fi0[j] = filter_input[i][0][j];
        flar_fi1[j] = filter_input[i][1][j];
      }
      temp.real = flar_fi0;
      temp.imag = flar_fi1;
      filters[i] = new mosseFilter();
      filters[i].load(temp);
    }
    
    num_Patches = numPatches;
  }
  
  this.getResponses = function(patches) {
    for (var i = 0;i < num_Patches;i++) {
      responses[i] = filters[i].getResponse(patches[i]);
      //responses[i] = logisticResponse(responses[i]);
      responses[i] = normalizeFilterMatrix(responses[i]);
    }
    
    return responses;
  }
  
  var logisticResponse = function(response) {
    // create probability by doing logistic transformation
    var filter_size = response.length;
    for (var j = 0;j < filter_size;j++) {
      response[j] = 1.0/(1.0 + Math.exp(- (response[j]-1.0) ));
    }
    return response;
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
      console.log("a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.")
      response = response.map(function() {return 1});
    } else {
      for (var i = 0;i < msize;i++) {
        response[i] = (response[i]-min)/dist;
      }
    }
    
    return response
  }
}

var left_eye_filter = {"real": [1.5419219943717721, 0.40010880110578706, -0.79043641265342957, -1.2685464969238938, 0.39878117336167285, -1.0673489992245377, -0.079880838229404019, -0.45374680224191505, -0.043474097938900787, -0.31125662385352687, 0.17092430376098702, -0.29613086164846153, 0.5616469648110296, -1.559786848789493, 0.64513037997492662, -1.2899747976234162, 1.1761667998175334, -1.2899747976233551, 0.64513037997490474, -1.5597868487894897, 0.56164696481102505, -0.29613086164845964, 0.17092430376099094, -0.31125662385352959, -0.043474097938900787, -0.45374680224191177, -0.079880838229404658, -1.0673489992245357, 0.39878117336167307, -1.2685464969238942, -0.79043641265343012, 0.40010880110578717, -1.3820969331049027, 0.069560471269205768, -1.9786339579213206, -1.9807415717551982, -0.78667274410450883, -1.2217002325587256, -0.19150029104902774, -0.35131617290773243, -0.17646388464205803, -0.16672095020503441, -0.092298612924566523, -0.028899376452253527, -0.1314555696102146, -0.32892265898101813, -0.40987148655061206, 0.11741827111366547, -0.67254330182605138, -0.46007833291519956, -0.67215259521101001, -0.44871907432473013, -0.034749316729184583, 0.0055639281302433969, -0.17675902360981591, -0.26196208085032191, -0.36301254306387037, -0.33546767337818123, -0.6458889740799838, -1.1981932989987978, 0.12372650763830917, -1.4996172161865935, -2.4084298023013888, -2.0505291279591722, -1.7249706159518585, -2.277646289702639, -3.1259631743419591, -2.9656385065342015, -2.8480835086962011, -1.4260964500310189, -0.61792590829173544, -0.2611655301498782, -0.38519889843539723, -0.17511899827006483, -0.32808050503227176, 0.0076800871037463036, -0.18710828510427668, 0.1976534820339281, -0.55444453100465052, 0.14583567590328381, -0.69844971117515287, -0.90188577233526623, -0.53500016384583371, -0.044420751861669799, 0.014727914354086128, -0.28084584584371913, -0.29890408748685848, -0.39431380149336548, -0.39569215798819307, -0.74351999988258299, -0.82502198370631752, -1.851491897104155, -0.74302378668934244, 0.21156442062863762, -3.3061472495599986, -1.7990472945779568, -2.2193764251732282, -2.3438802466919251, -3.3615971067123311, -3.5383249085863708, -2.2639673745086588, -2.0271757806780748, -0.75242583405872232, -0.30143411016839378, -0.3625272253546275, -0.25489431004647689, -0.18928491561467081, -0.1179891518538482, 0.027920290231533224, -0.035472107498143821, -0.29008721857562259, -0.3604588674139817, -0.39156143807433802, -0.82222257402876564, -0.44979914971695928, -0.098136330355476253, 0.065628582466229365, -0.33607304327303128, -0.32161201323497779, -0.41856090178723965, -0.64028425429629054, -0.7766428172010218, -1.3946448661671447, -2.2603422126144683, -0.38769722219534525, -0.95341593939478653, -1.412952994959813, -2.3602336858020432, -1.2756392437278019, -2.0983496132652038, -2.5682454610054268, -2.8791053946930378, -2.1809972632688095, -0.84281293847776861, -0.75998936793718697, -0.18584599820380068, -0.30105748355308259, -0.16098142942852958, -0.13792125740417191, -0.089790022871128708, -0.12321821342876504, -0.1128661923016878, -0.3924098378001975, -0.5780902167586397, -0.48685989567066695, -0.53565359443296234, -0.051036689850526382, -0.0068547033925117689, -0.18963405157839419, -0.22514761090777807, -0.35555823460888908, -0.46670603976585517, -0.56179541485257889, -0.7495095888115163, -1.4772075422260349, -1.5836466114968029, -2.3846549454186694, -1.4884613952536236, -1.8237453905245253, -1.6712324532934877, -1.5169157844507295, -1.6930052820597281, -2.1023566589276004, -2.2062031109308458, -1.7945281756942255, -0.26457398838912649, 0.22038139379151148, -0.43479836723775234, -0.19830827357221226, -0.18018565146479498, -0.097060879184795737, -0.10088329756370379, -0.063069709957272527, -0.17970932516041177, -0.1943040732581543, -0.37970560392277619, -0.47302301606251812, -0.30366967948052181, -0.064732391018915397, -0.08902516330269715, -0.082000200083027344, -0.22965854401457736, -0.32035624605031326, -0.31836783196552437, -0.40132058236311119, -0.65601747033470859, -0.59040483751417483, -1.8503084663080034, -1.8694842425148914, -1.9326778896298584, -1.6301578422923519, -1.4332006785118301, -1.305707665299106, -1.364200787821644, -1.5357935460809622, -1.6161992336951241, -0.74003518668370516, -0.29423824173210689, 0.025934598230976654, -0.043349004411304674, -0.25408021803022468, -0.066965686484977499, -0.075717498698635255, 0.007057189465364498, -0.042171356658338113, -0.036938315661768008, -0.34221561581756049, -0.20400167508805764, -0.37417116097079772, -0.25039909487805356, -0.070874531394524931, -0.0569972852039487, -0.067238206950403182, -0.17397285212300442, -0.20428337307808273, -0.23651154356493315, -0.33356498933276568, -0.07339749754226077, -0.70367959806681601, -0.82403680021595049, -1.6058616381755235, -1.6192427030685497, -1.5705638815427956, -1.4659201063980019, -0.95504179549951018, -0.97237526162739873, -1.0460191987834688, -0.91465668941265721, -0.60548232361398524, 0.01898438364933451, -0.19419044456729498, -0.039627851124307223, 0.0012357796666701798, -0.078110822445325079, 0.0048626364920250518, -0.040449089662379589, -0.0035054269587873454, -0.13387544724730729, -0.10031131456276647, -0.25968674675684189, -0.20555329767005767, -0.26509289948725284, -0.038788452621647145, -0.076999891872251258, -0.071661433038976499, -0.14182240789719938, -0.1654673053291095, -0.19859450279267193, -0.053382326365810369, -0.2156585383674445, -0.045097357284793499, -0.62449818579949512, -0.92624906744917224, -1.0411254782363617, -1.122035196738675, -1.0607692164246043, -0.57723811773534028, -0.63187735896388075, -0.54813311204421922, -0.55320252101738743, -0.30197299587482401, -0.047213249757838388, 0.082808930467383288, -0.067715134483222431, -0.01022881748368659, 0.042038311258956552, -0.063371767399980669, 0.029161890169972702, -0.091396316586836127, -0.0034600735070754811, -0.12424052925006424, -0.24432996418012101, -0.26521664175359499, -0.22745980283820413, -0.14361316535317664, -0.00075904203100577935, -0.020936168457862139, -0.14205665196423617, -0.19024248288823023, -0.079686122362245204, -0.15016133237735926, 0.049598910651295514, -0.11760486834511712, -0.1837522251545049, -0.38594205494114608, -0.53542516436999843, -0.57340991730807989, -0.52753621424018138, -0.23151163972118355, -0.22295096919949259, -0.33704349161770436, -0.26165852514054583, -0.13898866968588663, 0.034596483191139484, -0.012631210076789067, 0.047371310076345617, -0.038651839330751551, -0.0019970761454430018, 0.063048845258375494, -0.1124891762554399, 0.08556992539656616, -0.21043659051868199, -0.19223333969456, -0.39082994830035861, -0.19294368007162721, -0.41025595439938572, -0.17178084419175166, -0.010933041190555012, -0.089512936152074493, -0.21569610281495066, -0.09144756671688016, -0.19525258909505316, -0.029753598134641936, -0.021307245660079924, 0.029087127940551009, 0.037511290653097842, -0.20600990120705839, -0.26967580750352926, -0.21000923681194664, -0.28209018858285628, -0.11925518789339556, -0.24869348141289982, -0.21025892926356746, -0.15567029136726124, -0.040546729108395907, -0.0050266153100547101, 0.030710887069787196, -0.0061104340245858278, 0.0369376092260571, -0.054862661367900321, 0.013297880203253048, 0.19659447375886394, -0.2499491329142558, -0.062959699002865757, -0.53055029095956008, -0.38784811281629444, -0.53891285075962392, -0.41886712861154285, -0.099230097260325875, -0.16474199810952628, -0.28693665642627014, -0.0095667980850221105, -0.32619954993450928, -0.08627491478166284, -0.073253161755714766, 0.015634174038690329, 0.082440536547531792, 0.025411878261881942, -0.11318909242737961, -0.1270560226842935, -0.21657212936164139, -0.13993873549611191, -0.37510275237622831, -0.26472923111076219, -0.24460131567533192, -0.14127652303494026, -0.050428686591045178, 0.041347840374190772, -0.0061780445153000636, 0.0073990345210250153, -0.014062739037014381, 0.14348925152561878, -0.015321787554403667, 0.0017746672356015968, 0.25165135427361052, -0.626463828190993, -0.48167134330805639, -1.045863293770664, -0.69512591788493194, -0.44532127384388254, -0.28479724025368391, -0.39470955087317983, 0.20227228344720469, -0.53909912073488953, -0.12025629051789474, -0.1899243750597305, -0.048474806721595133, 0.060764771353227762, 0.090648151782516159, 0.091608208912697275, 0.0036582478916540977, -0.22492530005263131, -0.27295314658024766, -0.35559738025257359, -0.62902925014412947, -0.57166411974881004, -0.37258895173129181, -0.22157638610464933, 0.022494427132080854, 0.014769425415166171, 0.003526808789406817, -0.011346909674078769, 0.050921170848348289, 0.090308541799219627, 0.37260817254533324, -0.25909871392159911, -0.42379280974334355, -0.095380647808568128, -1.1906083748893519, -0.78599914414892469, -0.95277914352730275, -0.63659778359422337, -0.98026015008952749, 0.48173198285916102, -0.60092009018055192, -0.10265418316164113, -0.39913639006279306, -0.17310908908773887, -0.0194191171632387, 0.054047965289179878, 0.1388529643463832, 0.15661099050145999, -0.10898263774416243, -0.33291231456737602, -0.59569027865888713, -0.69353081584948972, -1.0999707493347484, -0.74392084753736687, -0.49074781214158159, -0.065190556733852961, 0.012289768389229717, 0.024577513704595676, 0.0040302804696096322, 0.036047756292976456, 0.058236765637246286, 0.13893846256790621, 0.036944676036934632, 0.41686279554239464, -0.85232286388185818, -1.2988315127624981, -0.47352779677305168, -0.81763632541546793, -0.77384457803621831, -1.4256240004519281, 0.52588993532360684, -0.89821724022902683, 0.1591911967653899, -0.55046596772346867, -0.30980016041271019, -0.16709614007114884, -0.046029700131955266, 0.044793268150423983, 0.1689242242845459, 0.14412365934528507, -0.0088250071313367359, -0.36778545124666312, -0.79393844517732104, -1.1610479066529615, -0.76523210008850662, -0.63009858032048405, -0.13947023057344932, -0.017173105577524262, 0.039030007688455846, 0.014491273083805401, 0.039792542943837252, 0.054072846696920814, 0.11729310469925348, 0.053609281522667675, 0.0081549498718087084, -0.30910813452845548, 0.25944224899607843, -1.3584842180322938, -1.5885570490138659, -0.65759582794618221, -1.139869490652734, 0.70928264080594694, -1.9674198903133462, 0.37712664425406606, -0.84336038390578949, -0.47788074719428036, -0.18342000086663721, -0.18811394573901796, -0.055050027645985648, 0.045043056834335606, 0.11486303559854361, 0.22023958716404868, 0.14735402009444676, -0.27894427087197998, -0.73080536953129638, -0.76794305693297227, -0.37355919765840223, 0.12353986794322802, 0.090505348376311842, 0.14069908672094206, 0.087373214380278855, 0.023353946735568523, 0.031400559920396587, 0.079550230446202241, 0.084927161382185437, 0.040777158255349423, -0.16274954314482293, -0.41184413435479567, -0.71871288822574875, 0.55302907456342854, -1.5309493464500674, -2.9026104205694736, 0.42043303599508353, -1.7138106264793671, 0.29513888249127102, -1.2517216433630918, -0.66769942176516839, -0.28576739334390183, -0.24127777006787937, -0.10778095858902549, -0.036092425009198861, 0.021519213385077923, 0.13414694961717147, 0.16917378957839613, 0.17307922682581758, 0.076246758829015673, -0.047904835134272621, -0.27544262702406924, 0.61826249566563185, 0.26987423123693399, 0.2085883517320696, 0.26073426210721973, 0.12070625812911842, 0.062945582093309679, 0.083649573916505432, 0.049688095345785867, 0.019564357607843069, -0.046035817476596949, -0.13409074070830324, -0.49027201814294552, -0.47756457321420159, -0.74403675135427549, -0.3080068432033089, -0.043712438842705037, -4.735594317158907, -0.043712438842706695, -0.30800684320330962, -0.74403675135427572, -0.47756457321420304, -0.49027201814294813, -0.13409074070830412, -0.046035817476598156, 0.019564357607843069, 0.049688095345786006, 0.083649573916506056, 0.062945582093310845, 0.12070625812911921, 0.26073426210722073, 0.20858835173207019, 0.26987423123693399, -0.37355919765836759, -0.27544262702403433, -0.047904835134273127, 0.076246758829012523, 0.17307922682581853, 0.16917378957839499, 0.13414694961716844, 0.02151921338507657, -0.036092425009199861, -0.1077809585890261, -0.24127777006787943, -0.2857673933439015, -0.66769942176516905, -1.2517216433630949, 0.29513888249127429, -1.7138106264793713, 0.42043303599507681, -2.902610420569474, -1.5309493464500692, 0.55302907456342232, -0.71871288822575019, -0.41184413435479833, -0.16274954314482265, 0.04077715825534866, 0.08492716138218645, 0.079550230446203143, 0.031400559920398419, 0.023353946735571576, 0.08737321438028138, 0.14069908672095732, 0.090505348376334033, 0.1235398679432393, -0.76523210008847808, -0.76794305693296139, -0.73080536953128505, -0.27894427087197604, 0.1473540200944477, 0.22023958716404682, 0.11486303559854165, 0.045043056834333829, -0.055050027645986453, -0.18811394573901843, -0.18342000086663854, -0.47788074719428042, -0.84336038390579149, 0.37712664425406617, -1.9674198903133469, 0.70928264080593695, -1.1398694906527307, -0.65759582794619398, -1.588557049013867, -1.3584842180322987, 0.25944224899607732, -0.30910813452845781, 0.0081549498718086911, 0.053609281522667279, 0.11729310469925426, 0.054072846696921202, 0.039792542943838709, 0.014491273083807311, 0.039030007688458185, -0.017173105577517028, -0.13947023057343994, -0.63009858032045107, -1.0999707493347308, -1.1610479066529467, -0.79393844517731305, -0.3677854512466584, -0.0088250071313340107, 0.14412365934528559, 0.16892422428454401, 0.044793268150420118, -0.046029700131956147, -0.16709614007115095, -0.30980016041271097, -0.55046596772347045, 0.15919119676539073, -0.8982172402290286, 0.52588993532360329, -1.4256240004519327, -0.77384457803621687, -0.8176363254154656, -0.47352779677305679, -1.2988315127625027, -0.85232286388185829, 0.41686279554239525, 0.036944676036935756, 0.13893846256790574, 0.058236765637246675, 0.036047756292977066, 0.0040302804696111128, 0.02457751370459911, 0.012289768389232913, -0.065190556733844662, -0.49074781214156804, -0.74392084753735632, -0.62902925014412903, -0.69353081584948562, -0.59569027865888302, -0.33291231456737491, -0.10898263774416028, 0.15661099050145985, 0.13885296434638142, 0.054047965289177706, -0.019419117163239467, -0.17310908908773912, -0.39913639006279433, -0.10265418316163986, -0.60092009018055315, 0.48173198285915786, -0.98026015008952594, -0.63659778359422126, -0.9527791435273002, -0.78599914414892458, -1.190608374889349, -0.095380647808570002, -0.42379280974334488, -0.25909871392159683, 0.37260817254533357, 0.09030854179921953, 0.050921170848348372, -0.011346909674079158, 0.0035268087894081549, 0.014769425415168456, 0.022494427132082863, -0.22157638610464575, -0.37258895173129003, -0.5716641197488066, -0.37510275237622537, -0.35559738025257059, -0.27295314658024672, -0.22492530005262792, 0.0036582478916564426, 0.091608208912696387, 0.090648151782514966, 0.060764771353224882, -0.048474806721595647, -0.18992437505973167, -0.12025629051789351, -0.53909912073488875, 0.20227228344720258, -0.39470955087317799, -0.28479724025368247, -0.44532127384387832, -0.69512591788493272, -1.04586329377066, -0.48167134330805861, -0.62646382819099156, 0.25165135427361029, 0.0017746672356018336, -0.0153217875544032, 0.14348925152561842, -0.01406273903701487, 0.0073990345210243587, -0.0061780445152985596, 0.04134784037419488, -0.050428686591041855, -0.1412765230349349, -0.2446013156753272, -0.26472923111076024, -0.11925518789339257, -0.13993873549610955, -0.21657212936163839, -0.1270560226842922, -0.11318909242737903, 0.025411878261882927, 0.082440536547530169, 0.015634174038688685, -0.073253161755715501, -0.086274914781661965, -0.326199549934509, -0.0095667980850238903, -0.28693665642627003, -0.16474199810952764, -0.099230097260324029, -0.41886712861154318, -0.53891285075962314, -0.38784811281629461, -0.53055029095956219, -0.062959699002866631, -0.24994913291425488, 0.1965944737588636, 0.013297880203252755, -0.054862661367901897, 0.036937609226056677, -0.0061104340245862225, 0.030710887069788338, -0.005026615310052167, -0.040546729108393256, -0.15567029136725916, -0.21025892926356554, -0.24869348141289621, -0.23151163972117689, -0.28209018858284918, -0.21000923681193823, -0.26967580750352416, -0.20600990120705304, 0.037511290653099091, 0.029087127940549885, -0.02130724566008323, -0.029753598134642099, -0.19525258909505444, -0.091447566716882075, -0.21569610281495041, -0.089512936152075853, -0.010933041190555782, -0.17178084419175305, -0.41025595439938806, -0.19294368007162768, -0.39082994830036216, -0.19223333969456258, -0.21043659051868269, 0.085569925396567076, -0.11248917625543933, 0.063048845258374231, -0.0019970761454456269, -0.038651839330752197, 0.047371310076345617, -0.012631210076786959, 0.034596483191142599, -0.13898866968588444, -0.26165852514053983, -0.33704349161769737, -0.22295096919948695, -0.57723811773534028, -0.52753621424018138, -0.57340991730807944, -0.53542516436999865, -0.38594205494114614, -0.1837522251545064, -0.11760486834511884, 0.049598910651293758, -0.15016133237735926, -0.07968612236224891, -0.1902424828882312, -0.14205665196423831, -0.020936168457862579, -0.00075904203100844866, -0.14361316535317845, -0.2274598028382093, -0.26521664175359499, -0.24432996418012529, -0.12424052925006639, -0.0034600735070760831, -0.09139631658683596, 0.029161890169972428, -0.063371767399980516, 0.042038311258955005, -0.01022881748368659, -0.067715134483221959, 0.082808930467383746, -0.047213249757837236, -0.3019729958748239, -0.55320252101738743, -0.548133112044219, -0.63187735896388053, -0.95504179549950285, -1.060769216424599, -1.1220351967386673, -1.0411254782363524, -0.92624906744916458, -0.62449818579949246, -0.045097357284792555, -0.21565853836744897, -0.053382326365811708, -0.19859450279267432, -0.16546730532911214, -0.14182240789720132, -0.07166143303897729, -0.076999891872253062, -0.038788452621649434, -0.2650928994872585, -0.20555329767005678, -0.25968674675684078, -0.10031131456276626, -0.13387544724730568, -0.0035054269587865765, -0.040449089662379971, 0.0048626364920241281, -0.078110822445325467, 0.0012357796666695618, -0.039627851124306598, -0.19419044456729473, 0.018984383649339364, -0.60548232361397991, -0.91465668941264988, -1.0460191987834631, -0.97237526162739263, -1.3057076652991049, -1.4659201063979992, -1.5705638815427927, -1.6192427030685486, -1.6058616381755215, -0.82403680021595249, -0.70367959806681868, -0.073397497542269388, -0.33356498933276529, -0.23651154356493967, -0.2042833730780847, -0.17397285212300875, -0.067238206950403417, -0.056997285203952995, -0.070874531394526111, -0.25039909487805306, -0.37417116097079761, -0.20400167508805389, -0.34221561581755761, -0.036938315661763657, -0.042171356658337315, 0.0070571894653653896, -0.075717498698634964, -0.066965686484977194, -0.25408021803022474, -0.043349004411301621, 0.025934598230977574, -0.29423824173210122, -0.74003518668370272, -1.6161992336951192, -1.5357935460809593, -1.3642007878216427, -1.5169157844507262, -1.4332006785118279, -1.6301578422923491, -1.932677889629856, -1.8694842425148879, -1.8503084663080056, -0.59040483751417916, -0.65601747033471336, -0.40132058236311047, -0.31836783196552787, -0.32035624605031593, -0.22965854401457814, -0.082000200083028219, -0.089025163302698024, -0.064732391018913552, -0.30366967948051288, -0.4730230160625184, -0.37970560392275871, -0.19430407325814622, -0.1797093251603995, -0.063069709957271444, -0.10088329756370083, -0.097060879184794432, -0.18018565146479387, -0.19830827357221226, -0.43479836723774673, 0.22038139379151372, -0.26457398838911428, -1.79452817569422, -2.2062031109308391, -2.102356658927595, -1.6930052820597257, -1.2756392437278008, -1.6712324532934884, -1.8237453905245253, -1.4884613952536252, -2.384654945418673, -1.5836466114968115, -1.4772075422260404, -0.74950958881152596, -0.561795414852579, -0.46670603976586306, -0.35555823460889052, -0.22514761090777982, -0.18963405157839525, -0.0068547033925124142, -0.051036689850529192, -0.53565359443295624, -0.48685989567066656, -0.57809021675862349, -0.39240983780018618, -0.11286619230167973, -0.12321821342876334, -0.089790022871127112, -0.13792125740417074, -0.16098142942852883, -0.30105748355308298, -0.18584599820379807, -0.75998936793718352, -0.8428129384777584, -2.1809972632688073, -2.8791053946930352, -2.5682454610054237, -2.0983496132652038, -2.219376425173226, -2.3602336858020396, -1.4129529949598048, -0.95341593939478875, -0.38769722219534936, -2.2603422126144772, -1.394644866167148, -0.77664281720103345, -0.64028425429629032, -0.41856090178724664, -0.3216120132349809, -0.33607304327303461, 0.065628582466230781, -0.098136330355478765, -0.44979914971695495, -0.82222257402878096, -0.39156143807433802, -0.36045886741397631, -0.29008721857562392, -0.035472107498135542, 0.027920290231535812, -0.117989151853845, -0.1892849156146684, -0.25489431004647656, -0.3625272253546275, -0.30143411016838906, -0.75242583405872021, -2.0271757806780628, -2.2639673745086539, -3.5383249085863659, -3.361597106712324, -2.3438802466919229, -1.7249706159518579, -1.7990472945779559, -3.3061472495599995, 0.21156442062862166, -0.74302378668934399, -1.8514918971041745, -0.82502198370632651, -0.74351999988260331, -0.39569215798819279, -0.3943138014933833, -0.29890408748686254, -0.28084584584372846, 0.01472791435408881, -0.04442075186168376, -0.53500016384583715, -0.90188577233528688, -0.69844971117515353, 0.14583567590324595, -0.5544445310046473, 0.1976534820339324, -0.18710828510427244, 0.0076800871037496377, -0.32808050503226982, -0.17511899827005836, -0.38519889843539723, -0.26116553014987143, -0.61792590829173255, -1.4260964500310052, -2.8480835086962002, -2.9656385065341997, -3.1259631743419583, -2.2776462897026373, -1.3820969331049018, -2.0505291279591713, -2.4084298023013879, -1.4996172161865962, 0.12372650763830863, -1.1981932989988076, -0.64588897407998824, -0.33546767337818667, -0.36301254306387043, -0.26196208085033179, -0.17675902360982099, 0.0055639281302357606, -0.034749316729180774, -0.44871907432473696, -0.67215259521100923, -0.46007833291523831, -0.67254330182605182, 0.11741827111366224, -0.409871486550618, -0.32892265898101625, -0.13145556961021479, -0.028899376452251727, -0.092298612924564649, -0.16672095020503341, -0.17646388464205828, -0.35131617290772521, -0.19150029104902661, -1.2217002325587201, -0.7866727441045076, -1.9807415717551959, -1.978633957921319, 0.069560471269209931], "bottom": {"real": [4103.3252596935745, 31959.928439656338, 10854.934870050551, 5174.7646941682715, 2670.3793024702013, 1512.8812431609856, 751.72119813508266, 487.34157279751093, 286.27976884850017, 202.21445228809756, 139.363320073941, 96.326676625874271, 67.416513392704019, 55.036039361563731, 42.617455049491909, 37.327841235406673, 35.198800209060273, 37.327841235406588, 42.617455049491802, 55.036039361563766, 67.416513392704019, 96.326676625874285, 139.36332007394108, 202.21445228809804, 286.27976884850017, 487.34157279751093, 751.72119813508289, 1512.8812431609856, 2670.3793024702018, 5174.7646941682751, 10854.934870050551, 31959.928439656363, 12454.694619943468, 7821.5833902765553, 5473.1790170642225, 2925.2286142376206, 1403.2127508507554, 917.05530556073552, 556.73350878905819, 335.58154911349368, 222.7562369115075, 161.71079893305554, 119.4497628246793, 75.609007514321249, 55.496087080936569, 43.998829489125107, 34.725029965122339, 29.983374804996487, 29.187336608781969, 30.714909872552553, 33.135728528562289, 38.780040560556557, 50.11926248444739, 62.426609296740132, 93.916765363567279, 123.96413175241418, 177.16967383039952, 250.50030243800805, 399.94920918463373, 596.1485322845399, 914.24633406931139, 1871.6210271277439, 4518.4223121248042, 13565.815861293135, 16084.742683461694, 10028.519769850123, 2736.2851033168113, 1377.4551350842332, 614.08174831750455, 382.39730464420114, 237.0105878631189, 156.24359018004319, 129.95938769710136, 95.53783206710068, 72.004092864891931, 47.804301653843083, 38.41781199466849, 32.452048622414502, 26.753427300507923, 23.772936248165699, 23.138404805980134, 23.598476471031617, 24.755859033283485, 28.713323989162731, 33.395537201677122, 40.850586549891439, 58.649881806718739, 74.872968711973769, 93.465129226367807, 123.19419955144703, 174.75706127058839, 262.71291650117263, 321.82068054258934, 657.05253635266399, 2163.5932265202309, 10212.960963472207, 3792.0213246064613, 2759.3366542985627, 1627.1011647050395, 788.44977202016776, 362.8509317865861, 253.90720770691448, 163.04342130809295, 117.95146004773997, 90.766106703902594, 66.207745096840526, 48.204553381452804, 35.429206551568903, 28.049881805648454, 23.25027473117818, 20.778936642061399, 19.004228801577, 17.585642163629327, 17.698181326434501, 18.806836162280465, 20.329571180523736, 23.456998427374465, 27.472702254518477, 37.193120035742723, 49.117252584083957, 59.574829012615233, 73.59994664128709, 112.97176733843995, 181.91972084309376, 284.0343016488693, 486.29648203694052, 857.05287855361007, 2037.977143592303, 2057.7285052573056, 2152.5952706253152, 1395.1090523951752, 736.25297680000074, 343.98700964912916, 189.68478304615005, 127.37774106216496, 91.12789293157843, 73.667255133763959, 55.964360327653644, 39.482567042532949, 28.14219415335706, 21.278934963706885, 18.193385040510105, 16.473354788100497, 15.086583853495943, 14.403945056404867, 14.533202056236952, 15.306988390608382, 16.092687824041843, 18.097466979870337, 20.289280537832838, 26.201109009342694, 34.023571220637564, 41.620492531599325, 50.685682074964014, 82.560701981631325, 127.19888958323958, 237.89761616945128, 410.06312322518994, 1062.2303232610248, 1612.0404058137353, 2295.1409914972487, 1787.3905923922546, 1192.2295048012345, 614.24882525880628, 260.82616895243024, 126.23242010647614, 84.151076288810984, 67.709414992782712, 60.122571559472298, 42.830591238304876, 28.733103940874788, 20.365121706656215, 16.211112474155353, 14.072758594539286, 13.483685068827034, 12.859628868618824, 11.868475605254234, 11.903201306554562, 12.816996745648828, 13.065794209061782, 14.312459824747068, 16.296926608708432, 20.657711991677495, 27.653390070235432, 33.593110413967857, 40.672720076575544, 55.856624618502167, 103.16047843117397, 179.23484372919035, 453.7513605151255, 902.26285048256875, 1683.7179352249004, 1784.083505146898, 1577.4265763170067, 936.66309122894188, 441.20892337587179, 176.43359667751182, 98.093971741535682, 65.995944695036641, 53.240295707495449, 45.85166507919449, 32.407485359783081, 22.496202298890402, 16.965613714417799, 13.731573445856062, 12.28236966845588, 11.747735381447885, 11.213125876643861, 10.807806034266576, 10.840341477375139, 11.066245600125107, 11.119452781179984, 12.424913044930788, 14.610027556462221, 16.887227742677396, 22.29462391228396, 27.978157381323118, 31.991250392971789, 44.052658881876532, 76.109568327798371, 159.13944268405785, 318.39207128278571, 686.00323178071869, 1336.6568589814267, 1043.0649603599104, 984.78746182807288, 618.0395600950327, 289.32426118556657, 132.98214831862998, 71.520048430881175, 52.57629039600819, 41.525598741467476, 33.633912722813989, 27.031255662449681, 19.489513580793098, 14.812681614273632, 12.21137274400836, 10.893625186679536, 10.482989068673637, 10.105487112246305, 9.7116899243817354, 9.9488069804828818, 10.161347795217756, 10.513274350469635, 11.294150924355744, 12.914730156139361, 14.478861048855546, 19.021661277112585, 23.907610167423496, 27.23573455134931, 38.742976413983023, 62.869223125902629, 111.92503010834605, 213.21062569137553, 433.61330953226366, 726.70269845820769, 562.21861410525219, 539.290746631297, 362.9446461846826, 211.46669660189423, 107.79772661917396, 60.676668375567573, 41.743276533116536, 34.42369696468284, 31.26708433258414, 25.313030406949355, 18.452791878453507, 13.956644256748325, 11.043046695375654, 10.111101652499672, 9.4550360444473061, 9.2429593469396529, 8.9521911222399257, 9.2321124164010211, 9.1908041669169815, 9.5662827353227868, 10.506864865879585, 11.825879962774797, 13.52813582962821, 15.999059082232355, 19.794027285196304, 22.478845287715099, 32.163272384867753, 47.951984523863096, 81.309242866655126, 150.93744536633105, 263.84630525991662, 439.27951033199258, 309.78719189559973, 306.63257299287005, 244.19689848939953, 151.98869039704036, 91.330482004276163, 57.352753322013349, 38.709160706067429, 33.087086449001383, 29.826655436967027, 23.987725148031473, 16.89184234468231, 13.041947079601194, 10.535809562752126, 9.5359471621683909, 8.9374483496855426, 8.5672643122912326, 8.4572735895659434, 8.4602804400971099, 8.5543457526330293, 9.3823763336699937, 10.377230417708629, 11.799673812944503, 12.361176855966248, 14.786798351390814, 16.739331260686697, 20.590148031359199, 25.753571174908508, 38.382893067866803, 63.23013835373618, 105.62437218489313, 175.68359833526657, 256.05708017959813, 224.80770864957879, 228.96753655549054, 172.54966320095522, 116.26546821946491, 80.46990226978717, 53.844334876610212, 37.483534347728245, 31.075867215997253, 26.284110672635684, 20.645029038002825, 15.025930043703783, 11.738276392873866, 10.201438772469425, 8.8203445227279982, 8.4945509856995365, 8.1968728409344909, 8.1790777304419588, 7.908377156922052, 8.6785477295074038, 9.6473715156890378, 10.339435934253908, 10.728035799158873, 11.983620383388951, 13.687783504221503, 16.433795900162693, 18.309125668572698, 22.485075799802843, 31.764684376383052, 49.270779583367755, 77.729851956279916, 123.50059845139852, 179.77043479461938, 159.57955623939222, 154.53329869380329, 127.96891349444883, 91.96783661678981, 66.921744646417011, 48.818287955893446, 36.591061653526779, 27.834277069623926, 23.041902152181589, 18.358713894245302, 13.567338391039883, 11.281789781943191, 9.6038295455647198, 8.7333277563772516, 8.338123667351045, 7.8992088061869676, 7.7104322426775909, 7.9983869985641034, 8.5436438490902269, 9.4432158844800043, 9.7870344410918424, 10.423788563863184, 11.253535869282553, 13.327727562718719, 15.178485067808285, 17.337135682195893, 21.668177060872456, 29.324826866357235, 41.492753218636352, 60.042185621507166, 90.318134916215342, 130.8805255687621, 104.24612656823803, 108.8028525877362, 89.555745560157249, 69.785669575239666, 56.746590174428, 44.171218069814536, 32.39838941697333, 25.182252759236459, 21.610810960419155, 16.972539478480535, 13.238333358456819, 10.614964653675054, 9.3806927946307859, 8.7492814832421075, 8.135298996501481, 7.634398163782139, 7.4852869889479292, 7.8316993609624435, 8.3853957415274643, 9.2315742305614634, 9.8389441629514209, 10.386606043801919, 11.354519695989614, 13.10793814373473, 15.240615596988986, 17.58979203150065, 21.461767190818147, 26.661531554562984, 36.123337028978582, 50.320731869274383, 69.065609719997795, 87.442260857354199, 76.677643276575949, 72.689147636714068, 65.50785565022079, 58.297658406357961, 49.558598741148941, 39.192968104211104, 29.954333787314212, 23.9458863540046, 19.99476230299754, 16.735792774046942, 12.901585018586687, 10.289994503175569, 9.3977765822679924, 8.6704196287803228, 7.8828328381527575, 7.5341478839019471, 7.4312929200041102, 7.7985159766257679, 8.1869060912414078, 9.142118635191288, 9.830770779287306, 10.340589147553608, 11.281171263863113, 12.984935438318947, 15.107976469419242, 18.19770629294057, 20.668346478928893, 25.761496132514932, 32.251424266134499, 39.95295492008993, 51.369438078625848, 65.569299562662465, 55.468849468729523, 55.86542097352833, 55.101138782973663, 50.531848368674723, 43.701470536966781, 35.866381440857431, 28.439396618258566, 22.42863349391925, 18.945030358082761, 15.596977301337661, 12.206893729284205, 10.050508612628137, 9.5473939482322869, 8.6463276280830179, 7.910774588392556, 7.2755514089661562, 7.5482842032534565, 7.5443879419641391, 7.9636187803325598, 8.9922015447577355, 9.6351997849989068, 10.0402383165895, 11.706284711344862, 13.336194440398188, 15.639893484781382, 17.791763978880017, 20.83131115781941, 25.6818621728188, 31.387128642704546, 37.143166049555219, 44.038766386780296, 52.904506937405849, 48.059389050535145, 50.072432061388255, 49.115010315515249, 44.160423468831148, 38.440953181308423, 32.35840674752788, 26.251111011761232, 21.081922570464979, 17.515101530242855, 14.10261621013594, 11.438986216849498, 10.422223192105227, 9.76207732514108, 8.7746319169344158, 7.5939451837729885, 7.2605949806802883, 7.4766919496025244, 7.419064545103371, 8.0187357229163059, 8.9266344512172839, 9.6975643711848107, 10.35689464395745, 11.451859179394017, 13.10319606057651, 15.460942026724263, 18.825021564083144, 22.773162722002358, 26.906890973713775, 31.803276985208164, 37.82902190094245, 42.388788349798304, 44.620105679799558, 48.859277780953818, 46.136110292205181, 41.908361785717766, 39.521022744549988, 35.052294933716347, 29.181343166303421, 23.802472557875006, 19.955746539759069, 16.22743816874862, 13.086681034223906, 11.226909204888067, 10.616041272149978, 9.8385492842648201, 8.6563005846195669, 7.8952678290472065, 7.2762149925656852, 7.3094046208482961, 7.276214992565686, 7.8952678290472083, 8.656300584619574, 9.8385492842648219, 10.616041272149991, 11.226909204888063, 13.086681034223904, 16.22743816874862, 19.955746539759044, 23.802472557874989, 29.181343166303414, 35.052294933716361, 39.521022744549988, 41.908361785717766, 46.136110292205167, 48.059389050535053, 44.620105679799494, 42.388788349798226, 37.82902190094255, 31.803276985208186, 26.906890973713825, 22.773162722002368, 18.825021564083173, 15.460942026724263, 13.103196060576508, 11.451859179394024, 10.356894643957441, 9.6975643711848196, 8.9266344512172751, 8.0187357229163005, 7.4190645451033674, 7.4766919496025226, 7.2605949806802919, 7.5939451837729957, 8.7746319169344176, 9.7620773251410817, 10.42222319210523, 11.438986216849505, 14.102616210135931, 17.515101530242841, 21.081922570464972, 26.251111011761218, 32.358406747527845, 38.44095318130838, 44.160423468831063, 49.115010315515164, 50.072432061388007, 55.468849468729232, 52.90450693740565, 44.038766386780239, 37.143166049555234, 31.387128642704567, 25.681862172818796, 20.831311157819414, 17.79176397888002, 15.639893484781378, 13.336194440398192, 11.706284711344864, 10.040238316589498, 9.6351997849989122, 8.9922015447577301, 7.9636187803325527, 7.5443879419641293, 7.5482842032534556, 7.2755514089661553, 7.9107745883925578, 8.6463276280830215, 9.5473939482322869, 10.050508612628141, 12.206893729284213, 15.596977301337661, 18.945030358082761, 22.428633493919239, 28.439396618258566, 35.866381440857388, 43.70147053696676, 50.531848368674645, 55.101138782973628, 55.865420973528131, 76.677643276575822, 65.569299562662337, 51.369438078625784, 39.952954920090001, 32.251424266134514, 25.761496132514949, 20.668346478928903, 18.19770629294057, 15.107976469419244, 12.984935438318949, 11.28117126386311, 10.340589147553612, 9.8307707792872954, 9.1421186351912915, 8.1869060912414096, 7.7985159766257599, 7.4312929200041067, 7.5341478839019365, 7.882832838152761, 8.6704196287803175, 9.3977765822679942, 10.289994503175565, 12.901585018586699, 16.735792774046928, 19.99476230299755, 23.945886354004585, 29.954333787314209, 39.19296810421109, 49.55859874114887, 58.297658406357826, 65.507855650220662, 72.689147636713827, 104.24612656823788, 87.442260857354213, 69.065609719997752, 50.320731869274368, 36.123337028978554, 26.661531554563005, 21.461767190818144, 17.589792031500657, 15.24061559698899, 13.107938143734735, 11.354519695989618, 10.386606043801917, 9.8389441629514156, 9.2315742305614528, 8.3853957415274607, 7.8316993609624372, 7.4852869889479301, 7.6343981637821372, 8.1352989965014757, 8.7492814832421111, 9.3806927946307876, 10.614964653675061, 13.238333358456819, 16.972539478480542, 21.610810960419155, 25.182252759236462, 32.39838941697333, 44.171218069814465, 56.746590174427972, 69.785669575239609, 89.555745560157177, 108.80285258773613, 159.57955623939191, 130.8805255687619, 90.318134916215342, 60.042185621507201, 41.492753218636388, 29.32482686635726, 21.668177060872456, 17.3371356821959, 15.178485067808294, 13.327727562718726, 11.253535869282553, 10.423788563863175, 9.7870344410918424, 9.4432158844800025, 8.5436438490902233, 7.9983869985640963, 7.7104322426775891, 7.8992088061869596, 8.338123667351045, 8.7333277563772533, 9.6038295455647145, 11.281789781943189, 13.567338391039879, 18.358713894245295, 23.041902152181613, 27.834277069623926, 36.591061653526772, 48.818287955893382, 66.921744646416983, 91.967836616789555, 127.96891349444861, 154.53329869380323, 224.80770864957842, 179.77043479461912, 123.50059845139837, 77.729851956279958, 49.27077958336772, 31.764684376383045, 22.485075799802868, 18.309125668572733, 16.433795900162703, 13.6877835042215, 11.983620383388955, 10.72803579915886, 10.339435934253904, 9.6473715156890254, 8.6785477295073967, 7.9083771569220449, 8.1790777304419553, 8.1968728409344909, 8.4945509856995347, 8.8203445227279857, 10.201438772469423, 11.738276392873869, 15.02593004370379, 20.645029038002832, 26.284110672635681, 31.075867215997253, 37.483534347728245, 53.844334876610127, 80.46990226978717, 116.26546821946468, 172.54966320095502, 228.96753655548974, 309.78719189559854, 256.05708017959756, 175.68359833526623, 105.62437218489312, 63.230138353736059, 38.382893067866796, 25.753571174908501, 20.590148031359234, 16.739331260686704, 14.786798351390802, 12.361176855966249, 11.799673812944498, 10.377230417708628, 9.3823763336699955, 8.5543457526330204, 8.4602804400970992, 8.4572735895659417, 8.567264312291222, 8.9374483496855408, 9.5359471621683873, 10.535809562752121, 13.041947079601202, 16.891842344682331, 23.987725148031487, 29.826655436966995, 33.087086449001376, 38.709160706067401, 57.35275332201325, 91.33048200427605, 151.98869039703999, 244.19689848939879, 306.63257299286875, 562.21861410525219, 439.27951033199281, 263.84630525991662, 150.93744536633125, 81.309242866655154, 47.951984523863167, 32.163272384867753, 22.478845287715121, 19.794027285196304, 15.999059082232373, 13.528135829628219, 11.82587996277479, 10.50686486587958, 9.5662827353227868, 9.190804166916978, 9.2321124164010211, 8.9521911222399257, 9.2429593469396458, 9.4550360444473043, 10.111101652499666, 11.043046695375656, 13.956644256748326, 18.452791878453535, 25.313030406949409, 31.26708433258414, 34.423696964682833, 41.743276533116529, 60.676668375567509, 107.79772661917396, 211.46669660189391, 362.94464618468271, 539.29074663129688, 1043.064960359907, 726.70269845820621, 433.61330953226252, 213.21062569137561, 111.92503010834588, 62.869223125902685, 38.742976413983001, 27.235734551349324, 23.907610167423499, 19.021661277112635, 14.478861048855546, 12.91473015613936, 11.29415092435573, 10.513274350469626, 10.161347795217747, 9.9488069804828747, 9.7116899243817407, 10.105487112246305, 10.48298906867363, 10.893625186679529, 12.211372744008356, 14.812681614273641, 19.489513580793123, 27.031255662449723, 33.633912722813953, 41.52559874146754, 52.576290396008147, 71.520048430881019, 132.98214831862995, 289.32426118556515, 618.03956009503077, 984.7874618280689, 1784.0835051468939, 1336.6568589814267, 686.00323178071801, 318.39207128278593, 159.1394426840578, 76.109568327798428, 44.05265888187656, 31.991250392971853, 27.978157381323108, 22.294623912284028, 16.8872277426774, 14.610027556462216, 12.424913044930783, 11.119452781179971, 11.066245600125109, 10.840341477375128, 10.807806034266576, 11.213125876643852, 11.747735381447885, 12.282369668455885, 13.731573445856061, 16.965613714417831, 22.496202298890438, 32.407485359783124, 45.85166507919449, 53.24029570749542, 65.99594469503667, 98.093971741535469, 176.43359667751176, 441.208923375871, 936.66309122894131, 1577.4265763170022, 2295.1409914972464, 1683.7179352249, 902.26285048256875, 453.75136051512629, 179.23484372919026, 103.16047843117398, 55.856624618502074, 40.67272007657553, 33.59311041396785, 27.653390070235471, 20.65771199167747, 16.296926608708404, 14.312459824747059, 13.065794209061782, 12.816996745648842, 11.903201306554557, 11.868475605254236, 12.859628868618808, 13.483685068827036, 14.072758594539282, 16.21111247415536, 20.365121706656261, 28.733103940874788, 42.830591238304876, 60.122571559472277, 67.709414992782641, 84.151076288810941, 126.23242010647601, 260.82616895243001, 614.24882525880446, 1192.2295048012338, 1787.3905923922525, 2057.7285052573061, 1612.0404058137367, 1062.2303232610243, 410.06312322519068, 237.89761616945111, 127.19888958323968, 82.560701981631325, 50.685682074964042, 41.620492531599325, 34.023571220637599, 26.201109009342684, 20.289280537832823, 18.09746697987034, 16.092687824041846, 15.306988390608383, 14.533202056236956, 14.40394505640487, 15.086583853495929, 16.47335478810048, 18.193385040510101, 21.278934963706888, 28.142194153357082, 39.482567042533006, 55.964360327653644, 73.667255133763959, 91.127892931578486, 127.37774106216477, 189.68478304614982, 343.98700964912905, 736.25297679999983, 1395.1090523951759, 2152.5952706253147, 3792.0213246064582, 2037.977143592301, 857.05287855360893, 486.29648203694126, 284.03430164886925, 181.91972084309373, 112.97176733844, 73.599946641287119, 59.574829012615204, 49.117252584084042, 37.193120035742723, 27.472702254518456, 23.456998427374469, 20.329571180523697, 18.806836162280476, 17.698181326434526, 17.585642163629327, 19.004228801577025, 20.77893664206141, 23.250274731178163, 28.04988180564844, 35.429206551568925, 48.204553381452833, 66.207745096840512, 90.766106703902594, 117.95146004774013, 163.04342130809306, 253.90720770691391, 362.85093178658599, 788.44977202016594, 1627.1011647050398, 2759.3366542985605, 16084.742683461691, 10212.960963472211, 2163.5932265202296, 657.05253635266433, 321.82068054258923, 262.71291650117286, 174.75706127058845, 123.19419955144714, 93.465129226367836, 74.87296871197394, 58.649881806718732, 40.850586549891496, 33.395537201677108, 28.713323989162756, 24.755859033283489, 23.59847647103166, 23.138404805980127, 23.77293624816566, 26.753427300507926, 32.452048622414509, 38.417811994668497, 47.804301653843162, 72.004092864892002, 95.537832067100879, 129.95938769710136, 156.24359018004333, 237.01058786311899, 382.39730464420052, 614.08174831750478, 1377.4551350842296, 2736.2851033168104, 10028.51976985012, 12454.694619943462, 13565.815861293133, 4518.4223121248006, 1871.6210271277425, 914.24633406931184, 596.14853228454001, 399.94920918463339, 250.50030243800833, 177.16967383039946, 123.96413175241405, 93.91676536356745, 62.426609296740118, 50.119262484447404, 38.780040560556571, 33.135728528562332, 30.714909872552628, 29.187336608781973, 29.983374804996448, 34.725029965122346, 43.998829489125086, 55.496087080936618, 75.609007514321277, 119.44976282467937, 161.7107989330556, 222.75623691150756, 335.58154911349339, 556.73350878905831, 917.05530556073529, 1403.2127508507556, 2925.2286142376206, 5473.1790170642225, 7821.5833902765453], "imag": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, "imag": [0.0, 1.2364045139795619, -1.417097645517349, -0.56399113050110405, 0.44612204752934542, -0.22857396866743723, -0.12233724123958876, -0.22175063370253431, -0.081847654478992143, -0.082144528420219798, -0.016319194688300422, -0.22782269816808789, -0.31553723757062957, -0.34627737582788948, 1.2289873494343531, 0.23714731979244019, 0.0, -0.23714731979246589, -1.2289873494343397, 0.34627737582788415, 0.31553723757063035, 0.2278226981680877, 0.01631919468830003, 0.082144528420217924, 0.081847654478992143, 0.22175063370252604, 0.12233724123959235, 0.22857396866743548, -0.4461220475293452, 0.56399113050110417, 1.4170976455173483, -1.2364045139795607, 1.3606412475160685, 1.4411916754234004, -0.33639737707180672, -0.40059289889731398, 0.077869450453742065, 0.11933000679564014, -0.37020072416606636, -0.12600631152419078, -0.29291723130603109, -0.050683391882318747, -0.2123566898458312, -0.16517843825664788, -0.37438920939802389, 0.24863698196016526, 0.3956337373323095, 0.63308466384569972, -0.42643957249522613, -0.28749903579454716, -0.42301568032975828, 0.27971623762695519, 0.37808744252168813, 0.15296582393741001, 0.10558339349178479, 0.088135922159038316, 0.16705843369252343, 0.44192642538865684, 0.29915225181421068, -0.047649837051313024, -0.91046045526812369, 2.0668504460046964, 1.0236217631956717, 0.095775027331123308, 0.74327388531916394, -0.034562131954570079, 0.86361123753821734, 0.14688086473050152, -0.48040960369964397, 0.72553155869356101, -0.87865322612555963, -0.4386512882244577, -0.48573187070285989, -0.32017932907631874, -0.2830749986604581, -0.26990648117356825, 0.02650526715020339, 0.070485969368273949, 0.54370981678608021, -0.067938463814195812, 0.03776680164210118, -0.5672727000472354, 0.70317235961628322, -0.021133356436660328, 0.38778419454697149, -0.021237638325898506, 0.22272564345822293, -0.028529610958430065, 0.38389116569218468, 0.54015192300945081, 0.60890484535134959, -0.58851980199844323, -0.44257995715790832, -0.31054780878465871, 0.29452233058397137, 1.4203730522981519, 0.55042973333868539, 1.0635980503754177, -0.33223282315731317, 0.34444399750770754, 1.4812987845914649, 0.28203486996779931, -0.50640814880924923, -0.63152717265524028, -0.629560147993244, -0.38143909917016799, -0.26024004494781861, -0.19111004712688462, -0.0074049384249937926, 0.048145781851611912, 0.1775395908501845, -0.010766872600667563, -0.30557543109041663, -0.12790446781789691, 0.26497067484017217, 0.49718432792560446, 0.10847923107858634, 0.18362798444374009, 0.075793328725680426, 0.060915837919074359, 0.19352569039908019, 0.72484131059467238, 0.44769863619371908, 0.65727781823869813, -2.7705947112358387, -1.8545676439466314, 1.2534804319442749, -0.33201271766592177, -0.88179245186152566, -0.49191197464707653, 0.25061306335904526, 0.36434278784189494, 1.1059306596340364, 0.5161851914994362, -0.56259342729123341, -0.79818525223670012, -0.50614275317411506, -0.33450192958610769, -0.24220392956380443, -0.29363217368931438, -0.008975674876931768, 0.015702042401798396, 0.0043287613096781611, 0.11769772351691517, 0.14847974776257877, -0.070404205425901215, 0.15992872157624444, 0.12352253662711207, 0.22950814490677918, 0.17507285276321555, 0.13747467644410705, -0.011397583876153142, 0.1103654239578318, 0.24404333943616102, 1.1434123865171273, -0.60566186502851971, -1.3365493166122577, -1.1338250057068979, -1.8155993788135394, -1.0591766157447149, -0.69010247253016843, -0.46612626705839738, -0.20353239335470796, 0.51663850122650623, 0.2659701357580968, -0.10657803670351902, -0.69730138646619577, -0.60106744735123385, -0.33542800506255299, -0.30564289179799509, -0.26672832014897818, -0.19552146385155383, -0.070818253229817343, 0.24716165183995586, 0.24871001073085031, 0.087906272749659523, -0.066392764898707615, 0.11463822183417864, -0.097387114548935505, 0.203187369929616, 0.14292038700683402, 0.20736891492564988, 0.013685307258532538, -0.024371827339033106, -0.07233585420108031, 0.49291782026044645, 0.048917424863979306, 0.2263553910474094, 0.10839909954530379, -1.1957108711898896, -1.074771579950879, -1.1175823884839045, -0.23901621577757862, -0.12877260884038383, 0.23623861798820722, 0.17111767858219226, 0.46208690546761177, -0.23679486634226876, -0.53178862418379869, -0.34021620386249335, -0.23658495591681358, -0.22275756609764855, -0.20572232603908905, -0.14494733127712761, 0.042569660785953992, 0.077733070765218709, 0.20558191361174552, 0.031719635617051253, -0.031869312918226284, 0.0092446671691655703, 0.056907964871031128, 0.0068114015565279478, 0.20608602740074444, 0.023033890597355222, 0.064944626021130644, -0.089711600875347838, 0.051877409795188095, -0.22207073062343796, 0.24248276916859771, 0.82932242402320455, -0.0013771249391571728, -0.059139017637367147, -0.22067925957890699, -0.040068334939789486, 0.30677005181582889, 0.34490255487274107, 0.22221882441368751, 0.30818729583490434, 0.041656273199714877, 0.04389672863559383, -0.30032839778423015, -0.24246246508902861, -0.17809969240648099, -0.14126159805998126, -0.15553553242798068, -0.11510723629505018, 0.024565620333015844, 0.099324973861333238, 0.11314958231531824, 0.10920785431914558, -0.0027478677647314949, -0.002197448124614016, -0.069940656171551499, 0.14656177726216579, 0.078770189361086429, 0.10241250352514093, -0.021646587213105684, 0.019211160710794505, -0.26281055993233693, -0.055146819140509458, 0.35375643597531126, 0.22411418130009836, 0.3875457037899096, 0.25423051056794166, 0.30581901500051717, 0.085911462662090085, 0.39107302489471046, 0.24699055160378858, 0.29387297238508514, 0.15373347718831995, -0.002984582508486469, -0.017260207971634491, 0.072827255227658205, -0.062089135258224848, -0.11929552297831249, -0.14454547803177953, -0.095268264507681985, -0.10089614648569632, 0.048761457544342565, 0.09609679137720413, 0.10529567915508571, 0.012037966893712556, 0.076343100588013288, -0.20546035926376272, 0.13459193297647368, 0.025732944443593878, 0.15693288111736364, 0.043761647973181578, 0.083479299629684353, -0.15251693996845875, -0.10657849804550705, 0.038641816351084518, -0.10231188547257405, 0.21957600785318707, 0.264774413728534, 0.35686545376814865, 0.29428485506006929, 0.49934135456989248, 0.18285689914378289, 0.23675610338446562, 0.10353827744151976, 0.063107949091445251, 0.03846587209017991, 0.13661717625997899, 0.098122755354854277, -0.0039906200934582137, -0.09009439550122221, -0.16647262366552176, -0.19237746196753996, 0.026647194061990561, 0.048509383863749068, 0.1654716680991411, -0.035914795239386038, 0.099273739290017232, -0.11264804489487676, 0.14792025966567318, -0.093791526140286519, 0.10021979561886898, 0.0080449729031339574, 0.25399382193861558, -0.10951630791791984, -0.033548122233328939, -0.024275040060575473, -0.15565555640319179, -0.035231426543991154, 0.063485173587351437, 0.22814734409189322, 0.10684932807789806, 0.2355554633425988, 0.10407564589066262, 0.015525256999673128, -0.011913601076877599, -0.0039295920043740662, 0.044701593559494178, 0.16547295203389104, 0.15907957720488838, 0.074664630781852093, -0.04243906497357175, -0.11246394245376416, -0.16972427151540079, -0.21530290856092721, -0.15195750501651711, 0.40553391670472716, -0.05507255324487554, 0.18488826926783164, -0.03502892748503212, 0.37001387119519319, -0.15325722566754377, 0.27168929694820759, -0.059701216280397838, 0.48380889084850937, -0.13017614472815323, 0.047673421329954523, 0.071305584258896423, -0.15733503034508964, -0.079593566039511662, -0.075485820780456944, 0.06620598245919887, -0.0012785498169557718, 0.063827392702187832, -0.0024400785561267094, 0.050942947866382475, -0.16777938781510365, -0.092876049496631446, -0.010378011540593757, 0.1594331602275639, 0.19523869399925237, 0.18464942262011227, 0.063285121982437173, -0.077221658551183223, -0.1491589460269398, -0.17912928670535297, -0.11490318286815507, -0.11189726101521234, -0.20487634390276235, 0.76498375921014561, -0.23800865420158882, 0.74083633586796371, -0.083846748018997191, 0.46372256496136105, -0.16834959597660981, 0.80860396068755125, -0.25555617026146216, 0.19463913189037535, 0.30097794980879877, -0.16302064990391507, -0.10690976611351247, -0.13497815095121923, -0.049772201284698683, -0.079451894223961053, 0.040407604161455342, -0.09194915079908203, -0.14731166465811507, -0.26499922119954511, -0.15560186997026057, -0.069496867327866862, 0.16530147427910571, 0.22324548858793775, 0.26718136730575287, 0.19166419512490171, 0.060936923017668537, -0.088626234734120041, -0.14658150061325853, -0.15988330961851069, -0.18630425886720847, 0.20226767939927157, -0.086142375874942501, -0.77766311308853053, 1.4090153420038987, -0.3038886842868152, 1.2888777036518135, -0.074470758177722834, 1.2482190395029773, -0.80123249840125998, 0.47055407026248136, 0.75577199845745624, -0.11037652666672627, -0.054458369263067428, -0.19918107544886771, -0.10899116252888987, -0.15332031366957818, 0.043654486348545761, -0.022154991408309264, -0.15239448898889341, -0.19372066727324463, -0.023769000914279543, 0.24457812931884867, 0.51722906014629877, 0.4305115414461968, 0.4535547011943914, 0.3493581281435395, 0.21076698450898321, 0.063956239140020049, -0.043266091906771414, -0.12307402164456596, -0.17588597679362872, -0.057160335632577657, -0.0784239564694168, 0.14377922602722487, -0.28610637890062873, -1.3867297053543113, 2.7076021550256568, -0.12637740826622818, 2.0046873376560574, -1.442367549533677, 0.41516901914168741, 1.1637625900709285, 0.02187980097289477, 0.2119473473283329, -0.11885489129254229, -0.1159981397382039, -0.18899269201836941, -0.058557612447368006, -0.043992442353640485, -0.026879139952720663, -0.019298206516241499, 0.11425895652801879, -0.10204731321609413, 1.0976378910336273, 0.78083298172023552, 0.7892686156410853, 0.65585934200294294, 0.41963945036689065, 0.21495377802894691, 0.10907743849349204, -0.015003791740152183, -0.098108352324783196, -0.099983461914768382, -0.14074054993981988, -0.049222223079134708, -0.60225632344115432, 0.14822551318702831, -0.20779420774376042, -2.0348285284263787, 4.0251861194251459, -2.7420491756647705, 0.80731183784095106, 1.5256751899261984, -0.015777910213839702, 0.56481061463658855, 0.094745090670160736, 0.032305441115857418, -0.15496296066268442, -0.070898922964745509, -0.15185446167165811, -0.099013378837322308, -0.033839198092832617, 0.20170996830425536, 0.31395340400286426, 0.91804938453195506, 0.90837480132247328, 1.1212514777242175, 0.8111537128593137, 0.65461252772518352, 0.44812018281017851, 0.29338133862249183, 0.080942506543726575, -0.0048876135806632161, -0.042045437831147933, -0.048279129686226115, -0.13748504879992385, -0.22867414001479172, -0.65206021382919555, -1.4414816884582773, 0.34900156189272369, 1.2663064061664397, -4.4519069438082273, 2.3981175155339263, 2.0366543361516527, -0.11574419894478603, 0.9938253915317542, 0.40253059498542132, 0.29787521235777348, -0.048128801274628118, 0.001834095689180732, -0.063578064856638958, -0.14911097873279969, -0.17638384543546787, -0.15061392115741762, -0.13638564133647435, 0.3292785127452974, 0.45646012748955422, 0.0, -0.098380202448484569, 0.2889284460821786, 0.39981384774007012, 0.45036129401449043, 0.22152969012775114, 0.11073273247786891, 0.024117745722163984, -0.0050280129542534689, -0.087106960572679512, -0.083653448195472838, -0.51521426958714545, -0.86124262209127589, -0.78939851130304506, -0.74423225935393889, -1.1208061672534484, 0.0, 1.1208061672534437, 0.74423225935393811, 0.78939851130304262, 0.86124262209127656, 0.51521426958714389, 0.083653448195472491, 0.087106960572679332, 0.0050280129542534689, -0.02411774572216354, -0.11073273247786865, -0.22152969012775126, -0.45036129401449054, -0.39981384774006945, -0.28892844608217838, 0.098380202448484236, -0.90837480132241455, -0.45646012748952708, -0.32927851274525655, 0.13638564133649642, 0.15061392115742586, 0.17638384543547556, 0.14911097873280194, 0.063578064856639402, -0.0018340956891807838, 0.048128801274627743, -0.29787521235777292, -0.40253059498541977, -0.99382539153175098, 0.11574419894478773, -2.0366543361516536, -2.3981175155339201, 4.4519069438082379, -1.2663064061664346, -0.34900156189272086, 1.4414816884582851, 0.65206021382919577, 0.22867414001479336, 0.13748504879992438, 0.048279129686228148, 0.042045437831148495, 0.0048876135806649092, -0.08094250654372509, -0.29338133862248678, -0.44812018281017507, -0.6546125277251672, -0.81115371285928006, -1.1212514777241573, -1.0976378910335893, -0.91804938453192131, -0.31395340400284805, -0.20170996830424559, 0.033839198092839695, 0.09901337883732686, 0.1518544616716597, 0.070898922964745925, 0.15496296066268506, -0.032305441115858605, -0.09474509067015939, -0.56481061463658966, 0.01577791021384034, -1.5256751899262007, -0.80731183784094906, 2.7420491756647789, -4.0251861194251397, 2.0348285284263898, 0.20779420774376556, -0.14822551318702654, 0.60225632344115643, 0.049222223079135728, 0.14074054993982205, 0.099983461914770352, 0.098108352324783543, 0.015003791740153814, -0.10907743849349157, -0.21495377802894516, -0.41963945036688871, -0.65585934200294016, -0.78926861564106898, -0.78083298172020743, -0.24457812931883488, 0.10204731321609872, -0.11425895652800651, 0.019298206516249916, 0.026879139952725083, 0.04399244235364598, 0.058557612447370934, 0.18899269201837068, 0.11599813973820411, 0.11885489129254237, -0.21194734732833254, -0.021879800972893531, -1.1637625900709312, -0.4151690191416833, 1.4423675495336772, -2.0046873376560503, 0.12637740826623128, -2.7076021550256497, 1.3867297053543175, 0.28610637890063351, -0.14377922602722601, 0.078423956469422018, 0.057160335632579121, 0.1758859767936318, 0.12307402164456618, 0.043266091906772712, -0.063956239140019508, -0.21076698450897985, -0.34935812814353639, -0.45355470119438479, -0.43051154144618958, -0.51722906014628944, 0.15560186997025938, 0.02376900091427863, 0.1937206672732435, 0.15239448898889424, 0.022154991408310377, -0.043654486348542874, 0.15332031366958024, 0.10899116252889095, 0.19918107544886771, 0.054458369263066984, 0.11037652666672665, -0.75577199845745313, -0.47055407026248092, 0.80123249840126409, -1.2482190395029769, 0.074470758177730939, -1.2888777036518133, 0.30388868428682408, -1.4090153420038936, 0.77766311308853486, 0.086142375874943722, -0.20226767939927121, 0.18630425886721039, 0.15988330961851227, 0.14658150061325848, 0.088626234734120221, -0.060936923017667906, -0.19166419512490052, -0.26718136730575015, -0.22324548858793625, -0.16530147427910316, 0.069496867327866751, 0.16777938781510615, 0.26499922119954472, 0.14731166465811715, 0.091949150799082904, -0.040407604161453857, 0.07945189422396387, 0.049772201284700876, 0.13497815095122018, 0.10690976611351256, 0.16302064990391416, -0.30097794980879811, -0.19463913189037552, 0.25555617026146371, -0.80860396068754914, 0.16834959597661264, -0.46372256496135666, 0.083846748018996636, -0.7408363358679616, 0.23800865420159104, -0.76498375921014483, 0.20487634390276502, 0.11189726101521356, 0.11490318286815554, 0.17912928670535491, 0.14915894602693922, 0.077221658551183403, -0.063285121982436826, -0.18464942262010964, -0.19523869399925042, -0.15943316022756263, 0.01037801154059478, 0.09287604949663382, -0.015525256999667727, -0.050942947866376924, 0.0024400785561280833, -0.06382739270218761, 0.0012785498169549515, -0.066205982459197621, 0.075485820780458165, 0.079593566039511995, 0.1573350303450895, -0.071305584258897034, -0.047673421329955501, 0.13017614472815361, -0.48380889084851014, 0.05970121628039872, -0.27168929694820493, 0.1532572256675474, -0.37001387119519286, 0.035028927485033709, -0.18488826926782947, 0.055072553244877351, -0.40553391670472771, 0.15195750501651886, 0.21530290856092779, 0.16972427151540051, 0.11246394245376394, 0.042439064973571687, -0.07466463078185144, -0.159079577204887, -0.16547295203388937, -0.044701593559494109, 0.0039295920043746378, 0.011913601076878649, -0.18285689914377989, -0.10407564589066134, -0.23555546334259525, -0.10684932807789757, -0.22814734409189458, -0.063485173587349258, 0.035231426543992535, 0.15565555640319106, 0.02427504006057549, 0.033548122233326282, 0.10951630791791925, -0.25399382193861758, -0.0080449729031340841, -0.10021979561886678, 0.093791526140288961, -0.14792025966567243, 0.11264804489487593, -0.099273739290014568, 0.035914795239387988, -0.16547166809914049, -0.048509383863747611, -0.02664719406199113, 0.19237746196754099, 0.16647262366552154, 0.090094395501221572, 0.0039906200934596249, -0.098122755354853444, -0.1366171762599761, -0.038465872090180271, -0.063107949091446403, -0.10353827744152141, -0.2367561033844639, -0.39107302489471046, -0.4993413545698927, -0.29428485506006952, -0.35686545376814893, -0.26477441372853422, -0.21957600785318782, 0.10231188547257429, -0.038641816351083894, 0.10657849804550705, 0.15251693996845886, -0.083479299629685352, -0.043761647973181945, -0.1569328811173635, -0.025732944443590922, -0.13459193297646968, 0.2054603592637699, -0.076343100588013288, -0.012037966893705314, -0.10529567915508307, -0.096096791377201118, -0.048761457544342024, 0.10089614648569795, 0.095268264507682401, 0.14454547803178003, 0.11929552297831249, 0.062089135258225896, -0.072827255227657164, 0.017260207971635976, 0.0029845825084865137, -0.15373347718831978, -0.29387297238508514, -0.24699055160378858, -0.30677005181583211, -0.085911462662095636, -0.30581901500051922, -0.25423051056794305, -0.38754570378991321, -0.22411418130009808, -0.35375643597530843, 0.055146819140510416, 0.26281055993233621, -0.019211160710796572, 0.021646587213104713, -0.10241250352514171, -0.078770189361086179, -0.14656177726216213, 0.06994065617155748, 0.0021974481246246629, 0.002747867764730379, -0.10920785431913549, -0.11314958231531393, -0.099324973861329505, -0.024565620333014605, 0.11510723629505158, 0.15553553242798132, 0.1412615980599799, 0.1780996924064803, 0.24246246508902825, 0.30032839778422799, -0.043896728635589223, -0.041656273199716369, -0.30818729583490695, -0.22221882441369153, -0.34490255487274391, 0.23901621577757753, 0.040068334939789153, 0.22067925957890494, 0.059139017637366363, 0.0013771249391555669, -0.82932242402320322, -0.24248276916859551, 0.2220707306234446, -0.05187740979518813, 0.089711600875348782, -0.064944626021130755, -0.023033890597353883, -0.20608602740074317, -0.0068114015565181726, -0.056907964871020393, -0.0092446671691440163, 0.031869312918224861, -0.031719635617030659, -0.20558191361174011, -0.077733070765210868, -0.042569660785952056, 0.14494733127713064, 0.20572232603908944, 0.22275756609764749, 0.23658495591681294, 0.34021620386249518, 0.53178862418379935, 0.23679486634227159, -0.4620869054676115, -0.17111767858219507, -0.23623861798820922, 0.12877260884038208, 0.69010247253016699, 1.1175823884838996, 1.0747715799508781, 1.1957108711898863, -0.10839909954530788, -0.2263553910474054, -0.048917424863973949, -0.49291782026044056, 0.072335854201080588, 0.024371827339035017, -0.013685307258532536, -0.20736891492564691, -0.14292038700683141, -0.20318736992960118, 0.09738711454895084, -0.11463822183415182, 0.066392764898706338, -0.087906272749638886, -0.24871001073084317, -0.24716165183994784, 0.070818253229819925, 0.19552146385155544, 0.2667283201489779, 0.30564289179799059, 0.3354280050625516, 0.60106744735122974, 0.69730138646619622, 0.10657803670352221, -0.26597013575809814, -0.51663850122651156, 0.20353239335470413, 0.46612626705839527, 0.88179245186152555, 1.0591766157447142, 1.8155993788135401, 1.1338250057069015, 1.3365493166122591, 0.60566186502852537, -1.1434123865171271, -0.24404333943615408, -0.11036542395783162, 0.011397583876152599, -0.13747467644410641, -0.17507285276321236, -0.22950814490677501, -0.12352253662709391, -0.15992872157622667, 0.070404205425931329, -0.14847974776257952, -0.11769772351689584, -0.0043287613096706385, -0.01570204240179094, 0.0089756748769360926, 0.29363217368931571, 0.24220392956380496, 0.33450192958610414, 0.50614275317411517, 0.79818525223669812, 0.56259342729123207, -0.51618519149943443, -1.105930659634037, -0.36434278784189805, -0.25061306335904537, 0.49191197464707476, -0.55042973333868328, 0.33201271766592266, -1.2534804319442674, 1.8545676439466376, 2.7705947112358396, -0.65727781823869658, -0.44769863619371453, -0.72484131059467105, -0.19352569039907921, -0.06091583791907576, -0.075793328725683354, -0.18362798444373865, -0.10847923107858076, -0.4971843279255892, -0.26497067484015147, 0.12790446781792791, 0.30557543109041591, 0.010766872600696822, -0.17753959085018053, -0.04814578185160058, 0.007404938424997059, 0.19111004712688809, 0.26024004494781799, 0.38143909917016594, 0.62956014799324345, 0.63152717265524094, 0.50640814880924734, -0.28203486996779953, -1.481298784591466, -0.34444399750770943, 0.33223282315731145, -1.0635980503754157, -0.74327388531916416, -1.4203730522981521, -0.29452233058397176, 0.31054780878465887, 0.4425799571579111, 0.58851980199844711, -0.60890484535134703, -0.54015192300944392, -0.38389116569218451, 0.02852961095842791, -0.22272564345822352, 0.021237638325899852, -0.387784194546968, 0.021133356436687064, -0.70317235961626101, 0.56727270004732033, -0.037766801642101375, 0.06793846381428098, -0.54370981678606445, -0.070485969368251203, -0.02650526715020184, 0.26990648117357596, 0.2830749986604576, 0.32017932907632218, 0.48573187070286, 0.4386512882244713, 0.87865322612556163, -0.72553155869355623, 0.48040960369964236, -0.14688086473049994, -0.86361123753821778, 0.034562131954569586, -1.3606412475160699, -0.095775027331125598, -1.0236217631956728, -2.0668504460046977, 0.91046045526812469, 0.047649837051312177, -0.29915225181420874, -0.44192642538866073, -0.16705843369252318, -0.088135922159037081, -0.1055833934917848, -0.15296582393740679, -0.37808744252168625, -0.27971623762694292, 0.42301568032977793, 0.28749903579456393, 0.42643957249522774, -0.63308466384561635, -0.39563373733228846, -0.2486369819601503, 0.37438920939802361, 0.1651784382566546, 0.21235668984583117, 0.050683391882319691, 0.2929172313060312, 0.12600631152419559, 0.37020072416606892, -0.11933000679563836, -0.077869450453742661, 0.40059289889731342, 0.33639737707180556, -1.4411916754233998], "height": 32, "width": 32, "top": {"real": [6327.0074679827858, 12787.448651417644, -8580.1357782693794, -6564.4296251926053, 1064.8969915597925, -1614.7722808334524, -60.048119421842308, -221.12968025641604, -12.445754708845573, -62.940587713583341, 23.820578453457951, -28.525301748952867, 37.864280125154337, -85.844490405627994, 27.493814969643068, -48.151974443362732, 41.399660199307142, -48.151974443360345, 27.493814969642067, -85.844490405627866, 37.864280125154032, -28.525301748952689, 23.820578453458509, -62.940587713584037, -12.445754708845573, -221.12968025641445, -60.048119421842806, -1614.7722808334493, 1064.8969915597932, -6564.4296251926125, -8580.1357782693849, 12787.448651417659, -17213.595236981997, 544.07302669902936, -10829.417860945707, -5794.1219231083051, -1103.8692252742003, -1120.3666800727638, -106.61462896985108, -117.89522553300084, -39.308430893651213, -26.960578056534288, -11.025047422886351, -2.1850531713376355, -7.2952697383625882, -14.472211987615461, -14.23279965231924, 3.5205960317557246, -19.629747734378611, -14.131264529804582, -22.272265924680642, -17.401343902608428, -1.7416101263052008, 0.34733716754184635, -16.60063574625633, -32.47390190466588, -64.31481385096977, -84.034753639409303, -258.32278440436392, -714.30117659130417, 113.1165060355228, -2806.7151144576005, -10882.302955904926, -27817.100568112113, -27745.708494118069, -22841.420845008699, -8553.5264674688351, -4085.033989629072, -1748.9561003744159, -545.33543865452532, -146.45498278007594, -40.805440061891282, -50.060212982262158, -16.730489448484349, -23.62313915150434, 0.36714120063527833, -7.1882909197809326, 6.4142604093545659, -14.833291452397129, 3.4669422259569203, -16.161012153790594, -21.283130178011959, -13.244388638951026, -1.2754674400463286, 0.49184661171499777, -11.472717532816317, -17.530689402649369, -29.523444921912187, -36.983418680226812, -91.597351236026796, -144.17841735614732, -486.41083616652156, -239.1204206916959, 139.00893917602821, -7153.1577950265046, -18373.599790964956, -8415.9227315857388, -6467.5546779833858, -5469.6585676007244, -2789.7914675082047, -821.48267137489768, -514.71454200305448, -122.6780822655295, -35.554593402553344, -32.905184819607861, -16.875977506192179, -9.1243948190511883, -4.1802620318744168, 0.78316084097390803, -0.82473624462572925, -6.0277039354546789, -6.8502427898926159, -6.8858593350514123, -14.551844205848809, -8.4592989146599145, -1.9950695133570466, 1.5394495557011576, -9.232834653609892, -11.961754213185406, -20.558561534905806, -38.14482496917136, -57.160869905334089, -157.5554953403844, -411.20082432868492, -110.11930975746139, -463.64281724562954, -1210.9754315912523, -4810.1023052011815, -2624.9192342435695, -4516.897453633137, -3582.982491421491, -2119.7499173636902, -750.23472664477231, -159.86878938364373, -96.805728919101412, -16.935754226078267, -22.178078450833883, -9.0092227225989756, -5.4454852920506616, -2.5268882566736748, -2.6219523498948458, -2.0534180946008638, -6.4643064804236223, -8.7214065300148622, -7.0127031874072934, -7.7847619200438425, -0.78121801903708787, -0.11031060182209244, -3.4318959866990171, -4.568083040130742, -9.3160200641569482, -15.87900618307528, -23.382201868158507, -37.989404730637517, -121.95929165874175, -201.43809047465336, -567.30372690179445, -610.36312853782476, -1937.2376557226703, -2694.0942422163184, -3481.5355975420739, -3026.0617140239533, -2506.4916383888308, -1355.1576691715959, -468.0599091435185, -33.397814851582154, 18.545331481583979, -29.439943085485314, -11.922803368680741, -7.7174579848963036, -2.7888603302094261, -2.0545006330536419, -1.0224301618297007, -2.5290059501700397, -2.6199349314032494, -4.8828731457816774, -5.6140621268617785, -3.6146413255535528, -0.82967484502750599, -1.1631844631411599, -1.1736245693095499, -3.7427284368684024, -6.6178270656423424, -8.8039498431578131, -13.481606634721874, -26.682014936266803, -32.97802138197703, -190.87870662958537, -335.07671606134033, -876.95522185704988, -1470.8308615232113, -2413.1056871868645, -2329.4915082040015, -2151.9265781424592, -1438.5211303616525, -713.08152385953474, -130.56706965452, -28.862997769748436, 1.7115783105395304, -2.3079138134833856, -11.650001060370576, -2.1701895043697319, -1.7033561682904694, 0.11972955037883275, -0.5790790812653609, -0.45369004788794814, -4.0202584980239315, -2.2874964618085922, -4.0439693313887179, -2.7144116941037568, -0.78431497120559035, -0.63377862148075614, -0.83542887465582039, -2.5417481635934216, -3.4497798452119151, -5.2729359146939485, -9.3325337684514853, -2.3480777220919959, -30.998957295773433, -62.717085150656224, -255.55592612696097, -515.55403813953228, -1077.4118984564277, -1959.4321649356718, -996.17063256475433, -957.58296584245431, -646.48124546709357, -264.63237090275328, -80.518340163143705, 1.3577640380308329, -10.209813205700032, -1.6455702447745924, 0.041564105453412992, -2.1114336115237951, 0.094770419949782359, -0.59915948675603581, -0.042806075220647909, -1.4583889440112543, -1.0515624140257636, -2.6242610725724358, -1.9962698899057387, -2.6373580888952279, -0.39414295752688255, -0.80952098820947449, -0.80935504019781357, -1.8315981280862581, -2.3957781219887311, -3.7775973636187947, -1.2762438485839673, -5.8736187047076989, -1.7472058496177199, -39.26171578474986, -103.67045476207605, -221.97901463800724, -486.52939506954147, -770.86385201715837, -324.53401456188737, -340.76561269504333, -198.94197841299811, -116.98390967138687, -32.552002455687223, -2.8647426984892057, 3.4567160839115996, -2.3310052693731942, -0.31982529888503974, 1.0641170511547686, -1.1693860348016079, 0.40700212695667498, -1.0092937918537679, -0.034985154955161232, -1.1746986822405672, -2.2583319261560799, -2.3742700657768196, -2.0999344700147127, -1.3199204785521135, -0.00726121067659493, -0.21997349279604808, -1.6799449140427334, -2.573626149077699, -1.2749029797075579, -2.9722975102588807, 1.1149262389696764, -3.7825574143705132, -8.8112838568342262, -31.380656277665604, -80.815706494855448, -151.2920880811312, -231.73584987382, -71.719340760371239, -68.364029336894532, -82.304975309081328, -39.769136567332623, -12.693902195545155, 1.9842035662706039, -0.48894354077452629, 1.5673786316984977, -1.1528550937233366, -0.047905313676576855, 1.0650111541187497, -1.4670778837513785, 0.90154843827712783, -2.0067122081730173, -1.7180755446076856, -3.3483434682482898, -1.6317774897434336, -3.4708804264384949, -1.4694727348954268, -0.10257790692130256, -0.92889636381571783, -2.5451436559397584, -1.1303995452351283, -2.8871606625355204, -0.49805533537312113, -0.43871934228158149, 0.74909741969065036, 1.439791857975526, -13.026034555561825, -28.484337861014367, -36.895178406765908, -72.231190035838395, -26.809485534889244, -56.942733796520351, -36.280107429422031, -18.099079313675151, -3.2627913287121535, -0.27065475805048173, 1.1511525903295716, -0.1898870363801409, 0.97087220888025372, -1.1326412370424179, 0.19981301766363371, 2.3076802702931332, -2.5498407756566031, -0.55532623625253041, -4.5067864970337066, -3.1791416623515807, -4.4078100962970312, -3.3125592316970569, -0.86117313527739991, -1.5893272599995412, -2.966763176308445, -0.10263295233944174, -3.9090515756474873, -1.1809123553765606, -1.2038275093350199, 0.28624805719871799, 1.8536817132476691, 0.80720029239974944, -5.5769148242308617, -9.876045833403877, -26.746787584056406, -25.156847324745275, -59.858730768373135, -40.909481344220289, -31.301364606284917, -12.992896188265529, -3.3747756869001191, 2.018530777741562, -0.2260612077575776, 0.20594677690592242, -0.32403225688254972, 2.6342781156582364, -0.20787587650621794, 0.020021422684959465, 2.4168167113542753, -5.4711139391067487, -4.0162352275216753, -8.2614925402209369, -5.3597212899808353, -3.5618518868969149, -2.433206189931258, -3.72732750056158, 1.9796458045760836, -5.6194552495050392, -1.3533084788499918, -2.5312603283156982, -0.73577412998862446, 1.0534870856485199, 1.9641802030644013, 2.686394865901931, 0.15179077698100021, -13.505006616733285, -24.652619118640313, -46.540772018331772, -65.573862825648789, -62.198686950727755, -33.367481359773279, -15.462856466374785, 1.2764820376727259, 0.65238351077916601, 0.11426292455840634, -0.28574074744887695, 1.1004477970868616, 1.5327652909312646, 4.9327112002405187, -2.7503236900904393, -3.9754701567757182, -0.83451213573114202, -9.6859551174636032, -6.0006304228248837, -7.1318253263858704, -4.9856424549649869, -8.2198692881497966, 4.4471445589999075, -5.9124192136821829, -1.066228559248251, -4.5320020023541776, -2.2691032318803455, -0.29595929991781228, 0.95069246916243999, 2.9800299945570488, 4.1754888650460389, -3.9368165535394617, -16.752391317324459, -41.141712299851427, -60.643902512124761, -84.343164732157774, -54.074972316673119, -32.147836838432397, -3.8004568077904586, 0.60906370022349199, 0.96326571070502964, 0.12072436644318041, 0.86319547550846831, 1.1644302862121145, 2.3252453178811563, 0.4766448788746544, 4.2895158747096422, -8.0099398507205173, -11.261414242734404, -3.7327404661807368, -6.1601929909303115, -5.7506657339441176, -11.117751544185502, 4.3054115148233869, -8.2116085703478756, 1.5649721654809707, -5.692142411938895, -3.4949086671880489, -2.1697325908161673, -0.69541562648800392, 0.8151347377023449, 3.4913843961972888, 3.7128410928274609, -0.28461904914440361, -14.694115553922865, -40.784171797776878, -76.129097997930202, -42.446544168449137, -35.20072244442639, -7.6849685309209681, -0.86778876706269803, 1.7056687310546401, 0.51974952798739482, 1.1316759112288739, 1.2127800605381189, 2.2221214293211369, 0.8361427470500683, 0.099546606452808759, -3.1066939683116588, 2.4769973579809332, -11.745899626687375, -12.566716735550758, -4.7843722525441121, -8.6040588700645948, 5.3511034027408684, -15.66778198729919, 3.3911987930307137, -8.1259457896856588, -4.7980365887404339, -2.1471667518999764, -2.5087241573260579, -0.86097656871748596, 0.80139543608377761, 2.3927476350849495, 5.6561627225456101, 4.6250195847240709, -10.360873371570062, -32.183766942993408, -40.627648783043085, -17.95302681367091, 6.1859416444601623, 4.4451711191118575, 6.2133312512745986, 3.3587096432927233, 0.75569650762962715, 0.82429958430179118, 1.6770717987294814, 1.4875078542842981, 0.57506461301516931, -1.8616897808321822, -4.2923314886050523, -7.0161307894352376, 4.8526265686569623, -11.625945416074895, -21.07467865045702, 3.1434482955713889, -12.714871655934468, 2.366640700254353, -11.173661544979289, -6.475058123170597, -2.9596627855411386, -2.7630790459355632, -1.4122750319888788, -0.55802289067111654, 0.40509965601619852, 3.0549503122920978, 4.5519407117959005, 5.5044865911271481, 2.884340309618711, -2.0306279174386637, -12.290279126535582, 30.20785901727286, 12.45094729737127, 8.7415961086741216, 10.304484703022888, 4.2310313602871554, 1.8368366318675937, 1.9910666876255574, 0.99156303676388502, 0.31747940339256003, -0.60245605946597469, -1.5054245711483085, -5.2047479791857727, -4.6985425899868174, -6.4406057657264579, -2.4317965202694718, -0.31806110286890099, -34.614374984304249, -0.31806110286891309, -2.431796520269478, -6.440605765726465, -4.6985425899868325, -5.2047479791858065, -1.5054245711483181, -0.60245605946599046, 0.31747940339256003, 0.99156303676388657, 1.9910666876255707, 1.836836631867627, 4.2310313602871847, 10.304484703022929, 8.7415961086741465, 12.450947297371268, -17.953026813669211, -12.290279126534008, -2.0306279174386814, 2.8843403096185996, 5.504486591127181, 4.5519407117958774, 3.0549503122920303, 0.40509965601617365, -0.55802289067113198, -1.4122750319888866, -2.7630790459355659, -2.9596627855411324, -6.4750581231706095, -11.173661544979305, 2.3666407002543774, -12.714871655934495, 3.1434482955713383, -21.074678650457034, -11.625945416074918, 4.852626568656909, -7.0161307894352527, -4.2923314886050816, -1.8616897808321802, 0.57506461301515821, 1.4875078542843145, 1.6770717987294999, 0.82429958430183881, 0.75569650762972518, 3.3587096432928165, 6.2133312512752603, 4.4451711191129393, 6.1859416444606961, -42.446544168447332, -40.62764878304236, -32.183766942992868, -10.360873371569919, 4.6250195847241038, 5.6561627225455622, 2.3927476350849091, 0.80139543608374619, -0.86097656871749839, -2.5087241573260646, -2.1471667518999924, -4.7980365887404339, -8.1259457896856819, 3.3911987930307128, -15.667781987299181, 5.3511034027407849, -8.6040588700645682, -4.7843722525441965, -12.56671673555077, -11.745899626687423, 2.4769973579809226, -3.1066939683116837, 0.09954660645280862, 0.83614274705006209, 2.2221214293211515, 1.2127800605381271, 1.1316759112289154, 0.51974952798746277, 1.7056687310547414, -0.86778876706233132, -7.6849685309204459, -35.200722444424429, -84.343164732156296, -76.12909799792908, -40.784171797776416, -14.694115553922703, -0.28461904914431585, 3.7128410928274769, 3.4913843961972515, 0.81513473770227451, -0.69541562648801736, -2.1697325908161949, -3.4949086671880569, -5.6921424119389146, 1.5649721654809774, -8.2116085703478952, 4.3054115148233585, -11.117751544185527, -5.7506657339441034, -6.1601929909302848, -3.732740466180779, -11.261414242734437, -8.0099398507205191, 4.2895158747096467, 0.47664487887466933, 2.3252453178811465, 1.164430286212123, 0.86319547550848241, 0.12072436644322472, 0.96326571070516387, 0.60906370022364942, -3.8004568077899661, -32.147836838431445, -54.074972316672174, -65.573862825648646, -60.643902512124406, -41.141712299851122, -16.752391317324395, -3.9368165535393809, 4.1754888650460389, 2.9800299945570101, 0.95069246916240213, -0.29595929991782399, -2.2691032318803499, -4.5320020023541936, -1.0662285592482377, -5.9124192136821918, 4.4471445589998728, -8.2198692881497806, -4.9856424549649665, -7.1318253263858535, -6.0006304228248819, -9.685955117463573, -0.83451213573115879, -3.9754701567757316, -2.7503236900904167, 4.9327112002405231, 1.5327652909312637, 1.1004477970868634, -0.28574074744888678, 0.11426292455844969, 0.65238351077926593, 1.2764820376728394, -15.462856466374522, -33.367481359773095, -62.198686950727343, -59.858730768372553, -46.54077201833131, -24.652619118640224, -13.505006616733089, 0.15179077698109764, 2.6863948659019079, 1.9641802030643754, 1.0534870856484704, -0.73577412998863267, -2.5312603283157147, -1.3533084788499781, -5.6194552495050258, 1.979645804576063, -3.7273275005615623, -2.4332061899312447, -3.561851886896878, -5.3597212899808397, -8.2614925402208979, -4.0162352275216939, -5.471113939106738, 2.4168167113542718, 0.020021422684962133, -0.20787587650621156, 2.6342781156582284, -0.32403225688256132, 0.20594677690590416, -0.22606120775752253, 2.01853077774176, -3.3747756868998953, -12.992896188264998, -31.301364606284267, -40.909481344219984, -26.809485534888527, -25.156847324744817, -26.746787584056008, -9.8760458334037811, -5.576914824230828, 0.80720029239978064, 1.8536817132476344, 0.28624805719868845, -1.2038275093350328, -1.1809123553765482, -3.9090515756474855, -0.10263295233946072, -2.9667631763084428, -1.589327259999552, -0.86117313527738326, -3.3125592316970565, -4.4078100962970224, -3.179141662351582, -4.5067864970337235, -0.5553262362525373, -2.5498407756565928, 2.3076802702931301, 0.19981301766362944, -1.1326412370424508, 0.9708722088802425, -0.18988703638015317, 1.1511525903296145, -0.27065475805034434, -3.2627913287119403, -18.099079313674878, -36.280107429421662, -56.94273379651932, -71.719340760368894, -72.231190035836406, -36.895178406764366, -28.484337861013827, -13.026034555561461, 1.4397918579755735, 0.74909741969062116, -0.43871934228165033, -0.49805533537312413, -2.8871606625355364, -1.1303995452351521, -2.5451436559397544, -0.92889636381573193, -0.10257790692130979, -1.4694727348954375, -3.47088042643851, -1.6317774897434372, -3.3483434682483164, -1.7180755446077083, -2.0067122081730231, 0.90154843827713693, -1.4670778837513718, 1.0650111541187297, -0.047905313676639846, -1.1528550937233546, 1.5673786316984972, -0.48894354077444435, 1.9842035662707793, -12.693902195544938, -39.769136567331614, -82.304975309079381, -68.364029336892514, -324.53401456188737, -231.73584987382011, -151.29208808113108, -80.815706494855604, -31.380656277665615, -8.8112838568343115, -3.7825574143705687, 1.1149262389696379, -2.9722975102588807, -1.2749029797076188, -2.5736261490777141, -1.679944914042758, -0.21997349279605258, -0.0072612106766204651, -1.3199204785521295, -2.0999344700147602, -2.3742700657768196, -2.2583319261561177, -1.1746986822405874, -0.034985154955167297, -1.0092937918537661, 0.40700212695667115, -1.169386034801607, 1.0641170511547318, -0.31982529888503974, -2.3310052693731778, 3.4567160839116182, -2.8647426984891329, -32.552002455687209, -116.98390967138671, -198.94197841299811, -340.7656126950431, -996.17063256474341, -770.86385201715302, -486.52939506953675, -221.97901463800537, -103.67045476207502, -39.261715784749725, -1.7472058496176825, -5.8736187047078241, -1.2762438485839995, -3.7775973636188502, -2.3957781219887693, -1.831598128086283, -0.80935504019782156, -0.8095209882094927, -0.39414295752690548, -2.6373580888952821, -1.9962698899057312, -2.6242610725724247, -1.0515624140257607, -1.4583889440112356, -0.042806075220638506, -0.5991594867560418, 0.09477041994976447, -2.1114336115238088, 0.041564105453392161, -1.645570244774569, -10.20981320570001, 1.3577640380311771, -80.518340163142966, -264.63237090274987, -646.481245467088, -957.58296584244442, -2329.4915082039943, -1959.4321649356684, -1077.4118984564245, -515.55403813953228, -255.55592612696057, -62.717085150656416, -30.998957295773568, -2.3480777220922762, -9.3325337684514711, -5.2729359146941093, -3.4497798452119492, -2.5417481635934838, -0.83542887465582316, -0.63377862148080311, -0.78431497120560345, -2.7144116941037488, -4.0439693313887171, -2.2874964618085487, -4.0202584980238978, -0.45369004788789491, -0.57907908126534979, 0.1197295503788481, -1.7033561682904654, -2.1701895043697248, -11.650001060370577, -2.3079138134832218, 1.7115783105395919, -28.862997769747821, -130.56706965451951, -713.08152385953122, -1438.5211303616488, -2151.926578142451, -3481.535597542063, -2413.1056871868605, -1470.8308615232086, -876.95522185705033, -335.07671606133954, -190.87870662958562, -32.978021381977214, -26.682014936266988, -13.48160663472185, -8.8039498431579233, -6.6178270656423885, -3.7427284368684086, -1.1736245693095617, -1.1631844631411714, -0.82967484502748334, -3.6146413255534444, -5.614062126861783, -4.8828731457814465, -2.619934931403141, -2.5290059501698665, -1.0224301618296836, -2.0545006330535864, -2.7888603302093888, -7.7174579848962566, -11.922803368680738, -29.439943085484902, 18.545331481584157, -33.397814851580577, -468.05990914351662, -1355.1576691715877, -2506.491638388823, -3026.0617140239456, -2624.9192342435676, -2694.0942422163216, -1937.2376557226696, -610.36312853782647, -567.30372690179479, -201.43809047465459, -121.9592916587422, -37.989404730638029, -23.382201868158511, -15.879006183075566, -9.3160200641569801, -4.5680830401307739, -3.4318959866990371, -0.11031060182210287, -0.78121801903713106, -7.7847619200437572, -7.0127031874072889, -8.7214065300146117, -6.4643064804234296, -2.0534180946007163, -2.6219523498948099, -2.5268882566736321, -5.4454852920506234, -9.009222722598933, -22.178078450833912, -16.935754226078039, -96.805728919100844, -159.86878938364163, -750.23472664477129, -2119.7499173636857, -3582.9824914214882, -4516.8974536331352, -8415.9227315857224, -4810.1023052011687, -1210.9754315912437, -463.64281724563131, -110.11930975746253, -411.20082432868651, -157.55549534038482, -57.160869905334962, -38.144824969171324, -20.558561534906186, -11.961754213185522, -9.2328346536099755, 1.5394495557011911, -1.9950695133570939, -8.4592989146598381, -14.5518442058491, -6.8858593350514123, -6.8502427898925227, -6.0277039354547091, -0.82473624462553619, 0.7831608409739802, -4.1802620318743058, -9.1243948190510782, -16.875977506192154, -32.905184819607861, -35.55459340255284, -122.67808226552924, -514.71454200305027, -821.48267137489574, -2789.7914675081947, -5469.6585676007144, -6467.5546779833749, -27745.708494118051, -18373.599790964956, -7153.1577950265018, 139.0089391760178, -239.12042069169632, -486.41083616652713, -144.17841735614894, -91.597351236029397, -36.983418680226798, -29.52344492191359, -17.530689402649603, -11.472717532816715, 0.49184661171508715, -1.2754674400467305, -13.244388638951113, -21.283130178012485, -16.161012153790605, 3.4669422259560148, -14.833291452397047, 6.4142604093547071, -7.1882909197807701, 0.36714120063543831, -23.623139151504223, -16.730489448483766, -50.060212982262158, -40.805440061890259, -146.45498278007528, -545.33543865451918, -1748.9561003744159, -4085.0339896290589, -8553.5264674688297, -22841.42084500867, -17213.595236981979, -27817.100568112102, -10882.302955904914, -2806.7151144576032, 113.11650603552238, -714.3011765913102, -258.32278440436545, -84.034753639410752, -64.314813850969756, -32.473901904667066, -16.600635746256838, 0.34733716754136956, -1.7416101263050108, -17.401343902608698, -22.272265924680649, -14.131264529805808, -19.629747734378629, 3.5205960317556224, -14.232799652319452, -14.472211987615371, -7.2952697383626051, -2.1850531713375001, -11.025047422886134, -26.960578056534132, -39.308430893651277, -117.89522553299832, -106.61462896985047, -1120.3666800727583, -1103.8692252741987, -5794.1219231082987, -10829.417860945698, 544.0730266990613], "imag": [0.0, 39515.399789254865, -15382.502646592808, -2918.5213899411633, 1191.3150820979915, -345.80526987183259, -91.963497561164317, -108.06830259743781, -23.431327605037779, -16.610810822958811, -2.2742971526945692, -21.94540337447156, -21.272420402577179, -19.057735286082725, 52.376313120912741, 8.8521975026144215, 0.0, -8.8521975026153612, -52.376313120912044, 19.05773528608244, 21.272420402577229, 21.945403374471546, 2.2742971526945164, 16.61081082295847, 23.431327605037779, 108.06830259743377, 91.963497561167031, 345.80526987182992, -1191.3150820979911, 2918.521389941166, 15382.502646592799, -39515.399789254865, 16946.371225111547, 11272.400870696511, -1841.1630655848539, -1171.8258105148211, 109.26740577843201, 109.43221584454041, -206.10314812122442, -42.28539321936541, -65.249140172269108, -8.1960517939269071, -25.365956236318521, -12.488977779350739, -20.777136166915728, 10.939736173955986, 13.738393384077787, 18.982014759380824, -12.446635345723246, -8.8305069728752752, -14.016932746731957, 10.847407040619599, 18.949463773817904, 9.5491377266946333, 9.9160507928571491, 10.925693066643545, 29.597688207921799, 110.70270321520633, 119.64570653889598, -28.406380421737747, -832.38513354395832, 3868.3607546707453, 4625.1554139998552, 1299.266384884336, 11955.369188695569, -346.60702359457474, 2363.0865643328202, 202.32180136864196, -295.01076934839688, 277.44131247872377, -208.25011765184485, -68.536452109290167, -63.125416501511275, -30.589238972650307, -20.382558491276782, -12.902690844348575, 1.0182743702449766, 2.2874141051372452, 14.546101055958879, -1.6150967690531892, 0.87386354462209248, -13.38677146472326, 17.407635810762027, -0.60680891034428563, 12.950261495215786, -0.86756998254740836, 13.062832664150156, -2.136096668655318, 35.880437410281019, 66.543583791324124, 106.41042136702396, -154.61175360170367, -142.43138300706809, -204.04622542072212, 637.22651951043269, 14506.214536688893, 2087.241286517743, 2934.8250858413794, -540.57641351250766, 271.57679130866723, 537.49064424335052, 71.610686309506633, -82.566517160157844, -74.489552074506761, -57.142723569279489, -25.254222647826957, -12.544755138678799, -6.7708773337384622, -0.20770764759918051, 1.1194026551973495, 3.6890839097334878, -0.20461611038051675, -5.3737401851528386, -2.2636764639022453, 4.9832600695280087, 10.10754418440443, 2.544597152813191, 5.0447569422202223, 2.8189903732027388, 2.9920185974422946, 11.529259915073515, 53.348281783168488, 50.577306165813702, 119.57179721034167, -786.94393395792213, -901.86972095078409, 1074.2990124084633, -676.63432998511303, -1814.4894639161919, -1058.8873901892568, 349.63255334068964, 268.2484621242063, 380.4257804868011, 97.912476061205965, -71.661879904778658, -72.736940205390937, -37.286147332183255, -18.720186517652358, -9.5628328869678363, -8.2634536416369482, -0.19099280196160889, 0.28567330333833429, 0.071309220847330906, 1.775656575203522, 2.1386941287610384, -1.0231985430634363, 2.4480270844924137, 1.9878096211738889, 4.1535160740617325, 3.5521022242715827, 3.6019889835361654, -0.38778650675348675, 4.5934633035837322, 12.369503115173785, 94.400929285346393, -77.039516694541632, -317.96189631496532, -464.93982303098932, -1928.5847150696218, -1707.435501473529, -1583.8824730375936, -833.14970460709912, -242.65732454029356, 317.34459246185179, 69.371971565542168, -13.453603503282046, -58.678662168810511, -40.697925231357267, -20.16679423742438, -13.090865763493372, -7.6639325468155182, -3.9818184076004774, -1.1480426683317846, 3.4782462601712654, 3.3535274581593777, 1.1304420427842019, -0.78798091056569108, 1.3645618319176871, -1.2482103302418359, 2.6548043613808709, 2.0455422971726147, 3.3794759874708187, 0.28270713586437868, -0.67396364813071052, -2.4299863370655714, 20.048308524208917, 2.7323622379310764, 23.350930435926227, 19.428895667387472, -542.55543458513841, -969.72646934413422, -1881.6935115818319, -426.42488803140975, -203.12933548649579, 221.2759941924873, 75.498646737827528, 81.527654709232138, -23.228148927519229, -35.095892631083601, -18.113211298120703, -10.84781416147373, -7.2190125620904571, -4.6279710639736358, -2.4591204313834965, 0.58454842364750592, 0.9547463106026568, 2.4151219203224654, 0.35567626693527166, -0.34443735246553581, 0.10021534895843374, 0.62975751586612239, 0.075739057981468363, 2.5606009702294732, 0.33652577635989583, 1.0967346902818458, -2.0000864020848046, 1.451434335785166, -7.1043203483245918, 10.682010714917043, 63.119371696969466, -0.21915489532378957, -18.829394319190524, -151.38668525810633, -53.557614725234764, 319.98109193688521, 339.65571159114432, 137.34002448547074, 89.166061674211335, 5.5395407010458548, 3.1394961579749192, -15.790153056071546, -10.068399035154068, -5.990189510359599, -3.8184783724455591, -3.031311871551015, -1.7050468427375405, 0.29997994657424726, 1.0820090369221063, 1.186145834536469, 1.1035985643781976, -0.026686439684276222, -0.021861987241408937, -0.71069133238487792, 1.5408441736495715, 0.88964240698419117, 1.3226298476418663, -0.31341792844049038, 0.36542819178090652, -6.2831724147446009, -1.5019641274621849, 13.705577255286178, 14.08988446983488, 43.376064565045787, 54.204646228028693, 132.60719521227122, 62.432091745032416, 219.86853407025291, 133.19971898528303, 106.65962198554573, 32.509510578136677, -0.32173120932219279, -1.0472919151881939, 3.0400482541159932, -2.1373375769283371, -3.730023177462626, -3.658884080605465, -1.7579654575817145, -1.4081716233776314, 0.53847505259675243, 0.97164442609396484, 0.99557444173589427, 0.11126643861839057, 0.68343802732828218, -1.8968331338371984, 1.2370080984335854, 0.24616862215977253, 1.6488725749132858, 0.51751999590405229, 1.1293193043526017, -2.440127533596657, -2.1096176983280075, 0.86862341139232857, -3.2906850406637935, 10.529105330387669, 21.528607110729595, 53.864359931260871, 77.645971701549371, 219.35042572397626, 56.646725304489401, 72.597133152544643, 25.283726226154087, 9.5916945360518575, 3.5131066387109651, 7.8353712095885921, 3.7982495059531924, -0.1320379920173739, -2.6872144914167868, -3.9932995411602161, -3.2496097582258021, 0.34753129477634409, 0.51108563039490118, 1.5779290838292743, -0.32098662744154494, 0.85050436376706862, -0.95269533500567993, 1.2514468795435802, -0.80232514327112969, 0.9402998385797201, 0.083484537520043395, 2.9970442493787712, -1.3537504507858638, -0.49606931853204594, -0.40634793694041288, -3.2049709482452999, -0.90733505109423707, 2.4367446291982722, 14.425788131967863, 11.285893196605656, 41.383231407558732, 26.64930600456886, 3.49019745229235, -2.7278278900775033, -0.6780497768719117, 5.1972517053508067, 13.315592278460391, 8.5655340270495799, 2.7986942524720009, -1.3188307478897958, -2.9560147101356677, -3.5039625138893249, -3.2351264422424442, -1.7837191938553947, 4.1370294214229899, -0.48575889336608402, 1.5705428299533406, -0.28712766434912351, 3.0263722138472242, -1.2120159426024513, 2.3578685311613299, -0.57595981339550073, 5.0023110313506045, -1.3965343408401143, 0.57130018359553258, 0.97601539997779907, -2.5856117766371072, -1.4572886030272592, -1.697304402058907, 2.1030121366448049, -0.062995146217583023, 4.9612937854964008, -0.301351161950073, 9.1580358876592047, -26.774160253651129, -14.352442298363405, -1.3280628610826339, 14.662722831107068, 13.065714024917916, 9.0142686843581057, 2.3156698002103213, -2.1494090398895276, -3.4369058394752816, -3.2885833247038136, -1.5589303641797962, -1.2624013759488528, -1.9675974847606275, 6.6808538974877765, -1.9845455926326385, 5.8520209082315047, -0.64649466936933941, 3.7090325345277484, -1.4383189901623872, 7.6358217658181289, -2.5011370399824613, 2.0288771570791524, 3.3870661540364408, -2.1726948090167277, -1.6227282885568253, -2.3401345171732064, -1.0784728701462345, -2.329913042321786, 1.6766227476276099, -5.5208279800184377, -13.304914803323909, -34.683237345909106, -16.2208922311743, -7.5614574111833583, 14.803696771258469, 15.579335900760762, 15.161631552742872, 8.4660409590375174, 1.9742581617985508, -2.2318082441740326, -3.1677451000476946, -2.7136257844502993, -2.4663578849843404, 2.1470642674041458, -0.80807516468244989, -6.8039934755458935, 11.462761097859508, -2.3200072133134322, 9.6476195054900042, -0.58323258923086063, 10.466810618341768, -7.3966372849294499, 4.629755222962074, 7.8499060069144697, -1.2532724460122664, -0.71383693570895246, -3.0356422051110532, -1.9171318821546592, -3.2905248775997, 1.1638954652799918, -0.8003122215164804, -7.6686022187651925, -13.379436000591463, -2.0784151782651223, 18.753674553162941, 37.596939514973172, 28.201887912811504, 26.441177038828386, 17.313699289624569, 8.2605837012813357, 1.9157665349814501, -1.0360449197814665, -2.460835808457071, -2.9435912594789997, -0.73745892985465067, -0.80698208098757895, 1.3512050433752707, -2.4806623635392717, -10.931358459008862, 20.399475046734903, -0.93914753929729033, 15.633586230850142, -11.808527677086202, 3.7955244266493091, 11.440683264496995, 0.22625003249094855, 2.3910143241324033, -1.5433230899620785, -1.7524971656611896, -3.4392335008624597, -1.2102890230410412, -1.1333111335531953, -0.86689054652400055, -0.77102037498198228, 5.869418392294464, -6.6911708498309208, 60.884710948918006, 43.621563233816303, 43.489599527484941, 33.141784821271486, 18.338861076357606, 7.7096141949396086, 3.1020965354201242, -0.33651454595896624, -1.858665713174497, -1.559439785993801, -1.7180049365163976, -0.49470837698954684, -5.7499783777067037, 1.2816063498557866, -1.6438131382345036, -14.80449956699707, 30.383248800411909, -20.687082737158004, 6.4291237134149926, 13.719178799652914, -0.15202331710012035, 5.6708331746907428, 1.1091130063870853, 0.43083164420390874, -2.423604198850752, -1.2614169037455483, -3.1633275417854718, -2.5428479485652051, -1.0621152637056996, 7.4921468465754772, 13.826120615216594, 48.56895003285198, 43.655937980459306, 56.143788442077074, 39.839822974553684, 28.907966432346075, 17.226166967005444, 9.493352687280801, 2.1248307248495761, -0.10304029106189502, -0.73643011249607004, -0.68086203692422775, -1.5726895782452095, -2.38329292549688, -6.3654622280486342, -12.648471231222512, 2.6502987300644998, 9.1941379366153466, -33.28553680715055, 17.791788634489134, 16.331392880531652, -1.0332061538290742, 9.6376857080971341, 4.1689669632335153, 3.4112249849533103, -0.63064111926197686, 0.028356847121888181, -1.1968584419289052, -3.3957285823190806, -4.7459408986565164, -4.790016252397657, -5.1593354130915703, 13.957717184896774, 20.367299127198688, 0.0, -4.5388798707327584, 12.10851784859719, 15.801052170121357, 15.786196904526062, 6.4645339091427658, 2.6357128260629894, 0.48128762074186304, -0.08159176932681525, -1.1399410088753745, -0.93916966756638098, -5.4695359499377414, -8.4733779831544798, -6.8332707948903639, -5.8759130146162706, -8.1552266379296245, 0.0, 8.1552266379295908, 5.8759130146162661, 6.8332707948903488, 8.4733779831544869, 5.4695359499377316, 0.93916966756637688, 1.1399410088753721, 0.08159176932681525, -0.48128762074185355, -2.635712826062981, -6.4645339091427667, -15.786196904526072, -15.801052170121331, -12.108517848597181, 4.5388798707327416, -43.655937980456407, -20.367299127197448, -13.957717184895019, 5.1593354130924185, 4.7900162523979217, 4.7459408986567313, 3.395728582319133, 1.1968584419289152, -0.028356847121888983, 0.63064111926197186, -3.4112249849533058, -4.1689669632334958, -9.6376857080971128, 1.0332061538290882, -16.331392880531649, -17.79178863448908, 33.285536807150621, -9.1941379366153146, -2.6502987300644807, 12.648471231222583, 6.3654622280486368, 2.3832929254968982, 1.5726895782452166, 0.68086203692425595, 0.73643011249607926, 0.10304029106193066, -2.1248307248495362, -9.4933526872806269, -17.226166967005291, -28.907966432345297, -39.839822974551964, -56.143788442073777, -60.884710948915583, -48.568950032850019, -13.82612061521586, -7.4921468465751166, 1.0621152637059226, 2.5428479485653219, 3.1633275417855051, 1.261416903745556, 2.4236041988507613, -0.43083164420392467, -1.1091130063870698, -5.6708331746907525, 0.15202331710012659, -13.719178799652925, -6.4291237134149712, 20.687082737158036, -30.383248800411859, 14.804499566997148, 1.6438131382345447, -1.2816063498557717, 5.7499783777067242, 0.49470837698955733, 1.7180049365164254, 1.5594397859938318, 1.8586657131745035, 0.33651454595900265, -3.1020965354201109, -7.7096141949395376, -18.33886107635751, -33.141784821271301, -43.48959952748401, -43.621563233814584, -18.753674553161851, 6.6911708498312086, -5.8694183922938254, 0.77102037498232001, 0.86689054652414343, 1.1333111335533379, 1.2102890230411023, 3.4392335008624833, 1.7524971656611932, 1.5433230899620796, -2.3910143241323984, -0.22625003249093578, -11.440683264497009, -3.7955244266492731, 11.808527677086207, -15.633586230850069, 0.93914753929731287, -20.399475046734821, 10.931358459008916, 2.4806623635393117, -1.3512050433752816, 0.80698208098763236, 0.73745892985467021, 2.943591259479049, 2.4608358084570767, 1.0360449197814969, -1.9157665349814337, -8.2605837012812007, -17.313699289624388, -26.441177038827938, -28.201887912810978, -37.59693951497237, 16.220892231174151, 2.0784151782650429, 13.379436000591378, 7.6686022187652325, 0.80031222151651993, -1.1638954652799158, 3.290524877599744, 1.917131882154679, 3.0356422051110541, 0.71383693570894702, 1.2532724460122711, -7.849906006914436, -4.6297552229620669, 7.3966372849294793, -10.466810618341761, 0.58323258923092358, -9.6476195054900042, 2.3200072133134997, -11.462761097859458, 6.8039934755459344, 0.80807516468246154, -2.1470642674041431, 2.4663578849843657, 2.7136257844503273, 3.1677451000476933, 2.2318082441740374, -1.9742581617985304, -8.4660409590374517, -15.161631552742712, -15.579335900760643, -14.803696771258227, 7.5614574111833424, 26.774160253651477, 34.683237345909006, 13.304914803324097, 5.5208279800184936, -1.6766227476275497, 2.3299130423218712, 1.078472870146282, 2.3401345171732237, 1.6227282885568277, 2.1726948090167162, -3.3870661540364337, -2.0288771570791524, 2.5011370399824764, -7.6358217658181076, 1.4383189901624107, -3.7090325345277098, 0.64649466936933497, -5.8520209082314834, 1.984545592632657, -6.6808538974877711, 1.9675974847606519, 1.2624013759488664, 1.5589303641798025, 3.2885833247038483, 3.4369058394752718, 2.1494090398895325, -2.3156698002103084, -9.0142686843579654, -13.06571402491778, -14.662722831106908, 1.3280628610827629, 14.352442298363767, -3.4901974522911301, -9.158035887658194, 0.30135116195024236, -4.9612937854963866, 0.062995146217542555, -2.1030121366447649, 1.6973044020589363, 1.457288603027268, 2.5856117766371063, -0.97601539997780706, -0.57130018359554458, 1.396534340840117, -5.0023110313506107, 0.5759598133955085, -2.357868531161305, 1.2120159426024788, -3.0263722138472202, 0.28712766434913656, -1.5705428299533217, 0.48575889336609929, -4.1370294214229943, 1.783719193855416, 3.2351264422424548, 3.5039625138893205, 2.956014710135662, 1.3188307478897938, -2.7986942524719765, -8.5655340270494893, -13.315592278460258, -5.1972517053507898, 0.67804977687200962, 2.7278278900777337, -56.64672530448825, -26.649306004568476, -41.383231407558036, -11.285893196605604, -14.42578813196792, -2.4367446291981882, 0.90733505109427226, 3.2049709482452906, 0.40634793694041338, 0.49606931853200609, 1.3537504507858567, -2.9970442493787934, -0.083484537520044699, -0.94029983857969957, 0.80232514327114979, -1.2514468795435723, 0.95269533500567261, -0.85050436376704475, 0.32098662744156231, -1.5779290838292681, -0.51108563039488553, -0.3475312947763517, 3.2496097582258234, 3.9932995411602126, 2.687214491416765, 0.13203799201742056, -3.7982495059531578, -7.8353712095884145, -3.5131066387109939, -9.5916945360520085, -25.283726226154414, -72.597133152543805, -219.86853407025291, -219.35042572397646, -77.645971701549428, -53.864359931260992, -21.52860711072962, -10.52910533038772, 3.2906850406638011, -0.86862341139231525, 2.1096176983280075, 2.4401275335966615, -1.1293193043526162, -0.5175199959040564, -1.6488725749132833, -0.24616862215974425, -1.2370080984335483, 1.8968331338372646, -0.68343802732828218, -0.11126643861832354, -0.99557444173586918, -0.97164442609393387, -0.53847505259674655, 1.408171623377654, 1.757965457581725, 3.6588840806054854, 3.730023177462626, 2.1373375769283727, -3.0400482541159497, 1.0472919151882829, 0.32173120932219762, -32.509510578136599, -106.65962198554575, -133.199718985283, -319.98109193688754, -62.432091745036324, -132.60719521227176, -54.20464622802902, -43.376064565046129, -14.089884469834875, -13.705577255286062, 1.5019641274622118, 6.2831724147445849, -0.36542819178094682, 0.31341792844047633, -1.3226298476418761, -0.88964240698418728, -1.5408441736495315, 0.71069133238493809, 0.021861987241514846, 0.026686439684265401, -1.1035985643780957, -1.1861458345364231, -1.082009036922065, -0.29997994657423205, 1.7050468427375622, 3.031311871551031, 3.818478372445528, 5.9901895103595697, 10.068399035154069, 15.790153056071421, -3.1394961579745826, -5.539540701046052, -89.166061674211662, -137.34002448547281, -339.65571159114575, 426.42488803140685, 53.557614725234323, 151.38668525810476, 18.829394319190286, 0.21915489532353391, -63.119371696969402, -10.682010714916952, 7.1043203483248183, -1.4514343357851665, 2.0000864020848317, -1.0967346902818478, -0.33652577635987613, -2.5606009702294568, -0.075739057981359575, -0.6297575158660037, -0.1002153489582, 0.34443735246552043, -0.35567626693504045, -2.4151219203224019, -0.95474631060256077, -0.58454842364747928, 2.4591204313835524, 4.6279710639736518, 7.2190125620904322, 10.8478141614737, 18.113211298120788, 35.095892631083665, 23.228148927519459, -81.527654709232067, -75.498646737828622, -221.27599419248904, 203.12933548649244, 1583.8824730375886, 1881.6935115818235, 969.72646934413342, 542.55543458513785, -19.428895667388197, -23.350930435925818, -2.7323622379307726, -20.048308524208668, 2.4299863370655803, 0.67396364813076437, -0.28270713586437829, -3.3794759874707645, -2.0455422971725761, -2.6548043613806773, 1.248210330242034, -1.3645618319173671, 0.78798091056567598, -1.130442042783935, -3.3535274581592822, -3.4782462601711517, 1.1480426683318268, 3.9818184076005192, 7.6639325468155093, 13.090865763493181, 20.166794237424288, 40.697925231356948, 58.678662168810511, 13.453603503282435, -69.371971565542452, -317.34459246185418, 242.65732454028887, 833.14970460709446, 1814.4894639161919, 1707.4355014735293, 1928.5847150696216, 464.93982303099165, 317.96189631496537, 77.039516694542414, -94.400929285346379, -12.369503115173439, -4.5934633035837242, 0.38778650675346871, -3.6019889835361467, -3.5521022242715157, -4.1535160740616579, -1.9878096211735969, -2.4480270844921423, 1.0231985430638744, -2.1386941287610495, -1.7756565752032289, -0.071309220847206922, -0.28567330333819857, 0.19099280196170093, 8.2634536416369908, 9.5628328869678718, 18.720186517652159, 37.286147332183262, 72.736940205390795, 71.661879904778388, -97.912476061205524, -380.4257804868011, -268.24846212420823, -349.63255334068992, 1058.8873901892528, -2087.2412865177334, 676.63432998511416, -1074.2990124084556, 901.86972095078841, 786.94393395792224, -119.57179721034137, -50.577306165813212, -53.34828178316841, -11.529259915073451, -2.9920185974423683, -2.8189903732028481, -5.0447569422201788, -2.5445971528130604, -10.107544184404102, -4.9832600695276215, 2.2636764639027978, 5.3737401851528253, 0.20461611038107308, -3.6890839097334074, -1.119402655197085, 0.20770764759927202, 6.7708773337385892, 12.544755138678777, 25.254222647826818, 57.14272356927944, 74.489552074506946, 82.566517160157588, -71.610686309506534, -537.49064424335074, -271.57679130866813, 540.57641351250493, -2934.8250858413712, -11955.369188695569, -14506.214536688902, -637.22651951043315, 204.04622542072232, 142.43138300706894, 154.61175360170483, -106.41042136702355, -66.543583791323357, -35.880437410281012, 2.1360966686551617, -13.062832664150186, 0.86756998254746465, -12.950261495215663, 0.6068089103450538, -17.40763581076148, 13.386771464725289, -0.87386354462209659, 1.6150967690552112, -14.546101055958459, -2.2874141051365076, -1.0182743702449171, 12.902690844348964, 20.382558491276768, 30.589238972650698, 63.125416501511289, 68.536452109292355, 208.25011765184536, -277.4413124787215, 295.01076934839597, -202.32180136863923, -2363.0865643328207, 346.60702359456963, -16946.371225111554, -1299.266384884367, -4625.1554139998561, -3868.3607546707453, 832.38513354395968, 28.40638042173725, -119.64570653889508, -110.70270321520742, -29.597688207921749, -10.925693066643381, -9.9160507928571686, -9.5491377266944308, -18.949463773817815, -10.847407040619126, 14.016932746732627, 8.8305069728758134, 12.446635345723296, -18.982014759378298, -13.73839338407706, -10.939736173955323, 20.777136166915735, 12.48897777935125, 25.365956236318532, 8.1960517939270616, 65.249140172269151, 42.285393219366988, 206.10314812122587, -109.43221584453876, -109.26740577843286, 1171.8258105148195, 1841.1630655848473, -11272.400870696491]}};

var right_eye_filter = {"real": [1.8229079259010603, 0.097810498648582461, -0.55840092137248587, -0.99970462842356578, 0.18757876559092043, -0.81617694863114465, -0.17457078018685562, 0.13822182613073089, -0.13830885288846723, 0.14088352498892928, -0.11242245121823281, 0.64190144530750459, -0.18081689212517704, 0.84757445423403044, -2.0367692338977426, 2.3450651085370726, -3.2227771639624789, 2.3450651085369945, -2.0367692338977101, 0.84757445423403854, -0.18081689212518082, 0.64190144530750537, -0.11242245121823641, 0.14088352498894169, -0.13830885288846723, 0.13822182613073003, -0.17457078018684977, -0.81617694863113666, 0.18757876559091927, -0.99970462842356267, -0.55840092137248509, 0.097810498648581323, -1.0960340088024074, -1.9119199452512301, -2.4095702665388479, -0.55162352804953241, -0.50505328014021589, -1.3068505835636635, -0.49437187443975589, 0.23127204626685588, 0.043350794355039197, -0.071787406154436809, 0.070854449837205943, 0.17221102955705567, 0.44602122877614064, 0.074148205591900399, -0.015611524439393865, -0.79896724757292914, -0.67155961060437053, -0.49375784053273253, 0.028167647498424205, -0.3191891366103583, 0.53582170174408217, 0.20183473326047152, 0.12710467168208503, -0.071400263514192097, 0.074560604448017792, -0.12451483496807111, 0.051358135035693774, -0.77282398113384443, -0.42160647485916963, -1.1311602269894514, -1.4625324828122208, 0.3337871926813808, -1.77011609742092, -1.4302518240108344, -3.0428213499062116, -0.60631225487833551, -0.13683987631719793, -2.3692424085383519, -0.34105804260324346, -0.25969657852057992, 0.23610433796619859, -0.27685168579683173, 0.16986402710259335, 0.073612585805396993, 0.57704911017089122, -0.059867579479423505, 0.51680219528175042, -1.0666730319382505, -0.0056552094670686828, -0.75809939060745812, -0.064701320250044561, -0.19832318162238824, 0.37732986906226257, 0.20777584924550152, 0.15451076629740465, -0.013991851319401968, 0.13841341497021592, -0.086883949036928532, 0.25261423619714551, -1.1991385138433881, -1.9867837341406693, -2.275962773839737, -3.1079150257531585, -2.3995834975192567, -2.2937692343725722, -2.1872173045799017, -1.3988308931230504, -2.2073957813815439, -1.9667941627002865, -1.6748432066635108, -0.92650772032564344, -0.054197237250052703, -0.18026213329416008, -0.15232124758200827, -0.016555865669302033, 0.31289505893870784, 0.42806542968274719, 0.23124350286854853, -0.087366382764607439, -0.12552132795947232, -0.12547908619418155, -0.38089709050886683, -0.28444076658768042, -0.11944152964360766, 0.12056136689031131, 0.21189843131890268, 0.16581598470154221, 0.13376820638890077, 0.22205496825803592, 0.15120100371874812, 0.12380717441381156, -1.2352596784322261, -2.2222759116849211, -3.3782856809989692, -3.415076124219421, -2.7282977589825297, -1.5576789730825455, -1.6478932560447972, -2.067907968457094, -2.004290735457912, -2.6427456963680775, -1.9309474679838621, -0.16864865516592528, -0.22642518078963322, -0.13321306237600422, -0.20692971953006373, -0.068344938706069003, 0.14493915164047724, 0.46505832474821007, 0.32789481378377733, 0.1852927254910838, -0.1356638873386786, -0.03585192275390197, -0.17203834043005656, -0.15954900238018477, -0.059581713131431221, -0.017228379259871054, 0.12402459410135006, 0.1524538089963059, 0.16271946063708961, 0.24817234249560427, 0.26053579199245763, -0.04769125558408293, -0.81459993031792033, -2.0740660744196808, -3.1544891598954261, -2.9212825837849565, -2.4077419089431702, -1.6411629435379833, -1.576697627349694, -1.8007661961299257, -2.0001525985045707, -1.5107274159689619, -1.2839600483166405, -0.63196711118813498, 0.31071574064221719, -0.21370915949408273, -0.16393119308951343, -0.097481704015839707, 0.13767875532837173, 0.23695911940697742, 0.36158774196732829, -0.00066962538216623158, -0.1026502520704758, -0.18176695609293456, -0.059120611570986824, -0.13702193484548986, -0.0073538013157976303, 0.020281061259494877, 0.15032501459982059, 0.16876214929227357, 0.18780650312631536, 0.16308782299108887, 0.18583847563977354, 0.44907296092364651, -0.020800703185301726, -1.4297091333258476, -2.4259185480628145, -2.4053777736536626, -1.9821307035675879, -1.3372642957080885, -1.2989780551182537, -1.3731553217305845, -1.2269479805895975, -1.2032897772528626, -0.22206046284030648, -0.18415548289109809, -0.23415549346578352, 0.077936325689887453, -0.21633114300736572, -0.059449722211757941, 0.013318329293146938, 0.28775437300883983, 0.21963323601338633, 0.16649953505637383, -0.18974789333555656, -0.17379580132236441, -0.23416130311702454, -0.13478379984035538, -0.13112027137758978, 0.033522350558872455, 0.13259733500476653, 0.171380750952115, 0.15737574769334917, 0.092477162456605533, 0.13448891283084696, 0.41161459424817443, 0.35609170613712005, -0.55316449850203808, -1.540851379466907, -1.8005308506197801, -1.5719273509294158, -0.94040577953316051, -0.83899416480020927, -0.77336416928591256, -0.66562411016509948, -0.4002640865408188, -0.36750971454913561, 0.2414435659750552, -0.065677188177910795, -0.21341675767851534, -0.041369287190265204, -0.087404842921799578, 0.04460147122538273, 0.19583406251777175, 0.1830038441344799, 0.012219595757112435, -0.034593273427629516, -0.14926486203867156, -0.036665364495836496, -0.1306231554860259, -0.086381624051938985, 0.0020423415177502173, 0.11239290885891448, 0.12781279786870536, 0.10954302373131516, 0.14996801038078228, 0.12965098836907352, 0.14324588543961356, 0.0567484033157079, -0.079342592710004017, -0.77145898274835445, -1.1705860882999304, -1.0754357146023004, -0.5463269466436349, -0.37732231762995055, -0.29185933229148303, -0.18245346577017127, -0.15252816280103537, 0.00992307413149221, -0.05334501261419268, 0.10643359538216246, -0.069635302172998029, -0.20853401810086858, 0.03871878333386132, 0.024793140061197767, 0.22215018549405419, 0.19090693911120471, 0.13458393481660855, -0.1226120784542075, -0.060406846903667433, -0.083060383383984554, -0.016925686434954099, -0.033352938195947639, 0.056474422402956102, 0.04873056181166667, 0.092098221631537164, 0.13542751998114597, 0.15540921558269635, 0.10094403852684107, 0.050489672984821425, 0.089840998861512991, 0.018445471788018183, -0.44280324831238471, -0.64559475523542686, -0.64123409592876457, -0.25571779361731761, -0.15489773108593363, -0.10375533557574061, -0.074018947274090252, -0.074021765277388324, -0.034956034554291714, 0.065331700319185393, -0.079299771847393716, 0.035725424097319471, -0.10252861876832269, -0.22107962950770202, 0.20755848067860572, 0.23157078159613848, 0.34013467279608883, 0.23914223921471078, 0.043591964044129211, -0.20186945117744784, 0.04172561648658725, -0.038784831650455286, 0.076248151019351923, -0.060276976551200255, -0.0031358189183758184, 0.13054366560393055, 0.16122088670834137, 0.12129220423951546, 0.064624059271846127, 0.063153633523538075, 0.036868954773612302, 0.025599835220982206, -0.17872521900892777, -0.37577894082797708, -0.33418602947194676, -0.19087946198724595, -0.092430559387195924, -0.15843188266072916, -0.1178716916883143, -0.079832857077627425, -0.038455609721169501, 0.015668866487966955, 0.0052656766019976578, -0.12562510603385343, -0.024970842739258457, -0.030574039563531708, -0.21969881405226874, 0.59191387060782108, 0.39965780946904006, 0.22624763974293111, 0.077143781693831126, 0.078542516808501189, -0.0038697706067326258, 0.21858339549378, 0.083870133819846815, -0.19149929243694441, 0.053673315180545467, 0.22201073984886688, 0.16922848807595664, 0.12511007660405662, 0.095111634855384236, 0.084195395624058103, 0.010195208908487656, -0.016589123313681509, -0.076922699191606214, -0.17023093289451025, -0.26998035777525775, -0.28189859430096476, -0.32749231042729776, -0.240474927671948, -0.27714150491609907, -0.16857870511456344, -0.056310601502134704, 0.026291702225850644, -0.028508101274427284, -0.066307061790986485, -0.27187198211633473, 0.066219257946262391, 0.32825613521942065, -0.21303758386300781, 1.1917438602108388, 0.2807025357241521, 0.23447144426975483, 0.16511529097099151, 0.40113579910945835, 0.017152293761210982, -0.1450564820739392, 0.13833755803781553, 0.2361453534315919, 0.16033553043559276, 0.17442784770871, 0.17988585089228476, 0.143828955260344, 0.07132752932041192, 0.016244981011210025, -0.013577696292888701, -0.024045478766089354, -0.11184718959298219, -0.14666574805211069, 0.0080876817899727359, -0.16842937423863127, -0.43543489379944889, -0.29975286851709859, -0.22753898297998879, -0.069036032916579657, 0.014977564920406507, -0.020420877906038481, -0.13524198281736244, -0.15739380447541659, -0.33401260934309918, 0.46815371603113554, 0.93021729768322625, -0.26733338229185055, 1.4166629302294091, 0.60664479608939625, 0.3603810860087503, 0.4605515924213614, -0.32071019427276837, 0.15916183580862947, 0.55142162929110317, 0.20419289205444074, 0.051921880744650394, 0.24679565757292954, 0.2706202610478356, 0.15590400545853877, 0.081859441051199647, -0.040980231107277158, 0.019596242113331462, 0.094364427919773758, 0.065278519802775978, -0.029838116901178446, 0.15435464918439099, -0.039386655944161249, -0.077796406016204619, -0.28158476919608044, -0.033110128022724852, 0.018548614200940845, 0.067046107651037731, -0.044245201781535203, -0.1525534131110006, -0.22929077462390082, -0.085302710712388263, -0.049302393683332231, 0.85744368463127818, 1.8448220136861364, -0.68744556831975012, 1.810880881855317, 0.66237039703709566, -0.31567909429421087, 0.83232217907475781, 0.97863365922997991, 0.30994126117364096, -0.041163173314528664, 0.17633262873947353, 0.23161854928102976, 0.25087704366702551, 0.19868696729956423, 0.093569887694449877, 0.032362324442922322, -0.0027507666136337679, 0.21821986304703905, 0.37880734420129814, 0.40096920495374666, 0.93561154193066887, 0.61326039955329958, 0.1638528799494483, 0.1009541629695858, 0.056601591846233243, 0.15561323638878471, 0.073000292812971587, -0.0055442453146915694, -0.13817032242271102, -0.11770539535407684, -0.096245457278447424, 0.46442894125096079, 0.45587938373772585, 1.3432089992982477, 2.5059762271849535, -0.93829697639874887, 0.6492144586928863, 0.79782532034480069, 1.3025551528258366, 0.94174974100732234, 0.099921469659734793, 0.28712786266920659, 0.15853680079243379, 0.130717335778855, 0.17042917695126153, 0.19986443324539621, 0.13725228360851915, 0.029723534135295476, 0.098143635338702864, 0.19047751491127032, 0.49335337973435212, 0.78089304412099847, 0.29686587428828037, 0.38131881367791126, 0.40643673652514922, 0.19053499830856688, 0.20033550195294633, 0.13287714582076576, 0.10533560171090586, 0.00085487897815388243, -0.049435296406267468, -0.081309411554979491, 0.057947716710403128, 0.48207744089263116, 1.2538880645112733, 1.2204758332556214, 0.87793863876625744, 1.8796293205435075, -0.45285576127696725, 2.2308682152547137, 1.2382711063757557, -0.035207936995945072, 0.88175672278306683, 0.55916031961657942, 0.18718127547628777, 0.011764673252274449, 0.068747006595272581, 0.13522321325294812, 0.17151312046437606, 0.10315365060575944, -0.0011239413451798732, 0.14993420590829387, 0.38959110228082239, 0.38845862264410469, 0.69716730182732189, 0.51802064393501412, 0.28837767001488912, 0.10583936790249669, 0.071016615670210698, 0.16001776032454934, 0.10715662619027419, 0.056188716759532117, -0.035825840155044517, 0.0061469917968111389, 0.089493877843564718, 0.68387295014825999, 1.1016519558168656, 1.3182383587355595, 0.66635382040153524, 0.54863673632097221, 4.0547770231048936, 0.54863673632097154, 0.66635382040153424, 1.3182383587355562, 1.101651955816866, 0.68387295014826144, 0.089493877843565064, 0.0061469917968106288, -0.035825840155044517, 0.056188716759530431, 0.10715662619027334, 0.16001776032454856, 0.071016615670210018, 0.10583936790249539, 0.28837767001488829, 0.51802064393501346, 0.29686587428830957, 0.38845862264412201, 0.38959110228084654, 0.14993420590831108, -0.001123941345177488, 0.10315365060576374, 0.17151312046437922, 0.13522321325295103, 0.068747006595273552, 0.011764673252275867, 0.18718127547628682, 0.55916031961657786, 0.88175672278306205, -0.035207936995944634, 1.2382711063757532, 2.2308682152547106, -0.45285576127696175, 1.8796293205435037, 0.87793863876625788, 1.2204758332556218, 1.2538880645112747, 0.48207744089263271, 0.057947716710405182, -0.081309411554979547, -0.04943529640626837, 0.0008548789781495917, 0.10533560171090255, 0.13287714582076099, 0.20033550195294705, 0.19053499830855697, 0.40643673652513423, 0.38131881367792647, 0.93561154193068807, 0.78089304412102711, 0.49335337973436416, 0.19047751491127579, 0.098143635338704696, 0.029723534135297252, 0.13725228360851988, 0.19986443324539666, 0.17042917695126267, 0.13071733577885625, 0.15853680079243432, 0.28712786266920415, 0.099921469659732226, 0.94174974100731712, 1.3025551528258328, 0.79782532034480169, 0.64921445869288619, -0.93829697639874454, 2.5059762271849548, 1.3432089992982421, 0.45587938373772335, 0.46442894125096235, -0.096245457278447175, -0.11770539535407766, -0.13817032242271216, -0.0055442453146945193, 0.073000292812969436, 0.15561323638878341, 0.056601591846233208, 0.10095416296958591, 0.16385287994945555, 0.61326039955330602, 0.15435464918440814, 0.40096920495376004, 0.37880734420131523, 0.21821986304704899, -0.0027507666136306606, 0.032362324442924542, 0.093569887694449336, 0.19868696729956342, 0.25087704366702551, 0.23161854928103059, 0.17633262873947425, -0.041163173314526666, 0.30994126117363635, 0.97863365922997669, 0.83232217907475781, -0.31567909429421592, 0.66237039703709699, 1.810880881855315, -0.68744556831974835, 1.8448220136861369, 0.85744368463127729, -0.049302393683329601, -0.085302710712388291, -0.22929077462390179, -0.15255341311100187, -0.044245201781539331, 0.067046107651035594, 0.018548614200940294, -0.033110128022720509, -0.28158476919607578, -0.077796406016200054, -0.039386655944146774, 0.0080876817899761342, -0.029838116901173162, 0.065278519802780599, 0.09436442791977824, 0.019596242113334855, -0.040980231107276652, 0.08185944105119855, 0.15590400545853553, 0.2706202610478356, 0.24679565757292907, 0.051921880744651677, 0.20419289205443947, 0.55142162929109917, 0.15916183580862936, -0.32071019427277109, 0.46055159242136157, 0.36038108600874996, 0.60664479608939803, 1.4166629302294049, -0.2673333822918505, 0.93021729768322625, 0.46815371603113448, -0.3340126093430979, -0.15739380447541709, -0.13524198281736302, -0.020420877906040739, 0.014977564920405115, -0.069036032916578546, -0.22753898297998629, -0.29975286851709465, -0.43543489379944272, -0.16842937423862725, -0.28189859430095965, -0.14666574805210519, -0.1118471895929717, -0.024045478766081197, -0.013577696292883971, 0.016244981011212731, 0.071327529320409616, 0.14382895526034015, 0.1798858508922834, 0.17442784770871067, 0.16033553043559218, 0.23614535343159299, 0.13833755803781575, -0.14505648207394348, 0.017152293761212339, 0.40113579910946001, 0.16511529097099317, 0.23447144426975502, 0.28070253572415382, 1.1917438602108381, -0.21303758386300692, 0.32825613521942509, 0.066219257946262128, -0.27187198211633529, -0.066307061790987387, -0.028508101274429792, 0.026291702225849072, -0.056310601502134003, -0.1685787051145608, -0.27714150491609657, -0.24047492767194564, -0.32749231042729471, -0.19087946198724157, -0.26998035777525148, -0.17023093289450597, -0.076922699191602426, -0.016589123313677224, 0.01019520890848913, 0.084195395624057118, 0.095111634855381461, 0.12511007660405613, 0.1692284880759545, 0.22201073984886613, 0.053673315180547035, -0.19149929243694658, 0.083870133819847037, 0.21858339549377928, -0.0038697706067313737, 0.078542516808500606, 0.077143781693833222, 0.22624763974292969, 0.39965780946904111, 0.59191387060782175, -0.21969881405226852, -0.030574039563530501, -0.024970842739259945, -0.12562510603385371, 0.0052656766019967011, 0.015668866487966397, -0.038455609721168606, -0.079832857077626357, -0.11787169168831266, -0.15843188266072725, -0.092430559387194119, -0.25571779361731145, -0.33418602947193976, -0.37577894082796842, -0.17872521900892277, 0.025599835220983438, 0.036868954773614897, 0.063153633523536298, 0.064624059271844544, 0.12129220423951372, 0.1612208867083402, 0.13054366560393141, -0.0031358189183775171, -0.060276976551199596, 0.076248151019351354, -0.038784831650455917, 0.041725616486585258, -0.20186945117744698, 0.043591964044126186, 0.23914223921470992, 0.34013467279608772, 0.23157078159613964, 0.20755848067860685, -0.22107962950770196, -0.1025286187683216, 0.035725424097318138, -0.079299771847393549, 0.065331700319185088, -0.034956034554288036, -0.074021765277385243, -0.074018947274086908, -0.10375533557573537, -0.15489773108592672, -0.5463269466436349, -0.64123409592876401, -0.64559475523542686, -0.44280324831238616, 0.018445471788017767, 0.089840998861510313, 0.050489672984820259, 0.10094403852683856, 0.15540921558269635, 0.13542751998114574, 0.092098221631536192, 0.04873056181166658, 0.056474422402955825, -0.033352938195949811, -0.016925686434957218, -0.083060383383990535, -0.060406846903667433, -0.1226120784542102, 0.13458393481660744, 0.19090693911120402, 0.22215018549405383, 0.024793140061198739, 0.038718783333861403, -0.20853401810086775, -0.069635302172998029, 0.1064335953821621, -0.053345012614191917, 0.0099230741314936706, -0.15252816280103507, -0.18245346577017069, -0.29185933229148298, -0.3773223176299495, -0.94040577953315341, -1.0754357146022915, -1.1705860882999222, -0.77145898274834734, -0.079342592710004489, 0.05674840331570382, 0.14324588543960787, 0.12965098836906794, 0.14996801038078103, 0.10954302373131049, 0.12781279786870486, 0.11239290885891398, 0.0020423415177517252, -0.086381624051941733, -0.13062315548602732, -0.036665364495841957, -0.14926486203867201, -0.034593273427634852, 0.012219595757109143, 0.18300384413447635, 0.1958340625177715, 0.044601471225384402, -0.08740484292179726, -0.041369287190262873, -0.21341675767851565, -0.065677188177906895, 0.24144356597505726, -0.36750971454912684, -0.40026408654081164, -0.66562411016509171, -0.77336416928590557, -0.83899416480020206, -1.3372642957080858, -1.5719273509294132, -1.8005308506197781, -1.5408513794669076, -0.5531644985020383, 0.35609170613711127, 0.4116145942481696, 0.13448891283084041, 0.092477162456604783, 0.15737574769334561, 0.17138075095211283, 0.13259733500476581, 0.033522350558872378, -0.13112027137759125, -0.13478379984035713, -0.23416130311702804, -0.17379580132236419, -0.18974789333555622, 0.16649953505637016, 0.21963323601338522, 0.28775437300883855, 0.013318329293152468, -0.059449722211754409, -0.21633114300735917, 0.077936325689886982, -0.23415549346577649, -0.18415548289109335, -0.2220604628402971, -1.203289777252857, -1.2269479805895944, -1.3731553217305832, -1.2989780551182524, -1.64116294353798, -1.9821307035675826, -2.4053777736536572, -2.4259185480628087, -1.4297091333258454, -0.020800703185307738, 0.44907296092364196, 0.1858384756397643, 0.16308782299108818, 0.18780650312631128, 0.16876214929227229, 0.15032501459981898, 0.020281061259495404, -0.0073538013157964212, -0.13702193484548603, -0.059120611570979247, -0.18176695609293375, -0.10265025207047287, -0.00066962538216844151, 0.36158774196732713, 0.23695911940697681, 0.13767875532837892, -0.097481704015832935, -0.1639311930895031, -0.2137091594940832, 0.3107157406422304, -0.63196711118812665, -1.2839600483166249, -1.5107274159689543, -2.0001525985045667, -1.8007661961299211, -1.5766976273496909, -1.557678973082544, -2.4077419089431698, -2.9212825837849552, -3.1544891598954266, -2.0740660744196822, -0.81459993031793032, -0.047691255584086843, 0.26053579199245114, 0.24817234249560396, 0.16271946063708537, 0.1524538089963024, 0.12402459410134949, -0.017228379259870755, -0.059581713131428661, -0.15954900238017725, -0.17203834043004781, -0.03585192275390215, -0.13566388733867182, 0.18529272549107739, 0.32789481378377766, 0.46505832474821029, 0.14493915164048596, -0.068344938706060163, -0.20692971953005135, -0.13321306237600447, -0.22642518078961474, -0.16864865516591712, -1.9309474679838514, -2.6427456963680735, -2.0042907354579094, -2.0679079684570931, -1.6478932560447959, -2.2937692343725686, -2.7282977589825266, -3.4150761242194196, -3.3782856809989692, -2.2222759116849189, -1.2352596784322352, 0.12380717441380877, 0.15120100371873849, 0.22205496825803597, 0.13376820638889372, 0.16581598470153905, 0.21189843131890046, 0.1205613668903118, -0.11944152964360859, -0.2844407665876833, -0.38089709050886922, -0.1254790861941803, -0.12552132795948226, -0.087366382764606329, 0.2312435028685467, 0.42806542968274885, 0.31289505893871633, -0.016555865669298363, -0.15232124758199445, -0.18026213329416074, -0.054197237250036744, -0.92650772032563145, -1.6748432066634964, -1.9667941627002861, -2.2073957813815404, -1.3988308931230464, -2.1872173045798968, -1.7701160974209202, -2.3995834975192563, -3.1079150257531603, -2.2759627738397428, -1.9867837341406698, -1.1991385138433992, 0.25261423619714135, -0.086883949036943853, 0.13841341497021575, -0.013991851319408978, 0.15451076629740176, 0.20777584924549833, 0.37732986906226257, -0.1983231816223982, -0.064701320250050057, -0.75809939060751286, -0.0056552094670687921, -1.0666730319382445, 0.5168021952817432, -0.059867579479421881, 0.57704911017089144, 0.073612585805408248, 0.16986402710259826, -0.27685168579682035, 0.23610433796619823, -0.25969657852055217, -0.34105804260323408, -2.3692424085383355, -0.13683987631719494, -0.60631225487832385, -3.0428213499062093, -1.4302518240108333, -1.0960340088024054, 0.33378719268138207, -1.462532482812221, -1.1311602269894572, -0.42160647485916913, -0.77282398113385165, 0.051358135035689895, -0.12451483496808054, 0.074560604448017945, -0.071400263514198398, 0.12710467168208478, 0.20183473326046905, 0.53582170174408328, -0.31918913661037829, 0.028167647498415164, -0.49375784053277538, -0.67155961060436897, -0.79896724757297555, -0.015611524439376389, 0.074148205591898927, 0.44602122877613765, 0.17221102955705844, 0.070854449837212188, -0.07178740615442418, 0.043350794355038948, 0.23127204626685866, -0.49437187443974595, -1.3068505835636535, -0.50505328014021666, -0.55162352804951942, -2.4095702665388461, -1.9119199452512292], "bottom": {"real": [4046.4241509904555, 28839.973879370627, 10191.512031448321, 4381.7455807293263, 2071.9548413915022, 1107.2825867454831, 528.34551553150311, 308.97330299730402, 185.81061014820492, 130.08108126519062, 87.65365496956511, 59.09553675383809, 43.506269455848397, 34.575222311681692, 30.798097452905282, 25.523061346781674, 25.951646714603154, 25.523061346781589, 30.798097452905211, 34.575222311681706, 43.506269455848361, 59.095536753838019, 87.653654969565167, 130.08108126519082, 185.81061014820492, 308.97330299730379, 528.34551553150334, 1107.2825867454833, 2071.9548413915031, 4381.7455807293318, 10191.512031448323, 28839.973879370646, 13493.224333003216, 16392.519951456161, 4871.8353885075767, 1571.3337493887716, 567.64096645153529, 373.6703632325283, 244.24660034041102, 142.79597933554817, 106.13580676907658, 75.858347622710525, 58.061322564190078, 39.854538473021549, 30.569538051408923, 25.094853941398878, 22.929757321982656, 20.808198660555963, 22.04888812794956, 23.622937534136554, 29.138683151168021, 38.082981074479143, 54.39934218587819, 77.826621492011839, 121.61239162949131, 177.63416429867098, 254.30481783551838, 430.47940336338183, 664.63525809030591, 1219.8477670929726, 1935.5749815900854, 3358.0980991141623, 6069.2175916916049, 8660.7039291701622, 15384.540019848313, 11175.737303494303, 2216.879483433469, 674.59859086304846, 308.57630937597452, 222.13525013253278, 144.62539946693838, 89.351789047349641, 66.200283784031811, 53.5880480094322, 38.104835124660447, 27.112718944065261, 23.033311632510099, 19.858863545618888, 19.408029208714364, 19.38215160767772, 20.871711217553681, 24.726583224259915, 29.429525390428484, 39.688247321389746, 54.026677368666491, 71.509189516532359, 109.7023701375058, 145.30501493896827, 183.84765943176359, 237.38092500065588, 353.51632607639505, 531.27032729465725, 719.17065363709912, 1367.9024800373177, 2857.7322293344778, 8625.0149751320059, 3926.8611003464521, 2444.3527388912776, 1028.8359899085383, 585.56166554295623, 332.57137935864864, 180.13762334769046, 105.17011934648103, 61.496465220623485, 45.504823244138095, 37.908460172493172, 30.846794326525625, 21.222093075405745, 18.718288087401802, 17.679602313216872, 17.409839934721365, 17.794391783727544, 18.55926931620083, 21.348409581307148, 24.39997962145603, 30.767061622228063, 38.103008129659166, 44.81251118779231, 61.577721983396955, 87.291485884508717, 106.99539363556599, 134.18057082099898, 169.37440787587974, 269.99443977661724, 364.56248799144709, 674.50839546082307, 1521.022929703357, 2634.7343778883865, 2024.849011267275, 1926.8239401023498, 1491.9478410044117, 611.83627918833088, 316.43682665588744, 155.99683828976544, 89.673379864750018, 50.013293885425135, 37.528009707176061, 32.4518569290627, 24.966146703509182, 19.130394004010906, 15.343856356149828, 14.843145659919786, 14.861292517540363, 15.036681247518919, 15.312687438440403, 17.138709070816013, 19.446116741695533, 22.850207518010233, 27.457440287528367, 32.581505308308088, 40.985078317877736, 54.412823137373685, 75.762437644265901, 94.733059346587339, 128.1467358253355, 182.29010502256293, 275.95498073034861, 527.32695578494008, 1132.5867520786235, 1866.4077869841046, 2266.7457843108455, 2107.8495169895441, 1234.6592459623437, 673.70954283481137, 266.6693703142779, 139.47307538222793, 66.730731113314476, 41.917980917088734, 31.2590211502941, 26.786032891522527, 20.35627716260823, 16.038244917569745, 14.499148814332527, 13.459360575928033, 13.144389071570192, 12.567429683418663, 12.230328765246885, 13.466217665828195, 15.328301760046475, 17.144679989146308, 19.68756841312236, 24.074134044705339, 32.727034612286992, 43.026808431859955, 59.193773906230369, 78.820512563030945, 92.994178544706472, 120.42091176843755, 201.9182636225097, 436.44287859172925, 877.8862778005057, 1439.6385876846614, 1649.9911393469049, 1699.4038415062987, 1039.2560892520887, 495.51661529689329, 235.81544665938515, 113.58561083280532, 60.921551746622946, 35.725757089786256, 25.024628727892466, 22.490208420088727, 17.759823791920191, 13.543616943986416, 12.938917306264663, 11.956364952724385, 11.421939705679733, 11.372903592178258, 11.032765269291394, 11.408229662647633, 12.481114469594095, 13.60275540468532, 14.91952183557448, 17.807798659576154, 23.277252196336814, 31.950664619854233, 47.416694920263616, 56.268912844615144, 64.959626427586997, 85.070809539003847, 130.82046704446515, 257.67571935262418, 627.46721263291658, 1194.1693767483914, 1080.1900960853154, 902.674413704009, 600.04055794738224, 333.23520458041997, 166.5653831933532, 88.946305547571072, 50.717393232745131, 31.897215329939996, 23.060146944216118, 18.579780562544169, 15.572760511446054, 13.273777678692191, 12.224913053368416, 11.68638594982613, 10.392044527527784, 10.422646713468026, 10.113778554702394, 10.225072497047023, 10.776357928238541, 11.471833692521193, 12.585855785544721, 15.131216045563061, 18.856673404562432, 24.858220592568422, 32.657280802203381, 43.037987696240307, 51.536289299313914, 59.08796022590117, 89.418846636419246, 168.44781050782811, 382.14409608709747, 745.27797642375992, 544.93523795856049, 546.63255943550212, 362.69938393816085, 222.4530875664997, 124.78354764285341, 69.155510455436144, 41.705006957506896, 26.952419331651356, 19.778555411511121, 17.743475715717917, 14.471795156401063, 12.765457059373892, 11.899864411981882, 10.85818004844373, 10.325708895539952, 9.5692719745256376, 9.2113480820863387, 9.0680329970290128, 9.7734432560375648, 10.488284402310196, 11.323197140674811, 13.23504887223821, 16.674941080314653, 20.771526297093455, 27.569417818570489, 33.103068481604573, 37.131609415223707, 47.964760430736398, 66.578817507603574, 119.94841754808802, 237.18123750436465, 409.50889944188236, 307.01018353832677, 288.66760291638428, 236.16368190120858, 154.32114106589464, 92.303773245695481, 55.468009578938926, 35.843010431224045, 24.692475316947775, 18.31794202383119, 15.815283825259236, 13.970653175227845, 12.321358605090134, 11.166914874614415, 10.532076178141539, 9.9738582700772742, 9.4395532322934326, 8.8524779931517035, 8.5960325163902862, 8.7233150994257276, 9.4605465058702709, 10.702831587760128, 11.695932356823622, 14.406509159208545, 18.893042693096525, 25.031728162044864, 28.237283751764654, 32.110154440123182, 38.781782989560291, 56.412497176172337, 92.156945934839797, 157.10935603874077, 244.24925081430641, 218.20881258539546, 216.90433976724358, 169.13695858014356, 112.62100758585845, 68.70530790516959, 45.982361360241526, 32.725373555898344, 22.221572369125827, 17.135020014783734, 14.281960328360523, 12.699315839076576, 11.637118060780567, 11.004832210877884, 10.393382355075936, 9.6259602166533789, 8.7652013845999868, 8.5650120282362092, 8.5854476517741869, 8.1851518687531879, 9.1270872628576836, 9.8383977040140795, 11.693775785132257, 13.287947209975439, 16.544820171465705, 21.672488061489553, 27.762087103200397, 29.796933806391767, 35.443188490268959, 51.826480475596476, 76.692136618774384, 114.45554660359461, 172.81330831382385, 150.57608076172548, 151.07347077509746, 119.80060075739978, 84.307080896770614, 55.59442395215062, 38.954469932130465, 27.855857296119165, 21.147383956632677, 16.095362427853317, 14.114901096214023, 12.272097096384341, 11.279235865220194, 10.45557624282722, 10.213403723979409, 9.4352540004127672, 9.059591562681522, 8.4285767722705884, 8.2189676916739653, 8.2483207838968369, 8.7272133851807343, 9.3853417746714207, 10.436988631508529, 11.737372360714012, 14.160123000210135, 18.959579242189701, 24.173697545434358, 27.627516713161423, 33.591531685669004, 45.297439984936638, 65.897539497304663, 88.610944521286385, 122.51594786973793, 102.9717174324289, 99.771144171631136, 80.254636450195903, 61.009511746666256, 44.786912913807178, 34.271024500323179, 25.65851527359349, 19.184194785168692, 15.369946764821139, 13.384363033613122, 12.174478369218367, 10.905903566261369, 10.181651684644681, 9.6408755952533021, 9.5742304532973588, 8.9763308722412098, 8.5493032688069821, 8.0644484764315241, 8.4276865988342919, 8.7290963397981098, 9.1210926884309647, 10.045472592156113, 11.11042960636143, 12.723025686618605, 15.905816263070196, 20.709491045732211, 24.791931496783274, 30.865414152401307, 41.835714185623189, 55.039237230593706, 69.947250962386548, 86.305656954166849, 70.246077416654487, 70.163190967731239, 60.941585457354833, 46.988681551588527, 37.426824229947464, 31.105114654781531, 24.717957509432271, 18.191253220796238, 14.40571246916225, 12.794837539520243, 11.45286867499644, 10.555233980671725, 10.138595490186225, 9.7617939326440126, 9.2945114067021706, 9.1978519570897408, 8.5269704682665068, 8.3694986947733589, 8.3900850601209847, 8.9019780110016278, 9.4993224051697176, 10.189344774081897, 11.038776504531349, 12.139382952337165, 14.296273465828305, 18.325479810549485, 23.245126565071455, 28.658675517780157, 38.227879918404916, 50.026933715114076, 59.913478880080035, 67.142884728350211, 60.673200962744687, 56.006191670589025, 48.952383807337931, 41.877289900824969, 35.563752651142892, 29.381802508508265, 23.323351055229473, 18.385352146785323, 14.386074460433269, 11.971768024825069, 10.813949701229296, 10.104265892334448, 9.7873837360525968, 9.7876450182420029, 9.3736570608881902, 9.1067837514303367, 8.8782742568475381, 8.4108604655156682, 8.5145291911969032, 9.1444638633789328, 9.5903687569174956, 10.050972673200409, 10.599176866339214, 12.029252308933458, 13.885440808202775, 16.933840539652099, 20.949767453800131, 27.392670445772911, 35.246945811604022, 44.844494652773115, 53.527920120533459, 59.511906457614529, 53.048047546158344, 49.667902591204928, 45.161821768091436, 39.040569237963304, 34.153498370766542, 28.812200326236947, 21.915204963295363, 16.825291542921018, 14.194481006765461, 11.707465544074312, 10.438485210613589, 9.8151245648739263, 9.283134914271356, 9.0821420578250009, 9.0470505359514384, 9.0924089982219396, 8.8466944313913523, 8.4735939062643375, 8.4980702515250215, 9.0113445722107848, 9.2898710387292809, 9.9575678866832327, 10.426127211956205, 11.223822839806607, 13.773251241827809, 16.688691177663944, 20.746301629724744, 26.320191080313084, 33.979868436578599, 42.0788988747152, 48.313306819757166, 52.481937610238305, 49.070609767594966, 49.766679020971374, 46.544272852318578, 40.229040805773415, 32.559839637781486, 26.180018982316895, 21.083604085221012, 16.97356072267355, 13.301411970448795, 11.264067238060504, 10.255581922596939, 9.6536213682112031, 9.2153774062964224, 8.7468110761182825, 8.4356542718375049, 8.7323846743622209, 8.8458599078381237, 8.7323846743622209, 8.4356542718375032, 8.7468110761182896, 9.2153774062964278, 9.653621368211212, 10.255581922596944, 11.264067238060512, 13.301411970448795, 16.973560722673543, 21.083604085221008, 26.180018982316881, 32.559839637781494, 40.22904080577338, 46.544272852318571, 49.766679020971331, 53.048047546158344, 52.481937610238198, 48.313306819757052, 42.0788988747152, 33.979868436578577, 26.320191080313087, 20.746301629724748, 16.688691177663959, 13.773251241827822, 11.22382283980661, 10.426127211956201, 9.9575678866832256, 9.2898710387292791, 9.0113445722107919, 8.4980702515250268, 8.4735939062643357, 8.8466944313913505, 9.0924089982219467, 9.0470505359514402, 9.0821420578250063, 9.2831349142713542, 9.815124564873928, 10.438485210613582, 11.707465544074319, 14.194481006765461, 16.825291542921001, 21.915204963295366, 28.81220032623694, 34.153498370766549, 39.040569237963282, 45.161821768091372, 49.667902591204999, 60.67320096274468, 59.511906457614572, 53.527920120533402, 44.844494652773072, 35.246945811604007, 27.392670445772925, 20.949767453800135, 16.933840539652095, 13.88544080820277, 12.029252308933453, 10.599176866339214, 10.0509726732004, 9.5903687569174885, 9.1444638633789328, 8.5145291911969103, 8.4108604655156718, 8.8782742568475346, 9.1067837514303385, 9.3736570608881848, 9.7876450182420012, 9.7873837360526021, 10.10426589233446, 10.813949701229305, 11.971768024825074, 14.386074460433273, 18.385352146785333, 23.323351055229466, 29.381802508508251, 35.563752651142877, 41.877289900824962, 48.952383807337874, 56.006191670588912, 70.246077416654288, 67.142884728350111, 59.913478880079964, 50.026933715114062, 38.227879918404923, 28.658675517780189, 23.24512656507147, 18.325479810549478, 14.296273465828305, 12.139382952337156, 11.038776504531354, 10.189344774081897, 9.4993224051697229, 8.9019780110016331, 8.3900850601209847, 8.3694986947733572, 8.5269704682665068, 9.197851957089739, 9.2945114067021635, 9.7617939326440144, 10.138595490186225, 10.555233980671725, 11.452868674996441, 12.794837539520238, 14.405712469162248, 18.191253220796249, 24.717957509432264, 31.105114654781534, 37.426824229947449, 46.988681551588492, 60.941585457354776, 70.163190967731154, 102.97171743242879, 86.305656954166807, 69.947250962386491, 55.039237230593784, 41.835714185623168, 30.865414152401325, 24.791931496783278, 20.709491045732214, 15.905816263070196, 12.723025686618604, 11.110429606361421, 10.045472592156107, 9.1210926884309735, 8.7290963397981116, 8.4276865988342937, 8.0644484764315258, 8.5493032688069821, 8.9763308722412134, 9.5742304532973623, 9.6408755952533056, 10.181651684644681, 10.905903566261369, 12.174478369218367, 13.384363033613125, 15.369946764821144, 19.184194785168689, 25.658515273593498, 34.271024500323144, 44.78691291380715, 61.009511746666178, 80.25463645019579, 99.771144171631008, 150.57608076172505, 122.51594786973776, 88.610944521286228, 65.897539497304635, 45.297439984936574, 33.591531685669025, 27.627516713161402, 24.173697545434354, 18.959579242189701, 14.160123000210131, 11.737372360714007, 10.436988631508525, 9.3853417746714172, 8.7272133851807325, 8.2483207838968262, 8.2189676916739689, 8.4285767722705902, 9.0595915626815202, 9.4352540004127654, 10.213403723979411, 10.455576242827226, 11.279235865220191, 12.272097096384336, 14.114901096214016, 16.095362427853331, 21.147383956632702, 27.855857296119162, 38.954469932130408, 55.594423952150585, 84.307080896770429, 119.80060075739968, 151.07347077509721, 218.20881258539524, 172.81330831382382, 114.45554660359451, 76.692136618774427, 51.826480475596433, 35.443188490268959, 29.796933806391745, 27.762087103200372, 21.672488061489549, 16.544820171465698, 13.287947209975444, 11.693775785132251, 9.8383977040140866, 9.1270872628576836, 8.1851518687531897, 8.5854476517741798, 8.5650120282362074, 8.7652013845999939, 9.6259602166533789, 10.393382355075921, 11.004832210877884, 11.637118060780562, 12.699315839076572, 14.281960328360505, 17.135020014783741, 22.221572369125806, 32.725373555898322, 45.982361360241477, 68.705307905169562, 112.62100758585837, 169.13695858014344, 216.90433976724287, 307.01018353832563, 244.24925081430612, 157.1093560387404, 92.156945934839754, 56.412497176172316, 38.781782989560298, 32.110154440123168, 28.237283751764654, 25.031728162044857, 18.893042693096504, 14.406509159208543, 11.695932356823613, 10.702831587760132, 9.4605465058702674, 8.7233150994257205, 8.5960325163902844, 8.8524779931517052, 9.4395532322934272, 9.9738582700772813, 10.53207617814155, 11.166914874614413, 12.321358605090143, 13.970653175227838, 15.815283825259229, 18.317942023831186, 24.692475316947778, 35.843010431224023, 55.468009578938812, 92.303773245695368, 154.32114106589427, 236.16368190120789, 288.66760291638354, 544.93523795856049, 409.50889944188282, 237.18123750436465, 119.9484175480882, 66.578817507603574, 47.964760430736455, 37.131609415223693, 33.103068481604559, 27.569417818570489, 20.771526297093434, 16.674941080314646, 13.235048872238199, 11.323197140674804, 10.488284402310184, 9.7734432560375648, 9.068032997029011, 9.2113480820863387, 9.5692719745256465, 10.325708895539957, 10.85818004844373, 11.899864411981884, 12.765457059373892, 14.471795156401072, 17.743475715717903, 19.778555411511121, 26.952419331651342, 41.705006957506882, 69.155510455435987, 124.78354764285341, 222.45308756649936, 362.69938393816079, 546.63255943550212, 1080.190096085312, 745.27797642375788, 382.14409608709587, 168.44781050782805, 89.418846636419218, 59.087960225901135, 51.536289299313893, 43.037987696240307, 32.657280802203374, 24.858220592568383, 18.8566734045624, 15.131216045563052, 12.585855785544725, 11.471833692521184, 10.776357928238532, 10.225072497047019, 10.113778554702403, 10.422646713468026, 10.39204452752778, 11.686385949826136, 12.224913053368411, 13.273777678692191, 15.572760511446033, 18.579780562544141, 23.060146944216122, 31.897215329939918, 50.717393232745025, 88.946305547570788, 166.56538319335294, 333.23520458041884, 600.04055794738088, 902.6744137040065, 1649.9911393469013, 1194.1693767483896, 627.46721263291568, 257.67571935262464, 130.82046704446515, 85.070809539003903, 64.95962642758694, 56.268912844615201, 47.416694920263595, 31.950664619854226, 23.277252196336782, 17.807798659576147, 14.919521835574468, 13.60275540468532, 12.481114469594081, 11.408229662647638, 11.032765269291394, 11.37290359217827, 11.421939705679749, 11.956364952724407, 12.938917306264656, 13.543616943986398, 17.759823791920141, 22.490208420088635, 25.024628727892473, 35.725757089786228, 60.92155174662291, 113.58561083280503, 235.81544665938503, 495.51661529689198, 1039.2560892520876, 1699.403841506296, 2266.7457843108423, 1439.6385876846614, 877.88627780050422, 436.44287859172971, 201.91826362250967, 120.42091176843758, 92.994178544706443, 78.820512563031016, 59.193773906230312, 43.026808431859891, 32.72703461228695, 24.074134044705318, 19.687568413122339, 17.144679989146301, 15.32830176004647, 13.466217665828177, 12.230328765246888, 12.567429683418707, 13.144389071570203, 13.459360575928054, 14.499148814332536, 16.038244917569731, 20.356277162608194, 26.786032891522385, 31.259021150294068, 41.917980917088691, 66.730731113314405, 139.47307538222745, 266.66937031427761, 673.70954283481024, 1234.659245962343, 2107.8495169895405, 2024.849011267276, 1866.4077869841062, 1132.5867520786235, 527.32695578494054, 275.9549807303485, 182.29010502256293, 128.14673582533564, 94.733059346587495, 75.762437644265887, 54.412823137373678, 40.985078317877729, 32.581505308308067, 27.457440287528364, 22.850207518010212, 19.446116741695526, 17.138709070816009, 15.312687438440404, 15.036681247518965, 14.861292517540395, 14.843145659919811, 15.343856356149836, 19.130394004010878, 24.966146703509185, 32.451856929062608, 37.528009707176075, 50.013293885425036, 89.673379864749862, 155.99683828976509, 316.43682665588739, 611.83627918833008, 1491.9478410044112, 1926.8239401023493, 3926.8611003464503, 2634.734377888386, 1521.0229297033557, 674.50839546082398, 364.56248799144686, 269.99443977661741, 169.37440787587963, 134.18057082099915, 106.99539363556599, 87.291485884508702, 61.577721983396906, 44.812511187792317, 38.103008129659166, 30.767061622228066, 24.399979621456033, 21.348409581307155, 18.55926931620083, 17.794391783727576, 17.409839934721397, 17.67960231321689, 18.718288087401788, 21.222093075405763, 30.846794326525586, 37.908460172493079, 45.504823244138095, 61.496465220623335, 105.17011934648093, 180.13762334768973, 332.57137935864853, 585.56166554295532, 1028.8359899085369, 2444.3527388912739, 15384.540019848309, 8625.0149751320096, 2857.732229334476, 1367.9024800373199, 719.17065363709889, 531.27032729465782, 353.51632607639499, 237.38092500065625, 183.84765943176362, 145.3050149389683, 109.70237013750581, 71.509189516532217, 54.026677368666469, 39.688247321389689, 29.429525390428473, 24.726583224260004, 20.871711217553678, 19.382151607677731, 19.408029208714364, 19.858863545618949, 23.033311632510113, 27.11271894406525, 38.10483512466044, 53.588048009431979, 66.200283784031811, 89.351789047349371, 144.62539946693829, 222.13525013253218, 308.57630937597452, 674.59859086304834, 2216.8794834334685, 11175.737303494294, 13493.224333003214, 8660.703929170164, 6069.217591691604, 3358.0980991141591, 1935.5749815900851, 1219.8477670929722, 664.63525809030602, 430.47940336338263, 254.30481783551838, 177.63416429867087, 121.61239162949126, 77.826621492011839, 54.399342185878204, 38.082981074479164, 29.138683151168006, 23.622937534136597, 22.04888812794956, 20.808198660555991, 22.92975732198266, 25.094853941398917, 30.569538051408909, 39.854538473021513, 58.061322564190043, 75.858347622710383, 106.13580676907651, 142.795979335548, 244.24660034041099, 373.67036323252819, 567.64096645153575, 1571.3337493887695, 4871.835388507573, 16392.519951456165], "imag": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, "imag": [0.0, -1.343901621327011, 1.7218333249394571, 1.3937858919791224, -0.90478486773452105, 0.99560837280124848, 0.13374486301187097, 0.79781790598886937, 0.081235430122062774, 0.47338465820086184, -0.078534664439216209, 0.29500119331967589, -0.62479162763734286, 0.29517170974740459, -1.8712573358384692, 0.14096520515743682, 0.0, -0.14096520515741828, 1.8712573358384816, -0.29517170974737811, 0.62479162763733731, -0.29500119331966146, 0.078534664439223051, -0.47338465820085962, -0.081235430122062774, -0.79781790598886004, -0.13374486301187369, -0.99560837280124848, 0.90478486773452083, -1.3937858919791226, -1.7218333249394557, 1.3439016213270112, 1.5871307700950068, 0.56370110925886319, 1.2603198064423482, 2.9317402483699255, -0.84090846562801591, 0.28395441021550888, 1.0393523450590114, 0.83104911233799517, 0.57334259680943678, 0.30043904784975134, 0.29266755430531527, 0.082930104064024152, -0.36502835420313379, -0.76213554623438207, -0.54316838647957832, -0.66028179246101926, 0.20250391456576916, -0.059393044635645806, 0.9730252602081445, -0.01782962412134724, 0.18933865423621773, 0.011846708470834844, 0.014269719747672759, -0.16281560265463554, -0.24410560657714067, -0.34963232308526365, -0.32234943120999277, -0.6691127915152365, 0.14162995430492473, -0.96909874957635578, -1.0446282489208472, 1.3612809629252947, 0.70464907960903156, 1.644656498115789, 1.457626495523533, -0.6840361527612806, -0.68147470897196205, 0.73118573872942638, 1.2398159392319572, 1.2736175855463174, 0.65550663858936109, 0.45897332874861019, 0.24948129249627171, 0.22959418801218681, -0.67355449112520815, -0.78917379441405089, -0.83366315887398401, -0.022795944141155908, 0.15753196960053301, 0.50836195321582456, 0.16292155588892338, 0.14047533322536612, -0.029004083819245446, 0.042892010327960466, -0.084780328674450303, -0.13217055914433321, -0.30971168558439593, -0.44061350110378666, -0.59380988893293396, -0.46719027726302742, -1.2068396187631398, -0.94969408752542273, 0.28960616895899027, -0.40038203546320883, 0.55620249274438927, 0.14065731983379595, 1.390858139437829, -0.67833342052765033, -2.0926714021100303, 1.5120873610545131, 1.7251253398279642, 1.0499886203647506, 0.6413261932720733, 0.43315407205712014, 0.32584150801434064, 0.17968537990933128, -0.27887473465086693, -0.72646685328097149, -0.95200274282553721, -0.38747270801185796, 0.31172559850957127, 0.3167238274644899, 0.25925756638025876, 0.20404659434537992, 0.19201780058118942, 0.0025274006864131854, -0.20254250312507205, -0.30146681989018398, -0.45156398868756126, -0.71583959398165176, -0.78636343748881254, -1.3089140888320636, 0.41393919421363429, -0.92308378001125058, -0.96086768695066838, 0.99489496750979101, -0.39381772120025876, -0.59264738582705545, -1.1161873115341066, -0.38954314931717815, 0.15776662530084892, 0.39428020304678685, 1.7089894615168126, 0.92868491010646026, 0.350973234838865, 0.32125288047075679, 0.23669588269358602, 0.19592809956108717, -0.32646969875093118, -0.5617213523679232, -0.76596137582739887, -0.4848448097808058, 0.0035984972735877987, 0.36234031342999445, 0.02584515038201295, 0.38132345964413561, 0.29172235696712762, 0.10460468586389386, -0.14942504329756173, -0.23406484190568319, -0.59702583417138988, -0.64791716309340031, -1.1972380302487249, -0.37081056871520385, 0.31028937719203986, -0.11181138224479611, 0.050224150041196312, -0.59242819035998506, -0.4703174557375312, -0.55567065990801068, -0.18073895078348806, -0.07286620180284567, 0.9080135456055427, 1.3448759333642528, 0.56391746646770735, 0.68824563709365816, 0.27054447491545203, 0.18555506064923469, 0.21628766080629094, 0.089227696371118645, -0.1330464412169812, -0.51044075389216448, -0.57482956329522528, -0.49992842173388902, -0.0088882269741786465, -0.026119953804469152, 0.16649307508904443, 0.18021017930701566, 0.24129551707937363, 0.051414766220490717, -0.054409237913703348, -0.21354052285348152, -0.34065461819243753, -0.83215952955458006, -0.65455151443418824, -0.51711748874533736, -0.34193161983249143, 0.25991967791421511, -0.31922838969311168, -0.25096233231074472, -0.019538834540300996, 0.23386297177233037, 0.39736821511961345, 0.72603823847118543, 0.8168297825020705, 1.0763583447726135, 0.80611388127401196, -0.023349196100613815, 0.09911450222934183, 0.11136104528246665, 0.11912886110102348, 0.1510441091194204, -0.10294245821664937, -0.34339582552474845, -0.46071517851602106, -0.27454988477746128, -0.021291149638310566, 0.035818571775215378, 0.08165077597052256, 0.19593247304969552, 0.16021249355536668, 0.073411101811457824, -0.094413661460628892, -0.14774431296435828, -0.3788677957633147, -0.38953231590523063, -0.62892221563014938, -0.80219655994233319, 0.19764870263775744, -0.28620163890053252, 0.041738811618785537, -0.047502831671303453, 0.42759360787810174, 0.44353340330745888, 0.61232760985457502, 0.70207820166905466, 0.72192603666395361, 0.63928791416871467, 0.32589261306510392, 0.17268517789923485, -0.11832252230060911, 0.049391580905339259, 0.12347320818582375, 0.064697981434016727, -0.0067102091687047593, -0.27370819742589092, -0.26939466873920287, -0.25176222060998998, 0.0081310228529853402, 0.013089808045657288, 0.185105158015618, 0.11713478908555769, 0.20974664985683705, 0.048745778368168922, -0.0069359534603168627, -0.18083894269546222, -0.15907734702039189, -0.30490326972545728, -0.49525043110711681, -0.20633977154543071, -0.28929558460304977, 0.017255660868048289, 0.089250558912840242, 0.34349774583430781, 0.38277205043588519, 0.56772102070325814, 0.48811475177501612, 0.43659704639067964, 0.38858095560717315, 0.35469216713625445, 0.11095078972724678, -0.050354098986562471, 0.010426074631633324, -0.023032548262322892, 0.068558023511485697, 0.1199625535177883, -0.0044347552049944488, -0.14323623171577524, -0.2742488887234425, -0.15059247824174998, -0.18695930341705305, 0.068136880204083733, 0.049162011568203073, 0.16362595169702937, 0.078325720975578078, 0.11153460033842358, -0.10346450949096449, -0.039250261082814021, -0.12027128091306288, -0.23113762225192683, -0.07885965728132803, -0.25982980909746989, -0.11219308206942701, -0.19372413199270677, 0.1843821378255476, 0.22544536191476028, 0.22838873802709003, 0.11856768480607156, 0.20011010622694389, 0.16443500436191505, 0.13441036107213369, 0.15600995797323181, 0.029984980792926118, -0.017950425588896139, -0.095997178609832007, 0.045815559155676679, 0.050600401411128038, 0.062555070152262096, 0.095362679574156284, -0.25884721081868178, -0.19352665675464722, -0.35758313393883445, 0.057241994101256972, -0.074489797691212556, 0.15760636776309309, 0.077110294424560424, 0.21793581355625327, -0.085384156728057797, 0.054278095118519401, -0.063513516552257218, -0.13925876777247995, -0.053090482149624127, -0.088550664314818353, 0.097053178020475875, 0.021982682491266895, -0.032057043772553598, 0.011606155190390474, 0.20297398899260674, 0.087137707937781944, 0.14235450986309658, -0.0066344290522644652, 0.069443551113326479, 0.037974760198136455, 0.02723212012252952, -0.030406089689394462, -0.062891439826322393, -0.033462702565164315, -0.13091343010011852, 0.18662643567548856, 0.14580799486104234, -0.055902479292620696, -0.12813917376287484, -0.43706544906097528, -0.097926711565598012, -0.17130544619864244, 0.30727693305691806, -0.2953735036217468, 0.22056405751825581, -0.16529753202987246, 0.34047455883436423, -0.0020829156573361745, -0.14775664020163451, -0.074518551996532553, -0.14850320507536691, -0.057283952975165524, 0.064665418468806371, 0.23511892751228192, 0.11870789775778592, 0.0625549454189092, -0.031298024311669038, 0.0064402102662661845, -0.079082063664820082, 0.0553675188881671, 0.025089334126077875, 0.063800091749494747, 0.015988824093570549, -0.048545356942180602, -0.11558972094679192, -0.097487343329346582, -0.00049492928166827938, -0.18157168721539482, 0.33942288969170681, 0.30522386319855388, -0.54639589537101751, 0.013830571830100688, -0.38850732372171304, 0.223075916963486, -0.31865288748874443, 0.25745374979707547, -0.81885347149607268, 0.50855168910609783, 0.20195166784419227, -0.069785743671890318, 0.050884179954965268, -0.18377202523753672, -0.17303565260659334, -0.15924096907201105, 0.092146932205194695, 0.26379835890783954, 0.35934781180257075, 0.24329566614587647, 0.18564601960109514, 0.24223475523068938, 0.28967646870392649, 0.15967526911861643, 0.22365254766781259, 0.12185893624650075, 0.070105150765122456, -0.045495054958742086, -0.13016925945874308, -0.12331107551299306, -0.096585989040329473, 0.22919036999990153, -0.39865574246315477, 0.3675558227752026, 0.45597314967635888, -0.93227329954448013, 0.32110784375924512, -0.59160024067907024, 0.39469180639349893, -1.5810043928614619, 0.71926168426275527, -0.095920343398476549, -0.084877145937512208, 0.40399573862898452, -0.041743397742653787, -0.10303861084867461, -0.21220001866353477, -0.17318676777643013, 0.015792824133968684, 0.37994714934071439, 0.4573031649737595, 0.48491463219236358, 0.278045596991639, 0.33099751793182702, 0.073461819418744675, 0.30445717300203839, 0.34815332433352569, 0.27336662525128713, 0.050113511094951838, -0.059405844632707393, -0.16261143474742229, -0.10891422081227681, -0.075070297171157971, 0.15348950297637756, 0.63132994551184463, -0.93604459898777692, 0.42325697055405742, 1.1782712312395338, -1.6700526640562563, 0.67307548012867768, -2.0365876460209691, 0.69665186973993298, -0.3990336398782886, -0.77751104168880314, 0.5628167675812451, 0.044468991083659551, 0.081278132745909976, -0.070993394358423723, -0.091852860366211056, -0.10956166698887557, 0.072941493094735049, 0.35167328896657912, 0.63488067579895824, 0.56199684605205058, 0.60627707274502263, -0.018970197624819329, 0.37668655809680845, 0.22317003729918369, 0.17596357423928666, 0.18202109891962201, 0.038230887694388455, -0.10475039750351335, -0.15123171993121909, -0.11278544639649389, 0.054110895720760221, 0.23852479976301316, 0.30802887330694756, 0.9587111952302938, -1.4503952980081434, -0.061458391503499687, 2.1389978076683529, -3.5142346308062073, 1.6764175122847729, -0.63252660757872936, -1.1970497991660358, 0.46027494328631885, -0.38868529039399474, 0.063066133758830412, 0.085024034480974467, 0.047630852690686495, -0.022775175196356016, 0.011847923870825085, 0.13095780969018636, 0.39870906011464996, 0.54100133271840178, 0.62859642257501991, 0.12859489098386362, -0.17371304193339687, -0.77135662792152559, -0.35696964892202537, -0.15350152251372465, -0.081108427029901481, -0.12890403553126065, -0.14048874638904316, -0.18368108831833277, -0.13143015225876853, 0.037318187986402011, 0.39815350666418903, 0.44685069559938601, 0.29854169122462937, 1.0193515924216456, -2.1696976832709978, -1.5476694374721527, 4.1115013597100623, -1.6066056279080296, -1.1129988398741206, 0.94641351902713078, -0.66464487182371434, -0.22405122030787003, -0.16845135046504239, 0.10281428632219461, 0.12776142939872864, 0.10903691838775739, 0.097477148948306502, 0.19622512970077355, 0.37581890109000032, 0.58152174984392169, 0.57588244896775898, 0.56835821427965394, 0.0, -0.075660591040101741, -0.5319087481065502, -0.49066506338765509, -0.25388521774435541, -0.17759782612549996, -0.17342105765811164, -0.17172379067488869, -0.17163733194043357, -0.0096917261033460274, 0.23320936684713756, 0.61544163247884909, 0.23031783223894131, 0.4504667163321609, -0.25422186154383031, -1.2338903723029762, 0.0, 1.2338903723029679, 0.25422186154382892, -0.45046671633216256, -0.23031783223894184, -0.61544163247884998, -0.2332093668471385, 0.0096917261033450074, 0.17163733194043357, 0.17172379067488885, 0.17342105765811225, 0.17759782612550043, 0.25388521774435546, 0.49066506338765642, 0.53190874810654976, 0.075660591040102837, 0.17371304193345374, -0.56835821427958444, -0.57588244896772234, -0.58152174984391114, -0.37581890108999716, -0.19622512970076919, -0.097477148948305531, -0.10903691838775691, -0.12776142939873006, -0.10281428632219579, 0.16845135046504131, 0.22405122030786861, 0.66464487182371257, -0.94641351902713855, 1.1129988398741191, 1.6066056279080259, -4.1115013597100578, 1.5476694374721462, 2.1696976832709982, -1.0193515924216472, -0.29854169122463242, -0.44685069559938834, -0.39815350666419097, -0.037318187986404149, 0.13143015225876697, 0.18368108831833294, 0.14048874638904385, 0.12890403553126736, 0.081108427029909794, 0.15350152251374438, 0.35696964892206617, 0.77135662792155557, 0.018970197624865574, -0.12859489098383081, -0.62859642257500314, -0.54100133271839856, -0.39870906011464824, -0.13095780969018578, -0.011847923870824886, 0.022775175196357359, -0.047630852690687051, -0.085024034480974689, -0.063066133758832563, 0.38868529039399197, -0.46027494328632251, 1.1970497991660347, 0.63252660757872814, -1.6764175122847713, 3.5142346308062082, -2.1389978076683618, 0.061458391503496308, 1.4503952980081398, -0.95871119523029713, -0.30802887330694667, -0.2385247997630148, -0.054110895720761498, 0.11278544639649313, 0.15123171993121864, 0.1047503975035141, -0.038230887694384028, -0.18202109891961496, -0.17596357423927469, -0.2231700372991606, -0.37668655809676688, -0.33099751793180054, -0.60627707274499654, -0.56199684605203359, -0.63488067579895024, -0.35167328896657507, -0.072941493094733675, 0.1095616669888755, 0.091852860366211805, 0.07099339435842339, -0.08127813274590924, -0.044468991083660724, -0.56281676758124888, 0.77751104168880336, 0.39903363987828711, -0.69665186973992854, 2.0365876460209633, -0.67307548012867546, 1.6700526640562572, -1.17827123123954, -0.42325697055405653, 0.93604459898777859, -0.63132994551184651, -0.15348950297637845, 0.075070297171156097, 0.10891422081227575, 0.16261143474742193, 0.059405844632708517, -0.050113511094947612, -0.27336662525128258, -0.34815332433351742, -0.30445717300202041, -0.073461819418725136, -0.24223475523068566, -0.278045596991634, -0.48491463219235625, -0.45730316497375445, -0.37994714934071183, -0.015792824133968708, 0.17318676777642958, 0.21220001866353638, 0.10303861084867472, 0.04174339774265514, -0.40399573862898636, 0.084877145937511819, 0.095920343398476396, -0.71926168426275106, 1.5810043928614652, -0.39469180639349549, 0.59160024067907002, -0.32110784375924395, 0.93227329954447902, -0.4559731496763611, -0.36755582277520088, 0.39865574246315755, -0.22919036999990147, 0.096585989040328876, 0.12331107551299252, 0.1301692594587433, 0.045495054958742745, -0.070105150765121721, -0.12185893624649984, -0.22365254766781101, -0.15967526911861513, -0.28967646870392255, -0.0064402102662611729, -0.18564601960108706, -0.24329566614587189, -0.35934781180256659, -0.2637983589078372, -0.092146932205193599, 0.15924096907201093, 0.17303565260659445, 0.18377202523753683, -0.05088417995496395, 0.069785743671890568, -0.20195166784419286, -0.50855168910609727, 0.81885347149607191, -0.25745374979707175, 0.3186528874887527, -0.22307591696348483, 0.38850732372171648, -0.013830571830101068, 0.54639589537102196, -0.30522386319855499, -0.33942288969170481, 0.18157168721539607, 0.0004949292816654663, 0.097487343329345874, 0.11558972094679136, 0.04854535694218131, -0.015988824093569717, -0.06380009174949508, -0.025089334126078205, -0.055367518888164394, 0.079082063664820748, -0.087137707937777115, 0.031298024311669809, -0.062554945418906743, -0.1187078977577844, -0.23511892751228033, -0.064665418468805788, 0.057283952975166051, 0.14850320507536743, 0.074518551996532748, 0.14775664020163706, 0.0020829156573371555, -0.34047455883436323, 0.16529753202987046, -0.22056405751825445, 0.29537350362174863, -0.30727693305691484, 0.17130544619864224, 0.097926711565602301, 0.43706544906097711, 0.12813917376287753, 0.055902479292622584, -0.14580799486104357, -0.18662643567548898, 0.13091343010011883, 0.03346270256516358, 0.062891439826322837, 0.030406089689394671, -0.027232120122529239, -0.037974760198136892, -0.069443551113326632, 0.0066344290522653612, -0.14235450986309153, -0.22838873802708684, -0.2029739889926063, -0.011606155190393872, 0.032057043772551565, -0.021982682491265355, -0.097053178020474681, 0.088550664314818187, 0.053090482149624633, 0.13925876777248047, 0.063513516552260146, -0.054278095118517826, 0.08538415672805795, -0.21793581355625447, -0.077110294424558939, -0.15760636776309189, 0.074489797691215873, -0.057241994101257694, 0.35758313393883939, 0.19352665675464989, 0.25884721081868634, -0.095362679574154757, -0.062555070152259029, -0.050600401411129904, -0.04581555915567774, 0.095997178609832215, 0.0179504255888953, -0.029984980792924155, -0.15600995797322956, -0.13441036107213314, -0.16443500436191227, -0.20011010622694014, -0.1185676848060674, -0.38277205043588519, -0.22544536191476122, -0.18438213782554774, 0.19372413199270522, 0.11219308206942717, 0.25982980909746933, 0.078859657281328127, 0.23113762225192824, 0.12027128091306288, 0.03925026108281611, 0.10346450949096535, -0.11153460033842286, -0.078325720975578092, -0.16362595169702793, -0.049162011568201006, -0.06813688020407746, 0.18695930341705305, 0.15059247824175628, 0.27424888872344572, 0.14323623171577923, 0.0044347552049957819, -0.11996255351778637, -0.068558023511484198, 0.023032548262321723, -0.010426074631633324, 0.05035409898656467, -0.11095078972724633, -0.35469216713625396, -0.38858095560717304, -0.43659704639067964, -0.48811475177501651, -0.56772102070325836, -0.42759360787810491, -0.34349774583431025, -0.089250558912845224, -0.017255660868053985, 0.28929558460304938, 0.20633977154543065, 0.49525043110711575, 0.30490326972545867, 0.15907734702039197, 0.18083894269546641, 0.0069359534603191542, -0.048745778368166057, -0.20974664985683683, -0.11713478908555536, -0.18510515801561447, -0.013089808045645399, -0.0081310228529865667, 0.25176222061000281, 0.2693946687392072, 0.27370819742589902, 0.0067102091687078758, -0.064697981434012536, -0.12347320818582255, -0.049391580905336498, 0.11832252230060807, -0.17268517789922963, -0.32589261306509937, -0.639287914168712, -0.72192603666395261, -0.70207820166905555, -0.61232760985457479, -0.44353340330746277, 0.019538834540299969, 0.047502831671299713, -0.041738811618788646, 0.28620163890052958, -0.19764870263775797, 0.80219655994233652, 0.62892221563015416, 0.38953231590523568, 0.37886779576331475, 0.14774431296436413, 0.094413661460632833, -0.073411101811453133, -0.16021249355536557, -0.19593247304968892, -0.081650775970513151, -0.035818571775193361, 0.021291149638310077, 0.27454988477748077, 0.46071517851603194, 0.34339582552475939, 0.10294245821665464, -0.15104410911941546, -0.11912886110102075, -0.1113610452824651, -0.099114502229341719, 0.023349196100618572, -0.80611388127400851, -1.0763583447726097, -0.8168297825020715, -0.72603823847118742, -0.39736821511961523, -0.23386297177233123, 0.4703174557375292, 0.25096233231074278, 0.31922838969310624, -0.25991967791422005, 0.34193161983249054, 0.51711748874534202, 0.65455151443419168, 0.83215952955458772, 0.34065461819243747, 0.21354052285348993, 0.054409237913707532, -0.051414766220484194, -0.24129551707937194, -0.18021017930700467, -0.16649307508903324, 0.02611995380449984, 0.0088882269741778191, 0.49992842173391339, 0.57482956329523827, 0.51044075389218102, 0.13304644121698647, -0.089227696371110179, -0.21628766080628836, -0.18555506064923225, -0.27054447491545125, -0.68824563709364495, -0.56391746646770291, -1.3448759333642519, -0.90801354560554559, 0.072866201802841646, 0.18073895078348665, 0.55567065990800857, 0.39381772120025821, 0.59242819035998318, -0.05022415004119711, 0.11181138224479518, -0.31028937719203997, 0.3708105687152074, 1.1972380302487284, 0.6479171630934083, 0.5970258341713901, 0.23406484190569171, 0.1494250432975659, -0.10460468586388726, -0.29172235696712562, -0.38132345964412656, -0.02584515038199919, -0.36234031342997136, -0.0035984972735878563, 0.484844809780835, 0.76596137582741564, 0.56172135236794307, 0.32646969875093934, -0.19592809956107882, -0.23669588269358457, -0.32125288047075495, -0.35097323483886461, -0.92868491010645182, -1.7089894615168075, -0.39428020304678818, -0.15776662530084981, 0.38954314931717637, 1.1161873115341048, 0.59264738582705467, -0.5562024927443876, -0.99489496750979078, 0.9608676869506646, 0.92308378001124303, -0.4139391942136339, 1.3089140888320632, 0.78636343748881332, 0.71583959398165975, 0.45156398868756092, 0.30146681989019147, 0.20254250312507779, -0.0025274006864056797, -0.19201780058118689, -0.20404659434537004, -0.25925756638025355, -0.31672382746446043, -0.31172559850957088, 0.38747270801190198, 0.95200274282556019, 0.72646685328099436, 0.2788747346508717, -0.17968537990931685, -0.3258415080143367, -0.43315407205711348, -0.64132619327207141, -1.0499886203647348, -1.7251253398279605, -1.5120873610545145, 2.092671402110029, 0.678333420527655, -1.3908581394378237, -0.14065731983379445, -0.704649079609032, 0.40038203546320611, -0.28960616895899161, 0.94969408752542028, 1.2068396187631385, 0.4671902772630277, 0.59380988893293563, 0.44061350110379105, 0.30971168558439566, 0.13217055914433939, 0.084780328674454217, -0.042892010327952576, 0.029004083819250629, -0.14047533322535441, -0.16292155588890503, -0.50836195321575572, -0.15753196960053273, 0.022795944141262008, 0.83366315887402342, 0.78917379441409119, 0.67355449112521315, -0.2295941880121686, -0.24948129249626083, -0.45897332874859181, -0.65550663858936065, -1.2736175855462932, -1.2398159392319503, -0.73118573872941783, 0.68147470897196305, 0.68403615276129381, -1.4576264955235338, -1.6446564981157907, -1.587130770095007, -1.361280962925294, 1.044628248920846, 0.96909874957635356, -0.14162995430492481, 0.66911279151523617, 0.32234943120999365, 0.34963232308526704, 0.244105606577141, 0.16281560265463482, -0.014269719747667968, -0.011846708470826806, -0.18933865423621368, 0.017829624121361812, -0.97302526020811309, 0.059393044635753067, -0.20250391456577049, 0.6602817924610479, 0.54316838647960974, 0.76213554623441593, 0.36502835420313695, -0.082930104064008442, -0.29266755430530383, -0.30043904784973247, -0.57334259680943656, -0.83104911233797829, -1.039352345059011, -0.28395441021550044, 0.84090846562801602, -2.9317402483699175, -1.2603198064423498, -0.5637011092588653], "height": 32, "width": 32, "top": {"real": [7376.2586563979694, 2820.8522261533344, -5690.9497085395169, -4380.4513376296127, 388.6547315083493, -903.73852292232903, -92.233688854560938, 42.706854165930984, -25.699252344104405, 18.326281263011424, -9.8542387499157442, 37.933510453511431, -7.8666684309670254, 29.305075180843883, -62.728617354661907, 59.853240627388928, -83.636374399044939, 59.853240627386725, -62.728617354660756, 29.305075180844177, -7.8666684309671826, 37.933510453511431, -9.8542387499160657, 18.326281263013065, -25.699252344104405, 42.706854165930693, -92.233688854557883, -903.7385229223205, 388.65473150834708, -4380.4513376296045, -5690.9497085395114, 2820.8522261533035, -14789.032757371706, -31341.185848117759, -11739.029695619593, -866.78466658113393, -286.68893204831016, -488.33133225087568, -120.74864963582691, 33.024718339611894, 4.6010715329524166, -5.4456740109959769, 4.1139030671062393, 6.8633911029603265, 13.634662924808397, 1.860738389345556, -0.35796846682150263, -16.625069210775109, -14.807142725465134, -11.664010623894898, 0.82076815557037341, -12.155673848711613, 29.148348103795893, 15.708115389403888, 15.457503110539639, -12.683126140048403, 18.961120931859309, -53.601071866945205, 34.13442733448511, -942.72760774202175, -816.05094481379797, -3798.5470080468217, -8876.4278731043305, 2890.8320511623124, -27232.421940549859, -15984.118662988652, -6745.5682223604126, -409.0173927639226, -42.225544009425754, -526.29225504527119, -49.325655652906178, -23.20435390028933, 15.630174176003303, -14.835941429972856, 6.4726407463551725, 1.9958373496876167, 13.291351981828791, -1.1889020916883648, 10.030112101155916, -20.674418420848429, -0.1180338988714332, -18.745207674116038, -1.9041291470929314, -7.8710994817942428, 20.385879097388031, 14.857882580655026, 16.950197274587541, -2.033086164989526, 25.446982376231624, -20.624592190095939, 89.303256695009566, -637.06671072120571, -1428.8365567175017, -3113.2951228079887, -8881.5889351276946, -20696.443600183225, -9007.3131796291163, -5346.3306090002807, -1439.1675666408983, -1292.5663502582722, -654.09944760377277, -301.70227472838963, -97.440927522083996, -3.3329385156017457, -8.2027965131620153, -5.7742639473870323, -0.51069538319854635, 6.6402880636318233, 8.0126520330591013, 4.0882931682311625, -1.5210347396074151, -2.2335756869246048, -2.3288001542285928, -8.1315470965115093, -6.9403489082507335, -3.6748649027980562, 4.5937507427443531, 9.4957008241539675, 10.21057060635477, 11.676825499792781, 23.758858737501662, 20.288236987689611, 20.969766857125101, -333.51324485695324, -810.1584353673162, -2278.6820540988888, -5194.4090916202094, -7188.339898707115, -3154.0647285180166, -3175.2001764803263, -3085.2108289353805, -1226.2977859942123, -836.26206181721784, -301.22169990911055, -15.123294918373254, -11.324269109892445, -4.9992210979693361, -6.7152536525607003, -1.7063097661780617, 2.7727430774894133, 7.1357881321682148, 4.8669904821248817, 2.7536893948953045, -2.039934630711028, -0.54898928719761031, -2.9485150656567436, -3.1026085263061316, -1.3614545093317587, -0.47304719477880164, 4.0409079710738931, 6.2483313015723709, 8.8540252326547968, 18.802141623354618, 24.681352644731625, -6.1114787305120313, -148.49350684902615, -572.34886359995278, -1663.447165744248, -3308.6059534728529, -4493.8282478995061, -3720.0991836319017, -3323.4413322476125, -2223.3326338682518, -1347.5218927383742, -402.86472873295924, -179.07785660663581, -42.171627369153548, 13.024576486879551, -6.680339136637107, -4.3910663300422375, -1.9843645852297727, 2.2081255979025864, 3.4357055351949564, 4.8667397989738959, -0.0088018165553918268, -1.2900498248809058, -2.2230696316747851, -0.79613102395178958, -2.1003135650570957, -0.12607857026311342, 0.39928478103703097, 3.6189445517483683, 5.5230847011321824, 8.0807144322734796, 9.6537837209938342, 14.647883903859292, 41.761171107733581, -2.5048396429986766, -288.68438568639834, -1058.7748743456032, -2111.6481404168812, -2853.5518466904464, -2206.4742388833251, -2207.4882969003411, -1427.0600295974209, -607.97311048711572, -283.75431628355585, -25.222873313531675, -11.21903778037437, -8.3653822807976148, 1.9503276148055415, -4.8653324939916756, -1.0558165909594248, 0.18037835028025548, 3.7232300368774149, 2.6260151255238955, 1.901747650437611, -2.157984497724208, -1.9174482807780493, -2.671365924063863, -1.6822520344543337, -1.7835969801453151, 0.50013744114287995, 2.3612666445612516, 3.9892729615099713, 5.0282597338489978, 4.3849613992965208, 7.5675449146459712, 26.738330274504211, 30.293009711209869, -72.365238046453968, -397.03998761961856, -1129.7740740979677, -1877.1475049531305, -1015.8170093531105, -757.33856581211364, -464.04986763483271, -221.8093865245269, -66.670140953208971, -32.688631361988044, 12.245388279073117, -2.09491941357581, -4.9214217924247272, -0.76863227802399692, -1.3611346863617457, 0.59203001318831738, 2.3940543871676745, 2.1386535528573569, 0.12698658321630202, -0.36055346759858353, -1.509631760657328, -0.37490601010058217, -1.4076418772333708, -0.9909556252137327, 0.025704615807234758, 1.7006413859335137, 2.4101241863335301, 2.7230446482899899, 4.8975474263529586, 5.5799176422335819, 7.3823613929523022, 3.3531473980019464, -7.0947231292717241, -129.95057654055665, -447.33256260550814, -801.49855315264267, -297.71280467242298, -206.25666421819534, -105.85720001872387, -40.587336797783266, -19.033005269759897, 0.68623525685047748, -2.224754122223199, 2.8686428937153541, -1.3772856826259623, -3.7001182860738422, 0.56033030111271653, 0.31649576481826269, 2.6435570864758695, 2.0729019173667451, 1.3896745329326239, -1.1733083260901866, -0.5564284933709801, -0.75319429727185228, -0.16542223594150865, -0.34981510145177358, 0.63947101827441399, 0.64495136714903334, 1.5357324193076427, 2.8130362926385226, 4.2845415972556511, 3.3415574201637508, 1.8747628167747612, 4.3092019872505318, 1.2280777000161127, -53.113548920223614, -153.12296297306594, -262.59106890839882, -78.507966752468647, -44.713956729763169, -24.503242066462356, -11.422688403833899, -6.8324882374101472, -1.9389416594991726, 2.3416848160301655, -1.9581076589813595, 0.65441624739147974, -1.6215192060328243, -3.0886268279599727, 2.5574024719687731, 2.5859312055320047, 3.582324284715654, 2.3851708003164411, 0.41148866509477899, -1.7870448740379692, 0.35867475608513461, -0.33833230756510146, 0.72134917870519866, -0.6451343286468626, -0.036676325952571386, 1.8806785141996829, 3.0459530955995717, 3.0361534846987768, 1.8248078988499765, 2.0278729258957484, 1.4298438030821454, 1.444150632114136, -16.470770345398158, -59.038387406403551, -81.624687331130716, -41.65158074717607, -20.04858945819673, -26.796686775361913, -13.274828683787613, -5.4849410264677898, -1.7682797425272327, 0.51276950901621554, 0.11701161368370352, -2.1525887062494076, -0.35663258536781867, -0.3882693848937121, -2.5566610369397282, 6.5139128293303541, 4.1537964250038213, 2.1778507792771808, 0.67618078211604771, 0.67271760119275736, -0.033223712968477391, 1.7891382881043305, 0.76549003012129324, -1.8840461990319548, 0.62764371336603419, 2.9500669911593334, 2.7998549031057309, 2.7114466415734606, 2.640497491382968, 2.5087646302130264, 0.3613507110411972, -0.85975587552377708, -5.8993661554875496, -19.483874473281006, -46.656198806892085, -42.447185502078959, -49.475399988407517, -28.809040802191632, -23.364991274814315, -9.3720360014436217, -2.193549633075087, 0.73237790527535407, -0.60287176352488314, -1.0672361910519921, -3.8374461384037328, 0.81264916316705227, 3.7024783733454592, -2.2274307006673761, 12.171761179906976, 2.6484997231173142, 2.1242155181960216, 1.3916869062247985, 3.2969221728544564, 0.14147762112210061, -1.2659388719629114, 1.2983452624583423, 2.4646463691490883, 1.8819178233751466, 2.4699197782172551, 3.4105600445409938, 3.476877662739366, 1.970602508408192, 0.54569379437115284, -0.61503488296082265, -1.584537886719974, -9.9108851118855448, -17.968893142628509, 0.83280248396047341, -16.804391379900096, -34.94566909960443, -18.287776152890835, -10.19076861522101, -2.3659355754892188, 0.38430207827148843, -0.39175809943359002, -2.0786620762717161, -2.1066158183404973, -4.0664292874937455, 5.1056392812224729, 9.4711485160420441, -2.577327881134023, 13.563457368659881, 5.4454444116217209, 3.0810071966308192, 3.71409458782056, -2.7028450063821525, 1.3893389983926554, 5.0295677911697716, 2.0512141006459763, 0.57687440104333121, 3.1399874906463121, 4.3044361492909653, 3.2286926050373941, 2.0294536549063102, -1.2648718051872285, 0.81982278416560761, 5.1937461344056883, 4.5660530070978913, -2.5751982814314349, 10.842808636227273, -2.7634934625905125, -4.7410363255116081, -13.231297049532179, -1.2392069417375811, 0.57695677140757395, 1.6572428400911725, -0.80487566941313171, -2.1976406054664008, -2.9337382106235617, -0.97696074341019468, -0.52039830113476337, 8.6932746740913363, 18.008772340009436, -6.3894706762347742, 16.656214263229323, 5.6480128145892756, -2.6420757676626341, 6.9832538798624686, 8.7117753152913409, 2.944231966553327, -0.4194257648970201, 1.9464964791115495, 2.8117062685871983, 3.586606822562346, 3.6410340078674697, 2.1750438821370093, 0.92746135521083628, -0.10515597578954901, 10.916870623975495, 22.695665816423684, 26.922229107827643, 56.766547106622895, 34.346379481364103, 8.0209890672230593, 4.2276867493724737, 2.0129650120803855, 4.5721973792850843, 1.7026114564114816, -0.10193290249876912, -1.9877285465951933, -1.4091416884493304, -1.0407935339809435, 4.6927135104950821, 4.461866465996299, 13.14685287043932, 23.490161756370188, -8.5448676586843408, 5.7639040157862622, 6.7103974452754578, 11.090643871879529, 8.6117964749879281, 0.95828374077, 2.8859143014026349, 1.6803595914225926, 1.5724318132354216, 2.3664842485474589, 3.3844724421254813, 2.8754034241014996, 0.81420697505183059, 3.4592633965370854, 8.5418678989119741, 26.408180301615612, 46.472433795130719, 15.748155014076561, 18.939305693948317, 18.355423454953527, 7.4385947937208252, 6.8421582395566496, 3.8284829441665011, 2.3084513014265475, 0.014383588041353481, -0.70170837590258439, -0.9519271341888802, 0.60488638387036897, 4.7316501322768225, 11.640012070252736, 11.084534895769892, 7.9427552323827459, 17.090358547431578, -4.0062765415124373, 18.903471314461139, 10.522914852414784, -0.31727085194714888, 8.1914062421872558, 5.5678768421215832, 1.9515757898119943, 0.13204460835173984, 0.94686979396028259, 2.2566984460298456, 3.5582629306092612, 2.7150237945754423, -0.038191379039643265, 6.3090662882758242, 18.822434458740766, 20.387061197767007, 34.210424610695611, 25.78016711295075, 13.422328957688888, 4.2578162502068047, 2.3122896178400203, 4.1892680028045364, 2.2592478817037658, 0.95372259584702335, -0.47653425908969432, 0.069240128911087023, 0.91781179579556105, 6.6018505246928774, 10.152138543237006, 11.530381877152177, 5.6211304516254526, 4.7909070280413646, 35.867989503906792, 4.7909070280413584, 5.6211304516254428, 11.530381877152157, 10.152138543237017, 6.6018505246928978, 0.91781179579556504, 0.069240128911081333, -0.47653425908969432, 0.95372259584699437, 2.2592478817037476, 4.1892680028045142, 2.3122896178399985, 4.2578162502067478, 13.422328957688848, 25.780167112950689, 15.74815501407811, 20.387061197767874, 18.822434458741888, 6.309066288276548, -0.038191379039562191, 2.7150237945755555, 3.5582629306093274, 2.2566984460298958, 0.94686979396029702, 0.1320446083517558, 1.9515757898119841, 5.5678768421215636, 8.1914062421872096, -0.31727085194714522, 10.522914852414772, 18.903471314461111, -4.0062765415123884, 17.090358547431556, 7.9427552323827513, 11.084534895769902, 11.640012070252746, 4.7316501322768385, 0.60488638387038995, -0.9519271341888812, -0.70170837590259716, 0.014383588041281273, 2.3084513014264751, 3.8284829441663635, 6.8421582395566762, 7.4385947937204353, 18.355423454952824, 18.939305693949098, 56.76654710662406, 46.472433795132453, 26.408180301616227, 8.5418678989122103, 3.4592633965371484, 0.81420697505187978, 2.8754034241015152, 3.3844724421254879, 2.3664842485474735, 1.572431813235436, 1.6803595914225982, 2.8859143014026083, 0.95828374076997469, 8.6117964749878801, 11.090643871879506, 6.7103974452754693, 5.7639040157862587, -8.5448676586843018, 23.490161756370185, 13.146852870439265, 4.4618664659962768, 4.6927135104951034, -1.0407935339809418, -1.4091416884493408, -1.9877285465952099, -0.10193290249882341, 1.7026114564114307, 4.5721973792850434, 2.0129650120803833, 4.2276867493724772, 8.0209890672234039, 34.346379481364394, 10.842808636228447, 26.9222291078285, 22.695665816424683, 10.916870623975989, -0.10515597578943024, 0.92746135521090112, 2.1750438821369982, 3.6410340078674537, 3.5866068225623455, 2.8117062685872063, 1.9464964791115584, -0.41942576489699973, 2.9442319665532848, 8.7117753152913178, 6.9832538798624677, -2.6420757676626758, 5.6480128145892872, 16.656214263229302, -6.389470676234752, 18.008772340009443, 8.6932746740913274, -0.52039830113473562, -0.97696074341019512, -2.9337382106235728, -2.1976406054664195, -0.80487566941320732, 1.657242840091119, 0.57695677140755675, -1.2392069417374179, -13.231297049531948, -4.7410363255113248, -2.7634934625894938, 0.83280248396082224, -2.5751982814309775, 4.566053007098211, 5.1937461344059415, 0.81982278416574916, -1.2648718051872136, 2.0294536549062832, 3.2286926050373275, 4.3044361492909653, 3.1399874906463054, 0.57687440104334498, 2.0512141006459621, 5.0295677911697396, 1.3893389983926545, -2.7028450063821756, 3.7140945878205622, 3.0810071966308161, 5.4454444116217395, 13.563457368659845, -2.5773278811340234, 9.4711485160420441, 5.1056392812224614, -4.0664292874937296, -2.1066158183405044, -2.0786620762717254, -0.39175809943363327, 0.38430207827145285, -2.3659355754891784, -10.190768615220893, -18.287776152890569, -34.94566909960389, -16.804391379899673, -42.447185502078064, -17.968893142627813, -9.910885111884598, -1.5845378867194357, -0.61503488296060749, 0.5456937943712441, 1.9706025084081269, 3.4768776627392723, 3.4105600445409681, 2.4699197782172639, 1.881917823375139, 2.4646463691490985, 1.2983452624583438, -1.2659388719629485, 0.1414776211221116, 3.2969221728544715, 1.3916869062248129, 2.1242155181960229, 2.6484997231173293, 12.171761179906973, -2.2274307006673681, 3.702478373345508, 0.81264916316704883, -3.8374461384037382, -1.0672361910520076, -0.60287176352493688, 0.73237790527531033, -2.1935496330750563, -9.3720360014434689, -23.364991274814056, -28.809040802191326, -49.47539998840697, -41.651580747175068, -46.656198806890998, -19.483874473280501, -5.8993661554872627, -0.85975587552355426, 0.36135071104124949, 2.5087646302129953, 2.6404974913828889, 2.7114466415734491, 2.7998549031056945, 2.9500669911593249, 0.62764371336605218, -1.8840461990319775, 0.76549003012129535, 1.789138288104325, -0.033223712968466615, 0.67271760119275215, 0.67618078211606658, 2.1778507792771671, 4.1537964250038266, 6.5139128293303603, -2.5566610369397242, -0.38826938489369667, -0.35663258536783948, -2.1525887062494129, 0.11701161368368215, 0.51276950901619689, -1.7682797425271899, -5.4849410264677143, -13.274828683787417, -26.796686775361572, -20.048589458196275, -78.507966752466459, -81.624687331128925, -59.038387406402052, -16.470770345397689, 1.4441506321142048, 1.429843803082246, 2.0278729258956907, 1.8248078988499317, 3.0361534846987319, 3.045953095599546, 1.8806785141996951, -0.036676325952591224, -0.64513432864685583, 0.72134917870519299, -0.33833230756510668, 0.35867475608511745, -1.7870448740379619, 0.41148866509475018, 2.385170800316434, 3.582324284715646, 2.5859312055320172, 2.5574024719687887, -3.08862682795997, -1.6215192060328063, 0.65441624739145521, -1.9581076589813555, 2.3416848160301531, -1.9389416594989646, -6.8324882374098541, -11.422688403833355, -24.503242066461048, -44.713956729761051, -297.71280467242298, -262.59106890839888, -153.12296297306594, -53.113548920223877, 1.2280777000160852, 4.3092019872504084, 1.8747628167747172, 3.341557420163666, 4.2845415972556511, 2.8130362926385151, 1.5357324193076258, 0.64495136714903156, 0.63947101827441044, -0.34981510145179595, -0.16542223594153913, -0.75319429727190634, -0.5564284933709801, -1.1733083260902135, 1.3896745329326132, 2.0729019173667376, 2.6435570864758655, 0.31649576481827513, 0.56033030111271809, -3.7001182860738244, -1.3772856826259623, 2.8686428937153425, -2.2247541222231662, 0.68623525685057685, -19.033005269759858, -40.587336797783074, -105.85720001872383, -206.25666421819477, -1015.8170093530997, -801.4985531526338, -447.3325626055032, -129.9505765405554, -7.094723129271765, 3.3531473980017035, 7.3823613929520073, 5.5799176422333412, 4.8975474263529168, 2.7230446482898696, 2.4101241863335168, 1.700641385933505, 0.025704615807253747, -0.99095562521376346, -1.4076418772333852, -0.37490601010063784, -1.5096317606573337, -0.36055346759863915, 0.12698658321626774, 2.1386535528573165, 2.3940543871676701, 0.59203001318833959, -1.3611346863617075, -0.7686322780239524, -4.9214217924247352, -2.0949194135756803, 12.245388279073195, -32.688631361987156, -66.670140953207678, -221.8093865245236, -464.04986763482748, -757.33856581210489, -2206.474238883316, -1877.1475049531243, -1129.7740740979648, -397.03998761961947, -72.365238046453996, 30.293009711209137, 26.738330274503873, 7.5675449146456106, 4.3849613992964835, 5.0282597338488833, 3.9892729615099154, 2.3612666445612378, 0.50013744114287839, -1.7835969801453351, -1.6822520344543537, -2.6713659240639043, -1.9174482807780469, -2.1579844977242066, 1.9017476504375719, 2.6260151255238871, 3.7232300368773963, 0.18037835028033014, -1.0558165909593593, -4.8653324939915077, 1.9503276148055297, -8.3653822807973572, -11.219037780374073, -25.222873313530549, -283.75431628355437, -607.97311048711254, -1427.0600295974182, -2207.4882969003352, -3720.099183631889, -2853.5518466904387, -2111.648140416873, -1058.7748743456016, -288.68438568639789, -2.5048396429994018, 41.761171107733141, 14.647883903858578, 9.6537837209937827, 8.0807144322732913, 5.5230847011321327, 3.6189445517483261, 0.3992847810370409, -0.12607857026309263, -2.1003135650570361, -0.79613102395168645, -2.2230696316747758, -1.2900498248808734, -0.0088018165554208834, 4.8667397989738888, 3.4357055351949497, 2.2081255979026997, -1.9843645852296314, -4.3910663300419372, -6.6803391366371141, 13.024576486880093, -42.171627369152944, -179.07785660663302, -402.86472873295685, -1347.5218927383694, -2223.332633868245, -3323.4413322476003, -3154.0647285180148, -4493.8282478995088, -3308.6059534728515, -1663.4471657442498, -572.34886359995289, -148.49350684902797, -6.1114787305125393, 24.681352644731046, 18.802141623354593, 8.8540252326545641, 6.2483313015722262, 4.0409079710738718, -0.47304719477879342, -1.361454509331699, -3.1026085263059842, -2.948515065656593, -0.5489892871976132, -2.0399346307109321, 2.7536893948952148, 4.8669904821248942, 7.1357881321682211, 2.7727430774895763, -1.706309766177841, -6.7152536525602784, -4.9992210979693477, -11.324269109891496, -15.123294918372496, -301.22169990910822, -836.26206181721659, -1226.2977859942089, -3085.2108289353782, -3175.2001764803231, -9007.3131796291, -7188.3398987071041, -5194.4090916202022, -2278.6820540988915, -810.15843536731495, -333.51324485695591, 20.969766857124615, 20.288236987688347, 23.758858737501669, 11.676825499792164, 10.210570606354567, 9.4957008241538698, 4.5937507427443718, -3.6748649027980855, -6.9403489082508045, -8.1315470965115626, -2.3288001542285697, -2.2335756869247856, -1.5210347396073984, 4.088293168231135, 8.0126520330591262, 6.6402880636320099, -0.51069538319843244, -5.7742639473864941, -8.2027965131620455, -3.3329385156007563, -97.440927522082632, -301.70227472838576, -654.09944760377232, -1292.5663502582681, -1439.1675666408923, -5346.3306090002598, -27232.421940549855, -20696.443600183229, -8881.5889351276946, -3113.2951228080024, -1428.8365567175017, -637.0667107212123, 89.303256695008088, -20.624592190099609, 25.446982376231595, -2.0330861649905447, 16.950197274587225, 14.857882580654767, 20.385879097388024, -7.8710994817946256, -1.9041291470930921, -18.745207674117459, -0.11803389887143546, -20.674418420848326, 10.030112101155776, -1.1889020916883362, 13.291351981828802, 1.9958373496879209, 6.472640746355359, -14.835941429972186, 15.63017417600328, -23.204353900286776, -49.325655652904786, -526.29225504526619, -42.22554400942483, -409.0173927639147, -6745.5682223604063, -15984.118662988625, -14789.032757371675, 2890.8320511623242, -8876.4278731043305, -3798.5470080468367, -816.05094481379706, -942.72760774203016, 34.134427334482538, -53.601071866949361, 18.961120931859348, -12.683126140049515, 15.4575031105396, 15.708115389403696, 29.148348103795954, -12.155673848712381, 0.82076815557010951, -11.664010623895932, -14.807142725465098, -16.625069210776097, -0.357968466821102, 1.860738389345522, 13.634662924808296, 6.8633911029604313, 4.1139030671065999, -5.4456740109950088, 4.6010715329523864, 33.024718339612249, -120.74864963582446, -488.33133225087175, -286.68893204831085, -866.78466658111245, -11739.029695619576, -31341.185848117751], "imag": [0.0, -38758.087655514835, 17548.085047269142, 6107.2151726624015, -1874.673387120311, 1102.4198144208276, 70.663498597697227, 246.50443360377355, 15.094404836632346, 61.578388193120794, -6.883850379905633, 17.433253862249003, -27.182352905748338, 10.205627484635695, -57.631165788617082, 3.5978635789949243, 0.0, -3.5978635789944389, 57.631165788617331, -10.205627484634785, 27.182352905748068, -17.433253862248129, 6.883850379906236, -61.578388193120603, -15.094404836632346, -246.50443360377048, -70.663498597698691, -1102.4198144208278, 1874.6733871203112, -6107.2151726624097, -17548.085047269135, 38758.087655514864, 21415.511526704078, 9240.4816801838842, 6140.0706338628515, 4606.7423967050836, -477.33409412636456, 106.10534760670754, 253.85827683649731, 118.670471872242, 60.852179067446961, 22.790809731222598, 16.992665274593527, 3.3051410229913314, -11.158748163655874, -19.125680216300069, -12.454719286949617, -13.739274709476872, 4.4649861577324996, -1.4030381833900449, 28.352674755287939, -0.67900523797834378, 10.299898240809684, 0.92198929608587377, 1.7353747463970652, -28.921613512340663, -62.077231813228479, -150.50951383829744, -214.24479740751684, -816.21574466320703, 274.1353961963593, -3254.3286688062722, -6340.0761451284025, 11789.651384311641, 10840.701965194425, 18380.248977426931, 3231.3822724351476, -461.44982475214084, -210.28695062763435, 162.42212696600188, 179.30887547689923, 113.80000983072934, 43.394725496932487, 24.59548477602943, 9.5064435172576207, 6.2249226907652995, -15.514190495563678, -15.67209469704693, -16.179758937655365, -0.44183444538403649, 3.28796177703477, 12.57005414423841, 4.7947040656811843, 5.5752197676029693, -1.566994278876134, 3.0671728952871886, -9.3006029966239439, -19.205045070959127, -56.939768493357462, -104.59324045979437, -209.92149032340299, -248.20433151041027, -867.92363746103479, -1299.0888976028032, 827.61688284819274, -3453.3010516440099, 2184.1299326736716, 343.81610498084558, 1430.9649107108667, -397.20604751762147, -695.96261474413006, 272.38382351444113, 181.43163787734565, 64.570588674311324, 29.183435066681643, 16.420203889130576, 10.051165980763317, 3.8132998567254699, -5.2200576234926608, -12.843645059741645, -16.574215370008311, -6.8948411718648668, 5.7853993354930253, 6.7615499928711884, 6.3258793363865982, 6.2779141420300757, 7.316455816584333, 0.11325917153592484, -12.472105947256997, -26.315486653091728, -48.315266721271882, -96.05176533673017, -133.18984159990899, -353.39952612993443, 150.90670251969735, -622.62775933130001, -1461.5017842629939, 2621.2839732861958, -797.42142339187535, -1141.9271710506443, -1665.2932495998289, -238.3366310615267, 49.923170262409073, 61.506465075545492, 153.25086116745172, 46.446591336114025, 13.171326963991907, 10.42525251508628, 5.9093841314446687, 3.7481817410606739, -5.009304162269796, -8.3377118534842136, -11.383176063308646, -7.2904568591879206, 0.055102663998529922, 6.2100452165049633, 0.50258781153530074, 8.7133201843540995, 8.0099491969619425, 3.408178127748358, -6.1241971022028396, -12.736128845291272, -45.232132533425755, -61.3791750629996, -153.42214558232837, -67.595097514570796, 85.625899103861244, -58.961155821254636, 56.883206971067992, -1105.7125877167775, -1066.0901100808514, -1171.2701320923616, -223.1510166903665, -49.090655504704266, 242.13940044346495, 187.57398243385657, 37.630624824958112, 28.849867481961542, 8.4569554634773265, 4.9702839577388582, 4.4028115702250554, 1.4310556478305501, -1.929060150422355, -6.8702061592831818, -7.5557834297932249, -6.2828152868831211, -0.10870593803434039, -0.35173698335235887, 2.5520560959229499, 3.0896458550054593, 4.7505220002799033, 1.2377659738692817, -1.7806530124299274, -9.1879671692559626, -20.164632449396379, -65.591240653702698, -60.869480400000917, -62.271759486118263, -69.04223895420877, 113.44009243151515, -280.24622279593513, -361.29505764988915, -32.2389038646619, 397.42763261597736, 412.9673372382922, 359.76401050336045, 192.6210800054142, 122.25882006598455, 49.109708531705785, -0.83416770813231356, 2.4803036198391495, 2.5045331182816137, 2.1157075816863125, 2.0456835555591151, -1.3319639541688315, -4.1057658132159602, -5.2622609905014679, -3.122429370817716, -0.23490025627283884, 0.40862649299968534, 1.0190926814192749, 2.6652215067301066, 2.3902937959311288, 1.3072901204360872, -2.1976906085986245, -4.7205289930149936, -17.964658686821839, -21.918559933832515, -40.854552179344822, -68.24351076369831, 25.856495589804037, -73.74721318359471, 26.189735785049756, -56.726426890704197, 461.88238037931336, 400.36625478870423, 367.42140066372616, 233.95717316464078, 120.24788693419019, 56.86229814651989, 16.528423808469729, 5.5081763037408891, -2.7285347510623343, 0.91768473485835012, 1.922818700657754, 0.85878662181529319, -0.082031723657331246, -3.1986596327501684, -2.7995613930163934, -2.624028681216124, 0.082235364558318216, 0.13384423623927519, 1.9947594371394528, 1.3437508199980643, 2.6398410865992954, 0.73758290379789837, -0.13078900915043976, -4.4953343292506398, -5.1950335909144894, -13.12242317098767, -25.523369493146305, -12.192196214097944, -25.868477512213353, 2.9066782920883432, 34.106574161015573, 256.00130492151595, 208.58597838816522, 310.3347945923577, 177.03891975992664, 97.122360992021001, 48.488510187113192, 24.528917872852549, 4.6272034575157122, -1.3571647909533122, 0.20621269482631011, -0.40867746076362715, 0.99215767258594889, 1.5313768256641691, -0.052772985639764861, -1.5552847934304939, -2.8318141898835965, -1.44106038161314, -1.7221472209588695, 0.61786747800524422, 0.48048213041449506, 1.7161555169971148, 0.88689757979195871, 1.4761658864245926, -1.7252645996654892, -0.81528783024945528, -3.3158091950668926, -7.6513645380807835, -2.9281859927886749, -12.462674546124116, -7.4696827367159715, -23.236903073402107, 43.7319836231637, 92.321882041990349, 70.117668379783737, 34.226649356314077, 47.258739472197036, 25.375797504306092, 12.406583490274283, 8.6535618432690811, 1.0747519793409035, -0.44324044078252567, -1.7584707522262704, 0.72458607165998101, 0.70692065864217946, 0.77076345191259166, 1.0649069250197341, -2.7261985428418192, -1.9302074459527443, -3.3754250277859401, 0.50673349306549698, -0.64031672309299725, 1.374850007673434, 0.72950552648490297, 2.3325303094340692, -0.9986473214357906, 0.78195787446934228, -1.1999635798104886, -3.4858876190620522, -1.4991310089769327, -2.8433755069243221, 3.7638952884372583, 1.2400980139632869, -2.9542792497780144, 1.8234355680479359, 49.576244746235467, 19.014215780516388, 30.877310974744532, -1.1221271518157558, 7.8208026967228905, 2.609067592037944, 1.2521971880796572, -0.99505064345958283, -1.3975466814991453, -0.57338407820284554, -1.8697004151394911, 2.3700280505641378, 1.696784850403636, -0.61519740478736618, -1.3317994275810727, -4.2071746247346917, -0.85834734780410316, -1.4672332071937433, 2.6381100233578914, -2.4176769851497171, 2.0131074000190821, -1.6262628596018911, 3.9814331515508763, -0.027677673297514374, -2.444607041276003, -1.6150024285043403, -4.1227589144067673, -1.7068861549694658, 2.2919486155920237, 12.185386506158462, 9.103962312567619, 7.1597604706792799, -5.4087151249860161, 0.96974162117579055, -11.947201833901593, 6.6330620252490986, 2.1152085218133552, 3.546929348907518, 0.62283616740311709, -1.3522725353705505, -2.4444202103018355, -1.5690941230144015, -0.0069858778603680161, -2.2282653754616524, 3.8284308308873771, 3.1912913728027457, -5.5805618725494135, 0.13049495818795365, -3.5197176720292105, 1.8802124921714005, -2.6189977871286096, 2.1235611153433935, -7.1463089769422359, 4.7729314123471731, 2.1077672614040215, -0.81910125894631813, 0.72052624692713507, -3.4842402749887631, -4.1829115306886377, -4.3994325344570067, 3.0953565929079909, 11.949390330752637, 23.680136621529911, 21.558658775121675, 22.744598059672118, 24.943328767928129, 28.901352722188449, 12.814680673201753, 13.644932734111247, 5.45768556544121, 2.4025753394703617, -1.1673355625318591, -2.4971924284976894, -1.8952846661475429, -1.2927419412763479, 2.7902732019969552, -4.3477010834394934, 3.7423253621601038, 4.3959804108055893, -8.9257994152947724, 2.882370251254919, -5.0577698714645729, 3.1829717367300585, -13.324209534416688, 6.2785045354550411, -0.87489834284363188, -0.85263104321571348, 4.4885662153073236, -0.53110232172652128, -1.6389132121610095, -4.3945543864166812, -4.2936344828625694, 0.48745205753098192, 15.895360345460416, 25.169617383292081, 33.918445473292621, 23.996907911576923, 23.251277269359601, 5.1543156647143702, 18.554102826608389, 16.359265688234956, 10.231244633613841, 1.5587865083621433, -1.4683911434431964, -2.958105786087343, -1.5689869488245065, -0.96051225634847215, 1.7578951205789273, 6.6638352938822507, -9.4901775499106478, 4.1317473271038834, 10.951455398944859, -15.360897164532773, 5.739294741971535, -17.045217645164048, 5.8449684444103616, -3.5521886878464675, -7.3858280585812945, 5.7347340895196259, 0.4908832539545151, 0.98666637905349674, -1.0149409800154179, -1.6832477381822215, -2.5467748158366237, 2.0904065823844138, 13.443724261124899, 31.761133485201317, 33.671186166611122, 40.707191608760652, -1.1509826127936451, 21.096779572504321, 10.924705320167561, 7.3688776104039446, 6.4733533392666498, 1.1232923919614803, -2.4431302941492752, -2.7804484266994756, -1.6225398299131664, 0.6478030911844409, 2.5793951871330139, 3.1124056384095988, 9.3832743597685244, -14.195954313031031, -0.57608988546761064, 19.479390479219269, -31.200338855208866, 14.100113777774123, -5.3856662644378392, -10.946378631138824, 4.4142064356850845, -3.9066652322250062, 0.66844910598604979, 1.0227755630951001, 0.66137538568075327, -0.38567118503773246, 0.24821124990411303, 3.5872841231435211, 14.053276636456639, 24.260931372233497, 33.647459095648763, 7.6529271231588307, -9.2151377078706354, -38.311665858686638, -16.121399661234683, -5.9927868178298498, -2.7701365304211771, -3.7140088945870469, -3.0788396721523017, -3.0904878618769733, -1.8655827999533823, 0.43690140001808947, 4.1561194908680763, 4.3858952392085353, 2.7714027971729753, 9.257895969243517, -19.629364588289473, -14.072043519544888, 36.373196183604975, -13.61372365841147, -9.4583423311161248, 8.5284583277520429, -6.1744651457950592, -2.2310052363098367, -1.75629520897435, 1.1539693350814639, 1.7596902661237352, 1.8196834579374304, 2.022290334087173, 5.1646829084835781, 12.770276815017757, 24.46979490512981, 27.822785449092486, 29.828540342091252, 0.0, -3.7653763488297254, -24.757305904406468, -19.738984856989376, -8.2664619761594462, -4.6495144591838038, -3.656340919703911, -2.9147641885479052, -2.2830188616483764, -0.10916825448095578, 2.3916977668177819, 5.941240494184604, 2.1224657474819093, 3.9401472638367774, -2.1445277323266949, -10.774805376941604, 0.0, 10.774805376941531, 2.1445277323266829, -3.9401472638367951, -2.1224657474819155, -5.9412404941846173, -2.3916977668177926, 0.10916825448094437, 2.2830188616483764, 2.9147641885479065, 3.6563409197039234, 4.6495144591838136, 8.2664619761594498, 19.738984856989411, 24.75730590440644, 3.7653763488297773, 9.2151377078736516, -29.828540342087546, -27.822785449090652, -24.469794905129369, -12.770276815017642, -5.1646829084834636, -2.022290334087153, -1.8196834579374237, -1.7596902661237566, -1.1539693350814775, 1.7562952089743382, 2.2310052363098212, 6.1744651457950415, -8.5284583277521193, 9.4583423311161194, 13.613723658411436, -36.373196183604932, 14.07204351954484, 19.629364588289484, -9.2578959692435348, -2.7714027971730033, -4.3858952392085584, -4.1561194908680932, -0.43690140001811473, 1.8655827999533605, 3.0904878618769729, 3.0788396721523177, 3.7140088945872396, 2.7701365304214618, 5.9927868178306163, 16.121399661236502, 38.31166585868818, 1.1509826127964509, -7.6529271231568821, -33.647459095647825, -24.260931372233326, -14.053276636456571, -3.5872841231435073, -0.24821124990410889, 0.3856711850377551, -0.66137538568076071, -1.0227755630951023, -0.66844910598607266, 3.9066652322249751, -4.4142064356851165, 10.946378631138813, 5.3856662644378339, -14.100113777774116, 31.200338855208855, -19.479390479219351, 0.57608988546757856, 14.195954313030994, -9.3832743597685617, -3.1124056384095939, -2.5793951871330338, -0.64780309118445645, 1.6225398299131557, 2.7804484266994689, 2.4431302941492916, -1.1232923919613498, -6.4733533392663958, -7.368877610403441, -10.924705320166419, -21.096779572501951, -23.251277269357679, -40.707191608758833, -33.671186166610063, -31.761133485200901, -13.443724261124748, -2.0904065823843769, 2.5467748158366237, 1.6832477381822348, 1.0149409800154132, -0.98666637905348697, -0.49088325395452825, -5.734734089519665, 7.3858280585813008, 3.5521886878464564, -5.8449684444103243, 17.045217645163998, -5.7392947419715155, 15.360897164532778, -10.951455398944907, -4.1317473271038754, 9.4901775499106655, -6.6638352938822702, -1.7578951205789377, 0.96051225634844795, 1.5689869488244914, 2.9581057860873381, 1.4683911434432237, -1.5587865083620118, -10.231244633613667, -16.359265688234554, -18.554102826607274, -5.1543156647129926, -24.943328767927714, -23.996907911576475, -33.918445473292081, -25.169617383291836, -15.895360345460304, -0.48745205753098297, 4.2936344828625552, 4.3945543864167149, 1.6389132121610113, 0.53110232172653837, -4.4885662153073405, 0.85263104321570904, 0.87489834284363122, -6.2785045354550046, 13.324209534416719, -3.1829717367300314, 5.0577698714645702, -2.8823702512549096, 8.9257994152947653, -4.3959804108056124, -3.742325362160086, 4.3477010834395236, -2.7902732019969543, 1.2927419412763401, 1.8952846661475349, 2.4971924284976934, 1.1673355625318764, -2.4025753394703346, -5.4576855654411656, -13.644932734111133, -12.814680673201632, -28.901352722188015, -0.96974162117503315, -22.744598059671098, -21.558658775121231, -23.680136621529623, -11.949390330752514, -3.0953565929079563, 4.3994325344570004, 4.1829115306886644, 3.4842402749887649, -0.72052624692711631, 0.81910125894632069, -2.1077672614040264, -4.772931412347166, 7.1463089769422288, -2.1235611153433602, 2.6189977871286789, -1.8802124921713912, 3.5197176720292407, -0.1304949581879572, 5.5805618725494606, -3.191291372802759, -3.8284308308873531, 2.228265375461667, 0.0069858778603283057, 1.5690941230143916, 2.4444202103018267, 1.3522725353705702, -0.62283616740308378, -3.546929348907534, -2.1152085218133783, -6.6330620252487682, 11.947201833901675, -19.014215780515311, 5.4087151249861485, -7.159760470678993, -9.1039623125675071, -12.18538650615837, -2.2919486155920032, 1.7068861549694805, 4.1227589144067789, 1.6150024285043441, 2.4446070412760443, 0.027677673297527423, -3.981433151550863, 1.6262628596018724, -2.0131074000190696, 2.4176769851497326, -2.6381100233578612, 1.4672332071937413, 0.85834734780414157, 4.2071746247347095, 1.3317994275810989, 0.61519740478738694, -1.6967848504036496, -2.3700280505641422, 1.8697004151394931, 0.5733840782028331, 1.397546681499154, 0.99505064345958893, -1.2521971880796428, -2.6090675920379729, -7.8208026967229012, 1.1221271518159068, -30.877310974743335, -70.117668379782501, -49.576244746235311, -1.8234355680484657, 2.9542792497778261, -1.2400980139631994, -3.7638952884372125, 2.8433755069243154, 1.499131008976947, 3.4858876190620642, 1.1999635798105424, -0.78195787446931964, 0.9986473214357916, -2.332530309434083, -0.72950552648488864, -1.3748500076734225, 0.64031672309302556, -0.50673349306550342, 3.375425027785985, 1.930207445952772, 2.7261985428418698, -1.0649069250197167, -0.77076345191255435, -0.7069206586422051, -0.72458607165999755, 1.7584707522262737, 0.44324044078250502, -1.0747519793408324, -8.6535618432689372, -12.406583490274217, -25.375797504305602, -47.258739472196019, -34.226649356312784, -208.58597838816522, -92.321882041990833, -43.731983623163735, 23.236903073401958, 7.4696827367159822, 12.462674546124104, 2.9281859927886771, 7.6513645380808271, 3.3158091950668926, 0.8152878302494978, 1.7252645996655027, -1.4761658864245817, -0.88689757979195827, -1.7161555169970975, -0.48048213041447491, -0.61786747800518715, 1.7221472209588695, 1.4410603816132017, 2.8318141898836315, 1.5552847934305372, 0.05277298563978073, -1.5313768256641445, -0.9921576725859278, 0.40867746076360606, -0.20621269482631011, 1.3571647909533708, -4.6272034575156908, -24.528917872852457, -48.488510187113178, -97.122360992020845, -177.03891975992676, -310.33479459235781, -461.88238037931535, -256.00130492151709, -34.106574161017335, -2.9066782920893011, 25.868477512213314, 12.192196214097933, 25.523369493146241, 13.12242317098773, 5.1950335909144911, 4.4953343292507366, 0.13078900915048275, -0.73758290379785463, -2.6398410865992936, -1.3437508199980366, -1.9947594371394133, -0.13384423623915356, -0.082235364558330679, 2.6240286812162581, 2.7995613930164374, 3.1986596327502643, 0.082031723657369299, -0.85878662181523746, -1.9228187006577324, -0.91768473485829727, 2.7285347510623108, -5.5081763037407097, -16.528423808469462, -56.862298146519464, -120.24788693418984, -233.9571731646403, -367.42140066372525, -400.36625478870656, 32.238903864660138, 56.726426890699642, -26.189735785051667, 73.747213183594084, -25.856495589804108, 68.243510763698637, 40.854552179345092, 21.918559933832821, 17.964658686821831, 4.7205289930151793, 2.1976906085987129, -1.3072901204360032, -2.3902937959311101, -2.6652215067300165, -1.0190926814191563, -0.40862649299943438, 0.23490025627283342, 3.1224293708179411, 5.2622609905015993, 4.1057658132160988, 1.331963954168899, -2.0456835555590454, -2.1157075816862583, -2.504533118281568, -2.4803036198391473, 0.83416770813248287, -49.109708531705543, -122.25882006598383, -192.62108000541431, -359.76401050336045, -412.96733723829362, -397.42763261597815, 1066.0901100808455, 361.29505764988636, 280.2462227959299, -113.44009243151743, 69.042238954208585, 62.271759486118853, 60.869480400001216, 65.591240653703366, 20.164632449396354, 9.1879671692563107, 1.7806530124300617, -1.2377659738691233, -4.7505220002798652, -3.0896458550052697, -2.5520560959227772, 0.35173698335277165, 0.10870593803433028, 6.2828152868834488, 7.5557834297934026, 6.8702061592834163, 1.9290601504224327, -1.4310556478304128, -4.402811570224995, -4.9702839577387659, -8.4569554634772928, -28.849867481960956, -37.630624824957771, -187.57398243385583, -242.13940044346546, 49.090655504701473, 223.15101669036463, 1171.270132092355, 797.42142339187455, 1105.712587716775, -56.883206971068894, 58.961155821254195, -85.625899103861229, 67.595097514571449, 153.42214558232899, 61.379175063000453, 45.232132533425769, 12.736128845291734, 6.1241971022030093, -3.4081781277481404, -8.0099491969618857, -8.7133201843538846, -0.50258781153503296, -6.2100452165045663, -0.055102663998530817, 7.2904568591883807, 11.383176063308918, 8.3377118534845209, 5.009304162269923, -3.7481817410605087, -5.9093841314446331, -10.425252515086191, -13.171326963991898, -46.446591336113507, -153.25086116745101, -61.506465075545563, -49.92317026240935, 238.33663106152528, 1665.2932495998259, 1141.9271710506428, -2184.1299326736644, -2621.2839732861944, 1461.5017842629866, 622.62775933129569, -150.90670251969709, 353.39952612993454, 133.18984159990904, 96.051765336731378, 48.315266721271847, 26.315486653092378, 12.47210594725734, -0.11325917153558851, -7.3164558165842362, -6.2779141420297728, -6.325879336386472, -6.7615499928705614, -5.7853993354930182, 6.8948411718656626, 16.574215370008741, 12.843645059742064, 5.220057623492746, -3.813299856725167, -10.051165980763182, -16.420203889130285, -29.183435066681557, -64.570588674310201, -181.43163787734505, -272.38382351444022, 695.96261474412938, 397.20604751762363, -1430.9649107108592, -343.81610498084137, -10840.70196519443, 3453.3010516439876, -827.61688284819604, 1299.088897602802, 867.92363746103354, 248.20433151041067, 209.92149032340356, 104.5932404597956, 56.939768493357427, 19.205045070960029, 9.3006029966243737, -3.0671728952866175, 1.5669942788764133, -5.5752197676024959, -4.7947040656806426, -12.570054144236753, -3.2879617770347638, 0.44183444538609318, 16.179758937656128, 15.672094697047777, 15.514190495563801, -6.2249226907648021, -9.506443517257205, -24.595484776028343, -43.394725496932452, -113.80000983072682, -179.30887547689812, -162.42212696599955, 210.28695062763467, 461.44982475214965, -3231.3822724351485, -18380.248977426934, -21415.511526704078, -11789.651384311637, 6340.0761451283952, 3254.3286688062612, -274.13539619635947, 816.21574466320624, 214.24479740751747, 150.50951383829917, 62.077231813228565, 28.921613512340514, -1.7353747463964819, -0.92198929608524804, -10.299898240809465, 0.679005237978899, -28.352674755287008, 1.4030381833925814, -4.4649861577325289, 13.739274709477487, 12.454719286950338, 19.12568021630095, 11.158748163655963, -3.3051410229907021, -16.992665274592852, -22.790809731221124, -60.852179067446897, -118.67047187223945, -253.85827683649717, -106.10534760670436, 477.33409412636507, -4606.7423967050645, -6140.0706338628543, -9240.4816801839206]}};

var nose_filter = {"real": [3.0408379415611857, 0.37822261363137938, 1.1747173276627942, 0.87791572866957501, 0.27429578524536991, -0.26015038230887205, -0.04589694540462394, 0.18194913988848616, -0.19035992880648842, -0.050064472206627651, -0.016354960260506177, 0.073013903204916478, 0.15053197342937383, -0.47378158951033317, -0.33333577915227425, -0.30434085338299055, 0.26372284406307878, -0.30434085338299449, -0.33333577915227836, -0.47378158951033028, 0.15053197342937505, 0.073013903204915576, -0.016354960260505706, -0.050064472206628768, -0.19035992880648842, 0.18194913988848541, -0.045896945404624932, -0.2601503823088735, 0.27429578524536952, 0.87791572866957457, 1.1747173276627942, 0.3782226136313796, 1.5638480130127614, -1.0502161590016221, -0.28144140655806793, 0.66007876209063499, 0.25575295594566605, -0.11471557758718753, 0.22751157962613866, -0.29100003355994186, 0.29242994865853178, 0.14950061311078749, -0.02214037256893809, 0.034503572794411914, -0.17412889405901469, 0.060297128338254534, -0.49717312826403576, -0.10585857030501357, 0.082780177430083327, -0.11056584899062352, -0.72499416196392785, -0.15612863751308359, -0.36586758725441215, -0.19718551141955715, -0.043434328775718851, -0.060442748560529294, 0.1383233017601592, -0.071714399164362744, 0.31725814504248484, 0.16801904131477274, 0.58675004339129755, 1.5155025161072819, 0.99232752650122635, -1.0129349058609656, -1.6725083187014358, -0.66077224091757791, 2.2217962254133496, 0.98390492917013095, -0.16319024909703178, 0.67271651276019184, 0.16380771710007475, 0.21989373735679621, -0.12296413131763459, -0.18224776807046256, -0.011168077595514938, 0.10364961453623285, 0.1359801303077682, -0.26339593722715721, 0.039396952561874886, -0.33866696253626483, -0.0043101841409601402, 0.17740469505625342, -0.10830858623812685, -0.33186051211466239, 0.12060535717228799, 0.084703474730635803, -0.086262426064740078, -0.1297477439897661, -0.28426475792188655, 0.21965961682575447, 0.24955999100642301, 0.46804088594167448, -0.46946217188769823, -0.31724338121882617, 1.7716903674840139, 0.56577881124465312, 0.71512595047089234, 0.4910867492003736, 0.92310119575914651, -1.0939934957719639, 0.18102793706565934, 0.35484595526186652, 0.48676804973393767, -0.097923922056211768, -7.8200063684236368e-05, -0.016095705887683227, -0.084263209727683391, 0.0080331646697531477, -0.044853168166469837, 0.0062170702455977547, -0.15890945005602863, 0.22815908460498888, 0.17110516400703918, -0.29052216822588434, -0.18224178792918044, -0.018421573557969844, 0.032535239746583795, -0.08121081614537759, -0.029175915115094032, 0.092929862040843081, 0.16824405476099044, -0.24679227625704581, 0.00030176170486908785, 0.023943010115746295, -0.43281940369171445, -0.57166836623813821, 1.2482467538506352, 0.74385331371246666, 0.54787745216002559, 0.92898941300038174, 0.54473301805335128, -0.4877077226718825, -0.40370312470951597, 0.39096527623299848, 0.12644205243911047, 0.0424168022733412, 0.049219693586430309, -0.077748186566067037, -0.043614575991677208, -0.0041635061718157101, 0.07897222152752259, 0.053944973270972973, 0.11232623315704468, -0.070445625631503542, -0.072057201177642652, 0.12706805404252855, -0.13757198545862262, -0.31639223989443899, -0.029240486244660611, 0.10973045554309782, 0.04102655016489614, 0.01891084247324943, -0.056922369420476487, 0.030805208457151782, -0.1512000567754726, -0.058122663069861934, -0.30506910575481916, -0.79033868813295383, 0.026308878650124864, 1.7521122327973906, 0.64566220514793682, 0.40042515638418263, -0.33461208848408402, -0.55369429571770823, -0.030140868003777126, 0.27230784814155334, 0.008194429787378334, -0.030645866021173498, 0.023928227792328404, 0.090706510487419734, -0.021569631345745226, -0.0064127815367693367, -0.026343653609827326, 0.0079365608693974798, -0.18791450411345401, 0.14319521757460194, 0.17003550871316017, -0.24661106106331201, -0.21167593067635404, -0.098206750946426211, -0.048592759447365123, 0.0028129255950110439, 0.049868939755416952, 0.039023237394102607, -0.0036267311717218496, -0.12930461632538415, -0.067411320951837855, 0.010801205387038189, 0.074931291451028059, 0.0657071501308885, -0.22930577035731903, -0.34894672904561169, 0.37107135970432037, 0.49136473213042475, 0.10645611255411164, -0.27741376367012294, -0.015474080963115254, 0.13417678111594294, -0.00046677178847574514, -0.064152245495751381, 0.038292630722847554, 0.017696001575747813, 0.0052157631901947011, -0.023883294001795573, 0.00060518848897931612, -0.088395093514235953, -0.11757758027246699, -0.019555697374511325, 0.025650630137194243, -0.043746931993730094, -0.08846438922222799, -0.074442233736365312, -0.0099104052825428144, 0.064871673248634965, 0.03524540038796406, -0.0065163090427015676, -0.065286881502250399, 0.023130909787941502, -0.028452880392855154, -0.098638256951482817, -0.23382813425811513, -0.077933249962103476, 0.069575788839216612, 0.86635491100693507, 0.33668792873864722, 0.34253529215971401, -0.019122146869457048, -0.016521364547457899, 0.098161425946316522, -0.011130337821465129, -0.092543539386829035, -0.0046973120649884213, -0.0086885172531639618, 0.041915865623990751, 0.004043623749362439, -0.037907101781325205, -0.03902357631677511, -0.037690963997300889, -0.11419210623041842, -0.076732726183972602, 0.12573571356763677, 0.012633232837618465, -0.098942005293094951, -0.061888601165125671, -0.0090012900880717846, -0.014247675993027883, 0.027023228844142398, -0.0050282774165556453, -0.031667371615742887, -0.033818612421867206, -0.01148054478950835, 0.002410610826462259, 0.084402385062259444, 0.069824545029654947, 0.043930202648898894, 0.078872059827297045, 0.055053359776043469, 0.0133618922103162, 0.068453834537836655, 0.081701490128889276, -0.023034972903402569, -0.042406682273060693, 0.035148729718712253, -0.037569795749083354, -0.0041828274568459145, 0.012961699649937737, -0.015475583648045266, -0.027857672235467607, 0.011169867646763214, -0.052725469767566695, -0.0070148730067698924, 0.053188940442296044, -0.095289423698368528, 0.016377276604450954, -0.015232023889946105, -0.077682201861009353, -0.048207693014304241, 0.014396940553395536, 0.021075876818278678, -0.0024069559974744716, -0.016758259649352489, 0.031865879706987753, 0.032824623168326933, -0.074182732572421473, -0.080036865363297197, 0.055652716802045685, 0.11130862970905192, 0.072952972838101024, 0.035098396968263772, 0.070242479721294437, 0.040747896169172335, -0.022459274904977577, 0.014257984954846613, 0.019887732012022184, -0.040811004345713692, 0.0088191619820225595, 0.012579825412037176, 0.0010521231241174983, 0.006430724322423028, -0.0107697788586532, -0.065313535534231901, -0.062769573271966861, 0.11452944472315789, 0.040092306883550234, 0.12236977033818973, 0.13624403896529741, 0.027359429809611827, -0.0259369720098807, -0.021419082531908291, -0.02797214703330133, 0.0080320798501609496, 0.0013100558771969473, -0.0052561866064455981, -0.0097375289490629017, 0.028455863356249757, 0.056425701741664613, -0.0014166564531825356, -0.032528933027316242, -0.019619112630340951, 0.065700689463789611, -0.027574493612476551, -0.032810042710946841, -0.021337213333243964, 0.019592990080295629, -0.019759750092557201, -0.030135000999364853, 0.031698488621630008, 0.023522054710864718, -0.021766544372719836, -0.0029827434432790761, 0.0078373508878122541, -0.0028057606437618767, -0.033232550782043478, -0.070830818977105636, 0.059080455669844649, 0.36951320092865036, 0.091411128665944044, 0.16404233417264649, 0.062598270207454218, 0.044986403519809939, -0.081510896741027974, -0.042341618007147169, -0.019737461729971092, -0.0048974001169366831, 0.01654258945512711, 0.0088012977267579451, -0.013696041496545335, -0.0020268100930256635, 0.015845305908809133, 0.016458612619479462, 0.00059474001274681517, -0.078529066256991353, -0.013402803779002156, -0.0036013299995064228, 0.037573591582025043, -0.0076636494655425947, -0.0021144243161753088, 0.032140257903604356, 0.0203960516207237, -0.013294834505552555, 0.0021732333351695376, 0.015223562704661827, -0.01560949450401263, 0.015060609589555461, 0.099241248305729085, -0.19173261262163863, 0.2450151004066432, 0.26383919410790774, 0.33963060937430933, 0.34586223870135113, 0.1799263961717959, 0.13137374163691276, -0.036493858716968031, -0.041981936023927684, -0.060258932305754208, -0.0099079010405710068, -0.010024579181759819, 0.016637976300589805, 0.031629106230310143, -0.011441752446030828, -0.010113005489671789, -0.027490450999295705, -0.018258671855796018, 0.061296263921274532, 0.10180457538537498, 0.097647442588889397, 0.063246459754600901, 0.0018146068702083835, -0.012539173975430248, -0.008345642545901014, 0.020763503842579332, 0.035996162307667869, -0.015675575409266479, -0.029662747125522582, 0.021312770128680401, 0.08274253615935892, -0.10682054954404013, -0.033729052458342255, -0.12662108035586636, 0.18053751384272299, 0.70940599530087467, -0.21337246968950196, 0.09763637724258438, 0.15345959634626791, -0.11690110107335749, 0.030098075283007168, 0.034186076429122053, -0.021402902917041923, 0.0046349096026103891, 0.011776544425663813, -0.0089673934626629329, 0.025701507905379287, -0.029333520957471802, 0.0051215906448292756, -9.6207917787476249e-05, -0.027988212781709806, 0.093695296747695708, 0.10174972110658871, 0.1657472048824174, 0.097917377942654052, 0.0075790364114050627, 0.031051039914505564, 0.015692601846482499, -0.015197940538429407, 0.026154427883432859, 0.038178927860434288, 0.07391023019778703, 0.0062621769249755375, 0.14786394857966245, -0.24133125024592791, 0.31944609486771169, 0.50715637065140251, 0.3628700351691248, 1.3060302403517683, -0.51467834215318464, -0.23826829822593068, -0.0093468845104923252, 0.16015106243072727, 0.02805086056222977, -0.0035899004350242064, -0.0055435803643170143, -0.026732920924661639, -0.010772744125100176, -0.024390993371714485, 0.012627431507206405, -0.0087203787935273471, 0.017390054254405412, 0.12375300108425245, 0.24623032442851775, 0.32220608581938742, 0.14931322115435106, -0.030460558403777184, 0.017160725214558258, -0.0094839753900354314, -0.040937983123008212, 0.021461267056132466, 0.052605242734811673, 0.082583919156513017, 0.060994911510484674, 0.19212485777257118, 0.030490864773295545, 0.17592824054988929, 0.2988664581273196, -0.3440942993173739, 1.3210609586317499, 0.059297195869563314, 0.66450815180221101, -0.18837022473974091, 0.34310308052214539, 0.41335445680081984, 0.18880277280085617, -0.0025387166530481263, -0.027056960212979365, -0.014992917126512385, -0.065825070508153574, 0.0059924400617866053, 0.016058491686644451, 0.096318218677932735, 0.11648035795043635, 0.083105673180647682, 0.32671351148908129, 0.21720864831216413, 0.29412476851730657, 0.14853499339159934, 0.018970981025801619, -0.0063191769526926853, -0.01818046447321878, -0.033255690840711123, 0.042088106914846247, 0.11598544099090941, 0.13859711550376397, -0.01128947190338442, 0.086653868373652729, 0.5559656926196751, 0.37680556076813249, 0.88488899182435787, -1.1019114795534606, 1.0072326177374842, -0.25696609936231163, 0.22595605005106759, 0.20124367647634336, 0.25020565573121023, 0.15194580903531119, 0.054188690404279045, -0.031588472981034531, -0.072469006265441768, -0.025752917138555535, -0.0471489088024604, 0.079293556902861981, 0.16429646904639772, 0.19502382977708654, 0.21468548823411557, 0.26510515339808566, 0.4121092848683291, 0.24137856300787749, 0.039245063706659693, 0.01328784964664576, 0.063938494959291103, -0.056808918829038463, -0.039159725623022792, 0.028736615885883075, 0.075674855852819814, 0.10715382591497533, 0.29941075202959372, 0.24751816611443234, -0.14971685838339424, 0.92725020353048182, -0.58030730777587713, 1.4965419400087898, -0.58030730777587824, 0.92725020353048315, -0.1497168583833956, 0.24751816611443228, 0.29941075202959394, 0.10715382591497527, 0.075674855852819592, 0.028736615885883075, -0.039159725623022994, -0.05680891882903847, 0.063938494959291325, 0.013287849646645755, 0.039245063706659422, 0.24137856300787733, 0.41210928486832976, 0.32671351148908501, 0.21468548823411696, 0.19502382977708541, 0.16429646904639761, 0.079293556902861592, -0.047148908802460261, -0.025752917138555938, -0.072469006265442351, -0.031588472981034621, 0.054188690404279302, 0.15194580903531177, 0.25020565573121106, 0.20124367647634361, 0.22595605005106811, -0.25696609936231096, 1.007232617737484, -1.1019114795534584, 0.88488899182435599, 0.37680556076813088, 0.55596569261967799, 0.086653868373653353, -0.011289471903383815, 0.13859711550376486, 0.11598544099090989, 0.042088106914846005, -0.033255690840710936, -0.018180464473218153, -0.0063191769526928735, 0.018970981025801085, 0.14853499339159881, 0.2941247685173069, 0.21720864831216452, 0.24623032442851894, 0.083105673180647877, 0.11648035795043606, 0.096318218677932249, 0.016058491686644111, 0.0059924400617863598, -0.065825070508153657, -0.014992917126512623, -0.027056960212979545, -0.0025387166530475712, 0.18880277280085703, 0.41335445680082034, 0.34310308052214539, -0.18837022473973905, 0.66450815180220979, 0.059297195869561843, 1.3210609586317497, -0.34409429931737384, 0.2988664581273206, 0.17592824054988909, 0.030490864773296485, 0.19212485777257224, 0.06099491151048541, 0.082583919156513072, 0.052605242734811722, 0.021461267056132546, -0.040937983123008163, -0.0094839753900352007, 0.017160725214558203, -0.030460558403777944, 0.14931322115435064, 0.32220608581938709, 0.093695296747697082, 0.1237530010842532, 0.017390054254405134, -0.0087203787935278762, 0.012627431507206037, -0.024390993371714388, -0.010772744125100247, -0.026732920924662094, -0.0055435803643169856, -0.0035899004350241582, 0.028050860562230339, 0.16015106243072863, -0.0093468845104915464, -0.23826829822592965, -0.51467834215318398, 1.3060302403517683, 0.3628700351691248, 0.5071563706514024, 0.31944609486771208, -0.24133125024592691, 0.14786394857966259, 0.0062621769249767058, 0.073910230197787502, 0.038178927860434683, 0.026154427883432977, -0.015197940538429256, 0.015692601846482655, 0.031051039914505491, 0.0075790364114048233, 0.097917377942654316, 0.16574720488241704, 0.10174972110658852, 0.10180457538537414, -0.027988212781709362, -9.6207917787273747e-05, 0.0051215906448293215, -0.029333520957471861, 0.025701507905379183, -0.0089673934626627958, 0.01177654442566399, 0.0046349096026104134, -0.021402902917041607, 0.03418607642912222, 0.030098075283007585, -0.1169011010733573, 0.15345959634626896, 0.097636377242582312, -0.21337246968950468, 0.70940599530087356, 0.18053751384272154, -0.12662108035586525, -0.033729052458342325, -0.1068205495440398, 0.082742536159359198, 0.021312770128680637, -0.02966274712552271, -0.015675575409266458, 0.03599616230766798, 0.020763503842579301, -0.0083456425459010053, -0.012539173975430275, 0.001814606870208027, 0.063246459754600651, 0.097647442588889521, -0.013402803779000846, 0.061296263921273318, -0.018258671855796015, -0.027490450999295389, -0.010113005489671634, -0.011441752446030721, 0.031629106230310081, 0.016637976300589583, -0.010024579181759814, -0.0099079010405707553, -0.060258932305754062, -0.041981936023927795, -0.036493858716967892, 0.13137374163691345, 0.17992639617179582, 0.34586223870135102, 0.3396306093743105, 0.26383919410790552, 0.24501510040664315, -0.19173261262163757, 0.099241248305728891, 0.015060609589555881, -0.015609494504012908, 0.015223562704661735, 0.0021732333351695459, -0.013294834505552418, 0.02039605162072387, 0.032140257903604141, -0.0021144243161752658, -0.0076636494655425496, 0.037573591582025015, -0.0036013299995066783, -0.027574493612477238, -0.078529066256990701, 0.00059474001274650921, 0.016458612619479316, 0.015845305908809258, -0.0020268100930258088, -0.013696041496545191, 0.0088012977267581446, 0.016542589455127047, -0.0048974001169365296, -0.019737461729970797, -0.042341618007146614, -0.081510896741027863, 0.044986403519810654, 0.062598270207453191, 0.16404233417264671, 0.091411128665944058, 0.36951320092864953, 0.059080455669845232, -0.070830818977105567, -0.033232550782043346, -0.0028057606437617323, 0.0078373508878121518, -0.0029827434432791976, -0.021766544372719729, 0.023522054710864829, 0.031698488621629953, -0.030135000999364715, -0.019759750092557201, 0.019592990080295483, -0.021337213333244086, -0.032810042710946043, 0.035098396968265562, 0.065700689463788264, -0.019619112630341284, -0.032528933027315937, -0.001416656453182401, 0.056425701741664544, 0.028455863356249663, -0.0097375289490629503, -0.0052561866064455114, 0.0013100558771970501, 0.0080320798501610554, -0.027972147033301621, -0.021419082531907708, -0.025936972009880058, 0.027359429809612409, 0.13624403896529816, 0.12236977033819048, 0.040092306883550463, 0.11452944472315811, -0.062769573271967, -0.065313535534231998, -0.010769778858653123, 0.0064307243224230766, 0.0010521231241173155, 0.012579825412037126, 0.0088191619820226965, -0.040811004345713484, 0.019887732012021993, 0.014257984954846677, -0.022459274904977573, 0.040747896169171828, 0.070242479721293896, 0.055053359776043469, 0.072952972838101038, 0.11130862970905199, 0.055652716802045553, -0.080036865363297377, -0.074182732572421722, 0.032824623168326919, 0.031865879706987767, -0.016758259649352489, -0.0024069559974744621, 0.021075876818278654, 0.014396940553395566, -0.048207693014304456, -0.077682201861009201, -0.015232023889946197, 0.016377276604451121, -0.095289423698368528, 0.053188940442295753, -0.007014873006770006, -0.052725469767566646, 0.011169867646763318, -0.027857672235467666, -0.015475583648045285, 0.012961699649937817, -0.0041828274568459145, -0.037569795749083312, 0.035148729718712211, -0.042406682273060839, -0.023034972903402475, 0.081701490128889429, 0.068453834537836641, 0.013361892210316275, 0.33668792873864783, 0.078872059827295921, 0.043930202648898596, 0.069824545029654933, 0.084402385062259569, 0.0024106108264623457, -0.011480544789508253, -0.033818612421867275, -0.031667371615742741, -0.00502827741655568, 0.027023228844142416, -0.014247675993027993, -0.0090012900880712381, -0.061888601165125297, -0.098942005293094756, 0.01263323283761883, 0.12573571356763666, -0.076732726183971922, -0.11419210623041808, -0.037690963997301666, -0.039023576316775249, -0.037907101781325052, 0.0040436237493623132, 0.041915865623990772, -0.0086885172531639947, -0.0046973120649884377, -0.092543539386828647, -0.011130337821465186, 0.098161425946316272, -0.016521364547457937, -0.019122146869457519, 0.34253529215971379, 0.37107135970432042, 0.86635491100693385, 0.069575788839216376, -0.077933249962103462, -0.23382813425811527, -0.09863825695148333, -0.028452880392854651, 0.02313090978794155, -0.065286881502250455, -0.0065163090427018816, 0.035245400387964018, 0.064871673248635159, -0.009910405282542379, -0.074442233736364327, -0.088464389222228296, -0.043746931993730226, 0.025650630137194336, -0.019555697374510583, -0.11757758027246727, -0.088395093514236453, 0.00060518848897882292, -0.023883294001795972, 0.0052157631901946491, 0.017696001575747886, 0.038292630722847568, -0.064152245495751242, -0.00046677178847581285, 0.13417678111594261, -0.015474080963115183, -0.27741376367012288, 0.10645611255411111, 0.49136473213042497, 0.64566220514793704, -0.34894672904561197, -0.22930577035731828, 0.065707150130888653, 0.074931291451028587, 0.010801205387038527, -0.067411320951837883, -0.12930461632538409, -0.0036267311717218076, 0.039023237394102503, 0.049868939755416578, 0.0028129255950110478, -0.048592759447364742, -0.098206750946425017, -0.21167593067635238, -0.24661106106331246, 0.17003550871316073, 0.14319521757460188, -0.18791450411345456, 0.0079365608693964147, -0.026343653609827478, -0.0064127815367698198, -0.021569631345745438, 0.090706510487419845, 0.023928227792328456, -0.030645866021173446, 0.0081944297873781744, 0.27230784814155318, -0.030140868003777244, -0.55369429571770845, -0.33461208848408441, 0.40042515638418175, 0.54787745216002603, 1.7521122327973908, 0.02630887865012415, -0.7903386881329546, -0.30506910575481899, -0.058122663069862322, -0.1512000567754718, 0.030805208457151921, -0.056922369420476453, 0.018910842473248434, 0.041026550164896362, 0.10973045554309671, -0.029240486244660805, -0.3163922398944391, -0.13757198545862218, 0.12706805404252669, -0.072057201177643526, -0.070445625631504319, 0.11232623315704422, 0.053944973270972273, 0.078972221527521536, -0.00416350617181658, -0.043614575991677791, -0.077748186566067592, 0.049219693586430378, 0.042416802273341318, 0.12644205243911066, 0.39096527623299887, -0.40370312470951558, -0.48770772267188234, 0.54473301805335006, 0.92898941300038107, 0.71512595047089234, 0.74385331371246577, 1.2482467538506357, -0.57166836623813888, -0.43281940369171401, 0.023943010115746077, 0.00030176170486927514, -0.24679227625704431, 0.16824405476099039, 0.092929862040844205, -0.029175915115095267, -0.081210816145377299, 0.032535239746583337, -0.018421573557966798, -0.18224178792918325, -0.29052216822588622, 0.17110516400703951, 0.22815908460498383, -0.15890945005602886, 0.0062170702455973175, -0.044853168166470177, 0.0080331646697508995, -0.0842632097276826, -0.016095705887684088, -7.820006368420517e-05, -0.097923922056211893, 0.48676804973393806, 0.35484595526186524, 0.18102793706565945, -1.0939934957719641, 0.9231011957591454, 0.49108674920037432, -1.6725083187014362, 0.56577881124465212, 1.7716903674840145, -0.31724338121882639, -0.46946217188769945, 0.46804088594167398, 0.2495599910064237, 0.21965961682575341, -0.28426475792188616, -0.12974774398976752, -0.086262426064739897, 0.084703474730633041, 0.12060535717228706, -0.33186051211466344, -0.10830858623812538, 0.17740469505624151, -0.0043101841409599814, -0.33866696253627121, 0.0393969525618715, -0.26339593722715154, 0.13598013030776659, 0.10364961453623139, -0.011168077595515221, -0.18224776807046214, -0.12296413131763473, 0.21989373735679868, 0.16380771710007386, 0.67271651276019118, -0.16319024909703136, 0.98390492917013062, 2.2217962254133496, -0.66077224091757802, 1.5638480130127614, -1.0129349058609651, 0.99232752650122769, 1.5155025161072819, 0.58675004339129688, 0.16801904131477496, 0.31725814504248523, -0.071714399164359219, 0.13832330176015895, -0.060442748560523389, -0.043434328775720558, -0.19718551141955432, -0.36586758725441343, -0.15612863751307984, -0.72499416196393318, -0.11056584899062098, 0.082780177430082938, -0.10585857030502534, -0.49717312826403764, 0.060297128338251849, -0.1741288940590143, 0.034503572794408806, -0.022140372568936324, 0.14950061311078439, 0.29242994865853189, -0.29100003355994419, 0.22751157962613949, -0.11471557758718894, 0.25575295594566499, 0.66007876209063321, -0.28144140655806837, -1.050216159001623], "bottom": {"real": [5837.242024063461, 16421.190569053251, 8527.2275182117319, 3938.41845467232, 969.21000043764093, 269.87346390668256, 193.53992011404216, 135.14631901451781, 58.523829955616236, 29.311524302602624, 20.863281166035843, 14.808244655557742, 11.401757474781434, 9.6386493936653288, 8.0538591094176777, 7.4714398062124641, 7.2711776175262477, 7.4714398062124241, 8.0538591094176653, 9.638649393665327, 11.401757474781416, 14.808244655557761, 20.863281166035858, 29.311524302602734, 58.523829955616236, 135.14631901451779, 193.53992011404216, 269.87346390668273, 969.21000043764059, 3938.4184546723186, 8527.2275182117301, 16421.190569053259, 14982.474506443205, 4755.1985786680198, 3083.8642011691682, 1502.8058566448747, 596.26851293203197, 262.35520536881347, 201.52671733896852, 108.13603201941604, 57.999252895466505, 28.992297197792773, 19.233600066307524, 13.871180161098291, 10.761421361485239, 8.8506828512506655, 7.8125584926262617, 7.1234501229356377, 6.6618841738486525, 7.1762884709606958, 8.1416889804019927, 8.8663632607744063, 11.328672883224474, 14.773155762512141, 20.058056528572401, 28.688344099071909, 57.995376649220766, 126.94075502116374, 211.30702516808583, 262.54448595536269, 752.96024761972672, 2141.5368023327851, 3891.3661973657881, 5583.5460990085166, 9500.6528086677681, 4225.3848573287933, 1708.1934289355966, 562.9450489663958, 526.41420210821127, 290.26438540271903, 151.34732218304231, 89.304597077577796, 50.161041673414204, 27.913649363990057, 19.253723232543383, 14.310501425792925, 10.622538104620221, 8.6520747616015274, 7.6796432895344049, 6.9311889072641995, 6.5442165971044535, 6.8352472963781299, 7.2531587818267287, 8.8125700197174819, 10.256278085266516, 13.962876672312342, 18.680978544191657, 27.641611280844351, 51.70590284164593, 92.383035108628647, 175.23999552329536, 238.6238705323166, 446.2130896894364, 589.84950644143407, 1459.800700420838, 2983.6421667816849, 6639.6980670311659, 3615.9230624353354, 1045.9683029029898, 541.95698246733389, 551.2292078912385, 269.01748102372522, 121.81702466062521, 78.619304368994548, 46.266475261147157, 26.895860643920866, 18.966790938143397, 14.141216689343503, 10.624588727031869, 8.5120674584469107, 7.5084078178458054, 7.127942028564302, 6.5505226441602824, 6.5976614570745582, 7.5671769464103935, 8.4293753366132105, 10.446992271959362, 13.755112829750516, 18.515159362788847, 25.868606557166544, 43.593025933023839, 76.847860046206122, 125.03651194390009, 241.77755195523335, 523.08499299660048, 633.36555854289827, 855.69476790230692, 2473.0066947172836, 2803.9496552692158, 1850.6397468329815, 645.77314636959022, 447.00563976857688, 381.6372474267435, 196.7755433249871, 93.586816816137713, 72.612365740746853, 43.711482900929944, 27.782222145016455, 18.647268962291367, 13.844325818813939, 10.350953817244738, 8.7140342291788997, 7.3516641952705921, 6.886474159529909, 6.1226840838853605, 6.7186507131414492, 7.0148796268173141, 7.9020046012660439, 9.6536048586555445, 13.744699851401226, 17.872279960942961, 25.025889030797831, 38.775468597625853, 68.228616115874232, 97.806667851841851, 178.52012140604563, 362.80035213144055, 506.96305099601244, 508.39785887877764, 1271.7627571663681, 1970.0713813014595, 1211.9708052435174, 448.82610946730466, 388.91483700385658, 294.26646682781944, 147.68595350614527, 76.563818921604266, 57.53338100926053, 37.094893522892434, 25.168957624787666, 17.405368368923924, 12.931100650610183, 10.211256703797424, 8.3783082007702649, 6.8845647672161352, 6.2954630311895761, 6.1711525716900475, 6.2852637654745882, 7.1318002457288037, 7.8126464614985656, 9.3562580058102185, 12.456351018938477, 17.12365920646733, 22.851195586808732, 34.22708028584924, 57.273965595477058, 82.126801795962663, 114.0182287640494, 248.69630423497404, 396.25330408596494, 382.77662480077174, 646.86127060880278, 1066.1959535576209, 770.74521634546068, 344.77171326148743, 328.21152488387929, 223.6732358532505, 124.39997807543456, 66.042124980023658, 49.833499278727395, 36.080812354502541, 24.010000461569632, 16.241873057458708, 11.94801819096188, 9.663620823552634, 8.1734814269274114, 6.8263151892657215, 6.1094005534487641, 5.5560364128866686, 6.0895161225652297, 6.7371553850929349, 7.5614980494024309, 9.3035460399083956, 11.904089289656403, 15.753116847312391, 22.116186781807659, 33.401673308199264, 45.775772578430164, 60.542807604917392, 91.412240599681226, 164.58197575163385, 289.97996145010069, 367.09809735899603, 433.51729458889321, 538.99714678707312, 436.9735160339028, 214.45700003213321, 226.27353228768246, 190.72528880002966, 104.24105887361294, 63.096291028318397, 44.895692115573439, 34.378771692376986, 22.743616534455565, 14.943130070579064, 11.306154107321404, 9.490726177732828, 8.057539867813432, 6.6249917656310648, 5.80530835366394, 5.6727594087192097, 5.6803994552580717, 6.4512234217242401, 7.3091820969038412, 8.9053168517088785, 11.147020821533502, 13.978280110937801, 19.308076160210224, 27.380246565144393, 39.45337202937899, 53.823906577074318, 70.344986096434994, 125.18491111128095, 220.27402210338624, 239.02095964975129, 253.52805611791149, 337.90351459641266, 301.96557555412664, 163.14862666039761, 176.1964602722083, 155.2640157056349, 91.402170635856464, 57.9223709952281, 44.132906015519367, 31.477476121951593, 19.994026988386896, 14.669539930994308, 10.52228821405147, 8.9499245076841039, 7.4160858759377133, 6.0652235290917815, 5.5336836347972449, 5.3453345499720717, 5.5251998880295652, 6.049359333226362, 7.4989544637904162, 8.6441457110666775, 10.17217183387787, 12.922358897606493, 16.497840703497662, 24.473525343022239, 36.893225695749557, 46.488420049270196, 60.008949681859932, 101.55324863109907, 172.57213254679178, 186.99209153386954, 177.69050334329205, 224.71115213821386, 219.35401697656656, 131.98987419085088, 139.58612051732652, 133.97255559613592, 83.717667146076977, 52.48988401855501, 39.948149447326976, 29.992243301791387, 19.792293237610107, 13.550533437330548, 10.191515700191694, 8.7251285565258279, 7.0130966099588283, 6.0412905020539274, 5.1695609453883042, 5.2335464217738297, 5.2748937601741233, 5.9151189141762242, 6.6404924993384702, 8.2525633194354864, 9.86505065394757, 11.857247340871206, 15.821441128307118, 21.809722446513074, 33.204842604441232, 45.292100524348356, 53.588237099743992, 86.046570865333351, 138.34470513151498, 144.28338780744133, 137.47961624149846, 184.36194031651988, 195.75342807998493, 111.69248606852764, 108.22432115257934, 110.31477193613289, 70.737952758714727, 45.316084199866118, 35.221165127538775, 25.132656234698317, 17.808708328889288, 13.205138111505544, 9.8319756958291311, 8.3851869715659166, 6.6979147774695411, 5.5577949826186863, 5.0640136062157843, 4.915924071395172, 5.1747375238921238, 5.621849200313962, 6.3926522738336189, 7.8838556827841355, 9.5880495388903295, 11.2702077240605, 14.675031755421674, 19.961485602094882, 28.340133288775018, 38.071807302485439, 46.211805942072893, 70.896288173632385, 115.46517347822648, 119.73189294718243, 110.21299532100423, 137.49078682179271, 145.21239507063157, 90.736475500202033, 84.820134662984671, 87.028887735269166, 58.592780926558561, 39.140475702708279, 30.030590426700837, 23.561748558876776, 17.07171690133644, 12.059539847971436, 9.1303469917064106, 7.8824915696444915, 6.5574577163884546, 5.4508327472027114, 4.6319487889521351, 4.8556936697031698, 5.0952139895658046, 5.3626061482546099, 6.0231504777267588, 7.3929124540904754, 8.8697139145390977, 10.839689611405708, 13.484775054138881, 18.125520870461621, 25.51326521840107, 33.814821758723511, 43.000660712963757, 59.289067569698993, 89.523100480006264, 96.002834501641331, 88.248403951505949, 114.22377327036351, 136.39754784126862, 85.133937137494058, 70.757613581150778, 72.941136445046624, 50.417131404641061, 35.29732697897046, 27.270836476117434, 20.337170678571646, 15.328073102208975, 11.227626008224609, 8.8894703399043831, 7.3066732490559749, 6.1117570010797007, 5.2240396800825373, 4.6927944754047974, 4.4462691956761722, 4.8997943453914665, 5.030120352863074, 5.8504059853355139, 6.9021030785339024, 8.6936481048417669, 10.014300550884089, 12.783989981758968, 16.915182814597465, 23.037925870268026, 29.781501040175886, 36.137786628939317, 46.659613314406016, 72.542096884400053, 80.098309124706176, 76.772569487169434, 87.96931966600367, 100.35528192635977, 66.390743706670477, 56.297220892995071, 57.449776912893739, 42.311817158751616, 30.3950154293208, 24.633397231162515, 18.800000277532757, 13.906344724637259, 10.685561612392762, 8.7726633553262534, 7.3765778222647702, 6.0238006847876884, 5.2082074327498917, 4.5886719664559399, 4.6166410353959701, 4.6115148269570563, 4.9920772309179853, 5.6445124723515159, 6.7980683034043956, 7.4975962005340229, 9.6066586869080837, 12.06652313244911, 16.251663826300451, 21.003823940677318, 27.622933573071261, 33.836788273738165, 41.877854423447168, 58.532835369703733, 68.239320296767971, 62.326149739008642, 71.583830955139561, 81.611391861674278, 61.290606958364052, 49.963935495056425, 49.744027747354536, 38.555965139546075, 28.313241453630848, 22.404376161451804, 17.69007471646092, 13.348861556136391, 10.46674681432143, 8.1507364253451371, 6.7804933363165576, 5.7303247473649055, 4.9737961537641029, 4.6591603640435002, 4.3691273542401863, 4.5876440856805019, 4.8775279419312838, 5.5643157824289267, 6.5106424514611394, 7.6517258774697732, 9.3978270674435489, 12.208065227131403, 15.316902707006163, 20.982661610009878, 26.317869298031479, 32.313264556083666, 38.062860113666069, 50.105167484748129, 58.219661124046276, 59.861157680092496, 61.019020366450448, 69.642805168222978, 57.634563900118565, 46.967523109902231, 43.85887184820529, 34.957640557066753, 26.505039983640941, 21.71545909517323, 16.236775633007202, 12.875777825381771, 9.9834739610997616, 7.8930497070448, 6.6829715569803687, 5.6672193704500122, 5.0794246915518446, 4.681264046528276, 4.6186632297539427, 4.6211092904368423, 5.0897865773957403, 5.5691810710797096, 6.5035187317991277, 7.2014326998593701, 9.1703620437944053, 12.475838151915292, 15.6623920524898, 19.528853843427893, 25.077028045409381, 30.833489267211636, 39.113296407386123, 46.770991427918581, 55.937485646720255, 57.570245344451905, 57.637817589030213, 59.860391566610879, 52.776195682327739, 43.839296967027423, 39.833012832197511, 33.598213438023812, 26.072700016436322, 20.276975501521189, 16.195672497297036, 12.239826665071417, 9.6875213487165066, 7.7379543893933631, 6.3786782254807157, 5.5950259327840586, 5.1232903828272987, 4.5786592309737237, 4.5178908466513406, 4.578659230973722, 5.1232903828272978, 5.5950259327840595, 6.3786782254807148, 7.7379543893933684, 9.6875213487165084, 12.239826665071421, 16.195672497297036, 20.276975501521182, 26.072700016436304, 33.598213438023805, 39.833012832197511, 43.839296967027416, 52.776195682327739, 59.860391566610836, 61.019020366450441, 57.570245344451919, 55.937485646720255, 46.770991427918517, 39.113296407386116, 30.833489267211636, 25.077028045409371, 19.528853843427893, 15.662392052489803, 12.475838151915296, 9.1703620437944, 7.2014326998593674, 6.5035187317991294, 5.5691810710797114, 5.0897865773957385, 4.6211092904368405, 4.6186632297539418, 4.6812640465282733, 5.0794246915518446, 5.6672193704500087, 6.6829715569803696, 7.8930497070447982, 9.9834739610997651, 12.875777825381764, 16.236775633007209, 21.715459095173252, 26.505039983640959, 34.957640557066775, 43.858871848205247, 46.967523109902238, 57.634563900118593, 69.642805168222949, 71.583830955139618, 59.861157680092482, 58.21966112404624, 50.105167484748137, 38.062860113666083, 32.313264556083666, 26.317869298031482, 20.982661610009867, 15.316902707006161, 12.208065227131399, 9.3978270674435436, 7.6517258774697732, 6.5106424514611394, 5.5643157824289284, 4.8775279419312856, 4.587644085680501, 4.3691273542401872, 4.6591603640434966, 4.973796153764102, 5.7303247473649037, 6.7804933363165594, 8.1507364253451371, 10.466746814321423, 13.348861556136386, 17.69007471646092, 22.404376161451779, 28.313241453630852, 38.555965139546039, 49.744027747354508, 49.963935495056418, 61.290606958364073, 81.611391861674235, 87.969319666003628, 62.326149739008663, 68.239320296767957, 58.532835369703712, 41.877854423447168, 33.836788273738158, 27.622933573071283, 21.003823940677311, 16.251663826300444, 12.066523132449106, 9.6066586869080837, 7.4975962005340202, 6.7980683034043929, 5.6445124723515194, 4.9920772309179844, 4.6115148269570598, 4.6166410353959719, 4.5886719664559372, 5.2082074327498908, 6.0238006847876866, 7.3765778222647711, 8.7726633553262587, 10.685561612392764, 13.906344724637258, 18.80000027753276, 24.633397231162533, 30.395015429320793, 42.311817158751609, 57.449776912893732, 56.297220892995028, 66.390743706670506, 100.35528192635975, 114.22377327036349, 76.772569487169449, 80.098309124706219, 72.542096884400053, 46.659613314406002, 36.137786628939317, 29.781501040175893, 23.037925870268033, 16.915182814597468, 12.783989981758975, 10.014300550884082, 8.6936481048417651, 6.9021030785339006, 5.8504059853355228, 5.0301203528630678, 4.8997943453914612, 4.4462691956761731, 4.6927944754048019, 5.2240396800825328, 6.1117570010797033, 7.3066732490559732, 8.8894703399043813, 11.227626008224609, 15.32807310220897, 20.337170678571649, 27.270836476117438, 35.297326978970467, 50.417131404641083, 72.941136445046638, 70.757613581150821, 85.133937137494073, 136.39754784126856, 137.49078682179265, 88.248403951506006, 96.002834501641317, 89.523100480006221, 59.289067569698972, 43.000660712963743, 33.814821758723504, 25.51326521840107, 18.125520870461624, 13.484775054138884, 10.839689611405714, 8.8697139145391013, 7.3929124540904745, 6.0231504777267606, 5.362606148254609, 5.0952139895658073, 4.8556936697031734, 4.6319487889521325, 5.4508327472027096, 6.5574577163884502, 7.8824915696444933, 9.1303469917064124, 12.059539847971436, 17.071716901336433, 23.561748558876776, 30.030590426700826, 39.140475702708279, 58.592780926558504, 87.028887735269151, 84.820134662984628, 90.73647550020199, 145.21239507063157, 184.36194031652002, 110.2129953210042, 119.73189294718252, 115.4651734782265, 70.896288173632399, 46.211805942072878, 38.071807302485439, 28.340133288775018, 19.961485602094882, 14.675031755421674, 11.270207724060505, 9.5880495388903277, 7.8838556827841364, 6.3926522738336233, 5.6218492003139611, 5.1747375238921247, 4.9159240713951746, 5.0640136062157861, 5.5577949826186863, 6.697914777469534, 8.3851869715659184, 9.8319756958291471, 13.205138111505548, 17.808708328889303, 25.13265623469832, 35.221165127538775, 45.316084199866097, 70.737952758714712, 110.3147719361329, 108.22432115257941, 111.69248606852767, 195.75342807998487, 224.71115213821398, 137.47961624149852, 144.28338780744139, 138.3447051315149, 86.046570865333337, 53.588237099743992, 45.292100524348342, 33.204842604441232, 21.809722446513071, 15.82144112830712, 11.857247340871208, 9.86505065394757, 8.2525633194354864, 6.640492499338464, 5.9151189141762242, 5.2748937601741233, 5.2335464217738297, 5.1695609453883096, 6.0412905020539247, 7.0130966099588328, 8.725128556525835, 10.191515700191694, 13.550533437330545, 19.792293237610103, 29.992243301791408, 39.948149447327012, 52.489884018555003, 83.717667146076963, 133.97255559613592, 139.58612051732649, 131.98987419085091, 219.35401697656653, 337.90351459641266, 177.69050334329216, 186.99209153386957, 172.57213254679186, 101.55324863109904, 60.008949681859939, 46.488420049270182, 36.893225695749557, 24.473525343022239, 16.497840703497666, 12.922358897606497, 10.172171833877869, 8.6441457110666775, 7.4989544637904126, 6.0493593332263602, 5.5251998880295634, 5.3453345499720717, 5.5336836347972458, 6.0652235290917806, 7.4160858759377133, 8.949924507684111, 10.522288214051471, 14.669539930994311, 19.994026988386889, 31.477476121951593, 44.13290601551936, 57.922370995228107, 91.402170635856407, 155.2640157056349, 176.1964602722083, 163.14862666039767, 301.96557555412653, 538.99714678707278, 253.5280561179116, 239.02095964975135, 220.2740221033861, 125.18491111128091, 70.344986096434923, 53.823906577074332, 39.453372029378947, 27.380246565144397, 19.308076160210238, 13.978280110937792, 11.147020821533506, 8.9053168517088821, 7.3091820969038439, 6.4512234217242392, 5.6803994552580681, 5.6727594087192079, 5.8053083536639392, 6.6249917656310666, 8.0575398678134302, 9.4907261777328298, 11.306154107321403, 14.943130070579077, 22.743616534455558, 34.378771692377001, 44.895692115573446, 63.09629102831839, 104.241058873613, 190.72528880002952, 226.27353228768243, 214.45700003213315, 436.97351603390234, 1066.1959535576209, 433.51729458889332, 367.09809735899597, 289.97996145010086, 164.58197575163376, 91.412240599681198, 60.542807604917364, 45.775772578430157, 33.401673308199264, 22.116186781807642, 15.753116847312398, 11.904089289656401, 9.303546039908392, 7.5614980494024362, 6.7371553850929322, 6.0895161225652297, 5.5560364128866677, 6.1094005534487605, 6.8263151892657268, 8.1734814269274096, 9.6636208235526269, 11.948018190961893, 16.241873057458722, 24.010000461569639, 36.080812354502534, 49.833499278727416, 66.042124980023658, 124.39997807543459, 223.67323585325056, 328.21152488387946, 344.77171326148749, 770.74521634546022, 1970.0713813014588, 646.86127060880301, 382.77662480077169, 396.2533040859646, 248.69630423497387, 114.01822876404945, 82.126801795962692, 57.273965595476994, 34.227080285849233, 22.851195586808739, 17.12365920646733, 12.456351018938475, 9.3562580058102167, 7.8126464614985656, 7.1318002457288037, 6.2852637654745882, 6.1711525716900475, 6.2954630311895787, 6.8845647672161387, 8.3783082007702561, 10.211256703797417, 12.931100650610194, 17.405368368923916, 25.168957624787648, 37.094893522892434, 57.533381009260523, 76.563818921604266, 147.68595350614524, 294.26646682781944, 388.91483700385658, 448.82610946730478, 1211.9708052435162, 2803.9496552692158, 1271.7627571663686, 508.39785887877753, 506.96305099601244, 362.80035213144049, 178.52012140604566, 97.806667851841823, 68.228616115874146, 38.775468597625853, 25.025889030797838, 17.872279960942965, 13.74469985140124, 9.6536048586555427, 7.9020046012660456, 7.0148796268173159, 6.7186507131414386, 6.1226840838853631, 6.8864741595299055, 7.3516641952705939, 8.7140342291789032, 10.350953817244736, 13.844325818813926, 18.64726896229136, 27.782222145016434, 43.711482900929937, 72.612365740746867, 93.586816816137656, 196.7755433249871, 381.63724742674339, 447.00563976857711, 645.77314636959011, 1850.6397468329806, 6639.698067031165, 2473.0066947172841, 855.69476790230704, 633.36555854289816, 523.0849929966007, 241.77755195523349, 125.03651194390012, 76.847860046206094, 43.593025933023853, 25.868606557166554, 18.515159362788847, 13.755112829750496, 10.44699227195936, 8.4293753366132176, 7.5671769464103917, 6.5976614570745609, 6.5505226441602797, 7.1279420285643011, 7.5084078178458018, 8.5120674584469089, 10.624588727031863, 14.141216689343507, 18.966790938143404, 26.895860643920869, 46.266475261147157, 78.619304368994577, 121.81702466062518, 269.01748102372505, 551.22920789123827, 541.95698246733446, 1045.9683029029893, 3615.9230624353345, 9500.6528086677681, 2983.6421667816867, 1459.8007004208382, 589.84950644143373, 446.2130896894364, 238.62387053231652, 175.23999552329542, 92.383035108628548, 51.705902841645944, 27.641611280844312, 18.680978544191667, 13.962876672312348, 10.256278085266519, 8.8125700197175014, 7.2531587818267402, 6.8352472963781254, 6.5442165971044517, 6.9311889072641986, 7.6796432895344058, 8.6520747616015186, 10.622538104620222, 14.310501425792921, 19.253723232543386, 27.913649363990046, 50.161041673414196, 89.304597077577839, 151.34732218304231, 290.26438540271909, 526.41420210821127, 562.94504896639603, 1708.1934289355961, 4225.3848573287914, 14982.474506443201, 5583.5460990085194, 3891.3661973657854, 2141.5368023327842, 752.9602476197266, 262.54448595536314, 211.30702516808589, 126.94075502116375, 57.995376649220773, 28.688344099071898, 20.058056528572401, 14.773155762512143, 11.328672883224492, 8.8663632607744169, 8.1416889804020052, 7.1762884709607109, 6.6618841738486552, 7.1234501229356262, 7.8125584926262608, 8.8506828512506655, 10.761421361485221, 13.871180161098282, 19.233600066307527, 28.992297197792769, 57.999252895466505, 108.13603201941609, 201.52671733896852, 262.35520536881353, 596.26851293203185, 1502.8058566448742, 3083.8642011691672, 4755.1985786680179], "imag": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, "imag": [0.0, -0.13871927675086548, -0.55827571496092165, -0.85383263724083092, -0.46701193383615952, 0.047809537946767316, -0.3426920474802197, 0.15831589374137059, 0.20589774809377467, -0.21944378817656213, -0.12654666933988368, -0.19544466487854822, -0.096446269450190733, -0.2265628704083768, 0.10840537931386109, -0.055457499657790024, 0.0, 0.055457499657784036, -0.10840537931384955, 0.22656287040837902, 0.096446269450189581, 0.19544466487854698, 0.1265466693398799, 0.21944378817656165, -0.20589774809377467, -0.15831589374137045, 0.34269204748022014, -0.047809537946767663, 0.46701193383615958, 0.85383263724082981, 0.55827571496092188, 0.13871927675086551, -1.279550584464304, -0.80695903119669909, -0.5066002982483182, -1.4561178571661, -0.13414323846201096, -0.25055297311154978, -0.64027404819684697, -0.19951473048424853, 0.33137912963317689, -0.047498243759329373, -0.24064979253244873, 0.1395435157196932, 0.17789090203786645, -0.044374683774468404, -0.36806933922228663, 0.19528135326361654, 0.32999646269142918, -0.13015393203002068, 0.1959681784626765, 0.28571835588367483, 0.22705801052414751, -0.14508578389537852, 0.23457709892354506, 0.21163882649779878, -0.20150153854805303, 0.060946718066345364, 0.31004041285864392, -0.26258387187077481, 0.13955938857309563, -0.2449447076348287, 0.068510051090992716, -1.1758121586434307, -1.3682011726001491, -1.6899938281650106, -0.75019268688655194, 0.68911282851752931, 0.44086751164178845, -0.23155959034339019, -0.54510315877084203, 0.26585764990704713, 0.2624959582084973, -0.045146091846171167, 0.011517786678353055, -0.057907844320373422, 0.076854124253010955, 0.0035092112107953564, -0.23588030048916386, -0.048277840380711062, -0.21489820671404902, 0.26817287178185911, -0.027468417557455105, 0.29873174456545087, 0.072899207479288589, 0.023829771073631354, 0.0066491170121291982, 0.0048298736467919092, -0.13057930804766787, -0.26946209593128789, 0.39240115086507854, 0.85317721061416651, -0.11752989424299426, 0.99998707389668307, 0.55787436167296656, -1.0314907737559058, -1.5043263157657751, -1.3830149751777026, -0.20446853801739454, -0.28364033101911829, 0.68303147161529054, 0.03902721990981517, -0.18273076418691683, -0.013466871818696777, 0.028891896132398606, -0.068621242214500722, -0.026792869524675422, 0.18842973810535804, 0.068392282535087487, -0.17097210848547728, -0.23288712244734344, -0.31416612358708279, 0.27284169211913839, -0.072438990225411531, 0.1202436118466384, 0.067443262257117789, -0.0040452326400191837, -0.15231382612170316, 0.020438609938141238, 0.10773227995476008, -0.1646961578394999, -0.034526357222664054, 0.34045699686317549, 0.3582562672018631, 0.50203282214727851, -0.90050578840564854, 1.1141651843864948, -1.0547694551235964, 0.15113295529591572, -0.42011997382980693, -1.1306687523196075, 0.59904899563122349, 0.12782476119556421, -0.26469201284371935, -0.0028564673365734715, 0.25878855026626812, 0.00078645201341040851, -0.17181320075266995, 0.013143422355015144, 0.050238396740841747, 0.072156922858987921, 0.041604673982592558, -0.0112715468909201, -0.016101389694688766, -0.31269363662275068, 0.094974334151300396, -0.067869387966968731, 0.14811752790953986, 0.010724993842567176, -0.071495661604357427, -0.037154193600070516, 0.05779163337935353, 0.094937931062917286, -0.1272883333159909, -0.017408469796284951, 0.62003762386212868, 0.074417004563971767, -0.12904495334293709, -1.4636084587609224, -0.036743730082780299, -0.77470195486186, -0.60686918965332004, 0.62872636200456622, -0.46894237369717673, -0.18368617361305525, 0.26816260026693856, 0.13447215443394772, 0.041008154934057489, -0.017107991149299187, -0.0680369475996533, -0.01240601732954733, 0.059356927461924558, 0.12615835751185028, 0.053672308552617974, 0.037749743592242362, -0.10597138547480481, 0.022172258797336697, -0.058681894939636686, 0.12792055300624736, -0.060154323470105943, 0.032298308999269151, -0.015701976416949814, 0.0007367446231764113, 0.07712597645572114, 0.055644084876628218, -0.065981740018919821, -0.14567847379973994, 0.039540586335773356, 0.24917063867041539, -0.6432675339043068, 0.36754664039982821, -0.069649876974948205, 0.2952063874323651, -0.011417685121456805, -1.0450941798657123, -0.17374602795438093, 0.19930768818407893, -0.1556206112582631, 0.055817128538076119, 0.17900936499027481, 0.0052374024484263433, -0.042596098823283683, -0.020834261928533829, 0.036686593450547957, 0.00076880260958178746, 0.061719692477810909, 0.027700064466885506, 0.10661165850403151, -0.13945785164947772, -0.050544258155460708, -0.11937846922953205, 0.034399558837173369, 0.02521789889898339, -0.0047468618190418272, -0.028256857448446509, 0.029239764799856467, 0.046539635523190424, -0.033316423902485758, -0.15842957449053302, 0.14741145606308881, 0.1752194519195491, -0.075629282019522973, -0.91968716746889456, -0.45504702044764878, -0.12177128356890227, -0.10462001130261206, 0.1804016527269576, -0.047005347129142915, -0.19200701617075638, 0.1133507155810588, 0.14315736053243958, -0.067545520465919412, -0.0022066601859379556, 0.0084001340729883862, 0.0049338006978182497, 0.016196236751475263, 0.066419574069183293, -0.012240235365542763, 0.014433128838537781, 0.023674926554515096, -0.017114978073907605, -0.14379898617791109, 0.084303725772385638, 0.0056779923265117766, 0.048321020114327755, 0.022849248448590556, 0.012395399307186424, 0.023865191872962568, -0.0060375152678495622, -0.021897165586045375, -0.092173327221095822, -0.17587869380621102, 0.091992925322100877, 0.014980171338571357, 0.0086063527817299007, 0.098130362913165001, 0.130338354535672, 0.094203305677599436, -0.19730813245926024, -0.10041080102702246, 0.075932042745803577, 0.021134287621881041, -0.059575532325244519, 0.0082992125422013917, 0.015272751179096435, 0.0090248536114060952, -0.019251438395354602, -0.015788972597383853, 0.003250945662643008, 0.032968975479440332, -0.015058032290062954, 0.017643635220692015, -0.068857003169531444, 0.030070570566209365, -0.12793975070954616, -0.043534834876013277, 0.0017423598025251587, 0.028592812647173952, 0.0036967486190428137, 0.0090453167838253218, -0.016259252101333162, -0.0052076817297827191, 0.018046490267575967, -0.026168992315517926, -0.083077034868877908, 0.032798990212655046, -0.027546636940189637, -0.12409557570816734, -0.0088356071822206419, -0.034333260533212014, 0.056177295075322577, 0.057994393113440637, -0.032566171451921935, -0.010465798382464413, 0.020079425427313443, -0.017396471253638312, -0.017870360905994483, 0.010919865739753416, 0.012363627500514531, -0.014313215415815719, 0.034201557920817111, 0.10796975372596095, 0.21064351401911494, -0.046016401957029461, -0.083803307026608195, -0.15128619138644472, -0.070745415327885988, -0.06951503451164355, 0.070638418064815661, -0.0080019035626016123, 0.02654845111893607, 0.019862790828204354, -0.0038698355762961699, -0.020866890117519677, -0.0048084432749845906, -0.0036451640963975299, -0.0026341276434262537, -0.0023099396961399252, 0.043908100663757832, 0.096772033742844044, 0.04918905117971889, 0.070600771528734244, -0.026773628633603139, -0.05973226936448265, 0.013552286601009756, 0.0032289712715403013, -0.023046337829044056, -0.010405151723131599, 0.019427096259341661, 0.018256051683800133, -0.056759870222018213, -0.032960560091393032, -0.017399004273096801, 0.12185282737260948, 0.037984234337547759, -0.12270039394272629, -0.055947445070994896, 0.04621845892784341, -0.1857137618225036, -0.063810892086947277, 0.07053672326645595, 0.03906718120738633, 0.01577517031413787, -0.0047456570119124219, -0.019801411975277249, 0.024637199092487311, 0.019457415538283646, -0.0067710516924983765, -0.028346017012320548, -0.0065098329229195291, -0.0013352981567293612, -0.023151697315085987, -0.06519950655383687, -0.1042256990435297, -0.023224778901070588, 0.0072963544312097499, -0.034708386400909293, -0.0023809918522597964, 0.016167541663812009, -0.0054340236987463174, 0.0027208109883995382, 0.0066841964011560553, -0.023932611476630623, -0.039005803210705049, 0.016396390506166533, 0.058067983872416511, 0.15112210885707694, 0.091953539825334094, -0.30878981646904202, 0.23773808023122034, -0.11496619447060334, -0.15410637679425881, 0.16398523681965352, -0.042160123131605852, -0.0057838899052936009, 0.021159112170606904, -0.00089194024953267548, -0.043632529695206895, 0.025594278576058532, 0.04827977682701759, 0.008656372871138766, -0.025123508323013247, -0.021222695138498154, 0.012727342518660707, 0.0027327273075534982, 0.016727222825244921, -0.062351999721191126, -0.095086485134993501, 0.009789315676701096, 0.020526538310240195, -0.0017729650965628652, 0.0042367342426148282, 0.0096918161659650896, -0.014518052118829046, -0.0045810353743419062, -0.069451504639382974, -0.036289590762181183, 0.22708969966074724, 0.3665342673930222, -0.59625762091850343, 0.36774772058649913, 0.2627894192943519, -0.10703484463450727, -0.13950623949894417, -0.080806137493453625, -0.018761372187400087, -0.14807112835842171, 0.015533054933596769, -0.013156170331620851, 0.033053795814988449, 0.012901798319361097, 0.00377273404398293, -0.0087050725477761486, -0.026200884444946355, -0.069983120002794802, -0.045851102565308888, -0.0080215915074989454, -0.048063119818621811, -0.055909188302307122, 0.015567334119352533, -0.069113723063989038, -0.0072474243428621863, 0.042526562044161449, 0.039806904605142766, -0.027582495505149615, 0.0049261433754253488, -0.04754279011762625, -0.074122351214929355, 0.22223052073714245, -0.0077417273670608014, -0.080825286420915882, 0.62390518555370134, -0.18988737165776209, -0.013310214088767865, -0.13538196037364333, -0.38370084094404533, -0.32853605615928166, -0.23095349816725769, 0.0024573680892555631, -0.060148519255615616, 0.030360401196840837, -0.0084214745912351684, 0.033027304906486414, 0.030739864021722699, 0.0097039926999964547, -0.060698557392036313, -0.04746519227839175, -0.053530300822357386, -0.046174234082698239, 0.0052112044634668644, 0.020596160958060543, -0.15878906772803425, -0.02193329327109424, 0.0079342423803948989, 0.027919450264392773, 0.0070430807013025284, 0.020085011074200358, -0.075935043112063852, -0.041652272773957957, 0.16889959178925665, 0.18647295463827135, 0.014154913289577202, 0.060521107576468509, -0.63637312336837226, 0.90191957381577925, 0.27922666018549325, -0.48673882008850766, -0.33370031357259311, -0.3513559741645329, -0.2622745616761345, -0.17197492325356289, -0.0057315027931728107, -0.04469091820910031, 0.053283308768024513, 0.039465599946438043, -0.014174493080162808, -0.04946368453866596, -0.022857667990750397, -0.16881304090125779, -0.080812526909499666, 0.00058234446045326367, 0.084338425974376116, -0.12518231680881081, -0.040506406338050263, -0.08374410530080309, -0.001350704062862842, 0.030200114248426838, 0.042299889130269351, 0.011310849850443803, 0.081759178445920677, 0.035289563063634108, -0.12268118169149339, 0.058157249768416797, 0.19453030659557169, -0.32019695139960919, 1.2621857571853639, 0.24947885513851267, -0.86827138009746829, -0.26346017480935824, -0.40312000272498982, -0.22820944406317872, 0.0813299683314989, 0.0038871367161771898, -0.10891614598123801, 0.00043800911235750134, -0.043289297526922212, 0.012593343881888895, -0.0042686636670869688, -0.008833188472774969, -0.024141666831206311, -0.03923653338217388, -0.25435751708422283, 0.0, 0.13106581900761921, 0.22441332039311174, -0.094768447018961463, 0.0024812204443016074, 0.036286268392045858, 0.0056698536823630895, -0.010705526315732189, 0.082044918512892334, 0.0010532751265284197, 0.015101105715477795, 0.090653505450224631, 0.27613653682211553, 0.17494988738073455, 0.30956821269851931, 0.034540476362705901, 0.0, -0.034540476362701966, -0.30956821269851847, -0.17494988738073411, -0.27613653682211575, -0.090653505450227129, -0.015101105715478274, -0.0010532751265287738, -0.082044918512892334, 0.010705526315732109, -0.0056698536823632534, -0.036286268392045982, -0.0024812204443015839, 0.094768447018961921, -0.22441332039311199, -0.13106581900761924, -0.00058234446045052997, 0.25435751708422299, 0.039236533382173068, 0.024141666831206714, 0.0088331884727747174, 0.0042686636670865837, -0.012593343881888968, 0.04328929752692242, -0.00043800911235742604, 0.1089161459812379, -0.0038871367161778295, -0.081329968331500163, 0.22820944406317803, 0.40312000272499021, 0.26346017480936107, 0.86827138009747229, -0.24947885513851314, -1.2621857571853632, 0.32019695139960752, -0.19453030659557399, -0.058157249768416672, 0.12268118169149329, -0.035289563063634503, -0.081759178445920649, -0.011310849850443772, -0.042299889130268831, -0.030200114248426623, 0.0013507040628626303, 0.083744105300802674, 0.040506406338050895, 0.12518231680881103, -0.084338425974375644, 0.046174234082696289, 0.08081252690949875, 0.16881304090125832, 0.022857667990750261, 0.04946368453866589, 0.014174493080162815, -0.039465599946438334, -0.053283308768024965, 0.044690918209100129, 0.0057315027931725314, 0.1719749232535624, 0.26227456167613394, 0.3513559741645329, 0.33370031357259339, 0.48673882008850966, -0.27922666018549408, -0.90191957381578058, 0.63637312336837637, -0.060521107576467635, -0.014154913289577215, -0.18647295463827063, -0.16889959178925743, 0.041652272773958374, 0.075935043112064518, -0.020085011074200423, -0.0070430807013025249, -0.02791945026439252, -0.0079342423803949891, 0.021933293271094115, 0.15878906772803389, -0.02059616095805995, -0.0052112044634673969, 0.008021591507499164, 0.05353030082235833, 0.047465192278390397, 0.060698557392036064, -0.0097039926999963107, -0.030739864021722747, -0.033027304906486449, 0.008421474591235165, -0.03036040119684083, 0.060148519255615165, -0.0024573680892551967, 0.23095349816725866, 0.32853605615928261, 0.38370084094404566, 0.13538196037364647, 0.01331021408876858, 0.18988737165776171, -0.62390518555369734, 0.080825286420916714, 0.0077417273670609497, -0.22223052073714233, 0.074122351214929452, 0.047542790117626306, -0.0049261433754249724, 0.027582495505149737, -0.039806904605142586, -0.042526562044161428, 0.0072474243428622999, 0.069113723063988594, -0.015567334119352927, 0.055909188302305984, 0.048063119818621922, -0.0027327273075535455, 0.045851102565308353, 0.069983120002795093, 0.026200884444945973, 0.0087050725477759057, -0.0037727340439829347, -0.012901798319361227, -0.033053795814988754, 0.013156170331620811, -0.015533054933596691, 0.14807112835842184, 0.018761372187400604, 0.080806137493454139, 0.13950623949894433, 0.10703484463450816, -0.26278941929435246, -0.36774772058649907, 0.59625762091850687, -0.36653426739302131, -0.22708969966074771, 0.036289590762181079, 0.069451504639383266, 0.0045810353743422905, 0.014518052118829408, -0.0096918161659651156, -0.0042367342426148941, 0.0017729650965627938, -0.020526538310240098, -0.0097893156767009607, 0.095086485134992974, 0.062351999721190751, -0.01672722282524514, 0.065199506553837216, -0.012727342518659639, 0.021222695138497675, 0.025123508323013449, -0.0086563728711387469, -0.048279776827017555, -0.025594278576058595, 0.04363252969520693, 0.00089194024953267417, -0.021159112170607019, 0.0057838899052937813, 0.042160123131606046, -0.16398523681965352, 0.15410637679425865, 0.11496619447060411, -0.23773808023122042, 0.30878981646904247, -0.091953539825333053, -0.15112210885707561, -0.058067983872415838, -0.016396390506166002, 0.039005803210705077, 0.023932611476630637, -0.006684196401155907, -0.0027208109883995707, 0.0054340236987462662, -0.016167541663812044, 0.0023809918522598732, 0.034708386400909043, -0.0072963544312100569, 0.023224778901069786, 0.10422569904352962, -0.049189051179719279, 0.023151697315084544, 0.0013352981567299367, 0.0065098329229194736, 0.02834601701232067, 0.0067710516924982707, -0.019457415538283694, -0.024637199092487305, 0.019801411975277225, 0.0047456570119127489, -0.015775170314137797, -0.03906718120738624, -0.070536723266456214, 0.063810892086946833, 0.18571376182250354, -0.046218458927843417, 0.055947445070994993, 0.1227003939427266, -0.037984234337547038, -0.12185282737260873, 0.01739900427309727, 0.032960560091392817, 0.056759870222018227, -0.018256051683799891, -0.019427096259341584, 0.010405151723131587, 0.023046337829044038, -0.0032289712715403299, -0.013552286601009622, 0.059732269364482463, 0.026773628633603028, -0.070600771528734632, 0.008835607182220432, -0.096772033742842878, -0.043908100663758234, 0.0023099396961402566, 0.0026341276434262424, 0.0036451640963978101, 0.004808443274984581, 0.020866890117519587, 0.0038698355762962896, -0.019862790828204292, -0.026548451118935872, 0.0080019035626017337, -0.0706384180648158, 0.069515034511643051, 0.070745415327885364, 0.15128619138644453, 0.083803307026608292, 0.046016401957029593, -0.21064351401911435, -0.10796975372596025, -0.034201557920816583, 0.01431321541581594, -0.012363627500514576, -0.010919865739753575, 0.017870360905994497, 0.017396471253638378, -0.02007942542731345, 0.01046579838246435, 0.032566171451921762, -0.05799439311344054, -0.056177295075322847, 0.034333260533212451, -0.130338354535672, 0.12409557570816747, 0.027546636940189641, -0.032798990212655137, 0.083077034868877978, 0.026168992315517999, -0.018046490267575984, 0.005207681729782772, 0.016259252101333162, -0.009045316783825131, -0.0036967486190426836, -0.028592812647173869, -0.0017423598025251908, 0.043534834876013409, 0.12793975070954686, -0.030070570566209584, 0.068857003169531444, -0.017643635220692303, 0.015058032290063167, -0.032968975479440069, -0.0032509456626426633, 0.015788972597383943, 0.019251438395354474, -0.0090248536114061126, -0.015272751179096435, -0.0082992125422014836, 0.059575532325244533, -0.021134287621881149, -0.075932042745803632, 0.10041080102702246, 0.19730813245926007, -0.094203305677599283, 0.12177128356890216, -0.098130362913164085, -0.0086063527817301141, -0.014980171338571052, -0.091992925322100572, 0.17587869380621152, 0.092173327221095724, 0.021897165586045479, 0.0060375152678496429, -0.02386519187296254, -0.012395399307185896, -0.022849248448590303, -0.048321020114327692, -0.0056779923265118269, -0.084303725772385735, 0.1437989861779114, 0.017114978073907855, -0.023674926554513819, -0.014433128838537597, 0.012240235365543381, -0.066419574069182918, -0.016196236751474646, -0.0049338006978182185, -0.0084001340729884851, 0.0022066601859378697, 0.067545520465919467, -0.1431573605324398, -0.11335071558105927, 0.19200701617075636, 0.047005347129143109, -0.18040165272695829, 0.10462001130261284, -0.2952063874323646, 0.45504702044764878, 0.91968716746889512, 0.075629282019523555, -0.17521945191954877, -0.14741145606308895, 0.15842957449053352, 0.033316423902485841, -0.046539635523190549, -0.029239764799856394, 0.028256857448446561, 0.0047468618190423398, -0.025217898898983061, -0.034399558837173737, 0.11937846922953181, 0.050544258155460854, 0.1394578516494776, -0.10661165850403005, -0.027700064466883865, -0.061719692477810832, -0.00076880260958129979, -0.036686593450547596, 0.020834261928534131, 0.042596098823283357, -0.0052374024484263546, -0.17900936499027453, -0.055817128538076279, 0.1556206112582631, -0.19930768818407885, 0.17374602795438115, 1.0450941798657116, 0.011417685121456564, 0.77470195486186011, 0.069649876974948843, -0.36754664039982721, 0.64326753390430735, -0.24917063867041553, -0.039540586335773169, 0.14567847379973992, 0.065981740018920196, -0.055644084876628218, -0.07712597645572139, -0.00073674462317611055, 0.015701976416950342, -0.032298308999268784, 0.060154323470106401, -0.12792055300624763, 0.058681894939637193, -0.022172258797336572, 0.105971385474807, -0.037749743592242278, -0.053672308552617946, -0.12615835751184992, -0.059356927461923663, 0.012406017329547145, 0.068036947599653702, 0.017107991149299284, -0.041008154934057878, -0.13447215443394772, -0.26816260026693878, 0.18368617361305564, 0.46894237369717728, -0.62872636200456533, 0.60686918965331993, -0.15113295529591578, 0.036743730082780514, 1.463608458760923, 0.12904495334293731, -0.074417004563971351, -0.62003762386212891, 0.017408469796284774, 0.1272883333159911, -0.094937931062917147, -0.057791633379354237, 0.037154193600070003, 0.071495661604358593, -0.010724993842566222, -0.14811752790953939, 0.067869387966970118, -0.094974334151300438, 0.31269363662275046, 0.01610138969469304, 0.011271546890921311, -0.0416046739825933, -0.072156922858987255, -0.050238396740841143, -0.013143422355014638, 0.17181320075266959, -0.00078645201341033847, -0.2587885502662679, 0.0028564673365731054, 0.26469201284371935, -0.12782476119556443, -0.59904899563122194, 1.130668752319608, 0.4201199738298067, 1.5043263157657745, 1.0547694551235969, -1.1141651843864939, 0.90050578840564954, -0.50203282214727829, -0.35825626720186354, -0.34045699686317654, 0.034526357222664658, 0.16469615783949973, -0.10773227995476188, -0.020438609938140839, 0.15231382612170388, 0.0040452326400209297, -0.067443262257115694, -0.120243611846639, 0.072438990225412794, -0.27284169211913811, 0.31416612358708496, 0.23288712244734694, 0.17097210848547523, -0.068392282535087695, -0.1884297381053561, 0.026792869524674363, 0.068621242214500999, -0.028891896132398506, 0.01346687181869554, 0.18273076418691644, -0.039027219909814823, -0.68303147161529132, 0.28364033101911773, 0.20446853801739462, 1.3830149751777017, 1.3682011726001495, 1.0314907737559074, -0.55787436167296611, -0.99998707389668362, 0.11752989424299397, -0.85317721061416762, -0.39240115086507799, 0.26946209593128989, 0.13057930804766765, -0.0048298736467919578, -0.0066491170121302217, -0.023829771073630532, -0.072899207479286424, -0.29873174456544577, 0.027468417557455698, -0.26817287178186705, 0.21489820671404794, 0.048277840380711166, 0.23588030048916517, -0.0035092112107956196, -0.076854124253010928, 0.057907844320372041, -0.011517786678353752, 0.045146091846170348, -0.2624959582084973, -0.26585764990704652, 0.54510315877084026, 0.23155959034339049, -0.44086751164178861, -0.68911282851752897, 0.75019268688655205, 1.6899938281650098, 1.2795505844643045, 1.1758121586434311, -0.068510051090992494, 0.24494470763482898, -0.13955938857309558, 0.26258387187077409, -0.31004041285864437, -0.060946718066343622, 0.20150153854805294, -0.21163882649779717, -0.23457709892354675, 0.14508578389537688, -0.22705801052414679, -0.2857183558836735, -0.19596817846267664, 0.13015393203001066, -0.32999646269142935, -0.19528135326362211, 0.36806933922229168, 0.044374683774465004, -0.17789090203786717, -0.13954351571969323, 0.24064979253244634, 0.047498243759330878, -0.33137912963317673, 0.19951473048424839, 0.64027404819684597, 0.25055297311154973, 0.13414323846201037, 1.4561178571660989, 0.50660029824831798, 0.80695903119669876], "height": 32, "width": 32, "top": {"real": [17750.107020847583, 6210.8656159662787, 10017.081922566327, 3457.5995074393513, 265.850218137708, -70.207684810343054, -8.882891147089472, 24.589756503786479, -11.140592103834141, -1.46746599378154, -0.34121813437428317, 1.0812077419156148, 1.7163290532419633, -4.5666146304635689, -2.6846394014203829, -2.2738643666223468, 1.917575640981823, -2.2738643666223641, -2.6846394014204122, -4.5666146304635395, 1.7163290532419744, 1.0812077419156028, -0.34121813437427356, -1.4674659937815784, -11.140592103834141, 24.589756503786372, -8.8828911470896639, -70.20768481034348, 265.8502181377076, 3457.5995074393486, 10017.081922566325, 6210.8656159662851, 23430.31298691556, -4993.9863865787011, -867.92707841112326, 991.9702295167051, 152.49743471969379, -30.096228916888641, 45.849661798659078, -31.467588946689009, 16.960718546454469, 4.334366206560186, -0.42583907131000093, 0.47860527443285728, -1.8738744001784808, 0.5336707597630499, -3.8841941455247579, -0.75407824565303971, 0.55147195392985549, -0.79345242739339283, -5.9026769793174889, -1.384293215600769, -4.1447942145798233, -2.9130522743117342, -0.87120822186396762, -1.7340023689981476, 8.0221119849442548, -9.1034799758133218, 67.038874839272566, 44.112472832699858, 441.79945796279685, 3245.5044122716786, 3861.5097933424754, -5655.7687421695537, -15889.920855591003, -2792.0170209163471, 3795.2577126849956, 553.88440852995757, -85.905664770254234, 195.26564512659746, 24.791859336013662, 19.637521614531401, -6.1680089153590449, -5.0872002952886746, -0.21502707506361318, 1.4832779566036471, 1.4444541156654906, -2.278921340791467, 0.30255454236990775, -2.3473646939882187, -0.02820677859184775, 1.2126049622480428, -0.78557937342030726, -2.9245439997897638, 1.236962081731878, 1.1827041713801927, -1.6114665304853286, -3.5864367039316236, -14.698165954413065, 20.29282209316159, 43.732891706759197, 111.68572777077691, -209.48016621032312, -187.12585183373633, 2586.3148393820152, 1688.0815183011628, 4748.2203910254093, 1775.7319020900284, 965.53459113591498, -592.89741380746352, 99.787886354888371, 95.459765036004825, 59.296635518443523, -7.6987106331430111, -0.0036180413118668542, -0.43290786252066471, -1.598202682681902, 0.11359872229615781, -0.47654646487314006, 0.052920121324431191, -1.1931569571302629, 1.6263047283546588, 1.1208282513608689, -1.916766911729648, -1.3790558562903064, -0.15528235781115707, 0.33989539819890602, -1.1170639390757935, -0.54019671791116586, 2.4039660385463359, 7.3342674422729397, -18.965458306286095, 0.037731231015075362, 5.7888823722245277, -226.40133474887324, -362.07505408372458, 1068.1182163210274, 1839.5542246985658, 1536.2207931138801, 1719.2247320855465, 351.77395499971544, -218.00810259302051, -154.06814929171503, 76.932404651951941, 11.83330919947551, 3.0799843602247936, 2.1514657945922595, -2.1600173905506583, -0.81329272919310069, -0.057640935991259419, 0.81743781787660708, 0.47007834357539929, 0.82578474649026179, -0.485121980563267, -0.44118347877967812, 0.85372587191033056, -0.96505091801449949, -2.5001329354507269, -0.28227610008090631, 1.5082121759974054, 0.73323799037869442, 0.47326064521443867, -2.2071915479661444, 2.1017967421924912, -14.788373732218277, -10.376064867674437, -110.67917899227207, -400.67251265606825, 13.375377575225066, 2228.271084047331, 1272.0006323499422, 485.30359922269918, -150.18264185504094, -215.33992676901764, -8.8694467345951651, 40.216044199991913, 0.62739683840663496, -1.763160286154926, 0.88761506214793773, 2.2829883187502253, -0.37542737915458396, -0.082924323502338942, -0.2690018095258665, 0.066494953017985289, -1.2937095742683771, 0.90148019848405447, 1.049315066873844, -1.5500155662664761, -1.5096304544124948, -0.76725462527686761, -0.45464639460382006, 0.03503878860161394, 0.85393872935961013, 0.89172763012310696, -0.12413241898971583, -7.4057881467564117, -5.5362761946156125, 1.231534306746803, 18.635135255424384, 26.03667534143711, -87.772888824715437, -225.72012452523001, 395.63478219787083, 378.71701677039351, 36.703056312438811, -91.050394397947159, -3.4611377608751717, 16.691588629055683, -0.030826600791664327, -3.1969308796412692, 1.3816292236713215, 0.42488100600164197, 0.084713763632908196, -0.28535803119354425, 0.0058483120842748731, -0.72249565507011926, -0.80262162213105104, -0.11947358836291604, 0.14251583505573939, -0.26639764768858409, -0.59599833623749066, -0.56289480519068613, -0.092201911820288446, 0.77223819072116529, 0.55522491064190738, -0.14411590791637013, -2.1806910872492855, 1.0588352659849944, -1.722617263430356, -9.0167440767821265, -38.483896322518888, -22.599080819691828, 25.541139705127691, 375.57983717352784, 181.47383294778024, 149.67885098073037, -4.1008782517976057, -3.7383475143657869, 18.721866312633942, -1.1602382001305473, -5.8391540939420112, -0.21088907624048864, -0.29870055099180226, 0.95331837446181378, 0.060424395643205529, -0.42858353450158043, -0.37036207729837234, -0.30369644506457261, -0.75652176347658984, -0.44545713631522388, 0.71326845215283508, 0.071761808928956311, -0.63829698193917794, -0.45235505563855871, -0.080159340308425753, -0.15881914095274483, 0.37773826228539636, -0.097086363313521495, -0.86706044290909523, -1.3342582973973045, -0.61792777020441481, 0.16957438507140327, 10.565905071599056, 15.380533375221102, 10.500239194747859, 19.996280010030237, 18.60272375866586, 4.0348314717803326, 11.168149094486143, 14.395513359675057, -3.576502394652771, -3.8760628092228431, 2.0358977627782506, -1.6580642648165547, -0.13166485139511078, 0.2591565726162201, -0.22701969228044258, -0.29312645623416966, 0.099969472199353857, -0.39101661164543228, -0.042546772814251566, 0.29433076927773894, -0.50935384874181677, 0.090487726861141621, -0.092143985882572388, -0.58253529440268437, -0.41671432281001702, 0.14644815319126417, 0.2723500443275419, -0.039709576626662156, -0.41013369223337515, 1.1756350920235048, 1.5259648698081887, -4.4516278662013118, -8.128003687892722, 9.6041080205516938, 20.813833475064627, 12.963050463991683, 7.8870012209429445, 15.407970089260955, 5.3783096889109112, -3.1350030536179672, 1.9101786820520572, 1.6649545288728531, -2.1421648847872561, 0.35230920085802159, 0.37729718445187704, 0.020823929394603979, 0.087139744957248078, -0.10976037032555666, -0.56986899401739111, -0.44020908152219312, 0.69190564661152409, 0.20725962387572397, 0.64042787368671816, 0.71867283099896739, 0.16183428074791173, -0.17223426808716463, -0.17676233483878773, -0.27594664738318686, 0.095238357444986099, 0.020726971935864241, -0.11463597101365795, -0.32333311610982368, 1.2888258236383847, 3.0237538834517612, -0.12189842989060286, -4.5002056479068653, -2.8307120360813536, 9.0325055742836575, -5.0836871456416608, -6.4226783361185662, -2.3832064029645541, 2.1204380507892151, -2.1797923249754283, -2.1316882770768921, 1.4364513793862836, 0.82847417311036753, -0.54705107713687473, -0.053118808001263999, 0.10349330090189141, -0.027586170457780668, -0.27866115184949369, -0.4744187891270259, 0.32835706009268834, 1.8712198771790323, 0.44937016780231559, 0.8488760221500451, 0.35191803530681381, 0.28758243475250972, -0.6426201464805843, -0.40597353100929789, -0.22244529364246871, -0.071869502235051647, 0.33021466122988641, 0.24942995069051274, -0.5214330526633183, -0.093662554700336673, 1.1233733739102922, 1.9003965613191238, 0.07120934753760759, -8.6549236119446178, -1.8427620371931031, -0.52295775466804406, 3.4092952720370122, -0.65003177967723336, -0.1840159964371442, 1.8831870902689825, 0.7983111628921199, -0.39925172982701851, 0.051205177403033819, 0.2598923527237304, -0.18824332097783142, 0.13750859145926242, 0.78226830313090512, -1.2572785001190825, 1.3355363328556911, 1.2220896356262305, 1.6491421999762639, 1.7622421170936717, 0.96487439834416744, 0.79128381470112286, -0.26979590260649117, -0.37236776211072153, -0.65318812250908376, -0.13360581679076858, -0.1817007191765827, 0.4244891020544192, 1.0695325895656689, -0.49200291489349485, -0.59959066580988751, -2.4610304070506381, -1.7528842523917616, 5.4092974592427554, 11.628502736704702, 13.318871722095572, 5.3844201289172418, 0.12839725172390623, -0.91462159985003533, -0.4207633568928546, 0.73289618436063231, 0.98164545606019249, -0.31879685258307294, -0.45467275635234927, 0.23929181220408463, 0.73553732103708702, -0.78050285180289647, -0.20614377250205779, -0.66147354811396541, 0.84722544756444751, 3.1542100241342745, -1.0454812204468338, 0.49112272834774079, 0.89780094097137775, -0.80686344960242362, 0.26166207514350009, 0.34234964401672252, -0.27361449647202374, 0.078400343257288022, 0.27130715748636108, -0.26706243773596255, 0.92879560872659372, -1.3686907450256593, 0.37153092475944227, -0.0077060915391855903, -2.1487270106054996, 8.2423115107991176, 10.211121947580189, 11.004080199445575, 5.512476255300478, 0.43541395104991965, 1.3138259234516578, 0.47697687525002369, -0.37437690637871951, 0.49170325146724819, 0.53092933204425685, 0.78977231856458541, 0.054935970034302541, 1.0907299238052364, -1.4537313504920895, 1.6637415256529431, 2.3271742206176285, 1.6752406948773604, 6.0227778178364675, -2.5693140331095292, -1.3449083811022362, -0.063540759326359386, 1.2007479971921078, 0.26947504329539196, -0.043317616642408706, -0.090092404474960294, -0.56149356452164179, -0.29757479536733583, -0.82531287850485402, 0.52880973840083989, -0.51042849628299192, 1.1866854822445434, 7.7130480762288167, 17.626109919920186, 26.29568713002228, 9.1514979514586194, -1.521929375229722, 0.85364359123771261, -0.365663824522519, -1.1590870007863963, 0.48082630002696519, 0.93059067445638155, 1.1024013035834523, 0.63841829574218301, 1.5659590764611497, 0.20674310541385973, 1.0081259505833966, 1.4865008399227626, -1.6031905208728288, 5.7718835709767413, 0.2720344299284404, 3.2411570780563994, -1.0481514144590243, 2.2338214812745694, 3.1628749936702945, 1.7743358086362808, -0.03099281849361625, -0.41442882712954171, -0.31459130661253104, -1.7323756021672925, 0.19363530105298493, 0.61123212270521721, 4.8260404786904152, 6.7814469674820135, 4.974801806376985, 19.935738411546794, 15.127019575257112, 16.951752765718286, 6.976320734749117, 0.83204582564536733, -0.22090351652873133, -0.48187393778382737, -0.72216259413318939, 0.68337514879437744, 1.4934027691778771, 1.3836806937153638, -0.089108362899698895, 0.57910533764344185, 3.1507795425198797, 1.9139554692796912, 4.1423990225960203, -5.0893580330573318, 4.6545320074577088, -1.3079026033800338, 1.2583901568403451, 1.3087920196200225, 1.8018391908724931, 1.3933980798910508, 0.67604933114803067, -0.49475104816844406, -1.4152366315362725, -0.64580662533466104, -1.4537653735214031, 3.1014323943375786, 7.6843087454063586, 10.909142678924193, 12.359496229531475, 15.280082473470735, 24.669023160454167, 12.739042274822818, 1.7204760023261643, 0.52929508548715176, 2.1482192005482723, -1.481161898887601, -0.79404079710432474, 0.46540881956838565, 0.92624711854277975, 1.0380549761479756, 2.3168267428989626, 1.5788387366050483, -0.83766970523004924, 4.7505720502223729, -2.6570294115495297, 6.7612131323955511, -2.6570294115495336, 4.7505720502223783, -0.83766970523005702, 1.5788387366050476, 2.3168267428989657, 1.0380549761479752, 0.9262471185427773, 0.46540881956838565, -0.79404079710432862, -1.4811618988876001, 2.148219200548279, 0.52929508548715154, 1.7204760023261518, 12.739042274822809, 24.669023160454188, 19.935738411547018, 12.35949622953156, 10.909142678924129, 7.6843087454063417, 3.101432394337563, -1.4537653735213987, -0.64580662533467093, -1.4152366315362839, -0.49475104816844562, 0.67604933114803401, 1.3933980798910552, 1.8018391908724987, 1.3087920196200247, 1.2583901568403484, -1.3079026033800301, 4.6545320074577061, -5.0893580330573203, 4.1423990225960088, 1.913955469279683, 3.1507795425198943, 0.57910533764344618, -0.089108362899694094, 1.3836806937153729, 1.4934027691778826, 0.68337514879437378, -0.72216259413318618, -0.4818739377838111, -0.22090351652873805, 0.83204582564534302, 6.976320734749093, 16.951752765718314, 15.127019575257133, 17.626109919920285, 4.9748018063769956, 6.7814469674819922, 4.8260404786903921, 0.61123212270520444, 0.19363530105297699, -1.7323756021672947, -0.31459130661253587, -0.41442882712954437, -0.030992818493609464, 1.7743358086362877, 3.162874993670298, 2.2338214812745694, -1.0481514144590143, 3.241157078056395, 0.27203442992843363, 5.7718835709767413, -1.6031905208728274, 1.4865008399227673, 1.0081259505833953, 0.20674310541386615, 1.5659590764611582, 0.63841829574219022, 1.1024013035834526, 0.93059067445638244, 0.48082630002696636, -1.1590870007863949, -0.36566382452250978, 0.85364359123770939, -1.5219293752297598, 9.1514979514585963, 26.295687130022237, 8.2423115107992349, 7.7130480762288673, 1.1866854822445243, -0.51042849628302267, 0.52880973840082446, -0.82531287850485058, -0.297574795367338, -0.56149356452165111, -0.090092404474959795, -0.043317616642408116, 0.26947504329539745, 1.2007479971921178, -0.063540759326354057, -1.3449083811022313, -2.5693140331095257, 6.0227778178364719, 1.675240694877361, 2.3271742206176271, 1.6637415256529449, -1.4537313504920832, 1.0907299238052377, 0.054935970034312817, 0.78977231856459063, 0.53092933204426229, 0.49170325146725052, -0.37437690637871607, 0.4769768752500283, 1.3138259234516545, 0.43541395104990577, 5.5124762553004887, 11.004080199445557, 10.211121947580168, 11.628502736704604, -2.1487270106054663, -0.0077060915391693741, 0.37153092475944555, -1.3686907450256616, 0.92879560872659006, -0.26706243773595856, 0.27130715748636525, 0.078400343257288438, -0.27361449647201985, 0.34234964401672396, 0.26166207514350365, -0.80686344960242218, 0.8978009409713853, 0.49112272834772974, -1.045481220446846, 3.1542100241342701, 0.84722544756444151, -0.66147354811395898, -0.20614377250205829, -0.7805028518028938, 0.73553732103708935, 0.23929181220408727, -0.45467275635235105, -0.31879685258307255, 0.98164545606019549, 0.73289618436063131, -0.42076335689285432, -0.91462159985003766, 0.12839725172388108, 5.3844201289172204, 13.318871722095585, -1.8427620371929223, 5.4092974592426515, -1.7528842523917609, -2.4610304070506088, -0.59959066580987819, -0.49200291489348996, 1.0695325895656664, 0.42448910205441354, -0.18170071917658262, -0.13360581679076522, -0.65318812250908243, -0.3723677621107227, -0.26979590260649017, 0.7912838147011273, 0.96487439834416677, 1.7622421170936722, 1.6491421999762705, 1.2220896356262196, 1.3355363328556904, -1.2572785001190749, 0.78226830313090379, 0.13750859145926628, -0.18824332097783478, 0.25989235272372874, 0.051205177403034013, -0.39925172982701423, 0.79831116289212645, 1.8831870902689678, -0.18401599643714045, -0.65003177967722925, 3.4092952720370082, -0.52295775466808114, -5.0836871456417914, -8.654923611944545, 0.071209347537571008, 1.9003965613191072, 1.1233733739103013, -0.093662554700343362, -0.52143305266331286, 0.24942995069051838, 0.33021466122988513, -0.071869502235049398, -0.22244529364246549, -0.40597353100929245, -0.64262014648058341, 0.28758243475251449, 0.35191803530680799, 0.84887602215004632, 0.44937016780231587, 1.8712198771790289, 0.32835706009269156, -0.47441878912702484, -0.27866115184949264, -0.027586170457779294, 0.1034933009018901, -0.053118808001266206, -0.54705107713687218, 0.82847417311037153, 1.4364513793862805, -2.1316882770768819, -2.1797923249754287, 2.1204380507892009, -2.3832064029645679, -6.4226783361184081, 7.8870012209433504, 9.0325055742834763, -2.8307120360814029, -4.5002056479068191, -0.12189842989059124, 3.0237538834517572, 1.2888258236383801, -0.32333311610982529, -0.11463597101365604, 0.020726971935865871, 0.09523835744498739, -0.27594664738318975, -0.17676233483878293, -0.17223426808716019, 0.16183428074791517, 0.71867283099897139, 0.64042787368672205, 0.20725962387572536, 0.69190564661152509, -0.44020908152219435, -0.56986899401739244, -0.10976037032555588, 0.087139744957248716, 0.020823929394600361, 0.37729718445187588, 0.35230920085802736, -2.142164884787245, 1.6649545288728365, 1.9101786820520659, -3.1350030536179658, 5.3783096889108455, 15.40797008926083, 18.60272375866586, 12.963050463991696, 20.813833475064644, 9.604108020551676, -8.1280036878927362, -4.4516278662013278, 1.5259648698081876, 1.1756350920235052, -0.41013369223337515, -0.039709576626662003, 0.27235004432754167, 0.14644815319126442, -0.41671432281001886, -0.58253529440268281, -0.092143985882572915, 0.090487726861142509, -0.50935384874181677, 0.29433076927773738, -0.042546772814252246, -0.39101661164543189, 0.09996947219935487, -0.29312645623417033, -0.22701969228044291, 0.2591565726162216, -0.13166485139511078, -1.6580642648165524, 2.0358977627782484, -3.8760628092228537, -3.5765023946527563, 14.395513359675084, 11.168149094486147, 4.034831471780354, 181.47383294778047, 19.996280010029963, 10.50023919474779, 15.38053337522109, 10.565905071599067, 0.16957438507140921, -0.6179277702044097, -1.3342582973973056, -0.86706044290909146, -0.097086363313522245, 0.37773826228539636, -0.1588191409527461, -0.080159340308420923, -0.45235505563855616, -0.63829698193917661, 0.071761808928958337, 0.71326845215283419, -0.44545713631521988, -0.75652176347658784, -0.30369644506457882, -0.37036207729837373, -0.42858353450157866, 0.060424395643203697, 0.953318374461814, -0.29870055099180354, -0.21088907624048941, -5.8391540939419864, -1.1602382001305538, 18.721866312633882, -3.7383475143657954, -4.100878251797706, 149.67885098073012, 395.63478219787089, 375.57983717352744, 25.541139705127598, -22.599080819691835, -38.483896322518888, -9.0167440767821709, -1.7226172634303247, 1.0588352659849967, -2.1806910872492873, -0.14411590791637696, 0.55522491064190693, 0.77223819072116739, -0.092201911820284366, -0.56289480519067914, -0.59599833623749254, -0.26639764768858487, 0.14251583505573989, -0.11947358836291141, -0.80262162213105359, -0.72249565507012325, 0.0058483120842701027, -0.28535803119354936, 0.084713763632907418, 0.42488100600164375, 1.3816292236713219, -3.1969308796412634, -0.030826600791668796, 16.691588629055648, -3.4611377608751566, -91.050394397947187, 36.703056312438633, 378.71701677039351, 1272.0006323499422, -225.72012452523026, -87.772888824715139, 26.036675341437149, 18.635135255424501, 1.2315343067468421, -5.536276194615616, -7.4057881467564011, -0.12413241898971436, 0.89172763012310496, 0.8539387293596038, 0.035038788601613982, -0.4546463946038164, -0.76725462527685828, -1.509630454412483, -1.550015566266479, 1.0493150668738473, 0.90148019848405436, -1.2937095742683813, 0.066494953017976297, -0.26900180952586789, -0.082924323502345257, -0.37542737915458746, 2.2829883187502267, 0.88761506214793962, -1.7631602861549227, 0.62739683840662275, 40.216044199991877, -8.8694467345952006, -215.33992676901772, -150.18264185504114, 485.3035992226977, 1536.2207931138814, 2228.2710840473319, 13.3753775752247, -400.67251265606865, -110.679178992272, -10.376064867674508, -14.788373732218197, 2.1017967421924983, -2.207191547966143, 0.47326064521441386, 0.73323799037869852, 1.5082121759973919, -0.28227610008090809, -2.5001329354507282, -0.96505091801449672, 0.8537258719103169, -0.44118347877968372, -0.48512198056327205, 0.82578474649025857, 0.47007834357539341, 0.81743781787659597, -0.057640935991271409, -0.81329272919311113, -2.1600173905506721, 2.1514657945922622, 3.0799843602248034, 11.833309199475522, 76.932404651952027, -154.06814929171483, -218.00810259302057, 351.77395499971459, 1719.2247320855445, 4748.2203910254084, 1839.554224698564, 1068.1182163210278, -362.07505408372486, -226.40133474887313, 5.7888823722244789, 0.037731231015098787, -18.965458306285971, 7.3342674422729388, 2.4039660385463657, -0.54019671791118873, -1.1170639390757879, 0.33989539819890108, -0.15528235781113153, -1.3790558562903272, -1.9167669117296613, 1.1208282513608707, 1.6263047283546224, -1.193156957130264, 0.052920121324427458, -0.4765464648731435, 0.11359872229612604, -1.5982026826818876, -0.43290786252068791, -0.0036180413118654109, -7.6987106331430244, 59.296635518443559, 95.459765036004427, 99.787886354888386, -592.89741380746432, 965.53459113591339, 1775.7319020900306, -15889.920855591006, 1688.0815183011607, 2586.3148393820165, -187.12585183373636, -209.48016621032366, 111.68572777077674, 43.732891706759332, 20.292822093161472, -14.698165954413049, -3.5864367039316574, -1.611466530485326, 1.1827041713801545, 1.2369620817318689, -2.9245439997897797, -0.78557937342029793, 1.2126049622479607, -0.028206778591846705, -2.3473646939882631, 0.30255454236988177, -2.278921340791416, 1.4444541156654738, 1.483277956603626, -0.21502707506361871, -5.0872002952886604, -6.168008915359052, 19.637521614531632, 24.791859336013527, 195.26564512659732, -85.905664770254006, 553.88440852995757, 3795.2577126849942, -2792.0170209163462, 23430.312986915553, -5655.7687421695537, 3861.5097933424781, 3245.5044122716772, 441.79945796279628, 44.112472832700512, 67.038874839272665, -9.1034799758128742, 8.0221119849442406, -1.7340023689979776, -0.87120822186400182, -2.9130522743116929, -4.1447942145798446, -1.3842932156007377, -5.9026769793175413, -0.7934524273933764, 0.55147195392985315, -0.75407824565312243, -3.8841941455247726, 0.53367075976302614, -1.8738744001784735, 0.47860527443281381, -0.42583907130996707, 4.3343662065600954, 16.960718546454476, -31.467588946689279, 45.849661798659248, -30.096228916889014, 152.4974347196931, 991.97022951670215, -867.92707841112428, -4993.9863865787029], "imag": [0.0, -2277.9356791272003, -4760.5440393641002, -3362.7502157108247, -452.6326365977277, 12.902525613472081, -66.324591493039264, 21.395810280639775, 12.049924797684376, -6.4322319301924828, -2.6401787430633603, -2.8941924141450355, -1.0996569736184965, -2.1837600734887777, 0.87308165169681873, -0.4143473704962265, 0.0, 0.41434737049617953, -0.87308165169672458, 2.1837600734887985, 1.0996569736184816, 2.8941924141450213, 2.6401787430632839, 6.4322319301924944, -12.049924797684376, -21.395810280639751, 66.32459149303935, -12.902525613472182, 452.63263659772764, 3362.7502157108197, 4760.5440393641011, 2277.9356791272016, -19170.834011440937, -3837.2504381898657, -1562.2865240696121, -2188.2624437144, -79.985389317630222, -65.73387671644744, -129.03232713044309, -21.574731283989859, 19.219741943874205, -1.3770831994436832, -4.6285618656089973, 1.9356332468609161, 1.914358953204174, -0.39274625271235863, -2.8755632420164114, 1.3910769799127469, 2.1983982122300691, -0.93402216187723952, 1.5955119590990243, 2.5332827335358816, 2.5722659267438073, -2.1433748844126024, 4.7051607105169868, 6.0715674792926295, -11.686157623491814, 7.736622407403881, 65.513717323045199, -68.939947660481394, 105.08267177765576, -524.55810593663, 266.59769699529215, -6565.2013915603111, -12998.804313286138, -7140.8743305075541, -1281.4742181751476, 387.93265499317204, 232.0789193763446, -67.213502175129548, -82.499903393484686, 23.742310304940581, 13.167070698799225, -1.2601921779485123, 0.22176027695648487, -0.82869028871129857, 0.81638586337482588, 0.030361957750051641, -1.8114765667849664, -0.33462283171345614, -1.4063404110660631, 1.8330278968089115, -0.19923279402973898, 2.6325944160953934, 0.74767454410312451, 0.33273215463055067, 0.1242120122414053, 0.13350548988021607, -6.7517210150420688, -24.893726268864828, 68.764375920932309, 203.58844824671789, -52.443377241039173, 589.84188198577237, 814.38538391702434, -3077.5993672243876, -9988.2725309741327, -5000.8757444384873, -213.86760970710961, -153.72085790515706, 376.50689706328353, 10.499004391497454, -22.259718007212538, -1.0587560944123571, 1.3367261976572533, -1.845627367813951, -0.50817475490747221, 2.6646257572641145, 0.72663987403826902, -1.4553261209412862, -1.7486114908592469, -2.2393579162674944, 1.7872556824974239, -0.4779279337995983, 0.90990468751900189, 0.56850457149088529, -0.042260514128558184, -2.0950938638350283, 0.37842412015856491, 2.7868839636562086, -7.1796038797667068, -2.6532766677526052, 42.569555354666804, 86.618323256686423, 262.60583525697274, -570.34935164465662, 953.3853188584327, -2608.4519239039555, 423.76919790180068, -777.49072200787271, -730.15551768721184, 267.77827954485849, 48.782690015655945, -52.084914641107339, -0.26732768536918228, 18.791248861451912, 0.034376983736590996, -4.7733525107569852, 0.2450889317389604, 0.69551673309505335, 0.74689297610787508, 0.36254455317814038, -0.082864627703290866, -0.11088180406499534, -1.9145243520823481, 0.63809937787576865, -0.47609558693404908, 1.1704253870693357, 0.1035348526676573, -0.98268640942924379, -0.66403014974353547, 1.4462870038602533, 3.6812627646537175, -8.6847068298461885, -1.7026644231740631, 110.68919188818327, 26.998515460375977, -65.421023262373453, -744.09540667092062, -46.729307478653496, -1526.2181503116453, -735.50774046161519, 282.18880697804167, -182.37864683063907, -54.052681314235201, 39.603849315110104, 10.295701682078779, 2.359337802307917, -0.63461911007383953, -1.7124190510555728, -0.21593130161202514, 0.76755040332111413, 1.2882353738829533, 0.44968314290067118, 0.25989055470659483, -0.66713893962057369, 0.13682839189736165, -0.36883118795348507, 0.91230383136371929, -0.46996446240256329, 0.30219131214854422, -0.19558932994062075, 0.012615763849470061, 1.762420772813289, 1.9045345605049635, -3.7790359077733253, -11.964107143689581, 4.5083476182968525, 61.96781696120042, -254.89688572081204, 140.68826246910922, -45.053807917761794, 314.74785574475101, -8.8001261891015723, -360.31891091191073, -57.025448777424479, 44.579795546563595, -19.359200628613657, 3.6862817789376683, 8.9206630611283089, 0.18896973496668307, -1.0227323524081069, -0.33838743748909134, 0.43833208591156991, 0.0074294169071561673, 0.50446476014305897, 0.1890893708139402, 0.65133332546862077, -0.77483290182744535, -0.30779007494077709, -0.80427129683389309, 0.26011219704759059, 0.23461588343644721, -0.056507066939534729, -0.44513357712322738, 0.64667209976975049, 1.5545017016282721, -1.5250850436867627, -9.5917712473092696, 13.475211488788412, 28.837963587037802, -21.930976284520067, -337.61540934331555, -197.27075321520147, -65.634374404237917, -45.716174186409042, 38.688397244662028, -10.636065931319841, -36.620593610799475, 11.815798616251305, 9.0326984830007095, -3.0325028906240825, -0.075862266735019127, 0.19104942819406223, 0.073726425569811863, 0.18311714867084192, 0.63036999033226249, -0.098626184449280729, 0.095619359707605045, 0.13744024889930673, -0.097089152898782338, -0.81683568275166918, 0.54386217024143169, 0.041501479859297266, 0.43031399471588661, 0.25470104821283102, 0.17326636360277617, 0.46079094226119127, -0.16530865667454481, -0.86391702025516282, -4.9611285532463629, -12.372184270457062, 11.516126179313906, 3.2997425925449804, 2.0570987009733872, 24.878800155519915, 44.041788084316849, 28.446155418037641, -32.190550839656126, -17.692027714058376, 11.789513877445387, 1.9317197634824401, -3.4507560855810171, 0.36626836712779354, 0.48074766055651602, 0.18044316667269442, -0.28240974426973131, -0.16613612027343372, 0.029095718259237995, 0.24450075339721422, -0.091330331747513635, 0.097634295479075664, -0.36806371804963284, 0.16614591312540511, -0.77395352504544712, -0.3264657443238585, 0.015061212014132833, 0.29085100346092985, 0.047770712409502418, 0.14922819541222396, -0.39792121836056477, -0.1921281774085053, 0.83895281997413806, -1.5703737430868945, -8.4367427775736452, 5.6601916863792336, -5.1510032561702124, -22.050605310259861, -1.9854594697574779, -7.5311386138630718, 7.4148341093741337, 8.0952123464619365, -4.3629732153959058, -0.87617222540110651, 1.0539667118389069, -0.69495683349647086, -0.53597221218340763, 0.21612918483653182, 0.16753374785242167, -0.14587335963051168, 0.2984129896925935, 0.75720231383362635, 1.2725586605629424, -0.23788459440434945, -0.43858849762191893, -0.79801858694486538, -0.41846754429723099, -0.46161406526582383, 0.58294801786464678, -0.078939183973068427, 0.31479155143425391, 0.31425797573231379, -0.084400039832661436, -0.69288180199641103, -0.21778449617622889, -0.19533791786522492, -0.22665765093841067, -0.31956792613405932, 6.3352095159571435, 13.304182061875538, 9.0685889178215753, 13.820343051841526, -2.9904131431626517, -6.4644843028741459, 1.4950174056035006, 0.22841081726546483, -1.0443697855595202, -0.36648156701751267, 0.48825453192442753, 0.32511669967392476, -0.74952192547288132, -0.32406742573949188, -0.14589390394899102, 0.81615985313544603, 0.21110858701983537, -0.62135646441400272, -0.27503339195756299, 0.23916839371037824, -1.0440547633891395, -0.40792084439497517, 0.55610134656922072, 0.37457806876122546, 0.17778944632276653, -0.069642667350154328, -0.39526560004564598, 0.69822150614317735, 0.740778974977921, -0.31290252683747921, -2.0096273906801638, -0.75165898775917361, -0.15987777595408992, -2.5516179078608778, -8.9643314564796608, -15.134863386021781, -2.1073345817546003, 0.61887776540407591, -3.0206322635570779, -0.13950893398737912, 0.63280527166495781, -0.16318694006603662, 0.06410706438489891, 0.11411070867346804, -0.2886162817684455, -0.35613651800395313, 0.12924440993745664, 0.38077834891929796, 0.82374133978448782, 0.42592408743381815, -1.4993887570975311, 1.2113263922466315, -0.61651842130949264, -0.92820589700907985, 1.2123284995709924, -0.37394823077908607, -0.062695571319925397, 0.28532586796592635, -0.016166881608109255, -1.1132083022635739, 0.86546596809253507, 2.0760623026361933, 0.51322827606545562, -2.2491343600113889, -2.0374388890600263, 1.1231676638159473, 0.31214242438772172, 2.2815521755579042, -5.3082712246609329, -6.7280927719717081, 0.71404381047788856, 1.0348891792697787, -0.062580928735681385, 0.11553928672311624, 0.19710411955257187, -0.22253376417909151, -0.05143415191355815, -0.61738709055352659, -0.26515618204121805, 1.3879170617746586, 1.9147895569711311, -2.7981144693643607, 1.6351053618238791, 1.2876141106871726, -0.53839815046157202, -0.81616813855627257, -0.55773229035800009, -0.16310476776122179, -1.4828287822897712, 0.19857441865721279, -0.22253902629935005, 0.76149089771667944, 0.38423492006819204, 0.13633825788919049, -0.40617531895308634, -1.9006670978622686, -5.6055295794952666, -3.5201069577585091, -0.70565394755327504, -4.8233879396582013, -3.7118525914264513, 0.87639764763224848, -3.9705679716456896, -0.3066516936670704, 1.2925955094882555, 0.98057929368147423, -0.5185509231518588, 0.068504647941653091, -0.50802141302695303, -0.65025043431383311, 1.6393007306999561, -0.046634622615140439, -0.42095485749155304, 2.8628962346767608, -0.87664183209871016, -0.061380249620325719, -0.67583720185830587, -2.165804182360429, -2.2334105499018988, -1.7315960703588726, 0.023607096501577673, -0.72578349898044647, 0.49340703388266716, -0.17688316963519096, 0.91231104952944475, 1.0401382704665323, 0.40638239361664558, -3.5528586670065754, -3.2389924588328527, -3.336337544628424, -3.3053285670589143, 0.42529364953930032, 1.2623512061316915, -7.9337267372836493, -1.0910503490681764, 0.30591237262721482, 0.79049013658839029, 0.15779582936744363, 0.35530534658354934, -1.0136463777621887, -0.43596379336607172, 1.3766560550226186, 1.2643786263280588, 0.081112249920068616, 0.30101965208538256, -2.9649644331404845, 3.940601481283172, 1.2809925361642973, -2.3740821954043603, -1.8568139214134618, -2.2875531209700912, -2.0068530505793203, -1.6161905886738603, -0.069970559948539499, -0.68452644609555957, 1.1180256373411306, 1.0386505011587546, -0.45802414484767806, -1.8827293053017495, -1.1452872829873149, -9.8282380345909921, -4.8375314158562777, 0.035534088492687292, 5.8735645683280717, -7.2148282372822923, -1.9024855757814658, -3.6729219824305321, -0.04721742712852893, 0.80045523566507781, 0.91856151213872683, 0.18365173124028911, 1.0527130168554166, 0.35231243394437906, -0.96832866520995176, 0.38866324603453256, 1.1024459216780036, -1.6264163010988006, 5.9086248051519128, 1.1522588148293589, -4.0123769411888297, -1.3409560614230069, -2.2450482885496141, -1.484164394238348, 0.58569229342098317, 0.035646451001070931, -1.3588202093923047, 0.006860270440306241, -0.84539036438792836, 0.31580363771161246, -0.13161779536446233, -0.34549511895795371, -1.1291296924180168, -2.1947930228924113, -14.643424663744323, 0.0, 7.8456512467946364, 11.843681310787776, -4.1545820919682566, 0.098834485797376734, 1.2191537903053737, 0.14782839419733967, -0.21707569483499198, 1.3287726303022269, 0.012891904979339022, 0.146292284007916, 0.70147269041246085, 1.7613861146868823, 0.97884915683286033, 1.5860078469473597, 0.15814907094033309, 0.0, -0.158149070940315, -1.5860078469473551, -0.978849156832858, -1.7613861146868832, -0.70147269041248062, -0.14629228400792066, -0.012891904979343359, -1.3287726303022269, 0.21707569483499028, -0.14782839419734384, -1.2191537903053777, -0.098834485797375804, 4.1545820919682752, -11.843681310787789, -7.8456512467946329, -0.035534088492520481, 14.643424663744339, 2.194793022892366, 1.1291296924180338, 0.34549511895794383, 0.13161779536445045, -0.31580363771161418, 0.84539036438793247, -0.0068602704403050632, 1.3588202093923036, -0.035646451001076773, -0.58569229342099216, 1.484164394238344, 2.2450482885496172, 1.3409560614230207, 4.0123769411888466, -1.1522588148293609, -5.9086248051519057, 1.6264163010987922, -1.102445921678016, -0.38866324603453178, 0.96832866520995076, -0.35231243394438311, -1.0527130168554157, -0.18365173124028866, -0.91856151213871651, -0.80045523566507271, 0.047217427128521561, 3.6729219824305099, 1.9024855757814958, 7.21482823728231, -5.8735645683280371, 3.3053285670587775, 4.8375314158562217, 9.828238034591017, 1.1452872829873084, 1.8827293053017478, 0.45802414484767828, -1.0386505011587623, -1.1180256373411397, 0.68452644609555668, 0.069970559948536071, 1.6161905886738546, 2.0068530505793158, 2.2875531209700912, 1.856813921413464, 2.3740821954043709, -1.2809925361643011, -3.9406014812831787, 2.9649644331405014, -0.30101965208537818, -0.081112249920068671, -1.2643786263280541, -1.3766560550226248, 0.43596379336607571, 1.0136463777621971, -0.35530534658355051, -0.15779582936744335, -0.79049013658838307, -0.30591237262721799, 1.0910503490681696, 7.9337267372836298, -1.2623512061316555, -0.42529364953934357, 0.70565394755329391, 3.336337544628484, 3.2389924588327594, 3.5528586670065589, -0.40638239361663953, -1.0401382704665338, -0.91231104952944642, 0.17688316963519082, -0.49340703388266682, 0.72578349898044092, -0.023607096501574151, 1.7315960703588795, 2.2334105499019046, 2.1658041823604322, 0.67583720185832141, 0.061380249620329064, 0.87664183209870872, -2.8628962346767408, 0.4209548574915572, 0.046634622615141327, -1.6393007306999556, 0.65025043431383445, 0.50802141302695369, -0.068504647941647845, 0.51855092315186113, -0.98057929368147057, -1.2925955094882546, 0.30665169366707518, 3.9705679716456634, -0.8763976476322699, 3.7118525914263776, 4.8233879396582111, -0.31214242438772705, 3.5201069577584692, 5.6055295794952924, 1.9006670978622411, 0.40617531895307485, -0.13633825788919066, -0.38423492006819598, -0.76149089771668665, 0.22253902629934941, -0.1985744186572119, 1.4828287822897714, 0.16310476776122626, 0.55773229035800354, 0.8161681385562749, 0.5383981504615758, -1.2876141106871739, -1.6351053618238793, 2.7981144693643794, -1.9147895569711246, -1.3879170617746619, 0.26515618204121727, 0.61738709055352903, 0.051434151913562459, 0.22253376417909698, -0.19710411955257243, -0.11553928672311804, 0.062580928735678887, -1.0348891792697743, -0.71404381047787879, 6.7280927719716743, 5.3082712246609018, -2.2815521755579335, 8.9643314564797052, -1.1231676638158536, 2.0374388890599802, 2.2491343600114058, -0.51322827606545429, -2.0760623026361911, -0.86546596809253695, 1.1132083022635748, 0.016166881608109234, -0.28532586796592802, 0.062695571319927382, 0.37394823077908801, -1.2123284995709924, 0.92820589700907918, 0.61651842130949663, -1.2113263922466326, 1.4993887570975342, -0.42592408743381316, -0.82374133978448039, -0.38077834891929335, -0.12924440993745251, 0.35613651800395341, 0.28861628176844567, -0.11411070867346547, -0.064107064384899673, 0.16318694006603501, -0.63280527166495915, 0.13950893398738348, 3.0206322635570557, -0.61887776540410167, 2.1073345817545266, 15.134863386021769, -9.0685889178216552, 2.5516179078607184, 0.15987777595415892, 0.75165898775916729, 2.0096273906801727, 0.31290252683747422, -0.74077897497792289, -0.69822150614317713, 0.39526560004564548, 0.069642667350159129, -0.1777894463227658, -0.37457806876122451, -0.55610134656922283, 0.40792084439497267, 1.0440547633891391, -0.23916839371037832, 0.2750333919575636, 0.62135646441400449, -0.21110858701983135, -0.81615985313544015, 0.14589390394899499, 0.32406742573949032, 0.74952192547288177, -0.32511669967392071, -0.48825453192442569, 0.36648156701751228, 1.0443697855595189, -0.22841081726546683, -1.4950174056034862, 6.4644843028741308, 2.9904131431626397, -13.820343051841599, 1.9854594697574315, -13.304182061875386, -6.3352095159572039, 0.3195679261341049, 0.22665765093840967, 0.19533791786523994, 0.21778449617622839, 0.69288180199640803, 0.084400039832664031, -0.31425797573231284, -0.31479155143425158, 0.078939183973069621, -0.58294801786464789, 0.46161406526582016, 0.41846754429722727, 0.79801858694486438, 0.43858849762191948, 0.23788459440435036, -1.2725586605629382, -0.75720231383362191, -0.29841298969258917, 0.14587335963051393, -0.16753374785242225, -0.21612918483653495, 0.53597221218340851, 0.69495683349647419, -1.0539667118389073, 0.87617222540110107, 4.3629732153958827, -8.0952123464619206, -7.414834109374171, 7.5311386138631651, -44.041788084316849, 22.050605310259897, 5.1510032561702142, -5.6601916863792514, 8.4367427775736488, 1.5703737430868991, -0.83895281997413862, 0.19212817740850727, 0.39792121836056477, -0.14922819541222082, -0.047770712409500753, -0.2908510034609289, -0.015061212014133111, 0.32646574432385927, 0.77395352504545112, -0.16614591312540627, 0.36806371804963284, -0.097634295479077288, 0.091330331747514912, -0.24450075339721231, -0.029095718259234932, 0.16613612027343472, 0.28240974426972948, -0.18044316667269469, -0.48074766055651602, -0.36626836712779753, 3.4507560855810189, -1.9317197634824488, -11.789513877445398, 17.692027714058376, 32.190550839656112, -28.446155418037584, 65.634374404237818, -24.878800155519698, -2.0570987009734387, -3.2997425925449106, -11.516126179313863, 12.372184270457085, 4.9611285532463585, 0.86391702025516592, 0.16530865667454706, -0.46079094226119116, -0.1732663636027687, -0.25470104821282824, -0.43031399471588622, -0.041501479859297648, -0.54386217024143224, 0.81683568275167029, 0.097089152898783726, -0.13744024889929929, -0.095619359707603865, 0.098626184449285684, -0.63036999033225893, -0.18311714867083492, -0.07372642556981146, -0.19104942819406445, 0.075862266735016198, 3.0325028906240856, -9.0326984830007238, -11.815798616251362, 36.620593610799446, 10.636065931319884, -38.68839724466217, 45.716174186409333, -314.74785574475044, 197.27075321520152, 337.61540934331566, 21.930976284520245, -28.837963587037731, -13.475211488788421, 9.5917712473092944, 1.5250850436867665, -1.5545017016282761, -0.64667209976974838, 0.44513357712322837, 0.056507066939540815, -0.23461588343644407, -0.26011219704759359, 0.8042712968338912, 0.30779007494077798, 0.77483290182744469, -0.65133332546861145, -0.18908937081392915, -0.50446476014305819, -0.0074294169071514489, -0.43833208591156614, 0.3383874374890965, 1.0227323524080993, -0.18896973496668346, -8.9206630611282982, -3.6862817789376785, 19.35920062861366, -44.579795546563581, 57.025448777424586, 360.31891091191056, 8.8001261891013822, 1526.2181503116453, 45.053807917762221, -140.68826246910882, 254.89688572081207, -61.967816961200398, -4.5083476182968329, 11.964107143689583, 3.7790359077733431, -1.9045345605049631, -1.7624207728132955, -0.01261576384946491, 0.19558932994062728, -0.30219131214854072, 0.4699644624025669, -0.91230383136372128, 0.36883118795348824, -0.13682839189736087, 0.66713893962058779, -0.25989055470659439, -0.44968314290067052, -1.2882353738829488, -0.76755040332110314, 0.21593130161202181, 1.7124190510555817, 0.6346191100738432, -2.3593378023079388, -10.295701682078779, -39.603849315110132, 54.052681314235315, 182.3786468306393, -282.18880697804133, 735.50774046161428, -423.76919790180091, 46.72930747865378, 744.09540667092074, 65.421023262373566, -26.998515460375824, -110.68919188818332, 1.7026644231740455, 8.6847068298461902, -3.6812627646537122, -1.4462870038602713, 0.66403014974352648, 0.98268640942926089, -0.10353485266764806, -1.1704253870693322, 0.47609558693405896, -0.63809937787576798, 1.9145243520823476, 0.11088180406502474, 0.082864627703299776, -0.36254455317814699, -0.74689297610786809, -0.69551673309504447, -0.24508893173895085, 4.773352510756971, -0.034376983736587929, -18.791248861451901, 0.26732768536914786, 52.084914641107339, -48.782690015656009, -267.77827954485792, 730.15551768721207, 777.49072200787191, 9988.2725309741272, 2608.4519239039569, -953.38531885843202, 570.34935164465708, -262.60583525697274, -86.61832325668658, -42.569555354666946, 2.6532766677526505, 7.1796038797667014, -2.7868839636562561, -0.37842412015855753, 2.0950938638350354, 0.042260514128576405, -0.56850457149086808, -0.90990468751900622, 0.4779279337996068, -1.7872556824974217, 2.239357916267509, 1.748611490859272, 1.4553261209412682, -0.72663987403827102, -2.6646257572640879, 0.50817475490745234, 1.8456273678139588, -1.3367261976572486, 1.0587560944122603, 22.259718007212488, -10.499004391497355, -376.50689706328382, 153.72085790515695, 213.86760970710961, 5000.8757444384837, 12998.804313286144, 3077.5993672243935, -814.38538391702389, -589.84188198577237, 52.443377241039045, -203.58844824671809, -68.764375920932238, 24.893726268864992, 6.7517210150420599, -0.13350548988021724, -0.12421201224142447, -0.33273215463053935, -0.74767454410310241, -2.6325944160953538, 0.19923279402974362, -1.8330278968089646, 1.4063404110660558, 0.33462283171345686, 1.8114765667849764, -0.030361957750053889, -0.81638586337482577, 0.82869028871127859, -0.22176027695649836, 1.260192177948489, -13.167070698799225, -23.742310304940538, 82.499903393484416, 67.213502175129648, -232.07891937634469, -387.93265499317204, 1281.4742181751474, 7140.8743305075477, 19170.83401144094, 6565.2013915603166, -266.59769699529113, 524.55810593663034, -105.08267177765569, 68.939947660481323, -65.513717323045299, -7.7366224074036598, 11.686157623491809, -6.0715674792925807, -4.7051607105170206, 2.1433748844125784, -2.5722659267438033, -2.5332827335358727, -1.5955119590990279, 0.93402216187716958, -2.1983982122300714, -1.3910769799127845, 2.875563242016451, 0.39274625271232855, -1.9143589532041787, -1.9356332468609148, 4.6285618656089529, 1.3770831994437265, -19.219741943874194, 21.574731283989856, 129.03232713044289, 65.73387671644744, 79.985389317629867, 2188.2624437143977, 1562.286524069611, 3837.2504381898625]}};

/**
 * Numeric Javascript
 * Copyright (C) 2011 by Sébastien Loisel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

"use strict";

var numeric = (typeof exports === "undefined")?(function numeric() {}):(exports);
if(typeof global !== "undefined") { global.numeric = numeric; }

numeric.version = "1.2.6";

// 1. Utility functions
numeric.bench = function bench (f,interval) {
    var t1,t2,n,i;
    if(typeof interval === "undefined") { interval = 15; }
    n = 0.5;
    t1 = new Date();
    while(1) {
        n*=2;
        for(i=n;i>3;i-=4) { f(); f(); f(); f(); }
        while(i>0) { f(); i--; }
        t2 = new Date();
        if(t2-t1 > interval) break;
    }
    for(i=n;i>3;i-=4) { f(); f(); f(); f(); }
    while(i>0) { f(); i--; }
    t2 = new Date();
    return 1000*(3*n-1)/(t2-t1);
}

numeric._myIndexOf = (function _myIndexOf(w) {
    var n = this.length,k;
    for(k=0;k<n;++k) if(this[k]===w) return k;
    return -1;
});
numeric.myIndexOf = (Array.prototype.indexOf)?Array.prototype.indexOf:numeric._myIndexOf;

numeric.Function = Function;
numeric.precision = 4;
numeric.largeArray = 50;

numeric.prettyPrint = function prettyPrint(x) {
    function fmtnum(x) {
        if(x === 0) { return '0'; }
        if(isNaN(x)) { return 'NaN'; }
        if(x<0) { return '-'+fmtnum(-x); }
        if(isFinite(x)) {
            var scale = Math.floor(Math.log(x) / Math.log(10));
            var normalized = x / Math.pow(10,scale);
            var basic = normalized.toPrecision(numeric.precision);
            if(parseFloat(basic) === 10) { scale++; normalized = 1; basic = normalized.toPrecision(numeric.precision); }
            return parseFloat(basic).toString()+'e'+scale.toString();
        }
        return 'Infinity';
    }
    var ret = [];
    function foo(x) {
        var k;
        if(typeof x === "undefined") { ret.push(Array(numeric.precision+8).join(' ')); return false; }
        if(typeof x === "string") { ret.push('"'+x+'"'); return false; }
        if(typeof x === "boolean") { ret.push(x.toString()); return false; }
        if(typeof x === "number") {
            var a = fmtnum(x);
            var b = x.toPrecision(numeric.precision);
            var c = parseFloat(x.toString()).toString();
            var d = [a,b,c,parseFloat(b).toString(),parseFloat(c).toString()];
            for(k=1;k<d.length;k++) { if(d[k].length < a.length) a = d[k]; }
            ret.push(Array(numeric.precision+8-a.length).join(' ')+a);
            return false;
        }
        if(x === null) { ret.push("null"); return false; }
        if(typeof x === "function") { 
            ret.push(x.toString());
            var flag = false;
            for(k in x) { if(x.hasOwnProperty(k)) { 
                if(flag) ret.push(',\n');
                else ret.push('\n{');
                flag = true; 
                ret.push(k); 
                ret.push(': \n'); 
                foo(x[k]); 
            } }
            if(flag) ret.push('}\n');
            return true;
        }
        if(x instanceof Array) {
            if(x.length > numeric.largeArray) { ret.push('...Large Array...'); return true; }
            var flag = false;
            ret.push('[');
            for(k=0;k<x.length;k++) { if(k>0) { ret.push(','); if(flag) ret.push('\n '); } flag = foo(x[k]); }
            ret.push(']');
            return true;
        }
        ret.push('{');
        var flag = false;
        for(k in x) { if(x.hasOwnProperty(k)) { if(flag) ret.push(',\n'); flag = true; ret.push(k); ret.push(': \n'); foo(x[k]); } }
        ret.push('}');
        return true;
    }
    foo(x);
    return ret.join('');
}

numeric.parseDate = function parseDate(d) {
    function foo(d) {
        if(typeof d === 'string') { return Date.parse(d.replace(/-/g,'/')); }
        if(!(d instanceof Array)) { throw new Error("parseDate: parameter must be arrays of strings"); }
        var ret = [],k;
        for(k=0;k<d.length;k++) { ret[k] = foo(d[k]); }
        return ret;
    }
    return foo(d);
}

numeric.parseFloat = function parseFloat_(d) {
    function foo(d) {
        if(typeof d === 'string') { return parseFloat(d); }
        if(!(d instanceof Array)) { throw new Error("parseFloat: parameter must be arrays of strings"); }
        var ret = [],k;
        for(k=0;k<d.length;k++) { ret[k] = foo(d[k]); }
        return ret;
    }
    return foo(d);
}

numeric.parseCSV = function parseCSV(t) {
    var foo = t.split('\n');
    var j,k;
    var ret = [];
    var pat = /(([^'",]*)|('[^']*')|("[^"]*")),/g;
    var patnum = /^\s*(([+-]?[0-9]+(\.[0-9]*)?(e[+-]?[0-9]+)?)|([+-]?[0-9]*(\.[0-9]+)?(e[+-]?[0-9]+)?))\s*$/;
    var stripper = function(n) { return n.substr(0,n.length-1); }
    var count = 0;
    for(k=0;k<foo.length;k++) {
      var bar = (foo[k]+",").match(pat),baz;
      if(bar.length>0) {
          ret[count] = [];
          for(j=0;j<bar.length;j++) {
              baz = stripper(bar[j]);
              if(patnum.test(baz)) { ret[count][j] = parseFloat(baz); }
              else ret[count][j] = baz;
          }
          count++;
      }
    }
    return ret;
}

numeric.toCSV = function toCSV(A) {
    var s = numeric.dim(A);
    var i,j,m,n,row,ret;
    m = s[0];
    n = s[1];
    ret = [];
    for(i=0;i<m;i++) {
        row = [];
        for(j=0;j<m;j++) { row[j] = A[i][j].toString(); }
        ret[i] = row.join(', ');
    }
    return ret.join('\n')+'\n';
}

numeric.getURL = function getURL(url) {
    var client = new XMLHttpRequest();
    client.open("GET",url,false);
    client.send();
    return client;
}

numeric.imageURL = function imageURL(img) {
    function base64(A) {
        var n = A.length, i,x,y,z,p,q,r,s;
        var key = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";
        var ret = "";
        for(i=0;i<n;i+=3) {
            x = A[i];
            y = A[i+1];
            z = A[i+2];
            p = x >> 2;
            q = ((x & 3) << 4) + (y >> 4);
            r = ((y & 15) << 2) + (z >> 6);
            s = z & 63;
            if(i+1>=n) { r = s = 64; }
            else if(i+2>=n) { s = 64; }
            ret += key.charAt(p) + key.charAt(q) + key.charAt(r) + key.charAt(s);
            }
        return ret;
    }
    function crc32Array (a,from,to) {
        if(typeof from === "undefined") { from = 0; }
        if(typeof to === "undefined") { to = a.length; }
        var table = [0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA, 0x076DC419, 0x706AF48F, 0xE963A535, 0x9E6495A3,
                     0x0EDB8832, 0x79DCB8A4, 0xE0D5E91E, 0x97D2D988, 0x09B64C2B, 0x7EB17CBD, 0xE7B82D07, 0x90BF1D91, 
                     0x1DB71064, 0x6AB020F2, 0xF3B97148, 0x84BE41DE, 0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
                     0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC, 0x14015C4F, 0x63066CD9, 0xFA0F3D63, 0x8D080DF5, 
                     0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172, 0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B, 
                     0x35B5A8FA, 0x42B2986C, 0xDBBBC9D6, 0xACBCF940, 0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59, 
                     0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116, 0x21B4F4B5, 0x56B3C423, 0xCFBA9599, 0xB8BDA50F,
                     0x2802B89E, 0x5F058808, 0xC60CD9B2, 0xB10BE924, 0x2F6F7C87, 0x58684C11, 0xC1611DAB, 0xB6662D3D,
                     0x76DC4190, 0x01DB7106, 0x98D220BC, 0xEFD5102A, 0x71B18589, 0x06B6B51F, 0x9FBFE4A5, 0xE8B8D433,
                     0x7807C9A2, 0x0F00F934, 0x9609A88E, 0xE10E9818, 0x7F6A0DBB, 0x086D3D2D, 0x91646C97, 0xE6635C01, 
                     0x6B6B51F4, 0x1C6C6162, 0x856530D8, 0xF262004E, 0x6C0695ED, 0x1B01A57B, 0x8208F4C1, 0xF50FC457, 
                     0x65B0D9C6, 0x12B7E950, 0x8BBEB8EA, 0xFCB9887C, 0x62DD1DDF, 0x15DA2D49, 0x8CD37CF3, 0xFBD44C65, 
                     0x4DB26158, 0x3AB551CE, 0xA3BC0074, 0xD4BB30E2, 0x4ADFA541, 0x3DD895D7, 0xA4D1C46D, 0xD3D6F4FB, 
                     0x4369E96A, 0x346ED9FC, 0xAD678846, 0xDA60B8D0, 0x44042D73, 0x33031DE5, 0xAA0A4C5F, 0xDD0D7CC9, 
                     0x5005713C, 0x270241AA, 0xBE0B1010, 0xC90C2086, 0x5768B525, 0x206F85B3, 0xB966D409, 0xCE61E49F, 
                     0x5EDEF90E, 0x29D9C998, 0xB0D09822, 0xC7D7A8B4, 0x59B33D17, 0x2EB40D81, 0xB7BD5C3B, 0xC0BA6CAD, 
                     0xEDB88320, 0x9ABFB3B6, 0x03B6E20C, 0x74B1D29A, 0xEAD54739, 0x9DD277AF, 0x04DB2615, 0x73DC1683, 
                     0xE3630B12, 0x94643B84, 0x0D6D6A3E, 0x7A6A5AA8, 0xE40ECF0B, 0x9309FF9D, 0x0A00AE27, 0x7D079EB1, 
                     0xF00F9344, 0x8708A3D2, 0x1E01F268, 0x6906C2FE, 0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7, 
                     0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC, 0xF9B9DF6F, 0x8EBEEFF9, 0x17B7BE43, 0x60B08ED5, 
                     0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252, 0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B, 
                     0xD80D2BDA, 0xAF0A1B4C, 0x36034AF6, 0x41047A60, 0xDF60EFC3, 0xA867DF55, 0x316E8EEF, 0x4669BE79, 
                     0xCB61B38C, 0xBC66831A, 0x256FD2A0, 0x5268E236, 0xCC0C7795, 0xBB0B4703, 0x220216B9, 0x5505262F, 
                     0xC5BA3BBE, 0xB2BD0B28, 0x2BB45A92, 0x5CB36A04, 0xC2D7FFA7, 0xB5D0CF31, 0x2CD99E8B, 0x5BDEAE1D, 
                     0x9B64C2B0, 0xEC63F226, 0x756AA39C, 0x026D930A, 0x9C0906A9, 0xEB0E363F, 0x72076785, 0x05005713, 
                     0x95BF4A82, 0xE2B87A14, 0x7BB12BAE, 0x0CB61B38, 0x92D28E9B, 0xE5D5BE0D, 0x7CDCEFB7, 0x0BDBDF21, 
                     0x86D3D2D4, 0xF1D4E242, 0x68DDB3F8, 0x1FDA836E, 0x81BE16CD, 0xF6B9265B, 0x6FB077E1, 0x18B74777, 
                     0x88085AE6, 0xFF0F6A70, 0x66063BCA, 0x11010B5C, 0x8F659EFF, 0xF862AE69, 0x616BFFD3, 0x166CCF45, 
                     0xA00AE278, 0xD70DD2EE, 0x4E048354, 0x3903B3C2, 0xA7672661, 0xD06016F7, 0x4969474D, 0x3E6E77DB, 
                     0xAED16A4A, 0xD9D65ADC, 0x40DF0B66, 0x37D83BF0, 0xA9BCAE53, 0xDEBB9EC5, 0x47B2CF7F, 0x30B5FFE9, 
                     0xBDBDF21C, 0xCABAC28A, 0x53B39330, 0x24B4A3A6, 0xBAD03605, 0xCDD70693, 0x54DE5729, 0x23D967BF, 
                     0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94, 0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D];
     
        var crc = -1, y = 0, n = a.length,i;

        for (i = from; i < to; i++) {
            y = (crc ^ a[i]) & 0xFF;
            crc = (crc >>> 8) ^ table[y];
        }
     
        return crc ^ (-1);
    }

    var h = img[0].length, w = img[0][0].length, s1, s2, next,k,length,a,b,i,j,adler32,crc32;
    var stream = [
                  137, 80, 78, 71, 13, 10, 26, 10,                           //  0: PNG signature
                  0,0,0,13,                                                  //  8: IHDR Chunk length
                  73, 72, 68, 82,                                            // 12: "IHDR" 
                  (w >> 24) & 255, (w >> 16) & 255, (w >> 8) & 255, w&255,   // 16: Width
                  (h >> 24) & 255, (h >> 16) & 255, (h >> 8) & 255, h&255,   // 20: Height
                  8,                                                         // 24: bit depth
                  2,                                                         // 25: RGB
                  0,                                                         // 26: deflate
                  0,                                                         // 27: no filter
                  0,                                                         // 28: no interlace
                  -1,-2,-3,-4,                                               // 29: CRC
                  -5,-6,-7,-8,                                               // 33: IDAT Chunk length
                  73, 68, 65, 84,                                            // 37: "IDAT"
                  // RFC 1950 header starts here
                  8,                                                         // 41: RFC1950 CMF
                  29                                                         // 42: RFC1950 FLG
                  ];
    crc32 = crc32Array(stream,12,29);
    stream[29] = (crc32>>24)&255;
    stream[30] = (crc32>>16)&255;
    stream[31] = (crc32>>8)&255;
    stream[32] = (crc32)&255;
    s1 = 1;
    s2 = 0;
    for(i=0;i<h;i++) {
        if(i<h-1) { stream.push(0); }
        else { stream.push(1); }
        a = (3*w+1+(i===0))&255; b = ((3*w+1+(i===0))>>8)&255;
        stream.push(a); stream.push(b);
        stream.push((~a)&255); stream.push((~b)&255);
        if(i===0) stream.push(0);
        for(j=0;j<w;j++) {
            for(k=0;k<3;k++) {
                a = img[k][i][j];
                if(a>255) a = 255;
                else if(a<0) a=0;
                else a = Math.round(a);
                s1 = (s1 + a )%65521;
                s2 = (s2 + s1)%65521;
                stream.push(a);
            }
        }
        stream.push(0);
    }
    adler32 = (s2<<16)+s1;
    stream.push((adler32>>24)&255);
    stream.push((adler32>>16)&255);
    stream.push((adler32>>8)&255);
    stream.push((adler32)&255);
    length = stream.length - 41;
    stream[33] = (length>>24)&255;
    stream[34] = (length>>16)&255;
    stream[35] = (length>>8)&255;
    stream[36] = (length)&255;
    crc32 = crc32Array(stream,37);
    stream.push((crc32>>24)&255);
    stream.push((crc32>>16)&255);
    stream.push((crc32>>8)&255);
    stream.push((crc32)&255);
    stream.push(0);
    stream.push(0);
    stream.push(0);
    stream.push(0);
//    a = stream.length;
    stream.push(73);  // I
    stream.push(69);  // E
    stream.push(78);  // N
    stream.push(68);  // D
    stream.push(174); // CRC1
    stream.push(66);  // CRC2
    stream.push(96);  // CRC3
    stream.push(130); // CRC4
    return 'data:image/png;base64,'+base64(stream);
}

// 2. Linear algebra with Arrays.
numeric._dim = function _dim(x) {
    var ret = [];
    while(typeof x === "object") { ret.push(x.length); x = x[0]; }
    return ret;
}

numeric.dim = function dim(x) {
    var y,z;
    if(typeof x === "object") {
        y = x[0];
        if(typeof y === "object") {
            z = y[0];
            if(typeof z === "object") {
                return numeric._dim(x);
            }
            return [x.length,y.length];
        }
        return [x.length];
    }
    return [];
}

numeric.mapreduce = function mapreduce(body,init) {
    return Function('x','accum','_s','_k',
            'if(typeof accum === "undefined") accum = '+init+';\n'+
            'if(typeof x === "number") { var xi = x; '+body+'; return accum; }\n'+
            'if(typeof _s === "undefined") _s = numeric.dim(x);\n'+
            'if(typeof _k === "undefined") _k = 0;\n'+
            'var _n = _s[_k];\n'+
            'var i,xi;\n'+
            'if(_k < _s.length-1) {\n'+
            '    for(i=_n-1;i>=0;i--) {\n'+
            '        accum = arguments.callee(x[i],accum,_s,_k+1);\n'+
            '    }'+
            '    return accum;\n'+
            '}\n'+
            'for(i=_n-1;i>=1;i-=2) { \n'+
            '    xi = x[i];\n'+
            '    '+body+';\n'+
            '    xi = x[i-1];\n'+
            '    '+body+';\n'+
            '}\n'+
            'if(i === 0) {\n'+
            '    xi = x[i];\n'+
            '    '+body+'\n'+
            '}\n'+
            'return accum;'
            );
}
numeric.mapreduce2 = function mapreduce2(body,setup) {
    return Function('x',
            'var n = x.length;\n'+
            'var i,xi;\n'+setup+';\n'+
            'for(i=n-1;i!==-1;--i) { \n'+
            '    xi = x[i];\n'+
            '    '+body+';\n'+
            '}\n'+
            'return accum;'
            );
}


numeric.same = function same(x,y) {
    var i,n;
    if(!(x instanceof Array) || !(y instanceof Array)) { return false; }
    n = x.length;
    if(n !== y.length) { return false; }
    for(i=0;i<n;i++) {
        if(x[i] === y[i]) { continue; }
        if(typeof x[i] === "object") { if(!same(x[i],y[i])) return false; }
        else { return false; }
    }
    return true;
}

numeric.rep = function rep(s,v,k) {
    if(typeof k === "undefined") { k=0; }
    var n = s[k], ret = Array(n), i;
    if(k === s.length-1) {
        for(i=n-2;i>=0;i-=2) { ret[i+1] = v; ret[i] = v; }
        if(i===-1) { ret[0] = v; }
        return ret;
    }
    for(i=n-1;i>=0;i--) { ret[i] = numeric.rep(s,v,k+1); }
    return ret;
}


numeric.dotMMsmall = function dotMMsmall(x,y) {
    var i,j,k,p,q,r,ret,foo,bar,woo,i0,k0,p0,r0;
    p = x.length; q = y.length; r = y[0].length;
    ret = Array(p);
    for(i=p-1;i>=0;i--) {
        foo = Array(r);
        bar = x[i];
        for(k=r-1;k>=0;k--) {
            woo = bar[q-1]*y[q-1][k];
            for(j=q-2;j>=1;j-=2) {
                i0 = j-1;
                woo += bar[j]*y[j][k] + bar[i0]*y[i0][k];
            }
            if(j===0) { woo += bar[0]*y[0][k]; }
            foo[k] = woo;
        }
        ret[i] = foo;
    }
    return ret;
}
numeric._getCol = function _getCol(A,j,x) {
    var n = A.length, i;
    for(i=n-1;i>0;--i) {
        x[i] = A[i][j];
        --i;
        x[i] = A[i][j];
    }
    if(i===0) x[0] = A[0][j];
}
numeric.dotMMbig = function dotMMbig(x,y){
    var gc = numeric._getCol, p = y.length, v = Array(p);
    var m = x.length, n = y[0].length, A = new Array(m), xj;
    var VV = numeric.dotVV;
    var i,j,k,z;
    --p;
    --m;
    for(i=m;i!==-1;--i) A[i] = Array(n);
    --n;
    for(i=n;i!==-1;--i) {
        gc(y,i,v);
        for(j=m;j!==-1;--j) {
            z=0;
            xj = x[j];
            A[j][i] = VV(xj,v);
        }
    }
    return A;
}

numeric.dotMV = function dotMV(x,y) {
    var p = x.length, q = y.length,i;
    var ret = Array(p), dotVV = numeric.dotVV;
    for(i=p-1;i>=0;i--) { ret[i] = dotVV(x[i],y); }
    return ret;
}

numeric.dotVM = function dotVM(x,y) {
    var i,j,k,p,q,r,ret,foo,bar,woo,i0,k0,p0,r0,s1,s2,s3,baz,accum;
    p = x.length; q = y[0].length;
    ret = Array(q);
    for(k=q-1;k>=0;k--) {
        woo = x[p-1]*y[p-1][k];
        for(j=p-2;j>=1;j-=2) {
            i0 = j-1;
            woo += x[j]*y[j][k] + x[i0]*y[i0][k];
        }
        if(j===0) { woo += x[0]*y[0][k]; }
        ret[k] = woo;
    }
    return ret;
}

numeric.dotVV = function dotVV(x,y) {
    var i,n=x.length,i1,ret = x[n-1]*y[n-1];
    for(i=n-2;i>=1;i-=2) {
        i1 = i-1;
        ret += x[i]*y[i] + x[i1]*y[i1];
    }
    if(i===0) { ret += x[0]*y[0]; }
    return ret;
}

numeric.dot = function dot(x,y) {
    var d = numeric.dim;
    switch(d(x).length*1000+d(y).length) {
    case 2002:
        if(y.length < 10) return numeric.dotMMsmall(x,y);
        else return numeric.dotMMbig(x,y);
    case 2001: return numeric.dotMV(x,y);
    case 1002: return numeric.dotVM(x,y);
    case 1001: return numeric.dotVV(x,y);
    case 1000: return numeric.mulVS(x,y);
    case 1: return numeric.mulSV(x,y);
    case 0: return x*y;
    default: throw new Error('numeric.dot only works on vectors and matrices');
    }
}

numeric.diag = function diag(d) {
    var i,i1,j,n = d.length, A = Array(n), Ai;
    for(i=n-1;i>=0;i--) {
        Ai = Array(n);
        i1 = i+2;
        for(j=n-1;j>=i1;j-=2) {
            Ai[j] = 0;
            Ai[j-1] = 0;
        }
        if(j>i) { Ai[j] = 0; }
        Ai[i] = d[i];
        for(j=i-1;j>=1;j-=2) {
            Ai[j] = 0;
            Ai[j-1] = 0;
        }
        if(j===0) { Ai[0] = 0; }
        A[i] = Ai;
    }
    return A;
}
numeric.getDiag = function(A) {
    var n = Math.min(A.length,A[0].length),i,ret = Array(n);
    for(i=n-1;i>=1;--i) {
        ret[i] = A[i][i];
        --i;
        ret[i] = A[i][i];
    }
    if(i===0) {
        ret[0] = A[0][0];
    }
    return ret;
}

numeric.identity = function identity(n) { return numeric.diag(numeric.rep([n],1)); }
numeric.pointwise = function pointwise(params,body,setup) {
    if(typeof setup === "undefined") { setup = ""; }
    var fun = [];
    var k;
    var avec = /\[i\]$/,p,thevec = '';
    var haveret = false;
    for(k=0;k<params.length;k++) {
        if(avec.test(params[k])) {
            p = params[k].substring(0,params[k].length-3);
            thevec = p;
        } else { p = params[k]; }
        if(p==='ret') haveret = true;
        fun.push(p);
    }
    fun[params.length] = '_s';
    fun[params.length+1] = '_k';
    fun[params.length+2] = (
            'if(typeof _s === "undefined") _s = numeric.dim('+thevec+');\n'+
            'if(typeof _k === "undefined") _k = 0;\n'+
            'var _n = _s[_k];\n'+
            'var i'+(haveret?'':', ret = Array(_n)')+';\n'+
            'if(_k < _s.length-1) {\n'+
            '    for(i=_n-1;i>=0;i--) ret[i] = arguments.callee('+params.join(',')+',_s,_k+1);\n'+
            '    return ret;\n'+
            '}\n'+
            setup+'\n'+
            'for(i=_n-1;i!==-1;--i) {\n'+
            '    '+body+'\n'+
            '}\n'+
            'return ret;'
            );
    return Function.apply(null,fun);
}
numeric.pointwise2 = function pointwise2(params,body,setup) {
    if(typeof setup === "undefined") { setup = ""; }
    var fun = [];
    var k;
    var avec = /\[i\]$/,p,thevec = '';
    var haveret = false;
    for(k=0;k<params.length;k++) {
        if(avec.test(params[k])) {
            p = params[k].substring(0,params[k].length-3);
            thevec = p;
        } else { p = params[k]; }
        if(p==='ret') haveret = true;
        fun.push(p);
    }
    fun[params.length] = (
            'var _n = '+thevec+'.length;\n'+
            'var i'+(haveret?'':', ret = Array(_n)')+';\n'+
            setup+'\n'+
            'for(i=_n-1;i!==-1;--i) {\n'+
            body+'\n'+
            '}\n'+
            'return ret;'
            );
    return Function.apply(null,fun);
}
numeric._biforeach = (function _biforeach(x,y,s,k,f) {
    if(k === s.length-1) { f(x,y); return; }
    var i,n=s[k];
    for(i=n-1;i>=0;i--) { _biforeach(typeof x==="object"?x[i]:x,typeof y==="object"?y[i]:y,s,k+1,f); }
});
numeric._biforeach2 = (function _biforeach2(x,y,s,k,f) {
    if(k === s.length-1) { return f(x,y); }
    var i,n=s[k],ret = Array(n);
    for(i=n-1;i>=0;--i) { ret[i] = _biforeach2(typeof x==="object"?x[i]:x,typeof y==="object"?y[i]:y,s,k+1,f); }
    return ret;
});
numeric._foreach = (function _foreach(x,s,k,f) {
    if(k === s.length-1) { f(x); return; }
    var i,n=s[k];
    for(i=n-1;i>=0;i--) { _foreach(x[i],s,k+1,f); }
});
numeric._foreach2 = (function _foreach2(x,s,k,f) {
    if(k === s.length-1) { return f(x); }
    var i,n=s[k], ret = Array(n);
    for(i=n-1;i>=0;i--) { ret[i] = _foreach2(x[i],s,k+1,f); }
    return ret;
});

/*numeric.anyV = numeric.mapreduce('if(xi) return true;','false');
numeric.allV = numeric.mapreduce('if(!xi) return false;','true');
numeric.any = function(x) { if(typeof x.length === "undefined") return x; return numeric.anyV(x); }
numeric.all = function(x) { if(typeof x.length === "undefined") return x; return numeric.allV(x); }*/

numeric.ops2 = {
        add: '+',
        sub: '-',
        mul: '*',
        div: '/',
        mod: '%',
        and: '&&',
        or:  '||',
        eq:  '===',
        neq: '!==',
        lt:  '<',
        gt:  '>',
        leq: '<=',
        geq: '>=',
        band: '&',
        bor: '|',
        bxor: '^',
        lshift: '<<',
        rshift: '>>',
        rrshift: '>>>'
};
numeric.opseq = {
        addeq: '+=',
        subeq: '-=',
        muleq: '*=',
        diveq: '/=',
        modeq: '%=',
        lshifteq: '<<=',
        rshifteq: '>>=',
        rrshifteq: '>>>=',
        bandeq: '&=',
        boreq: '|=',
        bxoreq: '^='
};
numeric.mathfuns = ['abs','acos','asin','atan','ceil','cos',
                    'exp','floor','log','round','sin','sqrt','tan',
                    'isNaN','isFinite'];
numeric.mathfuns2 = ['atan2','pow','max','min'];
numeric.ops1 = {
        neg: '-',
        not: '!',
        bnot: '~',
        clone: ''
};
numeric.mapreducers = {
        any: ['if(xi) return true;','var accum = false;'],
        all: ['if(!xi) return false;','var accum = true;'],
        sum: ['accum += xi;','var accum = 0;'],
        prod: ['accum *= xi;','var accum = 1;'],
        norm2Squared: ['accum += xi*xi;','var accum = 0;'],
        norminf: ['accum = max(accum,abs(xi));','var accum = 0, max = Math.max, abs = Math.abs;'],
        norm1: ['accum += abs(xi)','var accum = 0, abs = Math.abs;'],
        sup: ['accum = max(accum,xi);','var accum = -Infinity, max = Math.max;'],
        inf: ['accum = min(accum,xi);','var accum = Infinity, min = Math.min;']
};

(function () {
    var i,o;
    for(i=0;i<numeric.mathfuns2.length;++i) {
        o = numeric.mathfuns2[i];
        numeric.ops2[o] = o;
    }
    for(i in numeric.ops2) {
        if(numeric.ops2.hasOwnProperty(i)) {
            o = numeric.ops2[i];
            var code, codeeq, setup = '';
            if(numeric.myIndexOf.call(numeric.mathfuns2,i)!==-1) {
                setup = 'var '+o+' = Math.'+o+';\n';
                code = function(r,x,y) { return r+' = '+o+'('+x+','+y+')'; };
                codeeq = function(x,y) { return x+' = '+o+'('+x+','+y+')'; };
            } else {
                code = function(r,x,y) { return r+' = '+x+' '+o+' '+y; };
                if(numeric.opseq.hasOwnProperty(i+'eq')) {
                    codeeq = function(x,y) { return x+' '+o+'= '+y; };
                } else {
                    codeeq = function(x,y) { return x+' = '+x+' '+o+' '+y; };                    
                }
            }
            numeric[i+'VV'] = numeric.pointwise2(['x[i]','y[i]'],code('ret[i]','x[i]','y[i]'),setup);
            numeric[i+'SV'] = numeric.pointwise2(['x','y[i]'],code('ret[i]','x','y[i]'),setup);
            numeric[i+'VS'] = numeric.pointwise2(['x[i]','y'],code('ret[i]','x[i]','y'),setup);
            numeric[i] = Function(
                    'var n = arguments.length, i, x = arguments[0], y;\n'+
                    'var VV = numeric.'+i+'VV, VS = numeric.'+i+'VS, SV = numeric.'+i+'SV;\n'+
                    'var dim = numeric.dim;\n'+
                    'for(i=1;i!==n;++i) { \n'+
                    '  y = arguments[i];\n'+
                    '  if(typeof x === "object") {\n'+
                    '      if(typeof y === "object") x = numeric._biforeach2(x,y,dim(x),0,VV);\n'+
                    '      else x = numeric._biforeach2(x,y,dim(x),0,VS);\n'+
                    '  } else if(typeof y === "object") x = numeric._biforeach2(x,y,dim(y),0,SV);\n'+
                    '  else '+codeeq('x','y')+'\n'+
                    '}\nreturn x;\n');
            numeric[o] = numeric[i];
            numeric[i+'eqV'] = numeric.pointwise2(['ret[i]','x[i]'], codeeq('ret[i]','x[i]'),setup);
            numeric[i+'eqS'] = numeric.pointwise2(['ret[i]','x'], codeeq('ret[i]','x'),setup);
            numeric[i+'eq'] = Function(
                    'var n = arguments.length, i, x = arguments[0], y;\n'+
                    'var V = numeric.'+i+'eqV, S = numeric.'+i+'eqS\n'+
                    'var s = numeric.dim(x);\n'+
                    'for(i=1;i!==n;++i) { \n'+
                    '  y = arguments[i];\n'+
                    '  if(typeof y === "object") numeric._biforeach(x,y,s,0,V);\n'+
                    '  else numeric._biforeach(x,y,s,0,S);\n'+
                    '}\nreturn x;\n');
        }
    }
    for(i=0;i<numeric.mathfuns2.length;++i) {
        o = numeric.mathfuns2[i];
        delete numeric.ops2[o];
    }
    for(i=0;i<numeric.mathfuns.length;++i) {
        o = numeric.mathfuns[i];
        numeric.ops1[o] = o;
    }
    for(i in numeric.ops1) {
        if(numeric.ops1.hasOwnProperty(i)) {
            setup = '';
            o = numeric.ops1[i];
            if(numeric.myIndexOf.call(numeric.mathfuns,i)!==-1) {
                if(Math.hasOwnProperty(o)) setup = 'var '+o+' = Math.'+o+';\n';
            }
            numeric[i+'eqV'] = numeric.pointwise2(['ret[i]'],'ret[i] = '+o+'(ret[i]);',setup);
            numeric[i+'eq'] = Function('x',
                    'if(typeof x !== "object") return '+o+'x\n'+
                    'var i;\n'+
                    'var V = numeric.'+i+'eqV;\n'+
                    'var s = numeric.dim(x);\n'+
                    'numeric._foreach(x,s,0,V);\n'+
                    'return x;\n');
            numeric[i+'V'] = numeric.pointwise2(['x[i]'],'ret[i] = '+o+'(x[i]);',setup);
            numeric[i] = Function('x',
                    'if(typeof x !== "object") return '+o+'(x)\n'+
                    'var i;\n'+
                    'var V = numeric.'+i+'V;\n'+
                    'var s = numeric.dim(x);\n'+
                    'return numeric._foreach2(x,s,0,V);\n');
        }
    }
    for(i=0;i<numeric.mathfuns.length;++i) {
        o = numeric.mathfuns[i];
        delete numeric.ops1[o];
    }
    for(i in numeric.mapreducers) {
        if(numeric.mapreducers.hasOwnProperty(i)) {
            o = numeric.mapreducers[i];
            numeric[i+'V'] = numeric.mapreduce2(o[0],o[1]);
            numeric[i] = Function('x','s','k',
                    o[1]+
                    'if(typeof x !== "object") {'+
                    '    xi = x;\n'+
                    o[0]+';\n'+
                    '    return accum;\n'+
                    '}'+
                    'if(typeof s === "undefined") s = numeric.dim(x);\n'+
                    'if(typeof k === "undefined") k = 0;\n'+
                    'if(k === s.length-1) return numeric.'+i+'V(x);\n'+
                    'var xi;\n'+
                    'var n = x.length, i;\n'+
                    'for(i=n-1;i!==-1;--i) {\n'+
                    '   xi = arguments.callee(x[i]);\n'+
                    o[0]+';\n'+
                    '}\n'+
                    'return accum;\n');
        }
    }
}());

numeric.truncVV = numeric.pointwise(['x[i]','y[i]'],'ret[i] = round(x[i]/y[i])*y[i];','var round = Math.round;');
numeric.truncVS = numeric.pointwise(['x[i]','y'],'ret[i] = round(x[i]/y)*y;','var round = Math.round;');
numeric.truncSV = numeric.pointwise(['x','y[i]'],'ret[i] = round(x/y[i])*y[i];','var round = Math.round;');
numeric.trunc = function trunc(x,y) {
    if(typeof x === "object") {
        if(typeof y === "object") return numeric.truncVV(x,y);
        return numeric.truncVS(x,y);
    }
    if (typeof y === "object") return numeric.truncSV(x,y);
    return Math.round(x/y)*y;
}

numeric.inv = function inv(x) {
    var s = numeric.dim(x), abs = Math.abs, m = s[0], n = s[1];
    var A = numeric.clone(x), Ai, Aj;
    var I = numeric.identity(m), Ii, Ij;
    var i,j,k,x;
    for(j=0;j<n;++j) {
        var i0 = -1;
        var v0 = -1;
        for(i=j;i!==m;++i) { k = abs(A[i][j]); if(k>v0) { i0 = i; v0 = k; } }
        Aj = A[i0]; A[i0] = A[j]; A[j] = Aj;
        Ij = I[i0]; I[i0] = I[j]; I[j] = Ij;
        x = Aj[j];
        for(k=j;k!==n;++k)    Aj[k] /= x; 
        for(k=n-1;k!==-1;--k) Ij[k] /= x;
        for(i=m-1;i!==-1;--i) {
            if(i!==j) {
                Ai = A[i];
                Ii = I[i];
                x = Ai[j];
                for(k=j+1;k!==n;++k)  Ai[k] -= Aj[k]*x;
                for(k=n-1;k>0;--k) { Ii[k] -= Ij[k]*x; --k; Ii[k] -= Ij[k]*x; }
                if(k===0) Ii[0] -= Ij[0]*x;
            }
        }
    }
    return I;
}

numeric.det = function det(x) {
    var s = numeric.dim(x);
    if(s.length !== 2 || s[0] !== s[1]) { throw new Error('numeric: det() only works on square matrices'); }
    var n = s[0], ret = 1,i,j,k,A = numeric.clone(x),Aj,Ai,alpha,temp,k1,k2,k3;
    for(j=0;j<n-1;j++) {
        k=j;
        for(i=j+1;i<n;i++) { if(Math.abs(A[i][j]) > Math.abs(A[k][j])) { k = i; } }
        if(k !== j) {
            temp = A[k]; A[k] = A[j]; A[j] = temp;
            ret *= -1;
        }
        Aj = A[j];
        for(i=j+1;i<n;i++) {
            Ai = A[i];
            alpha = Ai[j]/Aj[j];
            for(k=j+1;k<n-1;k+=2) {
                k1 = k+1;
                Ai[k] -= Aj[k]*alpha;
                Ai[k1] -= Aj[k1]*alpha;
            }
            if(k!==n) { Ai[k] -= Aj[k]*alpha; }
        }
        if(Aj[j] === 0) { return 0; }
        ret *= Aj[j];
    }
    return ret*A[j][j];
}

numeric.transpose = function transpose(x) {
    var i,j,m = x.length,n = x[0].length, ret=Array(n),A0,A1,Bj;
    for(j=0;j<n;j++) ret[j] = Array(m);
    for(i=m-1;i>=1;i-=2) {
        A1 = x[i];
        A0 = x[i-1];
        for(j=n-1;j>=1;--j) {
            Bj = ret[j]; Bj[i] = A1[j]; Bj[i-1] = A0[j];
            --j;
            Bj = ret[j]; Bj[i] = A1[j]; Bj[i-1] = A0[j];
        }
        if(j===0) {
            Bj = ret[0]; Bj[i] = A1[0]; Bj[i-1] = A0[0];
        }
    }
    if(i===0) {
        A0 = x[0];
        for(j=n-1;j>=1;--j) {
            ret[j][0] = A0[j];
            --j;
            ret[j][0] = A0[j];
        }
        if(j===0) { ret[0][0] = A0[0]; }
    }
    return ret;
}
numeric.negtranspose = function negtranspose(x) {
    var i,j,m = x.length,n = x[0].length, ret=Array(n),A0,A1,Bj;
    for(j=0;j<n;j++) ret[j] = Array(m);
    for(i=m-1;i>=1;i-=2) {
        A1 = x[i];
        A0 = x[i-1];
        for(j=n-1;j>=1;--j) {
            Bj = ret[j]; Bj[i] = -A1[j]; Bj[i-1] = -A0[j];
            --j;
            Bj = ret[j]; Bj[i] = -A1[j]; Bj[i-1] = -A0[j];
        }
        if(j===0) {
            Bj = ret[0]; Bj[i] = -A1[0]; Bj[i-1] = -A0[0];
        }
    }
    if(i===0) {
        A0 = x[0];
        for(j=n-1;j>=1;--j) {
            ret[j][0] = -A0[j];
            --j;
            ret[j][0] = -A0[j];
        }
        if(j===0) { ret[0][0] = -A0[0]; }
    }
    return ret;
}

numeric._random = function _random(s,k) {
    var i,n=s[k],ret=Array(n), rnd;
    if(k === s.length-1) {
        rnd = Math.random;
        for(i=n-1;i>=1;i-=2) {
            ret[i] = rnd();
            ret[i-1] = rnd();
        }
        if(i===0) { ret[0] = rnd(); }
        return ret;
    }
    for(i=n-1;i>=0;i--) ret[i] = _random(s,k+1);
    return ret;
}
numeric.random = function random(s) { return numeric._random(s,0); }

numeric.norm2 = function norm2(x) { return Math.sqrt(numeric.norm2Squared(x)); }

numeric.linspace = function linspace(a,b,n) {
    if(typeof n === "undefined") n = Math.max(Math.round(b-a)+1,1);
    if(n<2) { return n===1?[a]:[]; }
    var i,ret = Array(n);
    n--;
    for(i=n;i>=0;i--) { ret[i] = (i*b+(n-i)*a)/n; }
    return ret;
}

numeric.getBlock = function getBlock(x,from,to) {
    var s = numeric.dim(x);
    function foo(x,k) {
        var i,a = from[k], n = to[k]-a, ret = Array(n);
        if(k === s.length-1) {
            for(i=n;i>=0;i--) { ret[i] = x[i+a]; }
            return ret;
        }
        for(i=n;i>=0;i--) { ret[i] = foo(x[i+a],k+1); }
        return ret;
    }
    return foo(x,0);
}

numeric.setBlock = function setBlock(x,from,to,B) {
    var s = numeric.dim(x);
    function foo(x,y,k) {
        var i,a = from[k], n = to[k]-a;
        if(k === s.length-1) { for(i=n;i>=0;i--) { x[i+a] = y[i]; } }
        for(i=n;i>=0;i--) { foo(x[i+a],y[i],k+1); }
    }
    foo(x,B,0);
    return x;
}

numeric.getRange = function getRange(A,I,J) {
    var m = I.length, n = J.length;
    var i,j;
    var B = Array(m), Bi, AI;
    for(i=m-1;i!==-1;--i) {
        B[i] = Array(n);
        Bi = B[i];
        AI = A[I[i]];
        for(j=n-1;j!==-1;--j) Bi[j] = AI[J[j]];
    }
    return B;
}

numeric.blockMatrix = function blockMatrix(X) {
    var s = numeric.dim(X);
    if(s.length<4) return numeric.blockMatrix([X]);
    var m=s[0],n=s[1],M,N,i,j,Xij;
    M = 0; N = 0;
    for(i=0;i<m;++i) M+=X[i][0].length;
    for(j=0;j<n;++j) N+=X[0][j][0].length;
    var Z = Array(M);
    for(i=0;i<M;++i) Z[i] = Array(N);
    var I=0,J,ZI,k,l,Xijk;
    for(i=0;i<m;++i) {
        J=N;
        for(j=n-1;j!==-1;--j) {
            Xij = X[i][j];
            J -= Xij[0].length;
            for(k=Xij.length-1;k!==-1;--k) {
                Xijk = Xij[k];
                ZI = Z[I+k];
                for(l = Xijk.length-1;l!==-1;--l) ZI[J+l] = Xijk[l];
            }
        }
        I += X[i][0].length;
    }
    return Z;
}

numeric.tensor = function tensor(x,y) {
    if(typeof x === "number" || typeof y === "number") return numeric.mul(x,y);
    var s1 = numeric.dim(x), s2 = numeric.dim(y);
    if(s1.length !== 1 || s2.length !== 1) {
        throw new Error('numeric: tensor product is only defined for vectors');
    }
    var m = s1[0], n = s2[0], A = Array(m), Ai, i,j,xi;
    for(i=m-1;i>=0;i--) {
        Ai = Array(n);
        xi = x[i];
        for(j=n-1;j>=3;--j) {
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j];
            --j;
            Ai[j] = xi * y[j];
        }
        while(j>=0) { Ai[j] = xi * y[j]; --j; }
        A[i] = Ai;
    }
    return A;
}

// 3. The Tensor type T
numeric.T = function T(x,y) { this.x = x; this.y = y; }
numeric.t = function t(x,y) { return new numeric.T(x,y); }

numeric.Tbinop = function Tbinop(rr,rc,cr,cc,setup) {
    var io = numeric.indexOf;
    if(typeof setup !== "string") {
        var k;
        setup = '';
        for(k in numeric) {
            if(numeric.hasOwnProperty(k) && (rr.indexOf(k)>=0 || rc.indexOf(k)>=0 || cr.indexOf(k)>=0 || cc.indexOf(k)>=0) && k.length>1) {
                setup += 'var '+k+' = numeric.'+k+';\n';
            }
        }
    }
    return Function(['y'],
            'var x = this;\n'+
            'if(!(y instanceof numeric.T)) { y = new numeric.T(y); }\n'+
            setup+'\n'+
            'if(x.y) {'+
            '  if(y.y) {'+
            '    return new numeric.T('+cc+');\n'+
            '  }\n'+
            '  return new numeric.T('+cr+');\n'+
            '}\n'+
            'if(y.y) {\n'+
            '  return new numeric.T('+rc+');\n'+
            '}\n'+
            'return new numeric.T('+rr+');\n'
    );
}

numeric.T.prototype.add = numeric.Tbinop(
        'add(x.x,y.x)',
        'add(x.x,y.x),y.y',
        'add(x.x,y.x),x.y',
        'add(x.x,y.x),add(x.y,y.y)');
numeric.T.prototype.sub = numeric.Tbinop(
        'sub(x.x,y.x)',
        'sub(x.x,y.x),neg(y.y)',
        'sub(x.x,y.x),x.y',
        'sub(x.x,y.x),sub(x.y,y.y)');
numeric.T.prototype.mul = numeric.Tbinop(
        'mul(x.x,y.x)',
        'mul(x.x,y.x),mul(x.x,y.y)',
        'mul(x.x,y.x),mul(x.y,y.x)',
        'sub(mul(x.x,y.x),mul(x.y,y.y)),add(mul(x.x,y.y),mul(x.y,y.x))');

numeric.T.prototype.reciprocal = function reciprocal() {
    var mul = numeric.mul, div = numeric.div;
    if(this.y) {
        var d = numeric.add(mul(this.x,this.x),mul(this.y,this.y));
        return new numeric.T(div(this.x,d),div(numeric.neg(this.y),d));
    }
    return new T(div(1,this.x));
}
numeric.T.prototype.div = function div(y) {
    if(!(y instanceof numeric.T)) y = new numeric.T(y);
    if(y.y) { return this.mul(y.reciprocal()); }
    var div = numeric.div;
    if(this.y) { return new numeric.T(div(this.x,y.x),div(this.y,y.x)); }
    return new numeric.T(div(this.x,y.x));
}
numeric.T.prototype.dot = numeric.Tbinop(
        'dot(x.x,y.x)',
        'dot(x.x,y.x),dot(x.x,y.y)',
        'dot(x.x,y.x),dot(x.y,y.x)',
        'sub(dot(x.x,y.x),dot(x.y,y.y)),add(dot(x.x,y.y),dot(x.y,y.x))'
        );
numeric.T.prototype.transpose = function transpose() {
    var t = numeric.transpose, x = this.x, y = this.y;
    if(y) { return new numeric.T(t(x),t(y)); }
    return new numeric.T(t(x));
}
numeric.T.prototype.transjugate = function transjugate() {
    var t = numeric.transpose, x = this.x, y = this.y;
    if(y) { return new numeric.T(t(x),numeric.negtranspose(y)); }
    return new numeric.T(t(x));
}
numeric.Tunop = function Tunop(r,c,s) {
    if(typeof s !== "string") { s = ''; }
    return Function(
            'var x = this;\n'+
            s+'\n'+
            'if(x.y) {'+
            '  '+c+';\n'+
            '}\n'+
            r+';\n'
    );
}

numeric.T.prototype.exp = numeric.Tunop(
        'return new numeric.T(ex)',
        'return new numeric.T(mul(cos(x.y),ex),mul(sin(x.y),ex))',
        'var ex = numeric.exp(x.x), cos = numeric.cos, sin = numeric.sin, mul = numeric.mul;');
numeric.T.prototype.conj = numeric.Tunop(
        'return new numeric.T(x.x);',
        'return new numeric.T(x.x,numeric.neg(x.y));');
numeric.T.prototype.neg = numeric.Tunop(
        'return new numeric.T(neg(x.x));',
        'return new numeric.T(neg(x.x),neg(x.y));',
        'var neg = numeric.neg;');
numeric.T.prototype.sin = numeric.Tunop(
        'return new numeric.T(numeric.sin(x.x))',
        'return x.exp().sub(x.neg().exp()).div(new numeric.T(0,2));');
numeric.T.prototype.cos = numeric.Tunop(
        'return new numeric.T(numeric.cos(x.x))',
        'return x.exp().add(x.neg().exp()).div(2);');
numeric.T.prototype.abs = numeric.Tunop(
        'return new numeric.T(numeric.abs(x.x));',
        'return new numeric.T(numeric.sqrt(numeric.add(mul(x.x,x.x),mul(x.y,x.y))));',
        'var mul = numeric.mul;');
numeric.T.prototype.log = numeric.Tunop(
        'return new numeric.T(numeric.log(x.x));',
        'var theta = new numeric.T(numeric.atan2(x.y,x.x)), r = x.abs();\n'+
        'return new numeric.T(numeric.log(r.x),theta.x);');
numeric.T.prototype.norm2 = numeric.Tunop(
        'return numeric.norm2(x.x);',
        'var f = numeric.norm2Squared;\n'+
        'return Math.sqrt(f(x.x)+f(x.y));');
numeric.T.prototype.inv = function inv() {
    var A = this;
    if(typeof A.y === "undefined") { return new numeric.T(numeric.inv(A.x)); }
    var n = A.x.length, i, j, k;
    var Rx = numeric.identity(n),Ry = numeric.rep([n,n],0);
    var Ax = numeric.clone(A.x), Ay = numeric.clone(A.y);
    var Aix, Aiy, Ajx, Ajy, Rix, Riy, Rjx, Rjy;
    var i,j,k,d,d1,ax,ay,bx,by,temp;
    for(i=0;i<n;i++) {
        ax = Ax[i][i]; ay = Ay[i][i];
        d = ax*ax+ay*ay;
        k = i;
        for(j=i+1;j<n;j++) {
            ax = Ax[j][i]; ay = Ay[j][i];
            d1 = ax*ax+ay*ay;
            if(d1 > d) { k=j; d = d1; }
        }
        if(k!==i) {
            temp = Ax[i]; Ax[i] = Ax[k]; Ax[k] = temp;
            temp = Ay[i]; Ay[i] = Ay[k]; Ay[k] = temp;
            temp = Rx[i]; Rx[i] = Rx[k]; Rx[k] = temp;
            temp = Ry[i]; Ry[i] = Ry[k]; Ry[k] = temp;
        }
        Aix = Ax[i]; Aiy = Ay[i];
        Rix = Rx[i]; Riy = Ry[i];
        ax = Aix[i]; ay = Aiy[i];
        for(j=i+1;j<n;j++) {
            bx = Aix[j]; by = Aiy[j];
            Aix[j] = (bx*ax+by*ay)/d;
            Aiy[j] = (by*ax-bx*ay)/d;
        }
        for(j=0;j<n;j++) {
            bx = Rix[j]; by = Riy[j];
            Rix[j] = (bx*ax+by*ay)/d;
            Riy[j] = (by*ax-bx*ay)/d;
        }
        for(j=i+1;j<n;j++) {
            Ajx = Ax[j]; Ajy = Ay[j];
            Rjx = Rx[j]; Rjy = Ry[j];
            ax = Ajx[i]; ay = Ajy[i];
            for(k=i+1;k<n;k++) {
                bx = Aix[k]; by = Aiy[k];
                Ajx[k] -= bx*ax-by*ay;
                Ajy[k] -= by*ax+bx*ay;
            }
            for(k=0;k<n;k++) {
                bx = Rix[k]; by = Riy[k];
                Rjx[k] -= bx*ax-by*ay;
                Rjy[k] -= by*ax+bx*ay;
            }
        }
    }
    for(i=n-1;i>0;i--) {
        Rix = Rx[i]; Riy = Ry[i];
        for(j=i-1;j>=0;j--) {
            Rjx = Rx[j]; Rjy = Ry[j];
            ax = Ax[j][i]; ay = Ay[j][i];
            for(k=n-1;k>=0;k--) {
                bx = Rix[k]; by = Riy[k];
                Rjx[k] -= ax*bx - ay*by;
                Rjy[k] -= ax*by + ay*bx;
            }
        }
    }
    return new numeric.T(Rx,Ry);
}
numeric.T.prototype.get = function get(i) {
    var x = this.x, y = this.y, k = 0, ik, n = i.length;
    if(y) {
        while(k<n) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++;
        }
        return new numeric.T(x,y);
    }
    while(k<n) {
        ik = i[k];
        x = x[ik];
        k++;
    }
    return new numeric.T(x);
}
numeric.T.prototype.set = function set(i,v) {
    var x = this.x, y = this.y, k = 0, ik, n = i.length, vx = v.x, vy = v.y;
    if(n===0) {
        if(vy) { this.y = vy; }
        else if(y) { this.y = undefined; }
        this.x = x;
        return this;
    }
    if(vy) {
        if(y) { /* ok */ }
        else {
            y = numeric.rep(numeric.dim(x),0);
            this.y = y;
        }
        while(k<n-1) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++;
        }
        ik = i[k];
        x[ik] = vx;
        y[ik] = vy;
        return this;
    }
    if(y) {
        while(k<n-1) {
            ik = i[k];
            x = x[ik];
            y = y[ik];
            k++;
        }
        ik = i[k];
        x[ik] = vx;
        if(vx instanceof Array) y[ik] = numeric.rep(numeric.dim(vx),0);
        else y[ik] = 0;
        return this;
    }
    while(k<n-1) {
        ik = i[k];
        x = x[ik];
        k++;
    }
    ik = i[k];
    x[ik] = vx;
    return this;
}
numeric.T.prototype.getRows = function getRows(i0,i1) {
    var n = i1-i0+1, j;
    var rx = Array(n), ry, x = this.x, y = this.y;
    for(j=i0;j<=i1;j++) { rx[j-i0] = x[j]; }
    if(y) {
        ry = Array(n);
        for(j=i0;j<=i1;j++) { ry[j-i0] = y[j]; }
        return new numeric.T(rx,ry);
    }
    return new numeric.T(rx);
}
numeric.T.prototype.setRows = function setRows(i0,i1,A) {
    var j;
    var rx = this.x, ry = this.y, x = A.x, y = A.y;
    for(j=i0;j<=i1;j++) { rx[j] = x[j-i0]; }
    if(y) {
        if(!ry) { ry = numeric.rep(numeric.dim(rx),0); this.y = ry; }
        for(j=i0;j<=i1;j++) { ry[j] = y[j-i0]; }
    } else if(ry) {
        for(j=i0;j<=i1;j++) { ry[j] = numeric.rep([x[j-i0].length],0); }
    }
    return this;
}
numeric.T.prototype.getRow = function getRow(k) {
    var x = this.x, y = this.y;
    if(y) { return new numeric.T(x[k],y[k]); }
    return new numeric.T(x[k]);
}
numeric.T.prototype.setRow = function setRow(i,v) {
    var rx = this.x, ry = this.y, x = v.x, y = v.y;
    rx[i] = x;
    if(y) {
        if(!ry) { ry = numeric.rep(numeric.dim(rx),0); this.y = ry; }
        ry[i] = y;
    } else if(ry) {
        ry = numeric.rep([x.length],0);
    }
    return this;
}

numeric.T.prototype.getBlock = function getBlock(from,to) {
    var x = this.x, y = this.y, b = numeric.getBlock;
    if(y) { return new numeric.T(b(x,from,to),b(y,from,to)); }
    return new numeric.T(b(x,from,to));
}
numeric.T.prototype.setBlock = function setBlock(from,to,A) {
    if(!(A instanceof numeric.T)) A = new numeric.T(A);
    var x = this.x, y = this.y, b = numeric.setBlock, Ax = A.x, Ay = A.y;
    if(Ay) {
        if(!y) { this.y = numeric.rep(numeric.dim(this),0); y = this.y; }
        b(x,from,to,Ax);
        b(y,from,to,Ay);
        return this;
    }
    b(x,from,to,Ax);
    if(y) b(y,from,to,numeric.rep(numeric.dim(Ax),0));
}
numeric.T.rep = function rep(s,v) {
    var T = numeric.T;
    if(!(v instanceof T)) v = new T(v);
    var x = v.x, y = v.y, r = numeric.rep;
    if(y) return new T(r(s,x),r(s,y));
    return new T(r(s,x));
}
numeric.T.diag = function diag(d) {
    if(!(d instanceof numeric.T)) d = new numeric.T(d);
    var x = d.x, y = d.y, diag = numeric.diag;
    if(y) return new numeric.T(diag(x),diag(y));
    return new numeric.T(diag(x));
}
numeric.T.eig = function eig() {
    if(this.y) { throw new Error('eig: not implemented for complex matrices.'); }
    return numeric.eig(this.x);
}
numeric.T.identity = function identity(n) { return new numeric.T(numeric.identity(n)); }
numeric.T.prototype.getDiag = function getDiag() {
    var n = numeric;
    var x = this.x, y = this.y;
    if(y) { return new n.T(n.getDiag(x),n.getDiag(y)); }
    return new n.T(n.getDiag(x));
}

// 4. Eigenvalues of real matrices

numeric.house = function house(x) {
    var v = numeric.clone(x);
    var s = x[0] >= 0 ? 1 : -1;
    var alpha = s*numeric.norm2(x);
    v[0] += alpha;
    var foo = numeric.norm2(v);
    if(foo === 0) { /* this should not happen */ throw new Error('eig: internal error'); }
    return numeric.div(v,foo);
}

numeric.toUpperHessenberg = function toUpperHessenberg(me) {
    var s = numeric.dim(me);
    if(s.length !== 2 || s[0] !== s[1]) { throw new Error('numeric: toUpperHessenberg() only works on square matrices'); }
    var m = s[0], i,j,k,x,v,A = numeric.clone(me),B,C,Ai,Ci,Q = numeric.identity(m),Qi;
    for(j=0;j<m-2;j++) {
        x = Array(m-j-1);
        for(i=j+1;i<m;i++) { x[i-j-1] = A[i][j]; }
        if(numeric.norm2(x)>0) {
            v = numeric.house(x);
            B = numeric.getBlock(A,[j+1,j],[m-1,m-1]);
            C = numeric.tensor(v,numeric.dot(v,B));
            for(i=j+1;i<m;i++) { Ai = A[i]; Ci = C[i-j-1]; for(k=j;k<m;k++) Ai[k] -= 2*Ci[k-j]; }
            B = numeric.getBlock(A,[0,j+1],[m-1,m-1]);
            C = numeric.tensor(numeric.dot(B,v),v);
            for(i=0;i<m;i++) { Ai = A[i]; Ci = C[i]; for(k=j+1;k<m;k++) Ai[k] -= 2*Ci[k-j-1]; }
            B = Array(m-j-1);
            for(i=j+1;i<m;i++) B[i-j-1] = Q[i];
            C = numeric.tensor(v,numeric.dot(v,B));
            for(i=j+1;i<m;i++) { Qi = Q[i]; Ci = C[i-j-1]; for(k=0;k<m;k++) Qi[k] -= 2*Ci[k]; }
        }
    }
    return {H:A, Q:Q};
}

numeric.epsilon = 2.220446049250313e-16;

numeric.QRFrancis = function(H,maxiter) {
    if(typeof maxiter === "undefined") { maxiter = 10000; }
    H = numeric.clone(H);
    var H0 = numeric.clone(H);
    var s = numeric.dim(H),m=s[0],x,v,a,b,c,d,det,tr, Hloc, Q = numeric.identity(m), Qi, Hi, B, C, Ci,i,j,k,iter;
    if(m<3) { return {Q:Q, B:[ [0,m-1] ]}; }
    var epsilon = numeric.epsilon;
    for(iter=0;iter<maxiter;iter++) {
        for(j=0;j<m-1;j++) {
            if(Math.abs(H[j+1][j]) < epsilon*(Math.abs(H[j][j])+Math.abs(H[j+1][j+1]))) {
                var QH1 = numeric.QRFrancis(numeric.getBlock(H,[0,0],[j,j]),maxiter);
                var QH2 = numeric.QRFrancis(numeric.getBlock(H,[j+1,j+1],[m-1,m-1]),maxiter);
                B = Array(j+1);
                for(i=0;i<=j;i++) { B[i] = Q[i]; }
                C = numeric.dot(QH1.Q,B);
                for(i=0;i<=j;i++) { Q[i] = C[i]; }
                B = Array(m-j-1);
                for(i=j+1;i<m;i++) { B[i-j-1] = Q[i]; }
                C = numeric.dot(QH2.Q,B);
                for(i=j+1;i<m;i++) { Q[i] = C[i-j-1]; }
                return {Q:Q,B:QH1.B.concat(numeric.add(QH2.B,j+1))};
            }
        }
        a = H[m-2][m-2]; b = H[m-2][m-1];
        c = H[m-1][m-2]; d = H[m-1][m-1];
        tr = a+d;
        det = (a*d-b*c);
        Hloc = numeric.getBlock(H, [0,0], [2,2]);
        if(tr*tr>=4*det) {
            var s1,s2;
            s1 = 0.5*(tr+Math.sqrt(tr*tr-4*det));
            s2 = 0.5*(tr-Math.sqrt(tr*tr-4*det));
            Hloc = numeric.add(numeric.sub(numeric.dot(Hloc,Hloc),
                                           numeric.mul(Hloc,s1+s2)),
                               numeric.diag(numeric.rep([3],s1*s2)));
        } else {
            Hloc = numeric.add(numeric.sub(numeric.dot(Hloc,Hloc),
                                           numeric.mul(Hloc,tr)),
                               numeric.diag(numeric.rep([3],det)));
        }
        x = [Hloc[0][0],Hloc[1][0],Hloc[2][0]];
        v = numeric.house(x);
        B = [H[0],H[1],H[2]];
        C = numeric.tensor(v,numeric.dot(v,B));
        for(i=0;i<3;i++) { Hi = H[i]; Ci = C[i]; for(k=0;k<m;k++) Hi[k] -= 2*Ci[k]; }
        B = numeric.getBlock(H, [0,0],[m-1,2]);
        C = numeric.tensor(numeric.dot(B,v),v);
        for(i=0;i<m;i++) { Hi = H[i]; Ci = C[i]; for(k=0;k<3;k++) Hi[k] -= 2*Ci[k]; }
        B = [Q[0],Q[1],Q[2]];
        C = numeric.tensor(v,numeric.dot(v,B));
        for(i=0;i<3;i++) { Qi = Q[i]; Ci = C[i]; for(k=0;k<m;k++) Qi[k] -= 2*Ci[k]; }
        var J;
        for(j=0;j<m-2;j++) {
            for(k=j;k<=j+1;k++) {
                if(Math.abs(H[k+1][k]) < epsilon*(Math.abs(H[k][k])+Math.abs(H[k+1][k+1]))) {
                    var QH1 = numeric.QRFrancis(numeric.getBlock(H,[0,0],[k,k]),maxiter);
                    var QH2 = numeric.QRFrancis(numeric.getBlock(H,[k+1,k+1],[m-1,m-1]),maxiter);
                    B = Array(k+1);
                    for(i=0;i<=k;i++) { B[i] = Q[i]; }
                    C = numeric.dot(QH1.Q,B);
                    for(i=0;i<=k;i++) { Q[i] = C[i]; }
                    B = Array(m-k-1);
                    for(i=k+1;i<m;i++) { B[i-k-1] = Q[i]; }
                    C = numeric.dot(QH2.Q,B);
                    for(i=k+1;i<m;i++) { Q[i] = C[i-k-1]; }
                    return {Q:Q,B:QH1.B.concat(numeric.add(QH2.B,k+1))};
                }
            }
            J = Math.min(m-1,j+3);
            x = Array(J-j);
            for(i=j+1;i<=J;i++) { x[i-j-1] = H[i][j]; }
            v = numeric.house(x);
            B = numeric.getBlock(H, [j+1,j],[J,m-1]);
            C = numeric.tensor(v,numeric.dot(v,B));
            for(i=j+1;i<=J;i++) { Hi = H[i]; Ci = C[i-j-1]; for(k=j;k<m;k++) Hi[k] -= 2*Ci[k-j]; }
            B = numeric.getBlock(H, [0,j+1],[m-1,J]);
            C = numeric.tensor(numeric.dot(B,v),v);
            for(i=0;i<m;i++) { Hi = H[i]; Ci = C[i]; for(k=j+1;k<=J;k++) Hi[k] -= 2*Ci[k-j-1]; }
            B = Array(J-j);
            for(i=j+1;i<=J;i++) B[i-j-1] = Q[i];
            C = numeric.tensor(v,numeric.dot(v,B));
            for(i=j+1;i<=J;i++) { Qi = Q[i]; Ci = C[i-j-1]; for(k=0;k<m;k++) Qi[k] -= 2*Ci[k]; }
        }
    }
    throw new Error('numeric: eigenvalue iteration does not converge -- increase maxiter?');
}

numeric.eig = function eig(A,maxiter) {
    var QH = numeric.toUpperHessenberg(A);
    var QB = numeric.QRFrancis(QH.H,maxiter);
    var T = numeric.T;
    var n = A.length,i,k,flag = false,B = QB.B,H = numeric.dot(QB.Q,numeric.dot(QH.H,numeric.transpose(QB.Q)));
    var Q = new T(numeric.dot(QB.Q,QH.Q)),Q0;
    var m = B.length,j;
    var a,b,c,d,p1,p2,disc,x,y,p,q,n1,n2;
    var sqrt = Math.sqrt;
    for(k=0;k<m;k++) {
        i = B[k][0];
        if(i === B[k][1]) {
            // nothing
        } else {
            j = i+1;
            a = H[i][i];
            b = H[i][j];
            c = H[j][i];
            d = H[j][j];
            if(b === 0 && c === 0) continue;
            p1 = -a-d;
            p2 = a*d-b*c;
            disc = p1*p1-4*p2;
            if(disc>=0) {
                if(p1<0) x = -0.5*(p1-sqrt(disc));
                else     x = -0.5*(p1+sqrt(disc));
                n1 = (a-x)*(a-x)+b*b;
                n2 = c*c+(d-x)*(d-x);
                if(n1>n2) {
                    n1 = sqrt(n1);
                    p = (a-x)/n1;
                    q = b/n1;
                } else {
                    n2 = sqrt(n2);
                    p = c/n2;
                    q = (d-x)/n2;
                }
                Q0 = new T([[q,-p],[p,q]]);
                Q.setRows(i,j,Q0.dot(Q.getRows(i,j)));
            } else {
                x = -0.5*p1;
                y = 0.5*sqrt(-disc);
                n1 = (a-x)*(a-x)+b*b;
                n2 = c*c+(d-x)*(d-x);
                if(n1>n2) {
                    n1 = sqrt(n1+y*y);
                    p = (a-x)/n1;
                    q = b/n1;
                    x = 0;
                    y /= n1;
                } else {
                    n2 = sqrt(n2+y*y);
                    p = c/n2;
                    q = (d-x)/n2;
                    x = y/n2;
                    y = 0;
                }
                Q0 = new T([[q,-p],[p,q]],[[x,y],[y,-x]]);
                Q.setRows(i,j,Q0.dot(Q.getRows(i,j)));
            }
        }
    }
    var R = Q.dot(A).dot(Q.transjugate()), n = A.length, E = numeric.T.identity(n);
    for(j=0;j<n;j++) {
        if(j>0) {
            for(k=j-1;k>=0;k--) {
                var Rk = R.get([k,k]), Rj = R.get([j,j]);
                if(numeric.neq(Rk.x,Rj.x) || numeric.neq(Rk.y,Rj.y)) {
                    x = R.getRow(k).getBlock([k],[j-1]);
                    y = E.getRow(j).getBlock([k],[j-1]);
                    E.set([j,k],(R.get([k,j]).neg().sub(x.dot(y))).div(Rk.sub(Rj)));
                } else {
                    E.setRow(j,E.getRow(k));
                    continue;
                }
            }
        }
    }
    for(j=0;j<n;j++) {
        x = E.getRow(j);
        E.setRow(j,x.div(x.norm2()));
    }
    E = E.transpose();
    E = Q.transjugate().dot(E);
    return { lambda:R.getDiag(), E:E };
};

// 5. Compressed Column Storage matrices
numeric.ccsSparse = function ccsSparse(A) {
    var m = A.length,n,foo, i,j, counts = [];
    for(i=m-1;i!==-1;--i) {
        foo = A[i];
        for(j in foo) {
            j = parseInt(j);
            while(j>=counts.length) counts[counts.length] = 0;
            if(foo[j]!==0) counts[j]++;
        }
    }
    var n = counts.length;
    var Ai = Array(n+1);
    Ai[0] = 0;
    for(i=0;i<n;++i) Ai[i+1] = Ai[i] + counts[i];
    var Aj = Array(Ai[n]), Av = Array(Ai[n]);
    for(i=m-1;i!==-1;--i) {
        foo = A[i];
        for(j in foo) {
            if(foo[j]!==0) {
                counts[j]--;
                Aj[Ai[j]+counts[j]] = i;
                Av[Ai[j]+counts[j]] = foo[j];
            }
        }
    }
    return [Ai,Aj,Av];
}
numeric.ccsFull = function ccsFull(A) {
    var Ai = A[0], Aj = A[1], Av = A[2], s = numeric.ccsDim(A), m = s[0], n = s[1], i,j,j0,j1,k;
    var B = numeric.rep([m,n],0);
    for(i=0;i<n;i++) {
        j0 = Ai[i];
        j1 = Ai[i+1];
        for(j=j0;j<j1;++j) { B[Aj[j]][i] = Av[j]; }
    }
    return B;
}
numeric.ccsTSolve = function ccsTSolve(A,b,x,bj,xj) {
    var Ai = A[0], Aj = A[1], Av = A[2],m = Ai.length-1, max = Math.max,n=0;
    if(typeof bj === "undefined") x = numeric.rep([m],0);
    if(typeof bj === "undefined") bj = numeric.linspace(0,x.length-1);
    if(typeof xj === "undefined") xj = [];
    function dfs(j) {
        var k;
        if(x[j] !== 0) return;
        x[j] = 1;
        for(k=Ai[j];k<Ai[j+1];++k) dfs(Aj[k]);
        xj[n] = j;
        ++n;
    }
    var i,j,j0,j1,k,l,l0,l1,a;
    for(i=bj.length-1;i!==-1;--i) { dfs(bj[i]); }
    xj.length = n;
    for(i=xj.length-1;i!==-1;--i) { x[xj[i]] = 0; }
    for(i=bj.length-1;i!==-1;--i) { j = bj[i]; x[j] = b[j]; }
    for(i=xj.length-1;i!==-1;--i) {
        j = xj[i];
        j0 = Ai[j];
        j1 = max(Ai[j+1],j0);
        for(k=j0;k!==j1;++k) { if(Aj[k] === j) { x[j] /= Av[k]; break; } }
        a = x[j];
        for(k=j0;k!==j1;++k) {
            l = Aj[k];
            if(l !== j) x[l] -= a*Av[k];
        }
    }
    return x;
}
numeric.ccsDFS = function ccsDFS(n) {
    this.k = Array(n);
    this.k1 = Array(n);
    this.j = Array(n);
}
numeric.ccsDFS.prototype.dfs = function dfs(J,Ai,Aj,x,xj,Pinv) {
    var m = 0,foo,n=xj.length;
    var k = this.k, k1 = this.k1, j = this.j,km,k11;
    if(x[J]!==0) return;
    x[J] = 1;
    j[0] = J;
    k[0] = km = Ai[J];
    k1[0] = k11 = Ai[J+1];
    while(1) {
        if(km >= k11) {
            xj[n] = j[m];
            if(m===0) return;
            ++n;
            --m;
            km = k[m];
            k11 = k1[m];
        } else {
            foo = Pinv[Aj[km]];
            if(x[foo] === 0) {
                x[foo] = 1;
                k[m] = km;
                ++m;
                j[m] = foo;
                km = Ai[foo];
                k1[m] = k11 = Ai[foo+1];
            } else ++km;
        }
    }
}
numeric.ccsLPSolve = function ccsLPSolve(A,B,x,xj,I,Pinv,dfs) {
    var Ai = A[0], Aj = A[1], Av = A[2],m = Ai.length-1, n=0;
    var Bi = B[0], Bj = B[1], Bv = B[2];
    
    var i,i0,i1,j,J,j0,j1,k,l,l0,l1,a;
    i0 = Bi[I];
    i1 = Bi[I+1];
    xj.length = 0;
    for(i=i0;i<i1;++i) { dfs.dfs(Pinv[Bj[i]],Ai,Aj,x,xj,Pinv); }
    for(i=xj.length-1;i!==-1;--i) { x[xj[i]] = 0; }
    for(i=i0;i!==i1;++i) { j = Pinv[Bj[i]]; x[j] = Bv[i]; }
    for(i=xj.length-1;i!==-1;--i) {
        j = xj[i];
        j0 = Ai[j];
        j1 = Ai[j+1];
        for(k=j0;k<j1;++k) { if(Pinv[Aj[k]] === j) { x[j] /= Av[k]; break; } }
        a = x[j];
        for(k=j0;k<j1;++k) {
            l = Pinv[Aj[k]];
            if(l !== j) x[l] -= a*Av[k];
        }
    }
    return x;
}
numeric.ccsLUP1 = function ccsLUP1(A,threshold) {
    var m = A[0].length-1;
    var L = [numeric.rep([m+1],0),[],[]], U = [numeric.rep([m+1], 0),[],[]];
    var Li = L[0], Lj = L[1], Lv = L[2], Ui = U[0], Uj = U[1], Uv = U[2];
    var x = numeric.rep([m],0), xj = numeric.rep([m],0);
    var i,j,k,j0,j1,a,e,c,d,K;
    var sol = numeric.ccsLPSolve, max = Math.max, abs = Math.abs;
    var P = numeric.linspace(0,m-1),Pinv = numeric.linspace(0,m-1);
    var dfs = new numeric.ccsDFS(m);
    if(typeof threshold === "undefined") { threshold = 1; }
    for(i=0;i<m;++i) {
        sol(L,A,x,xj,i,Pinv,dfs);
        a = -1;
        e = -1;
        for(j=xj.length-1;j!==-1;--j) {
            k = xj[j];
            if(k <= i) continue;
            c = abs(x[k]);
            if(c > a) { e = k; a = c; }
        }
        if(abs(x[i])<threshold*a) {
            j = P[i];
            a = P[e];
            P[i] = a; Pinv[a] = i;
            P[e] = j; Pinv[j] = e;
            a = x[i]; x[i] = x[e]; x[e] = a;
        }
        a = Li[i];
        e = Ui[i];
        d = x[i];
        Lj[a] = P[i];
        Lv[a] = 1;
        ++a;
        for(j=xj.length-1;j!==-1;--j) {
            k = xj[j];
            c = x[k];
            xj[j] = 0;
            x[k] = 0;
            if(k<=i) { Uj[e] = k; Uv[e] = c;   ++e; }
            else     { Lj[a] = P[k]; Lv[a] = c/d; ++a; }
        }
        Li[i+1] = a;
        Ui[i+1] = e;
    }
    for(j=Lj.length-1;j!==-1;--j) { Lj[j] = Pinv[Lj[j]]; }
    return {L:L, U:U, P:P, Pinv:Pinv};
}
numeric.ccsDFS0 = function ccsDFS0(n) {
    this.k = Array(n);
    this.k1 = Array(n);
    this.j = Array(n);
}
numeric.ccsDFS0.prototype.dfs = function dfs(J,Ai,Aj,x,xj,Pinv,P) {
    var m = 0,foo,n=xj.length;
    var k = this.k, k1 = this.k1, j = this.j,km,k11;
    if(x[J]!==0) return;
    x[J] = 1;
    j[0] = J;
    k[0] = km = Ai[Pinv[J]];
    k1[0] = k11 = Ai[Pinv[J]+1];
    while(1) {
        if(isNaN(km)) throw new Error("Ow!");
        if(km >= k11) {
            xj[n] = Pinv[j[m]];
            if(m===0) return;
            ++n;
            --m;
            km = k[m];
            k11 = k1[m];
        } else {
            foo = Aj[km];
            if(x[foo] === 0) {
                x[foo] = 1;
                k[m] = km;
                ++m;
                j[m] = foo;
                foo = Pinv[foo];
                km = Ai[foo];
                k1[m] = k11 = Ai[foo+1];
            } else ++km;
        }
    }
}
numeric.ccsLPSolve0 = function ccsLPSolve0(A,B,y,xj,I,Pinv,P,dfs) {
    var Ai = A[0], Aj = A[1], Av = A[2],m = Ai.length-1, n=0;
    var Bi = B[0], Bj = B[1], Bv = B[2];
    
    var i,i0,i1,j,J,j0,j1,k,l,l0,l1,a;
    i0 = Bi[I];
    i1 = Bi[I+1];
    xj.length = 0;
    for(i=i0;i<i1;++i) { dfs.dfs(Bj[i],Ai,Aj,y,xj,Pinv,P); }
    for(i=xj.length-1;i!==-1;--i) { j = xj[i]; y[P[j]] = 0; }
    for(i=i0;i!==i1;++i) { j = Bj[i]; y[j] = Bv[i]; }
    for(i=xj.length-1;i!==-1;--i) {
        j = xj[i];
        l = P[j];
        j0 = Ai[j];
        j1 = Ai[j+1];
        for(k=j0;k<j1;++k) { if(Aj[k] === l) { y[l] /= Av[k]; break; } }
        a = y[l];
        for(k=j0;k<j1;++k) y[Aj[k]] -= a*Av[k];
        y[l] = a;
    }
}
numeric.ccsLUP0 = function ccsLUP0(A,threshold) {
    var m = A[0].length-1;
    var L = [numeric.rep([m+1],0),[],[]], U = [numeric.rep([m+1], 0),[],[]];
    var Li = L[0], Lj = L[1], Lv = L[2], Ui = U[0], Uj = U[1], Uv = U[2];
    var y = numeric.rep([m],0), xj = numeric.rep([m],0);
    var i,j,k,j0,j1,a,e,c,d,K;
    var sol = numeric.ccsLPSolve0, max = Math.max, abs = Math.abs;
    var P = numeric.linspace(0,m-1),Pinv = numeric.linspace(0,m-1);
    var dfs = new numeric.ccsDFS0(m);
    if(typeof threshold === "undefined") { threshold = 1; }
    for(i=0;i<m;++i) {
        sol(L,A,y,xj,i,Pinv,P,dfs);
        a = -1;
        e = -1;
        for(j=xj.length-1;j!==-1;--j) {
            k = xj[j];
            if(k <= i) continue;
            c = abs(y[P[k]]);
            if(c > a) { e = k; a = c; }
        }
        if(abs(y[P[i]])<threshold*a) {
            j = P[i];
            a = P[e];
            P[i] = a; Pinv[a] = i;
            P[e] = j; Pinv[j] = e;
        }
        a = Li[i];
        e = Ui[i];
        d = y[P[i]];
        Lj[a] = P[i];
        Lv[a] = 1;
        ++a;
        for(j=xj.length-1;j!==-1;--j) {
            k = xj[j];
            c = y[P[k]];
            xj[j] = 0;
            y[P[k]] = 0;
            if(k<=i) { Uj[e] = k; Uv[e] = c;   ++e; }
            else     { Lj[a] = P[k]; Lv[a] = c/d; ++a; }
        }
        Li[i+1] = a;
        Ui[i+1] = e;
    }
    for(j=Lj.length-1;j!==-1;--j) { Lj[j] = Pinv[Lj[j]]; }
    return {L:L, U:U, P:P, Pinv:Pinv};
}
numeric.ccsLUP = numeric.ccsLUP0;

numeric.ccsDim = function ccsDim(A) { return [numeric.sup(A[1])+1,A[0].length-1]; }
numeric.ccsGetBlock = function ccsGetBlock(A,i,j) {
    var s = numeric.ccsDim(A),m=s[0],n=s[1];
    if(typeof i === "undefined") { i = numeric.linspace(0,m-1); }
    else if(typeof i === "number") { i = [i]; }
    if(typeof j === "undefined") { j = numeric.linspace(0,n-1); }
    else if(typeof j === "number") { j = [j]; }
    var p,p0,p1,P = i.length,q,Q = j.length,r,jq,ip;
    var Bi = numeric.rep([n],0), Bj=[], Bv=[], B = [Bi,Bj,Bv];
    var Ai = A[0], Aj = A[1], Av = A[2];
    var x = numeric.rep([m],0),count=0,flags = numeric.rep([m],0);
    for(q=0;q<Q;++q) {
        jq = j[q];
        var q0 = Ai[jq];
        var q1 = Ai[jq+1];
        for(p=q0;p<q1;++p) {
            r = Aj[p];
            flags[r] = 1;
            x[r] = Av[p];
        }
        for(p=0;p<P;++p) {
            ip = i[p];
            if(flags[ip]) {
                Bj[count] = p;
                Bv[count] = x[i[p]];
                ++count;
            }
        }
        for(p=q0;p<q1;++p) {
            r = Aj[p];
            flags[r] = 0;
        }
        Bi[q+1] = count;
    }
    return B;
}

numeric.ccsDot = function ccsDot(A,B) {
    var Ai = A[0], Aj = A[1], Av = A[2];
    var Bi = B[0], Bj = B[1], Bv = B[2];
    var sA = numeric.ccsDim(A), sB = numeric.ccsDim(B);
    var m = sA[0], n = sA[1], o = sB[1];
    var x = numeric.rep([m],0), flags = numeric.rep([m],0), xj = Array(m);
    var Ci = numeric.rep([o],0), Cj = [], Cv = [], C = [Ci,Cj,Cv];
    var i,j,k,j0,j1,i0,i1,l,p,a,b;
    for(k=0;k!==o;++k) {
        j0 = Bi[k];
        j1 = Bi[k+1];
        p = 0;
        for(j=j0;j<j1;++j) {
            a = Bj[j];
            b = Bv[j];
            i0 = Ai[a];
            i1 = Ai[a+1];
            for(i=i0;i<i1;++i) {
                l = Aj[i];
                if(flags[l]===0) {
                    xj[p] = l;
                    flags[l] = 1;
                    p = p+1;
                }
                x[l] = x[l] + Av[i]*b;
            }
        }
        j0 = Ci[k];
        j1 = j0+p;
        Ci[k+1] = j1;
        for(j=p-1;j!==-1;--j) {
            b = j0+j;
            i = xj[j];
            Cj[b] = i;
            Cv[b] = x[i];
            flags[i] = 0;
            x[i] = 0;
        }
        Ci[k+1] = Ci[k]+p;
    }
    return C;
}

numeric.ccsLUPSolve = function ccsLUPSolve(LUP,B) {
    var L = LUP.L, U = LUP.U, P = LUP.P;
    var Bi = B[0];
    var flag = false;
    if(typeof Bi !== "object") { B = [[0,B.length],numeric.linspace(0,B.length-1),B]; Bi = B[0]; flag = true; }
    var Bj = B[1], Bv = B[2];
    var n = L[0].length-1, m = Bi.length-1;
    var x = numeric.rep([n],0), xj = Array(n);
    var b = numeric.rep([n],0), bj = Array(n);
    var Xi = numeric.rep([m+1],0), Xj = [], Xv = [];
    var sol = numeric.ccsTSolve;
    var i,j,j0,j1,k,J,N=0;
    for(i=0;i<m;++i) {
        k = 0;
        j0 = Bi[i];
        j1 = Bi[i+1];
        for(j=j0;j<j1;++j) { 
            J = LUP.Pinv[Bj[j]];
            bj[k] = J;
            b[J] = Bv[j];
            ++k;
        }
        bj.length = k;
        sol(L,b,x,bj,xj);
        for(j=bj.length-1;j!==-1;--j) b[bj[j]] = 0;
        sol(U,x,b,xj,bj);
        if(flag) return b;
        for(j=xj.length-1;j!==-1;--j) x[xj[j]] = 0;
        for(j=bj.length-1;j!==-1;--j) {
            J = bj[j];
            Xj[N] = J;
            Xv[N] = b[J];
            b[J] = 0;
            ++N;
        }
        Xi[i+1] = N;
    }
    return [Xi,Xj,Xv];
}

numeric.ccsbinop = function ccsbinop(body,setup) {
    if(typeof setup === "undefined") setup='';
    return Function('X','Y',
            'var Xi = X[0], Xj = X[1], Xv = X[2];\n'+
            'var Yi = Y[0], Yj = Y[1], Yv = Y[2];\n'+
            'var n = Xi.length-1,m = Math.max(numeric.sup(Xj),numeric.sup(Yj))+1;\n'+
            'var Zi = numeric.rep([n+1],0), Zj = [], Zv = [];\n'+
            'var x = numeric.rep([m],0),y = numeric.rep([m],0);\n'+
            'var xk,yk,zk;\n'+
            'var i,j,j0,j1,k,p=0;\n'+
            setup+
            'for(i=0;i<n;++i) {\n'+
            '  j0 = Xi[i]; j1 = Xi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) {\n'+
            '    k = Xj[j];\n'+
            '    x[k] = 1;\n'+
            '    Zj[p] = k;\n'+
            '    ++p;\n'+
            '  }\n'+
            '  j0 = Yi[i]; j1 = Yi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) {\n'+
            '    k = Yj[j];\n'+
            '    y[k] = Yv[j];\n'+
            '    if(x[k] === 0) {\n'+
            '      Zj[p] = k;\n'+
            '      ++p;\n'+
            '    }\n'+
            '  }\n'+
            '  Zi[i+1] = p;\n'+
            '  j0 = Xi[i]; j1 = Xi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) x[Xj[j]] = Xv[j];\n'+
            '  j0 = Zi[i]; j1 = Zi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) {\n'+
            '    k = Zj[j];\n'+
            '    xk = x[k];\n'+
            '    yk = y[k];\n'+
            body+'\n'+
            '    Zv[j] = zk;\n'+
            '  }\n'+
            '  j0 = Xi[i]; j1 = Xi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) x[Xj[j]] = 0;\n'+
            '  j0 = Yi[i]; j1 = Yi[i+1];\n'+
            '  for(j=j0;j!==j1;++j) y[Yj[j]] = 0;\n'+
            '}\n'+
            'return [Zi,Zj,Zv];'
            );
};

(function() {
    var k,A,B,C;
    for(k in numeric.ops2) {
        if(isFinite(eval('1'+numeric.ops2[k]+'0'))) A = '[Y[0],Y[1],numeric.'+k+'(X,Y[2])]';
        else A = 'NaN';
        if(isFinite(eval('0'+numeric.ops2[k]+'1'))) B = '[X[0],X[1],numeric.'+k+'(X[2],Y)]';
        else B = 'NaN';
        if(isFinite(eval('1'+numeric.ops2[k]+'0')) && isFinite(eval('0'+numeric.ops2[k]+'1'))) C = 'numeric.ccs'+k+'MM(X,Y)';
        else C = 'NaN';
        numeric['ccs'+k+'MM'] = numeric.ccsbinop('zk = xk '+numeric.ops2[k]+'yk;');
        numeric['ccs'+k] = Function('X','Y',
                'if(typeof X === "number") return '+A+';\n'+
                'if(typeof Y === "number") return '+B+';\n'+
                'return '+C+';\n'
                );
    }
}());

numeric.ccsScatter = function ccsScatter(A) {
    var Ai = A[0], Aj = A[1], Av = A[2];
    var n = numeric.sup(Aj)+1,m=Ai.length;
    var Ri = numeric.rep([n],0),Rj=Array(m), Rv = Array(m);
    var counts = numeric.rep([n],0),i;
    for(i=0;i<m;++i) counts[Aj[i]]++;
    for(i=0;i<n;++i) Ri[i+1] = Ri[i] + counts[i];
    var ptr = Ri.slice(0),k,Aii;
    for(i=0;i<m;++i) {
        Aii = Aj[i];
        k = ptr[Aii];
        Rj[k] = Ai[i];
        Rv[k] = Av[i];
        ptr[Aii]=ptr[Aii]+1;
    }
    return [Ri,Rj,Rv];
}

numeric.ccsGather = function ccsGather(A) {
    var Ai = A[0], Aj = A[1], Av = A[2];
    var n = Ai.length-1,m = Aj.length;
    var Ri = Array(m), Rj = Array(m), Rv = Array(m);
    var i,j,j0,j1,p;
    p=0;
    for(i=0;i<n;++i) {
        j0 = Ai[i];
        j1 = Ai[i+1];
        for(j=j0;j!==j1;++j) {
            Rj[p] = i;
            Ri[p] = Aj[j];
            Rv[p] = Av[j];
            ++p;
        }
    }
    return [Ri,Rj,Rv];
}

// The following sparse linear algebra routines are deprecated.

numeric.sdim = function dim(A,ret,k) {
    if(typeof ret === "undefined") { ret = []; }
    if(typeof A !== "object") return ret;
    if(typeof k === "undefined") { k=0; }
    if(!(k in ret)) { ret[k] = 0; }
    if(A.length > ret[k]) ret[k] = A.length;
    var i;
    for(i in A) {
        if(A.hasOwnProperty(i)) dim(A[i],ret,k+1);
    }
    return ret;
};

numeric.sclone = function clone(A,k,n) {
    if(typeof k === "undefined") { k=0; }
    if(typeof n === "undefined") { n = numeric.sdim(A).length; }
    var i,ret = Array(A.length);
    if(k === n-1) {
        for(i in A) { if(A.hasOwnProperty(i)) ret[i] = A[i]; }
        return ret;
    }
    for(i in A) {
        if(A.hasOwnProperty(i)) ret[i] = clone(A[i],k+1,n);
    }
    return ret;
}

numeric.sdiag = function diag(d) {
    var n = d.length,i,ret = Array(n),i1,i2,i3;
    for(i=n-1;i>=1;i-=2) {
        i1 = i-1;
        ret[i] = []; ret[i][i] = d[i];
        ret[i1] = []; ret[i1][i1] = d[i1];
    }
    if(i===0) { ret[0] = []; ret[0][0] = d[i]; }
    return ret;
}

numeric.sidentity = function identity(n) { return numeric.sdiag(numeric.rep([n],1)); }

numeric.stranspose = function transpose(A) {
    var ret = [], n = A.length, i,j,Ai;
    for(i in A) {
        if(!(A.hasOwnProperty(i))) continue;
        Ai = A[i];
        for(j in Ai) {
            if(!(Ai.hasOwnProperty(j))) continue;
            if(typeof ret[j] !== "object") { ret[j] = []; }
            ret[j][i] = Ai[j];
        }
    }
    return ret;
}

numeric.sLUP = function LUP(A,tol) {
    throw new Error("The function numeric.sLUP had a bug in it and has been removed. Please use the new numeric.ccsLUP function instead.");
};

numeric.sdotMM = function dotMM(A,B) {
    var p = A.length, q = B.length, BT = numeric.stranspose(B), r = BT.length, Ai, BTk;
    var i,j,k,accum;
    var ret = Array(p),reti;
    for(i=p-1;i>=0;i--) {
        reti = [];
        Ai = A[i];
        for(k=r-1;k>=0;k--) {
            accum = 0;
            BTk = BT[k];
            for(j in Ai) {
                if(!(Ai.hasOwnProperty(j))) continue;
                if(j in BTk) { accum += Ai[j]*BTk[j]; }
            }
            if(accum) reti[k] = accum;
        }
        ret[i] = reti;
    }
    return ret;
}

numeric.sdotMV = function dotMV(A,x) {
    var p = A.length, Ai, i,j;
    var ret = Array(p), accum;
    for(i=p-1;i>=0;i--) {
        Ai = A[i];
        accum = 0;
        for(j in Ai) {
            if(!(Ai.hasOwnProperty(j))) continue;
            if(x[j]) accum += Ai[j]*x[j];
        }
        if(accum) ret[i] = accum;
    }
    return ret;
}

numeric.sdotVM = function dotMV(x,A) {
    var i,j,Ai,alpha;
    var ret = [], accum;
    for(i in x) {
        if(!x.hasOwnProperty(i)) continue;
        Ai = A[i];
        alpha = x[i];
        for(j in Ai) {
            if(!Ai.hasOwnProperty(j)) continue;
            if(!ret[j]) { ret[j] = 0; }
            ret[j] += alpha*Ai[j];
        }
    }
    return ret;
}

numeric.sdotVV = function dotVV(x,y) {
    var i,ret=0;
    for(i in x) { if(x[i] && y[i]) ret+= x[i]*y[i]; }
    return ret;
}

numeric.sdot = function dot(A,B) {
    var m = numeric.sdim(A).length, n = numeric.sdim(B).length;
    var k = m*1000+n;
    switch(k) {
    case 0: return A*B;
    case 1001: return numeric.sdotVV(A,B);
    case 2001: return numeric.sdotMV(A,B);
    case 1002: return numeric.sdotVM(A,B);
    case 2002: return numeric.sdotMM(A,B);
    default: throw new Error('numeric.sdot not implemented for tensors of order '+m+' and '+n);
    }
}

numeric.sscatter = function scatter(V) {
    var n = V[0].length, Vij, i, j, m = V.length, A = [], Aj;
    for(i=n-1;i>=0;--i) {
        if(!V[m-1][i]) continue;
        Aj = A;
        for(j=0;j<m-2;j++) {
            Vij = V[j][i];
            if(!Aj[Vij]) Aj[Vij] = [];
            Aj = Aj[Vij];
        }
        Aj[V[j][i]] = V[j+1][i];
    }
    return A;
}

numeric.sgather = function gather(A,ret,k) {
    if(typeof ret === "undefined") ret = [];
    if(typeof k === "undefined") k = [];
    var n,i,Ai;
    n = k.length;
    for(i in A) {
        if(A.hasOwnProperty(i)) {
            k[n] = parseInt(i);
            Ai = A[i];
            if(typeof Ai === "number") {
                if(Ai) {
                    if(ret.length === 0) {
                        for(i=n+1;i>=0;--i) ret[i] = [];
                    }
                    for(i=n;i>=0;--i) ret[i].push(k[i]);
                    ret[n+1].push(Ai);
                }
            } else gather(Ai,ret,k);
        }
    }
    if(k.length>n) k.pop();
    return ret;
}

// 6. Coordinate matrices
numeric.cLU = function LU(A) {
    var I = A[0], J = A[1], V = A[2];
    var p = I.length, m=0, i,j,k,a,b,c;
    for(i=0;i<p;i++) if(I[i]>m) m=I[i];
    m++;
    var L = Array(m), U = Array(m), left = numeric.rep([m],Infinity), right = numeric.rep([m],-Infinity);
    var Ui, Uj,alpha;
    for(k=0;k<p;k++) {
        i = I[k];
        j = J[k];
        if(j<left[i]) left[i] = j;
        if(j>right[i]) right[i] = j;
    }
    for(i=0;i<m-1;i++) { if(right[i] > right[i+1]) right[i+1] = right[i]; }
    for(i=m-1;i>=1;i--) { if(left[i]<left[i-1]) left[i-1] = left[i]; }
    var countL = 0, countU = 0;
    for(i=0;i<m;i++) {
        U[i] = numeric.rep([right[i]-left[i]+1],0);
        L[i] = numeric.rep([i-left[i]],0);
        countL += i-left[i]+1;
        countU += right[i]-i+1;
    }
    for(k=0;k<p;k++) { i = I[k]; U[i][J[k]-left[i]] = V[k]; }
    for(i=0;i<m-1;i++) {
        a = i-left[i];
        Ui = U[i];
        for(j=i+1;left[j]<=i && j<m;j++) {
            b = i-left[j];
            c = right[i]-i;
            Uj = U[j];
            alpha = Uj[b]/Ui[a];
            if(alpha) {
                for(k=1;k<=c;k++) { Uj[k+b] -= alpha*Ui[k+a]; }
                L[j][i-left[j]] = alpha;
            }
        }
    }
    var Ui = [], Uj = [], Uv = [], Li = [], Lj = [], Lv = [];
    var p,q,foo;
    p=0; q=0;
    for(i=0;i<m;i++) {
        a = left[i];
        b = right[i];
        foo = U[i];
        for(j=i;j<=b;j++) {
            if(foo[j-a]) {
                Ui[p] = i;
                Uj[p] = j;
                Uv[p] = foo[j-a];
                p++;
            }
        }
        foo = L[i];
        for(j=a;j<i;j++) {
            if(foo[j-a]) {
                Li[q] = i;
                Lj[q] = j;
                Lv[q] = foo[j-a];
                q++;
            }
        }
        Li[q] = i;
        Lj[q] = i;
        Lv[q] = 1;
        q++;
    }
    return {U:[Ui,Uj,Uv], L:[Li,Lj,Lv]};
};

numeric.cLUsolve = function LUsolve(lu,b) {
    var L = lu.L, U = lu.U, ret = numeric.clone(b);
    var Li = L[0], Lj = L[1], Lv = L[2];
    var Ui = U[0], Uj = U[1], Uv = U[2];
    var p = Ui.length, q = Li.length;
    var m = ret.length,i,j,k;
    k = 0;
    for(i=0;i<m;i++) {
        while(Lj[k] < i) {
            ret[i] -= Lv[k]*ret[Lj[k]];
            k++;
        }
        k++;
    }
    k = p-1;
    for(i=m-1;i>=0;i--) {
        while(Uj[k] > i) {
            ret[i] -= Uv[k]*ret[Uj[k]];
            k--;
        }
        ret[i] /= Uv[k];
        k--;
    }
    return ret;
};

numeric.cgrid = function grid(n,shape) {
    if(typeof n === "number") n = [n,n];
    var ret = numeric.rep(n,-1);
    var i,j,count;
    if(typeof shape !== "function") {
        switch(shape) {
        case 'L':
            shape = function(i,j) { return (i>=n[0]/2 || j<n[1]/2); }
            break;
        default:
            shape = function(i,j) { return true; };
            break;
        }
    }
    count=0;
    for(i=1;i<n[0]-1;i++) for(j=1;j<n[1]-1;j++) 
        if(shape(i,j)) {
            ret[i][j] = count;
            count++;
        }
    return ret;
}

numeric.cdelsq = function delsq(g) {
    var dir = [[-1,0],[0,-1],[0,1],[1,0]];
    var s = numeric.dim(g), m = s[0], n = s[1], i,j,k,p,q;
    var Li = [], Lj = [], Lv = [];
    for(i=1;i<m-1;i++) for(j=1;j<n-1;j++) {
        if(g[i][j]<0) continue;
        for(k=0;k<4;k++) {
            p = i+dir[k][0];
            q = j+dir[k][1];
            if(g[p][q]<0) continue;
            Li.push(g[i][j]);
            Lj.push(g[p][q]);
            Lv.push(-1);
        }
        Li.push(g[i][j]);
        Lj.push(g[i][j]);
        Lv.push(4);
    }
    return [Li,Lj,Lv];
}

numeric.cdotMV = function dotMV(A,x) {
    var ret, Ai = A[0], Aj = A[1], Av = A[2],k,p=Ai.length,N;
    N=0;
    for(k=0;k<p;k++) { if(Ai[k]>N) N = Ai[k]; }
    N++;
    ret = numeric.rep([N],0);
    for(k=0;k<p;k++) { ret[Ai[k]]+=Av[k]*x[Aj[k]]; }
    return ret;
}

// 7. Splines

numeric.Spline = function Spline(x,yl,yr,kl,kr) { this.x = x; this.yl = yl; this.yr = yr; this.kl = kl; this.kr = kr; }
numeric.Spline.prototype._at = function _at(x1,p) {
    var x = this.x;
    var yl = this.yl;
    var yr = this.yr;
    var kl = this.kl;
    var kr = this.kr;
    var x1,a,b,t;
    var add = numeric.add, sub = numeric.sub, mul = numeric.mul;
    a = sub(mul(kl[p],x[p+1]-x[p]),sub(yr[p+1],yl[p]));
    b = add(mul(kr[p+1],x[p]-x[p+1]),sub(yr[p+1],yl[p]));
    t = (x1-x[p])/(x[p+1]-x[p]);
    var s = t*(1-t);
    return add(add(add(mul(1-t,yl[p]),mul(t,yr[p+1])),mul(a,s*(1-t))),mul(b,s*t));
}
numeric.Spline.prototype.at = function at(x0) {
    if(typeof x0 === "number") {
        var x = this.x;
        var n = x.length;
        var p,q,mid,floor = Math.floor,a,b,t;
        p = 0;
        q = n-1;
        while(q-p>1) {
            mid = floor((p+q)/2);
            if(x[mid] <= x0) p = mid;
            else q = mid;
        }
        return this._at(x0,p);
    }
    var n = x0.length, i, ret = Array(n);
    for(i=n-1;i!==-1;--i) ret[i] = this.at(x0[i]);
    return ret;
}
numeric.Spline.prototype.diff = function diff() {
    var x = this.x;
    var yl = this.yl;
    var yr = this.yr;
    var kl = this.kl;
    var kr = this.kr;
    var n = yl.length;
    var i,dx,dy;
    var zl = kl, zr = kr, pl = Array(n), pr = Array(n);
    var add = numeric.add, mul = numeric.mul, div = numeric.div, sub = numeric.sub;
    for(i=n-1;i!==-1;--i) {
        dx = x[i+1]-x[i];
        dy = sub(yr[i+1],yl[i]);
        pl[i] = div(add(mul(dy, 6),mul(kl[i],-4*dx),mul(kr[i+1],-2*dx)),dx*dx);
        pr[i+1] = div(add(mul(dy,-6),mul(kl[i], 2*dx),mul(kr[i+1], 4*dx)),dx*dx);
    }
    return new numeric.Spline(x,zl,zr,pl,pr);
}
numeric.Spline.prototype.roots = function roots() {
    function sqr(x) { return x*x; }
    function heval(y0,y1,k0,k1,x) {
        var A = k0*2-(y1-y0);
        var B = -k1*2+(y1-y0);
        var t = (x+1)*0.5;
        var s = t*(1-t);
        return (1-t)*y0+t*y1+A*s*(1-t)+B*s*t;
    }
    var ret = [];
    var x = this.x, yl = this.yl, yr = this.yr, kl = this.kl, kr = this.kr;
    if(typeof yl[0] === "number") {
        yl = [yl];
        yr = [yr];
        kl = [kl];
        kr = [kr];
    }
    var m = yl.length,n=x.length-1,i,j,k,y,s,t;
    var ai,bi,ci,di, ret = Array(m),ri,k0,k1,y0,y1,A,B,D,dx,cx,stops,z0,z1,zm,t0,t1,tm;
    var sqrt = Math.sqrt;
    for(i=0;i!==m;++i) {
        ai = yl[i];
        bi = yr[i];
        ci = kl[i];
        di = kr[i];
        ri = [];
        for(j=0;j!==n;j++) {
            if(j>0 && bi[j]*ai[j]<0) ri.push(x[j]);
            dx = (x[j+1]-x[j]);
            cx = x[j];
            y0 = ai[j];
            y1 = bi[j+1];
            k0 = ci[j]/dx;
            k1 = di[j+1]/dx;
            D = sqr(k0-k1+3*(y0-y1)) + 12*k1*y0;
            A = k1+3*y0+2*k0-3*y1;
            B = 3*(k1+k0+2*(y0-y1));
            if(D<=0) {
                z0 = A/B;
                if(z0>x[j] && z0<x[j+1]) stops = [x[j],z0,x[j+1]];
                else stops = [x[j],x[j+1]];
            } else {
                z0 = (A-sqrt(D))/B;
                z1 = (A+sqrt(D))/B;
                stops = [x[j]];
                if(z0>x[j] && z0<x[j+1]) stops.push(z0);
                if(z1>x[j] && z1<x[j+1]) stops.push(z1);
                stops.push(x[j+1]);
            }
            t0 = stops[0];
            z0 = this._at(t0,j);
            for(k=0;k<stops.length-1;k++) {
                t1 = stops[k+1];
                z1 = this._at(t1,j);
                if(z0 === 0) {
                    ri.push(t0); 
                    t0 = t1;
                    z0 = z1;
                    continue;
                }
                if(z1 === 0 || z0*z1>0) {
                    t0 = t1;
                    z0 = z1;
                    continue;
                }
                var side = 0;
                while(1) {
                    tm = (z0*t1-z1*t0)/(z0-z1);
                    if(tm <= t0 || tm >= t1) { break; }
                    zm = this._at(tm,j);
                    if(zm*z1>0) {
                        t1 = tm;
                        z1 = zm;
                        if(side === -1) z0*=0.5;
                        side = -1;
                    } else if(zm*z0>0) {
                        t0 = tm;
                        z0 = zm;
                        if(side === 1) z1*=0.5;
                        side = 1;
                    } else break;
                }
                ri.push(tm);
                t0 = stops[k+1];
                z0 = this._at(t0, j);
            }
            if(z1 === 0) ri.push(t1);
        }
        ret[i] = ri;
    }
    if(typeof this.yl[0] === "number") return ret[0];
    return ret;
}
numeric.spline = function spline(x,y,k1,kn) {
    var n = x.length, b = [], dx = [], dy = [];
    var i;
    var sub = numeric.sub,mul = numeric.mul,add = numeric.add;
    for(i=n-2;i>=0;i--) { dx[i] = x[i+1]-x[i]; dy[i] = sub(y[i+1],y[i]); }
    if(typeof k1 === "string" || typeof kn === "string") { 
        k1 = kn = "periodic";
    }
    // Build sparse tridiagonal system
    var T = [[],[],[]];
    switch(typeof k1) {
    case "undefined":
        b[0] = mul(3/(dx[0]*dx[0]),dy[0]);
        T[0].push(0,0);
        T[1].push(0,1);
        T[2].push(2/dx[0],1/dx[0]);
        break;
    case "string":
        b[0] = add(mul(3/(dx[n-2]*dx[n-2]),dy[n-2]),mul(3/(dx[0]*dx[0]),dy[0]));
        T[0].push(0,0,0);
        T[1].push(n-2,0,1);
        T[2].push(1/dx[n-2],2/dx[n-2]+2/dx[0],1/dx[0]);
        break;
    default:
        b[0] = k1;
        T[0].push(0);
        T[1].push(0);
        T[2].push(1);
        break;
    }
    for(i=1;i<n-1;i++) {
        b[i] = add(mul(3/(dx[i-1]*dx[i-1]),dy[i-1]),mul(3/(dx[i]*dx[i]),dy[i]));
        T[0].push(i,i,i);
        T[1].push(i-1,i,i+1);
        T[2].push(1/dx[i-1],2/dx[i-1]+2/dx[i],1/dx[i]);
    }
    switch(typeof kn) {
    case "undefined":
        b[n-1] = mul(3/(dx[n-2]*dx[n-2]),dy[n-2]);
        T[0].push(n-1,n-1);
        T[1].push(n-2,n-1);
        T[2].push(1/dx[n-2],2/dx[n-2]);
        break;
    case "string":
        T[1][T[1].length-1] = 0;
        break;
    default:
        b[n-1] = kn;
        T[0].push(n-1);
        T[1].push(n-1);
        T[2].push(1);
        break;
    }
    if(typeof b[0] !== "number") b = numeric.transpose(b);
    else b = [b];
    var k = Array(b.length);
    if(typeof k1 === "string") {
        for(i=k.length-1;i!==-1;--i) {
            k[i] = numeric.ccsLUPSolve(numeric.ccsLUP(numeric.ccsScatter(T)),b[i]);
            k[i][n-1] = k[i][0];
        }
    } else {
        for(i=k.length-1;i!==-1;--i) {
            k[i] = numeric.cLUsolve(numeric.cLU(T),b[i]);
        }
    }
    if(typeof y[0] === "number") k = k[0];
    else k = numeric.transpose(k);
    return new numeric.Spline(x,y,y,k,k);
}

// 8. FFT
numeric.fftpow2 = function fftpow2(x,y) {
    var n = x.length;
    if(n === 1) return;
    var cos = Math.cos, sin = Math.sin, i,j;
    var xe = Array(n/2), ye = Array(n/2), xo = Array(n/2), yo = Array(n/2);
    j = n/2;
    for(i=n-1;i!==-1;--i) {
        --j;
        xo[j] = x[i];
        yo[j] = y[i];
        --i;
        xe[j] = x[i];
        ye[j] = y[i];
    }
    fftpow2(xe,ye);
    fftpow2(xo,yo);
    j = n/2;
    var t,k = (-6.2831853071795864769252867665590057683943387987502116419/n),ci,si;
    for(i=n-1;i!==-1;--i) {
        --j;
        if(j === -1) j = n/2-1;
        t = k*i;
        ci = cos(t);
        si = sin(t);
        x[i] = xe[j] + ci*xo[j] - si*yo[j];
        y[i] = ye[j] + ci*yo[j] + si*xo[j];
    }
}
numeric._ifftpow2 = function _ifftpow2(x,y) {
    var n = x.length;
    if(n === 1) return;
    var cos = Math.cos, sin = Math.sin, i,j;
    var xe = Array(n/2), ye = Array(n/2), xo = Array(n/2), yo = Array(n/2);
    j = n/2;
    for(i=n-1;i!==-1;--i) {
        --j;
        xo[j] = x[i];
        yo[j] = y[i];
        --i;
        xe[j] = x[i];
        ye[j] = y[i];
    }
    _ifftpow2(xe,ye);
    _ifftpow2(xo,yo);
    j = n/2;
    var t,k = (6.2831853071795864769252867665590057683943387987502116419/n),ci,si;
    for(i=n-1;i!==-1;--i) {
        --j;
        if(j === -1) j = n/2-1;
        t = k*i;
        ci = cos(t);
        si = sin(t);
        x[i] = xe[j] + ci*xo[j] - si*yo[j];
        y[i] = ye[j] + ci*yo[j] + si*xo[j];
    }
}
numeric.ifftpow2 = function ifftpow2(x,y) {
    numeric._ifftpow2(x,y);
    numeric.diveq(x,x.length);
    numeric.diveq(y,y.length);
}
numeric.convpow2 = function convpow2(ax,ay,bx,by) {
    numeric.fftpow2(ax,ay);
    numeric.fftpow2(bx,by);
    var i,n = ax.length,axi,bxi,ayi,byi;
    for(i=n-1;i!==-1;--i) {
        axi = ax[i]; ayi = ay[i]; bxi = bx[i]; byi = by[i];
        ax[i] = axi*bxi-ayi*byi;
        ay[i] = axi*byi+ayi*bxi;
    }
    numeric.ifftpow2(ax,ay);
}
numeric.T.prototype.fft = function fft() {
    var x = this.x, y = this.y;
    var n = x.length, log = Math.log, log2 = log(2),
        p = Math.ceil(log(2*n-1)/log2), m = Math.pow(2,p);
    var cx = numeric.rep([m],0), cy = numeric.rep([m],0), cos = Math.cos, sin = Math.sin;
    var k, c = (-3.141592653589793238462643383279502884197169399375105820/n),t;
    var a = numeric.rep([m],0), b = numeric.rep([m],0),nhalf = Math.floor(n/2);
    for(k=0;k<n;k++) a[k] = x[k];
    if(typeof y !== "undefined") for(k=0;k<n;k++) b[k] = y[k];
    cx[0] = 1;
    for(k=1;k<=m/2;k++) {
        t = c*k*k;
        cx[k] = cos(t);
        cy[k] = sin(t);
        cx[m-k] = cos(t);
        cy[m-k] = sin(t)
    }
    var X = new numeric.T(a,b), Y = new numeric.T(cx,cy);
    X = X.mul(Y);
    numeric.convpow2(X.x,X.y,numeric.clone(Y.x),numeric.neg(Y.y));
    X = X.mul(Y);
    X.x.length = n;
    X.y.length = n;
    return X;
}
numeric.T.prototype.ifft = function ifft() {
    var x = this.x, y = this.y;
    var n = x.length, log = Math.log, log2 = log(2),
        p = Math.ceil(log(2*n-1)/log2), m = Math.pow(2,p);
    var cx = numeric.rep([m],0), cy = numeric.rep([m],0), cos = Math.cos, sin = Math.sin;
    var k, c = (3.141592653589793238462643383279502884197169399375105820/n),t;
    var a = numeric.rep([m],0), b = numeric.rep([m],0),nhalf = Math.floor(n/2);
    for(k=0;k<n;k++) a[k] = x[k];
    if(typeof y !== "undefined") for(k=0;k<n;k++) b[k] = y[k];
    cx[0] = 1;
    for(k=1;k<=m/2;k++) {
        t = c*k*k;
        cx[k] = cos(t);
        cy[k] = sin(t);
        cx[m-k] = cos(t);
        cy[m-k] = sin(t)
    }
    var X = new numeric.T(a,b), Y = new numeric.T(cx,cy);
    X = X.mul(Y);
    numeric.convpow2(X.x,X.y,numeric.clone(Y.x),numeric.neg(Y.y));
    X = X.mul(Y);
    X.x.length = n;
    X.y.length = n;
    return X.div(n);
}

//9. Unconstrained optimization
numeric.gradient = function gradient(f,x) {
    var n = x.length;
    var f0 = f(x);
    if(isNaN(f0)) throw new Error('gradient: f(x) is a NaN!');
    var max = Math.max;
    var i,x0 = numeric.clone(x),f1,f2, J = Array(n);
    var div = numeric.div, sub = numeric.sub,errest,roundoff,max = Math.max,eps = 1e-3,abs = Math.abs, min = Math.min;
    var t0,t1,t2,it=0,d1,d2,N;
    for(i=0;i<n;i++) {
        var h = max(1e-6*f0,1e-8);
        while(1) {
            ++it;
            if(it>20) { throw new Error("Numerical gradient fails"); }
            x0[i] = x[i]+h;
            f1 = f(x0);
            x0[i] = x[i]-h;
            f2 = f(x0);
            x0[i] = x[i];
            if(isNaN(f1) || isNaN(f2)) { h/=16; continue; }
            J[i] = (f1-f2)/(2*h);
            t0 = x[i]-h;
            t1 = x[i];
            t2 = x[i]+h;
            d1 = (f1-f0)/h;
            d2 = (f0-f2)/h;
            N = max(abs(J[i]),abs(f0),abs(f1),abs(f2),abs(t0),abs(t1),abs(t2),1e-8);
            errest = min(max(abs(d1-J[i]),abs(d2-J[i]),abs(d1-d2))/N,h/N);
            if(errest>eps) { h/=16; }
            else break;
            }
    }
    return J;
}

numeric.uncmin = function uncmin(f,x0,tol,gradient,maxit,callback,options) {
    var grad = numeric.gradient;
    if(typeof options === "undefined") { options = {}; }
    if(typeof tol === "undefined") { tol = 1e-8; }
    if(typeof gradient === "undefined") { gradient = function(x) { return grad(f,x); }; }
    if(typeof maxit === "undefined") maxit = 1000;
    x0 = numeric.clone(x0);
    var n = x0.length;
    var f0 = f(x0),f1,df0;
    if(isNaN(f0)) throw new Error('uncmin: f(x0) is a NaN!');
    var max = Math.max, norm2 = numeric.norm2;
    tol = max(tol,numeric.epsilon);
    var step,g0,g1,H1 = options.Hinv || numeric.identity(n);
    var dot = numeric.dot, inv = numeric.inv, sub = numeric.sub, add = numeric.add, ten = numeric.tensor, div = numeric.div, mul = numeric.mul;
    var all = numeric.all, isfinite = numeric.isFinite, neg = numeric.neg;
    var it=0,i,s,x1,y,Hy,Hs,ys,i0,t,nstep,t1,t2;
    var msg = "";
    g0 = gradient(x0);
    while(it<maxit) {
        if(typeof callback === "function") { if(callback(it,x0,f0,g0,H1)) { msg = "Callback returned true"; break; } }
        if(!all(isfinite(g0))) { msg = "Gradient has Infinity or NaN"; break; }
        step = neg(dot(H1,g0));
        if(!all(isfinite(step))) { msg = "Search direction has Infinity or NaN"; break; }
        nstep = norm2(step);
        if(nstep < tol) { msg="Newton step smaller than tol"; break; }
        t = 1;
        df0 = dot(g0,step);
        // line search
        x1 = x0;
        while(it < maxit) {
            if(t*nstep < tol) { break; }
            s = mul(step,t);
            x1 = add(x0,s);
            f1 = f(x1);
            if(f1-f0 >= 0.1*t*df0 || isNaN(f1)) {
                t *= 0.5;
                ++it;
                continue;
            }
            break;
        }
        if(t*nstep < tol) { msg = "Line search step size smaller than tol"; break; }
        if(it === maxit) { msg = "maxit reached during line search"; break; }
        g1 = gradient(x1);
        y = sub(g1,g0);
        ys = dot(y,s);
        Hy = dot(H1,y);
        H1 = sub(add(H1,
                mul(
                        (ys+dot(y,Hy))/(ys*ys),
                        ten(s,s)    )),
                div(add(ten(Hy,s),ten(s,Hy)),ys));
        x0 = x1;
        f0 = f1;
        g0 = g1;
        ++it;
    }
    return {solution: x0, f: f0, gradient: g0, invHessian: H1, iterations:it, message: msg};
}

// 10. Ode solver (Dormand-Prince)
numeric.Dopri = function Dopri(x,y,f,ymid,iterations,msg,events) {
    this.x = x;
    this.y = y;
    this.f = f;
    this.ymid = ymid;
    this.iterations = iterations;
    this.events = events;
    this.message = msg;
}
numeric.Dopri.prototype._at = function _at(xi,j) {
    function sqr(x) { return x*x; }
    var sol = this;
    var xs = sol.x;
    var ys = sol.y;
    var k1 = sol.f;
    var ymid = sol.ymid;
    var n = xs.length;
    var x0,x1,xh,y0,y1,yh,xi;
    var floor = Math.floor,h;
    var c = 0.5;
    var add = numeric.add, mul = numeric.mul,sub = numeric.sub, p,q,w;
    x0 = xs[j];
    x1 = xs[j+1];
    y0 = ys[j];
    y1 = ys[j+1];
    h  = x1-x0;
    xh = x0+c*h;
    yh = ymid[j];
    p = sub(k1[j  ],mul(y0,1/(x0-xh)+2/(x0-x1)));
    q = sub(k1[j+1],mul(y1,1/(x1-xh)+2/(x1-x0)));
    w = [sqr(xi - x1) * (xi - xh) / sqr(x0 - x1) / (x0 - xh),
         sqr(xi - x0) * sqr(xi - x1) / sqr(x0 - xh) / sqr(x1 - xh),
         sqr(xi - x0) * (xi - xh) / sqr(x1 - x0) / (x1 - xh),
         (xi - x0) * sqr(xi - x1) * (xi - xh) / sqr(x0-x1) / (x0 - xh),
         (xi - x1) * sqr(xi - x0) * (xi - xh) / sqr(x0-x1) / (x1 - xh)];
    return add(add(add(add(mul(y0,w[0]),
                           mul(yh,w[1])),
                           mul(y1,w[2])),
                           mul( p,w[3])),
                           mul( q,w[4]));
}
numeric.Dopri.prototype.at = function at(x) {
    var i,j,k,floor = Math.floor;
    if(typeof x !== "number") {
        var n = x.length, ret = Array(n);
        for(i=n-1;i!==-1;--i) {
            ret[i] = this.at(x[i]);
        }
        return ret;
    }
    var x0 = this.x;
    i = 0; j = x0.length-1;
    while(j-i>1) {
        k = floor(0.5*(i+j));
        if(x0[k] <= x) i = k;
        else j = k;
    }
    return this._at(x,i);
}

numeric.dopri = function dopri(x0,x1,y0,f,tol,maxit,event) {
    if(typeof tol === "undefined") { tol = 1e-6; }
    if(typeof maxit === "undefined") { maxit = 1000; }
    var xs = [x0], ys = [y0], k1 = [f(x0,y0)], k2,k3,k4,k5,k6,k7, ymid = [];
    var A2 = 1/5;
    var A3 = [3/40,9/40];
    var A4 = [44/45,-56/15,32/9];
    var A5 = [19372/6561,-25360/2187,64448/6561,-212/729];
    var A6 = [9017/3168,-355/33,46732/5247,49/176,-5103/18656];
    var b = [35/384,0,500/1113,125/192,-2187/6784,11/84];
    var bm = [0.5*6025192743/30085553152,
              0,
              0.5*51252292925/65400821598,
              0.5*-2691868925/45128329728,
              0.5*187940372067/1594534317056,
              0.5*-1776094331/19743644256,
              0.5*11237099/235043384];
    var c = [1/5,3/10,4/5,8/9,1,1];
    var e = [-71/57600,0,71/16695,-71/1920,17253/339200,-22/525,1/40];
    var i = 0,er,j;
    var h = (x1-x0)/10;
    var it = 0;
    var add = numeric.add, mul = numeric.mul, y1,erinf;
    var max = Math.max, min = Math.min, abs = Math.abs, norminf = numeric.norminf,pow = Math.pow;
    var any = numeric.any, lt = numeric.lt, and = numeric.and, sub = numeric.sub;
    var e0, e1, ev;
    var ret = new numeric.Dopri(xs,ys,k1,ymid,-1,"");
    if(typeof event === "function") e0 = event(x0,y0);
    while(x0<x1 && it<maxit) {
        ++it;
        if(x0+h>x1) h = x1-x0;
        k2 = f(x0+c[0]*h,                add(y0,mul(   A2*h,k1[i])));
        k3 = f(x0+c[1]*h,            add(add(y0,mul(A3[0]*h,k1[i])),mul(A3[1]*h,k2)));
        k4 = f(x0+c[2]*h,        add(add(add(y0,mul(A4[0]*h,k1[i])),mul(A4[1]*h,k2)),mul(A4[2]*h,k3)));
        k5 = f(x0+c[3]*h,    add(add(add(add(y0,mul(A5[0]*h,k1[i])),mul(A5[1]*h,k2)),mul(A5[2]*h,k3)),mul(A5[3]*h,k4)));
        k6 = f(x0+c[4]*h,add(add(add(add(add(y0,mul(A6[0]*h,k1[i])),mul(A6[1]*h,k2)),mul(A6[2]*h,k3)),mul(A6[3]*h,k4)),mul(A6[4]*h,k5)));
        y1 = add(add(add(add(add(y0,mul(k1[i],h*b[0])),mul(k3,h*b[2])),mul(k4,h*b[3])),mul(k5,h*b[4])),mul(k6,h*b[5]));
        k7 = f(x0+h,y1);
        er = add(add(add(add(add(mul(k1[i],h*e[0]),mul(k3,h*e[2])),mul(k4,h*e[3])),mul(k5,h*e[4])),mul(k6,h*e[5])),mul(k7,h*e[6]));
        if(typeof er === "number") erinf = abs(er);
        else erinf = norminf(er);
        if(erinf > tol) { // reject
            h = 0.2*h*pow(tol/erinf,0.25);
            if(x0+h === x0) {
                ret.msg = "Step size became too small";
                break;
            }
            continue;
        }
        ymid[i] = add(add(add(add(add(add(y0,
                mul(k1[i],h*bm[0])),
                mul(k3   ,h*bm[2])),
                mul(k4   ,h*bm[3])),
                mul(k5   ,h*bm[4])),
                mul(k6   ,h*bm[5])),
                mul(k7   ,h*bm[6]));
        ++i;
        xs[i] = x0+h;
        ys[i] = y1;
        k1[i] = k7;
        if(typeof event === "function") {
            var yi,xl = x0,xr = x0+0.5*h,xi;
            e1 = event(xr,ymid[i-1]);
            ev = and(lt(e0,0),lt(0,e1));
            if(!any(ev)) { xl = xr; xr = x0+h; e0 = e1; e1 = event(xr,y1); ev = and(lt(e0,0),lt(0,e1)); }
            if(any(ev)) {
                var xc, yc, en,ei;
                var side=0, sl = 1.0, sr = 1.0;
                while(1) {
                    if(typeof e0 === "number") xi = (sr*e1*xl-sl*e0*xr)/(sr*e1-sl*e0);
                    else {
                        xi = xr;
                        for(j=e0.length-1;j!==-1;--j) {
                            if(e0[j]<0 && e1[j]>0) xi = min(xi,(sr*e1[j]*xl-sl*e0[j]*xr)/(sr*e1[j]-sl*e0[j]));
                        }
                    }
                    if(xi <= xl || xi >= xr) break;
                    yi = ret._at(xi, i-1);
                    ei = event(xi,yi);
                    en = and(lt(e0,0),lt(0,ei));
                    if(any(en)) {
                        xr = xi;
                        e1 = ei;
                        ev = en;
                        sr = 1.0;
                        if(side === -1) sl *= 0.5;
                        else sl = 1.0;
                        side = -1;
                    } else {
                        xl = xi;
                        e0 = ei;
                        sl = 1.0;
                        if(side === 1) sr *= 0.5;
                        else sr = 1.0;
                        side = 1;
                    }
                }
                y1 = ret._at(0.5*(x0+xi),i-1);
                ret.f[i] = f(xi,yi);
                ret.x[i] = xi;
                ret.y[i] = yi;
                ret.ymid[i-1] = y1;
                ret.events = ev;
                ret.iterations = it;
                return ret;
            }
        }
        x0 += h;
        y0 = y1;
        e0 = e1;
        h = min(0.8*h*pow(tol/erinf,0.25),4*h);
    }
    ret.iterations = it;
    return ret;
}

// 11. Ax = b
numeric.LU = function(A, fast) {
  fast = fast || false;

  var abs = Math.abs;
  var i, j, k, absAjk, Akk, Ak, Pk, Ai;
  var max;
  var n = A.length, n1 = n-1;
  var P = new Array(n);
  if(!fast) A = numeric.clone(A);

  for (k = 0; k < n; ++k) {
    Pk = k;
    Ak = A[k];
    max = abs(Ak[k]);
    for (j = k + 1; j < n; ++j) {
      absAjk = abs(A[j][k]);
      if (max < absAjk) {
        max = absAjk;
        Pk = j;
      }
    }
    P[k] = Pk;

    if (Pk != k) {
      A[k] = A[Pk];
      A[Pk] = Ak;
      Ak = A[k];
    }

    Akk = Ak[k];

    for (i = k + 1; i < n; ++i) {
      A[i][k] /= Akk;
    }

    for (i = k + 1; i < n; ++i) {
      Ai = A[i];
      for (j = k + 1; j < n1; ++j) {
        Ai[j] -= Ai[k] * Ak[j];
        ++j;
        Ai[j] -= Ai[k] * Ak[j];
      }
      if(j===n1) Ai[j] -= Ai[k] * Ak[j];
    }
  }

  return {
    LU: A,
    P:  P
  };
}

numeric.LUsolve = function LUsolve(LUP, b) {
  var i, j;
  var LU = LUP.LU;
  var n   = LU.length;
  var x = numeric.clone(b);
  var P   = LUP.P;
  var Pi, LUi, LUii, tmp;

  for (i=n-1;i!==-1;--i) x[i] = b[i];
  for (i = 0; i < n; ++i) {
    Pi = P[i];
    if (P[i] !== i) {
      tmp = x[i];
      x[i] = x[Pi];
      x[Pi] = tmp;
    }

    LUi = LU[i];
    for (j = 0; j < i; ++j) {
      x[i] -= x[j] * LUi[j];
    }
  }

  for (i = n - 1; i >= 0; --i) {
    LUi = LU[i];
    for (j = i + 1; j < n; ++j) {
      x[i] -= x[j] * LUi[j];
    }

    x[i] /= LUi[i];
  }

  return x;
}

numeric.solve = function solve(A,b,fast) { return numeric.LUsolve(numeric.LU(A,fast), b); }

// 12. Linear programming
numeric.echelonize = function echelonize(A) {
    var s = numeric.dim(A), m = s[0], n = s[1];
    var I = numeric.identity(m);
    var P = Array(m);
    var i,j,k,l,Ai,Ii,Z,a;
    var abs = Math.abs;
    var diveq = numeric.diveq;
    A = numeric.clone(A);
    for(i=0;i<m;++i) {
        k = 0;
        Ai = A[i];
        Ii = I[i];
        for(j=1;j<n;++j) if(abs(Ai[k])<abs(Ai[j])) k=j;
        P[i] = k;
        diveq(Ii,Ai[k]);
        diveq(Ai,Ai[k]);
        for(j=0;j<m;++j) if(j!==i) {
            Z = A[j]; a = Z[k];
            for(l=n-1;l!==-1;--l) Z[l] -= Ai[l]*a;
            Z = I[j];
            for(l=m-1;l!==-1;--l) Z[l] -= Ii[l]*a;
        }
    }
    return {I:I, A:A, P:P};
}

numeric.__solveLP = function __solveLP(c,A,b,tol,maxit,x,flag) {
    var sum = numeric.sum, log = numeric.log, mul = numeric.mul, sub = numeric.sub, dot = numeric.dot, div = numeric.div, add = numeric.add;
    var m = c.length, n = b.length,y;
    var unbounded = false, cb,i0=0;
    var alpha = 1.0;
    var f0,df0,AT = numeric.transpose(A), svd = numeric.svd,transpose = numeric.transpose,leq = numeric.leq, sqrt = Math.sqrt, abs = Math.abs;
    var muleq = numeric.muleq;
    var norm = numeric.norminf, any = numeric.any,min = Math.min;
    var all = numeric.all, gt = numeric.gt;
    var p = Array(m), A0 = Array(n),e=numeric.rep([n],1), H;
    var solve = numeric.solve, z = sub(b,dot(A,x)),count;
    var dotcc = dot(c,c);
    var g;
    for(count=i0;count<maxit;++count) {
        var i,j,d;
        for(i=n-1;i!==-1;--i) A0[i] = div(A[i],z[i]);
        var A1 = transpose(A0);
        for(i=m-1;i!==-1;--i) p[i] = (/*x[i]+*/sum(A1[i]));
        alpha = 0.25*abs(dotcc/dot(c,p));
        var a1 = 100*sqrt(dotcc/dot(p,p));
        if(!isFinite(alpha) || alpha>a1) alpha = a1;
        g = add(c,mul(alpha,p));
        H = dot(A1,A0);
        for(i=m-1;i!==-1;--i) H[i][i] += 1;
        d = solve(H,div(g,alpha),true);
        var t0 = div(z,dot(A,d));
        var t = 1.0;
        for(i=n-1;i!==-1;--i) if(t0[i]<0) t = min(t,-0.999*t0[i]);
        y = sub(x,mul(d,t));
        z = sub(b,dot(A,y));
        if(!all(gt(z,0))) return { solution: x, message: "", iterations: count };
        x = y;
        if(alpha<tol) return { solution: y, message: "", iterations: count };
        if(flag) {
            var s = dot(c,g), Ag = dot(A,g);
            unbounded = true;
            for(i=n-1;i!==-1;--i) if(s*Ag[i]<0) { unbounded = false; break; }
        } else {
            if(x[m-1]>=0) unbounded = false;
            else unbounded = true;
        }
        if(unbounded) return { solution: y, message: "Unbounded", iterations: count };
    }
    return { solution: x, message: "maximum iteration count exceeded", iterations:count };
}

numeric._solveLP = function _solveLP(c,A,b,tol,maxit) {
    var m = c.length, n = b.length,y;
    var sum = numeric.sum, log = numeric.log, mul = numeric.mul, sub = numeric.sub, dot = numeric.dot, div = numeric.div, add = numeric.add;
    var c0 = numeric.rep([m],0).concat([1]);
    var J = numeric.rep([n,1],-1);
    var A0 = numeric.blockMatrix([[A                   ,   J  ]]);
    var b0 = b;
    var y = numeric.rep([m],0).concat(Math.max(0,numeric.sup(numeric.neg(b)))+1);
    var x0 = numeric.__solveLP(c0,A0,b0,tol,maxit,y,false);
    var x = numeric.clone(x0.solution);
    x.length = m;
    var foo = numeric.inf(sub(b,dot(A,x)));
    if(foo<0) { return { solution: NaN, message: "Infeasible", iterations: x0.iterations }; }
    var ret = numeric.__solveLP(c, A, b, tol, maxit-x0.iterations, x, true);
    ret.iterations += x0.iterations;
    return ret;
};

numeric.solveLP = function solveLP(c,A,b,Aeq,beq,tol,maxit) {
    if(typeof maxit === "undefined") maxit = 1000;
    if(typeof tol === "undefined") tol = numeric.epsilon;
    if(typeof Aeq === "undefined") return numeric._solveLP(c,A,b,tol,maxit);
    var m = Aeq.length, n = Aeq[0].length, o = A.length;
    var B = numeric.echelonize(Aeq);
    var flags = numeric.rep([n],0);
    var P = B.P;
    var Q = [];
    var i;
    for(i=P.length-1;i!==-1;--i) flags[P[i]] = 1;
    for(i=n-1;i!==-1;--i) if(flags[i]===0) Q.push(i);
    var g = numeric.getRange;
    var I = numeric.linspace(0,m-1), J = numeric.linspace(0,o-1);
    var Aeq2 = g(Aeq,I,Q), A1 = g(A,J,P), A2 = g(A,J,Q), dot = numeric.dot, sub = numeric.sub;
    var A3 = dot(A1,B.I);
    var A4 = sub(A2,dot(A3,Aeq2)), b4 = sub(b,dot(A3,beq));
    var c1 = Array(P.length), c2 = Array(Q.length);
    for(i=P.length-1;i!==-1;--i) c1[i] = c[P[i]];
    for(i=Q.length-1;i!==-1;--i) c2[i] = c[Q[i]];
    var c4 = sub(c2,dot(c1,dot(B.I,Aeq2)));
    var S = numeric._solveLP(c4,A4,b4,tol,maxit);
    var x2 = S.solution;
    if(x2!==x2) return S;
    var x1 = dot(B.I,sub(beq,dot(Aeq2,x2)));
    var x = Array(c.length);
    for(i=P.length-1;i!==-1;--i) x[P[i]] = x1[i];
    for(i=Q.length-1;i!==-1;--i) x[Q[i]] = x2[i];
    return { solution: x, message:S.message, iterations: S.iterations };
}

numeric.MPStoLP = function MPStoLP(MPS) {
    if(MPS instanceof String) { MPS.split('\n'); }
    var state = 0;
    var states = ['Initial state','NAME','ROWS','COLUMNS','RHS','BOUNDS','ENDATA'];
    var n = MPS.length;
    var i,j,z,N=0,rows = {}, sign = [], rl = 0, vars = {}, nv = 0;
    var name;
    var c = [], A = [], b = [];
    function err(e) { throw new Error('MPStoLP: '+e+'\nLine '+i+': '+MPS[i]+'\nCurrent state: '+states[state]+'\n'); }
    for(i=0;i<n;++i) {
        z = MPS[i];
        var w0 = z.match(/\S*/g);
        var w = [];
        for(j=0;j<w0.length;++j) if(w0[j]!=="") w.push(w0[j]);
        if(w.length === 0) continue;
        for(j=0;j<states.length;++j) if(z.substr(0,states[j].length) === states[j]) break;
        if(j<states.length) {
            state = j;
            if(j===1) { name = w[1]; }
            if(j===6) return { name:name, c:c, A:numeric.transpose(A), b:b, rows:rows, vars:vars };
            continue;
        }
        switch(state) {
        case 0: case 1: err('Unexpected line');
        case 2: 
            switch(w[0]) {
            case 'N': if(N===0) N = w[1]; else err('Two or more N rows'); break;
            case 'L': rows[w[1]] = rl; sign[rl] = 1; b[rl] = 0; ++rl; break;
            case 'G': rows[w[1]] = rl; sign[rl] = -1;b[rl] = 0; ++rl; break;
            case 'E': rows[w[1]] = rl; sign[rl] = 0;b[rl] = 0; ++rl; break;
            default: err('Parse error '+numeric.prettyPrint(w));
            }
            break;
        case 3:
            if(!vars.hasOwnProperty(w[0])) { vars[w[0]] = nv; c[nv] = 0; A[nv] = numeric.rep([rl],0); ++nv; }
            var p = vars[w[0]];
            for(j=1;j<w.length;j+=2) {
                if(w[j] === N) { c[p] = parseFloat(w[j+1]); continue; }
                var q = rows[w[j]];
                A[p][q] = (sign[q]<0?-1:1)*parseFloat(w[j+1]);
            }
            break;
        case 4:
            for(j=1;j<w.length;j+=2) b[rows[w[j]]] = (sign[rows[w[j]]]<0?-1:1)*parseFloat(w[j+1]);
            break;
        case 5: /*FIXME*/ break;
        case 6: err('Internal error');
        }
    }
    err('Reached end of file without ENDATA');
}
// seedrandom.js version 2.0.
// Author: David Bau 4/2/2011
//
// Defines a method Math.seedrandom() that, when called, substitutes
// an explicitly seeded RC4-based algorithm for Math.random().  Also
// supports automatic seeding from local or network sources of entropy.
//
// Usage:
//
//   <script src=http://davidbau.com/encode/seedrandom-min.js></script>
//
//   Math.seedrandom('yipee'); Sets Math.random to a function that is
//                             initialized using the given explicit seed.
//
//   Math.seedrandom();        Sets Math.random to a function that is
//                             seeded using the current time, dom state,
//                             and other accumulated local entropy.
//                             The generated seed string is returned.
//
//   Math.seedrandom('yowza', true);
//                             Seeds using the given explicit seed mixed
//                             together with accumulated entropy.
//
//   <script src="http://bit.ly/srandom-512"></script>
//                             Seeds using physical random bits downloaded
//                             from random.org.
//
//   <script src="https://jsonlib.appspot.com/urandom?callback=Math.seedrandom">
//   </script>                 Seeds using urandom bits from call.jsonlib.com,
//                             which is faster than random.org.
//
// Examples:
//
//   Math.seedrandom("hello");            // Use "hello" as the seed.
//   document.write(Math.random());       // Always 0.5463663768140734
//   document.write(Math.random());       // Always 0.43973793770592234
//   var rng1 = Math.random;              // Remember the current prng.
//
//   var autoseed = Math.seedrandom();    // New prng with an automatic seed.
//   document.write(Math.random());       // Pretty much unpredictable.
//
//   Math.random = rng1;                  // Continue "hello" prng sequence.
//   document.write(Math.random());       // Always 0.554769432473455
//
//   Math.seedrandom(autoseed);           // Restart at the previous seed.
//   document.write(Math.random());       // Repeat the 'unpredictable' value.
//
// Notes:
//
// Each time seedrandom('arg') is called, entropy from the passed seed
// is accumulated in a pool to help generate future seeds for the
// zero-argument form of Math.seedrandom, so entropy can be injected over
// time by calling seedrandom with explicit data repeatedly.
//
// On speed - This javascript implementation of Math.random() is about
// 3-10x slower than the built-in Math.random() because it is not native
// code, but this is typically fast enough anyway.  Seeding is more expensive,
// especially if you use auto-seeding.  Some details (timings on Chrome 4):
//
// Our Math.random()            - avg less than 0.002 milliseconds per call
// seedrandom('explicit')       - avg less than 0.5 milliseconds per call
// seedrandom('explicit', true) - avg less than 2 milliseconds per call
// seedrandom()                 - avg about 38 milliseconds per call
//
// LICENSE (BSD):
//
// Copyright 2010 David Bau, all rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
// 
//   3. Neither the name of this module nor the names of its contributors may
//      be used to endorse or promote products derived from this software
//      without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
/**
 * All code is in an anonymous closure to keep the global namespace clean.
 *
 * @param {number=} overflow 
 * @param {number=} startdenom
 */

// Patched by Seb so that seedrandom.js does not pollute the Math object.
// My tests suggest that doing Math.trouble = 1 makes Math lookups about 5%
// slower.
numeric.seedrandom = { pow:Math.pow, random:Math.random };

(function (pool, math, width, chunks, significance, overflow, startdenom) {


//
// seedrandom()
// This is the seedrandom function described above.
//
math['seedrandom'] = function seedrandom(seed, use_entropy) {
  var key = [];
  var arc4;

  // Flatten the seed string or build one from local entropy if needed.
  seed = mixkey(flatten(
    use_entropy ? [seed, pool] :
    arguments.length ? seed :
    [new Date().getTime(), pool, window], 3), key);

  // Use the seed to initialize an ARC4 generator.
  arc4 = new ARC4(key);

  // Mix the randomness into accumulated entropy.
  mixkey(arc4.S, pool);

  // Override Math.random

  // This function returns a random double in [0, 1) that contains
  // randomness in every bit of the mantissa of the IEEE 754 value.

  math['random'] = function random() {  // Closure to return a random double:
    var n = arc4.g(chunks);             // Start with a numerator n < 2 ^ 48
    var d = startdenom;                 //   and denominator d = 2 ^ 48.
    var x = 0;                          //   and no 'extra last byte'.
    while (n < significance) {          // Fill up all significant digits by
      n = (n + x) * width;              //   shifting numerator and
      d *= width;                       //   denominator and generating a
      x = arc4.g(1);                    //   new least-significant-byte.
    }
    while (n >= overflow) {             // To avoid rounding up, before adding
      n /= 2;                           //   last byte, shift everything
      d /= 2;                           //   right using integer math until
      x >>>= 1;                         //   we have exactly the desired bits.
    }
    return (n + x) / d;                 // Form the number within [0, 1).
  };

  // Return the seed that was used
  return seed;
};

//
// ARC4
//
// An ARC4 implementation.  The constructor takes a key in the form of
// an array of at most (width) integers that should be 0 <= x < (width).
//
// The g(count) method returns a pseudorandom integer that concatenates
// the next (count) outputs from ARC4.  Its return value is a number x
// that is in the range 0 <= x < (width ^ count).
//
/** @constructor */
function ARC4(key) {
  var t, u, me = this, keylen = key.length;
  var i = 0, j = me.i = me.j = me.m = 0;
  me.S = [];
  me.c = [];

  // The empty key [] is treated as [0].
  if (!keylen) { key = [keylen++]; }

  // Set up S using the standard key scheduling algorithm.
  while (i < width) { me.S[i] = i++; }
  for (i = 0; i < width; i++) {
    t = me.S[i];
    j = lowbits(j + t + key[i % keylen]);
    u = me.S[j];
    me.S[i] = u;
    me.S[j] = t;
  }

  // The "g" method returns the next (count) outputs as one number.
  me.g = function getnext(count) {
    var s = me.S;
    var i = lowbits(me.i + 1); var t = s[i];
    var j = lowbits(me.j + t); var u = s[j];
    s[i] = u;
    s[j] = t;
    var r = s[lowbits(t + u)];
    while (--count) {
      i = lowbits(i + 1); t = s[i];
      j = lowbits(j + t); u = s[j];
      s[i] = u;
      s[j] = t;
      r = r * width + s[lowbits(t + u)];
    }
    me.i = i;
    me.j = j;
    return r;
  };
  // For robust unpredictability discard an initial batch of values.
  // See http://www.rsa.com/rsalabs/node.asp?id=2009
  me.g(width);
}

//
// flatten()
// Converts an object tree to nested arrays of strings.
//
/** @param {Object=} result 
  * @param {string=} prop
  * @param {string=} typ */
function flatten(obj, depth, result, prop, typ) {
  result = [];
  typ = typeof(obj);
  if (depth && typ == 'object') {
    for (prop in obj) {
      if (prop.indexOf('S') < 5) {    // Avoid FF3 bug (local/sessionStorage)
        try { result.push(flatten(obj[prop], depth - 1)); } catch (e) {}
      }
    }
  }
  return (result.length ? result : obj + (typ != 'string' ? '\0' : ''));
}

//
// mixkey()
// Mixes a string seed into a key that is an array of integers, and
// returns a shortened string seed that is equivalent to the result key.
//
/** @param {number=} smear 
  * @param {number=} j */
function mixkey(seed, key, smear, j) {
  seed += '';                         // Ensure the seed is a string
  smear = 0;
  for (j = 0; j < seed.length; j++) {
    key[lowbits(j)] =
      lowbits((smear ^= key[lowbits(j)] * 19) + seed.charCodeAt(j));
  }
  seed = '';
  for (j in key) { seed += String.fromCharCode(key[j]); }
  return seed;
}

//
// lowbits()
// A quick "n mod width" for width a power of 2.
//
function lowbits(n) { return n & (width - 1); }

//
// The following constants are related to IEEE 754 limits.
//
startdenom = math.pow(width, chunks);
significance = math.pow(2, significance);
overflow = significance * 2;

//
// When seedrandom.js is loaded, we immediately mix a few bits
// from the built-in RNG into the entropy pool.  Because we do
// not want to intefere with determinstic PRNG state later,
// seedrandom will not call math.random on its own again after
// initialization.
//
mixkey(math.random(), pool);

// End anonymous scope, and pass initial values.
}(
  [],   // pool: entropy pool starts empty
  numeric.seedrandom, // math: package containing random, pow, and seedrandom
  256,  // width: each RC4 output is 0 <= x < 256
  6,    // chunks: at least six RC4 outputs for each double
  52    // significance: there are 52 significant digits in a double
  ));
/* This file is a slightly modified version of quadprog.js from Alberto Santini.
 * It has been slightly modified by SÃ©bastien Loisel to make sure that it handles
 * 0-based Arrays instead of 1-based Arrays.
 * License is in resources/LICENSE.quadprog */
(function(exports) {

function base0to1(A) {
    if(typeof A !== "object") { return A; }
    var ret = [], i,n=A.length;
    for(i=0;i<n;i++) ret[i+1] = base0to1(A[i]);
    return ret;
}
function base1to0(A) {
    if(typeof A !== "object") { return A; }
    var ret = [], i,n=A.length;
    for(i=1;i<n;i++) ret[i-1] = base1to0(A[i]);
    return ret;
}

function dpori(a, lda, n) {
    var i, j, k, kp1, t;

    for (k = 1; k <= n; k = k + 1) {
        a[k][k] = 1 / a[k][k];
        t = -a[k][k];
        //~ dscal(k - 1, t, a[1][k], 1);
        for (i = 1; i < k; i = i + 1) {
            a[i][k] = t * a[i][k];
        }

        kp1 = k + 1;
        if (n < kp1) {
            break;
        }
        for (j = kp1; j <= n; j = j + 1) {
            t = a[k][j];
            a[k][j] = 0;
            //~ daxpy(k, t, a[1][k], 1, a[1][j], 1);
            for (i = 1; i <= k; i = i + 1) {
                a[i][j] = a[i][j] + (t * a[i][k]);
            }
        }
    }

}

function dposl(a, lda, n, b) {
    var i, k, kb, t;

    for (k = 1; k <= n; k = k + 1) {
        //~ t = ddot(k - 1, a[1][k], 1, b[1], 1);
        t = 0;
        for (i = 1; i < k; i = i + 1) {
            t = t + (a[i][k] * b[i]);
        }

        b[k] = (b[k] - t) / a[k][k];
    }

    for (kb = 1; kb <= n; kb = kb + 1) {
        k = n + 1 - kb;
        b[k] = b[k] / a[k][k];
        t = -b[k];
        //~ daxpy(k - 1, t, a[1][k], 1, b[1], 1);
        for (i = 1; i < k; i = i + 1) {
            b[i] = b[i] + (t * a[i][k]);
        }
    }
}

function dpofa(a, lda, n, info) {
    var i, j, jm1, k, t, s;

    for (j = 1; j <= n; j = j + 1) {
        info[1] = j;
        s = 0;
        jm1 = j - 1;
        if (jm1 < 1) {
            s = a[j][j] - s;
            if (s <= 0) {
                break;
            }
            a[j][j] = Math.sqrt(s);
        } else {
            for (k = 1; k <= jm1; k = k + 1) {
                //~ t = a[k][j] - ddot(k - 1, a[1][k], 1, a[1][j], 1);
                t = a[k][j];
                for (i = 1; i < k; i = i + 1) {
                    t = t - (a[i][j] * a[i][k]);
                }
                t = t / a[k][k];
                a[k][j] = t;
                s = s + t * t;
            }
            s = a[j][j] - s;
            if (s <= 0) {
                break;
            }
            a[j][j] = Math.sqrt(s);
        }
        info[1] = 0;
    }
}

function qpgen2(dmat, dvec, fddmat, n, sol, crval, amat,
    bvec, fdamat, q, meq, iact, nact, iter, work, ierr) {

    var i, j, l, l1, info, it1, iwzv, iwrv, iwrm, iwsv, iwuv, nvl, r, iwnbv,
        temp, sum, t1, tt, gc, gs, nu,
        t1inf, t2min,
        vsmall, tmpa, tmpb,
        go;

    r = Math.min(n, q);
    l = 2 * n + (r * (r + 5)) / 2 + 2 * q + 1;

    vsmall = 1.0e-60;
    do {
        vsmall = vsmall + vsmall;
        tmpa = 1 + 0.1 * vsmall;
        tmpb = 1 + 0.2 * vsmall;
    } while (tmpa <= 1 || tmpb <= 1);

    for (i = 1; i <= n; i = i + 1) {
        work[i] = dvec[i];
    }
    for (i = n + 1; i <= l; i = i + 1) {
        work[i] = 0;
    }
    for (i = 1; i <= q; i = i + 1) {
        iact[i] = 0;
    }

    info = [];

    if (ierr[1] === 0) {
        dpofa(dmat, fddmat, n, info);
        if (info[1] !== 0) {
            ierr[1] = 2;
            return;
        }
        dposl(dmat, fddmat, n, dvec);
        dpori(dmat, fddmat, n);
    } else {
        for (j = 1; j <= n; j = j + 1) {
            sol[j] = 0;
            for (i = 1; i <= j; i = i + 1) {
                sol[j] = sol[j] + dmat[i][j] * dvec[i];
            }
        }
        for (j = 1; j <= n; j = j + 1) {
            dvec[j] = 0;
            for (i = j; i <= n; i = i + 1) {
                dvec[j] = dvec[j] + dmat[j][i] * sol[i];
            }
        }
    }

    crval[1] = 0;
    for (j = 1; j <= n; j = j + 1) {
        sol[j] = dvec[j];
        crval[1] = crval[1] + work[j] * sol[j];
        work[j] = 0;
        for (i = j + 1; i <= n; i = i + 1) {
            dmat[i][j] = 0;
        }
    }
    crval[1] = -crval[1] / 2;
    ierr[1] = 0;

    iwzv = n;
    iwrv = iwzv + n;
    iwuv = iwrv + r;
    iwrm = iwuv + r + 1;
    iwsv = iwrm + (r * (r + 1)) / 2;
    iwnbv = iwsv + q;

    for (i = 1; i <= q; i = i + 1) {
        sum = 0;
        for (j = 1; j <= n; j = j + 1) {
            sum = sum + amat[j][i] * amat[j][i];
        }
        work[iwnbv + i] = Math.sqrt(sum);
    }
    nact = 0;
    iter[1] = 0;
    iter[2] = 0;

    function fn_goto_50() {
        iter[1] = iter[1] + 1;

        l = iwsv;
        for (i = 1; i <= q; i = i + 1) {
            l = l + 1;
            sum = -bvec[i];
            for (j = 1; j <= n; j = j + 1) {
                sum = sum + amat[j][i] * sol[j];
            }
            if (Math.abs(sum) < vsmall) {
                sum = 0;
            }
            if (i > meq) {
                work[l] = sum;
            } else {
                work[l] = -Math.abs(sum);
                if (sum > 0) {
                    for (j = 1; j <= n; j = j + 1) {
                        amat[j][i] = -amat[j][i];
                    }
                    bvec[i] = -bvec[i];
                }
            }
        }

        for (i = 1; i <= nact; i = i + 1) {
            work[iwsv + iact[i]] = 0;
        }

        nvl = 0;
        temp = 0;
        for (i = 1; i <= q; i = i + 1) {
            if (work[iwsv + i] < temp * work[iwnbv + i]) {
                nvl = i;
                temp = work[iwsv + i] / work[iwnbv + i];
            }
        }
        if (nvl === 0) {
            return 999;
        }

        return 0;
    }

    function fn_goto_55() {
        for (i = 1; i <= n; i = i + 1) {
            sum = 0;
            for (j = 1; j <= n; j = j + 1) {
                sum = sum + dmat[j][i] * amat[j][nvl];
            }
            work[i] = sum;
        }

        l1 = iwzv;
        for (i = 1; i <= n; i = i + 1) {
            work[l1 + i] = 0;
        }
        for (j = nact + 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                work[l1 + i] = work[l1 + i] + dmat[i][j] * work[j];
            }
        }

        t1inf = true;
        for (i = nact; i >= 1; i = i - 1) {
            sum = work[i];
            l = iwrm + (i * (i + 3)) / 2;
            l1 = l - i;
            for (j = i + 1; j <= nact; j = j + 1) {
                sum = sum - work[l] * work[iwrv + j];
                l = l + j;
            }
            sum = sum / work[l1];
            work[iwrv + i] = sum;
            if (iact[i] < meq) {
                // continue;
                break;
            }
            if (sum < 0) {
                // continue;
                break;
            }
            t1inf = false;
            it1 = i;
        }

        if (!t1inf) {
            t1 = work[iwuv + it1] / work[iwrv + it1];
            for (i = 1; i <= nact; i = i + 1) {
                if (iact[i] < meq) {
                    // continue;
                    break;
                }
                if (work[iwrv + i] < 0) {
                    // continue;
                    break;
                }
                temp = work[iwuv + i] / work[iwrv + i];
                if (temp < t1) {
                    t1 = temp;
                    it1 = i;
                }
            }
        }

        sum = 0;
        for (i = iwzv + 1; i <= iwzv + n; i = i + 1) {
            sum = sum + work[i] * work[i];
        }
        if (Math.abs(sum) <= vsmall) {
            if (t1inf) {
                ierr[1] = 1;
                // GOTO 999
                return 999;
            } else {
                for (i = 1; i <= nact; i = i + 1) {
                    work[iwuv + i] = work[iwuv + i] - t1 * work[iwrv + i];
                }
                work[iwuv + nact + 1] = work[iwuv + nact + 1] + t1;
                // GOTO 700
                return 700;
            }
        } else {
            sum = 0;
            for (i = 1; i <= n; i = i + 1) {
                sum = sum + work[iwzv + i] * amat[i][nvl];
            }
            tt = -work[iwsv + nvl] / sum;
            t2min = true;
            if (!t1inf) {
                if (t1 < tt) {
                    tt = t1;
                    t2min = false;
                }
            }

            for (i = 1; i <= n; i = i + 1) {
                sol[i] = sol[i] + tt * work[iwzv + i];
                if (Math.abs(sol[i]) < vsmall) {
                    sol[i] = 0;
                }
            }

            crval[1] = crval[1] + tt * sum * (tt / 2 + work[iwuv + nact + 1]);
            for (i = 1; i <= nact; i = i + 1) {
                work[iwuv + i] = work[iwuv + i] - tt * work[iwrv + i];
            }
            work[iwuv + nact + 1] = work[iwuv + nact + 1] + tt;

            if (t2min) {
                nact = nact + 1;
                iact[nact] = nvl;

                l = iwrm + ((nact - 1) * nact) / 2 + 1;
                for (i = 1; i <= nact - 1; i = i + 1) {
                    work[l] = work[i];
                    l = l + 1;
                }

                if (nact === n) {
                    work[l] = work[n];
                } else {
                    for (i = n; i >= nact + 1; i = i - 1) {
                        if (work[i] === 0) {
                            // continue;
                            break;
                        }
                        gc = Math.max(Math.abs(work[i - 1]), Math.abs(work[i]));
                        gs = Math.min(Math.abs(work[i - 1]), Math.abs(work[i]));
                        if (work[i - 1] >= 0) {
                            temp = Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)));
                        } else {
                            temp = -Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)));
                        }
                        gc = work[i - 1] / temp;
                        gs = work[i] / temp;

                        if (gc === 1) {
                            // continue;
                            break;
                        }
                        if (gc === 0) {
                            work[i - 1] = gs * temp;
                            for (j = 1; j <= n; j = j + 1) {
                                temp = dmat[j][i - 1];
                                dmat[j][i - 1] = dmat[j][i];
                                dmat[j][i] = temp;
                            }
                        } else {
                            work[i - 1] = temp;
                            nu = gs / (1 + gc);
                            for (j = 1; j <= n; j = j + 1) {
                                temp = gc * dmat[j][i - 1] + gs * dmat[j][i];
                                dmat[j][i] = nu * (dmat[j][i - 1] + temp) - dmat[j][i];
                                dmat[j][i - 1] = temp;

                            }
                        }
                    }
                    work[l] = work[nact];
                }
            } else {
                sum = -bvec[nvl];
                for (j = 1; j <= n; j = j + 1) {
                    sum = sum + sol[j] * amat[j][nvl];
                }
                if (nvl > meq) {
                    work[iwsv + nvl] = sum;
                } else {
                    work[iwsv + nvl] = -Math.abs(sum);
                    if (sum > 0) {
                        for (j = 1; j <= n; j = j + 1) {
                            amat[j][nvl] = -amat[j][nvl];
                        }
                        bvec[nvl] = -bvec[nvl];
                    }
                }
                // GOTO 700
                return 700;
            }
        }

        return 0;
    }

    function fn_goto_797() {
        l = iwrm + (it1 * (it1 + 1)) / 2 + 1;
        l1 = l + it1;
        if (work[l1] === 0) {
            // GOTO 798
            return 798;
        }
        gc = Math.max(Math.abs(work[l1 - 1]), Math.abs(work[l1]));
        gs = Math.min(Math.abs(work[l1 - 1]), Math.abs(work[l1]));
        if (work[l1 - 1] >= 0) {
            temp = Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)));
        } else {
            temp = -Math.abs(gc * Math.sqrt(1 + gs * gs / (gc * gc)));
        }
        gc = work[l1 - 1] / temp;
        gs = work[l1] / temp;

        if (gc === 1) {
            // GOTO 798
            return 798;
        }
        if (gc === 0) {
            for (i = it1 + 1; i <= nact; i = i + 1) {
                temp = work[l1 - 1];
                work[l1 - 1] = work[l1];
                work[l1] = temp;
                l1 = l1 + i;
            }
            for (i = 1; i <= n; i = i + 1) {
                temp = dmat[i][it1];
                dmat[i][it1] = dmat[i][it1 + 1];
                dmat[i][it1 + 1] = temp;
            }
        } else {
            nu = gs / (1 + gc);
            for (i = it1 + 1; i <= nact; i = i + 1) {
                temp = gc * work[l1 - 1] + gs * work[l1];
                work[l1] = nu * (work[l1 - 1] + temp) - work[l1];
                work[l1 - 1] = temp;
                l1 = l1 + i;
            }
            for (i = 1; i <= n; i = i + 1) {
                temp = gc * dmat[i][it1] + gs * dmat[i][it1 + 1];
                dmat[i][it1 + 1] = nu * (dmat[i][it1] + temp) - dmat[i][it1 + 1];
                dmat[i][it1] = temp;
            }
        }

        return 0;
    }

    function fn_goto_798() {
        l1 = l - it1;
        for (i = 1; i <= it1; i = i + 1) {
            work[l1] = work[l];
            l = l + 1;
            l1 = l1 + 1;
        }

        work[iwuv + it1] = work[iwuv + it1 + 1];
        iact[it1] = iact[it1 + 1];
        it1 = it1 + 1;
        if (it1 < nact) {
            // GOTO 797
            return 797;
        }

        return 0;
    }

    function fn_goto_799() {
        work[iwuv + nact] = work[iwuv + nact + 1];
        work[iwuv + nact + 1] = 0;
        iact[nact] = 0;
        nact = nact - 1;
        iter[2] = iter[2] + 1;

        return 0;
    }

    go = 0;
    while (true) {
        go = fn_goto_50();
        if (go === 999) {
            return;
        }
        while (true) {
            go = fn_goto_55();
            if (go === 0) {
                break;
            }
            if (go === 999) {
                return;
            }
            if (go === 700) {
                if (it1 === nact) {
                    fn_goto_799();
                } else {
                    while (true) {
                        fn_goto_797();
                        go = fn_goto_798();
                        if (go !== 797) {
                            break;
                        }
                    }
                    fn_goto_799();
                }
            }
        }
    }

}

function solveQP(Dmat, dvec, Amat, bvec, meq, factorized) {
    Dmat = base0to1(Dmat);
    dvec = base0to1(dvec);
    Amat = base0to1(Amat);
    var i, n, q,
        nact, r,
        crval = [], iact = [], sol = [], work = [], iter = [],
        message;

    meq = meq || 0;
    factorized = factorized ? base0to1(factorized) : [undefined, 0];
    bvec = bvec ? base0to1(bvec) : [];

    // In Fortran the array index starts from 1
    n = Dmat.length - 1;
    q = Amat[1].length - 1;

    if (!bvec) {
        for (i = 1; i <= q; i = i + 1) {
            bvec[i] = 0;
        }
    }
    for (i = 1; i <= q; i = i + 1) {
        iact[i] = 0;
    }
    nact = 0;
    r = Math.min(n, q);
    for (i = 1; i <= n; i = i + 1) {
        sol[i] = 0;
    }
    crval[1] = 0;
    for (i = 1; i <= (2 * n + (r * (r + 5)) / 2 + 2 * q + 1); i = i + 1) {
        work[i] = 0;
    }
    for (i = 1; i <= 2; i = i + 1) {
        iter[i] = 0;
    }

    qpgen2(Dmat, dvec, n, n, sol, crval, Amat,
        bvec, n, q, meq, iact, nact, iter, work, factorized);

    message = "";
    if (factorized[1] === 1) {
        message = "constraints are inconsistent, no solution!";
    }
    if (factorized[1] === 2) {
        message = "matrix D in quadratic function is not positive definite!";
    }

    return {
        solution: base1to0(sol),
        value: base1to0(crval),
        unconstrained_solution: base1to0(dvec),
        iterations: base1to0(iter),
        iact: base1to0(iact),
        message: message
    };
}
exports.solveQP = solveQP;
}(numeric));
/*
Shanti Rao sent me this routine by private email. I had to modify it
slightly to work on Arrays instead of using a Matrix object.
It is apparently translated from http://stitchpanorama.sourceforge.net/Python/svd.py
*/

numeric.svd= function svd(A) {
    var temp;
//Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
	var prec= numeric.epsilon; //Math.pow(2,-52) // assumes double prec
	var tolerance= 1.e-64/prec;
	var itmax= 50;
	var c=0;
	var i=0;
	var j=0;
	var k=0;
	var l=0;
	
	var u= numeric.clone(A);
	var m= u.length;
	
	var n= u[0].length;
	
	if (m < n) throw "Need more rows than columns"
	
	var e = new Array(n);
	var q = new Array(n);
	for (i=0; i<n; i++) e[i] = q[i] = 0.0;
	var v = numeric.rep([n,n],0);
//	v.zero();
	
 	function pythag(a,b)
 	{
		a = Math.abs(a)
		b = Math.abs(b)
		if (a > b)
			return a*Math.sqrt(1.0+(b*b/a/a))
		else if (b == 0.0) 
			return a
		return b*Math.sqrt(1.0+(a*a/b/b))
	}

	//Householder's reduction to bidiagonal form

	var f= 0.0;
	var g= 0.0;
	var h= 0.0;
	var x= 0.0;
	var y= 0.0;
	var z= 0.0;
	var s= 0.0;
	
	for (i=0; i < n; i++)
	{	
		e[i]= g;
		s= 0.0;
		l= i+1;
		for (j=i; j < m; j++) 
			s += (u[j][i]*u[j][i]);
		if (s <= tolerance)
			g= 0.0;
		else
		{	
			f= u[i][i];
			g= Math.sqrt(s);
			if (f >= 0.0) g= -g;
			h= f*g-s
			u[i][i]=f-g;
			for (j=l; j < n; j++)
			{
				s= 0.0
				for (k=i; k < m; k++) 
					s += u[k][i]*u[k][j]
				f= s/h
				for (k=i; k < m; k++) 
					u[k][j]+=f*u[k][i]
			}
		}
		q[i]= g
		s= 0.0
		for (j=l; j < n; j++) 
			s= s + u[i][j]*u[i][j]
		if (s <= tolerance)
			g= 0.0
		else
		{	
			f= u[i][i+1]
			g= Math.sqrt(s)
			if (f >= 0.0) g= -g
			h= f*g - s
			u[i][i+1] = f-g;
			for (j=l; j < n; j++) e[j]= u[i][j]/h
			for (j=l; j < m; j++)
			{	
				s=0.0
				for (k=l; k < n; k++) 
					s += (u[j][k]*u[i][k])
				for (k=l; k < n; k++) 
					u[j][k]+=s*e[k]
			}	
		}
		y= Math.abs(q[i])+Math.abs(e[i])
		if (y>x) 
			x=y
	}
	
	// accumulation of right hand gtransformations
	for (i=n-1; i != -1; i+= -1)
	{	
		if (g != 0.0)
		{
		 	h= g*u[i][i+1]
			for (j=l; j < n; j++) 
				v[j][i]=u[i][j]/h
			for (j=l; j < n; j++)
			{	
				s=0.0
				for (k=l; k < n; k++) 
					s += u[i][k]*v[k][j]
				for (k=l; k < n; k++) 
					v[k][j]+=(s*v[k][i])
			}	
		}
		for (j=l; j < n; j++)
		{
			v[i][j] = 0;
			v[j][i] = 0;
		}
		v[i][i] = 1;
		g= e[i]
		l= i
	}
	
	// accumulation of left hand transformations
	for (i=n-1; i != -1; i+= -1)
	{	
		l= i+1
		g= q[i]
		for (j=l; j < n; j++) 
			u[i][j] = 0;
		if (g != 0.0)
		{
			h= u[i][i]*g
			for (j=l; j < n; j++)
			{
				s=0.0
				for (k=l; k < m; k++) s += u[k][i]*u[k][j];
				f= s/h
				for (k=i; k < m; k++) u[k][j]+=f*u[k][i];
			}
			for (j=i; j < m; j++) u[j][i] = u[j][i]/g;
		}
		else
			for (j=i; j < m; j++) u[j][i] = 0;
		u[i][i] += 1;
	}
	
	// diagonalization of the bidiagonal form
	prec= prec*x
	for (k=n-1; k != -1; k+= -1)
	{
		for (var iteration=0; iteration < itmax; iteration++)
		{	// test f splitting
			var test_convergence = false
			for (l=k; l != -1; l+= -1)
			{	
				if (Math.abs(e[l]) <= prec)
				{	test_convergence= true
					break 
				}
				if (Math.abs(q[l-1]) <= prec)
					break 
			}
			if (!test_convergence)
			{	// cancellation of e[l] if l>0
				c= 0.0
				s= 1.0
				var l1= l-1
				for (i =l; i<k+1; i++)
				{	
					f= s*e[i]
					e[i]= c*e[i]
					if (Math.abs(f) <= prec)
						break
					g= q[i]
					h= pythag(f,g)
					q[i]= h
					c= g/h
					s= -f/h
					for (j=0; j < m; j++)
					{	
						y= u[j][l1]
						z= u[j][i]
						u[j][l1] =  y*c+(z*s)
						u[j][i] = -y*s+(z*c)
					} 
				}	
			}
			// test f convergence
			z= q[k]
			if (l== k)
			{	//convergence
				if (z<0.0)
				{	//q[k] is made non-negative
					q[k]= -z
					for (j=0; j < n; j++)
						v[j][k] = -v[j][k]
				}
				break  //break out of iteration loop and move on to next k value
			}
			if (iteration >= itmax-1)
				throw 'Error: no convergence.'
			// shift from bottom 2x2 minor
			x= q[l]
			y= q[k-1]
			g= e[k-1]
			h= e[k]
			f= ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
			g= pythag(f,1.0)
			if (f < 0.0)
				f= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
			else
				f= ((x-z)*(x+z)+h*(y/(f+g)-h))/x
			// next QR transformation
			c= 1.0
			s= 1.0
			for (i=l+1; i< k+1; i++)
			{	
				g= e[i]
				y= q[i]
				h= s*g
				g= c*g
				z= pythag(f,h)
				e[i-1]= z
				c= f/z
				s= h/z
				f= x*c+g*s
				g= -x*s+g*c
				h= y*s
				y= y*c
				for (j=0; j < n; j++)
				{	
					x= v[j][i-1]
					z= v[j][i]
					v[j][i-1] = x*c+z*s
					v[j][i] = -x*s+z*c
				}
				z= pythag(f,h)
				q[i-1]= z
				c= f/z
				s= h/z
				f= c*g+s*y
				x= -s*g+c*y
				for (j=0; j < m; j++)
				{
					y= u[j][i-1]
					z= u[j][i]
					u[j][i-1] = y*c+z*s
					u[j][i] = -y*s+z*c
				}
			}
			e[l]= 0.0
			e[k]= f
			q[k]= x
		} 
	}
		
	//vt= transpose(v)
	//return (u,q,vt)
	for (i=0;i<q.length; i++) 
	  if (q[i] < prec) q[i] = 0
	  
	//sort eigenvalues	
	for (i=0; i< n; i++)
	{	 
	//writeln(q)
	 for (j=i-1; j >= 0; j--)
	 {
	  if (q[j] < q[i])
	  {
	//  writeln(i,'-',j)
	   c = q[j]
	   q[j] = q[i]
	   q[i] = c
	   for(k=0;k<u.length;k++) { temp = u[k][i]; u[k][i] = u[k][j]; u[k][j] = temp; }
	   for(k=0;k<v.length;k++) { temp = v[k][i]; v[k][i] = v[k][j]; v[k][j] = temp; }
//	   u.swapCols(i,j)
//	   v.swapCols(i,j)
	   i = j	   
	  }
	 }	
	}
	
	return {U:u,S:q,V:v}
};

/**
 * Viola-Jones-like face detection algorithm
 * Some explanation here: http://liuliu.me/eyes/javascript-face-detection-explained/
 *
 * @author Liu Liu / github.com/liuliu
 *
 * Copyright (c) 2010, Liu Liu
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * * Neither the name of the authors nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

var ccv = {};

ccv.grayscale = function (canvas) {
  /* detect_objects requires gray-scale image */
  var ctx = canvas.getContext("2d");
  var imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
  var data = imageData.data;
  var pix1, pix2, pix = canvas.width * canvas.height * 4;
  while (pix > 0)
    data[pix -= 4] = data[pix1 = pix + 1] = data[pix2 = pix + 2] = (data[pix] * 0.3 + data[pix1] * 0.59 + data[pix2] * 0.11);
  ctx.putImageData(imageData, 0, 0);
  return canvas;
};

ccv.array_group = function (seq, gfunc) {
  var i, j;
  var node = new Array(seq.length);
  for (i = 0; i < seq.length; i++)
    node[i] = {"parent" : -1,
           "element" : seq[i],
           "rank" : 0};
  for (i = 0; i < seq.length; i++)
  {
    if (!node[i].element)
      continue;
    var root = i;
    while (node[root].parent != -1)
      root = node[root].parent;
    for (j = 0; j < seq.length; j++)
    {
      if( i != j && node[j].element && gfunc(node[i].element, node[j].element))
      {
        var root2 = j;

        while (node[root2].parent != -1)
          root2 = node[root2].parent;

        if(root2 != root)
        {
          if(node[root].rank > node[root2].rank)
            node[root2].parent = root;
          else
          {
            node[root].parent = root2;
            if (node[root].rank == node[root2].rank)
            node[root2].rank++;
            root = root2;
          }

          /* compress path from node2 to the root: */
          var temp, node2 = j;
          while (node[node2].parent != -1)
          {
            temp = node2;
            node2 = node[node2].parent;
            node[temp].parent = root;
          }

          /* compress path from node to the root: */
          node2 = i;
          while (node[node2].parent != -1)
          {
            temp = node2;
            node2 = node[node2].parent;
            node[temp].parent = root;
          }
        }
      }
    }
  }
  var idx = new Array(seq.length);
  var class_idx = 0;
  for(i = 0; i < seq.length; i++)
  {
    j = -1;
    var node1 = i;
    if(node[node1].element)
    {
      while (node[node1].parent != -1)
        node1 = node[node1].parent;
      if(node[node1].rank >= 0)
        node[node1].rank = ~class_idx++;
      j = ~node[node1].rank;
    }
    idx[i] = j;
  }
  return {"index" : idx, "cat" : class_idx};
};

ccv.detect_objects = function (canvas, cascade, interval, min_neighbors) {
  var scale = Math.pow(2, 1 / (interval + 1));
  var next = interval + 1;
  var scale_upto = Math.floor(Math.log(Math.min(cascade.width, cascade.height)) / Math.log(scale));
  var pyr = new Array((scale_upto + next * 2) * 4);
  pyr[0] = canvas;
  pyr[0].data = pyr[0].getContext("2d").getImageData(0, 0, pyr[0].width, pyr[0].height).data;
  var i, j, k, x, y, q;
  for (i = 1; i <= interval; i++) {
    pyr[i * 4] = document.createElement("canvas");
    pyr[i * 4].width = Math.floor(pyr[0].width / Math.pow(scale, i));
    pyr[i * 4].height = Math.floor(pyr[0].height / Math.pow(scale, i));
    pyr[i * 4].getContext("2d").drawImage(pyr[0], 0, 0, pyr[0].width, pyr[0].height, 0, 0, pyr[i * 4].width, pyr[i * 4].height);
    pyr[i * 4].data = pyr[i * 4].getContext("2d").getImageData(0, 0, pyr[i * 4].width, pyr[i * 4].height).data;
  }
  for (i = next; i < scale_upto + next * 2; i++) {
    pyr[i * 4] = document.createElement("canvas");
    pyr[i * 4].width = Math.floor(pyr[i * 4 - next * 4].width / 2);
    pyr[i * 4].height = Math.floor(pyr[i * 4 - next * 4].height / 2);
    pyr[i * 4].getContext("2d").drawImage(pyr[i * 4 - next * 4], 0, 0, pyr[i * 4 - next * 4].width, pyr[i * 4 - next * 4].height, 0, 0, pyr[i * 4].width, pyr[i * 4].height);
    pyr[i * 4].data = pyr[i * 4].getContext("2d").getImageData(0, 0, pyr[i * 4].width, pyr[i * 4].height).data;
  }
  for (i = next * 2; i < scale_upto + next * 2; i++) {
    pyr[i * 4 + 1] = document.createElement("canvas");
    pyr[i * 4 + 1].width = Math.floor(pyr[i * 4 - next * 4].width / 2);
    pyr[i * 4 + 1].height = Math.floor(pyr[i * 4 - next * 4].height / 2);
    pyr[i * 4 + 1].getContext("2d").drawImage(pyr[i * 4 - next * 4], 1, 0, pyr[i * 4 - next * 4].width - 1, pyr[i * 4 - next * 4].height, 0, 0, pyr[i * 4 + 1].width - 2, pyr[i * 4 + 1].height);
    pyr[i * 4 + 1].data = pyr[i * 4 + 1].getContext("2d").getImageData(0, 0, pyr[i * 4 + 1].width, pyr[i * 4 + 1].height).data;
    pyr[i * 4 + 2] = document.createElement("canvas");
    pyr[i * 4 + 2].width = Math.floor(pyr[i * 4 - next * 4].width / 2);
    pyr[i * 4 + 2].height = Math.floor(pyr[i * 4 - next * 4].height / 2);
    pyr[i * 4 + 2].getContext("2d").drawImage(pyr[i * 4 - next * 4], 0, 1, pyr[i * 4 - next * 4].width, pyr[i * 4 - next * 4].height - 1, 0, 0, pyr[i * 4 + 2].width, pyr[i * 4 + 2].height - 2);
    pyr[i * 4 + 2].data = pyr[i * 4 + 2].getContext("2d").getImageData(0, 0, pyr[i * 4 + 2].width, pyr[i * 4 + 2].height).data;
    pyr[i * 4 + 3] = document.createElement("canvas");
    pyr[i * 4 + 3].width = Math.floor(pyr[i * 4 - next * 4].width / 2);
    pyr[i * 4 + 3].height = Math.floor(pyr[i * 4 - next * 4].height / 2);
    pyr[i * 4 + 3].getContext("2d").drawImage(pyr[i * 4 - next * 4], 1, 1, pyr[i * 4 - next * 4].width - 1, pyr[i * 4 - next * 4].height - 1, 0, 0, pyr[i * 4 + 3].width - 2, pyr[i * 4 + 3].height - 2);
    pyr[i * 4 + 3].data = pyr[i * 4 + 3].getContext("2d").getImageData(0, 0, pyr[i * 4 + 3].width, pyr[i * 4 + 3].height).data;
  }
  for (j = 0; j < cascade.stage_classifier.length; j++)
    cascade.stage_classifier[j].orig_feature = cascade.stage_classifier[j].feature;
  var scale_x = 1, scale_y = 1;
  var dx = [0, 1, 0, 1];
  var dy = [0, 0, 1, 1];
  var seq = [];
  for (i = 0; i < scale_upto; i++) {
    var qw = pyr[i * 4 + next * 8].width - Math.floor(cascade.width / 4);
    var qh = pyr[i * 4 + next * 8].height - Math.floor(cascade.height / 4);
    var step = [pyr[i * 4].width * 4, pyr[i * 4 + next * 4].width * 4, pyr[i * 4 + next * 8].width * 4];
    var paddings = [pyr[i * 4].width * 16 - qw * 16,
            pyr[i * 4 + next * 4].width * 8 - qw * 8,
            pyr[i * 4 + next * 8].width * 4 - qw * 4];
    for (j = 0; j < cascade.stage_classifier.length; j++) {
      var orig_feature = cascade.stage_classifier[j].orig_feature;
      var feature = cascade.stage_classifier[j].feature = new Array(cascade.stage_classifier[j].count);
      for (k = 0; k < cascade.stage_classifier[j].count; k++) {
        feature[k] = {"size" : orig_feature[k].size,
                "px" : new Array(orig_feature[k].size),
                "pz" : new Array(orig_feature[k].size),
                "nx" : new Array(orig_feature[k].size),
                "nz" : new Array(orig_feature[k].size)};
        for (q = 0; q < orig_feature[k].size; q++) {
          feature[k].px[q] = orig_feature[k].px[q] * 4 + orig_feature[k].py[q] * step[orig_feature[k].pz[q]];
          feature[k].pz[q] = orig_feature[k].pz[q];
          feature[k].nx[q] = orig_feature[k].nx[q] * 4 + orig_feature[k].ny[q] * step[orig_feature[k].nz[q]];
          feature[k].nz[q] = orig_feature[k].nz[q];
        }
      }
    }
    for (q = 0; q < 4; q++) {
      var u8 = [pyr[i * 4].data, pyr[i * 4 + next * 4].data, pyr[i * 4 + next * 8 + q].data];
      var u8o = [dx[q] * 8 + dy[q] * pyr[i * 4].width * 8, dx[q] * 4 + dy[q] * pyr[i * 4 + next * 4].width * 4, 0];
      for (y = 0; y < qh; y++) {
        for (x = 0; x < qw; x++) {
          var sum = 0;
          var flag = true;
          for (j = 0; j < cascade.stage_classifier.length; j++) {
            sum = 0;
            var alpha = cascade.stage_classifier[j].alpha;
            var feature = cascade.stage_classifier[j].feature;
            for (k = 0; k < cascade.stage_classifier[j].count; k++) {
              var feature_k = feature[k];
              var p, pmin = u8[feature_k.pz[0]][u8o[feature_k.pz[0]] + feature_k.px[0]];
              var n, nmax = u8[feature_k.nz[0]][u8o[feature_k.nz[0]] + feature_k.nx[0]];
              if (pmin <= nmax) {
                sum += alpha[k * 2];
              } else {
                var f, shortcut = true;
                for (f = 0; f < feature_k.size; f++) {
                  if (feature_k.pz[f] >= 0) {
                    p = u8[feature_k.pz[f]][u8o[feature_k.pz[f]] + feature_k.px[f]];
                    if (p < pmin) {
                      if (p <= nmax) {
                        shortcut = false;
                        break;
                      }
                      pmin = p;
                    }
                  }
                  if (feature_k.nz[f] >= 0) {
                    n = u8[feature_k.nz[f]][u8o[feature_k.nz[f]] + feature_k.nx[f]];
                    if (n > nmax) {
                      if (pmin <= n) {
                        shortcut = false;
                        break;
                      }
                      nmax = n;
                    }
                  }
                }
                sum += (shortcut) ? alpha[k * 2 + 1] : alpha[k * 2];
              }
            }
            if (sum < cascade.stage_classifier[j].threshold) {
              flag = false;
              break;
            }
          }
          if (flag) {
            seq.push({"x" : (x * 4 + dx[q] * 2) * scale_x,
                  "y" : (y * 4 + dy[q] * 2) * scale_y,
                  "width" : cascade.width * scale_x,
                  "height" : cascade.height * scale_y,
                  "neighbor" : 1,
                  "confidence" : sum});
          }
          u8o[0] += 16;
          u8o[1] += 8;
          u8o[2] += 4;
        }
        u8o[0] += paddings[0];
        u8o[1] += paddings[1];
        u8o[2] += paddings[2];
      }
    }
    scale_x *= scale;
    scale_y *= scale;
  }
  for (j = 0; j < cascade.stage_classifier.length; j++)
    cascade.stage_classifier[j].feature = cascade.stage_classifier[j].orig_feature;
  if (!(min_neighbors > 0))
    return seq;
  else {
    var result = ccv.array_group(seq, function (r1, r2) {
      var distance = Math.floor(r1.width * 0.25 + 0.5);

      return r2.x <= r1.x + distance &&
           r2.x >= r1.x - distance &&
           r2.y <= r1.y + distance &&
           r2.y >= r1.y - distance &&
           r2.width <= Math.floor(r1.width * 1.5 + 0.5) &&
           Math.floor(r2.width * 1.5 + 0.5) >= r1.width;
    });
    var ncomp = result.cat;
    var idx_seq = result.index;
    var comps = new Array(ncomp + 1);
    for (i = 0; i < comps.length; i++)
      comps[i] = {"neighbors" : 0,
            "x" : 0,
            "y" : 0,
            "width" : 0,
            "height" : 0,
            "confidence" : 0};

    // count number of neighbors
    for(i = 0; i < seq.length; i++)
    {
      var r1 = seq[i];
      var idx = idx_seq[i];

      if (comps[idx].neighbors == 0)
        comps[idx].confidence = r1.confidence;

      ++comps[idx].neighbors;

      comps[idx].x += r1.x;
      comps[idx].y += r1.y;
      comps[idx].width += r1.width;
      comps[idx].height += r1.height;
      comps[idx].confidence = Math.max(comps[idx].confidence, r1.confidence);
    }

    var seq2 = [];
    // calculate average bounding box
    for(i = 0; i < ncomp; i++)
    {
      var n = comps[i].neighbors;
      if (n >= min_neighbors)
        seq2.push({"x" : (comps[i].x * 2 + n) / (2 * n),
               "y" : (comps[i].y * 2 + n) / (2 * n),
               "width" : (comps[i].width * 2 + n) / (2 * n),
               "height" : (comps[i].height * 2 + n) / (2 * n),
               "neighbors" : comps[i].neighbors,
               "confidence" : comps[i].confidence});
    }

    var result_seq = [];
    // filter out small face rectangles inside large face rectangles
    for(i = 0; i < seq2.length; i++)
    {
      var r1 = seq2[i];
      var flag = true;
      for(j = 0; j < seq2.length; j++)
      {
        var r2 = seq2[j];
        var distance = Math.floor(r2.width * 0.25 + 0.5);

        if(i != j &&
           r1.x >= r2.x - distance &&
           r1.y >= r2.y - distance &&
           r1.x + r1.width <= r2.x + r2.width + distance &&
           r1.y + r1.height <= r2.y + r2.height + distance &&
           (r2.neighbors > Math.max(3, r1.neighbors) || r1.neighbors < 3))
        {
          flag = false;
          break;
        }
      }

      if(flag)
        result_seq.push(r1);
    }
    return result_seq;
  }
};

/**
 * Data for ccv facedetection
 *
 * @author Liu Liu / github.com/liuliu
 *
 * Copyright (c) 2010, Liu Liu
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * * Neither the name of the authors nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

ccv.cascade = {"count" : 16, "width" : 24, "height" : 24, "stage_classifier" : [{"count":4,"threshold":-4.577530e+00,"feature":[{"size":4,"px":[3,5,8,11],"py":[2,2,6,3],"pz":[2,1,1,0],"nx":[8,4,0,0],"ny":[4,4,0,0],"nz":[1,1,-1,-1]},{"size":3,"px":[3,6,7],"py":[7,13,0],"pz":[1,0,-1],"nx":[2,3,4],"ny":[5,4,4],"nz":[2,1,1]},{"size":5,"px":[5,3,10,13,11],"py":[1,0,3,2,2],"pz":[1,2,0,0,0],"nx":[0,11,0,11,11],"ny":[0,2,3,1,1],"nz":[1,1,0,1,-1]},{"size":5,"px":[6,12,12,9,12],"py":[4,13,12,7,11],"pz":[1,0,0,1,0],"nx":[8,0,8,2,11],"ny":[4,0,8,5,1],"nz":[1,-1,-1,-1,-1]}],"alpha":[-2.879683e+00,2.879683e+00,-1.569341e+00,1.569341e+00,-1.286131e+00,1.286131e+00,-1.157626e+00,1.157626e+00]},{"count":4,"threshold":-4.339908e+00,"feature":[{"size":5,"px":[13,12,3,11,17],"py":[3,3,1,4,13],"pz":[0,0,2,0,0],"nx":[4,3,8,15,15],"ny":[4,5,4,8,8],"nz":[1,2,1,0,-1]},{"size":5,"px":[6,7,6,3,3],"py":[13,13,4,2,7],"pz":[0,0,1,2,1],"nx":[4,8,3,0,15],"ny":[4,4,4,3,8],"nz":[1,1,-1,-1,-1]},{"size":3,"px":[2,2,11],"py":[3,2,5],"pz":[2,2,0],"nx":[3,8,3],"ny":[4,4,4],"nz":[1,-1,-1]},{"size":5,"px":[15,13,9,11,7],"py":[2,1,2,1,0],"pz":[0,0,0,0,1],"nx":[23,11,23,22,23],"ny":[1,0,2,0,0],"nz":[0,1,0,0,0]}],"alpha":[-2.466029e+00,2.466029e+00,-1.839510e+00,1.839510e+00,-1.060559e+00,1.060559e+00,-1.094927e+00,1.094927e+00]},{"count":7,"threshold":-5.052474e+00,"feature":[{"size":5,"px":[17,13,3,11,10],"py":[13,2,1,4,3],"pz":[0,0,2,0,0],"nx":[4,8,8,3,7],"ny":[2,8,4,5,4],"nz":[2,0,1,2,1]},{"size":5,"px":[6,7,3,6,6],"py":[4,12,2,13,14],"pz":[1,0,2,0,0],"nx":[8,3,4,4,3],"ny":[4,4,2,0,2],"nz":[1,1,-1,-1,-1]},{"size":5,"px":[7,4,5,3,3],"py":[2,1,3,1,1],"pz":[0,1,0,1,-1],"nx":[1,0,1,1,0],"ny":[1,3,2,0,4],"nz":[0,0,0,0,0]},{"size":5,"px":[11,11,11,3,2],"py":[11,13,10,7,2],"pz":[0,0,0,1,2],"nx":[4,1,8,2,0],"ny":[4,1,12,0,4],"nz":[1,-1,-1,-1,-1]},{"size":3,"px":[9,13,1],"py":[7,19,4],"pz":[1,-1,-1],"nx":[4,7,4],"ny":[5,8,2],"nz":[2,1,2]},{"size":5,"px":[12,8,16,4,4],"py":[12,1,2,0,0],"pz":[0,1,0,2,-1],"nx":[11,22,11,23,23],"ny":[2,0,1,1,5],"nz":[1,0,1,0,0]},{"size":3,"px":[11,17,17],"py":[6,11,12],"pz":[0,0,0],"nx":[15,1,11],"ny":[9,1,1],"nz":[0,-1,-1]}],"alpha":[-2.156890e+00,2.156890e+00,-1.718246e+00,1.718246e+00,-9.651329e-01,9.651329e-01,-9.948090e-01,9.948090e-01,-8.802466e-01,8.802466e-01,-8.486741e-01,8.486741e-01,-8.141777e-01,8.141777e-01]},{"count":13,"threshold":-5.774400e+00,"feature":[{"size":5,"px":[6,10,3,12,14],"py":[5,3,1,2,2],"pz":[1,0,2,0,0],"nx":[3,4,14,8,4],"ny":[5,4,8,4,2],"nz":[2,1,0,1,2]},{"size":5,"px":[10,6,11,5,12],"py":[4,13,4,2,4],"pz":[0,0,0,1,0],"nx":[1,4,8,1,1],"ny":[2,4,4,4,3],"nz":[0,1,1,0,0]},{"size":3,"px":[18,6,12],"py":[12,4,8],"pz":[0,1,0],"nx":[7,4,8],"ny":[4,2,4],"nz":[1,-1,-1]},{"size":5,"px":[7,5,6,3,17],"py":[13,12,3,8,13],"pz":[0,0,1,1,0],"nx":[3,3,0,1,8],"ny":[4,5,5,10,4],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[16,7,16,7,7],"py":[1,1,2,0,0],"pz":[0,1,0,1,-1],"nx":[23,23,23,11,5],"ny":[2,14,1,2,1],"nz":[0,0,0,1,2]},{"size":3,"px":[9,18,16],"py":[7,14,2],"pz":[1,0,-1],"nx":[8,4,9],"ny":[10,2,4],"nz":[1,2,1]},{"size":4,"px":[3,16,1,22],"py":[7,4,5,11],"pz":[1,-1,-1,-1],"nx":[3,9,4,2],"ny":[4,9,7,5],"nz":[1,0,1,2]},{"size":5,"px":[4,7,8,8,9],"py":[0,2,2,1,1],"pz":[1,0,0,0,0],"nx":[0,0,1,0,0],"ny":[15,16,19,0,14],"nz":[0,0,0,1,0]},{"size":5,"px":[4,4,7,8,12],"py":[2,5,6,7,10],"pz":[2,2,1,1,0],"nx":[8,5,10,0,0],"ny":[4,2,5,3,14],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[11,0],"py":[13,4],"pz":[0,-1],"nx":[3,14],"ny":[4,16],"nz":[1,0]},{"size":5,"px":[17,8,18,4,4],"py":[3,1,3,0,0],"pz":[0,1,0,2,-1],"nx":[21,22,5,11,22],"ny":[0,1,0,1,2],"nz":[0,0,2,1,0]},{"size":4,"px":[7,8,2,11],"py":[13,12,2,7],"pz":[0,0,2,0],"nx":[4,0,23,3],"ny":[4,1,1,11],"nz":[1,-1,-1,-1]},{"size":5,"px":[4,18,8,9,15],"py":[4,16,7,7,23],"pz":[2,0,1,1,0],"nx":[0,1,1,1,1],"ny":[10,21,23,22,22],"nz":[1,0,0,0,-1]}],"alpha":[-1.956565e+00,1.956565e+00,-1.262438e+00,1.262438e+00,-1.056941e+00,1.056941e+00,-9.712509e-01,9.712509e-01,-8.261028e-01,8.261028e-01,-8.456506e-01,8.456506e-01,-6.652113e-01,6.652113e-01,-6.026287e-01,6.026287e-01,-6.915425e-01,6.915425e-01,-5.539286e-01,5.539286e-01,-5.515072e-01,5.515072e-01,-6.685884e-01,6.685884e-01,-4.656070e-01,4.656070e-01]},{"count":20,"threshold":-5.606853e+00,"feature":[{"size":5,"px":[17,11,6,14,9],"py":[13,4,4,3,3],"pz":[0,0,1,0,0],"nx":[14,4,8,7,8],"ny":[8,4,4,4,8],"nz":[0,1,1,1,0]},{"size":5,"px":[3,9,10,11,11],"py":[7,2,2,3,3],"pz":[1,0,0,0,-1],"nx":[3,8,4,2,5],"ny":[4,4,10,2,8],"nz":[1,1,1,2,1]},{"size":5,"px":[12,12,12,5,12],"py":[12,9,10,12,11],"pz":[0,0,0,0,0],"nx":[0,0,0,0,0],"ny":[2,1,3,0,0],"nz":[0,0,0,0,-1]},{"size":5,"px":[9,18,9,9,12],"py":[7,14,19,5,11],"pz":[1,-1,-1,-1,-1],"nx":[23,4,23,23,8],"ny":[13,5,14,16,4],"nz":[0,2,0,0,1]},{"size":5,"px":[12,12,12,6,1],"py":[13,11,12,6,5],"pz":[0,0,0,-1,-1],"nx":[4,6,8,4,9],"ny":[2,8,4,4,4],"nz":[2,1,1,1,1]},{"size":4,"px":[12,11,11,6],"py":[5,5,6,13],"pz":[0,0,0,0],"nx":[8,3,2,8],"ny":[4,4,17,2],"nz":[1,1,-1,-1]},{"size":5,"px":[3,14,12,15,13],"py":[0,2,2,2,2],"pz":[2,0,0,0,0],"nx":[22,23,22,23,7],"ny":[0,3,1,2,4],"nz":[0,0,0,0,1]},{"size":5,"px":[16,15,18,19,9],"py":[12,11,12,12,9],"pz":[0,0,0,0,1],"nx":[8,2,22,23,21],"ny":[4,1,1,2,20],"nz":[1,-1,-1,-1,-1]},{"size":3,"px":[4,7,7],"py":[0,2,2],"pz":[1,0,-1],"nx":[1,2,2],"ny":[2,0,2],"nz":[1,0,0]},{"size":3,"px":[4,11,11],"py":[6,9,8],"pz":[1,0,0],"nx":[9,2,8],"ny":[9,4,5],"nz":[0,-1,-1]},{"size":4,"px":[2,7,6,6],"py":[4,23,21,22],"pz":[2,0,0,0],"nx":[9,3,8,17],"ny":[21,2,5,1],"nz":[0,-1,-1,-1]},{"size":2,"px":[2,8],"py":[4,12],"pz":[2,0],"nx":[3,0],"ny":[4,4],"nz":[1,-1]},{"size":5,"px":[4,5,1,8,4],"py":[15,12,3,23,12],"pz":[0,0,2,0,0],"nx":[0,0,0,0,0],"ny":[23,10,22,21,11],"nz":[0,1,0,0,-1]},{"size":2,"px":[21,5],"py":[13,4],"pz":[0,2],"nx":[23,4],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[15,17],"py":[2,3],"pz":[0,0],"nx":[19,20],"ny":[2,1],"nz":[0,0]},{"size":5,"px":[12,1,8,17,4],"py":[14,2,13,6,12],"pz":[0,-1,-1,-1,-1],"nx":[8,13,15,15,7],"ny":[10,9,15,14,8],"nz":[1,0,0,0,1]},{"size":2,"px":[8,5],"py":[7,4],"pz":[1,-1],"nx":[4,13],"ny":[2,21],"nz":[2,0]},{"size":2,"px":[3,4],"py":[7,0],"pz":[1,-1],"nx":[4,2],"ny":[7,5],"nz":[1,2]},{"size":4,"px":[4,14,3,11],"py":[3,23,2,5],"pz":[2,0,2,0],"nx":[7,8,2,16],"ny":[8,0,1,15],"nz":[0,-1,-1,-1]},{"size":2,"px":[9,8],"py":[0,0],"pz":[0,0],"nx":[2,2],"ny":[3,5],"nz":[2,2]}],"alpha":[-1.957970e+00,1.957970e+00,-1.225984e+00,1.225984e+00,-8.310246e-01,8.310246e-01,-8.315741e-01,8.315741e-01,-7.973616e-01,7.973616e-01,-7.661959e-01,7.661959e-01,-6.042118e-01,6.042118e-01,-6.506833e-01,6.506833e-01,-4.808219e-01,4.808219e-01,-6.079504e-01,6.079504e-01,-5.163994e-01,5.163994e-01,-5.268142e-01,5.268142e-01,-4.935685e-01,4.935685e-01,-4.427544e-01,4.427544e-01,-4.053949e-01,4.053949e-01,-4.701274e-01,4.701274e-01,-4.387648e-01,4.387648e-01,-4.305499e-01,4.305499e-01,-4.042607e-01,4.042607e-01,-4.372088e-01,4.372088e-01]},{"count":22,"threshold":-5.679317e+00,"feature":[{"size":5,"px":[11,3,17,14,13],"py":[4,0,13,2,3],"pz":[0,2,0,0,0],"nx":[7,4,14,23,11],"ny":[8,4,8,4,0],"nz":[1,1,0,0,1]},{"size":5,"px":[7,12,6,12,12],"py":[12,8,3,10,9],"pz":[0,0,1,0,0],"nx":[4,9,8,15,15],"ny":[4,8,4,8,8],"nz":[1,0,1,0,-1]},{"size":3,"px":[4,2,10],"py":[1,4,1],"pz":[1,2,0],"nx":[2,3,8],"ny":[5,4,4],"nz":[2,1,-1]},{"size":5,"px":[3,17,6,6,16],"py":[2,12,4,14,12],"pz":[2,0,1,0,0],"nx":[8,3,7,5,15],"ny":[4,4,4,4,8],"nz":[1,1,-1,-1,-1]},{"size":5,"px":[5,6,7,4,8],"py":[3,3,3,1,3],"pz":[0,0,0,1,0],"nx":[0,0,0,0,1],"ny":[5,4,3,2,0],"nz":[0,0,0,0,0]},{"size":3,"px":[18,9,0],"py":[14,7,0],"pz":[0,1,-1],"nx":[8,14,8],"ny":[10,9,4],"nz":[1,0,1]},{"size":2,"px":[9,5],"py":[18,13],"pz":[0,0],"nx":[10,3],"ny":[16,4],"nz":[0,-1]},{"size":5,"px":[11,11,11,11,6],"py":[10,12,11,13,6],"pz":[0,0,0,0,-1],"nx":[5,21,22,22,22],"ny":[4,22,17,19,18],"nz":[2,0,0,0,0]},{"size":4,"px":[8,9,15,4],"py":[7,7,23,4],"pz":[1,1,0,2],"nx":[8,5,0,3],"ny":[4,18,4,9],"nz":[1,-1,-1,-1]},{"size":5,"px":[11,10,12,11,11],"py":[4,4,4,5,5],"pz":[0,0,0,0,-1],"nx":[4,6,8,2,8],"ny":[4,9,9,2,4],"nz":[1,1,0,2,1]},{"size":5,"px":[2,2,3,3,4],"py":[10,9,14,13,15],"pz":[1,1,0,0,0],"nx":[0,0,0,0,0],"ny":[5,9,10,19,18],"nz":[2,1,1,0,-1]},{"size":2,"px":[11,11],"py":[13,12],"pz":[0,0],"nx":[9,2],"ny":[15,2],"nz":[0,-1]},{"size":5,"px":[2,4,3,3,4],"py":[5,11,6,9,12],"pz":[1,0,1,0,0],"nx":[6,2,11,11,0],"ny":[9,1,5,20,18],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[18,9,17,19,16],"py":[2,0,2,2,1],"pz":[0,1,0,0,0],"nx":[22,23,11,23,23],"ny":[0,2,0,1,1],"nz":[0,0,1,0,-1]},{"size":5,"px":[5,5,6,7,6],"py":[17,16,15,23,22],"pz":[0,0,0,0,0],"nx":[7,6,2,5,23],"ny":[8,1,2,3,1],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[12,12,11,10,6],"py":[14,13,18,4,22],"pz":[0,-1,-1,-1,-1],"nx":[3,2,4,1,2],"ny":[19,4,23,13,16],"nz":[0,0,0,0,0]},{"size":4,"px":[11,16,11,17],"py":[7,11,8,12],"pz":[0,0,0,0],"nx":[7,14,10,4],"ny":[4,7,10,4],"nz":[1,0,-1,-1]},{"size":2,"px":[3,3],"py":[8,7],"pz":[1,1],"nx":[4,2],"ny":[10,2],"nz":[1,-1]},{"size":2,"px":[3,9],"py":[0,1],"pz":[1,0],"nx":[4,5],"ny":[1,0],"nz":[0,0]},{"size":2,"px":[14,16],"py":[3,3],"pz":[0,0],"nx":[9,14],"ny":[4,21],"nz":[1,0]},{"size":2,"px":[9,1],"py":[7,1],"pz":[1,-1],"nx":[8,9],"ny":[7,4],"nz":[1,1]},{"size":2,"px":[1,0],"py":[8,3],"pz":[0,2],"nx":[20,0],"ny":[3,3],"nz":[0,-1]}],"alpha":[-1.581077e+00,1.581077e+00,-1.389689e+00,1.389689e+00,-8.733094e-01,8.733094e-01,-8.525177e-01,8.525177e-01,-7.416304e-01,7.416304e-01,-6.609002e-01,6.609002e-01,-7.119043e-01,7.119043e-01,-6.204438e-01,6.204438e-01,-6.638519e-01,6.638519e-01,-5.518876e-01,5.518876e-01,-4.898991e-01,4.898991e-01,-5.508243e-01,5.508243e-01,-4.635525e-01,4.635525e-01,-5.163159e-01,5.163159e-01,-4.495338e-01,4.495338e-01,-4.515036e-01,4.515036e-01,-5.130473e-01,5.130473e-01,-4.694233e-01,4.694233e-01,-4.022514e-01,4.022514e-01,-4.055690e-01,4.055690e-01,-4.151817e-01,4.151817e-01,-3.352302e-01,3.352302e-01]},{"count":32,"threshold":-5.363782e+00,"feature":[{"size":5,"px":[12,9,6,8,14],"py":[4,2,13,3,3],"pz":[0,0,0,0,0],"nx":[0,15,0,9,5],"ny":[2,7,3,8,8],"nz":[0,0,0,0,1]},{"size":5,"px":[13,16,3,6,11],"py":[3,13,1,4,3],"pz":[0,0,2,1,0],"nx":[7,4,8,14,14],"ny":[4,4,4,8,8],"nz":[1,1,1,0,-1]},{"size":5,"px":[10,19,18,19,19],"py":[6,13,13,12,12],"pz":[1,0,0,0,-1],"nx":[23,5,23,23,11],"ny":[12,2,13,14,8],"nz":[0,2,0,0,1]},{"size":5,"px":[12,12,12,12,6],"py":[11,13,12,10,6],"pz":[0,0,0,0,1],"nx":[6,8,3,9,9],"ny":[8,4,4,4,4],"nz":[1,1,1,1,-1]},{"size":5,"px":[5,3,5,8,11],"py":[12,8,3,11,8],"pz":[0,1,1,0,0],"nx":[4,0,1,1,9],"ny":[4,3,4,3,4],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[13,3,12,14,12],"py":[1,0,1,2,3],"pz":[0,2,0,0,0],"nx":[7,9,8,4,4],"ny":[5,4,10,2,2],"nz":[1,1,1,2,-1]},{"size":5,"px":[18,16,12,15,8],"py":[12,23,7,11,8],"pz":[0,0,0,0,1],"nx":[8,6,10,12,4],"ny":[4,4,10,6,3],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[4,4,5,2,2],"py":[13,14,14,7,7],"pz":[0,0,0,1,-1],"nx":[0,0,0,0,1],"ny":[15,4,14,13,17],"nz":[0,2,0,0,0]},{"size":2,"px":[9,9],"py":[7,7],"pz":[1,-1],"nx":[4,7],"ny":[5,8],"nz":[2,1]},{"size":5,"px":[3,4,6,5,4],"py":[2,2,14,6,9],"pz":[1,1,0,1,1],"nx":[23,23,23,23,11],"ny":[0,3,2,1,0],"nz":[0,0,0,0,-1]},{"size":3,"px":[10,2,3],"py":[23,4,7],"pz":[0,2,1],"nx":[10,21,23],"ny":[21,9,2],"nz":[0,-1,-1]},{"size":5,"px":[20,21,21,10,12],"py":[13,12,8,8,12],"pz":[0,0,0,1,0],"nx":[8,16,3,3,11],"ny":[4,8,4,3,0],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[2,21],"py":[4,12],"pz":[2,-1],"nx":[2,3],"ny":[5,4],"nz":[2,1]},{"size":5,"px":[8,5,6,8,7],"py":[0,2,1,1,1],"pz":[0,0,0,0,0],"nx":[3,2,2,2,2],"ny":[0,0,1,2,2],"nz":[0,0,0,0,-1]},{"size":5,"px":[11,2,2,11,10],"py":[10,12,8,11,12],"pz":[0,0,0,0,0],"nx":[3,5,2,4,2],"ny":[4,1,4,2,2],"nz":[1,-1,-1,-1,-1]},{"size":4,"px":[15,16,8,17],"py":[2,1,0,2],"pz":[0,0,1,0],"nx":[19,20,0,8],"ny":[1,2,11,10],"nz":[0,0,-1,-1]},{"size":2,"px":[17,16],"py":[12,12],"pz":[0,0],"nx":[8,9],"ny":[5,1],"nz":[1,-1]},{"size":4,"px":[11,11,0,0],"py":[12,13,0,0],"pz":[0,0,-1,-1],"nx":[10,10,9,10],"ny":[10,12,13,11],"nz":[0,0,0,0]},{"size":3,"px":[11,10,8],"py":[5,2,6],"pz":[0,-1,-1],"nx":[8,12,4],"ny":[4,17,4],"nz":[1,0,1]},{"size":5,"px":[10,21,10,20,20],"py":[11,13,7,13,14],"pz":[1,0,1,0,0],"nx":[23,23,11,23,17],"ny":[23,22,11,21,21],"nz":[0,0,1,-1,-1]},{"size":2,"px":[4,7],"py":[3,9],"pz":[2,1],"nx":[9,23],"ny":[4,22],"nz":[1,-1]},{"size":4,"px":[3,2,2,5],"py":[11,5,4,20],"pz":[1,2,2,0],"nx":[4,23,11,23],"ny":[10,22,11,21],"nz":[1,-1,-1,-1]},{"size":2,"px":[7,5],"py":[13,4],"pz":[0,-1],"nx":[4,4],"ny":[8,6],"nz":[1,1]},{"size":2,"px":[2,5],"py":[4,9],"pz":[2,1],"nx":[10,10],"ny":[16,16],"nz":[0,-1]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,0],"ny":[4,0],"nz":[1,-1]},{"size":5,"px":[7,3,12,13,6],"py":[11,5,23,23,7],"pz":[1,2,0,0,1],"nx":[1,0,0,0,0],"ny":[23,20,19,21,21],"nz":[0,0,0,0,-1]},{"size":5,"px":[0,0,0,0,0],"py":[10,9,6,13,13],"pz":[0,0,1,0,-1],"nx":[8,8,4,4,9],"ny":[4,11,5,4,5],"nz":[1,1,2,2,1]},{"size":2,"px":[9,18],"py":[8,15],"pz":[1,0],"nx":[15,4],"ny":[15,2],"nz":[0,-1]},{"size":2,"px":[5,13],"py":[6,17],"pz":[1,-1],"nx":[1,2],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[19,10,20,18,18],"py":[2,0,2,2,2],"pz":[0,1,0,0,-1],"nx":[22,23,22,11,23],"ny":[1,3,0,1,2],"nz":[0,0,0,1,0]},{"size":5,"px":[4,2,2,2,6],"py":[7,2,5,4,14],"pz":[1,2,2,2,0],"nx":[16,7,9,15,23],"ny":[8,0,3,11,2],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[10,10,9,9,5],"py":[2,0,0,1,0],"pz":[0,0,0,0,1],"nx":[3,2,3,2,2],"ny":[11,3,9,5,5],"nz":[1,2,1,2,-1]}],"alpha":[-1.490426e+00,1.490426e+00,-1.214280e+00,1.214280e+00,-8.124863e-01,8.124863e-01,-7.307594e-01,7.307594e-01,-7.377259e-01,7.377259e-01,-5.982859e-01,5.982859e-01,-6.451736e-01,6.451736e-01,-6.117417e-01,6.117417e-01,-5.438949e-01,5.438949e-01,-4.563701e-01,4.563701e-01,-4.975362e-01,4.975362e-01,-4.707373e-01,4.707373e-01,-5.013868e-01,5.013868e-01,-5.139018e-01,5.139018e-01,-4.728007e-01,4.728007e-01,-4.839748e-01,4.839748e-01,-4.852528e-01,4.852528e-01,-5.768956e-01,5.768956e-01,-3.635091e-01,3.635091e-01,-4.190090e-01,4.190090e-01,-3.854715e-01,3.854715e-01,-3.409591e-01,3.409591e-01,-3.440222e-01,3.440222e-01,-3.375895e-01,3.375895e-01,-3.367032e-01,3.367032e-01,-3.708106e-01,3.708106e-01,-3.260956e-01,3.260956e-01,-3.657681e-01,3.657681e-01,-3.518800e-01,3.518800e-01,-3.845758e-01,3.845758e-01,-2.832236e-01,2.832236e-01,-2.865156e-01,2.865156e-01]},{"count":45,"threshold":-5.479836e+00,"feature":[{"size":5,"px":[15,6,17,6,9],"py":[2,13,13,4,3],"pz":[0,0,0,1,0],"nx":[3,9,4,8,14],"ny":[5,8,4,4,8],"nz":[2,0,1,1,0]},{"size":5,"px":[9,8,11,6,7],"py":[1,2,3,14,2],"pz":[0,0,0,0,0],"nx":[0,0,4,0,0],"ny":[4,2,4,1,0],"nz":[0,0,1,0,0]},{"size":5,"px":[2,2,11,11,11],"py":[2,4,10,8,6],"pz":[2,2,0,0,0],"nx":[8,4,3,23,23],"ny":[4,4,4,16,18],"nz":[1,1,-1,-1,-1]},{"size":5,"px":[18,16,17,15,9],"py":[2,2,2,2,1],"pz":[0,0,0,0,1],"nx":[22,22,21,23,23],"ny":[1,2,0,5,4],"nz":[0,0,0,0,0]},{"size":5,"px":[15,3,17,18,6],"py":[11,2,11,11,4],"pz":[0,2,0,0,1],"nx":[3,8,1,4,23],"ny":[4,4,3,9,4],"nz":[1,1,-1,-1,-1]},{"size":2,"px":[4,5],"py":[4,0],"pz":[2,-1],"nx":[7,4],"ny":[8,5],"nz":[1,2]},{"size":2,"px":[11,5],"py":[12,5],"pz":[0,-1],"nx":[4,9],"ny":[10,15],"nz":[1,0]},{"size":4,"px":[2,2,7,1],"py":[7,7,3,4],"pz":[1,-1,-1,-1],"nx":[0,2,1,2],"ny":[6,20,14,16],"nz":[1,0,0,0]},{"size":5,"px":[14,12,12,13,9],"py":[23,5,6,5,7],"pz":[0,0,0,0,1],"nx":[8,18,2,8,14],"ny":[4,9,0,12,7],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[3,10,13,11,9],"py":[0,3,2,3,2],"pz":[2,0,0,0,0],"nx":[3,11,22,22,22],"ny":[2,6,15,2,0],"nz":[2,1,0,0,0]},{"size":5,"px":[8,7,5,8,5],"py":[23,12,12,12,13],"pz":[0,0,0,0,0],"nx":[3,18,3,1,22],"ny":[4,4,4,2,0],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[22,22,22,21,22],"py":[9,11,10,14,12],"pz":[0,0,0,0,0],"nx":[23,23,11,1,22],"ny":[23,23,11,2,0],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[9,3],"py":[18,7],"pz":[0,1],"nx":[10,8],"ny":[16,19],"nz":[0,-1]},{"size":5,"px":[10,12,11,6,6],"py":[4,4,4,2,2],"pz":[0,0,0,1,-1],"nx":[3,8,7,8,4],"ny":[5,4,4,10,4],"nz":[2,1,1,0,1]},{"size":4,"px":[12,12,4,15],"py":[13,12,0,11],"pz":[0,0,-1,-1],"nx":[13,14,13,14],"ny":[9,12,10,13],"nz":[0,0,0,0]},{"size":2,"px":[4,4],"py":[3,3],"pz":[2,-1],"nx":[9,4],"ny":[4,2],"nz":[1,2]},{"size":3,"px":[9,7,0],"py":[7,5,5],"pz":[1,-1,-1],"nx":[4,15,9],"ny":[5,14,9],"nz":[2,0,1]},{"size":5,"px":[15,20,7,10,16],"py":[17,12,6,4,23],"pz":[0,0,1,1,0],"nx":[1,2,2,1,1],"ny":[3,0,1,2,2],"nz":[0,0,0,0,-1]},{"size":5,"px":[2,1,1,11,2],"py":[16,4,5,12,14],"pz":[0,1,1,0,0],"nx":[4,6,3,19,1],"ny":[4,2,5,19,2],"nz":[1,-1,-1,-1,-1]},{"size":3,"px":[15,14,14],"py":[1,1,0],"pz":[0,0,0],"nx":[4,8,4],"ny":[3,4,2],"nz":[2,1,2]},{"size":5,"px":[2,3,1,2,7],"py":[8,12,4,9,13],"pz":[1,0,2,1,0],"nx":[1,1,0,0,0],"ny":[21,20,18,17,9],"nz":[0,0,0,0,1]},{"size":5,"px":[17,15,17,16,16],"py":[12,12,22,23,12],"pz":[0,0,0,0,0],"nx":[7,3,16,1,0],"ny":[8,6,8,3,9],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[9,17,18,18,18],"py":[6,12,12,13,13],"pz":[1,0,0,0,-1],"nx":[23,23,20,11,11],"ny":[12,13,23,7,8],"nz":[0,0,0,1,1]},{"size":2,"px":[2,4],"py":[4,7],"pz":[2,1],"nx":[4,4],"ny":[10,5],"nz":[1,-1]},{"size":4,"px":[4,22,19,12],"py":[5,8,14,9],"pz":[2,0,0,0],"nx":[8,4,4,2],"ny":[4,4,1,2],"nz":[1,-1,-1,-1]},{"size":2,"px":[3,21],"py":[7,14],"pz":[1,-1],"nx":[4,2],"ny":[7,2],"nz":[1,2]},{"size":3,"px":[7,4,17],"py":[3,1,6],"pz":[0,1,-1],"nx":[3,4,5],"ny":[0,2,1],"nz":[1,0,0]},{"size":4,"px":[15,7,14,0],"py":[3,1,3,7],"pz":[0,1,0,-1],"nx":[8,18,17,18],"ny":[0,1,1,2],"nz":[1,0,0,0]},{"size":5,"px":[12,12,12,12,6],"py":[10,11,12,13,6],"pz":[0,0,0,0,-1],"nx":[8,15,15,4,8],"ny":[10,10,9,2,4],"nz":[0,0,0,2,1]},{"size":2,"px":[17,12],"py":[13,11],"pz":[0,-1],"nx":[9,8],"ny":[4,10],"nz":[1,1]},{"size":5,"px":[0,0,0,0,0],"py":[10,9,12,11,4],"pz":[0,0,0,0,1],"nx":[8,9,8,9,9],"ny":[10,4,4,5,5],"nz":[1,1,1,1,-1]},{"size":3,"px":[7,0,1],"py":[1,9,8],"pz":[0,-1,-1],"nx":[4,3,3],"ny":[7,15,16],"nz":[0,0,0]},{"size":2,"px":[4,7],"py":[15,23],"pz":[0,0],"nx":[9,18],"ny":[21,3],"nz":[0,-1]},{"size":5,"px":[17,4,19,18,8],"py":[12,3,12,17,6],"pz":[0,2,0,0,1],"nx":[23,23,11,22,16],"ny":[0,1,0,21,-1],"nz":[0,0,-1,-1,-1]},{"size":2,"px":[7,4],"py":[13,5],"pz":[0,-1],"nx":[4,2],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[21,20,10,10,21],"py":[13,14,10,7,11],"pz":[0,0,1,1,0],"nx":[4,4,4,5,5],"ny":[18,17,19,20,20],"nz":[0,0,0,0,-1]},{"size":2,"px":[2,3],"py":[11,13],"pz":[1,0],"nx":[12,4],"ny":[17,17],"nz":[0,-1]},{"size":2,"px":[11,5],"py":[13,1],"pz":[0,-1],"nx":[1,2],"ny":[1,4],"nz":[2,1]},{"size":2,"px":[15,7],"py":[17,7],"pz":[0,1],"nx":[14,4],"ny":[15,3],"nz":[0,-1]},{"size":2,"px":[3,11],"py":[3,8],"pz":[2,0],"nx":[13,13],"ny":[9,8],"nz":[0,0]},{"size":2,"px":[8,3],"py":[11,2],"pz":[0,-1],"nx":[8,4],"ny":[9,5],"nz":[0,1]},{"size":3,"px":[12,6,9],"py":[9,10,11],"pz":[0,-1,-1],"nx":[2,1,5],"ny":[2,1,6],"nz":[2,2,1]},{"size":4,"px":[4,5,5,1],"py":[11,11,11,3],"pz":[1,0,1,2],"nx":[0,0,5,4],"ny":[23,22,0,0],"nz":[0,0,-1,-1]},{"size":5,"px":[15,7,17,15,16],"py":[1,0,2,2,0],"pz":[0,1,0,0,0],"nx":[7,4,7,4,8],"ny":[5,2,4,3,4],"nz":[1,2,1,2,-1]},{"size":2,"px":[6,12],"py":[11,23],"pz":[1,0],"nx":[12,4],"ny":[21,2],"nz":[0,-1]}],"alpha":[-1.535800e+00,1.535800e+00,-8.580514e-01,8.580514e-01,-8.625210e-01,8.625210e-01,-7.177500e-01,7.177500e-01,-6.832222e-01,6.832222e-01,-5.736298e-01,5.736298e-01,-5.028217e-01,5.028217e-01,-5.091788e-01,5.091788e-01,-5.791940e-01,5.791940e-01,-4.924942e-01,4.924942e-01,-5.489055e-01,5.489055e-01,-4.528190e-01,4.528190e-01,-4.748324e-01,4.748324e-01,-4.150403e-01,4.150403e-01,-4.820464e-01,4.820464e-01,-4.840212e-01,4.840212e-01,-3.941872e-01,3.941872e-01,-3.663507e-01,3.663507e-01,-3.814835e-01,3.814835e-01,-3.936426e-01,3.936426e-01,-3.049970e-01,3.049970e-01,-3.604256e-01,3.604256e-01,-3.974041e-01,3.974041e-01,-4.203486e-01,4.203486e-01,-3.174435e-01,3.174435e-01,-3.426336e-01,3.426336e-01,-4.492150e-01,4.492150e-01,-3.538784e-01,3.538784e-01,-3.679703e-01,3.679703e-01,-3.985452e-01,3.985452e-01,-2.884028e-01,2.884028e-01,-2.797264e-01,2.797264e-01,-2.664214e-01,2.664214e-01,-2.484857e-01,2.484857e-01,-2.581492e-01,2.581492e-01,-2.943778e-01,2.943778e-01,-2.315507e-01,2.315507e-01,-2.979337e-01,2.979337e-01,-2.976173e-01,2.976173e-01,-2.847965e-01,2.847965e-01,-2.814763e-01,2.814763e-01,-2.489068e-01,2.489068e-01,-2.632427e-01,2.632427e-01,-3.308292e-01,3.308292e-01,-2.790170e-01,2.790170e-01]},{"count":61,"threshold":-5.239104e+00,"feature":[{"size":5,"px":[8,8,11,15,6],"py":[3,6,5,3,4],"pz":[0,1,0,0,1],"nx":[3,9,14,8,4],"ny":[4,8,8,7,2],"nz":[1,0,0,0,2]},{"size":5,"px":[11,12,10,6,9],"py":[3,3,2,13,2],"pz":[0,0,0,0,0],"nx":[0,0,5,2,2],"ny":[13,1,8,5,2],"nz":[0,1,1,2,2]},{"size":5,"px":[11,5,11,11,4],"py":[9,13,10,11,6],"pz":[0,0,0,0,1],"nx":[4,15,9,3,3],"ny":[5,8,9,4,4],"nz":[1,0,0,1,-1]},{"size":5,"px":[15,16,8,17,17],"py":[1,2,0,2,2],"pz":[0,0,1,0,-1],"nx":[23,23,23,23,23],"ny":[4,0,2,3,1],"nz":[0,0,0,0,0]},{"size":4,"px":[9,18,17,18],"py":[7,13,13,14],"pz":[1,0,0,0],"nx":[9,7,4,8],"ny":[4,10,2,4],"nz":[1,1,2,1]},{"size":5,"px":[12,11,12,12,6],"py":[6,5,14,5,3],"pz":[0,0,0,0,1],"nx":[13,8,14,7,7],"ny":[16,4,7,4,4],"nz":[0,1,0,1,-1]},{"size":5,"px":[12,6,3,7,12],"py":[7,12,7,11,8],"pz":[0,0,1,0,0],"nx":[16,4,4,4,7],"ny":[8,4,4,4,4],"nz":[0,1,-1,-1,-1]},{"size":5,"px":[6,4,5,3,3],"py":[2,3,2,0,0],"pz":[0,0,0,1,-1],"nx":[1,0,1,0,0],"ny":[0,3,1,1,2],"nz":[0,0,0,1,0]},{"size":2,"px":[15,9],"py":[11,6],"pz":[0,1],"nx":[14,5],"ny":[9,11],"nz":[0,-1]},{"size":5,"px":[10,19,19,10,20],"py":[7,20,14,6,12],"pz":[1,0,0,1,0],"nx":[23,22,11,23,23],"ny":[21,23,9,20,20],"nz":[0,0,1,0,-1]},{"size":5,"px":[1,1,5,1,1],"py":[8,6,6,9,4],"pz":[0,1,1,0,2],"nx":[3,3,3,2,5],"ny":[4,4,2,5,4],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[13,12,3,11,11],"py":[2,2,0,1,2],"pz":[0,0,2,0,0],"nx":[3,6,8,4,3],"ny":[2,9,4,4,5],"nz":[2,1,1,1,-1]},{"size":3,"px":[12,12,6],"py":[11,12,9],"pz":[0,0,-1],"nx":[2,1,9],"ny":[6,1,14],"nz":[0,2,0]},{"size":5,"px":[6,3,17,16,16],"py":[4,2,14,23,13],"pz":[1,2,0,0,0],"nx":[8,10,21,5,1],"ny":[4,10,11,0,0],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[5,6,1,3,3],"py":[15,14,4,7,7],"pz":[0,0,2,1,-1],"nx":[1,0,0,1,1],"ny":[5,8,7,18,17],"nz":[2,1,1,0,0]},{"size":4,"px":[6,12,5,3],"py":[6,12,2,7],"pz":[1,-1,-1,-1],"nx":[14,13,13,7],"ny":[12,10,9,8],"nz":[0,0,0,1]},{"size":2,"px":[3,6],"py":[7,15],"pz":[1,0],"nx":[3,3],"ny":[4,2],"nz":[1,-1]},{"size":4,"px":[11,10,12,2],"py":[18,18,18,3],"pz":[0,0,0,2],"nx":[11,17,4,16],"ny":[16,4,4,21],"nz":[0,-1,-1,-1]},{"size":5,"px":[9,8,8,5,2],"py":[4,4,4,2,3],"pz":[0,0,-1,-1,-1],"nx":[2,2,4,4,2],"ny":[1,2,10,5,4],"nz":[2,2,1,1,2]},{"size":4,"px":[8,18,14,18],"py":[7,16,23,15],"pz":[1,0,0,0],"nx":[14,3,1,0],"ny":[21,1,9,3],"nz":[0,-1,-1,-1]},{"size":2,"px":[12,3],"py":[9,5],"pz":[0,2],"nx":[8,1],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[9,9],"py":[1,1],"pz":[1,-1],"nx":[19,20],"ny":[1,2],"nz":[0,0]},{"size":3,"px":[10,10,10],"py":[6,6,8],"pz":[1,-1,-1],"nx":[22,21,22],"ny":[13,18,12],"nz":[0,0,0]},{"size":2,"px":[2,2],"py":[4,1],"pz":[2,-1],"nx":[2,4],"ny":[5,4],"nz":[2,1]},{"size":5,"px":[21,21,21,21,21],"py":[19,17,18,15,16],"pz":[0,0,0,0,0],"nx":[11,21,6,1,21],"ny":[17,1,10,0,2],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[7,3,4,4,4],"py":[23,13,14,16,13],"pz":[0,0,0,0,0],"nx":[21,22,22,22,22],"ny":[23,21,20,19,19],"nz":[0,0,0,0,-1]},{"size":2,"px":[11,8],"py":[6,6],"pz":[0,1],"nx":[8,4],"ny":[4,2],"nz":[1,-1]},{"size":5,"px":[23,23,11,23,23],"py":[8,12,6,11,10],"pz":[0,0,1,0,0],"nx":[4,4,3,8,8],"ny":[3,8,4,4,4],"nz":[1,1,1,1,-1]},{"size":5,"px":[8,9,4,7,10],"py":[2,1,0,2,1],"pz":[0,0,1,0,0],"nx":[5,5,6,4,4],"ny":[1,0,0,2,1],"nz":[0,0,0,0,-1]},{"size":2,"px":[12,2],"py":[13,6],"pz":[0,-1],"nx":[15,9],"ny":[15,4],"nz":[0,1]},{"size":2,"px":[2,4],"py":[4,9],"pz":[2,1],"nx":[3,13],"ny":[4,1],"nz":[1,-1]},{"size":3,"px":[3,6,2],"py":[10,22,4],"pz":[1,0,2],"nx":[4,2,1],"ny":[10,4,3],"nz":[1,-1,-1]},{"size":2,"px":[1,0],"py":[9,7],"pz":[0,1],"nx":[0,0],"ny":[23,22],"nz":[0,0]},{"size":2,"px":[8,7],"py":[0,1],"pz":[0,0],"nx":[4,4],"ny":[8,8],"nz":[1,-1]},{"size":5,"px":[7,4,4,6,3],"py":[8,4,5,5,3],"pz":[1,2,2,1,2],"nx":[1,0,2,0,0],"ny":[1,0,0,2,4],"nz":[0,2,0,1,-1]},{"size":3,"px":[10,4,4],"py":[6,1,5],"pz":[1,-1,-1],"nx":[5,23,22],"ny":[4,13,7],"nz":[2,0,0]},{"size":2,"px":[2,2],"py":[6,5],"pz":[1,1],"nx":[6,0],"ny":[9,2],"nz":[0,-1]},{"size":5,"px":[0,1,1,0,0],"py":[5,18,19,16,6],"pz":[2,0,0,0,1],"nx":[5,9,4,8,8],"ny":[8,7,3,7,7],"nz":[1,0,1,0,-1]},{"size":2,"px":[13,12],"py":[23,23],"pz":[0,0],"nx":[7,6],"ny":[8,10],"nz":[0,-1]},{"size":2,"px":[14,19],"py":[12,8],"pz":[0,0],"nx":[18,5],"ny":[8,11],"nz":[0,-1]},{"size":5,"px":[2,8,6,4,4],"py":[3,23,14,6,9],"pz":[2,0,0,1,1],"nx":[0,0,0,0,1],"ny":[21,20,5,19,23],"nz":[0,0,2,0,0]},{"size":2,"px":[11,22],"py":[4,14],"pz":[0,-1],"nx":[3,8],"ny":[1,4],"nz":[2,1]},{"size":5,"px":[1,1,0,1,1],"py":[6,8,3,12,7],"pz":[1,1,2,0,1],"nx":[21,21,19,10,10],"ny":[14,16,23,9,9],"nz":[0,0,0,1,-1]},{"size":2,"px":[10,3],"py":[23,2],"pz":[0,2],"nx":[10,3],"ny":[21,5],"nz":[0,-1]},{"size":2,"px":[9,9],"py":[7,0],"pz":[1,-1],"nx":[9,9],"ny":[11,10],"nz":[1,1]},{"size":5,"px":[23,11,23,23,23],"py":[18,10,19,20,16],"pz":[0,1,0,0,0],"nx":[3,3,2,3,2],"ny":[15,16,10,17,9],"nz":[0,0,1,0,-1]},{"size":2,"px":[9,14],"py":[7,18],"pz":[1,0],"nx":[7,10],"ny":[8,8],"nz":[1,-1]},{"size":2,"px":[12,5],"py":[6,4],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[4,5],"py":[13,4],"pz":[0,-1],"nx":[4,4],"ny":[17,19],"nz":[0,0]},{"size":3,"px":[2,3,3],"py":[11,17,19],"pz":[1,0,0],"nx":[7,7,4],"ny":[8,8,5],"nz":[1,-1,-1]},{"size":2,"px":[6,6],"py":[6,5],"pz":[1,-1],"nx":[2,9],"ny":[4,12],"nz":[1,0]},{"size":5,"px":[8,8,9,2,2],"py":[18,13,12,3,3],"pz":[0,0,0,2,-1],"nx":[23,11,23,11,11],"ny":[13,6,14,7,8],"nz":[0,1,0,1,1]},{"size":2,"px":[9,11],"py":[6,13],"pz":[1,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[8,10],"py":[0,6],"pz":[1,1],"nx":[9,4],"ny":[6,7],"nz":[1,-1]},{"size":3,"px":[3,10,9],"py":[8,6,0],"pz":[1,-1,-1],"nx":[2,2,2],"ny":[15,16,9],"nz":[0,0,1]},{"size":3,"px":[14,15,0],"py":[2,2,5],"pz":[0,0,-1],"nx":[17,17,18],"ny":[0,1,2],"nz":[0,0,0]},{"size":2,"px":[11,5],"py":[14,1],"pz":[0,-1],"nx":[10,9],"ny":[12,14],"nz":[0,0]},{"size":2,"px":[8,8],"py":[7,8],"pz":[1,1],"nx":[8,4],"ny":[4,4],"nz":[1,-1]},{"size":5,"px":[0,0,0,0,0],"py":[19,18,10,5,20],"pz":[0,0,1,2,0],"nx":[4,8,2,4,4],"ny":[4,15,5,10,10],"nz":[1,0,2,1,-1]},{"size":2,"px":[7,0],"py":[13,18],"pz":[0,-1],"nx":[4,3],"ny":[4,4],"nz":[1,1]},{"size":5,"px":[23,22,22,11,22],"py":[16,13,7,6,14],"pz":[0,0,0,1,0],"nx":[13,7,15,14,14],"ny":[6,3,7,6,6],"nz":[0,1,0,0,-1]}],"alpha":[-1.428861e+00,1.428861e+00,-8.591837e-01,8.591837e-01,-7.734305e-01,7.734305e-01,-6.534460e-01,6.534460e-01,-6.262547e-01,6.262547e-01,-5.231782e-01,5.231782e-01,-4.984303e-01,4.984303e-01,-4.913187e-01,4.913187e-01,-4.852198e-01,4.852198e-01,-4.906681e-01,4.906681e-01,-4.126248e-01,4.126248e-01,-4.590814e-01,4.590814e-01,-4.653825e-01,4.653825e-01,-4.179600e-01,4.179600e-01,-4.357392e-01,4.357392e-01,-4.087982e-01,4.087982e-01,-4.594812e-01,4.594812e-01,-4.858794e-01,4.858794e-01,-3.713580e-01,3.713580e-01,-3.894534e-01,3.894534e-01,-3.127168e-01,3.127168e-01,-4.012654e-01,4.012654e-01,-3.370552e-01,3.370552e-01,-3.534712e-01,3.534712e-01,-3.843450e-01,3.843450e-01,-2.688805e-01,2.688805e-01,-3.500203e-01,3.500203e-01,-2.827120e-01,2.827120e-01,-3.742119e-01,3.742119e-01,-3.219074e-01,3.219074e-01,-2.544953e-01,2.544953e-01,-3.355513e-01,3.355513e-01,-2.672670e-01,2.672670e-01,-2.932047e-01,2.932047e-01,-2.404618e-01,2.404618e-01,-2.354372e-01,2.354372e-01,-2.657955e-01,2.657955e-01,-2.293701e-01,2.293701e-01,-2.708918e-01,2.708918e-01,-2.340181e-01,2.340181e-01,-2.464815e-01,2.464815e-01,-2.944239e-01,2.944239e-01,-2.407960e-01,2.407960e-01,-3.029642e-01,3.029642e-01,-2.684602e-01,2.684602e-01,-2.495078e-01,2.495078e-01,-2.539708e-01,2.539708e-01,-2.989293e-01,2.989293e-01,-2.391309e-01,2.391309e-01,-2.531372e-01,2.531372e-01,-2.500390e-01,2.500390e-01,-2.295077e-01,2.295077e-01,-2.526125e-01,2.526125e-01,-2.337182e-01,2.337182e-01,-1.984756e-01,1.984756e-01,-3.089996e-01,3.089996e-01,-2.589053e-01,2.589053e-01,-2.962490e-01,2.962490e-01,-2.458660e-01,2.458660e-01,-2.515206e-01,2.515206e-01,-2.637299e-01,2.637299e-01]},{"count":80,"threshold":-5.185898e+00,"feature":[{"size":5,"px":[12,17,13,10,15],"py":[9,13,3,3,2],"pz":[0,0,0,0,0],"nx":[8,14,6,9,4],"ny":[10,9,8,8,2],"nz":[1,0,1,0,2]},{"size":5,"px":[3,11,8,10,9],"py":[7,4,3,3,3],"pz":[1,0,0,0,0],"nx":[2,1,5,0,0],"ny":[2,15,8,4,13],"nz":[2,0,1,0,0]},{"size":5,"px":[11,11,11,4,17],"py":[7,9,8,6,11],"pz":[0,0,0,1,0],"nx":[8,8,8,3,0],"ny":[4,8,8,8,13],"nz":[1,0,-1,-1,-1]},{"size":5,"px":[14,15,7,16,16],"py":[3,3,1,3,3],"pz":[0,0,1,0,-1],"nx":[23,22,23,22,22],"ny":[6,2,14,3,4],"nz":[0,0,0,0,0]},{"size":4,"px":[6,4,7,15],"py":[4,2,6,17],"pz":[1,2,1,0],"nx":[3,8,3,14],"ny":[4,4,10,22],"nz":[1,1,-1,-1]},{"size":3,"px":[3,5,22],"py":[7,7,5],"pz":[1,-1,-1],"nx":[2,2,4],"ny":[5,2,7],"nz":[2,2,1]},{"size":5,"px":[7,6,5,6,3],"py":[0,1,2,2,0],"pz":[0,0,0,0,1],"nx":[0,1,1,0,1],"ny":[0,2,1,2,0],"nz":[2,0,0,1,0]},{"size":5,"px":[11,11,11,11,5],"py":[11,10,13,12,6],"pz":[0,0,0,0,-1],"nx":[15,14,5,2,8],"ny":[9,8,10,2,10],"nz":[0,0,1,2,0]},{"size":5,"px":[8,5,6,8,7],"py":[12,12,12,23,12],"pz":[0,0,0,0,0],"nx":[3,17,5,2,8],"ny":[4,0,10,2,10],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[10,10,10,19,20],"py":[8,10,9,15,13],"pz":[1,1,1,0,0],"nx":[23,11,5,23,23],"ny":[20,10,5,19,19],"nz":[0,1,2,0,-1]},{"size":5,"px":[9,13,3,10,12],"py":[2,0,0,1,1],"pz":[0,0,2,0,0],"nx":[3,3,6,7,7],"ny":[5,2,11,4,4],"nz":[2,2,1,1,-1]},{"size":2,"px":[15,7],"py":[17,6],"pz":[0,1],"nx":[14,0],"ny":[16,10],"nz":[0,-1]},{"size":5,"px":[17,15,18,12,19],"py":[22,12,13,7,15],"pz":[0,0,0,0,0],"nx":[8,15,6,1,7],"ny":[4,8,22,5,4],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[10,9,18,19,8],"py":[2,1,3,3,1],"pz":[1,1,0,0,1],"nx":[23,23,23,11,11],"ny":[0,1,2,0,1],"nz":[0,0,0,1,-1]},{"size":5,"px":[12,23,0,1,8],"py":[14,5,0,17,1],"pz":[0,-1,-1,-1,-1],"nx":[8,14,15,18,14],"ny":[10,11,14,19,10],"nz":[1,0,0,0,0]},{"size":2,"px":[4,6],"py":[6,13],"pz":[1,0],"nx":[4,12],"ny":[10,14],"nz":[1,-1]},{"size":5,"px":[5,23,11,23,13],"py":[3,10,4,11,12],"pz":[2,0,1,0,0],"nx":[7,4,9,8,8],"ny":[4,2,4,4,4],"nz":[1,2,1,1,-1]},{"size":3,"px":[9,5,11],"py":[4,2,4],"pz":[0,1,-1],"nx":[5,2,4],"ny":[0,1,2],"nz":[0,2,0]},{"size":5,"px":[5,2,2,5,8],"py":[12,4,4,6,13],"pz":[0,2,1,1,0],"nx":[3,9,4,4,8],"ny":[4,0,2,2,4],"nz":[1,-1,-1,-1,-1]},{"size":3,"px":[9,5,22],"py":[7,4,20],"pz":[1,-1,-1],"nx":[8,19,4],"ny":[4,18,5],"nz":[1,0,2]},{"size":5,"px":[2,3,3,3,3],"py":[10,16,15,14,13],"pz":[1,0,0,0,0],"nx":[0,0,0,1,0],"ny":[10,20,5,23,21],"nz":[1,0,2,0,0]},{"size":2,"px":[12,11],"py":[4,18],"pz":[0,0],"nx":[11,23],"ny":[17,13],"nz":[0,-1]},{"size":2,"px":[17,8],"py":[16,7],"pz":[0,1],"nx":[8,3],"ny":[4,6],"nz":[1,-1]},{"size":5,"px":[13,5,14,12,3],"py":[4,7,4,5,3],"pz":[0,1,0,0,1],"nx":[21,20,21,21,21],"ny":[2,0,4,3,3],"nz":[0,0,0,0,-1]},{"size":4,"px":[20,20,20,10],"py":[21,19,20,8],"pz":[0,0,0,1],"nx":[8,11,0,2],"ny":[10,8,1,3],"nz":[1,-1,-1,-1]},{"size":4,"px":[6,7,12,8],"py":[12,12,8,11],"pz":[0,0,0,0],"nx":[9,5,5,18],"ny":[9,2,0,20],"nz":[0,-1,-1,-1]},{"size":3,"px":[11,5,9],"py":[0,0,0],"pz":[0,1,0],"nx":[2,6,3],"ny":[3,7,4],"nz":[2,0,1]},{"size":5,"px":[18,18,9,17,17],"py":[15,14,7,14,14],"pz":[0,0,1,0,-1],"nx":[21,21,21,22,20],"ny":[15,21,17,14,23],"nz":[0,0,0,0,0]},{"size":5,"px":[9,12,12,7,4],"py":[4,11,12,6,5],"pz":[1,0,0,1,2],"nx":[16,11,9,6,20],"ny":[8,4,11,10,23],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[12,11,10,11,11],"py":[23,4,4,5,23],"pz":[0,0,0,0,0],"nx":[11,11,7,3,20],"ny":[21,21,11,1,23],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[12,1],"py":[12,3],"pz":[0,-1],"nx":[10,10],"ny":[3,2],"nz":[1,1]},{"size":5,"px":[9,4,15,9,9],"py":[8,4,23,7,7],"pz":[1,2,0,1,-1],"nx":[5,3,3,3,2],"ny":[23,19,17,18,15],"nz":[0,0,0,0,0]},{"size":2,"px":[2,0],"py":[16,3],"pz":[0,2],"nx":[9,4],"ny":[15,2],"nz":[0,-1]},{"size":2,"px":[2,3],"py":[3,7],"pz":[2,1],"nx":[3,8],"ny":[4,10],"nz":[1,-1]},{"size":3,"px":[9,4,3],"py":[18,0,14],"pz":[0,-1,-1],"nx":[3,5,2],"ny":[5,8,5],"nz":[2,1,2]},{"size":3,"px":[1,1,10],"py":[2,1,7],"pz":[1,-1,-1],"nx":[0,0,0],"ny":[3,5,1],"nz":[0,0,1]},{"size":4,"px":[11,11,5,2],"py":[12,13,7,3],"pz":[0,0,-1,-1],"nx":[5,10,10,9],"ny":[6,9,10,13],"nz":[1,0,0,0]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[9,1],"ny":[4,3],"nz":[1,-1]},{"size":5,"px":[0,0,1,1,0],"py":[4,10,12,13,5],"pz":[1,0,0,0,1],"nx":[4,4,8,7,7],"ny":[3,2,10,4,4],"nz":[2,2,1,1,-1]},{"size":3,"px":[3,4,3],"py":[1,1,2],"pz":[1,-1,-1],"nx":[4,5,3],"ny":[1,0,2],"nz":[0,0,0]},{"size":2,"px":[9,2],"py":[6,4],"pz":[1,-1],"nx":[8,4],"ny":[6,2],"nz":[1,2]},{"size":5,"px":[12,13,15,16,7],"py":[1,1,2,2,1],"pz":[0,0,0,0,1],"nx":[4,4,4,3,7],"ny":[2,2,4,2,4],"nz":[2,-1,-1,-1,-1]},{"size":5,"px":[9,3,2,11,5],"py":[23,7,4,10,6],"pz":[0,1,2,0,1],"nx":[21,20,11,21,21],"ny":[21,23,8,20,20],"nz":[0,0,1,0,-1]},{"size":4,"px":[12,6,13,12],"py":[7,3,5,6],"pz":[0,1,0,0],"nx":[3,0,5,10],"ny":[4,6,5,1],"nz":[1,-1,-1,-1]},{"size":2,"px":[10,4],"py":[4,0],"pz":[0,-1],"nx":[12,11],"ny":[2,1],"nz":[0,0]},{"size":4,"px":[2,3,22,5],"py":[6,1,18,5],"pz":[1,-1,-1,-1],"nx":[0,0,0,3],"ny":[14,3,12,18],"nz":[0,2,0,0]},{"size":3,"px":[10,20,21],"py":[10,18,15],"pz":[1,0,0],"nx":[15,1,2],"ny":[7,0,8],"nz":[0,-1,-1]},{"size":5,"px":[0,0,0,0,0],"py":[4,7,13,4,6],"pz":[1,1,0,2,1],"nx":[5,9,8,4,4],"ny":[3,7,7,3,3],"nz":[1,0,0,1,-1]},{"size":3,"px":[13,12,14],"py":[2,2,2],"pz":[0,0,0],"nx":[4,4,4],"ny":[2,2,5],"nz":[2,-1,-1]},{"size":5,"px":[5,4,6,2,12],"py":[7,9,7,4,10],"pz":[0,1,0,2,0],"nx":[6,1,2,5,2],"ny":[9,2,4,13,4],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[11,1],"py":[12,5],"pz":[0,-1],"nx":[1,0],"ny":[7,2],"nz":[0,2]},{"size":5,"px":[8,8,1,16,6],"py":[6,6,4,8,11],"pz":[1,-1,-1,-1,-1],"nx":[13,5,4,4,13],"ny":[12,1,2,5,11],"nz":[0,2,2,2,0]},{"size":2,"px":[5,6],"py":[4,14],"pz":[1,0],"nx":[9,5],"ny":[7,1],"nz":[0,-1]},{"size":2,"px":[2,6],"py":[4,14],"pz":[2,0],"nx":[9,2],"ny":[15,1],"nz":[0,-1]},{"size":5,"px":[10,19,20,10,9],"py":[1,2,3,0,0],"pz":[1,0,0,1,-1],"nx":[11,23,23,11,23],"ny":[0,3,1,1,2],"nz":[1,0,0,1,0]},{"size":2,"px":[2,9],"py":[3,12],"pz":[2,0],"nx":[2,6],"ny":[4,6],"nz":[1,-1]},{"size":5,"px":[0,0,0,0,0],"py":[4,10,11,9,9],"pz":[1,0,0,0,-1],"nx":[16,2,17,8,4],"ny":[10,2,9,4,4],"nz":[0,2,0,1,1]},{"size":2,"px":[12,0],"py":[5,4],"pz":[0,-1],"nx":[7,8],"ny":[4,8],"nz":[1,1]},{"size":2,"px":[21,21],"py":[9,10],"pz":[0,0],"nx":[11,8],"ny":[18,8],"nz":[0,-1]},{"size":2,"px":[14,7],"py":[23,9],"pz":[0,1],"nx":[7,13],"ny":[10,4],"nz":[1,-1]},{"size":5,"px":[12,12,12,6,2],"py":[11,13,12,6,4],"pz":[0,0,0,-1,-1],"nx":[0,0,0,0,0],"ny":[14,13,6,12,11],"nz":[0,0,1,0,0]},{"size":2,"px":[8,9],"py":[6,11],"pz":[1,-1],"nx":[15,15],"ny":[11,10],"nz":[0,0]},{"size":4,"px":[4,6,7,2],"py":[8,4,23,7],"pz":[1,-1,-1,-1],"nx":[4,20,19,17],"ny":[0,3,1,1],"nz":[2,0,0,0]},{"size":2,"px":[7,0],"py":[6,0],"pz":[1,-1],"nx":[7,4],"ny":[8,2],"nz":[1,2]},{"size":2,"px":[10,0],"py":[7,0],"pz":[1,-1],"nx":[15,15],"ny":[15,14],"nz":[0,0]},{"size":5,"px":[6,2,5,2,4],"py":[23,7,21,8,16],"pz":[0,1,0,1,0],"nx":[18,2,10,0,11],"ny":[9,3,23,5,3],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[9,9,8,10,4],"py":[0,2,2,1,1],"pz":[0,0,0,0,1],"nx":[4,3,2,2,5],"ny":[7,3,4,2,17],"nz":[0,1,2,2,0]},{"size":2,"px":[10,7],"py":[5,6],"pz":[1,-1],"nx":[11,11],"ny":[6,5],"nz":[1,1]},{"size":5,"px":[11,11,5,6,11],"py":[8,10,5,5,9],"pz":[0,0,1,1,0],"nx":[13,16,11,14,4],"ny":[9,13,11,20,23],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[7,14],"py":[14,22],"pz":[0,-1],"nx":[3,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[4,11],"py":[4,5],"pz":[2,-1],"nx":[2,4],"ny":[5,7],"nz":[2,1]},{"size":2,"px":[1,0],"py":[0,0],"pz":[0,1],"nx":[0,4],"ny":[0,2],"nz":[0,-1]},{"size":5,"px":[11,11,11,4,9],"py":[5,5,2,9,23],"pz":[0,-1,-1,-1,-1],"nx":[11,12,10,9,5],"ny":[2,2,2,2,1],"nz":[0,0,0,0,1]},{"size":3,"px":[16,14,15],"py":[1,1,0],"pz":[0,0,0],"nx":[4,7,4],"ny":[2,4,4],"nz":[2,1,-1]},{"size":2,"px":[5,0],"py":[14,5],"pz":[0,-1],"nx":[2,4],"ny":[5,17],"nz":[2,0]},{"size":5,"px":[18,7,16,19,4],"py":[13,6,23,13,3],"pz":[0,1,0,0,2],"nx":[5,2,3,4,4],"ny":[1,1,4,1,3],"nz":[0,1,0,0,0]},{"size":2,"px":[8,8],"py":[7,6],"pz":[1,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[2,1],"py":[10,4],"pz":[1,2],"nx":[4,4],"ny":[3,3],"nz":[2,-1]},{"size":2,"px":[10,5],"py":[19,1],"pz":[0,-1],"nx":[4,12],"ny":[10,17],"nz":[1,0]},{"size":5,"px":[12,6,2,4,11],"py":[14,4,2,1,5],"pz":[0,-1,-1,-1,-1],"nx":[3,4,3,4,3],"ny":[13,17,14,16,15],"nz":[0,0,0,0,0]}],"alpha":[-1.368326e+00,1.368326e+00,-7.706897e-01,7.706897e-01,-8.378147e-01,8.378147e-01,-6.120624e-01,6.120624e-01,-5.139189e-01,5.139189e-01,-4.759130e-01,4.759130e-01,-5.161374e-01,5.161374e-01,-5.407743e-01,5.407743e-01,-4.216105e-01,4.216105e-01,-4.418693e-01,4.418693e-01,-4.435335e-01,4.435335e-01,-4.052076e-01,4.052076e-01,-4.293050e-01,4.293050e-01,-3.431154e-01,3.431154e-01,-4.231203e-01,4.231203e-01,-3.917100e-01,3.917100e-01,-3.623450e-01,3.623450e-01,-3.202670e-01,3.202670e-01,-3.331602e-01,3.331602e-01,-3.552034e-01,3.552034e-01,-3.784556e-01,3.784556e-01,-3.295428e-01,3.295428e-01,-3.587038e-01,3.587038e-01,-2.861332e-01,2.861332e-01,-3.403258e-01,3.403258e-01,-3.989002e-01,3.989002e-01,-2.631159e-01,2.631159e-01,-3.272156e-01,3.272156e-01,-2.816567e-01,2.816567e-01,-3.125926e-01,3.125926e-01,-3.146982e-01,3.146982e-01,-2.521825e-01,2.521825e-01,-2.434554e-01,2.434554e-01,-3.435378e-01,3.435378e-01,-3.161172e-01,3.161172e-01,-2.805027e-01,2.805027e-01,-3.303579e-01,3.303579e-01,-2.725089e-01,2.725089e-01,-2.575051e-01,2.575051e-01,-3.210646e-01,3.210646e-01,-2.986997e-01,2.986997e-01,-2.408925e-01,2.408925e-01,-2.456291e-01,2.456291e-01,-2.836550e-01,2.836550e-01,-2.469860e-01,2.469860e-01,-2.915900e-01,2.915900e-01,-2.513559e-01,2.513559e-01,-2.433728e-01,2.433728e-01,-2.377905e-01,2.377905e-01,-2.089327e-01,2.089327e-01,-1.978434e-01,1.978434e-01,-3.017699e-01,3.017699e-01,-2.339661e-01,2.339661e-01,-1.932560e-01,1.932560e-01,-2.278285e-01,2.278285e-01,-2.438200e-01,2.438200e-01,-2.216769e-01,2.216769e-01,-1.941995e-01,1.941995e-01,-2.129081e-01,2.129081e-01,-2.270319e-01,2.270319e-01,-2.393942e-01,2.393942e-01,-2.132518e-01,2.132518e-01,-1.867741e-01,1.867741e-01,-2.394237e-01,2.394237e-01,-2.005917e-01,2.005917e-01,-2.445217e-01,2.445217e-01,-2.229078e-01,2.229078e-01,-2.342967e-01,2.342967e-01,-2.481784e-01,2.481784e-01,-2.735603e-01,2.735603e-01,-2.187604e-01,2.187604e-01,-1.677239e-01,1.677239e-01,-2.248867e-01,2.248867e-01,-2.505358e-01,2.505358e-01,-1.867706e-01,1.867706e-01,-1.904305e-01,1.904305e-01,-1.939881e-01,1.939881e-01,-2.249474e-01,2.249474e-01,-1.762483e-01,1.762483e-01,-2.299974e-01,2.299974e-01]},{"count":115,"threshold":-5.151920e+00,"feature":[{"size":5,"px":[7,14,7,10,6],"py":[3,3,12,4,4],"pz":[0,0,0,0,1],"nx":[14,3,14,9,3],"ny":[7,4,8,8,5],"nz":[0,1,0,0,2]},{"size":5,"px":[13,18,16,17,15],"py":[1,13,1,2,0],"pz":[0,0,0,0,0],"nx":[23,23,8,11,22],"ny":[3,4,4,8,0],"nz":[0,0,1,1,0]},{"size":5,"px":[16,6,6,7,12],"py":[12,13,4,12,5],"pz":[0,0,1,0,0],"nx":[0,0,8,4,0],"ny":[0,2,4,4,2],"nz":[0,0,1,1,-1]},{"size":3,"px":[12,13,7],"py":[13,18,6],"pz":[0,0,1],"nx":[13,5,6],"ny":[16,3,8],"nz":[0,-1,-1]},{"size":5,"px":[10,12,9,13,11],"py":[3,3,3,3,3],"pz":[0,0,0,0,0],"nx":[3,4,15,4,4],"ny":[2,5,10,4,4],"nz":[2,1,0,1,-1]},{"size":5,"px":[12,12,12,3,12],"py":[7,9,8,3,10],"pz":[0,0,0,2,0],"nx":[4,8,15,9,9],"ny":[4,4,8,8,8],"nz":[1,1,0,0,-1]},{"size":5,"px":[6,3,4,4,2],"py":[22,12,13,14,7],"pz":[0,0,0,0,1],"nx":[2,0,1,1,1],"ny":[23,5,22,21,21],"nz":[0,2,0,0,-1]},{"size":2,"px":[3,3],"py":[8,8],"pz":[1,-1],"nx":[3,4],"ny":[4,10],"nz":[1,1]},{"size":5,"px":[11,11,11,11,0],"py":[10,12,11,13,2],"pz":[0,0,0,-1,-1],"nx":[8,13,13,13,13],"ny":[10,8,9,11,10],"nz":[1,0,0,0,0]},{"size":5,"px":[16,16,15,17,18],"py":[12,23,11,12,12],"pz":[0,0,0,0,0],"nx":[8,8,9,3,13],"ny":[4,4,12,3,10],"nz":[1,-1,-1,-1,-1]},{"size":4,"px":[17,16,6,5],"py":[14,13,4,5],"pz":[0,0,-1,-1],"nx":[8,15,4,7],"ny":[10,14,4,8],"nz":[1,0,2,1]},{"size":5,"px":[20,10,20,21,19],"py":[14,7,13,12,22],"pz":[0,1,0,0,0],"nx":[22,23,11,23,23],"ny":[23,22,11,21,20],"nz":[0,0,1,0,-1]},{"size":4,"px":[12,13,1,18],"py":[14,23,3,5],"pz":[0,-1,-1,-1],"nx":[2,10,5,9],"ny":[2,9,8,14],"nz":[2,0,1,0]},{"size":5,"px":[10,4,7,9,8],"py":[1,0,2,0,1],"pz":[0,1,0,0,0],"nx":[2,3,5,3,3],"ny":[2,4,8,3,3],"nz":[2,1,1,1,-1]},{"size":4,"px":[11,2,2,11],"py":[6,4,5,7],"pz":[0,2,2,0],"nx":[3,0,5,3],"ny":[4,9,8,3],"nz":[1,-1,-1,-1]},{"size":5,"px":[12,10,9,12,12],"py":[11,2,1,10,10],"pz":[0,1,1,0,-1],"nx":[22,11,5,22,23],"ny":[1,1,0,0,3],"nz":[0,1,2,0,0]},{"size":4,"px":[5,10,7,11],"py":[14,3,0,4],"pz":[0,-1,-1,-1],"nx":[4,4,4,4],"ny":[17,18,15,16],"nz":[0,0,0,0]},{"size":5,"px":[2,2,3,2,2],"py":[16,12,20,15,17],"pz":[0,0,0,0,0],"nx":[12,8,4,15,15],"ny":[17,4,4,8,8],"nz":[0,1,1,0,-1]},{"size":5,"px":[12,12,1,6,12],"py":[11,10,3,6,10],"pz":[0,0,-1,-1,-1],"nx":[0,0,1,0,2],"ny":[4,0,2,1,0],"nz":[0,2,0,1,0]},{"size":5,"px":[21,20,21,21,14],"py":[9,16,11,8,12],"pz":[0,0,0,0,0],"nx":[17,6,15,0,2],"ny":[8,23,13,2,0],"nz":[0,-1,-1,-1,-1]},{"size":4,"px":[6,9,9,5],"py":[14,18,23,14],"pz":[0,0,0,0],"nx":[9,5,5,12],"ny":[21,5,3,1],"nz":[0,-1,-1,-1]},{"size":2,"px":[12,13],"py":[4,4],"pz":[0,0],"nx":[4,3],"ny":[4,1],"nz":[1,2]},{"size":5,"px":[7,8,11,4,10],"py":[3,3,2,1,2],"pz":[0,0,0,1,0],"nx":[19,20,19,20,20],"ny":[0,3,1,2,2],"nz":[0,0,0,0,-1]},{"size":2,"px":[9,1],"py":[7,4],"pz":[1,-1],"nx":[4,7],"ny":[5,9],"nz":[2,1]},{"size":5,"px":[11,10,1,5,1],"py":[10,12,6,6,5],"pz":[0,0,1,1,1],"nx":[16,3,2,4,4],"ny":[10,4,2,4,4],"nz":[0,1,2,1,-1]},{"size":2,"px":[15,0],"py":[17,0],"pz":[0,-1],"nx":[7,4],"ny":[8,5],"nz":[1,2]},{"size":5,"px":[8,10,9,9,9],"py":[2,2,2,1,1],"pz":[0,0,0,0,-1],"nx":[4,2,3,3,2],"ny":[0,3,2,1,4],"nz":[0,0,0,0,0]},{"size":4,"px":[11,15,17,16],"py":[8,10,11,11],"pz":[0,0,0,0],"nx":[14,1,1,2],"ny":[9,5,7,0],"nz":[0,-1,-1,-1]},{"size":3,"px":[3,5,9],"py":[8,6,12],"pz":[0,1,0],"nx":[3,4,18],"ny":[4,2,22],"nz":[1,-1,-1]},{"size":5,"px":[6,1,7,3,3],"py":[13,4,13,7,7],"pz":[0,2,0,1,-1],"nx":[0,0,0,0,0],"ny":[16,15,8,13,14],"nz":[0,0,1,0,0]},{"size":2,"px":[5,16],"py":[13,10],"pz":[0,-1],"nx":[3,4],"ny":[4,5],"nz":[1,1]},{"size":5,"px":[5,23,11,23,23],"py":[5,12,4,16,15],"pz":[2,0,1,0,0],"nx":[3,2,4,5,5],"ny":[4,2,4,11,11],"nz":[1,2,1,1,-1]},{"size":4,"px":[10,10,3,23],"py":[7,7,3,16],"pz":[1,-1,-1,-1],"nx":[5,23,11,22],"ny":[4,13,7,16],"nz":[2,0,1,0]},{"size":5,"px":[15,14,13,15,16],"py":[1,0,0,0,1],"pz":[0,0,0,0,0],"nx":[4,9,8,8,8],"ny":[2,4,9,4,4],"nz":[2,1,1,1,-1]},{"size":2,"px":[10,4],"py":[5,5],"pz":[0,-1],"nx":[3,15],"ny":[1,8],"nz":[2,0]},{"size":2,"px":[6,12],"py":[6,9],"pz":[1,0],"nx":[10,10],"ny":[10,10],"nz":[0,-1]},{"size":5,"px":[1,0,0,0,0],"py":[5,4,11,9,12],"pz":[0,1,0,0,0],"nx":[9,8,2,4,7],"ny":[7,7,2,4,7],"nz":[0,0,2,1,0]},{"size":2,"px":[4,8],"py":[4,7],"pz":[2,1],"nx":[9,8],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[5,6],"py":[4,1],"pz":[2,-1],"nx":[8,6],"ny":[7,3],"nz":[1,1]},{"size":5,"px":[8,5,7,6,11],"py":[12,5,13,13,22],"pz":[0,1,0,0,0],"nx":[23,23,23,22,22],"ny":[20,19,21,23,23],"nz":[0,0,0,0,-1]},{"size":2,"px":[3,17],"py":[6,9],"pz":[1,-1],"nx":[3,3],"ny":[10,9],"nz":[1,1]},{"size":2,"px":[14,11],"py":[23,5],"pz":[0,0],"nx":[7,3],"ny":[10,20],"nz":[1,-1]},{"size":2,"px":[3,4],"py":[8,8],"pz":[1,1],"nx":[9,4],"ny":[15,4],"nz":[0,-1]},{"size":2,"px":[2,4],"py":[4,7],"pz":[2,1],"nx":[2,4],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[23,11],"py":[21,10],"pz":[0,1],"nx":[2,3],"ny":[11,14],"nz":[1,0]},{"size":4,"px":[11,11,11,3],"py":[13,12,11,4],"pz":[0,0,0,-1],"nx":[14,13,13,6],"ny":[13,11,10,5],"nz":[0,0,0,1]},{"size":2,"px":[4,7],"py":[3,6],"pz":[2,1],"nx":[9,19],"ny":[4,14],"nz":[1,-1]},{"size":3,"px":[10,5,7],"py":[5,0,6],"pz":[1,-1,-1],"nx":[10,21,5],"ny":[0,5,3],"nz":[1,0,2]},{"size":2,"px":[16,13],"py":[3,15],"pz":[0,-1],"nx":[17,7],"ny":[23,8],"nz":[0,1]},{"size":3,"px":[4,2,2],"py":[15,7,19],"pz":[0,1,-1],"nx":[2,8,4],"ny":[5,14,9],"nz":[2,0,1]},{"size":3,"px":[8,3,6],"py":[10,2,4],"pz":[0,2,1],"nx":[3,8,4],"ny":[4,14,9],"nz":[1,-1,-1]},{"size":2,"px":[14,3],"py":[18,3],"pz":[0,-1],"nx":[12,14],"ny":[17,9],"nz":[0,0]},{"size":3,"px":[7,1,10],"py":[14,10,10],"pz":[0,-1,-1],"nx":[9,6,2],"ny":[13,18,2],"nz":[0,0,2]},{"size":2,"px":[11,8],"py":[13,11],"pz":[0,-1],"nx":[2,4],"ny":[7,18],"nz":[1,0]},{"size":2,"px":[5,4],"py":[21,17],"pz":[0,0],"nx":[9,3],"ny":[5,1],"nz":[1,-1]},{"size":2,"px":[6,6],"py":[4,0],"pz":[0,-1],"nx":[4,3],"ny":[2,0],"nz":[0,1]},{"size":2,"px":[2,1],"py":[1,5],"pz":[0,-1],"nx":[0,1],"ny":[1,0],"nz":[1,0]},{"size":2,"px":[18,1],"py":[13,5],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[0,0,0,0,1],"py":[4,3,2,12,15],"pz":[1,1,2,0,0],"nx":[5,9,4,8,8],"ny":[3,6,3,6,6],"nz":[1,0,1,0,-1]},{"size":2,"px":[2,5],"py":[0,2],"pz":[1,-1],"nx":[2,1],"ny":[0,1],"nz":[0,1]},{"size":4,"px":[7,15,4,20],"py":[8,23,4,8],"pz":[1,0,2,0],"nx":[6,0,3,4],"ny":[9,2,13,6],"nz":[0,-1,-1,-1]},{"size":4,"px":[11,11,10,20],"py":[10,9,11,8],"pz":[0,0,0,-1],"nx":[21,20,21,21],"ny":[18,23,19,17],"nz":[0,0,0,0]},{"size":2,"px":[3,8],"py":[7,5],"pz":[1,-1],"nx":[3,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[5,11],"py":[3,4],"pz":[2,1],"nx":[8,7],"ny":[5,12],"nz":[1,0]},{"size":2,"px":[4,1],"py":[1,3],"pz":[1,-1],"nx":[3,6],"ny":[0,0],"nz":[1,0]},{"size":2,"px":[19,9],"py":[16,8],"pz":[0,1],"nx":[14,6],"ny":[15,1],"nz":[0,-1]},{"size":2,"px":[12,6],"py":[13,5],"pz":[0,-1],"nx":[5,5],"ny":[1,2],"nz":[2,2]},{"size":5,"px":[16,14,4,15,12],"py":[1,1,1,2,1],"pz":[0,0,2,0,0],"nx":[6,4,3,2,10],"ny":[22,8,2,1,7],"nz":[0,1,1,2,0]},{"size":5,"px":[6,8,6,5,5],"py":[1,0,0,1,0],"pz":[0,0,0,0,0],"nx":[4,4,4,4,8],"ny":[4,3,2,5,10],"nz":[2,2,2,2,1]},{"size":2,"px":[9,8],"py":[17,0],"pz":[0,-1],"nx":[2,5],"ny":[5,8],"nz":[2,1]},{"size":2,"px":[8,0],"py":[7,3],"pz":[1,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[10,21],"py":[11,20],"pz":[1,0],"nx":[11,4],"ny":[17,1],"nz":[0,-1]},{"size":5,"px":[5,10,4,17,10],"py":[3,6,3,11,5],"pz":[1,0,1,0,0],"nx":[21,20,9,19,10],"ny":[4,3,0,2,1],"nz":[0,0,1,0,-1]},{"size":2,"px":[23,23],"py":[10,10],"pz":[0,-1],"nx":[23,23],"ny":[21,22],"nz":[0,0]},{"size":5,"px":[9,20,19,20,20],"py":[0,3,1,2,2],"pz":[1,0,0,0,-1],"nx":[11,23,11,23,5],"ny":[1,2,0,1,0],"nz":[1,0,1,0,2]},{"size":3,"px":[6,8,7],"py":[4,10,11],"pz":[1,0,0],"nx":[8,3,4],"ny":[9,4,4],"nz":[0,-1,-1]},{"size":4,"px":[13,13,10,4],"py":[14,23,1,5],"pz":[0,-1,-1,-1],"nx":[15,14,8,8],"ny":[13,12,8,9],"nz":[0,0,1,1]},{"size":2,"px":[11,9],"py":[5,8],"pz":[0,-1],"nx":[7,8],"ny":[7,4],"nz":[0,1]},{"size":5,"px":[4,8,4,7,7],"py":[2,3,3,11,11],"pz":[2,1,2,1,-1],"nx":[0,0,1,0,0],"ny":[4,6,15,3,2],"nz":[1,1,0,2,2]},{"size":2,"px":[6,1],"py":[12,1],"pz":[0,-1],"nx":[1,10],"ny":[2,11],"nz":[2,0]},{"size":5,"px":[0,0,2,3,7],"py":[0,1,4,3,11],"pz":[0,-1,-1,-1,-1],"nx":[9,11,9,6,12],"ny":[2,1,1,0,2],"nz":[0,0,0,1,0]},{"size":2,"px":[10,11],"py":[4,4],"pz":[0,0],"nx":[8,4],"ny":[4,2],"nz":[1,-1]},{"size":5,"px":[1,1,1,1,1],"py":[15,10,19,16,18],"pz":[0,1,0,0,0],"nx":[4,5,3,5,6],"ny":[4,19,9,18,19],"nz":[1,0,1,0,-1]},{"size":5,"px":[12,12,12,12,20],"py":[11,12,13,13,18],"pz":[0,0,0,-1,-1],"nx":[0,0,0,0,0],"ny":[4,2,7,6,12],"nz":[1,2,1,1,0]},{"size":2,"px":[0,0],"py":[9,11],"pz":[0,0],"nx":[10,4],"ny":[5,3],"nz":[1,-1]},{"size":2,"px":[11,8],"py":[9,6],"pz":[0,1],"nx":[13,13],"ny":[10,10],"nz":[0,-1]},{"size":2,"px":[6,3],"py":[5,3],"pz":[1,2],"nx":[3,3],"ny":[5,5],"nz":[2,-1]},{"size":2,"px":[19,9],"py":[10,6],"pz":[0,1],"nx":[4,1],"ny":[2,2],"nz":[2,-1]},{"size":2,"px":[14,4],"py":[19,12],"pz":[0,-1],"nx":[14,8],"ny":[17,10],"nz":[0,1]},{"size":4,"px":[4,2,13,2],"py":[12,6,9,3],"pz":[0,1,-1,-1],"nx":[1,0,1,0],"ny":[16,14,11,15],"nz":[0,0,1,0]},{"size":2,"px":[3,3],"py":[8,7],"pz":[1,1],"nx":[4,4],"ny":[4,8],"nz":[1,-1]},{"size":5,"px":[9,11,12,6,10],"py":[2,1,2,1,2],"pz":[0,0,0,1,0],"nx":[4,6,4,6,2],"ny":[4,0,9,1,8],"nz":[0,0,1,0,1]},{"size":5,"px":[4,4,7,2,2],"py":[19,20,23,8,9],"pz":[0,0,0,1,1],"nx":[7,0,5,6,2],"ny":[10,5,4,1,8],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[18,18,17,18,18],"py":[15,16,14,20,17],"pz":[0,0,0,0,0],"nx":[15,2,2,5,2],"ny":[8,0,2,9,4],"nz":[0,-1,-1,-1,-1]},{"size":4,"px":[13,13,13,18],"py":[11,12,12,20],"pz":[0,0,-1,-1],"nx":[1,3,10,10],"ny":[1,6,12,11],"nz":[2,0,0,0]},{"size":2,"px":[8,9],"py":[0,1],"pz":[1,1],"nx":[19,4],"ny":[2,2],"nz":[0,-1]},{"size":2,"px":[6,3],"py":[4,2],"pz":[1,2],"nx":[8,4],"ny":[4,0],"nz":[1,-1]},{"size":5,"px":[23,11,22,13,13],"py":[8,3,3,12,12],"pz":[0,1,0,0,-1],"nx":[15,7,14,13,8],"ny":[7,3,6,6,3],"nz":[0,1,0,0,1]},{"size":3,"px":[9,11,19],"py":[7,3,0],"pz":[1,-1,-1],"nx":[23,23,11],"ny":[16,12,7],"nz":[0,0,1]},{"size":2,"px":[15,8],"py":[23,7],"pz":[0,-1],"nx":[4,3],"ny":[5,4],"nz":[2,2]},{"size":2,"px":[4,10],"py":[6,13],"pz":[1,-1],"nx":[2,3],"ny":[4,10],"nz":[2,1]},{"size":2,"px":[4,1],"py":[11,2],"pz":[1,2],"nx":[9,2],"ny":[5,2],"nz":[1,-1]},{"size":2,"px":[22,22],"py":[22,21],"pz":[0,0],"nx":[3,0],"ny":[5,3],"nz":[1,-1]},{"size":2,"px":[20,10],"py":[12,6],"pz":[0,1],"nx":[20,10],"ny":[23,11],"nz":[0,-1]},{"size":4,"px":[10,3,3,4],"py":[5,3,4,9],"pz":[0,-1,-1,-1],"nx":[14,4,3,11],"ny":[2,1,1,3],"nz":[0,2,2,0]},{"size":3,"px":[15,15,3],"py":[1,1,4],"pz":[0,-1,-1],"nx":[7,4,4],"ny":[8,2,3],"nz":[1,2,2]},{"size":3,"px":[0,0,0],"py":[3,4,6],"pz":[2,2,1],"nx":[0,21,4],"ny":[23,14,3],"nz":[0,-1,-1]},{"size":5,"px":[4,4,5,3,4],"py":[9,11,8,4,8],"pz":[1,1,1,2,1],"nx":[21,21,10,19,19],"ny":[3,4,1,0,0],"nz":[0,0,1,0,-1]},{"size":4,"px":[21,20,20,21],"py":[18,21,20,17],"pz":[0,0,0,0],"nx":[8,1,4,2],"ny":[10,0,2,4],"nz":[1,-1,-1,-1]},{"size":2,"px":[3,6],"py":[7,14],"pz":[1,0],"nx":[3,5],"ny":[4,5],"nz":[1,-1]},{"size":3,"px":[12,0,23],"py":[20,2,13],"pz":[0,-1,-1],"nx":[12,2,9],"ny":[19,2,7],"nz":[0,2,0]},{"size":2,"px":[0,6],"py":[22,11],"pz":[0,-1],"nx":[20,18],"ny":[12,23],"nz":[0,0]},{"size":5,"px":[9,15,15,16,8],"py":[2,1,2,2,1],"pz":[1,0,0,0,1],"nx":[1,1,1,1,1],"ny":[16,10,17,18,18],"nz":[0,1,0,0,-1]},{"size":5,"px":[10,5,3,5,8],"py":[14,2,1,4,1],"pz":[0,-1,-1,-1,-1],"nx":[23,23,23,23,23],"ny":[18,15,16,14,17],"nz":[0,0,0,0,0]},{"size":5,"px":[2,2,2,3,2],"py":[16,17,15,20,11],"pz":[0,0,0,0,1],"nx":[8,22,2,1,23],"ny":[20,11,5,0,17],"nz":[0,-1,-1,-1,-1]}],"alpha":[-1.299972e+00,1.299972e+00,-7.630804e-01,7.630804e-01,-5.530378e-01,5.530378e-01,-5.444703e-01,5.444703e-01,-5.207701e-01,5.207701e-01,-5.035143e-01,5.035143e-01,-4.514416e-01,4.514416e-01,-4.897723e-01,4.897723e-01,-5.006264e-01,5.006264e-01,-4.626049e-01,4.626049e-01,-4.375402e-01,4.375402e-01,-3.742565e-01,3.742565e-01,-3.873996e-01,3.873996e-01,-3.715484e-01,3.715484e-01,-3.562480e-01,3.562480e-01,-3.216189e-01,3.216189e-01,-3.983409e-01,3.983409e-01,-3.191891e-01,3.191891e-01,-3.242173e-01,3.242173e-01,-3.528040e-01,3.528040e-01,-3.562318e-01,3.562318e-01,-3.592398e-01,3.592398e-01,-2.557584e-01,2.557584e-01,-2.747951e-01,2.747951e-01,-2.747554e-01,2.747554e-01,-2.980481e-01,2.980481e-01,-2.887670e-01,2.887670e-01,-3.895318e-01,3.895318e-01,-2.786896e-01,2.786896e-01,-2.763841e-01,2.763841e-01,-2.704816e-01,2.704816e-01,-2.075489e-01,2.075489e-01,-3.104773e-01,3.104773e-01,-2.580337e-01,2.580337e-01,-2.448334e-01,2.448334e-01,-3.054279e-01,3.054279e-01,-2.335804e-01,2.335804e-01,-2.972322e-01,2.972322e-01,-2.270521e-01,2.270521e-01,-2.134621e-01,2.134621e-01,-2.261655e-01,2.261655e-01,-2.091024e-01,2.091024e-01,-2.478928e-01,2.478928e-01,-2.468972e-01,2.468972e-01,-1.919746e-01,1.919746e-01,-2.756623e-01,2.756623e-01,-2.629717e-01,2.629717e-01,-2.198653e-01,2.198653e-01,-2.174434e-01,2.174434e-01,-2.193626e-01,2.193626e-01,-1.956262e-01,1.956262e-01,-1.720459e-01,1.720459e-01,-1.781067e-01,1.781067e-01,-1.773484e-01,1.773484e-01,-1.793871e-01,1.793871e-01,-1.973396e-01,1.973396e-01,-2.397262e-01,2.397262e-01,-2.164685e-01,2.164685e-01,-2.214348e-01,2.214348e-01,-2.265941e-01,2.265941e-01,-2.075436e-01,2.075436e-01,-2.244070e-01,2.244070e-01,-2.291992e-01,2.291992e-01,-2.223506e-01,2.223506e-01,-1.639398e-01,1.639398e-01,-1.732374e-01,1.732374e-01,-1.808631e-01,1.808631e-01,-1.860962e-01,1.860962e-01,-1.781604e-01,1.781604e-01,-2.108322e-01,2.108322e-01,-2.386390e-01,2.386390e-01,-1.942083e-01,1.942083e-01,-1.949161e-01,1.949161e-01,-1.953729e-01,1.953729e-01,-2.317591e-01,2.317591e-01,-2.335136e-01,2.335136e-01,-2.282835e-01,2.282835e-01,-2.148716e-01,2.148716e-01,-1.588127e-01,1.588127e-01,-1.566765e-01,1.566765e-01,-1.644839e-01,1.644839e-01,-2.386947e-01,2.386947e-01,-1.704126e-01,1.704126e-01,-2.213945e-01,2.213945e-01,-1.740398e-01,1.740398e-01,-2.451678e-01,2.451678e-01,-2.120524e-01,2.120524e-01,-1.886646e-01,1.886646e-01,-2.824447e-01,2.824447e-01,-1.900364e-01,1.900364e-01,-2.179183e-01,2.179183e-01,-2.257696e-01,2.257696e-01,-2.023404e-01,2.023404e-01,-1.886901e-01,1.886901e-01,-1.850663e-01,1.850663e-01,-2.035414e-01,2.035414e-01,-1.930174e-01,1.930174e-01,-1.898282e-01,1.898282e-01,-1.666640e-01,1.666640e-01,-1.646143e-01,1.646143e-01,-1.543475e-01,1.543475e-01,-1.366289e-01,1.366289e-01,-1.636837e-01,1.636837e-01,-2.547716e-01,2.547716e-01,-1.281869e-01,1.281869e-01,-1.509159e-01,1.509159e-01,-1.447827e-01,1.447827e-01,-1.626126e-01,1.626126e-01,-2.387014e-01,2.387014e-01,-2.571160e-01,2.571160e-01,-1.719175e-01,1.719175e-01,-1.646742e-01,1.646742e-01,-1.717041e-01,1.717041e-01,-2.039217e-01,2.039217e-01,-1.796907e-01,1.796907e-01]},{"count":153,"threshold":-4.971032e+00,"feature":[{"size":5,"px":[14,13,18,10,16],"py":[2,2,13,3,12],"pz":[0,0,0,0,0],"nx":[21,7,14,23,23],"ny":[16,7,8,3,13],"nz":[0,1,0,0,0]},{"size":5,"px":[12,12,12,15,14],"py":[9,10,11,3,3],"pz":[0,0,0,0,0],"nx":[9,9,8,14,3],"ny":[9,8,5,9,5],"nz":[0,0,1,0,2]},{"size":5,"px":[5,11,7,6,8],"py":[12,8,12,12,11],"pz":[0,0,0,0,0],"nx":[8,4,3,9,9],"ny":[4,4,4,9,9],"nz":[1,1,1,0,-1]},{"size":5,"px":[9,8,4,10,6],"py":[2,2,1,3,13],"pz":[0,0,1,0,0],"nx":[1,1,5,1,1],"ny":[2,3,8,4,16],"nz":[0,0,1,0,0]},{"size":5,"px":[3,16,6,17,15],"py":[2,17,4,12,12],"pz":[2,0,1,0,0],"nx":[4,8,15,1,1],"ny":[4,4,8,16,16],"nz":[1,1,-1,-1,-1]},{"size":4,"px":[18,15,8,17],"py":[12,23,6,12],"pz":[0,0,1,0],"nx":[15,4,10,5],"ny":[21,8,14,3],"nz":[0,-1,-1,-1]},{"size":5,"px":[18,17,9,19,19],"py":[3,1,0,3,3],"pz":[0,0,1,0,-1],"nx":[22,11,23,23,23],"ny":[0,1,2,3,4],"nz":[0,1,0,0,0]},{"size":4,"px":[9,5,5,10],"py":[18,15,14,18],"pz":[0,0,0,0],"nx":[10,11,2,0],"ny":[16,7,12,7],"nz":[0,-1,-1,-1]},{"size":2,"px":[2,12],"py":[4,6],"pz":[2,0],"nx":[3,12],"ny":[4,19],"nz":[1,-1]},{"size":5,"px":[3,4,5,2,2],"py":[3,3,3,1,1],"pz":[0,0,0,1,-1],"nx":[0,0,1,0,0],"ny":[3,4,0,1,2],"nz":[0,0,0,1,0]},{"size":5,"px":[12,12,12,8,10],"py":[13,12,12,1,18],"pz":[0,0,-1,-1,-1],"nx":[13,8,7,14,9],"ny":[10,10,7,13,4],"nz":[0,1,1,0,1]},{"size":5,"px":[15,4,12,14,12],"py":[12,3,9,10,8],"pz":[0,2,0,0,0],"nx":[14,7,11,2,9],"ny":[8,4,7,5,4],"nz":[0,1,-1,-1,-1]},{"size":3,"px":[3,9,7],"py":[7,23,15],"pz":[1,-1,-1],"nx":[4,4,2],"ny":[9,7,5],"nz":[1,1,2]},{"size":3,"px":[5,17,5],"py":[3,23,4],"pz":[2,0,2],"nx":[23,2,4],"ny":[23,16,4],"nz":[0,-1,-1]},{"size":5,"px":[4,9,9,10,8],"py":[1,0,1,0,2],"pz":[1,0,0,0,0],"nx":[2,5,4,2,2],"ny":[2,19,11,4,1],"nz":[2,0,1,2,2]},{"size":5,"px":[8,3,8,4,7],"py":[23,9,13,8,16],"pz":[0,1,0,1,0],"nx":[8,2,5,3,2],"ny":[8,15,1,1,1],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[11,5],"py":[14,5],"pz":[0,-1],"nx":[1,9],"ny":[3,13],"nz":[2,0]},{"size":5,"px":[5,8,1,8,6],"py":[12,12,3,23,12],"pz":[0,0,2,0,0],"nx":[1,1,2,1,1],"ny":[22,21,23,20,20],"nz":[0,0,0,0,-1]},{"size":5,"px":[14,21,19,21,20],"py":[13,8,20,10,7],"pz":[0,0,0,0,0],"nx":[16,0,14,23,1],"ny":[8,1,23,10,20],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[15,16,13,14,14],"py":[3,3,3,3,3],"pz":[0,0,0,0,-1],"nx":[18,19,18,9,17],"ny":[2,2,1,1,0],"nz":[0,0,0,1,0]},{"size":2,"px":[17,9],"py":[14,4],"pz":[0,-1],"nx":[9,18],"ny":[4,18],"nz":[1,0]},{"size":2,"px":[21,20],"py":[17,21],"pz":[0,0],"nx":[12,3],"ny":[17,10],"nz":[0,-1]},{"size":2,"px":[2,1],"py":[10,4],"pz":[1,2],"nx":[4,1],"ny":[10,5],"nz":[1,-1]},{"size":5,"px":[7,8,4,9,9],"py":[2,2,0,2,2],"pz":[0,0,1,0,-1],"nx":[5,5,4,6,3],"ny":[0,1,2,0,0],"nz":[0,0,0,0,1]},{"size":2,"px":[2,5],"py":[3,5],"pz":[2,-1],"nx":[3,2],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[0,0,0,0,0],"py":[0,1,3,4,4],"pz":[2,2,1,1,-1],"nx":[20,20,19,20,19],"ny":[21,20,23,19,22],"nz":[0,0,0,0,0]},{"size":2,"px":[9,18],"py":[8,16],"pz":[1,0],"nx":[14,6],"ny":[15,16],"nz":[0,-1]},{"size":3,"px":[3,4,7],"py":[3,3,9],"pz":[2,2,1],"nx":[8,9,7],"ny":[4,11,4],"nz":[1,-1,-1]},{"size":5,"px":[6,14,4,7,7],"py":[4,23,3,6,6],"pz":[1,0,2,1,-1],"nx":[2,0,2,1,3],"ny":[20,4,21,10,23],"nz":[0,2,0,1,0]},{"size":5,"px":[2,4,8,9,10],"py":[3,8,13,23,23],"pz":[2,1,0,0,0],"nx":[10,4,0,3,3],"ny":[21,3,0,3,23],"nz":[0,-1,-1,-1,-1]},{"size":3,"px":[11,10,11],"py":[6,5,5],"pz":[0,0,0],"nx":[14,6,1],"ny":[7,9,5],"nz":[0,1,-1]},{"size":5,"px":[11,11,11,11,6],"py":[11,12,10,13,6],"pz":[0,0,0,0,1],"nx":[9,13,13,13,4],"ny":[4,9,10,11,2],"nz":[1,0,0,0,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,11],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[1,2],"py":[4,11],"pz":[2,0],"nx":[8,8],"ny":[15,15],"nz":[0,-1]},{"size":5,"px":[12,12,13,12,12],"py":[10,11,13,12,12],"pz":[0,0,0,0,-1],"nx":[0,0,0,1,0],"ny":[13,2,12,5,14],"nz":[0,2,0,0,0]},{"size":5,"px":[0,0,0,1,1],"py":[4,3,11,15,13],"pz":[1,2,0,0,0],"nx":[2,3,3,1,0],"ny":[2,4,4,5,14],"nz":[2,1,-1,-1,-1]},{"size":2,"px":[4,11],"py":[12,10],"pz":[0,-1],"nx":[1,2],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[18,8,9,9,9],"py":[15,7,8,10,7],"pz":[0,1,1,1,1],"nx":[22,23,21,22,11],"ny":[20,16,23,19,9],"nz":[0,0,0,0,1]},{"size":5,"px":[14,12,13,14,15],"py":[1,0,0,0,1],"pz":[0,0,0,0,0],"nx":[4,9,4,7,7],"ny":[2,3,1,8,8],"nz":[2,1,2,1,-1]},{"size":2,"px":[13,9],"py":[14,19],"pz":[0,-1],"nx":[6,10],"ny":[0,2],"nz":[1,0]},{"size":2,"px":[13,12],"py":[4,4],"pz":[0,0],"nx":[3,3],"ny":[1,1],"nz":[2,-1]},{"size":3,"px":[14,5,5],"py":[18,3,4],"pz":[0,-1,-1],"nx":[8,7,8],"ny":[4,8,10],"nz":[1,1,1]},{"size":2,"px":[8,18],"py":[6,11],"pz":[1,0],"nx":[9,1],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[16,11],"py":[9,7],"pz":[0,0],"nx":[7,7],"ny":[4,4],"nz":[1,-1]},{"size":5,"px":[23,11,23,11,23],"py":[13,4,12,7,10],"pz":[0,1,0,1,0],"nx":[7,4,8,15,15],"ny":[9,2,4,8,8],"nz":[0,2,1,0,-1]},{"size":2,"px":[6,3],"py":[1,0],"pz":[0,1],"nx":[4,1],"ny":[1,2],"nz":[0,-1]},{"size":2,"px":[5,5],"py":[7,6],"pz":[0,1],"nx":[6,4],"ny":[9,11],"nz":[0,-1]},{"size":4,"px":[5,6,5,5],"py":[8,6,11,6],"pz":[1,1,1,0],"nx":[23,0,4,5],"ny":[0,2,2,1],"nz":[0,-1,-1,-1]},{"size":2,"px":[18,4],"py":[13,3],"pz":[0,-1],"nx":[15,4],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[4,0],"py":[8,0],"pz":[1,-1],"nx":[9,2],"ny":[15,5],"nz":[0,2]},{"size":5,"px":[15,15,16,14,14],"py":[0,1,1,0,0],"pz":[0,0,0,0,-1],"nx":[4,4,8,8,15],"ny":[4,5,4,11,23],"nz":[2,2,1,1,0]},{"size":4,"px":[12,11,3,14],"py":[14,22,1,0],"pz":[0,-1,-1,-1],"nx":[8,15,7,16],"ny":[2,3,1,3],"nz":[1,0,1,0]},{"size":2,"px":[5,12],"py":[6,17],"pz":[1,-1],"nx":[2,1],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[13,12,12,7,7],"py":[5,6,5,14,14],"pz":[0,0,0,0,-1],"nx":[10,3,10,1,10],"ny":[13,8,11,3,10],"nz":[0,0,0,1,0]},{"size":2,"px":[4,4],"py":[15,0],"pz":[0,-1],"nx":[4,4],"ny":[16,17],"nz":[0,0]},{"size":5,"px":[1,4,2,1,2],"py":[4,0,1,1,0],"pz":[1,1,1,2,1],"nx":[4,9,1,5,1],"ny":[3,4,4,5,5],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[10,3],"py":[3,1],"pz":[0,2],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[16,0],"py":[21,0],"pz":[0,-1],"nx":[6,8],"ny":[8,4],"nz":[1,1]},{"size":2,"px":[7,11],"py":[4,18],"pz":[0,-1],"nx":[5,7],"ny":[0,2],"nz":[2,0]},{"size":2,"px":[9,7],"py":[0,3],"pz":[1,-1],"nx":[20,10],"ny":[0,1],"nz":[0,1]},{"size":4,"px":[10,4,1,5],"py":[0,6,8,4],"pz":[1,-1,-1,-1],"nx":[6,15,4,14],"ny":[3,5,1,5],"nz":[1,0,2,0]},{"size":2,"px":[4,4],"py":[3,4],"pz":[2,2],"nx":[9,2],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[8,4],"py":[3,4],"pz":[0,-1],"nx":[8,6],"ny":[2,1],"nz":[0,0]},{"size":2,"px":[2,0],"py":[6,3],"pz":[1,2],"nx":[0,7],"ny":[7,8],"nz":[1,-1]},{"size":2,"px":[10,0],"py":[7,3],"pz":[1,-1],"nx":[15,4],"ny":[14,4],"nz":[0,2]},{"size":4,"px":[3,1,2,2],"py":[20,7,18,17],"pz":[0,1,0,0],"nx":[9,5,5,4],"ny":[5,4,18,4],"nz":[1,-1,-1,-1]},{"size":2,"px":[5,4],"py":[3,1],"pz":[2,-1],"nx":[23,23],"ny":[14,13],"nz":[0,0]},{"size":2,"px":[12,4],"py":[6,1],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":5,"px":[22,22,11,11,11],"py":[12,13,4,6,6],"pz":[0,0,1,1,-1],"nx":[4,4,4,4,3],"ny":[16,15,18,14,11],"nz":[0,0,0,0,1]},{"size":2,"px":[4,10],"py":[0,1],"pz":[1,0],"nx":[2,2],"ny":[2,2],"nz":[2,-1]},{"size":2,"px":[15,6],"py":[4,4],"pz":[0,-1],"nx":[15,4],"ny":[2,1],"nz":[0,2]},{"size":2,"px":[11,2],"py":[10,20],"pz":[0,-1],"nx":[4,9],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[4,19],"py":[3,8],"pz":[2,0],"nx":[8,21],"ny":[4,20],"nz":[1,-1]},{"size":5,"px":[4,6,7,6,2],"py":[6,15,13,14,3],"pz":[1,0,0,0,-1],"nx":[21,22,19,21,10],"ny":[6,12,0,3,2],"nz":[0,0,0,0,1]},{"size":5,"px":[8,12,15,14,13],"py":[0,0,0,0,0],"pz":[1,0,0,0,0],"nx":[4,3,1,3,4],"ny":[19,16,3,15,4],"nz":[0,0,2,0,1]},{"size":2,"px":[3,3],"py":[2,3],"pz":[2,2],"nx":[8,4],"ny":[4,1],"nz":[1,-1]},{"size":4,"px":[0,0,0,5],"py":[10,9,11,21],"pz":[1,1,-1,-1],"nx":[12,4,3,11],"ny":[3,1,1,3],"nz":[0,2,2,0]},{"size":2,"px":[3,1],"py":[0,0],"pz":[1,2],"nx":[1,4],"ny":[2,1],"nz":[1,-1]},{"size":5,"px":[2,5,1,0,1],"py":[14,23,7,5,9],"pz":[0,0,1,1,1],"nx":[0,0,7,9,11],"ny":[23,22,4,9,3],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[8,9],"py":[7,1],"pz":[1,-1],"nx":[8,8],"ny":[8,9],"nz":[1,1]},{"size":2,"px":[11,9],"py":[11,3],"pz":[1,-1],"nx":[3,2],"ny":[14,10],"nz":[0,1]},{"size":4,"px":[2,4,5,4],"py":[8,20,22,16],"pz":[1,0,0,0],"nx":[8,2,11,3],"ny":[7,4,15,4],"nz":[0,-1,-1,-1]},{"size":3,"px":[1,2,3],"py":[2,1,0],"pz":[0,0,0],"nx":[0,0,15],"ny":[1,0,11],"nz":[0,0,-1]},{"size":2,"px":[12,22],"py":[6,7],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":3,"px":[13,0,5],"py":[19,10,2],"pz":[0,-1,-1],"nx":[3,4,6],"ny":[5,5,9],"nz":[2,2,1]},{"size":2,"px":[8,15],"py":[8,22],"pz":[1,0],"nx":[7,4],"ny":[10,7],"nz":[1,-1]},{"size":2,"px":[10,10],"py":[7,6],"pz":[1,1],"nx":[10,1],"ny":[9,0],"nz":[1,-1]},{"size":2,"px":[9,11],"py":[4,3],"pz":[0,-1],"nx":[5,9],"ny":[0,1],"nz":[1,0]},{"size":5,"px":[14,13,14,12,15],"py":[1,2,2,2,2],"pz":[0,0,0,0,0],"nx":[4,8,4,7,4],"ny":[2,4,3,4,4],"nz":[2,1,2,1,-1]},{"size":3,"px":[13,8,2],"py":[14,5,8],"pz":[0,-1,-1],"nx":[6,8,9],"ny":[3,2,2],"nz":[0,0,0]},{"size":3,"px":[3,6,8],"py":[7,4,12],"pz":[1,1,0],"nx":[3,8,9],"ny":[5,2,2],"nz":[1,-1,-1]},{"size":2,"px":[13,4],"py":[16,3],"pz":[0,2],"nx":[13,7],"ny":[15,5],"nz":[0,-1]},{"size":2,"px":[3,0],"py":[7,9],"pz":[1,-1],"nx":[2,8],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[3,6,8,7,7],"py":[0,1,0,0,0],"pz":[1,0,0,0,-1],"nx":[7,9,4,3,4],"ny":[9,7,4,2,2],"nz":[1,1,1,2,2]},{"size":3,"px":[3,4,16],"py":[4,4,6],"pz":[1,2,0],"nx":[2,2,2],"ny":[0,0,1],"nz":[0,-1,-1]},{"size":2,"px":[0,0],"py":[1,0],"pz":[2,2],"nx":[5,5],"ny":[2,2],"nz":[1,-1]},{"size":2,"px":[9,3],"py":[7,20],"pz":[1,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[8,21],"py":[10,18],"pz":[0,-1],"nx":[9,4],"ny":[10,4],"nz":[0,1]},{"size":2,"px":[6,13],"py":[6,23],"pz":[1,-1],"nx":[10,10],"ny":[11,12],"nz":[0,0]},{"size":5,"px":[10,9,5,10,10],"py":[9,13,6,10,10],"pz":[0,0,1,0,-1],"nx":[21,21,21,10,21],"ny":[18,20,19,11,17],"nz":[0,0,0,1,0]},{"size":2,"px":[8,8],"py":[7,6],"pz":[1,1],"nx":[8,1],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[11,4],"py":[14,7],"pz":[0,-1],"nx":[13,13],"ny":[13,11],"nz":[0,0]},{"size":2,"px":[4,4],"py":[4,5],"pz":[2,2],"nx":[12,5],"ny":[16,2],"nz":[0,-1]},{"size":3,"px":[1,3,20],"py":[3,9,2],"pz":[2,-1,-1],"nx":[0,0,0],"ny":[7,4,13],"nz":[1,2,0]},{"size":2,"px":[0,0],"py":[4,2],"pz":[1,2],"nx":[1,0],"ny":[4,4],"nz":[1,-1]},{"size":3,"px":[8,9,11],"py":[2,1,2],"pz":[0,0,0],"nx":[2,2,0],"ny":[2,2,13],"nz":[2,-1,-1]},{"size":2,"px":[1,10],"py":[23,5],"pz":[0,-1],"nx":[3,6],"ny":[1,1],"nz":[2,1]},{"size":4,"px":[13,6,3,4],"py":[8,6,4,2],"pz":[0,-1,-1,-1],"nx":[1,1,1,4],"ny":[9,7,8,20],"nz":[1,1,1,0]},{"size":5,"px":[11,4,4,10,3],"py":[9,16,13,12,7],"pz":[0,0,0,0,0],"nx":[7,11,3,17,4],"ny":[8,11,9,0,4],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[6,6],"py":[6,8],"pz":[1,-1],"nx":[0,0],"ny":[1,2],"nz":[2,2]},{"size":2,"px":[10,5],"py":[7,2],"pz":[0,-1],"nx":[4,13],"ny":[5,9],"nz":[2,0]},{"size":2,"px":[10,5],"py":[8,2],"pz":[1,-1],"nx":[16,4],"ny":[14,5],"nz":[0,2]},{"size":2,"px":[1,1],"py":[16,15],"pz":[0,0],"nx":[1,20],"ny":[23,1],"nz":[0,-1]},{"size":2,"px":[2,3],"py":[4,7],"pz":[2,1],"nx":[2,3],"ny":[5,4],"nz":[2,-1]},{"size":2,"px":[19,8],"py":[5,4],"pz":[0,-1],"nx":[10,10],"ny":[1,3],"nz":[1,1]},{"size":2,"px":[21,21],"py":[18,16],"pz":[0,0],"nx":[10,3],"ny":[17,5],"nz":[0,-1]},{"size":2,"px":[9,2],"py":[23,4],"pz":[0,2],"nx":[5,11],"ny":[3,7],"nz":[2,1]},{"size":2,"px":[7,0],"py":[3,2],"pz":[0,-1],"nx":[3,6],"ny":[1,1],"nz":[1,0]},{"size":4,"px":[5,9,8,9],"py":[8,12,13,18],"pz":[0,0,0,0],"nx":[6,5,2,5],"ny":[8,4,7,11],"nz":[0,-1,-1,-1]},{"size":2,"px":[7,2],"py":[0,0],"pz":[0,2],"nx":[5,5],"ny":[3,4],"nz":[1,-1]},{"size":2,"px":[11,11],"py":[12,13],"pz":[0,0],"nx":[9,1],"ny":[14,3],"nz":[0,-1]},{"size":5,"px":[8,16,9,4,15],"py":[11,13,8,4,12],"pz":[1,0,1,2,0],"nx":[3,3,3,3,4],"ny":[4,2,1,3,0],"nz":[0,0,0,0,0]},{"size":2,"px":[9,5],"py":[7,6],"pz":[1,-1],"nx":[19,8],"ny":[17,11],"nz":[0,1]},{"size":5,"px":[14,15,12,13,13],"py":[2,2,2,2,2],"pz":[0,0,0,0,-1],"nx":[20,9,19,20,4],"ny":[14,2,5,15,1],"nz":[0,1,0,0,2]},{"size":2,"px":[18,8],"py":[20,7],"pz":[0,1],"nx":[4,9],"ny":[2,2],"nz":[2,-1]},{"size":2,"px":[6,3],"py":[11,5],"pz":[1,2],"nx":[13,19],"ny":[20,20],"nz":[0,-1]},{"size":3,"px":[12,11,3],"py":[20,20,5],"pz":[0,0,-1],"nx":[11,12,6],"ny":[21,21,10],"nz":[0,0,1]},{"size":2,"px":[3,6],"py":[7,14],"pz":[1,0],"nx":[3,13],"ny":[4,8],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[5,9],"pz":[2,1],"nx":[2,11],"ny":[8,6],"nz":[1,-1]},{"size":2,"px":[2,2],"py":[5,5],"pz":[1,-1],"nx":[0,0],"ny":[6,3],"nz":[1,2]},{"size":2,"px":[11,23],"py":[5,9],"pz":[1,0],"nx":[8,2],"ny":[11,0],"nz":[0,-1]},{"size":2,"px":[11,23],"py":[12,9],"pz":[0,-1],"nx":[11,22],"ny":[10,21],"nz":[1,0]},{"size":2,"px":[12,12],"py":[7,7],"pz":[0,-1],"nx":[5,4],"ny":[7,10],"nz":[1,1]},{"size":2,"px":[9,8],"py":[18,1],"pz":[0,-1],"nx":[5,4],"ny":[8,10],"nz":[1,1]},{"size":2,"px":[16,17],"py":[11,11],"pz":[0,0],"nx":[15,2],"ny":[9,4],"nz":[0,-1]},{"size":2,"px":[0,1],"py":[3,0],"pz":[2,-1],"nx":[9,10],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[13,13],"py":[20,21],"pz":[0,-1],"nx":[2,2],"ny":[6,5],"nz":[1,1]},{"size":5,"px":[20,20,4,18,19],"py":[17,16,5,22,20],"pz":[0,0,2,0,0],"nx":[8,11,5,6,2],"ny":[10,15,11,10,1],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[11,11],"py":[4,4],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":3,"px":[6,5,6],"py":[8,10,10],"pz":[1,1,1],"nx":[11,8,22],"ny":[19,2,15],"nz":[0,-1,-1]},{"size":3,"px":[5,2,13],"py":[7,10,10],"pz":[1,-1,-1],"nx":[11,11,23],"ny":[8,9,14],"nz":[1,1,0]},{"size":5,"px":[3,6,1,5,10],"py":[7,14,1,9,2],"pz":[1,-1,-1,-1,-1],"nx":[11,0,1,5,1],"ny":[14,12,18,5,19],"nz":[0,0,0,1,0]},{"size":3,"px":[21,21,10],"py":[16,17,10],"pz":[0,0,1],"nx":[5,5,1],"ny":[9,9,18],"nz":[1,-1,-1]},{"size":2,"px":[6,21],"py":[6,17],"pz":[1,-1],"nx":[20,10],"ny":[7,4],"nz":[0,1]},{"size":2,"px":[10,11],"py":[0,0],"pz":[1,-1],"nx":[6,13],"ny":[2,4],"nz":[1,0]},{"size":4,"px":[4,4,7,9],"py":[3,4,10,3],"pz":[2,2,1,1],"nx":[21,2,15,5],"ny":[0,0,0,2],"nz":[0,-1,-1,-1]},{"size":3,"px":[11,11,11],"py":[7,6,9],"pz":[1,1,1],"nx":[23,4,9],"ny":[23,5,6],"nz":[0,-1,-1]},{"size":2,"px":[14,15],"py":[1,1],"pz":[0,0],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[11,23,11,23,23],"py":[11,22,10,21,20],"pz":[1,0,1,0,0],"nx":[10,9,19,10,10],"ny":[10,11,20,9,9],"nz":[1,1,0,1,-1]},{"size":2,"px":[7,23],"py":[13,22],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[12,1],"py":[19,0],"pz":[0,-1],"nx":[11,12],"ny":[22,17],"nz":[0,0]},{"size":2,"px":[10,8],"py":[4,3],"pz":[1,-1],"nx":[5,23],"ny":[2,7],"nz":[2,0]},{"size":2,"px":[9,10],"py":[6,20],"pz":[1,-1],"nx":[8,8],"ny":[4,6],"nz":[1,1]}],"alpha":[-1.135386e+00,1.135386e+00,-9.090800e-01,9.090800e-01,-5.913780e-01,5.913780e-01,-5.556534e-01,5.556534e-01,-5.084150e-01,5.084150e-01,-4.464489e-01,4.464489e-01,-4.463241e-01,4.463241e-01,-4.985226e-01,4.985226e-01,-4.424638e-01,4.424638e-01,-4.300093e-01,4.300093e-01,-4.231341e-01,4.231341e-01,-4.087428e-01,4.087428e-01,-3.374480e-01,3.374480e-01,-3.230151e-01,3.230151e-01,-3.084427e-01,3.084427e-01,-3.235494e-01,3.235494e-01,-2.589281e-01,2.589281e-01,-2.970292e-01,2.970292e-01,-2.957065e-01,2.957065e-01,-3.997619e-01,3.997619e-01,-3.535901e-01,3.535901e-01,-2.725396e-01,2.725396e-01,-2.649725e-01,2.649725e-01,-3.103888e-01,3.103888e-01,-3.117775e-01,3.117775e-01,-2.589620e-01,2.589620e-01,-2.689202e-01,2.689202e-01,-2.127024e-01,2.127024e-01,-2.436322e-01,2.436322e-01,-3.120574e-01,3.120574e-01,-2.786010e-01,2.786010e-01,-2.649072e-01,2.649072e-01,-2.766509e-01,2.766509e-01,-2.367237e-01,2.367237e-01,-2.658049e-01,2.658049e-01,-2.103463e-01,2.103463e-01,-1.911522e-01,1.911522e-01,-2.535425e-01,2.535425e-01,-2.434696e-01,2.434696e-01,-2.180788e-01,2.180788e-01,-2.496873e-01,2.496873e-01,-2.700969e-01,2.700969e-01,-2.565479e-01,2.565479e-01,-2.737741e-01,2.737741e-01,-1.675507e-01,1.675507e-01,-2.551417e-01,2.551417e-01,-2.067648e-01,2.067648e-01,-1.636834e-01,1.636834e-01,-2.129306e-01,2.129306e-01,-1.656758e-01,1.656758e-01,-1.919369e-01,1.919369e-01,-2.031763e-01,2.031763e-01,-2.062327e-01,2.062327e-01,-2.577950e-01,2.577950e-01,-2.951823e-01,2.951823e-01,-2.023160e-01,2.023160e-01,-2.022234e-01,2.022234e-01,-2.132906e-01,2.132906e-01,-1.653278e-01,1.653278e-01,-1.648474e-01,1.648474e-01,-1.593352e-01,1.593352e-01,-1.735650e-01,1.735650e-01,-1.688778e-01,1.688778e-01,-1.519705e-01,1.519705e-01,-1.812202e-01,1.812202e-01,-1.967481e-01,1.967481e-01,-1.852954e-01,1.852954e-01,-2.317780e-01,2.317780e-01,-2.036251e-01,2.036251e-01,-1.609324e-01,1.609324e-01,-2.160205e-01,2.160205e-01,-2.026190e-01,2.026190e-01,-1.854761e-01,1.854761e-01,-1.832038e-01,1.832038e-01,-2.001141e-01,2.001141e-01,-1.418333e-01,1.418333e-01,-1.704773e-01,1.704773e-01,-1.586261e-01,1.586261e-01,-1.587582e-01,1.587582e-01,-1.899489e-01,1.899489e-01,-1.477160e-01,1.477160e-01,-2.260467e-01,2.260467e-01,-2.393598e-01,2.393598e-01,-1.582373e-01,1.582373e-01,-1.702498e-01,1.702498e-01,-1.737398e-01,1.737398e-01,-1.462529e-01,1.462529e-01,-1.396517e-01,1.396517e-01,-1.629625e-01,1.629625e-01,-1.446933e-01,1.446933e-01,-1.811657e-01,1.811657e-01,-1.336427e-01,1.336427e-01,-1.924813e-01,1.924813e-01,-1.457520e-01,1.457520e-01,-1.600259e-01,1.600259e-01,-1.297000e-01,1.297000e-01,-2.076199e-01,2.076199e-01,-1.510060e-01,1.510060e-01,-1.914568e-01,1.914568e-01,-2.138162e-01,2.138162e-01,-1.856916e-01,1.856916e-01,-1.843047e-01,1.843047e-01,-1.526846e-01,1.526846e-01,-1.328320e-01,1.328320e-01,-1.751311e-01,1.751311e-01,-1.643908e-01,1.643908e-01,-1.482706e-01,1.482706e-01,-1.622298e-01,1.622298e-01,-1.884979e-01,1.884979e-01,-1.633604e-01,1.633604e-01,-1.554166e-01,1.554166e-01,-1.405332e-01,1.405332e-01,-1.772398e-01,1.772398e-01,-1.410008e-01,1.410008e-01,-1.362301e-01,1.362301e-01,-1.709087e-01,1.709087e-01,-1.584613e-01,1.584613e-01,-1.188814e-01,1.188814e-01,-1.423888e-01,1.423888e-01,-1.345565e-01,1.345565e-01,-1.835986e-01,1.835986e-01,-1.445329e-01,1.445329e-01,-1.385826e-01,1.385826e-01,-1.558917e-01,1.558917e-01,-1.476053e-01,1.476053e-01,-1.370722e-01,1.370722e-01,-2.362666e-01,2.362666e-01,-2.907774e-01,2.907774e-01,-1.656360e-01,1.656360e-01,-1.644407e-01,1.644407e-01,-1.443394e-01,1.443394e-01,-1.438823e-01,1.438823e-01,-1.476964e-01,1.476964e-01,-1.956593e-01,1.956593e-01,-2.417519e-01,2.417519e-01,-1.659315e-01,1.659315e-01,-1.466254e-01,1.466254e-01,-2.034909e-01,2.034909e-01,-2.128771e-01,2.128771e-01,-1.665429e-01,1.665429e-01,-1.387131e-01,1.387131e-01,-1.298823e-01,1.298823e-01,-1.329495e-01,1.329495e-01,-1.769587e-01,1.769587e-01,-1.366530e-01,1.366530e-01,-1.254359e-01,1.254359e-01,-1.673022e-01,1.673022e-01,-1.602519e-01,1.602519e-01,-1.897245e-01,1.897245e-01,-1.893579e-01,1.893579e-01,-1.579350e-01,1.579350e-01,-1.472589e-01,1.472589e-01,-1.614193e-01,1.614193e-01]},{"count":203,"threshold":-4.769677e+00,"feature":[{"size":5,"px":[12,5,14,9,7],"py":[9,13,3,1,3],"pz":[0,0,0,0,0],"nx":[1,0,5,14,9],"ny":[5,3,8,8,9],"nz":[2,0,1,0,0]},{"size":5,"px":[14,13,11,17,12],"py":[2,2,4,13,3],"pz":[0,0,0,0,0],"nx":[7,22,8,23,22],"ny":[8,15,11,12,3],"nz":[1,0,1,0,0]},{"size":5,"px":[9,11,11,11,16],"py":[4,8,7,9,12],"pz":[0,0,0,0,0],"nx":[4,8,14,9,9],"ny":[4,4,8,8,8],"nz":[1,1,0,0,-1]},{"size":5,"px":[6,12,12,8,3],"py":[11,7,8,10,2],"pz":[0,0,0,0,2],"nx":[8,4,4,4,0],"ny":[4,4,4,11,0],"nz":[1,1,-1,-1,-1]},{"size":5,"px":[19,17,18,9,9],"py":[3,2,3,1,1],"pz":[0,0,0,1,-1],"nx":[21,21,10,22,22],"ny":[1,2,0,4,3],"nz":[0,0,1,0,0]},{"size":2,"px":[4,7],"py":[4,6],"pz":[2,1],"nx":[8,7],"ny":[4,10],"nz":[1,1]},{"size":5,"px":[14,17,17,13,12],"py":[18,15,16,18,18],"pz":[0,0,0,0,0],"nx":[13,19,5,20,6],"ny":[16,4,1,19,0],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[6,7,4,5,5],"py":[15,23,6,12,16],"pz":[0,0,1,0,0],"nx":[3,14,14,6,6],"ny":[4,11,11,9,0],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[16,9,6,3,11],"py":[2,2,5,3,2],"pz":[0,0,1,2,0],"nx":[3,4,2,5,5],"ny":[4,11,2,8,8],"nz":[1,1,2,1,-1]},{"size":5,"px":[6,1,5,3,3],"py":[14,4,15,7,7],"pz":[0,2,0,1,-1],"nx":[0,0,1,1,1],"ny":[7,8,18,17,5],"nz":[1,1,0,0,2]},{"size":5,"px":[12,12,9,5,3],"py":[14,14,0,3,7],"pz":[0,-1,-1,-1,-1],"nx":[7,7,14,8,13],"ny":[7,8,13,10,10],"nz":[1,1,0,1,0]},{"size":2,"px":[3,4],"py":[7,9],"pz":[1,-1],"nx":[2,4],"ny":[5,4],"nz":[2,1]},{"size":3,"px":[10,21,17],"py":[7,11,23],"pz":[1,0,0],"nx":[21,9,3],"ny":[23,5,5],"nz":[0,-1,-1]},{"size":5,"px":[8,11,9,10,11],"py":[2,0,1,1,2],"pz":[0,0,0,0,0],"nx":[4,5,6,4,3],"ny":[8,4,18,7,4],"nz":[1,1,0,1,-1]},{"size":5,"px":[20,22,3,19,10],"py":[20,9,4,22,3],"pz":[0,0,2,0,1],"nx":[8,20,8,3,2],"ny":[4,3,6,4,3],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[4,4],"py":[8,7],"pz":[1,1],"nx":[9,2],"ny":[15,5],"nz":[0,-1]},{"size":2,"px":[11,13],"py":[13,4],"pz":[0,-1],"nx":[20,21],"ny":[1,4],"nz":[0,0]},{"size":5,"px":[1,2,7,6,8],"py":[0,2,3,3,3],"pz":[2,1,0,0,0],"nx":[1,2,1,1,1],"ny":[0,0,4,3,3],"nz":[1,0,0,0,-1]},{"size":2,"px":[3,10],"py":[9,11],"pz":[0,0],"nx":[6,3],"ny":[9,2],"nz":[0,-1]},{"size":5,"px":[12,12,12,12,6],"py":[10,11,13,12,6],"pz":[0,0,0,0,-1],"nx":[10,2,1,10,10],"ny":[10,4,2,11,9],"nz":[0,1,2,0,0]},{"size":5,"px":[16,18,11,17,15],"py":[11,12,8,12,11],"pz":[0,0,0,0,0],"nx":[14,0,19,0,10],"ny":[9,3,14,8,9],"nz":[0,-1,-1,-1,-1]},{"size":4,"px":[5,9,5,8],"py":[21,18,20,23],"pz":[0,0,0,0],"nx":[8,4,3,1],"ny":[20,3,4,3],"nz":[0,-1,-1,-1]},{"size":2,"px":[2,3],"py":[3,2],"pz":[2,2],"nx":[3,12],"ny":[4,23],"nz":[1,-1]},{"size":5,"px":[0,1,1,1,1],"py":[2,16,14,13,12],"pz":[2,0,0,0,0],"nx":[8,4,9,4,7],"ny":[9,3,4,2,9],"nz":[1,2,1,2,1]},{"size":2,"px":[4,9],"py":[3,7],"pz":[2,-1],"nx":[4,9],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[15,16,17,15,8],"py":[3,3,3,18,1],"pz":[0,0,0,0,1],"nx":[1,2,2,1,3],"ny":[5,3,2,6,0],"nz":[0,0,0,0,0]},{"size":2,"px":[4,17],"py":[4,14],"pz":[2,0],"nx":[15,7],"ny":[15,10],"nz":[0,-1]},{"size":3,"px":[14,12,3],"py":[3,13,3],"pz":[0,-1,-1],"nx":[4,17,4],"ny":[3,19,4],"nz":[2,0,2]},{"size":4,"px":[4,5,12,2],"py":[9,6,19,4],"pz":[1,1,0,2],"nx":[12,17,4,4],"ny":[18,19,4,4],"nz":[0,-1,-1,-1]},{"size":5,"px":[10,19,20,20,19],"py":[7,14,13,14,13],"pz":[1,0,0,0,-1],"nx":[11,23,23,23,23],"ny":[9,15,13,16,14],"nz":[1,0,0,0,0]},{"size":4,"px":[0,0,0,2],"py":[5,6,5,14],"pz":[1,1,2,0],"nx":[0,3,3,17],"ny":[23,5,5,9],"nz":[0,-1,-1,-1]},{"size":2,"px":[15,4],"py":[23,5],"pz":[0,2],"nx":[9,3],"ny":[4,4],"nz":[1,-1]},{"size":4,"px":[6,5,10,12],"py":[3,3,23,23],"pz":[1,1,0,0],"nx":[11,1,1,4],"ny":[21,3,5,5],"nz":[0,-1,-1,-1]},{"size":2,"px":[5,2],"py":[9,4],"pz":[1,2],"nx":[4,9],"ny":[4,2],"nz":[1,-1]},{"size":5,"px":[23,23,23,23,23],"py":[14,9,13,11,12],"pz":[0,0,0,0,0],"nx":[6,13,7,8,8],"ny":[9,6,3,3,3],"nz":[1,0,1,1,-1]},{"size":2,"px":[10,3],"py":[4,5],"pz":[0,-1],"nx":[3,8],"ny":[1,3],"nz":[2,1]},{"size":2,"px":[3,12],"py":[4,18],"pz":[2,0],"nx":[12,0],"ny":[16,3],"nz":[0,-1]},{"size":2,"px":[16,2],"py":[4,4],"pz":[0,-1],"nx":[16,4],"ny":[1,0],"nz":[0,2]},{"size":2,"px":[3,4],"py":[7,1],"pz":[1,-1],"nx":[5,3],"ny":[19,9],"nz":[0,1]},{"size":4,"px":[20,19,20,21],"py":[2,0,1,3],"pz":[0,0,0,0],"nx":[11,5,23,11],"ny":[0,0,1,1],"nz":[1,2,0,1]},{"size":2,"px":[12,13],"py":[7,5],"pz":[0,0],"nx":[8,5],"ny":[3,5],"nz":[1,-1]},{"size":5,"px":[22,21,22,22,22],"py":[20,22,18,19,16],"pz":[0,0,0,0,0],"nx":[2,3,3,15,15],"ny":[4,5,4,7,7],"nz":[1,2,1,0,-1]},{"size":3,"px":[15,14,14],"py":[1,1,1],"pz":[0,0,-1],"nx":[17,18,16],"ny":[1,2,1],"nz":[0,0,0]},{"size":4,"px":[17,16,16,15],"py":[2,1,0,0],"pz":[0,0,0,0],"nx":[7,4,2,11],"ny":[11,2,1,4],"nz":[1,2,-1,-1]},{"size":4,"px":[18,0,0,0],"py":[14,6,5,4],"pz":[0,-1,-1,-1],"nx":[19,19,19,19],"ny":[16,19,17,18],"nz":[0,0,0,0]},{"size":4,"px":[11,5,5,0],"py":[14,1,4,4],"pz":[0,-1,-1,-1],"nx":[11,8,2,15],"ny":[17,14,1,9],"nz":[0,0,2,0]},{"size":2,"px":[4,5],"py":[19,21],"pz":[0,0],"nx":[10,2],"ny":[15,4],"nz":[0,-1]},{"size":2,"px":[6,4],"py":[4,6],"pz":[1,1],"nx":[3,3],"ny":[4,5],"nz":[1,-1]},{"size":2,"px":[2,7],"py":[1,13],"pz":[2,0],"nx":[7,2],"ny":[1,4],"nz":[1,-1]},{"size":4,"px":[15,10,4,7],"py":[23,3,1,7],"pz":[0,1,2,1],"nx":[0,4,1,1],"ny":[0,2,0,-1900147915],"nz":[0,-1,-1,-1]},{"size":2,"px":[7,2],"py":[12,11],"pz":[0,-1],"nx":[2,4],"ny":[2,5],"nz":[2,1]},{"size":5,"px":[0,0,0,1,0],"py":[9,4,3,2,6],"pz":[0,1,2,1,1],"nx":[9,4,2,16,16],"ny":[7,4,2,8,8],"nz":[0,1,2,0,-1]},{"size":5,"px":[18,4,9,4,4],"py":[12,5,6,3,4],"pz":[0,2,1,2,-1],"nx":[4,3,3,2,3],"ny":[23,19,21,16,18],"nz":[0,0,0,0,0]},{"size":2,"px":[6,6],"py":[14,13],"pz":[0,0],"nx":[3,10],"ny":[4,7],"nz":[1,-1]},{"size":5,"px":[3,4,4,2,2],"py":[8,11,7,4,4],"pz":[1,1,1,2,-1],"nx":[20,18,19,20,19],"ny":[4,0,2,3,1],"nz":[0,0,0,0,0]},{"size":5,"px":[17,12,14,8,16],"py":[2,0,0,0,0],"pz":[0,0,0,1,0],"nx":[3,15,3,2,2],"ny":[2,9,7,2,2],"nz":[2,0,1,2,-1]},{"size":5,"px":[11,10,11,11,11],"py":[10,12,11,12,12],"pz":[0,0,0,0,-1],"nx":[13,13,20,10,13],"ny":[9,11,8,4,10],"nz":[0,0,0,1,0]},{"size":2,"px":[8,16],"py":[7,13],"pz":[1,0],"nx":[8,13],"ny":[4,11],"nz":[1,-1]},{"size":2,"px":[6,7],"py":[20,3],"pz":[0,-1],"nx":[3,4],"ny":[10,10],"nz":[1,1]},{"size":3,"px":[13,10,17],"py":[9,3,5],"pz":[0,-1,-1],"nx":[1,3,1],"ny":[5,16,6],"nz":[2,0,1]},{"size":2,"px":[0,0],"py":[5,5],"pz":[2,-1],"nx":[8,3],"ny":[14,10],"nz":[0,1]},{"size":4,"px":[11,9,12,10],"py":[2,2,2,2],"pz":[0,0,0,0],"nx":[4,4,4,10],"ny":[5,5,0,16],"nz":[1,-1,-1,-1]},{"size":3,"px":[7,9,12],"py":[2,2,2],"pz":[1,-1,-1],"nx":[4,7,2],"ny":[3,1,0],"nz":[0,0,2]},{"size":2,"px":[2,4],"py":[3,12],"pz":[2,0],"nx":[7,4],"ny":[6,5],"nz":[1,2]},{"size":4,"px":[12,12,6,3],"py":[12,11,21,7],"pz":[0,0,-1,-1],"nx":[1,0,0,0],"ny":[13,3,6,5],"nz":[0,2,1,1]},{"size":3,"px":[3,1,3],"py":[21,8,18],"pz":[0,1,0],"nx":[11,20,0],"ny":[17,17,6],"nz":[0,-1,-1]},{"size":2,"px":[2,8],"py":[3,12],"pz":[2,0],"nx":[2,20],"ny":[4,17],"nz":[1,-1]},{"size":5,"px":[2,3,4,3,2],"py":[10,14,14,15,13],"pz":[1,0,0,0,0],"nx":[0,0,1,0,0],"ny":[21,20,23,19,19],"nz":[0,0,0,0,-1]},{"size":2,"px":[2,15],"py":[7,4],"pz":[1,-1],"nx":[3,8],"ny":[4,14],"nz":[1,0]},{"size":5,"px":[19,14,12,15,4],"py":[8,12,10,16,2],"pz":[0,0,0,0,2],"nx":[8,0,12,4,0],"ny":[4,1,12,2,19],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[18,9],"py":[15,3],"pz":[0,-1],"nx":[8,15],"ny":[9,14],"nz":[1,0]},{"size":5,"px":[4,2,3,4,9],"py":[9,4,3,8,23],"pz":[1,2,1,1,0],"nx":[11,23,23,11,11],"ny":[0,2,3,1,1],"nz":[1,0,0,1,-1]},{"size":2,"px":[6,7],"py":[1,1],"pz":[0,0],"nx":[3,4],"ny":[10,5],"nz":[1,-1]},{"size":4,"px":[11,9,8,5],"py":[12,15,13,3],"pz":[0,-1,-1,-1],"nx":[3,12,14,13],"ny":[0,3,3,3],"nz":[2,0,0,0]},{"size":2,"px":[11,11],"py":[6,5],"pz":[0,0],"nx":[8,11],"ny":[4,20],"nz":[1,-1]},{"size":5,"px":[21,20,21,21,21],"py":[18,21,17,19,19],"pz":[0,0,0,0,-1],"nx":[2,5,4,4,5],"ny":[5,12,11,10,10],"nz":[1,0,0,0,0]},{"size":5,"px":[1,1,1,1,1],"py":[10,11,7,9,8],"pz":[0,0,0,0,0],"nx":[11,23,23,23,23],"ny":[10,20,21,19,19],"nz":[1,0,0,0,-1]},{"size":5,"px":[7,8,7,3,1],"py":[14,13,13,2,2],"pz":[0,0,-1,-1,-1],"nx":[1,10,2,2,10],"ny":[2,13,4,16,12],"nz":[2,0,1,0,0]},{"size":2,"px":[17,18],"py":[12,12],"pz":[0,0],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[17,0],"py":[5,20],"pz":[0,-1],"nx":[4,9],"ny":[0,2],"nz":[2,1]},{"size":5,"px":[22,22,22,11,23],"py":[16,15,14,6,13],"pz":[0,0,0,1,0],"nx":[16,15,7,9,9],"ny":[15,8,4,10,10],"nz":[0,0,1,1,-1]},{"size":2,"px":[13,3],"py":[3,1],"pz":[0,2],"nx":[8,3],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[5,6],"py":[4,1],"pz":[1,-1],"nx":[6,3],"ny":[4,2],"nz":[1,2]},{"size":3,"px":[4,2,6],"py":[6,3,4],"pz":[1,2,1],"nx":[10,0,4],"ny":[9,4,3],"nz":[0,-1,-1]},{"size":4,"px":[2,8,4,10],"py":[4,23,7,23],"pz":[2,0,1,0],"nx":[9,4,11,9],"ny":[21,5,16,0],"nz":[0,-1,-1,-1]},{"size":2,"px":[6,3],"py":[13,0],"pz":[0,-1],"nx":[8,2],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[3,3],"py":[1,4],"pz":[1,-1],"nx":[3,5],"ny":[0,1],"nz":[1,0]},{"size":2,"px":[7,2],"py":[0,0],"pz":[0,2],"nx":[2,10],"ny":[1,6],"nz":[2,0]},{"size":2,"px":[10,2],"py":[7,0],"pz":[1,-1],"nx":[21,5],"ny":[15,4],"nz":[0,2]},{"size":2,"px":[1,1],"py":[10,9],"pz":[0,0],"nx":[0,3],"ny":[13,11],"nz":[0,-1]},{"size":2,"px":[11,9],"py":[13,0],"pz":[0,-1],"nx":[3,3],"ny":[4,3],"nz":[1,1]},{"size":5,"px":[14,13,13,14,14],"py":[12,10,11,13,13],"pz":[0,0,0,0,-1],"nx":[9,8,4,5,7],"ny":[4,4,2,2,4],"nz":[0,0,1,1,0]},{"size":3,"px":[2,4,1],"py":[2,0,0],"pz":[0,0,1],"nx":[0,7,4],"ny":[0,3,2],"nz":[1,-1,-1]},{"size":2,"px":[11,4],"py":[5,0],"pz":[0,-1],"nx":[8,6],"ny":[4,9],"nz":[1,1]},{"size":3,"px":[0,0,0],"py":[20,2,4],"pz":[0,-1,-1],"nx":[12,3,10],"ny":[3,1,3],"nz":[0,2,0]},{"size":5,"px":[5,11,10,13,13],"py":[0,0,0,2,2],"pz":[1,0,0,0,-1],"nx":[4,5,5,4,5],"ny":[14,0,2,6,1],"nz":[0,0,0,0,0]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,11],"ny":[4,1],"nz":[1,-1]},{"size":2,"px":[14,-1715597992],"py":[19,9],"pz":[0,-1],"nx":[7,14],"ny":[10,17],"nz":[1,0]},{"size":2,"px":[11,1],"py":[9,0],"pz":[0,-1],"nx":[1,12],"ny":[2,10],"nz":[2,0]},{"size":2,"px":[17,9],"py":[13,17],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[0,7],"py":[1,9],"pz":[1,-1],"nx":[18,4],"ny":[14,2],"nz":[0,2]},{"size":2,"px":[14,7],"py":[23,9],"pz":[0,-1],"nx":[4,8],"ny":[5,10],"nz":[2,1]},{"size":2,"px":[8,7],"py":[17,9],"pz":[0,-1],"nx":[3,2],"ny":[0,3],"nz":[0,0]},{"size":2,"px":[13,4],"py":[20,1],"pz":[0,-1],"nx":[5,3],"ny":[21,17],"nz":[0,0]},{"size":3,"px":[0,0,1],"py":[3,6,15],"pz":[2,1,0],"nx":[10,8,3],"ny":[6,4,2],"nz":[0,-1,-1]},{"size":2,"px":[8,8],"py":[18,8],"pz":[0,-1],"nx":[5,4],"ny":[8,10],"nz":[1,1]},{"size":2,"px":[6,5],"py":[2,2],"pz":[1,1],"nx":[8,9],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[6,3],"py":[11,5],"pz":[1,2],"nx":[13,3],"ny":[19,2],"nz":[0,-1]},{"size":2,"px":[4,6],"py":[1,11],"pz":[2,-1],"nx":[3,2],"ny":[1,0],"nz":[1,2]},{"size":2,"px":[9,4],"py":[10,5],"pz":[1,2],"nx":[8,4],"ny":[10,4],"nz":[1,-1]},{"size":2,"px":[12,12],"py":[11,20],"pz":[0,-1],"nx":[0,0],"ny":[6,10],"nz":[1,0]},{"size":2,"px":[7,12],"py":[2,20],"pz":[0,-1],"nx":[2,2],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[0,15],"py":[5,21],"pz":[1,-1],"nx":[10,9],"ny":[3,3],"nz":[0,1]},{"size":2,"px":[15,9],"py":[1,0],"pz":[0,1],"nx":[19,3],"ny":[0,3],"nz":[0,-1]},{"size":2,"px":[21,5],"py":[13,5],"pz":[0,2],"nx":[23,6],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[5,8],"py":[3,1],"pz":[2,-1],"nx":[9,9],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[2,2],"py":[7,7],"pz":[1,-1],"nx":[5,3],"ny":[23,17],"nz":[0,0]},{"size":2,"px":[11,3],"py":[6,4],"pz":[0,-1],"nx":[2,4],"ny":[2,4],"nz":[2,1]},{"size":3,"px":[14,0,17],"py":[20,3,21],"pz":[0,-1,-1],"nx":[11,11,11],"ny":[7,9,10],"nz":[1,1,1]},{"size":5,"px":[11,11,23,23,12],"py":[10,11,21,20,12],"pz":[1,1,0,0,0],"nx":[8,3,6,7,7],"ny":[4,5,11,11,11],"nz":[1,2,1,1,-1]},{"size":2,"px":[11,11],"py":[11,10],"pz":[0,0],"nx":[9,3],"ny":[2,5],"nz":[1,-1]},{"size":2,"px":[12,14],"py":[19,19],"pz":[0,0],"nx":[12,13],"ny":[18,17],"nz":[0,-1]},{"size":5,"px":[13,14,12,15,14],"py":[0,0,1,1,1],"pz":[0,0,0,0,0],"nx":[4,8,4,7,7],"ny":[3,4,2,5,5],"nz":[2,1,2,1,-1]},{"size":2,"px":[17,5],"py":[10,2],"pz":[0,-1],"nx":[4,9],"ny":[2,3],"nz":[2,1]},{"size":2,"px":[18,10],"py":[6,10],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[8,18,8,4,16],"py":[6,12,9,4,13],"pz":[1,0,1,2,0],"nx":[3,4,3,5,5],"ny":[0,2,3,1,1],"nz":[1,0,0,0,-1]},{"size":2,"px":[3,6],"py":[2,4],"pz":[2,1],"nx":[8,0],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[4,5],"pz":[2,-1],"nx":[4,2],"ny":[14,7],"nz":[0,1]},{"size":4,"px":[3,4,4,3],"py":[11,12,12,2],"pz":[0,0,-1,-1],"nx":[1,2,1,2],"ny":[11,14,12,16],"nz":[0,0,0,0]},{"size":2,"px":[6,0],"py":[11,0],"pz":[0,-1],"nx":[3,4],"ny":[4,5],"nz":[1,1]},{"size":2,"px":[3,2],"py":[21,11],"pz":[0,1],"nx":[3,2],"ny":[10,0],"nz":[1,-1]},{"size":3,"px":[10,3,13],"py":[2,0,2],"pz":[0,2,0],"nx":[7,16,1],"ny":[10,4,1],"nz":[0,-1,-1]},{"size":2,"px":[6,12],"py":[2,5],"pz":[1,0],"nx":[6,18],"ny":[1,19],"nz":[1,-1]},{"size":2,"px":[3,16],"py":[0,16],"pz":[1,-1],"nx":[11,2],"ny":[5,1],"nz":[0,2]},{"size":2,"px":[11,10],"py":[13,1],"pz":[0,-1],"nx":[1,1],"ny":[22,21],"nz":[0,0]},{"size":2,"px":[11,10],"py":[18,18],"pz":[0,0],"nx":[5,8],"ny":[9,0],"nz":[1,-1]},{"size":2,"px":[3,2],"py":[20,18],"pz":[0,0],"nx":[8,3],"ny":[5,1],"nz":[1,-1]},{"size":2,"px":[14,2],"py":[17,1],"pz":[0,-1],"nx":[14,13],"ny":[15,15],"nz":[0,0]},{"size":2,"px":[3,4],"py":[2,3],"pz":[2,2],"nx":[8,3],"ny":[4,0],"nz":[1,-1]},{"size":5,"px":[8,18,18,8,7],"py":[6,11,11,7,9],"pz":[1,0,-1,-1,-1],"nx":[5,13,5,11,5],"ny":[3,11,0,8,2],"nz":[2,0,2,1,2]},{"size":5,"px":[12,0,5,4,7],"py":[15,0,4,0,9],"pz":[0,-1,-1,-1,-1],"nx":[8,7,4,16,6],"ny":[17,12,9,10,12],"nz":[0,0,1,0,0]},{"size":2,"px":[6,7],"py":[14,1],"pz":[0,-1],"nx":[5,4],"ny":[9,4],"nz":[1,1]},{"size":4,"px":[8,0,22,4],"py":[4,4,23,0],"pz":[0,-1,-1,-1],"nx":[2,4,2,5],"ny":[0,1,2,9],"nz":[2,1,2,1]},{"size":5,"px":[9,9,10,10,8],"py":[0,1,1,2,0],"pz":[1,1,1,1,1],"nx":[4,16,16,16,6],"ny":[2,11,11,11,12],"nz":[2,0,-1,-1,-1]},{"size":2,"px":[6,6],"py":[6,5],"pz":[1,1],"nx":[0,4],"ny":[3,2],"nz":[1,-1]},{"size":3,"px":[10,3,4],"py":[5,9,8],"pz":[1,-1,-1],"nx":[11,23,23],"ny":[7,12,11],"nz":[1,0,0]},{"size":3,"px":[13,12,7],"py":[19,19,10],"pz":[0,0,1],"nx":[13,5,19],"ny":[20,15,22],"nz":[0,-1,-1]},{"size":2,"px":[12,12],"py":[12,13],"pz":[0,0],"nx":[9,10],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[0,12],"py":[1,13],"pz":[2,-1],"nx":[2,7],"ny":[2,13],"nz":[2,0]},{"size":2,"px":[10,10],"py":[8,9],"pz":[1,1],"nx":[19,7],"ny":[23,13],"nz":[0,-1]},{"size":4,"px":[8,7,23,15],"py":[11,12,4,21],"pz":[0,0,-1,-1],"nx":[2,5,1,10],"ny":[6,6,2,13],"nz":[0,1,1,0]},{"size":2,"px":[10,9],"py":[3,3],"pz":[0,0],"nx":[2,3],"ny":[2,4],"nz":[2,-1]},{"size":2,"px":[5,2],"py":[3,4],"pz":[2,-1],"nx":[3,6],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[7,11],"py":[20,16],"pz":[0,-1],"nx":[2,4],"ny":[5,20],"nz":[2,0]},{"size":2,"px":[9,7],"py":[7,5],"pz":[1,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[4,2],"py":[11,3],"pz":[1,2],"nx":[5,5],"ny":[3,5],"nz":[2,-1]},{"size":2,"px":[11,3],"py":[11,5],"pz":[1,-1],"nx":[4,1],"ny":[12,3],"nz":[0,2]},{"size":2,"px":[9,11],"py":[6,4],"pz":[1,-1],"nx":[10,20],"ny":[9,18],"nz":[1,0]},{"size":5,"px":[2,2,2,2,1],"py":[15,13,16,14,7],"pz":[0,0,0,0,1],"nx":[15,8,9,8,4],"ny":[11,6,5,5,4],"nz":[0,1,1,1,-1]},{"size":2,"px":[12,2],"py":[5,5],"pz":[0,-1],"nx":[3,2],"ny":[7,2],"nz":[1,2]},{"size":2,"px":[5,11],"py":[1,3],"pz":[2,1],"nx":[10,10],"ny":[3,3],"nz":[1,-1]},{"size":2,"px":[17,11],"py":[13,18],"pz":[0,-1],"nx":[6,9],"ny":[9,4],"nz":[1,1]},{"size":5,"px":[5,1,2,5,6],"py":[14,4,9,15,23],"pz":[0,2,1,0,0],"nx":[4,9,18,16,17],"ny":[0,1,1,0,0],"nz":[2,1,0,0,0]},{"size":2,"px":[16,17],"py":[0,0],"pz":[0,0],"nx":[23,23],"ny":[5,4],"nz":[0,-1]},{"size":2,"px":[13,8],"py":[20,6],"pz":[0,-1],"nx":[5,6],"ny":[12,10],"nz":[0,1]},{"size":2,"px":[6,15],"py":[15,0],"pz":[0,-1],"nx":[6,3],"ny":[16,4],"nz":[0,1]},{"size":2,"px":[18,20],"py":[7,8],"pz":[0,0],"nx":[18,11],"ny":[9,14],"nz":[0,-1]},{"size":2,"px":[9,4],"py":[12,6],"pz":[0,1],"nx":[3,15],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[5,2],"pz":[1,2],"nx":[5,5],"ny":[2,2],"nz":[1,-1]},{"size":2,"px":[5,20],"py":[1,20],"pz":[1,-1],"nx":[15,17],"ny":[1,2],"nz":[0,0]},{"size":2,"px":[7,2],"py":[16,4],"pz":[0,2],"nx":[4,0],"ny":[10,6],"nz":[1,-1]},{"size":2,"px":[3,8],"py":[5,0],"pz":[1,-1],"nx":[1,1],"ny":[10,18],"nz":[1,0]},{"size":2,"px":[22,0],"py":[3,0],"pz":[0,-1],"nx":[23,11],"ny":[4,1],"nz":[0,1]},{"size":3,"px":[19,10,20],"py":[21,8,18],"pz":[0,1,0],"nx":[3,6,20],"ny":[5,11,14],"nz":[2,-1,-1]},{"size":4,"px":[2,1,6,5],"py":[7,4,23,22],"pz":[1,2,0,0],"nx":[9,19,20,4],"ny":[8,11,9,2],"nz":[0,-1,-1,-1]},{"size":2,"px":[3,6],"py":[2,11],"pz":[2,1],"nx":[12,10],"ny":[21,9],"nz":[0,-1]},{"size":4,"px":[6,0,2,2],"py":[6,1,4,1],"pz":[1,-1,-1,-1],"nx":[0,0,0,0],"ny":[5,8,9,4],"nz":[1,0,0,1]},{"size":5,"px":[3,13,6,11,9],"py":[0,3,1,1,2],"pz":[2,0,1,0,0],"nx":[7,20,16,4,7],"ny":[7,2,19,2,6],"nz":[1,0,0,2,1]},{"size":4,"px":[7,5,2,6],"py":[7,7,4,11],"pz":[0,0,2,1],"nx":[7,1,21,0],"ny":[8,4,11,3],"nz":[0,-1,-1,-1]},{"size":2,"px":[2,2],"py":[3,2],"pz":[2,2],"nx":[8,9],"ny":[3,11],"nz":[1,-1]},{"size":2,"px":[7,13],"py":[3,5],"pz":[1,0],"nx":[4,3],"ny":[2,2],"nz":[1,-1]},{"size":4,"px":[3,12,13,11],"py":[0,1,1,1],"pz":[2,0,0,0],"nx":[8,9,13,0],"ny":[4,1,16,3],"nz":[1,-1,-1,-1]},{"size":2,"px":[10,1],"py":[4,14],"pz":[0,-1],"nx":[5,10],"ny":[1,2],"nz":[1,0]},{"size":2,"px":[11,12],"py":[21,21],"pz":[0,0],"nx":[10,11],"ny":[19,19],"nz":[0,0]},{"size":2,"px":[8,12],"py":[6,21],"pz":[1,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[11,7],"py":[19,0],"pz":[0,-1],"nx":[6,5],"ny":[9,11],"nz":[1,1]},{"size":5,"px":[11,11,11,10,10],"py":[10,12,11,13,13],"pz":[0,0,0,0,-1],"nx":[7,13,6,12,7],"ny":[10,6,3,6,11],"nz":[0,0,1,0,0]},{"size":2,"px":[12,11],"py":[6,12],"pz":[0,-1],"nx":[4,8],"ny":[4,4],"nz":[1,1]},{"size":5,"px":[16,15,16,15,17],"py":[1,0,0,1,1],"pz":[0,0,0,0,0],"nx":[13,7,6,12,12],"ny":[5,4,3,6,6],"nz":[0,1,1,0,-1]},{"size":2,"px":[2,3],"py":[1,3],"pz":[2,1],"nx":[1,5],"ny":[1,3],"nz":[2,-1]},{"size":2,"px":[6,3],"py":[13,6],"pz":[0,1],"nx":[4,9],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[0,3],"py":[4,3],"pz":[1,-1],"nx":[4,8],"ny":[3,6],"nz":[2,1]},{"size":2,"px":[6,3],"py":[2,1],"pz":[0,1],"nx":[5,5],"ny":[7,21],"nz":[1,-1]},{"size":2,"px":[8,4],"py":[0,0],"pz":[1,-1],"nx":[19,17],"ny":[1,0],"nz":[0,0]},{"size":4,"px":[8,11,5,0],"py":[6,1,1,22],"pz":[1,-1,-1,-1],"nx":[0,10,10,1],"ny":[6,12,13,4],"nz":[1,0,0,1]},{"size":2,"px":[8,17],"py":[6,13],"pz":[1,0],"nx":[14,17],"ny":[9,3],"nz":[0,-1]},{"size":2,"px":[5,8],"py":[0,4],"pz":[2,-1],"nx":[9,8],"ny":[1,1],"nz":[0,0]},{"size":2,"px":[11,14],"py":[13,9],"pz":[0,-1],"nx":[23,23],"ny":[21,19],"nz":[0,0]},{"size":2,"px":[10,9],"py":[9,3],"pz":[0,-1],"nx":[6,3],"ny":[2,1],"nz":[1,2]},{"size":2,"px":[11,1],"py":[4,4],"pz":[0,-1],"nx":[2,4],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[5,9],"py":[3,3],"pz":[2,-1],"nx":[17,9],"ny":[12,5],"nz":[0,1]},{"size":2,"px":[9,7],"py":[18,16],"pz":[0,-1],"nx":[5,2],"ny":[9,5],"nz":[1,2]},{"size":2,"px":[3,6],"py":[0,1],"pz":[1,-1],"nx":[4,5],"ny":[1,0],"nz":[0,0]}],"alpha":[-1.149973e+00,1.149973e+00,-6.844773e-01,6.844773e-01,-6.635048e-01,6.635048e-01,-4.888349e-01,4.888349e-01,-4.267976e-01,4.267976e-01,-4.258100e-01,4.258100e-01,-4.815853e-01,4.815853e-01,-4.091859e-01,4.091859e-01,-3.137414e-01,3.137414e-01,-3.339860e-01,3.339860e-01,-3.891196e-01,3.891196e-01,-4.167691e-01,4.167691e-01,-3.186609e-01,3.186609e-01,-2.957171e-01,2.957171e-01,-3.210062e-01,3.210062e-01,-2.725684e-01,2.725684e-01,-2.452176e-01,2.452176e-01,-2.812662e-01,2.812662e-01,-3.029622e-01,3.029622e-01,-3.293745e-01,3.293745e-01,-3.441536e-01,3.441536e-01,-2.946918e-01,2.946918e-01,-2.890545e-01,2.890545e-01,-1.949205e-01,1.949205e-01,-2.176102e-01,2.176102e-01,-2.595190e-01,2.595190e-01,-2.690931e-01,2.690931e-01,-2.130294e-01,2.130294e-01,-2.316308e-01,2.316308e-01,-2.798562e-01,2.798562e-01,-2.146988e-01,2.146988e-01,-2.332089e-01,2.332089e-01,-2.470614e-01,2.470614e-01,-2.204300e-01,2.204300e-01,-2.272045e-01,2.272045e-01,-2.583686e-01,2.583686e-01,-2.072299e-01,2.072299e-01,-1.834971e-01,1.834971e-01,-2.332656e-01,2.332656e-01,-3.271297e-01,3.271297e-01,-2.401937e-01,2.401937e-01,-2.006316e-01,2.006316e-01,-2.401947e-01,2.401947e-01,-2.475346e-01,2.475346e-01,-2.579532e-01,2.579532e-01,-2.466235e-01,2.466235e-01,-1.787582e-01,1.787582e-01,-2.036892e-01,2.036892e-01,-1.665028e-01,1.665028e-01,-1.576510e-01,1.576510e-01,-2.036997e-01,2.036997e-01,-2.040734e-01,2.040734e-01,-1.792532e-01,1.792532e-01,-2.174767e-01,2.174767e-01,-1.876948e-01,1.876948e-01,-1.883137e-01,1.883137e-01,-1.923872e-01,1.923872e-01,-2.620218e-01,2.620218e-01,-1.659873e-01,1.659873e-01,-1.475948e-01,1.475948e-01,-1.731607e-01,1.731607e-01,-2.059256e-01,2.059256e-01,-1.586309e-01,1.586309e-01,-1.607668e-01,1.607668e-01,-1.975101e-01,1.975101e-01,-2.130745e-01,2.130745e-01,-1.898872e-01,1.898872e-01,-2.052598e-01,2.052598e-01,-1.599397e-01,1.599397e-01,-1.770134e-01,1.770134e-01,-1.888249e-01,1.888249e-01,-1.515406e-01,1.515406e-01,-1.907771e-01,1.907771e-01,-1.698406e-01,1.698406e-01,-2.079535e-01,2.079535e-01,-1.966967e-01,1.966967e-01,-1.631391e-01,1.631391e-01,-2.158666e-01,2.158666e-01,-2.891774e-01,2.891774e-01,-1.581556e-01,1.581556e-01,-1.475359e-01,1.475359e-01,-1.806169e-01,1.806169e-01,-1.782238e-01,1.782238e-01,-1.660440e-01,1.660440e-01,-1.576919e-01,1.576919e-01,-1.741775e-01,1.741775e-01,-1.427265e-01,1.427265e-01,-1.695880e-01,1.695880e-01,-1.486712e-01,1.486712e-01,-1.533565e-01,1.533565e-01,-1.601464e-01,1.601464e-01,-1.978414e-01,1.978414e-01,-1.746566e-01,1.746566e-01,-1.794736e-01,1.794736e-01,-1.896567e-01,1.896567e-01,-1.666197e-01,1.666197e-01,-1.969351e-01,1.969351e-01,-2.321735e-01,2.321735e-01,-1.592485e-01,1.592485e-01,-1.671464e-01,1.671464e-01,-1.688885e-01,1.688885e-01,-1.868042e-01,1.868042e-01,-1.301138e-01,1.301138e-01,-1.330094e-01,1.330094e-01,-1.268423e-01,1.268423e-01,-1.820868e-01,1.820868e-01,-1.881020e-01,1.881020e-01,-1.580814e-01,1.580814e-01,-1.302653e-01,1.302653e-01,-1.787262e-01,1.787262e-01,-1.658453e-01,1.658453e-01,-1.240772e-01,1.240772e-01,-1.315621e-01,1.315621e-01,-1.756341e-01,1.756341e-01,-1.429438e-01,1.429438e-01,-1.351775e-01,1.351775e-01,-2.035692e-01,2.035692e-01,-1.267670e-01,1.267670e-01,-1.288470e-01,1.288470e-01,-1.393648e-01,1.393648e-01,-1.755962e-01,1.755962e-01,-1.308445e-01,1.308445e-01,-1.703894e-01,1.703894e-01,-1.461334e-01,1.461334e-01,-1.368683e-01,1.368683e-01,-1.244085e-01,1.244085e-01,-1.718163e-01,1.718163e-01,-1.415624e-01,1.415624e-01,-1.752024e-01,1.752024e-01,-1.666463e-01,1.666463e-01,-1.407325e-01,1.407325e-01,-1.258317e-01,1.258317e-01,-1.416511e-01,1.416511e-01,-1.420816e-01,1.420816e-01,-1.562547e-01,1.562547e-01,-1.542952e-01,1.542952e-01,-1.158829e-01,1.158829e-01,-1.392875e-01,1.392875e-01,-1.610095e-01,1.610095e-01,-1.546440e-01,1.546440e-01,-1.416235e-01,1.416235e-01,-2.028817e-01,2.028817e-01,-1.106779e-01,1.106779e-01,-9.231660e-02,9.231660e-02,-1.164460e-01,1.164460e-01,-1.701578e-01,1.701578e-01,-1.277995e-01,1.277995e-01,-1.946177e-01,1.946177e-01,-1.394509e-01,1.394509e-01,-1.370145e-01,1.370145e-01,-1.446031e-01,1.446031e-01,-1.665215e-01,1.665215e-01,-1.435822e-01,1.435822e-01,-1.559354e-01,1.559354e-01,-1.591860e-01,1.591860e-01,-1.193338e-01,1.193338e-01,-1.236954e-01,1.236954e-01,-1.209139e-01,1.209139e-01,-1.267385e-01,1.267385e-01,-1.232397e-01,1.232397e-01,-1.299632e-01,1.299632e-01,-1.302020e-01,1.302020e-01,-1.202975e-01,1.202975e-01,-1.525378e-01,1.525378e-01,-1.123073e-01,1.123073e-01,-1.605678e-01,1.605678e-01,-1.406867e-01,1.406867e-01,-1.354273e-01,1.354273e-01,-1.393192e-01,1.393192e-01,-1.278263e-01,1.278263e-01,-1.172073e-01,1.172073e-01,-1.153493e-01,1.153493e-01,-1.356318e-01,1.356318e-01,-1.316614e-01,1.316614e-01,-1.374489e-01,1.374489e-01,-1.018254e-01,1.018254e-01,-1.473336e-01,1.473336e-01,-1.289687e-01,1.289687e-01,-1.299183e-01,1.299183e-01,-1.178391e-01,1.178391e-01,-1.619059e-01,1.619059e-01,-1.842569e-01,1.842569e-01,-1.829095e-01,1.829095e-01,-1.939918e-01,1.939918e-01,-1.395362e-01,1.395362e-01,-1.774673e-01,1.774673e-01,-1.688216e-01,1.688216e-01,-1.671747e-01,1.671747e-01,-1.850178e-01,1.850178e-01,-1.106695e-01,1.106695e-01,-1.258323e-01,1.258323e-01,-1.246819e-01,1.246819e-01,-9.892193e-02,9.892193e-02,-1.399638e-01,1.399638e-01,-1.228375e-01,1.228375e-01,-1.756236e-01,1.756236e-01,-1.360307e-01,1.360307e-01,-1.266574e-01,1.266574e-01,-1.372135e-01,1.372135e-01,-1.175947e-01,1.175947e-01,-1.330075e-01,1.330075e-01,-1.396152e-01,1.396152e-01,-2.088443e-01,2.088443e-01]},{"count":301,"threshold":-4.887516e+00,"feature":[{"size":5,"px":[8,11,8,14,10],"py":[6,9,3,3,4],"pz":[1,0,0,0,0],"nx":[8,7,19,7,13],"ny":[11,8,8,5,8],"nz":[1,1,0,1,0]},{"size":5,"px":[14,3,13,12,12],"py":[4,6,4,4,8],"pz":[0,1,0,0,0],"nx":[2,5,2,10,10],"ny":[2,8,5,8,8],"nz":[2,1,2,0,-1]},{"size":5,"px":[6,5,3,7,7],"py":[2,3,1,2,2],"pz":[0,0,1,0,-1],"nx":[2,2,1,2,1],"ny":[3,1,2,2,2],"nz":[0,0,2,0,1]},{"size":5,"px":[3,3,6,12,8],"py":[4,2,4,10,17],"pz":[2,2,1,0,0],"nx":[4,8,8,2,1],"ny":[4,4,4,2,2],"nz":[1,1,-1,-1,-1]},{"size":5,"px":[18,19,17,9,16],"py":[1,2,2,0,2],"pz":[0,0,0,1,0],"nx":[23,23,22,22,22],"ny":[4,3,1,0,2],"nz":[0,0,0,0,0]},{"size":3,"px":[15,4,14],"py":[23,4,18],"pz":[0,2,0],"nx":[7,0,5],"ny":[10,4,9],"nz":[1,-1,-1]},{"size":5,"px":[11,11,16,11,17],"py":[8,6,11,7,11],"pz":[0,0,0,0,0],"nx":[8,4,14,14,1],"ny":[4,4,8,8,5],"nz":[1,1,0,-1,-1]},{"size":5,"px":[12,12,12,12,12],"py":[13,10,11,12,12],"pz":[0,0,0,0,-1],"nx":[4,4,1,2,9],"ny":[8,10,2,4,15],"nz":[0,1,2,1,0]},{"size":2,"px":[19,0],"py":[14,17],"pz":[0,-1],"nx":[20,19],"ny":[15,22],"nz":[0,0]},{"size":5,"px":[3,3,1,3,5],"py":[13,15,6,14,22],"pz":[0,0,1,0,0],"nx":[0,0,1,0,0],"ny":[11,21,23,5,5],"nz":[1,0,0,2,-1]},{"size":5,"px":[4,2,10,4,3],"py":[19,4,13,16,13],"pz":[0,1,0,0,0],"nx":[3,20,7,4,0],"ny":[4,19,5,1,5],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[11,5],"py":[4,4],"pz":[0,-1],"nx":[15,3],"ny":[15,1],"nz":[0,2]},{"size":4,"px":[17,17,12,11],"py":[14,15,18,18],"pz":[0,0,0,0],"nx":[11,4,1,0],"ny":[17,20,8,5],"nz":[0,-1,-1,-1]},{"size":5,"px":[6,2,1,2,11],"py":[14,4,1,1,18],"pz":[0,-1,-1,-1,-1],"nx":[5,5,3,5,2],"ny":[18,17,7,9,2],"nz":[0,0,1,1,2]},{"size":5,"px":[20,19,20,15,20],"py":[17,20,12,12,8],"pz":[0,0,0,0,0],"nx":[17,0,5,2,2],"ny":[8,4,9,2,2],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[6,8],"py":[7,11],"pz":[1,-1],"nx":[7,8],"ny":[7,10],"nz":[1,1]},{"size":5,"px":[15,16,14,8,8],"py":[2,2,2,0,0],"pz":[0,0,0,1,-1],"nx":[20,11,21,18,19],"ny":[3,6,5,1,2],"nz":[0,1,0,0,0]},{"size":4,"px":[17,18,9,8],"py":[23,21,7,8],"pz":[0,0,1,1],"nx":[8,17,10,18],"ny":[4,12,2,1],"nz":[1,-1,-1,-1]},{"size":5,"px":[2,2,9,4,8],"py":[7,3,12,12,23],"pz":[1,1,0,0,0],"nx":[0,0,0,0,0],"ny":[3,1,2,4,4],"nz":[0,0,0,0,-1]},{"size":3,"px":[7,8,5],"py":[22,23,9],"pz":[0,0,1],"nx":[9,4,2],"ny":[21,4,0],"nz":[0,-1,-1]},{"size":2,"px":[3,3],"py":[7,7],"pz":[1,-1],"nx":[3,2],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[15,11,10,3,17],"py":[0,1,2,3,1],"pz":[0,0,0,2,0],"nx":[5,8,4,3,3],"ny":[9,4,7,10,10],"nz":[1,1,1,1,-1]},{"size":3,"px":[22,11,22],"py":[12,5,14],"pz":[0,1,0],"nx":[23,23,3],"ny":[22,23,8],"nz":[0,0,-1]},{"size":2,"px":[3,11],"py":[7,5],"pz":[1,-1],"nx":[8,2],"ny":[14,5],"nz":[0,2]},{"size":4,"px":[17,16,2,4],"py":[14,13,5,0],"pz":[0,0,-1,-1],"nx":[8,9,15,8],"ny":[8,9,14,7],"nz":[1,1,0,1]},{"size":2,"px":[5,16],"py":[6,13],"pz":[1,-1],"nx":[2,1],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[1,0,1,2,1],"py":[15,2,16,19,12],"pz":[0,2,0,0,0],"nx":[8,7,4,9,9],"ny":[5,11,4,5,5],"nz":[1,1,1,1,-1]},{"size":2,"px":[8,7],"py":[11,12],"pz":[0,0],"nx":[9,1],"ny":[10,16],"nz":[0,-1]},{"size":2,"px":[15,13],"py":[17,10],"pz":[0,-1],"nx":[7,4],"ny":[8,4],"nz":[1,2]},{"size":5,"px":[11,10,7,8,9],"py":[0,0,1,1,1],"pz":[0,0,0,0,0],"nx":[4,5,4,5,6],"ny":[1,0,2,1,0],"nz":[0,0,0,0,-1]},{"size":2,"px":[2,2],"py":[4,3],"pz":[2,2],"nx":[3,21],"ny":[4,20],"nz":[1,-1]},{"size":5,"px":[10,11,5,2,11],"py":[12,10,6,11,11],"pz":[0,0,1,0,0],"nx":[4,15,16,7,7],"ny":[5,10,11,10,10],"nz":[1,0,0,0,-1]},{"size":5,"px":[13,14,1,11,11],"py":[2,2,3,2,2],"pz":[0,0,2,0,-1],"nx":[3,0,0,1,0],"ny":[23,15,14,9,8],"nz":[0,0,0,1,1]},{"size":2,"px":[17,2],"py":[13,5],"pz":[0,-1],"nx":[4,9],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[10,5],"py":[4,1],"pz":[0,-1],"nx":[11,3],"ny":[3,0],"nz":[0,2]},{"size":2,"px":[5,3],"py":[3,3],"pz":[2,-1],"nx":[11,23],"ny":[8,14],"nz":[1,0]},{"size":3,"px":[22,22,22],"py":[16,18,9],"pz":[0,0,0],"nx":[13,2,0],"ny":[17,3,5],"nz":[0,-1,-1]},{"size":5,"px":[13,10,13,14,11],"py":[2,2,1,2,1],"pz":[0,0,0,0,0],"nx":[3,3,8,6,6],"ny":[2,5,4,11,11],"nz":[2,2,1,1,-1]},{"size":3,"px":[12,1,1],"py":[14,0,1],"pz":[0,-1,-1],"nx":[8,15,7],"ny":[1,2,0],"nz":[1,0,1]},{"size":2,"px":[4,5],"py":[20,23],"pz":[0,0],"nx":[3,3],"ny":[10,2],"nz":[1,-1]},{"size":2,"px":[2,4],"py":[7,2],"pz":[1,-1],"nx":[4,3],"ny":[23,16],"nz":[0,0]},{"size":3,"px":[3,3,6],"py":[5,2,4],"pz":[2,2,1],"nx":[3,1,2],"ny":[5,17,0],"nz":[1,-1,-1]},{"size":2,"px":[14,8],"py":[17,6],"pz":[0,1],"nx":[13,10],"ny":[16,9],"nz":[0,-1]},{"size":5,"px":[15,7,14,13,14],"py":[1,0,0,0,1],"pz":[0,1,0,0,0],"nx":[4,4,4,8,8],"ny":[5,3,2,10,10],"nz":[2,2,2,1,-1]},{"size":5,"px":[8,9,4,5,4],"py":[13,12,9,5,7],"pz":[0,0,1,1,1],"nx":[22,21,22,22,22],"ny":[4,0,3,2,2],"nz":[0,0,0,0,-1]},{"size":2,"px":[17,17],"py":[16,13],"pz":[0,0],"nx":[14,21],"ny":[8,0],"nz":[0,-1]},{"size":2,"px":[16,10],"py":[4,9],"pz":[0,-1],"nx":[16,10],"ny":[3,3],"nz":[0,1]},{"size":5,"px":[1,1,0,1,0],"py":[17,16,7,15,8],"pz":[0,0,1,0,0],"nx":[4,3,8,9,7],"ny":[3,3,6,6,6],"nz":[1,1,0,0,-1]},{"size":2,"px":[3,3],"py":[2,3],"pz":[2,2],"nx":[8,3],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[10,2],"py":[17,4],"pz":[0,2],"nx":[10,12],"ny":[15,14],"nz":[0,-1]},{"size":2,"px":[11,11],"py":[14,12],"pz":[0,0],"nx":[9,10],"ny":[13,11],"nz":[0,0]},{"size":2,"px":[12,13],"py":[5,5],"pz":[0,0],"nx":[3,4],"ny":[4,1],"nz":[1,-1]},{"size":5,"px":[7,10,8,11,11],"py":[13,2,12,2,2],"pz":[0,0,0,0,-1],"nx":[10,1,1,10,1],"ny":[12,5,3,13,1],"nz":[0,1,1,0,2]},{"size":2,"px":[6,10],"py":[4,2],"pz":[1,-1],"nx":[4,6],"ny":[4,9],"nz":[1,1]},{"size":2,"px":[20,20],"py":[21,22],"pz":[0,0],"nx":[15,8],"ny":[5,5],"nz":[0,-1]},{"size":2,"px":[4,3],"py":[3,3],"pz":[2,2],"nx":[9,17],"ny":[4,15],"nz":[1,-1]},{"size":3,"px":[2,2,4],"py":[3,3,7],"pz":[2,-1,-1],"nx":[7,4,4],"ny":[6,5,4],"nz":[1,2,2]},{"size":5,"px":[8,9,16,17,17],"py":[1,2,1,1,1],"pz":[1,1,0,0,-1],"nx":[2,2,4,2,4],"ny":[16,14,22,15,21],"nz":[0,0,0,0,0]},{"size":2,"px":[9,9],"py":[18,0],"pz":[0,-1],"nx":[2,5],"ny":[5,8],"nz":[2,1]},{"size":2,"px":[7,8],"py":[11,11],"pz":[0,0],"nx":[15,5],"ny":[8,8],"nz":[0,-1]},{"size":2,"px":[0,3],"py":[4,3],"pz":[2,-1],"nx":[1,6],"ny":[4,14],"nz":[2,0]},{"size":2,"px":[6,12],"py":[7,11],"pz":[1,-1],"nx":[0,0],"ny":[7,12],"nz":[1,0]},{"size":2,"px":[3,7],"py":[10,22],"pz":[1,0],"nx":[4,3],"ny":[10,0],"nz":[1,-1]},{"size":2,"px":[5,19],"py":[4,21],"pz":[2,-1],"nx":[11,11],"ny":[8,9],"nz":[1,1]},{"size":2,"px":[3,3],"py":[8,7],"pz":[1,1],"nx":[4,20],"ny":[4,5],"nz":[1,-1]},{"size":5,"px":[11,23,23,23,23],"py":[7,13,19,20,21],"pz":[1,0,0,0,0],"nx":[4,3,2,8,8],"ny":[11,5,5,23,23],"nz":[1,1,2,0,-1]},{"size":2,"px":[4,1],"py":[0,2],"pz":[0,0],"nx":[0,6],"ny":[0,11],"nz":[0,-1]},{"size":2,"px":[11,8],"py":[12,1],"pz":[0,-1],"nx":[23,23],"ny":[13,12],"nz":[0,0]},{"size":5,"px":[23,11,23,11,11],"py":[13,7,12,5,6],"pz":[0,1,0,1,1],"nx":[6,3,8,7,7],"ny":[12,4,4,11,11],"nz":[0,1,1,0,-1]},{"size":2,"px":[20,5],"py":[15,5],"pz":[0,-1],"nx":[10,10],"ny":[11,10],"nz":[1,1]},{"size":2,"px":[11,4],"py":[19,8],"pz":[0,1],"nx":[11,19],"ny":[18,2],"nz":[0,-1]},{"size":2,"px":[14,6],"py":[3,4],"pz":[0,-1],"nx":[8,15],"ny":[1,0],"nz":[1,0]},{"size":4,"px":[14,5,13,12],"py":[23,3,23,23],"pz":[0,1,0,0],"nx":[12,0,1,4],"ny":[21,3,2,4],"nz":[0,-1,-1,-1]},{"size":2,"px":[19,5],"py":[12,2],"pz":[0,-1],"nx":[4,7],"ny":[3,5],"nz":[2,1]},{"size":2,"px":[0,8],"py":[5,3],"pz":[2,-1],"nx":[5,22],"ny":[3,11],"nz":[2,0]},{"size":2,"px":[2,6],"py":[3,12],"pz":[2,0],"nx":[3,5],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[5,5],"py":[0,6],"pz":[2,-1],"nx":[14,6],"ny":[4,2],"nz":[0,1]},{"size":2,"px":[16,11],"py":[1,0],"pz":[0,-1],"nx":[4,8],"ny":[4,10],"nz":[2,1]},{"size":2,"px":[9,4],"py":[4,3],"pz":[1,1],"nx":[5,8],"ny":[0,10],"nz":[2,-1]},{"size":2,"px":[16,1],"py":[22,1],"pz":[0,-1],"nx":[2,2],"ny":[4,2],"nz":[2,2]},{"size":2,"px":[12,2],"py":[11,2],"pz":[0,-1],"nx":[5,5],"ny":[1,0],"nz":[2,2]},{"size":2,"px":[11,11],"py":[4,3],"pz":[1,1],"nx":[7,5],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[9,2],"py":[22,3],"pz":[0,2],"nx":[4,9],"ny":[10,11],"nz":[1,-1]},{"size":2,"px":[2,4],"py":[8,10],"pz":[1,-1],"nx":[5,3],"ny":[23,18],"nz":[0,0]},{"size":2,"px":[12,6],"py":[21,9],"pz":[0,-1],"nx":[11,23],"ny":[6,10],"nz":[1,0]},{"size":2,"px":[9,9],"py":[8,7],"pz":[1,1],"nx":[18,8],"ny":[18,6],"nz":[0,-1]},{"size":2,"px":[13,3],"py":[19,0],"pz":[0,-1],"nx":[6,5],"ny":[9,11],"nz":[1,1]},{"size":5,"px":[2,10,9,7,8],"py":[0,1,0,1,0],"pz":[2,0,0,0,0],"nx":[3,4,6,8,8],"ny":[2,4,9,4,4],"nz":[2,1,1,1,-1]},{"size":2,"px":[8,4],"py":[6,3],"pz":[1,2],"nx":[9,4],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[0,4],"py":[23,3],"pz":[0,-1],"nx":[12,9],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[4,2],"py":[10,3],"pz":[1,2],"nx":[0,2],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[12,14],"py":[18,0],"pz":[0,-1],"nx":[12,8],"ny":[16,10],"nz":[0,1]},{"size":4,"px":[10,18,7,5],"py":[14,8,0,3],"pz":[0,-1,-1,-1],"nx":[8,6,8,5],"ny":[11,12,5,5],"nz":[0,0,1,1]},{"size":2,"px":[6,5],"py":[2,2],"pz":[1,1],"nx":[8,8],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[12,10],"py":[20,20],"pz":[0,0],"nx":[11,10],"ny":[19,19],"nz":[0,0]},{"size":2,"px":[17,10],"py":[16,20],"pz":[0,-1],"nx":[8,7],"ny":[4,8],"nz":[1,1]},{"size":3,"px":[2,1,3],"py":[20,4,21],"pz":[0,2,0],"nx":[3,4,0],"ny":[10,1,0],"nz":[1,-1,-1]},{"size":5,"px":[6,7,3,6,6],"py":[15,14,7,16,19],"pz":[0,0,1,0,0],"nx":[0,0,0,0,0],"ny":[18,19,16,17,17],"nz":[0,0,0,0,-1]},{"size":2,"px":[8,16],"py":[6,12],"pz":[1,0],"nx":[8,15],"ny":[4,10],"nz":[1,-1]},{"size":5,"px":[0,0,0,0,0],"py":[1,3,2,0,4],"pz":[2,2,2,2,1],"nx":[13,8,14,4,7],"ny":[23,6,23,3,9],"nz":[0,1,0,2,-1]},{"size":2,"px":[3,6],"py":[3,5],"pz":[2,1],"nx":[10,8],"ny":[11,6],"nz":[0,-1]},{"size":2,"px":[11,10],"py":[4,4],"pz":[0,0],"nx":[8,5],"ny":[4,9],"nz":[1,-1]},{"size":5,"px":[15,18,9,16,4],"py":[12,13,6,23,3],"pz":[0,0,1,0,2],"nx":[6,3,6,2,7],"ny":[2,3,0,1,0],"nz":[0,0,0,1,0]},{"size":2,"px":[4,18],"py":[12,13],"pz":[0,-1],"nx":[2,8],"ny":[3,4],"nz":[2,1]},{"size":2,"px":[4,2],"py":[10,4],"pz":[1,2],"nx":[3,3],"ny":[5,0],"nz":[2,-1]},{"size":2,"px":[9,19],"py":[7,8],"pz":[1,0],"nx":[8,3],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[6,0],"py":[6,0],"pz":[0,-1],"nx":[0,0],"ny":[7,2],"nz":[1,2]},{"size":2,"px":[8,8],"py":[0,0],"pz":[1,-1],"nx":[17,18],"ny":[0,2],"nz":[0,0]},{"size":4,"px":[13,4,4,1],"py":[14,7,3,5],"pz":[0,-1,-1,-1],"nx":[3,16,3,7],"ny":[1,15,5,13],"nz":[2,0,2,0]},{"size":2,"px":[4,9],"py":[6,11],"pz":[1,0],"nx":[3,23],"ny":[4,8],"nz":[1,-1]},{"size":5,"px":[9,17,4,16,16],"py":[2,3,1,3,3],"pz":[1,0,2,0,-1],"nx":[2,3,3,2,3],"ny":[1,7,2,3,3],"nz":[2,1,1,1,1]},{"size":2,"px":[10,5],"py":[22,9],"pz":[0,1],"nx":[10,3],"ny":[21,2],"nz":[0,-1]},{"size":2,"px":[11,11],"py":[6,3],"pz":[0,-1],"nx":[8,5],"ny":[4,3],"nz":[1,1]},{"size":2,"px":[10,5],"py":[8,3],"pz":[0,-1],"nx":[14,5],"ny":[14,2],"nz":[0,2]},{"size":2,"px":[7,8],"py":[3,2],"pz":[0,-1],"nx":[8,2],"ny":[18,2],"nz":[0,2]},{"size":2,"px":[1,1],"py":[19,11],"pz":[0,1],"nx":[9,4],"ny":[5,1],"nz":[0,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,3],"ny":[4,4],"nz":[1,-1]},{"size":5,"px":[7,15,13,14,4],"py":[6,12,9,11,4],"pz":[1,0,0,0,2],"nx":[7,3,8,4,5],"ny":[0,3,0,2,1],"nz":[0,0,0,0,0]},{"size":5,"px":[10,13,7,8,9],"py":[0,1,1,0,1],"pz":[0,0,0,0,0],"nx":[7,4,4,4,8],"ny":[8,3,4,2,4],"nz":[1,2,2,2,1]},{"size":2,"px":[6,1],"py":[6,0],"pz":[1,-1],"nx":[11,7],"ny":[3,2],"nz":[0,1]},{"size":2,"px":[13,0],"py":[13,2],"pz":[0,-1],"nx":[0,1],"ny":[13,16],"nz":[0,0]},{"size":2,"px":[8,17],"py":[6,13],"pz":[1,0],"nx":[8,1],"ny":[4,16],"nz":[1,-1]},{"size":5,"px":[12,11,3,6,17],"py":[4,4,1,2,14],"pz":[0,0,2,1,0],"nx":[6,23,23,6,23],"ny":[5,7,6,6,14],"nz":[1,0,0,1,0]},{"size":2,"px":[5,22],"py":[4,17],"pz":[2,-1],"nx":[4,8],"ny":[5,7],"nz":[2,1]},{"size":2,"px":[15,14],"py":[1,1],"pz":[0,0],"nx":[4,7],"ny":[2,4],"nz":[2,-1]},{"size":2,"px":[15,17],"py":[12,7],"pz":[0,-1],"nx":[14,10],"ny":[11,4],"nz":[0,1]},{"size":4,"px":[10,2,9,15],"py":[5,11,1,13],"pz":[0,-1,-1,-1],"nx":[11,3,3,13],"ny":[1,1,0,1],"nz":[0,2,2,0]},{"size":2,"px":[7,21],"py":[15,22],"pz":[0,-1],"nx":[4,9],"ny":[8,14],"nz":[1,0]},{"size":2,"px":[6,5],"py":[21,2],"pz":[0,-1],"nx":[3,5],"ny":[11,21],"nz":[1,0]},{"size":2,"px":[17,7],"py":[2,0],"pz":[0,-1],"nx":[4,8],"ny":[5,11],"nz":[2,1]},{"size":2,"px":[11,8],"py":[10,4],"pz":[0,-1],"nx":[13,12],"ny":[3,3],"nz":[0,0]},{"size":2,"px":[6,5],"py":[2,2],"pz":[1,1],"nx":[7,1],"ny":[8,2],"nz":[0,-1]},{"size":5,"px":[0,0,1,0,0],"py":[12,4,14,0,2],"pz":[0,1,0,2,2],"nx":[9,5,8,4,4],"ny":[6,3,6,3,3],"nz":[0,1,0,1,-1]},{"size":5,"px":[8,0,0,3,2],"py":[6,5,0,8,2],"pz":[1,-1,-1,-1,-1],"nx":[23,7,22,11,4],"ny":[12,6,14,4,3],"nz":[0,1,0,1,2]},{"size":4,"px":[12,12,4,8],"py":[12,11,3,10],"pz":[0,0,-1,-1],"nx":[0,0,0,0],"ny":[2,1,0,3],"nz":[1,2,2,1]},{"size":2,"px":[10,6],"py":[7,6],"pz":[1,-1],"nx":[16,4],"ny":[12,2],"nz":[0,2]},{"size":5,"px":[2,1,3,3,3],"py":[14,8,20,21,21],"pz":[0,1,0,0,-1],"nx":[20,10,21,21,21],"ny":[23,11,21,23,20],"nz":[0,1,0,0,0]},{"size":2,"px":[6,13],"py":[2,4],"pz":[1,0],"nx":[7,21],"ny":[8,0],"nz":[0,-1]},{"size":2,"px":[12,3],"py":[17,4],"pz":[0,2],"nx":[11,10],"ny":[15,7],"nz":[0,-1]},{"size":4,"px":[11,0,19,2],"py":[15,2,23,10],"pz":[0,-1,-1,-1],"nx":[6,8,16,2],"ny":[13,11,10,2],"nz":[0,0,0,2]},{"size":2,"px":[6,3],"py":[14,7],"pz":[0,1],"nx":[3,1],"ny":[4,1],"nz":[1,-1]},{"size":4,"px":[12,17,5,10],"py":[19,15,14,3],"pz":[0,-1,-1,-1],"nx":[4,12,6,12],"ny":[4,18,9,22],"nz":[1,0,1,0]},{"size":2,"px":[8,3],"py":[13,5],"pz":[0,-1],"nx":[3,4],"ny":[4,9],"nz":[1,1]},{"size":5,"px":[6,5,4,5,3],"py":[2,1,2,2,0],"pz":[0,0,0,0,1],"nx":[7,4,9,18,18],"ny":[4,4,7,14,14],"nz":[1,1,1,0,-1]},{"size":4,"px":[8,3,20,1],"py":[6,3,18,0],"pz":[1,-1,-1,-1],"nx":[13,11,5,22],"ny":[12,6,2,17],"nz":[0,1,2,0]},{"size":2,"px":[6,3],"py":[6,3],"pz":[1,2],"nx":[8,5],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[21,7],"py":[14,7],"pz":[0,1],"nx":[16,11],"ny":[14,6],"nz":[0,-1]},{"size":2,"px":[10,4],"py":[3,1],"pz":[0,-1],"nx":[9,5],"ny":[0,0],"nz":[0,1]},{"size":2,"px":[4,10],"py":[5,8],"pz":[2,1],"nx":[5,14],"ny":[9,7],"nz":[1,-1]},{"size":2,"px":[9,2],"py":[23,4],"pz":[0,2],"nx":[2,2],"ny":[5,5],"nz":[2,-1]},{"size":5,"px":[10,9,11,10,10],"py":[2,2,1,1,1],"pz":[0,0,0,0,-1],"nx":[2,3,2,4,5],"ny":[4,10,2,4,3],"nz":[2,1,1,0,0]},{"size":2,"px":[11,4],"py":[13,4],"pz":[0,-1],"nx":[8,4],"ny":[4,1],"nz":[1,2]},{"size":2,"px":[17,5],"py":[15,1],"pz":[0,-1],"nx":[20,19],"ny":[14,14],"nz":[0,0]},{"size":2,"px":[2,2],"py":[20,18],"pz":[0,0],"nx":[2,1],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[10,1],"py":[18,3],"pz":[0,2],"nx":[11,3],"ny":[16,5],"nz":[0,-1]},{"size":2,"px":[3,8],"py":[6,10],"pz":[1,0],"nx":[9,0],"ny":[9,3],"nz":[0,-1]},{"size":2,"px":[20,10],"py":[21,7],"pz":[0,1],"nx":[7,2],"ny":[3,5],"nz":[1,-1]},{"size":2,"px":[10,6],"py":[4,7],"pz":[1,-1],"nx":[23,5],"ny":[9,2],"nz":[0,2]},{"size":5,"px":[2,4,5,3,4],"py":[0,1,1,2,2],"pz":[1,0,0,0,0],"nx":[1,0,1,1,1],"ny":[2,1,0,1,1],"nz":[0,1,0,0,-1]},{"size":2,"px":[8,16],"py":[7,13],"pz":[1,0],"nx":[8,3],"ny":[4,16],"nz":[1,-1]},{"size":2,"px":[17,15],"py":[7,19],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[4,3],"py":[11,5],"pz":[1,2],"nx":[7,8],"ny":[9,4],"nz":[1,-1]},{"size":2,"px":[23,11],"py":[9,6],"pz":[0,1],"nx":[22,22],"ny":[23,23],"nz":[0,-1]},{"size":2,"px":[23,23],"py":[21,20],"pz":[0,0],"nx":[2,2],"ny":[5,4],"nz":[1,-1]},{"size":2,"px":[17,4],"py":[12,2],"pz":[0,-1],"nx":[9,8],"ny":[4,5],"nz":[1,1]},{"size":2,"px":[6,14],"py":[2,4],"pz":[1,0],"nx":[7,18],"ny":[1,1],"nz":[1,-1]},{"size":2,"px":[20,22],"py":[1,2],"pz":[0,0],"nx":[23,23],"ny":[1,1],"nz":[0,-1]},{"size":2,"px":[0,1],"py":[9,10],"pz":[1,1],"nx":[8,0],"ny":[15,0],"nz":[0,-1]},{"size":3,"px":[11,11,6],"py":[10,11,11],"pz":[0,0,-1],"nx":[23,23,23],"ny":[19,21,20],"nz":[0,0,0]},{"size":5,"px":[23,23,23,6,6],"py":[21,22,22,3,6],"pz":[0,0,-1,-1,-1],"nx":[8,8,8,17,4],"ny":[7,10,8,16,5],"nz":[1,1,1,0,2]},{"size":2,"px":[10,23],"py":[1,22],"pz":[0,-1],"nx":[7,2],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[7,14],"py":[3,10],"pz":[1,-1],"nx":[5,3],"ny":[2,1],"nz":[0,1]},{"size":2,"px":[5,3],"py":[13,7],"pz":[0,1],"nx":[4,10],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[10,0],"py":[15,6],"pz":[0,-1],"nx":[3,6],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[13,4],"py":[18,17],"pz":[0,-1],"nx":[7,6],"ny":[10,7],"nz":[1,1]},{"size":2,"px":[12,11],"py":[3,8],"pz":[0,-1],"nx":[7,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[17,4],"py":[5,7],"pz":[0,1],"nx":[17,10],"ny":[4,0],"nz":[0,-1]},{"size":5,"px":[16,8,16,15,15],"py":[0,0,1,0,1],"pz":[0,1,0,0,0],"nx":[7,4,7,4,4],"ny":[7,5,8,1,1],"nz":[1,2,1,2,-1]},{"size":2,"px":[13,11],"py":[5,6],"pz":[0,-1],"nx":[4,5],"ny":[2,2],"nz":[1,1]},{"size":2,"px":[3,6],"py":[3,6],"pz":[2,1],"nx":[8,4],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[10,16],"py":[8,10],"pz":[0,0],"nx":[7,2],"ny":[3,3],"nz":[1,-1]},{"size":2,"px":[6,8],"py":[4,11],"pz":[1,0],"nx":[10,1],"ny":[9,20],"nz":[0,-1]},{"size":2,"px":[5,1],"py":[4,2],"pz":[2,-1],"nx":[23,23],"ny":[15,16],"nz":[0,0]},{"size":5,"px":[9,8,2,4,9],"py":[1,1,0,1,2],"pz":[0,0,2,1,0],"nx":[8,3,8,4,4],"ny":[6,2,4,2,2],"nz":[1,2,1,2,-1]},{"size":2,"px":[13,6],"py":[10,5],"pz":[0,-1],"nx":[13,7],"ny":[6,3],"nz":[0,1]},{"size":2,"px":[11,5],"py":[10,5],"pz":[1,2],"nx":[10,8],"ny":[10,9],"nz":[1,-1]},{"size":2,"px":[7,4],"py":[6,3],"pz":[1,2],"nx":[9,14],"ny":[4,9],"nz":[1,-1]},{"size":3,"px":[5,2,15],"py":[3,1,22],"pz":[1,-1,-1],"nx":[15,9,4],"ny":[0,1,0],"nz":[0,1,2]},{"size":2,"px":[10,19],"py":[9,21],"pz":[1,0],"nx":[2,17],"ny":[5,14],"nz":[2,-1]},{"size":3,"px":[16,2,1],"py":[2,10,4],"pz":[0,-1,-1],"nx":[4,4,9],"ny":[3,2,6],"nz":[2,2,1]},{"size":2,"px":[10,2],"py":[6,10],"pz":[1,-1],"nx":[21,22],"ny":[16,12],"nz":[0,0]},{"size":2,"px":[7,16],"py":[4,23],"pz":[0,-1],"nx":[7,3],"ny":[3,3],"nz":[0,1]},{"size":2,"px":[1,1],"py":[13,14],"pz":[0,0],"nx":[1,2],"ny":[18,3],"nz":[0,-1]},{"size":2,"px":[18,5],"py":[13,4],"pz":[0,-1],"nx":[4,13],"ny":[2,11],"nz":[2,0]},{"size":2,"px":[18,17],"py":[3,3],"pz":[0,0],"nx":[19,19],"ny":[1,1],"nz":[0,-1]},{"size":2,"px":[9,5],"py":[0,5],"pz":[1,-1],"nx":[12,3],"ny":[5,1],"nz":[0,2]},{"size":2,"px":[5,3],"py":[2,1],"pz":[1,2],"nx":[18,4],"ny":[4,1],"nz":[0,-1]},{"size":5,"px":[13,13,2,10,15],"py":[11,12,13,17,23],"pz":[0,-1,-1,-1,-1],"nx":[12,13,4,3,8],"ny":[4,4,1,0,3],"nz":[0,0,2,2,1]},{"size":2,"px":[9,3],"py":[2,2],"pz":[0,-1],"nx":[4,2],"ny":[7,2],"nz":[1,2]},{"size":2,"px":[13,4],"py":[5,1],"pz":[0,-1],"nx":[18,4],"ny":[12,2],"nz":[0,2]},{"size":2,"px":[19,4],"py":[11,1],"pz":[0,-1],"nx":[4,7],"ny":[2,2],"nz":[2,1]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,2],"ny":[4,5],"nz":[1,-1]},{"size":2,"px":[4,0],"py":[7,7],"pz":[0,-1],"nx":[4,9],"ny":[0,2],"nz":[2,1]},{"size":2,"px":[4,9],"py":[0,2],"pz":[2,1],"nx":[6,4],"ny":[3,4],"nz":[0,-1]},{"size":2,"px":[4,2],"py":[9,4],"pz":[1,2],"nx":[13,5],"ny":[18,2],"nz":[0,-1]},{"size":3,"px":[5,23,23],"py":[2,8,7],"pz":[2,0,0],"nx":[10,12,1],"ny":[4,1,0],"nz":[1,-1,-1]},{"size":2,"px":[13,0],"py":[3,3],"pz":[0,-1],"nx":[4,4],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[6,5],"py":[10,5],"pz":[0,-1],"nx":[0,0],"ny":[4,11],"nz":[1,0]},{"size":2,"px":[11,2],"py":[14,11],"pz":[0,-1],"nx":[10,11],"ny":[4,13],"nz":[1,0]},{"size":2,"px":[5,6],"py":[21,23],"pz":[0,0],"nx":[7,0],"ny":[21,3],"nz":[0,-1]},{"size":2,"px":[8,4],"py":[6,3],"pz":[1,2],"nx":[8,5],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[7,6],"py":[8,8],"pz":[0,0],"nx":[6,14],"ny":[9,15],"nz":[0,-1]},{"size":2,"px":[16,6],"py":[4,8],"pz":[0,-1],"nx":[16,8],"ny":[0,1],"nz":[0,1]},{"size":4,"px":[3,6,0,9],"py":[0,8,5,23],"pz":[1,-1,-1,-1],"nx":[12,2,6,10],"ny":[5,0,3,5],"nz":[0,2,1,0]},{"size":2,"px":[3,6],"py":[7,13],"pz":[1,0],"nx":[3,9],"ny":[4,9],"nz":[1,-1]},{"size":2,"px":[2,5],"py":[8,23],"pz":[1,0],"nx":[8,9],"ny":[15,0],"nz":[0,-1]},{"size":2,"px":[13,18],"py":[8,0],"pz":[0,-1],"nx":[1,1],"ny":[9,8],"nz":[1,1]},{"size":2,"px":[2,7],"py":[4,21],"pz":[2,0],"nx":[13,11],"ny":[8,9],"nz":[0,-1]},{"size":2,"px":[5,4],"py":[8,8],"pz":[0,0],"nx":[6,1],"ny":[8,5],"nz":[0,-1]},{"size":2,"px":[7,3],"py":[20,7],"pz":[0,-1],"nx":[4,3],"ny":[10,4],"nz":[1,1]},{"size":2,"px":[9,9],"py":[8,7],"pz":[1,-1],"nx":[1,2],"ny":[4,9],"nz":[2,1]},{"size":2,"px":[5,10],"py":[5,13],"pz":[1,-1],"nx":[3,6],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[12,5],"py":[6,3],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[10,10],"py":[4,4],"pz":[1,-1],"nx":[5,11],"ny":[2,5],"nz":[2,1]},{"size":5,"px":[11,23,11,23,11],"py":[4,9,5,10,6],"pz":[1,0,1,0,1],"nx":[7,14,13,7,3],"ny":[9,5,6,4,4],"nz":[0,0,0,1,-1]},{"size":2,"px":[8,5],"py":[0,0],"pz":[1,-1],"nx":[9,20],"ny":[1,4],"nz":[1,0]},{"size":2,"px":[19,20],"py":[0,3],"pz":[0,0],"nx":[4,6],"ny":[11,3],"nz":[1,-1]},{"size":4,"px":[13,5,20,5],"py":[14,3,23,4],"pz":[0,-1,-1,-1],"nx":[8,15,7,16],"ny":[8,14,6,15],"nz":[1,0,1,0]},{"size":2,"px":[10,20],"py":[5,17],"pz":[0,-1],"nx":[7,3],"ny":[10,1],"nz":[0,2]},{"size":3,"px":[1,12,7],"py":[3,7,10],"pz":[2,0,0],"nx":[2,2,3],"ny":[3,2,2],"nz":[1,-1,-1]},{"size":3,"px":[10,5,7],"py":[7,10,10],"pz":[1,-1,-1],"nx":[10,10,18],"ny":[10,9,23],"nz":[1,1,0]},{"size":3,"px":[14,14,4],"py":[3,3,4],"pz":[0,-1,-1],"nx":[4,4,8],"ny":[3,2,6],"nz":[2,2,1]},{"size":2,"px":[4,12],"py":[4,17],"pz":[2,0],"nx":[13,1],"ny":[15,4],"nz":[0,-1]},{"size":2,"px":[10,20],"py":[9,22],"pz":[0,-1],"nx":[9,4],"ny":[2,0],"nz":[1,2]},{"size":2,"px":[11,2],"py":[3,6],"pz":[0,-1],"nx":[2,4],"ny":[2,4],"nz":[2,1]},{"size":3,"px":[15,10,1],"py":[12,2,3],"pz":[0,-1,-1],"nx":[7,5,10],"ny":[2,1,1],"nz":[0,1,0]},{"size":5,"px":[9,11,10,12,12],"py":[0,0,0,0,0],"pz":[0,0,0,0,-1],"nx":[8,4,16,5,10],"ny":[4,4,10,3,6],"nz":[1,1,0,1,0]},{"size":2,"px":[0,10],"py":[3,5],"pz":[2,-1],"nx":[3,6],"ny":[0,1],"nz":[2,1]},{"size":5,"px":[7,8,7,2,12],"py":[14,13,13,16,0],"pz":[0,0,-1,-1,-1],"nx":[10,1,10,1,1],"ny":[13,2,12,4,9],"nz":[0,2,0,1,0]},{"size":3,"px":[6,14,13],"py":[1,2,1],"pz":[1,0,0],"nx":[8,21,10],"ny":[4,23,12],"nz":[1,-1,-1]},{"size":2,"px":[19,19],"py":[22,21],"pz":[0,0],"nx":[20,1],"ny":[22,5],"nz":[0,-1]},{"size":2,"px":[13,12],"py":[19,22],"pz":[0,-1],"nx":[2,3],"ny":[0,1],"nz":[2,1]},{"size":4,"px":[11,9,21,4],"py":[13,3,19,5],"pz":[0,-1,-1,-1],"nx":[9,9,9,5],"ny":[13,14,12,6],"nz":[0,0,0,1]},{"size":4,"px":[11,12,13,14],"py":[22,22,22,22],"pz":[0,0,0,0],"nx":[13,2,4,5],"ny":[20,0,0,6],"nz":[0,-1,-1,-1]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,1],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[0,1],"pz":[2,2],"nx":[9,4],"ny":[6,5],"nz":[1,-1]},{"size":2,"px":[17,0],"py":[10,1],"pz":[0,-1],"nx":[9,4],"ny":[3,2],"nz":[1,2]},{"size":2,"px":[10,4],"py":[3,1],"pz":[1,2],"nx":[12,18],"ny":[17,4],"nz":[0,-1]},{"size":3,"px":[2,3,4],"py":[4,3,9],"pz":[2,2,1],"nx":[0,3,17],"ny":[0,1,18],"nz":[0,-1,-1]},{"size":2,"px":[7,3],"py":[12,6],"pz":[0,1],"nx":[5,1],"ny":[11,1],"nz":[1,-1]},{"size":2,"px":[10,17],"py":[20,6],"pz":[0,-1],"nx":[5,2],"ny":[9,5],"nz":[1,2]},{"size":2,"px":[8,11],"py":[18,2],"pz":[0,-1],"nx":[5,4],"ny":[9,9],"nz":[1,1]},{"size":2,"px":[16,15],"py":[2,2],"pz":[0,0],"nx":[17,12],"ny":[2,2],"nz":[0,-1]},{"size":2,"px":[18,4],"py":[5,5],"pz":[0,-1],"nx":[7,5],"ny":[23,19],"nz":[0,0]},{"size":2,"px":[12,13],"py":[23,23],"pz":[0,0],"nx":[7,11],"ny":[10,20],"nz":[1,-1]},{"size":2,"px":[5,10],"py":[3,18],"pz":[2,-1],"nx":[9,9],"ny":[5,6],"nz":[1,1]},{"size":2,"px":[5,10],"py":[2,4],"pz":[1,0],"nx":[4,23],"ny":[4,20],"nz":[1,-1]},{"size":2,"px":[2,3],"py":[8,1],"pz":[1,-1],"nx":[15,12],"ny":[2,1],"nz":[0,0]},{"size":2,"px":[4,7],"py":[3,10],"pz":[2,1],"nx":[10,1],"ny":[20,4],"nz":[0,-1]},{"size":2,"px":[11,11],"py":[10,11],"pz":[0,0],"nx":[22,3],"ny":[5,4],"nz":[0,-1]},{"size":5,"px":[8,17,17,9,18],"py":[0,1,0,1,0],"pz":[1,0,0,1,0],"nx":[11,8,9,4,4],"ny":[23,4,6,2,2],"nz":[0,1,0,2,-1]},{"size":2,"px":[5,5],"py":[4,4],"pz":[1,-1],"nx":[13,4],"ny":[9,2],"nz":[0,2]},{"size":5,"px":[9,4,8,7,7],"py":[3,1,3,3,3],"pz":[0,1,0,0,-1],"nx":[4,2,5,3,2],"ny":[1,15,1,4,13],"nz":[0,0,0,0,0]},{"size":2,"px":[17,7],"py":[13,7],"pz":[0,-1],"nx":[4,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[1,2],"py":[1,12],"pz":[2,0],"nx":[9,21],"ny":[5,4],"nz":[0,-1]},{"size":2,"px":[12,0],"py":[14,1],"pz":[0,-1],"nx":[1,1],"ny":[19,10],"nz":[0,1]},{"size":2,"px":[16,1],"py":[5,9],"pz":[0,-1],"nx":[16,15],"ny":[3,3],"nz":[0,0]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[8,4],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[11,6],"py":[17,15],"pz":[0,0],"nx":[11,0],"ny":[16,4],"nz":[0,-1]},{"size":4,"px":[12,11,0,3],"py":[16,8,7,1],"pz":[0,-1,-1,-1],"nx":[10,5,10,5],"ny":[11,9,10,8],"nz":[0,1,0,1]},{"size":2,"px":[3,6],"py":[7,13],"pz":[1,0],"nx":[4,14],"ny":[4,16],"nz":[1,-1]},{"size":2,"px":[7,17],"py":[6,13],"pz":[0,-1],"nx":[4,8],"ny":[4,9],"nz":[2,1]},{"size":2,"px":[15,11],"py":[3,2],"pz":[0,-1],"nx":[4,15],"ny":[1,2],"nz":[2,0]},{"size":2,"px":[10,11],"py":[18,4],"pz":[0,-1],"nx":[5,5],"ny":[8,9],"nz":[1,1]},{"size":2,"px":[8,4],"py":[7,4],"pz":[1,2],"nx":[4,3],"ny":[5,7],"nz":[2,-1]},{"size":2,"px":[12,4],"py":[15,4],"pz":[0,-1],"nx":[11,8],"ny":[14,19],"nz":[0,0]},{"size":2,"px":[18,13],"py":[13,20],"pz":[0,0],"nx":[13,4],"ny":[18,2],"nz":[0,-1]},{"size":2,"px":[12,4],"py":[6,3],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[21,5,11,5,10],"py":[1,1,3,0,0],"pz":[0,2,1,2,1],"nx":[7,14,15,4,8],"ny":[3,6,11,3,4],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[10,6],"py":[15,10],"pz":[0,-1],"nx":[21,22],"ny":[14,12],"nz":[0,0]},{"size":2,"px":[18,0],"py":[20,0],"pz":[0,-1],"nx":[2,3],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[12,6,13,11,7],"py":[1,1,1,2,1],"pz":[0,1,0,0,1],"nx":[7,6,8,5,5],"ny":[4,15,4,16,16],"nz":[1,0,1,0,-1]},{"size":3,"px":[22,21,21],"py":[14,15,17],"pz":[0,0,0],"nx":[5,9,4],"ny":[0,5,0],"nz":[2,-1,-1]},{"size":2,"px":[10,2],"py":[14,1],"pz":[0,-1],"nx":[23,11],"ny":[16,8],"nz":[0,1]},{"size":4,"px":[21,21,0,18],"py":[14,15,5,4],"pz":[0,0,-1,-1],"nx":[8,8,9,4],"ny":[7,8,10,5],"nz":[1,1,1,2]},{"size":2,"px":[15,5],"py":[18,1],"pz":[0,-1],"nx":[23,23],"ny":[16,18],"nz":[0,0]},{"size":2,"px":[15,14],"py":[1,1],"pz":[0,0],"nx":[4,4],"ny":[2,3],"nz":[2,-1]},{"size":2,"px":[2,6],"py":[6,5],"pz":[1,-1],"nx":[14,11],"ny":[1,1],"nz":[0,0]},{"size":2,"px":[3,17],"py":[2,8],"pz":[2,0],"nx":[8,3],"ny":[4,9],"nz":[1,-1]},{"size":2,"px":[17,8],"py":[13,10],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[0,0],"py":[8,3],"pz":[0,1],"nx":[1,11],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[6,8],"py":[5,0],"pz":[1,-1],"nx":[0,0],"ny":[3,1],"nz":[1,2]},{"size":2,"px":[0,0],"py":[5,3],"pz":[1,2],"nx":[1,18],"ny":[5,7],"nz":[1,-1]},{"size":2,"px":[7,3],"py":[6,6],"pz":[0,1],"nx":[7,12],"ny":[5,20],"nz":[0,-1]},{"size":2,"px":[8,1],"py":[0,5],"pz":[0,-1],"nx":[4,2],"ny":[9,3],"nz":[1,2]},{"size":2,"px":[0,0],"py":[10,11],"pz":[0,0],"nx":[0,5],"ny":[5,9],"nz":[0,-1]},{"size":2,"px":[8,1],"py":[23,4],"pz":[0,2],"nx":[0,0],"ny":[13,2],"nz":[0,-1]},{"size":2,"px":[4,1],"py":[6,4],"pz":[0,-1],"nx":[4,4],"ny":[4,5],"nz":[2,2]},{"size":2,"px":[7,6],"py":[6,5],"pz":[1,1],"nx":[3,9],"ny":[4,16],"nz":[1,-1]},{"size":2,"px":[5,3],"py":[9,13],"pz":[0,-1],"nx":[4,10],"ny":[3,7],"nz":[1,0]},{"size":5,"px":[13,9,6,10,10],"py":[2,2,1,2,2],"pz":[0,0,1,0,-1],"nx":[7,5,6,5,6],"ny":[0,2,2,1,1],"nz":[0,0,0,0,0]}],"alpha":[-1.119615e+00,1.119615e+00,-8.169953e-01,8.169953e-01,-5.291213e-01,5.291213e-01,-4.904488e-01,4.904488e-01,-4.930982e-01,4.930982e-01,-4.106179e-01,4.106179e-01,-4.246842e-01,4.246842e-01,-3.802383e-01,3.802383e-01,-3.364358e-01,3.364358e-01,-3.214186e-01,3.214186e-01,-3.210798e-01,3.210798e-01,-2.993167e-01,2.993167e-01,-3.426336e-01,3.426336e-01,-3.199184e-01,3.199184e-01,-3.061071e-01,3.061071e-01,-2.758972e-01,2.758972e-01,-3.075590e-01,3.075590e-01,-3.009565e-01,3.009565e-01,-2.015739e-01,2.015739e-01,-2.603266e-01,2.603266e-01,-2.772993e-01,2.772993e-01,-2.184913e-01,2.184913e-01,-2.306681e-01,2.306681e-01,-1.983223e-01,1.983223e-01,-2.194760e-01,2.194760e-01,-2.528421e-01,2.528421e-01,-2.436416e-01,2.436416e-01,-3.032886e-01,3.032886e-01,-2.556071e-01,2.556071e-01,-2.562170e-01,2.562170e-01,-1.930298e-01,1.930298e-01,-2.735898e-01,2.735898e-01,-1.814703e-01,1.814703e-01,-2.054824e-01,2.054824e-01,-1.986146e-01,1.986146e-01,-1.769226e-01,1.769226e-01,-1.775257e-01,1.775257e-01,-2.167927e-01,2.167927e-01,-1.823633e-01,1.823633e-01,-1.584280e-01,1.584280e-01,-1.778321e-01,1.778321e-01,-1.826777e-01,1.826777e-01,-1.979903e-01,1.979903e-01,-1.898326e-01,1.898326e-01,-1.835506e-01,1.835506e-01,-1.967860e-01,1.967860e-01,-1.871528e-01,1.871528e-01,-1.772414e-01,1.772414e-01,-1.985514e-01,1.985514e-01,-2.144078e-01,2.144078e-01,-2.742303e-01,2.742303e-01,-2.240550e-01,2.240550e-01,-2.132534e-01,2.132534e-01,-1.552127e-01,1.552127e-01,-1.568276e-01,1.568276e-01,-1.630086e-01,1.630086e-01,-1.458232e-01,1.458232e-01,-1.559541e-01,1.559541e-01,-1.720131e-01,1.720131e-01,-1.708434e-01,1.708434e-01,-1.624431e-01,1.624431e-01,-1.814161e-01,1.814161e-01,-1.552639e-01,1.552639e-01,-1.242354e-01,1.242354e-01,-1.552139e-01,1.552139e-01,-1.694359e-01,1.694359e-01,-1.801481e-01,1.801481e-01,-1.387182e-01,1.387182e-01,-1.409679e-01,1.409679e-01,-1.486724e-01,1.486724e-01,-1.779553e-01,1.779553e-01,-1.524595e-01,1.524595e-01,-1.788086e-01,1.788086e-01,-1.671479e-01,1.671479e-01,-1.376197e-01,1.376197e-01,-1.511808e-01,1.511808e-01,-1.524632e-01,1.524632e-01,-1.198986e-01,1.198986e-01,-1.382641e-01,1.382641e-01,-1.148901e-01,1.148901e-01,-1.131803e-01,1.131803e-01,-1.273508e-01,1.273508e-01,-1.405125e-01,1.405125e-01,-1.322132e-01,1.322132e-01,-1.386966e-01,1.386966e-01,-1.275621e-01,1.275621e-01,-1.180573e-01,1.180573e-01,-1.238803e-01,1.238803e-01,-1.428389e-01,1.428389e-01,-1.694437e-01,1.694437e-01,-1.290855e-01,1.290855e-01,-1.520260e-01,1.520260e-01,-1.398282e-01,1.398282e-01,-1.890736e-01,1.890736e-01,-2.280428e-01,2.280428e-01,-1.325099e-01,1.325099e-01,-1.342873e-01,1.342873e-01,-1.463841e-01,1.463841e-01,-1.983567e-01,1.983567e-01,-1.585711e-01,1.585711e-01,-1.260154e-01,1.260154e-01,-1.426774e-01,1.426774e-01,-1.554278e-01,1.554278e-01,-1.361201e-01,1.361201e-01,-1.181856e-01,1.181856e-01,-1.255941e-01,1.255941e-01,-1.113275e-01,1.113275e-01,-1.506576e-01,1.506576e-01,-1.202859e-01,1.202859e-01,-2.159751e-01,2.159751e-01,-1.443150e-01,1.443150e-01,-1.379194e-01,1.379194e-01,-1.805758e-01,1.805758e-01,-1.465612e-01,1.465612e-01,-1.328856e-01,1.328856e-01,-1.532173e-01,1.532173e-01,-1.590635e-01,1.590635e-01,-1.462229e-01,1.462229e-01,-1.350012e-01,1.350012e-01,-1.195634e-01,1.195634e-01,-1.173221e-01,1.173221e-01,-1.192867e-01,1.192867e-01,-1.595013e-01,1.595013e-01,-1.209751e-01,1.209751e-01,-1.571290e-01,1.571290e-01,-1.527274e-01,1.527274e-01,-1.373708e-01,1.373708e-01,-1.318313e-01,1.318313e-01,-1.273391e-01,1.273391e-01,-1.271365e-01,1.271365e-01,-1.528693e-01,1.528693e-01,-1.590476e-01,1.590476e-01,-1.581911e-01,1.581911e-01,-1.183023e-01,1.183023e-01,-1.559822e-01,1.559822e-01,-1.214999e-01,1.214999e-01,-1.283378e-01,1.283378e-01,-1.542583e-01,1.542583e-01,-1.336377e-01,1.336377e-01,-1.800416e-01,1.800416e-01,-1.710931e-01,1.710931e-01,-1.621737e-01,1.621737e-01,-1.239002e-01,1.239002e-01,-1.432928e-01,1.432928e-01,-1.392447e-01,1.392447e-01,-1.383938e-01,1.383938e-01,-1.357633e-01,1.357633e-01,-1.175842e-01,1.175842e-01,-1.085318e-01,1.085318e-01,-1.148885e-01,1.148885e-01,-1.320396e-01,1.320396e-01,-1.351204e-01,1.351204e-01,-1.581518e-01,1.581518e-01,-1.459574e-01,1.459574e-01,-1.180068e-01,1.180068e-01,-1.464196e-01,1.464196e-01,-1.179543e-01,1.179543e-01,-1.004204e-01,1.004204e-01,-1.294660e-01,1.294660e-01,-1.534244e-01,1.534244e-01,-1.378970e-01,1.378970e-01,-1.226545e-01,1.226545e-01,-1.281182e-01,1.281182e-01,-1.201471e-01,1.201471e-01,-1.448701e-01,1.448701e-01,-1.290980e-01,1.290980e-01,-1.388764e-01,1.388764e-01,-9.605773e-02,9.605773e-02,-1.411021e-01,1.411021e-01,-1.295693e-01,1.295693e-01,-1.371739e-01,1.371739e-01,-1.167579e-01,1.167579e-01,-1.400486e-01,1.400486e-01,-1.214224e-01,1.214224e-01,-1.287835e-01,1.287835e-01,-1.197646e-01,1.197646e-01,-1.192358e-01,1.192358e-01,-1.218651e-01,1.218651e-01,-1.564816e-01,1.564816e-01,-1.172391e-01,1.172391e-01,-1.342268e-01,1.342268e-01,-1.492471e-01,1.492471e-01,-1.157299e-01,1.157299e-01,-1.046703e-01,1.046703e-01,-1.255571e-01,1.255571e-01,-1.100135e-01,1.100135e-01,-1.501592e-01,1.501592e-01,-1.155712e-01,1.155712e-01,-1.145563e-01,1.145563e-01,-1.013425e-01,1.013425e-01,-1.145783e-01,1.145783e-01,-1.328031e-01,1.328031e-01,-1.077413e-01,1.077413e-01,-1.064996e-01,1.064996e-01,-1.191170e-01,1.191170e-01,-1.213217e-01,1.213217e-01,-1.260969e-01,1.260969e-01,-1.156494e-01,1.156494e-01,-1.268126e-01,1.268126e-01,-1.070999e-01,1.070999e-01,-1.112365e-01,1.112365e-01,-1.243916e-01,1.243916e-01,-1.283152e-01,1.283152e-01,-1.166925e-01,1.166925e-01,-8.997633e-02,8.997633e-02,-1.583840e-01,1.583840e-01,-1.211178e-01,1.211178e-01,-1.090830e-01,1.090830e-01,-1.030818e-01,1.030818e-01,-1.440600e-01,1.440600e-01,-1.458713e-01,1.458713e-01,-1.559082e-01,1.559082e-01,-1.058868e-01,1.058868e-01,-1.010130e-01,1.010130e-01,-1.642301e-01,1.642301e-01,-1.236850e-01,1.236850e-01,-1.467589e-01,1.467589e-01,-1.109359e-01,1.109359e-01,-1.673655e-01,1.673655e-01,-1.239984e-01,1.239984e-01,-1.039509e-01,1.039509e-01,-1.089378e-01,1.089378e-01,-1.545085e-01,1.545085e-01,-1.200862e-01,1.200862e-01,-1.105608e-01,1.105608e-01,-1.235262e-01,1.235262e-01,-8.496153e-02,8.496153e-02,-1.181372e-01,1.181372e-01,-1.139467e-01,1.139467e-01,-1.189317e-01,1.189317e-01,-1.266519e-01,1.266519e-01,-9.470736e-02,9.470736e-02,-1.336735e-01,1.336735e-01,-8.726601e-02,8.726601e-02,-1.304782e-01,1.304782e-01,-1.186529e-01,1.186529e-01,-1.355944e-01,1.355944e-01,-9.568801e-02,9.568801e-02,-1.282618e-01,1.282618e-01,-1.625632e-01,1.625632e-01,-1.167652e-01,1.167652e-01,-1.001301e-01,1.001301e-01,-1.292419e-01,1.292419e-01,-1.904213e-01,1.904213e-01,-1.511542e-01,1.511542e-01,-9.814394e-02,9.814394e-02,-1.171564e-01,1.171564e-01,-9.806486e-02,9.806486e-02,-9.217615e-02,9.217615e-02,-8.505645e-02,8.505645e-02,-1.573637e-01,1.573637e-01,-1.419174e-01,1.419174e-01,-1.298601e-01,1.298601e-01,-1.120613e-01,1.120613e-01,-1.158363e-01,1.158363e-01,-1.090957e-01,1.090957e-01,-1.204516e-01,1.204516e-01,-1.139852e-01,1.139852e-01,-9.642479e-02,9.642479e-02,-1.410872e-01,1.410872e-01,-1.142779e-01,1.142779e-01,-1.043991e-01,1.043991e-01,-9.736463e-02,9.736463e-02,-1.451046e-01,1.451046e-01,-1.205668e-01,1.205668e-01,-9.881445e-02,9.881445e-02,-1.612822e-01,1.612822e-01,-1.175681e-01,1.175681e-01,-1.522528e-01,1.522528e-01,-1.617520e-01,1.617520e-01,-1.582938e-01,1.582938e-01,-1.208202e-01,1.208202e-01,-1.016003e-01,1.016003e-01,-1.232059e-01,1.232059e-01,-9.583025e-02,9.583025e-02,-1.013990e-01,1.013990e-01,-1.178752e-01,1.178752e-01,-1.215972e-01,1.215972e-01,-1.294932e-01,1.294932e-01,-1.158270e-01,1.158270e-01,-1.008645e-01,1.008645e-01,-9.699190e-02,9.699190e-02,-1.022144e-01,1.022144e-01,-9.878768e-02,9.878768e-02,-1.339052e-01,1.339052e-01,-9.279961e-02,9.279961e-02,-1.047606e-01,1.047606e-01,-1.141163e-01,1.141163e-01,-1.267600e-01,1.267600e-01,-1.252763e-01,1.252763e-01,-9.775003e-02,9.775003e-02,-9.169116e-02,9.169116e-02,-1.006496e-01,1.006496e-01,-9.493293e-02,9.493293e-02,-1.213694e-01,1.213694e-01,-1.109243e-01,1.109243e-01,-1.115973e-01,1.115973e-01,-7.979327e-02,7.979327e-02,-9.220953e-02,9.220953e-02,-1.028913e-01,1.028913e-01,-1.253510e-01,1.253510e-01]},{"count":391,"threshold":-4.665692e+00,"feature":[{"size":5,"px":[14,9,11,17,12],"py":[2,3,9,13,3],"pz":[0,0,0,0,0],"nx":[21,8,7,20,13],"ny":[16,10,7,7,9],"nz":[0,1,1,0,0]},{"size":5,"px":[12,10,6,11,13],"py":[9,3,13,3,4],"pz":[0,0,0,0,0],"nx":[10,4,5,10,2],"ny":[9,10,8,8,2],"nz":[0,1,1,0,2]},{"size":5,"px":[6,9,7,8,8],"py":[3,3,3,3,3],"pz":[0,0,0,0,-1],"nx":[0,0,0,4,9],"ny":[4,2,3,10,8],"nz":[0,0,0,1,0]},{"size":5,"px":[6,2,16,6,8],"py":[16,2,11,4,11],"pz":[0,2,0,1,0],"nx":[3,8,4,1,1],"ny":[4,4,4,5,13],"nz":[1,1,-1,-1,-1]},{"size":3,"px":[16,13,9],"py":[23,18,10],"pz":[0,0,1],"nx":[14,15,8],"ny":[21,22,3],"nz":[0,-1,-1]},{"size":5,"px":[9,16,19,17,17],"py":[1,2,3,2,2],"pz":[1,0,0,0,-1],"nx":[23,23,23,23,23],"ny":[6,2,1,3,5],"nz":[0,0,0,0,0]},{"size":5,"px":[12,12,12,12,12],"py":[10,11,12,13,13],"pz":[0,0,0,0,-1],"nx":[4,8,14,4,6],"ny":[2,4,7,4,8],"nz":[2,1,0,1,1]},{"size":5,"px":[1,2,3,6,4],"py":[6,10,12,23,13],"pz":[1,1,0,0,0],"nx":[2,0,0,1,1],"ny":[23,5,10,21,21],"nz":[0,2,1,0,-1]},{"size":5,"px":[12,16,12,4,12],"py":[6,17,7,2,8],"pz":[0,0,0,2,0],"nx":[8,8,12,0,6],"ny":[4,4,16,0,8],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[9,2],"py":[18,4],"pz":[0,-1],"nx":[4,9],"ny":[10,16],"nz":[1,0]},{"size":5,"px":[9,9,2,0,12],"py":[6,6,21,4,8],"pz":[1,-1,-1,-1,-1],"nx":[8,4,9,7,7],"ny":[10,2,4,5,8],"nz":[1,2,1,1,1]},{"size":5,"px":[10,10,10,18,19],"py":[10,8,7,14,14],"pz":[1,1,1,0,0],"nx":[21,23,22,22,11],"ny":[23,19,21,22,10],"nz":[0,0,0,0,-1]},{"size":5,"px":[12,3,15,4,19],"py":[14,0,5,5,14],"pz":[0,-1,-1,-1,-1],"nx":[12,17,15,3,8],"ny":[18,18,14,2,10],"nz":[0,0,0,2,0]},{"size":5,"px":[8,11,3,11,4],"py":[23,7,9,8,8],"pz":[0,0,1,0,1],"nx":[8,0,10,0,8],"ny":[8,2,8,4,10],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[10,11,12,8,4],"py":[3,0,0,1,1],"pz":[0,0,0,0,1],"nx":[2,3,4,3,3],"ny":[14,5,0,1,2],"nz":[0,0,0,0,0]},{"size":2,"px":[3,11],"py":[7,0],"pz":[1,-1],"nx":[5,2],"ny":[9,5],"nz":[1,2]},{"size":5,"px":[7,1,0,10,1],"py":[0,0,2,12,6],"pz":[0,2,2,0,1],"nx":[4,6,2,8,8],"ny":[4,11,2,4,4],"nz":[1,1,2,1,-1]},{"size":2,"px":[4,15],"py":[4,12],"pz":[2,0],"nx":[4,6],"ny":[5,11],"nz":[2,-1]},{"size":5,"px":[9,4,16,14,14],"py":[8,4,23,18,18],"pz":[1,2,0,0,-1],"nx":[0,2,1,1,0],"ny":[2,0,3,2,3],"nz":[1,0,0,0,1]},{"size":5,"px":[17,7,7,18,19],"py":[7,11,8,7,7],"pz":[0,1,1,0,0],"nx":[17,5,8,2,0],"ny":[8,0,7,5,3],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[5,14],"py":[12,3],"pz":[0,-1],"nx":[4,3],"ny":[5,4],"nz":[1,1]},{"size":5,"px":[10,8,16,11,11],"py":[5,6,12,4,4],"pz":[0,1,0,0,-1],"nx":[14,13,5,9,5],"ny":[13,10,1,4,2],"nz":[0,0,2,1,2]},{"size":5,"px":[15,14,16,8,8],"py":[2,2,2,0,0],"pz":[0,0,0,1,-1],"nx":[9,18,19,18,17],"ny":[0,0,2,1,0],"nz":[1,0,0,0,0]},{"size":2,"px":[17,15],"py":[12,11],"pz":[0,0],"nx":[14,4],"ny":[9,15],"nz":[0,-1]},{"size":3,"px":[5,11,11],"py":[3,4,5],"pz":[2,1,1],"nx":[14,3,18],"ny":[6,5,0],"nz":[0,1,-1]},{"size":5,"px":[16,14,17,15,9],"py":[2,2,2,2,1],"pz":[0,0,0,0,1],"nx":[21,20,11,21,21],"ny":[2,0,7,3,3],"nz":[0,0,1,0,-1]},{"size":5,"px":[2,1,1,1,5],"py":[12,9,7,3,6],"pz":[0,0,1,1,1],"nx":[4,8,3,4,17],"ny":[4,4,0,8,0],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[8,4],"py":[6,3],"pz":[1,2],"nx":[9,2],"ny":[4,17],"nz":[1,-1]},{"size":2,"px":[8,5],"py":[16,9],"pz":[0,1],"nx":[10,17],"ny":[16,8],"nz":[0,-1]},{"size":4,"px":[11,5,9,15],"py":[14,9,11,5],"pz":[0,-1,-1,-1],"nx":[10,1,9,4],"ny":[9,2,13,7],"nz":[0,2,0,1]},{"size":5,"px":[2,5,10,7,10],"py":[7,12,2,13,3],"pz":[1,-1,-1,-1,-1],"nx":[5,2,3,3,2],"ny":[23,15,17,16,14],"nz":[0,0,0,0,0]},{"size":2,"px":[11,7],"py":[8,10],"pz":[0,-1],"nx":[7,14],"ny":[5,8],"nz":[1,0]},{"size":2,"px":[9,16],"py":[7,23],"pz":[1,0],"nx":[4,4],"ny":[2,1],"nz":[2,-1]},{"size":5,"px":[16,14,18,4,17],"py":[0,0,4,0,1],"pz":[0,0,0,2,0],"nx":[8,8,16,9,9],"ny":[5,4,11,7,7],"nz":[1,1,0,0,-1]},{"size":5,"px":[12,13,7,8,4],"py":[9,12,6,11,5],"pz":[0,0,1,1,2],"nx":[23,23,16,9,9],"ny":[0,1,11,7,7],"nz":[0,-1,-1,-1,-1]},{"size":3,"px":[6,7,2],"py":[21,23,4],"pz":[0,0,2],"nx":[4,1,16],"ny":[10,5,11],"nz":[1,-1,-1]},{"size":2,"px":[2,2],"py":[3,4],"pz":[2,2],"nx":[3,1],"ny":[4,5],"nz":[1,-1]},{"size":5,"px":[1,2,1,0,1],"py":[7,13,12,4,13],"pz":[0,0,0,2,0],"nx":[18,9,9,19,19],"ny":[23,5,11,19,19],"nz":[0,1,1,0,-1]},{"size":3,"px":[4,10,12],"py":[6,2,5],"pz":[1,-1,-1],"nx":[10,0,0],"ny":[12,1,3],"nz":[0,2,2]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,0],"ny":[4,3],"nz":[1,-1]},{"size":5,"px":[19,17,10,14,18],"py":[2,1,7,0,1],"pz":[0,0,1,0,0],"nx":[3,3,3,7,5],"ny":[9,10,7,23,18],"nz":[1,1,1,0,0]},{"size":2,"px":[10,10],"py":[8,7],"pz":[1,1],"nx":[14,4],"ny":[15,6],"nz":[0,-1]},{"size":2,"px":[7,15],"py":[1,3],"pz":[1,0],"nx":[16,19],"ny":[1,3],"nz":[0,-1]},{"size":5,"px":[11,11,1,2,11],"py":[11,12,1,13,12],"pz":[0,0,-1,-1,-1],"nx":[12,17,8,16,8],"ny":[7,12,11,16,6],"nz":[0,0,0,0,1]},{"size":5,"px":[13,11,10,12,5],"py":[0,0,0,0,0],"pz":[0,0,0,0,1],"nx":[8,4,3,4,4],"ny":[4,5,2,4,4],"nz":[1,1,2,1,-1]},{"size":5,"px":[6,1,3,2,3],"py":[13,3,3,4,10],"pz":[0,2,1,1,1],"nx":[0,1,0,0,0],"ny":[2,0,5,4,4],"nz":[0,0,0,0,-1]},{"size":2,"px":[15,1],"py":[4,3],"pz":[0,-1],"nx":[16,15],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[3,7],"py":[7,13],"pz":[1,0],"nx":[3,0],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[14,15],"py":[18,14],"pz":[0,-1],"nx":[4,14],"ny":[4,16],"nz":[1,0]},{"size":2,"px":[4,6],"py":[3,4],"pz":[2,1],"nx":[9,5],"ny":[14,2],"nz":[0,-1]},{"size":2,"px":[16,6],"py":[1,5],"pz":[0,-1],"nx":[4,9],"ny":[0,4],"nz":[2,1]},{"size":2,"px":[9,0],"py":[4,2],"pz":[0,-1],"nx":[5,3],"ny":[1,0],"nz":[1,2]},{"size":5,"px":[1,1,1,0,0],"py":[16,15,17,6,9],"pz":[0,0,0,1,0],"nx":[9,5,4,9,8],"ny":[7,3,3,6,7],"nz":[0,1,1,0,-1]},{"size":2,"px":[9,1],"py":[8,15],"pz":[1,-1],"nx":[9,8],"ny":[9,4],"nz":[1,1]},{"size":2,"px":[20,19],"py":[19,22],"pz":[0,0],"nx":[7,0],"ny":[3,0],"nz":[1,-1]},{"size":5,"px":[8,4,2,5,5],"py":[12,6,3,5,5],"pz":[0,1,2,1,-1],"nx":[22,21,20,21,22],"ny":[17,20,22,19,16],"nz":[0,0,0,0,0]},{"size":2,"px":[6,12],"py":[2,6],"pz":[1,0],"nx":[8,3],"ny":[3,2],"nz":[1,-1]},{"size":2,"px":[11,11],"py":[9,4],"pz":[1,1],"nx":[12,4],"ny":[17,5],"nz":[0,-1]},{"size":3,"px":[0,1,0],"py":[5,13,3],"pz":[2,0,2],"nx":[0,4,11],"ny":[23,5,1],"nz":[0,-1,-1]},{"size":2,"px":[10,5],"py":[6,3],"pz":[0,1],"nx":[4,4],"ny":[3,0],"nz":[1,-1]},{"size":2,"px":[6,5],"py":[7,3],"pz":[0,-1],"nx":[0,1],"ny":[4,10],"nz":[2,1]},{"size":5,"px":[12,13,12,12,12],"py":[12,13,11,10,10],"pz":[0,0,0,0,-1],"nx":[10,8,8,16,15],"ny":[7,4,10,11,10],"nz":[0,1,0,0,0]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[4,2],"ny":[5,5],"nz":[2,-1]},{"size":2,"px":[9,17],"py":[17,7],"pz":[0,-1],"nx":[5,2],"ny":[9,4],"nz":[1,2]},{"size":2,"px":[4,4],"py":[3,5],"pz":[2,2],"nx":[12,8],"ny":[16,2],"nz":[0,-1]},{"size":2,"px":[1,1],"py":[2,0],"pz":[1,1],"nx":[0,4],"ny":[0,1],"nz":[2,-1]},{"size":2,"px":[11,1],"py":[5,0],"pz":[0,-1],"nx":[2,3],"ny":[2,4],"nz":[2,1]},{"size":4,"px":[0,6,4,22],"py":[23,2,4,12],"pz":[0,-1,-1,-1],"nx":[7,6,8,5],"ny":[1,1,2,1],"nz":[1,1,1,1]},{"size":2,"px":[4,10],"py":[0,9],"pz":[1,-1],"nx":[2,4],"ny":[3,10],"nz":[2,1]},{"size":2,"px":[11,8],"py":[15,13],"pz":[0,-1],"nx":[23,11],"ny":[13,5],"nz":[0,1]},{"size":2,"px":[18,4],"py":[5,4],"pz":[0,-1],"nx":[18,20],"ny":[4,7],"nz":[0,0]},{"size":5,"px":[21,20,20,10,20],"py":[17,22,19,10,21],"pz":[0,0,0,1,0],"nx":[5,5,3,14,7],"ny":[9,9,0,8,4],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[3,7,13,7,3],"py":[6,12,3,0,3],"pz":[1,-1,-1,-1,-1],"nx":[1,5,0,0,2],"ny":[16,6,13,5,4],"nz":[0,1,0,1,0]},{"size":2,"px":[7,4],"py":[6,3],"pz":[1,2],"nx":[9,5],"ny":[4,6],"nz":[1,-1]},{"size":3,"px":[14,9,13],"py":[19,22,8],"pz":[0,-1,-1],"nx":[13,4,4],"ny":[17,2,5],"nz":[0,2,2]},{"size":2,"px":[16,4],"py":[9,3],"pz":[0,2],"nx":[7,4],"ny":[4,5],"nz":[1,-1]},{"size":4,"px":[10,2,4,2],"py":[23,4,8,3],"pz":[0,2,1,2],"nx":[14,0,4,11],"ny":[19,3,5,3],"nz":[0,-1,-1,-1]},{"size":5,"px":[9,10,8,7,11],"py":[2,2,2,2,2],"pz":[0,0,0,0,0],"nx":[6,5,3,4,4],"ny":[0,1,0,2,2],"nz":[0,0,1,0,-1]},{"size":2,"px":[6,4],"py":[13,6],"pz":[0,-1],"nx":[15,4],"ny":[8,4],"nz":[0,1]},{"size":2,"px":[0,8],"py":[1,2],"pz":[2,-1],"nx":[5,4],"ny":[2,2],"nz":[1,1]},{"size":5,"px":[16,13,14,15,15],"py":[1,0,0,0,0],"pz":[0,0,0,0,-1],"nx":[4,9,4,18,8],"ny":[5,9,4,18,11],"nz":[2,1,2,0,1]},{"size":2,"px":[5,6],"py":[2,6],"pz":[2,1],"nx":[22,9],"ny":[23,9],"nz":[0,-1]},{"size":2,"px":[19,19],"py":[5,5],"pz":[0,-1],"nx":[21,22],"ny":[2,4],"nz":[0,0]},{"size":2,"px":[2,5],"py":[8,6],"pz":[0,1],"nx":[3,4],"ny":[4,9],"nz":[1,-1]},{"size":2,"px":[18,14],"py":[13,17],"pz":[0,0],"nx":[14,4],"ny":[16,3],"nz":[0,-1]},{"size":2,"px":[6,6],"py":[6,3],"pz":[1,-1],"nx":[1,0],"ny":[2,2],"nz":[1,2]},{"size":2,"px":[23,21],"py":[21,14],"pz":[0,-1],"nx":[7,5],"ny":[0,0],"nz":[0,1]},{"size":2,"px":[15,10],"py":[23,7],"pz":[0,-1],"nx":[9,4],"ny":[4,5],"nz":[1,2]},{"size":2,"px":[4,18],"py":[3,8],"pz":[2,0],"nx":[8,4],"ny":[4,5],"nz":[1,-1]},{"size":2,"px":[13,7],"py":[2,11],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":5,"px":[2,3,5,6,1],"py":[7,14,2,2,4],"pz":[1,0,0,0,2],"nx":[8,4,4,7,7],"ny":[7,5,4,9,9],"nz":[1,2,2,1,-1]},{"size":2,"px":[5,3],"py":[6,3],"pz":[1,-1],"nx":[1,2],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[7,20,4,10,10],"py":[9,16,4,10,8],"pz":[1,0,2,1,1],"nx":[4,2,3,5,3],"ny":[11,5,6,12,5],"nz":[0,1,1,0,-1]},{"size":2,"px":[6,11],"py":[4,18],"pz":[1,-1],"nx":[8,6],"ny":[4,9],"nz":[1,1]},{"size":2,"px":[2,8],"py":[5,23],"pz":[2,0],"nx":[9,4],"ny":[0,2],"nz":[1,-1]},{"size":5,"px":[3,1,2,2,2],"py":[12,6,12,11,11],"pz":[0,1,0,0,-1],"nx":[0,0,0,0,0],"ny":[13,12,11,14,7],"nz":[0,0,0,0,1]},{"size":2,"px":[3,6],"py":[1,2],"pz":[2,1],"nx":[8,4],"ny":[4,14],"nz":[1,-1]},{"size":5,"px":[11,23,23,22,22],"py":[8,12,6,13,14],"pz":[1,0,0,0,0],"nx":[13,8,7,6,6],"ny":[6,3,3,9,9],"nz":[0,1,1,0,-1]},{"size":4,"px":[9,23,23,22],"py":[7,12,6,13],"pz":[1,-1,-1,-1],"nx":[11,23,23,23],"ny":[6,13,17,10],"nz":[1,0,0,0]},{"size":5,"px":[0,0,0,0,0],"py":[19,5,9,16,10],"pz":[0,2,1,0,1],"nx":[5,2,1,2,2],"ny":[18,10,5,9,9],"nz":[0,1,2,1,-1]},{"size":2,"px":[11,5],"py":[10,4],"pz":[1,2],"nx":[23,14],"ny":[23,3],"nz":[0,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,1],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[8,10],"py":[4,8],"pz":[0,-1],"nx":[8,8],"ny":[2,3],"nz":[0,0]},{"size":3,"px":[7,10,11],"py":[1,6,13],"pz":[0,-1,-1],"nx":[4,4,2],"ny":[3,8,2],"nz":[1,1,2]},{"size":2,"px":[8,4],"py":[8,2],"pz":[1,2],"nx":[10,5],"ny":[10,0],"nz":[0,-1]},{"size":2,"px":[7,16],"py":[20,21],"pz":[0,-1],"nx":[2,4],"ny":[5,10],"nz":[2,1]},{"size":2,"px":[3,10],"py":[7,8],"pz":[1,-1],"nx":[7,4],"ny":[20,7],"nz":[0,1]},{"size":5,"px":[11,11,11,11,11],"py":[10,12,13,11,11],"pz":[0,0,0,0,-1],"nx":[11,12,16,3,8],"ny":[6,6,10,1,8],"nz":[0,0,0,2,0]},{"size":2,"px":[12,6],"py":[4,2],"pz":[0,1],"nx":[7,7],"ny":[8,1],"nz":[0,-1]},{"size":5,"px":[23,23,23,23,23],"py":[22,20,21,19,19],"pz":[0,0,0,0,-1],"nx":[4,6,3,4,3],"ny":[19,23,15,20,16],"nz":[0,0,0,0,0]},{"size":3,"px":[8,4,14],"py":[12,3,8],"pz":[0,-1,-1],"nx":[4,2,10],"ny":[10,3,13],"nz":[1,2,0]},{"size":2,"px":[11,18],"py":[13,23],"pz":[0,-1],"nx":[5,5],"ny":[1,2],"nz":[2,2]},{"size":3,"px":[11,2,10],"py":[17,4,17],"pz":[0,2,0],"nx":[11,0,22],"ny":[15,2,4],"nz":[0,-1,-1]},{"size":3,"px":[11,3,0],"py":[15,4,8],"pz":[0,-1,-1],"nx":[14,11,4],"ny":[9,17,7],"nz":[0,0,1]},{"size":2,"px":[17,16],"py":[2,1],"pz":[0,0],"nx":[9,11],"ny":[4,6],"nz":[1,-1]},{"size":2,"px":[3,4],"py":[21,23],"pz":[0,0],"nx":[4,0],"ny":[3,3],"nz":[1,-1]},{"size":2,"px":[18,2],"py":[20,0],"pz":[0,-1],"nx":[4,9],"ny":[5,10],"nz":[2,1]},{"size":2,"px":[9,1],"py":[19,3],"pz":[0,-1],"nx":[0,0],"ny":[9,21],"nz":[1,0]},{"size":2,"px":[19,19],"py":[21,22],"pz":[0,0],"nx":[19,0],"ny":[23,0],"nz":[0,-1]},{"size":4,"px":[11,2,3,2],"py":[6,6,9,4],"pz":[0,-1,-1,-1],"nx":[4,9,19,19],"ny":[5,10,17,18],"nz":[2,1,0,0]},{"size":2,"px":[2,4],"py":[4,8],"pz":[2,1],"nx":[4,9],"ny":[10,10],"nz":[1,-1]},{"size":2,"px":[23,22],"py":[8,12],"pz":[0,-1],"nx":[7,4],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[12,1],"py":[5,2],"pz":[0,-1],"nx":[9,11],"ny":[2,1],"nz":[0,0]},{"size":2,"px":[4,4],"py":[2,2],"pz":[0,-1],"nx":[3,2],"ny":[1,2],"nz":[0,0]},{"size":2,"px":[17,9],"py":[13,7],"pz":[0,1],"nx":[9,5],"ny":[4,0],"nz":[1,-1]},{"size":4,"px":[0,0,9,13],"py":[3,3,7,3],"pz":[2,-1,-1,-1],"nx":[2,4,4,11],"ny":[1,2,8,5],"nz":[2,1,1,0]},{"size":5,"px":[3,6,5,6,6],"py":[0,0,2,1,1],"pz":[1,0,0,0,-1],"nx":[2,2,2,1,1],"ny":[21,19,20,16,17],"nz":[0,0,0,0,0]},{"size":2,"px":[13,3],"py":[22,10],"pz":[0,-1],"nx":[7,4],"ny":[10,5],"nz":[1,2]},{"size":2,"px":[3,2],"py":[7,3],"pz":[1,2],"nx":[8,4],"ny":[4,5],"nz":[1,-1]},{"size":5,"px":[17,8,15,7,15],"py":[13,6,16,5,12],"pz":[0,1,0,1,0],"nx":[5,4,6,3,4],"ny":[1,2,1,0,3],"nz":[0,0,0,1,-1]},{"size":5,"px":[12,9,11,12,10],"py":[0,1,2,2,0],"pz":[0,0,0,0,0],"nx":[8,16,7,4,4],"ny":[9,23,9,3,2],"nz":[1,0,1,2,-1]},{"size":2,"px":[4,11],"py":[1,4],"pz":[2,-1],"nx":[8,7],"ny":[4,4],"nz":[0,0]},{"size":4,"px":[7,4,5,8],"py":[13,2,1,3],"pz":[0,-1,-1,-1],"nx":[9,4,9,9],"ny":[9,5,10,11],"nz":[0,1,0,0]},{"size":2,"px":[10,11],"py":[10,11],"pz":[0,-1],"nx":[2,6],"ny":[2,2],"nz":[2,1]},{"size":2,"px":[21,3],"py":[11,2],"pz":[0,-1],"nx":[22,22],"ny":[20,18],"nz":[0,0]},{"size":2,"px":[7,6],"py":[1,2],"pz":[0,0],"nx":[5,10],"ny":[1,0],"nz":[0,-1]},{"size":2,"px":[21,3],"py":[18,1],"pz":[0,-1],"nx":[16,15],"ny":[4,4],"nz":[0,0]},{"size":2,"px":[12,7],"py":[4,1],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[13,11],"py":[23,17],"pz":[0,0],"nx":[11,21],"ny":[16,0],"nz":[0,-1]},{"size":2,"px":[1,2],"py":[0,6],"pz":[1,-1],"nx":[16,16],"ny":[9,11],"nz":[0,0]},{"size":2,"px":[12,13],"py":[20,20],"pz":[0,0],"nx":[11,3],"ny":[21,7],"nz":[0,-1]},{"size":3,"px":[19,20,9],"py":[21,18,11],"pz":[0,0,1],"nx":[17,4,11],"ny":[19,2,0],"nz":[0,-1,-1]},{"size":2,"px":[12,5],"py":[5,2],"pz":[0,1],"nx":[7,9],"ny":[7,8],"nz":[0,-1]},{"size":5,"px":[8,4,4,8,4],"py":[4,4,5,10,3],"pz":[1,1,2,0,2],"nx":[11,22,11,23,23],"ny":[0,0,1,3,3],"nz":[1,0,1,0,-1]},{"size":2,"px":[8,14],"py":[10,23],"pz":[1,0],"nx":[7,2],"ny":[10,9],"nz":[1,-1]},{"size":2,"px":[5,14],"py":[6,23],"pz":[1,-1],"nx":[1,2],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[11,2],"py":[19,3],"pz":[0,-1],"nx":[10,12],"ny":[18,18],"nz":[0,0]},{"size":2,"px":[12,3],"py":[4,1],"pz":[0,2],"nx":[6,6],"ny":[11,11],"nz":[1,-1]},{"size":5,"px":[0,0,0,0,0],"py":[18,10,20,19,19],"pz":[0,1,0,0,-1],"nx":[11,10,14,12,13],"ny":[2,2,2,2,2],"nz":[0,0,0,0,0]},{"size":3,"px":[12,2,9],"py":[14,5,10],"pz":[0,-1,-1],"nx":[11,10,5],"ny":[10,13,5],"nz":[0,0,1]},{"size":2,"px":[2,3],"py":[3,7],"pz":[2,1],"nx":[3,10],"ny":[4,13],"nz":[1,-1]},{"size":2,"px":[9,3],"py":[21,7],"pz":[0,-1],"nx":[10,21],"ny":[7,15],"nz":[1,0]},{"size":2,"px":[21,10],"py":[16,8],"pz":[0,1],"nx":[8,2],"ny":[10,8],"nz":[1,-1]},{"size":2,"px":[8,8],"py":[6,7],"pz":[1,-1],"nx":[12,11],"ny":[11,7],"nz":[0,1]},{"size":2,"px":[3,11],"py":[4,20],"pz":[2,0],"nx":[11,10],"ny":[19,1],"nz":[0,-1]},{"size":2,"px":[17,5],"py":[13,3],"pz":[0,-1],"nx":[7,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[7,1],"py":[23,3],"pz":[0,2],"nx":[14,6],"ny":[12,9],"nz":[0,-1]},{"size":2,"px":[12,5],"py":[11,2],"pz":[0,-1],"nx":[11,7],"ny":[3,1],"nz":[0,1]},{"size":2,"px":[9,6],"py":[2,17],"pz":[0,-1],"nx":[4,6],"ny":[4,12],"nz":[1,0]},{"size":2,"px":[14,19],"py":[5,6],"pz":[0,-1],"nx":[9,3],"ny":[9,1],"nz":[0,2]},{"size":5,"px":[12,13,13,13,12],"py":[9,11,12,13,10],"pz":[0,0,0,0,0],"nx":[2,4,4,4,4],"ny":[7,18,17,14,14],"nz":[1,0,0,0,-1]},{"size":2,"px":[10,10],"py":[6,6],"pz":[1,-1],"nx":[20,18],"ny":[18,23],"nz":[0,0]},{"size":2,"px":[5,6],"py":[4,14],"pz":[1,-1],"nx":[9,4],"ny":[2,1],"nz":[1,2]},{"size":2,"px":[11,9],"py":[4,18],"pz":[0,-1],"nx":[4,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[15,0],"py":[18,4],"pz":[0,-1],"nx":[3,4],"ny":[5,4],"nz":[2,2]},{"size":4,"px":[7,3,6,6],"py":[8,4,6,5],"pz":[1,2,1,1],"nx":[10,4,13,0],"ny":[10,4,9,22],"nz":[0,-1,-1,-1]},{"size":2,"px":[10,8],"py":[18,11],"pz":[0,-1],"nx":[5,4],"ny":[8,10],"nz":[1,1]},{"size":4,"px":[17,2,10,2],"py":[14,1,10,3],"pz":[0,-1,-1,-1],"nx":[8,8,17,8],"ny":[4,5,12,6],"nz":[1,1,0,1]},{"size":5,"px":[9,11,9,4,10],"py":[1,1,0,0,1],"pz":[0,0,0,1,0],"nx":[8,4,7,15,15],"ny":[7,2,4,17,17],"nz":[1,2,1,0,-1]},{"size":2,"px":[4,3],"py":[11,8],"pz":[0,-1],"nx":[2,2],"ny":[1,2],"nz":[2,2]},{"size":2,"px":[11,3],"py":[13,8],"pz":[0,-1],"nx":[1,1],"ny":[5,2],"nz":[1,2]},{"size":2,"px":[6,2],"py":[8,3],"pz":[0,2],"nx":[3,1],"ny":[5,2],"nz":[1,-1]},{"size":5,"px":[10,5,7,8,6],"py":[9,7,7,7,7],"pz":[0,0,0,0,0],"nx":[7,3,0,2,15],"ny":[8,0,1,18,17],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[17,8],"py":[12,6],"pz":[0,1],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":5,"px":[3,11,8,10,12],"py":[0,2,10,2,3],"pz":[2,0,0,0,0],"nx":[3,2,10,2,2],"ny":[6,4,11,3,3],"nz":[0,1,0,1,-1]},{"size":2,"px":[3,6],"py":[2,4],"pz":[2,1],"nx":[8,19],"ny":[4,16],"nz":[1,-1]},{"size":2,"px":[2,2],"py":[1,1],"pz":[2,-1],"nx":[7,17],"ny":[1,2],"nz":[1,0]},{"size":5,"px":[16,15,14,13,7],"py":[0,0,0,0,0],"pz":[0,0,0,0,-1],"nx":[6,4,8,3,11],"ny":[3,4,4,1,6],"nz":[1,1,1,2,0]},{"size":2,"px":[11,1],"py":[8,5],"pz":[0,-1],"nx":[13,4],"ny":[10,2],"nz":[0,2]},{"size":2,"px":[4,9],"py":[0,2],"pz":[2,1],"nx":[4,11],"ny":[0,2],"nz":[0,-1]},{"size":2,"px":[15,15],"py":[2,2],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[8,17],"py":[9,22],"pz":[1,0],"nx":[8,20],"ny":[10,2],"nz":[1,-1]},{"size":2,"px":[10,10],"py":[14,22],"pz":[0,-1],"nx":[3,11],"ny":[3,3],"nz":[1,0]},{"size":2,"px":[4,2],"py":[1,0],"pz":[1,2],"nx":[5,8],"ny":[3,9],"nz":[0,-1]},{"size":2,"px":[2,3],"py":[4,8],"pz":[2,1],"nx":[9,5],"ny":[15,19],"nz":[0,-1]},{"size":2,"px":[5,2],"py":[1,1],"pz":[0,1],"nx":[10,10],"ny":[6,6],"nz":[0,-1]},{"size":2,"px":[17,6],"py":[10,2],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":3,"px":[13,7,3],"py":[5,2,6],"pz":[0,1,-1],"nx":[17,16,17],"ny":[1,1,2],"nz":[0,0,0]},{"size":2,"px":[11,10],"py":[3,3],"pz":[0,0],"nx":[8,4],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[4,8],"py":[0,8],"pz":[2,-1],"nx":[3,4],"ny":[0,0],"nz":[1,1]},{"size":5,"px":[9,2,4,1,2],"py":[13,3,9,2,5],"pz":[0,2,1,2,2],"nx":[9,5,10,4,10],"ny":[5,1,3,0,0],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[6,12],"py":[5,9],"pz":[1,0],"nx":[0,2],"ny":[23,9],"nz":[0,-1]},{"size":2,"px":[22,11],"py":[21,8],"pz":[0,1],"nx":[10,0],"ny":[17,2],"nz":[0,-1]},{"size":2,"px":[3,1],"py":[22,9],"pz":[0,1],"nx":[22,5],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[5,6],"ny":[10,9],"nz":[1,-1]},{"size":4,"px":[7,3,17,7],"py":[8,2,10,11],"pz":[0,2,0,1],"nx":[6,10,5,23],"ny":[9,21,1,23],"nz":[0,-1,-1,-1]},{"size":2,"px":[8,3],"py":[7,2],"pz":[1,2],"nx":[8,9],"ny":[4,9],"nz":[1,-1]},{"size":2,"px":[9,5],"py":[14,6],"pz":[0,1],"nx":[8,8],"ny":[13,13],"nz":[0,-1]},{"size":3,"px":[11,6,8],"py":[20,3,20],"pz":[0,-1,-1],"nx":[5,3,12],"ny":[9,5,18],"nz":[1,2,0]},{"size":2,"px":[3,9],"py":[1,3],"pz":[1,0],"nx":[2,8],"ny":[5,8],"nz":[0,-1]},{"size":2,"px":[15,9],"py":[21,3],"pz":[0,-1],"nx":[3,4],"ny":[5,5],"nz":[2,2]},{"size":2,"px":[2,9],"py":[7,11],"pz":[1,-1],"nx":[2,2],"ny":[8,9],"nz":[1,1]},{"size":4,"px":[3,4,3,1],"py":[14,21,19,6],"pz":[0,0,0,1],"nx":[10,16,4,5],"ny":[8,1,7,6],"nz":[0,-1,-1,-1]},{"size":4,"px":[10,4,3,1],"py":[5,21,19,6],"pz":[1,-1,-1,-1],"nx":[21,10,5,11],"ny":[4,2,3,4],"nz":[0,1,2,1]},{"size":2,"px":[4,17],"py":[3,8],"pz":[2,0],"nx":[17,2],"ny":[9,22],"nz":[0,-1]},{"size":2,"px":[17,12],"py":[14,20],"pz":[0,-1],"nx":[7,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[10,12],"py":[9,20],"pz":[0,-1],"nx":[11,23],"ny":[8,18],"nz":[1,0]},{"size":2,"px":[5,11],"py":[4,7],"pz":[2,1],"nx":[8,15],"ny":[7,5],"nz":[1,-1]},{"size":2,"px":[11,15],"py":[13,8],"pz":[0,-1],"nx":[11,11],"ny":[6,7],"nz":[1,1]},{"size":2,"px":[6,15],"py":[14,8],"pz":[0,-1],"nx":[4,4],"ny":[12,13],"nz":[0,0]},{"size":2,"px":[5,5],"py":[0,1],"pz":[2,2],"nx":[15,4],"ny":[5,5],"nz":[0,-1]},{"size":2,"px":[16,17],"py":[2,2],"pz":[0,0],"nx":[20,8],"ny":[3,7],"nz":[0,-1]},{"size":3,"px":[6,3,2],"py":[10,6,1],"pz":[0,-1,-1],"nx":[4,3,2],"ny":[3,4,2],"nz":[1,1,2]},{"size":2,"px":[10,6],"py":[4,6],"pz":[0,-1],"nx":[6,13],"ny":[0,1],"nz":[1,0]},{"size":2,"px":[10,10],"py":[8,7],"pz":[1,1],"nx":[8,2],"ny":[7,2],"nz":[1,-1]},{"size":2,"px":[7,1],"py":[12,4],"pz":[0,-1],"nx":[3,4],"ny":[5,5],"nz":[1,1]},{"size":2,"px":[11,15],"py":[15,14],"pz":[0,-1],"nx":[3,11],"ny":[4,13],"nz":[1,0]},{"size":5,"px":[13,9,11,14,12],"py":[0,2,0,0,2],"pz":[0,0,0,0,0],"nx":[5,4,4,3,4],"ny":[4,4,18,7,17],"nz":[1,1,0,1,0]},{"size":3,"px":[13,12,11],"py":[22,22,22],"pz":[0,0,0],"nx":[11,12,13],"ny":[20,20,20],"nz":[0,0,0]},{"size":2,"px":[6,13],"py":[2,4],"pz":[1,0],"nx":[7,6],"ny":[8,9],"nz":[0,-1]},{"size":2,"px":[0,0],"py":[23,4],"pz":[0,-1],"nx":[5,9],"ny":[1,1],"nz":[1,0]},{"size":2,"px":[14,14],"py":[19,19],"pz":[0,-1],"nx":[11,11],"ny":[10,9],"nz":[1,1]},{"size":2,"px":[23,23],"py":[11,9],"pz":[0,0],"nx":[23,23],"ny":[0,11],"nz":[0,-1]},{"size":2,"px":[23,3],"py":[23,5],"pz":[0,-1],"nx":[4,1],"ny":[23,10],"nz":[0,1]},{"size":2,"px":[9,1],"py":[7,4],"pz":[1,-1],"nx":[19,10],"ny":[20,9],"nz":[0,1]},{"size":2,"px":[16,1],"py":[9,4],"pz":[0,-1],"nx":[7,8],"ny":[3,3],"nz":[1,1]},{"size":2,"px":[7,6],"py":[13,13],"pz":[0,0],"nx":[4,5],"ny":[4,11],"nz":[1,-1]},{"size":5,"px":[19,20,20,10,10],"py":[0,0,2,0,1],"pz":[0,0,0,1,1],"nx":[7,7,15,4,4],"ny":[4,13,7,4,4],"nz":[1,0,0,1,-1]},{"size":2,"px":[12,23],"py":[6,5],"pz":[0,-1],"nx":[18,18],"ny":[17,16],"nz":[0,0]},{"size":2,"px":[6,3],"py":[9,2],"pz":[1,2],"nx":[14,18],"ny":[9,1],"nz":[0,-1]},{"size":2,"px":[9,13],"py":[16,5],"pz":[0,-1],"nx":[5,4],"ny":[7,9],"nz":[1,1]},{"size":2,"px":[10,10],"py":[8,10],"pz":[1,1],"nx":[4,1],"ny":[5,3],"nz":[2,-1]},{"size":2,"px":[12,11],"py":[13,4],"pz":[0,-1],"nx":[0,0],"ny":[14,15],"nz":[0,0]},{"size":2,"px":[2,1],"py":[20,17],"pz":[0,0],"nx":[12,12],"ny":[22,2],"nz":[0,-1]},{"size":2,"px":[2,3],"py":[6,7],"pz":[1,-1],"nx":[21,21],"ny":[13,12],"nz":[0,0]},{"size":2,"px":[3,10],"py":[4,23],"pz":[2,0],"nx":[10,2],"ny":[21,5],"nz":[0,-1]},{"size":2,"px":[6,12],"py":[3,6],"pz":[1,0],"nx":[11,0],"ny":[17,1],"nz":[0,-1]},{"size":2,"px":[11,4],"py":[21,9],"pz":[0,-1],"nx":[2,3],"ny":[18,22],"nz":[0,0]},{"size":2,"px":[13,5],"py":[18,9],"pz":[0,-1],"nx":[6,7],"ny":[8,9],"nz":[1,1]},{"size":2,"px":[21,4],"py":[16,3],"pz":[0,-1],"nx":[23,23],"ny":[16,15],"nz":[0,0]},{"size":2,"px":[2,0],"py":[7,4],"pz":[1,-1],"nx":[3,8],"ny":[7,4],"nz":[1,1]},{"size":2,"px":[15,16],"py":[11,12],"pz":[0,0],"nx":[8,5],"ny":[4,5],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[7,5],"pz":[0,0],"nx":[17,17],"ny":[11,10],"nz":[0,-1]},{"size":5,"px":[8,13,12,3,3],"py":[6,23,23,3,3],"pz":[1,0,0,2,-1],"nx":[0,1,0,0,0],"ny":[2,13,4,5,6],"nz":[2,0,1,1,1]},{"size":2,"px":[0,1],"py":[7,8],"pz":[1,-1],"nx":[0,0],"ny":[1,0],"nz":[2,2]},{"size":2,"px":[2,12],"py":[1,7],"pz":[1,-1],"nx":[0,0],"ny":[12,14],"nz":[0,0]},{"size":2,"px":[5,1],"py":[7,4],"pz":[1,2],"nx":[8,0],"ny":[15,14],"nz":[0,-1]},{"size":2,"px":[7,4],"py":[14,8],"pz":[0,-1],"nx":[2,4],"ny":[1,4],"nz":[2,1]},{"size":2,"px":[5,3],"py":[3,1],"pz":[2,-1],"nx":[9,9],"ny":[5,6],"nz":[1,1]},{"size":2,"px":[4,5],"py":[2,3],"pz":[1,-1],"nx":[11,12],"ny":[23,23],"nz":[0,0]},{"size":2,"px":[10,5],"py":[7,0],"pz":[1,-1],"nx":[22,22],"ny":[19,18],"nz":[0,0]},{"size":3,"px":[10,2,9],"py":[20,9,4],"pz":[0,-1,-1],"nx":[1,10,11],"ny":[2,11,9],"nz":[2,0,0]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[9,3],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[17,6],"py":[7,16],"pz":[0,-1],"nx":[17,17],"ny":[9,6],"nz":[0,0]},{"size":3,"px":[8,1,9],"py":[6,3,4],"pz":[1,-1,-1],"nx":[2,9,2],"ny":[5,13,3],"nz":[2,0,2]},{"size":4,"px":[10,10,9,2],"py":[12,11,2,10],"pz":[0,0,-1,-1],"nx":[6,11,3,13],"ny":[2,4,1,4],"nz":[1,0,2,0]},{"size":2,"px":[3,3],"py":[7,1],"pz":[1,-1],"nx":[4,3],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[0,0],"py":[4,8],"pz":[2,1],"nx":[4,4],"ny":[15,5],"nz":[0,-1]},{"size":2,"px":[5,0],"py":[4,8],"pz":[1,-1],"nx":[13,13],"ny":[9,10],"nz":[0,0]},{"size":2,"px":[6,3],"py":[2,1],"pz":[1,2],"nx":[8,17],"ny":[4,12],"nz":[1,-1]},{"size":2,"px":[15,16],"py":[11,6],"pz":[0,0],"nx":[16,17],"ny":[5,12],"nz":[0,-1]},{"size":2,"px":[13,11],"py":[9,7],"pz":[0,-1],"nx":[0,1],"ny":[9,20],"nz":[1,0]},{"size":3,"px":[16,11,20],"py":[4,7,23],"pz":[0,-1,-1],"nx":[8,9,4],"ny":[4,6,4],"nz":[1,1,2]},{"size":2,"px":[1,1],"py":[18,17],"pz":[0,0],"nx":[9,6],"ny":[7,11],"nz":[0,-1]},{"size":3,"px":[4,4,19],"py":[3,2,9],"pz":[2,2,0],"nx":[2,14,11],"ny":[5,3,9],"nz":[1,-1,-1]},{"size":2,"px":[11,19],"py":[13,9],"pz":[0,-1],"nx":[11,11],"ny":[4,5],"nz":[1,1]},{"size":2,"px":[13,7],"py":[19,2],"pz":[0,-1],"nx":[3,5],"ny":[6,12],"nz":[1,0]},{"size":4,"px":[9,4,4,2],"py":[13,9,8,4],"pz":[0,1,1,2],"nx":[13,0,0,14],"ny":[18,11,6,1],"nz":[0,-1,-1,-1]},{"size":2,"px":[11,15],"py":[8,10],"pz":[0,0],"nx":[14,11],"ny":[9,2],"nz":[0,-1]},{"size":2,"px":[3,2],"py":[8,5],"pz":[1,2],"nx":[4,4],"ny":[10,10],"nz":[1,-1]},{"size":4,"px":[4,6,16,14],"py":[1,1,1,7],"pz":[2,1,0,0],"nx":[10,1,1,2],"ny":[8,5,10,3],"nz":[0,-1,-1,-1]},{"size":4,"px":[2,3,1,2],"py":[3,1,0,2],"pz":[0,0,1,0],"nx":[0,0,0,0],"ny":[1,1,2,0],"nz":[0,1,0,1]},{"size":2,"px":[8,8],"py":[6,7],"pz":[1,1],"nx":[8,0],"ny":[4,1],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[3,0],"pz":[0,1],"nx":[2,2],"ny":[1,16],"nz":[1,-1]},{"size":2,"px":[6,6],"py":[19,18],"pz":[0,0],"nx":[2,10],"ny":[5,8],"nz":[2,-1]},{"size":2,"px":[8,5],"py":[21,11],"pz":[0,-1],"nx":[3,2],"ny":[11,5],"nz":[1,2]},{"size":2,"px":[4,9],"py":[4,7],"pz":[2,1],"nx":[8,7],"ny":[10,4],"nz":[1,-1]},{"size":5,"px":[4,18,19,16,19],"py":[3,12,12,23,13],"pz":[2,0,0,0,0],"nx":[2,8,3,2,2],"ny":[4,23,10,5,5],"nz":[2,0,1,2,-1]},{"size":2,"px":[4,8],"py":[6,11],"pz":[1,0],"nx":[8,3],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[3,12],"py":[4,13],"pz":[2,0],"nx":[10,5],"ny":[15,21],"nz":[0,-1]},{"size":2,"px":[2,9],"py":[4,23],"pz":[2,0],"nx":[19,4],"ny":[9,3],"nz":[0,2]},{"size":2,"px":[3,6],"py":[8,15],"pz":[1,0],"nx":[6,1],"ny":[18,5],"nz":[0,-1]},{"size":2,"px":[9,0],"py":[20,3],"pz":[0,-1],"nx":[2,10],"ny":[5,17],"nz":[2,0]},{"size":3,"px":[10,6,3],"py":[2,7,3],"pz":[0,-1,-1],"nx":[5,4,2],"ny":[9,7,2],"nz":[1,1,2]},{"size":2,"px":[14,6],"py":[12,7],"pz":[0,-1],"nx":[2,10],"ny":[0,1],"nz":[2,0]},{"size":3,"px":[10,5,1],"py":[15,5,4],"pz":[0,-1,-1],"nx":[9,4,18],"ny":[2,0,4],"nz":[1,2,0]},{"size":2,"px":[17,2],"py":[12,6],"pz":[0,-1],"nx":[8,16],"ny":[4,11],"nz":[1,0]},{"size":3,"px":[7,13,4],"py":[0,0,1],"pz":[1,0,-1],"nx":[18,4,4],"ny":[13,2,3],"nz":[0,2,2]},{"size":2,"px":[1,11],"py":[10,6],"pz":[0,-1],"nx":[0,1],"ny":[15,17],"nz":[0,0]},{"size":3,"px":[9,12,8],"py":[8,17,11],"pz":[1,0,1],"nx":[12,0,20],"ny":[16,9,13],"nz":[0,-1,-1]},{"size":2,"px":[11,4],"py":[5,8],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[16,3],"py":[9,8],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[6,3],"py":[11,5],"pz":[1,2],"nx":[11,5],"ny":[21,5],"nz":[0,-1]},{"size":2,"px":[11,13],"py":[1,1],"pz":[0,0],"nx":[4,4],"ny":[5,5],"nz":[1,-1]},{"size":2,"px":[14,4],"py":[4,3],"pz":[0,-1],"nx":[12,10],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[3,6],"py":[2,4],"pz":[2,1],"nx":[9,7],"ny":[9,7],"nz":[0,-1]},{"size":3,"px":[5,6,6],"py":[4,4,4],"pz":[1,-1,-1],"nx":[13,8,7],"ny":[8,3,4],"nz":[0,1,1]},{"size":2,"px":[5,5],"py":[2,11],"pz":[1,1],"nx":[10,11],"ny":[22,22],"nz":[0,0]},{"size":2,"px":[16,9],"py":[13,7],"pz":[0,1],"nx":[8,14],"ny":[4,12],"nz":[1,-1]},{"size":2,"px":[13,5],"py":[13,3],"pz":[0,2],"nx":[16,22],"ny":[13,6],"nz":[0,-1]},{"size":4,"px":[4,4,3,4],"py":[4,3,4,5],"pz":[2,2,2,2],"nx":[21,5,17,7],"ny":[0,2,5,23],"nz":[0,-1,-1,-1]},{"size":2,"px":[4,16],"py":[0,1],"pz":[2,0],"nx":[15,1],"ny":[23,10],"nz":[0,-1]},{"size":2,"px":[4,6],"py":[11,2],"pz":[0,-1],"nx":[15,6],"ny":[2,1],"nz":[0,1]},{"size":2,"px":[6,3],"py":[2,1],"pz":[1,2],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":3,"px":[13,14,5],"py":[9,15,2],"pz":[0,-1,-1],"nx":[11,1,11],"ny":[10,3,11],"nz":[0,1,0]},{"size":2,"px":[5,1],"py":[6,2],"pz":[1,-1],"nx":[1,1],"ny":[2,5],"nz":[2,1]},{"size":2,"px":[11,5],"py":[1,0],"pz":[1,2],"nx":[10,4],"ny":[2,3],"nz":[1,-1]},{"size":2,"px":[11,11],"py":[8,9],"pz":[1,1],"nx":[23,4],"ny":[23,2],"nz":[0,-1]},{"size":2,"px":[5,2],"py":[10,2],"pz":[0,-1],"nx":[18,10],"ny":[0,1],"nz":[0,1]},{"size":2,"px":[20,4],"py":[7,3],"pz":[0,2],"nx":[8,4],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[10,4],"py":[5,4],"pz":[1,-1],"nx":[11,11],"ny":[5,6],"nz":[1,1]},{"size":3,"px":[14,15,16],"py":[0,0,1],"pz":[0,0,0],"nx":[8,5,15],"ny":[7,2,10],"nz":[1,-1,-1]},{"size":2,"px":[2,2],"py":[1,1],"pz":[2,-1],"nx":[17,18],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[13,8],"py":[15,7],"pz":[0,-1],"nx":[9,4],"ny":[5,2],"nz":[0,1]},{"size":2,"px":[4,0],"py":[6,17],"pz":[1,-1],"nx":[3,2],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[14,8],"py":[17,9],"pz":[0,-1],"nx":[7,6],"ny":[8,8],"nz":[1,1]},{"size":2,"px":[10,4],"py":[7,1],"pz":[1,-1],"nx":[15,6],"ny":[14,4],"nz":[0,1]},{"size":2,"px":[3,12],"py":[8,19],"pz":[1,0],"nx":[13,10],"ny":[17,9],"nz":[0,-1]},{"size":2,"px":[7,12],"py":[2,4],"pz":[1,0],"nx":[6,11],"ny":[3,2],"nz":[0,-1]},{"size":4,"px":[2,1,6,1],"py":[10,3,23,8],"pz":[1,2,0,1],"nx":[17,10,23,0],"ny":[9,2,20,3],"nz":[0,-1,-1,-1]},{"size":2,"px":[9,9],"py":[2,8],"pz":[0,-1],"nx":[2,2],"ny":[4,2],"nz":[2,2]},{"size":2,"px":[3,16],"py":[1,6],"pz":[2,0],"nx":[8,4],"ny":[2,5],"nz":[1,-1]},{"size":2,"px":[3,6],"py":[1,2],"pz":[2,1],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[5,6],"py":[3,0],"pz":[2,-1],"nx":[9,5],"ny":[2,1],"nz":[0,1]},{"size":2,"px":[3,16],"py":[5,23],"pz":[1,-1],"nx":[0,0],"ny":[6,3],"nz":[1,2]},{"size":4,"px":[0,0,0,0],"py":[3,2,12,5],"pz":[2,2,0,1],"nx":[2,3,2,13],"ny":[5,5,2,19],"nz":[1,-1,-1,-1]},{"size":2,"px":[11,11],"py":[10,11],"pz":[0,0],"nx":[5,5],"ny":[1,1],"nz":[2,-1]},{"size":2,"px":[5,2],"py":[0,4],"pz":[2,-1],"nx":[2,2],"ny":[10,8],"nz":[1,1]},{"size":4,"px":[16,2,8,4],"py":[14,0,11,5],"pz":[0,-1,-1,-1],"nx":[18,14,7,7],"ny":[13,14,8,6],"nz":[0,0,1,1]},{"size":2,"px":[8,9],"py":[2,2],"pz":[0,0],"nx":[5,14],"ny":[4,14],"nz":[1,-1]},{"size":2,"px":[3,5],"py":[11,20],"pz":[1,0],"nx":[11,4],"ny":[0,2],"nz":[0,-1]},{"size":2,"px":[2,2],"py":[3,4],"pz":[2,2],"nx":[3,4],"ny":[4,2],"nz":[1,-1]},{"size":3,"px":[10,4,3],"py":[5,5,3],"pz":[0,-1,-1],"nx":[11,3,10],"ny":[2,0,2],"nz":[0,2,0]},{"size":2,"px":[15,15],"py":[1,1],"pz":[0,-1],"nx":[7,4],"ny":[5,2],"nz":[1,2]},{"size":4,"px":[9,5,2,6],"py":[22,8,4,19],"pz":[0,1,2,0],"nx":[9,5,0,3],"ny":[20,5,22,4],"nz":[0,-1,-1,-1]},{"size":3,"px":[1,4,10],"py":[3,9,12],"pz":[2,1,0],"nx":[0,10,0],"ny":[0,5,0],"nz":[0,-1,-1]},{"size":2,"px":[1,6],"py":[0,7],"pz":[0,-1],"nx":[20,19],"ny":[14,14],"nz":[0,0]},{"size":2,"px":[13,4],"py":[14,15],"pz":[0,-1],"nx":[2,1],"ny":[5,7],"nz":[0,0]},{"size":2,"px":[17,7],"py":[9,11],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[17,9],"py":[12,6],"pz":[0,1],"nx":[15,10],"ny":[9,8],"nz":[0,-1]},{"size":2,"px":[0,0],"py":[0,1],"pz":[2,2],"nx":[9,7],"ny":[6,17],"nz":[1,-1]},{"size":3,"px":[3,3,15],"py":[3,4,6],"pz":[2,1,0],"nx":[0,2,22],"ny":[5,8,9],"nz":[0,-1,-1]},{"size":4,"px":[15,15,15,1],"py":[12,6,6,1],"pz":[0,-1,-1,-1],"nx":[4,7,13,4],"ny":[4,7,12,2],"nz":[2,1,0,2]},{"size":2,"px":[3,15],"py":[12,6],"pz":[0,-1],"nx":[9,1],"ny":[14,2],"nz":[0,2]},{"size":2,"px":[12,12],"py":[11,12],"pz":[0,0],"nx":[9,5],"ny":[4,4],"nz":[1,-1]},{"size":3,"px":[23,6,7],"py":[23,3,4],"pz":[0,-1,-1],"nx":[19,16,17],"ny":[17,14,15],"nz":[0,0,0]},{"size":2,"px":[9,5],"py":[2,7],"pz":[1,-1],"nx":[11,23],"ny":[10,18],"nz":[1,0]},{"size":3,"px":[0,0,0],"py":[4,9,2],"pz":[1,0,2],"nx":[2,0,0],"ny":[9,2,1],"nz":[0,-1,-1]},{"size":2,"px":[12,0],"py":[11,9],"pz":[0,-1],"nx":[1,0],"ny":[18,5],"nz":[0,2]},{"size":2,"px":[5,4],"py":[10,6],"pz":[0,1],"nx":[10,6],"ny":[10,18],"nz":[0,-1]},{"size":2,"px":[13,12],"py":[13,13],"pz":[0,-1],"nx":[5,11],"ny":[1,3],"nz":[2,1]},{"size":2,"px":[10,19],"py":[5,22],"pz":[1,-1],"nx":[4,12],"ny":[1,5],"nz":[2,0]},{"size":2,"px":[8,6],"py":[0,0],"pz":[0,0],"nx":[3,12],"ny":[0,3],"nz":[0,-1]},{"size":2,"px":[9,6],"py":[7,0],"pz":[1,-1],"nx":[12,12],"ny":[10,11],"nz":[0,0]},{"size":4,"px":[3,1,3,2],"py":[20,9,21,19],"pz":[0,1,0,0],"nx":[20,20,5,12],"ny":[10,15,2,10],"nz":[0,-1,-1,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,1],"ny":[4,6],"nz":[1,-1]},{"size":3,"px":[5,11,11],"py":[1,3,4],"pz":[2,1,1],"nx":[3,3,7],"ny":[5,5,0],"nz":[1,-1,-1]},{"size":3,"px":[8,6,7],"py":[10,5,6],"pz":[1,1,1],"nx":[23,3,7],"ny":[0,5,0],"nz":[0,-1,-1]},{"size":2,"px":[2,7],"py":[2,14],"pz":[1,-1],"nx":[7,3],"ny":[12,4],"nz":[0,1]},{"size":2,"px":[5,3],"py":[6,3],"pz":[1,2],"nx":[13,3],"ny":[12,4],"nz":[0,-1]},{"size":2,"px":[11,18],"py":[11,4],"pz":[0,-1],"nx":[23,11],"ny":[19,10],"nz":[0,1]},{"size":2,"px":[7,2],"py":[12,3],"pz":[0,-1],"nx":[8,4],"ny":[11,5],"nz":[0,1]},{"size":2,"px":[11,11],"py":[0,11],"pz":[1,-1],"nx":[3,3],"ny":[19,18],"nz":[0,0]},{"size":2,"px":[11,1],"py":[11,11],"pz":[1,-1],"nx":[13,15],"ny":[6,5],"nz":[0,0]},{"size":2,"px":[8,8],"py":[9,9],"pz":[0,-1],"nx":[5,11],"ny":[1,3],"nz":[2,1]},{"size":4,"px":[6,4,8,3],"py":[6,2,4,3],"pz":[0,2,1,2],"nx":[7,0,15,8],"ny":[8,8,16,7],"nz":[0,-1,-1,-1]},{"size":2,"px":[4,3],"py":[22,20],"pz":[0,0],"nx":[2,8],"ny":[5,4],"nz":[2,-1]},{"size":2,"px":[12,6],"py":[11,0],"pz":[0,-1],"nx":[0,0],"ny":[3,1],"nz":[1,2]},{"size":2,"px":[0,0],"py":[12,7],"pz":[0,1],"nx":[3,1],"ny":[23,9],"nz":[0,-1]},{"size":2,"px":[7,0],"py":[11,5],"pz":[1,-1],"nx":[0,0],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[8,8],"py":[10,10],"pz":[0,-1],"nx":[4,3],"ny":[5,4],"nz":[2,2]},{"size":2,"px":[13,3],"py":[2,4],"pz":[0,-1],"nx":[4,3],"ny":[3,5],"nz":[2,2]},{"size":2,"px":[1,1],"py":[23,22],"pz":[0,0],"nx":[9,0],"ny":[7,3],"nz":[0,-1]},{"size":2,"px":[1,0],"py":[16,15],"pz":[0,0],"nx":[0,14],"ny":[23,12],"nz":[0,-1]},{"size":2,"px":[13,8],"py":[22,0],"pz":[0,-1],"nx":[5,3],"ny":[0,1],"nz":[1,1]},{"size":2,"px":[13,13],"py":[7,7],"pz":[0,-1],"nx":[3,2],"ny":[17,10],"nz":[0,1]},{"size":2,"px":[20,20],"py":[15,16],"pz":[0,0],"nx":[7,3],"ny":[9,17],"nz":[1,-1]},{"size":5,"px":[10,12,11,13,11],"py":[2,2,1,2,2],"pz":[0,0,0,0,0],"nx":[10,18,21,21,19],"ny":[3,1,13,11,2],"nz":[1,0,0,0,0]},{"size":2,"px":[16,3],"py":[6,1],"pz":[0,2],"nx":[15,18],"ny":[8,1],"nz":[0,-1]},{"size":2,"px":[19,3],"py":[8,1],"pz":[0,-1],"nx":[9,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[10,3],"py":[15,18],"pz":[0,-1],"nx":[3,3],"ny":[0,1],"nz":[2,2]},{"size":2,"px":[3,3],"py":[2,3],"pz":[2,2],"nx":[7,3],"ny":[11,1],"nz":[1,-1]},{"size":2,"px":[11,10],"py":[17,9],"pz":[0,-1],"nx":[11,10],"ny":[15,15],"nz":[0,0]},{"size":2,"px":[5,10],"py":[2,4],"pz":[1,0],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[9,10],"py":[3,4],"pz":[0,-1],"nx":[9,10],"ny":[2,1],"nz":[0,0]},{"size":2,"px":[23,11],"py":[13,10],"pz":[0,1],"nx":[14,7],"ny":[5,14],"nz":[0,-1]},{"size":2,"px":[4,4],"py":[5,4],"pz":[2,2],"nx":[9,8],"ny":[3,3],"nz":[1,-1]},{"size":3,"px":[12,4,15],"py":[5,4,7],"pz":[0,-1,-1],"nx":[3,4,2],"ny":[7,11,5],"nz":[1,1,2]},{"size":2,"px":[11,4],"py":[15,4],"pz":[0,-1],"nx":[5,9],"ny":[7,15],"nz":[1,0]},{"size":2,"px":[9,7],"py":[0,1],"pz":[1,-1],"nx":[11,11],"ny":[8,7],"nz":[1,1]},{"size":5,"px":[1,1,1,1,1],"py":[11,12,10,9,9],"pz":[0,0,0,0,-1],"nx":[4,5,8,16,11],"ny":[4,3,8,8,6],"nz":[1,1,0,0,0]}],"alpha":[-1.059083e+00,1.059083e+00,-7.846122e-01,7.846122e-01,-4.451160e-01,4.451160e-01,-4.483277e-01,4.483277e-01,-3.905999e-01,3.905999e-01,-3.789250e-01,3.789250e-01,-3.874610e-01,3.874610e-01,-3.110541e-01,3.110541e-01,-3.565056e-01,3.565056e-01,-3.812617e-01,3.812617e-01,-3.325142e-01,3.325142e-01,-2.787282e-01,2.787282e-01,-3.238869e-01,3.238869e-01,-2.993499e-01,2.993499e-01,-2.807737e-01,2.807737e-01,-2.855285e-01,2.855285e-01,-2.277550e-01,2.277550e-01,-2.031261e-01,2.031261e-01,-2.071574e-01,2.071574e-01,-2.534142e-01,2.534142e-01,-2.266871e-01,2.266871e-01,-2.229078e-01,2.229078e-01,-2.716325e-01,2.716325e-01,-3.046938e-01,3.046938e-01,-2.271601e-01,2.271601e-01,-1.987651e-01,1.987651e-01,-1.953664e-01,1.953664e-01,-2.178737e-01,2.178737e-01,-2.285148e-01,2.285148e-01,-1.891073e-01,1.891073e-01,-2.926469e-01,2.926469e-01,-2.094783e-01,2.094783e-01,-1.478037e-01,1.478037e-01,-1.707579e-01,1.707579e-01,-1.464390e-01,1.464390e-01,-2.462321e-01,2.462321e-01,-2.319978e-01,2.319978e-01,-1.781651e-01,1.781651e-01,-1.471349e-01,1.471349e-01,-1.953006e-01,1.953006e-01,-2.145108e-01,2.145108e-01,-1.567881e-01,1.567881e-01,-2.024617e-01,2.024617e-01,-1.883198e-01,1.883198e-01,-1.996976e-01,1.996976e-01,-1.292330e-01,1.292330e-01,-2.142242e-01,2.142242e-01,-2.473748e-01,2.473748e-01,-1.880902e-01,1.880902e-01,-1.874572e-01,1.874572e-01,-1.495984e-01,1.495984e-01,-1.608525e-01,1.608525e-01,-1.698402e-01,1.698402e-01,-1.898871e-01,1.898871e-01,-1.350238e-01,1.350238e-01,-1.727032e-01,1.727032e-01,-1.593352e-01,1.593352e-01,-1.476968e-01,1.476968e-01,-1.428431e-01,1.428431e-01,-1.766261e-01,1.766261e-01,-1.453226e-01,1.453226e-01,-1.929885e-01,1.929885e-01,-1.337582e-01,1.337582e-01,-1.629078e-01,1.629078e-01,-9.973085e-02,9.973085e-02,-1.172760e-01,1.172760e-01,-1.399242e-01,1.399242e-01,-1.613189e-01,1.613189e-01,-1.145695e-01,1.145695e-01,-1.191093e-01,1.191093e-01,-1.225900e-01,1.225900e-01,-1.641114e-01,1.641114e-01,-1.419878e-01,1.419878e-01,-2.183465e-01,2.183465e-01,-1.566968e-01,1.566968e-01,-1.288216e-01,1.288216e-01,-1.422831e-01,1.422831e-01,-2.000107e-01,2.000107e-01,-1.817265e-01,1.817265e-01,-1.793796e-01,1.793796e-01,-1.428926e-01,1.428926e-01,-1.182032e-01,1.182032e-01,-1.150421e-01,1.150421e-01,-1.336584e-01,1.336584e-01,-1.656178e-01,1.656178e-01,-1.386549e-01,1.386549e-01,-1.387461e-01,1.387461e-01,-1.313023e-01,1.313023e-01,-1.360391e-01,1.360391e-01,-1.305505e-01,1.305505e-01,-1.323399e-01,1.323399e-01,-1.502891e-01,1.502891e-01,-1.488859e-01,1.488859e-01,-1.126628e-01,1.126628e-01,-1.233623e-01,1.233623e-01,-1.702106e-01,1.702106e-01,-1.629639e-01,1.629639e-01,-1.337706e-01,1.337706e-01,-1.290384e-01,1.290384e-01,-1.165519e-01,1.165519e-01,-1.412778e-01,1.412778e-01,-1.470204e-01,1.470204e-01,-2.213780e-01,2.213780e-01,-1.472619e-01,1.472619e-01,-1.357071e-01,1.357071e-01,-1.416513e-01,1.416513e-01,-1.050208e-01,1.050208e-01,-1.480033e-01,1.480033e-01,-1.899871e-01,1.899871e-01,-1.466249e-01,1.466249e-01,-1.076952e-01,1.076952e-01,-1.035096e-01,1.035096e-01,-1.566970e-01,1.566970e-01,-1.364115e-01,1.364115e-01,-1.512889e-01,1.512889e-01,-1.252851e-01,1.252851e-01,-1.206300e-01,1.206300e-01,-1.059134e-01,1.059134e-01,-1.140398e-01,1.140398e-01,-1.359912e-01,1.359912e-01,-1.231201e-01,1.231201e-01,-1.231867e-01,1.231867e-01,-9.789923e-02,9.789923e-02,-1.590213e-01,1.590213e-01,-1.002206e-01,1.002206e-01,-1.518339e-01,1.518339e-01,-1.055203e-01,1.055203e-01,-1.012579e-01,1.012579e-01,-1.094956e-01,1.094956e-01,-1.429592e-01,1.429592e-01,-1.108838e-01,1.108838e-01,-1.116475e-01,1.116475e-01,-1.735371e-01,1.735371e-01,-1.067758e-01,1.067758e-01,-1.290406e-01,1.290406e-01,-1.156822e-01,1.156822e-01,-9.668217e-02,9.668217e-02,-1.170053e-01,1.170053e-01,-1.252092e-01,1.252092e-01,-1.135158e-01,1.135158e-01,-1.105896e-01,1.105896e-01,-1.038175e-01,1.038175e-01,-1.210459e-01,1.210459e-01,-1.078878e-01,1.078878e-01,-1.050808e-01,1.050808e-01,-1.428227e-01,1.428227e-01,-1.664600e-01,1.664600e-01,-1.013508e-01,1.013508e-01,-1.206930e-01,1.206930e-01,-1.088972e-01,1.088972e-01,-1.381026e-01,1.381026e-01,-1.109115e-01,1.109115e-01,-7.921549e-02,7.921549e-02,-1.057832e-01,1.057832e-01,-9.385827e-02,9.385827e-02,-1.486035e-01,1.486035e-01,-1.247401e-01,1.247401e-01,-9.451327e-02,9.451327e-02,-1.272805e-01,1.272805e-01,-9.616206e-02,9.616206e-02,-9.051084e-02,9.051084e-02,-1.138458e-01,1.138458e-01,-1.047581e-01,1.047581e-01,-1.382394e-01,1.382394e-01,-1.122203e-01,1.122203e-01,-1.052936e-01,1.052936e-01,-1.239318e-01,1.239318e-01,-1.241439e-01,1.241439e-01,-1.259012e-01,1.259012e-01,-1.211701e-01,1.211701e-01,-1.344131e-01,1.344131e-01,-1.127778e-01,1.127778e-01,-1.609745e-01,1.609745e-01,-1.901382e-01,1.901382e-01,-1.618962e-01,1.618962e-01,-1.230398e-01,1.230398e-01,-1.319311e-01,1.319311e-01,-1.431410e-01,1.431410e-01,-1.143306e-01,1.143306e-01,-9.390938e-02,9.390938e-02,-1.154161e-01,1.154161e-01,-1.141205e-01,1.141205e-01,-1.098048e-01,1.098048e-01,-8.870072e-02,8.870072e-02,-1.122444e-01,1.122444e-01,-1.114147e-01,1.114147e-01,-1.185710e-01,1.185710e-01,-1.107775e-01,1.107775e-01,-1.259167e-01,1.259167e-01,-1.105176e-01,1.105176e-01,-1.020691e-01,1.020691e-01,-9.607863e-02,9.607863e-02,-9.573700e-02,9.573700e-02,-1.054349e-01,1.054349e-01,-1.137856e-01,1.137856e-01,-1.192043e-01,1.192043e-01,-1.113264e-01,1.113264e-01,-1.093137e-01,1.093137e-01,-1.010919e-01,1.010919e-01,-9.625901e-02,9.625901e-02,-9.338459e-02,9.338459e-02,-1.142944e-01,1.142944e-01,-1.038877e-01,1.038877e-01,-9.772862e-02,9.772862e-02,-1.375298e-01,1.375298e-01,-1.394776e-01,1.394776e-01,-9.454765e-02,9.454765e-02,-1.203246e-01,1.203246e-01,-8.684943e-02,8.684943e-02,-1.135622e-01,1.135622e-01,-1.058181e-01,1.058181e-01,-1.082152e-01,1.082152e-01,-1.411355e-01,1.411355e-01,-9.978846e-02,9.978846e-02,-1.057874e-01,1.057874e-01,-1.415366e-01,1.415366e-01,-9.981014e-02,9.981014e-02,-9.261151e-02,9.261151e-02,-1.737173e-01,1.737173e-01,-1.580335e-01,1.580335e-01,-9.594668e-02,9.594668e-02,-9.336013e-02,9.336013e-02,-1.102373e-01,1.102373e-01,-8.546557e-02,8.546557e-02,-9.945057e-02,9.945057e-02,-1.146358e-01,1.146358e-01,-1.324734e-01,1.324734e-01,-1.422296e-01,1.422296e-01,-9.937990e-02,9.937990e-02,-8.381049e-02,8.381049e-02,-1.270714e-01,1.270714e-01,-1.091738e-01,1.091738e-01,-1.314881e-01,1.314881e-01,-1.085159e-01,1.085159e-01,-9.247554e-02,9.247554e-02,-8.121645e-02,8.121645e-02,-1.059589e-01,1.059589e-01,-8.307793e-02,8.307793e-02,-1.033103e-01,1.033103e-01,-1.056706e-01,1.056706e-01,-1.032803e-01,1.032803e-01,-1.266840e-01,1.266840e-01,-9.341601e-02,9.341601e-02,-7.683570e-02,7.683570e-02,-1.030530e-01,1.030530e-01,-1.051872e-01,1.051872e-01,-9.114946e-02,9.114946e-02,-1.329341e-01,1.329341e-01,-9.270830e-02,9.270830e-02,-1.141750e-01,1.141750e-01,-9.889318e-02,9.889318e-02,-8.856485e-02,8.856485e-02,-1.054210e-01,1.054210e-01,-1.092704e-01,1.092704e-01,-8.729085e-02,8.729085e-02,-1.141057e-01,1.141057e-01,-1.530774e-01,1.530774e-01,-8.129720e-02,8.129720e-02,-1.143335e-01,1.143335e-01,-1.175777e-01,1.175777e-01,-1.371729e-01,1.371729e-01,-1.394356e-01,1.394356e-01,-1.016308e-01,1.016308e-01,-1.125547e-01,1.125547e-01,-9.672600e-02,9.672600e-02,-1.036631e-01,1.036631e-01,-8.702514e-02,8.702514e-02,-1.264807e-01,1.264807e-01,-1.465688e-01,1.465688e-01,-8.781464e-02,8.781464e-02,-8.552605e-02,8.552605e-02,-1.145072e-01,1.145072e-01,-1.378489e-01,1.378489e-01,-1.013312e-01,1.013312e-01,-1.020083e-01,1.020083e-01,-1.015816e-01,1.015816e-01,-8.407101e-02,8.407101e-02,-8.296485e-02,8.296485e-02,-8.033655e-02,8.033655e-02,-9.003615e-02,9.003615e-02,-7.504954e-02,7.504954e-02,-1.224941e-01,1.224941e-01,-9.347814e-02,9.347814e-02,-9.555575e-02,9.555575e-02,-9.810025e-02,9.810025e-02,-1.237068e-01,1.237068e-01,-1.283586e-01,1.283586e-01,-1.082763e-01,1.082763e-01,-1.018145e-01,1.018145e-01,-1.175161e-01,1.175161e-01,-1.252279e-01,1.252279e-01,-1.370559e-01,1.370559e-01,-9.941339e-02,9.941339e-02,-8.506938e-02,8.506938e-02,-1.260902e-01,1.260902e-01,-1.014152e-01,1.014152e-01,-9.728694e-02,9.728694e-02,-9.374910e-02,9.374910e-02,-9.587429e-02,9.587429e-02,-9.516036e-02,9.516036e-02,-7.375173e-02,7.375173e-02,-9.332487e-02,9.332487e-02,-9.020733e-02,9.020733e-02,-1.133381e-01,1.133381e-01,-1.542180e-01,1.542180e-01,-9.692168e-02,9.692168e-02,-7.960904e-02,7.960904e-02,-8.947089e-02,8.947089e-02,-7.830286e-02,7.830286e-02,-9.900050e-02,9.900050e-02,-1.041293e-01,1.041293e-01,-9.572501e-02,9.572501e-02,-8.230575e-02,8.230575e-02,-9.194901e-02,9.194901e-02,-1.076971e-01,1.076971e-01,-1.027782e-01,1.027782e-01,-1.028538e-01,1.028538e-01,-1.013992e-01,1.013992e-01,-9.087585e-02,9.087585e-02,-1.100706e-01,1.100706e-01,-1.094934e-01,1.094934e-01,-1.107879e-01,1.107879e-01,-1.026915e-01,1.026915e-01,-1.017572e-01,1.017572e-01,-7.984776e-02,7.984776e-02,-9.015413e-02,9.015413e-02,-1.299870e-01,1.299870e-01,-9.164982e-02,9.164982e-02,-1.062788e-01,1.062788e-01,-1.160203e-01,1.160203e-01,-8.858603e-02,8.858603e-02,-9.762964e-02,9.762964e-02,-1.070694e-01,1.070694e-01,-9.549046e-02,9.549046e-02,-1.533034e-01,1.533034e-01,-8.663316e-02,8.663316e-02,-9.303018e-02,9.303018e-02,-9.853582e-02,9.853582e-02,-9.733371e-02,9.733371e-02,-1.048555e-01,1.048555e-01,-9.056041e-02,9.056041e-02,-7.552283e-02,7.552283e-02,-8.780631e-02,8.780631e-02,-1.123953e-01,1.123953e-01,-1.452948e-01,1.452948e-01,-1.156423e-01,1.156423e-01,-8.701142e-02,8.701142e-02,-9.713334e-02,9.713334e-02,-9.970888e-02,9.970888e-02,-8.614129e-02,8.614129e-02,-7.459861e-02,7.459861e-02,-9.253517e-02,9.253517e-02,-9.570092e-02,9.570092e-02,-9.485535e-02,9.485535e-02,-1.148365e-01,1.148365e-01,-1.063193e-01,1.063193e-01,-9.986686e-02,9.986686e-02,-7.523412e-02,7.523412e-02,-1.005881e-01,1.005881e-01,-8.249716e-02,8.249716e-02,-1.055866e-01,1.055866e-01,-1.343050e-01,1.343050e-01,-1.371056e-01,1.371056e-01,-9.604689e-02,9.604689e-02,-1.224268e-01,1.224268e-01,-9.211478e-02,9.211478e-02,-1.108371e-01,1.108371e-01,-1.100547e-01,1.100547e-01,-8.938970e-02,8.938970e-02,-8.655951e-02,8.655951e-02,-7.085816e-02,7.085816e-02,-8.101028e-02,8.101028e-02,-8.338046e-02,8.338046e-02,-8.309588e-02,8.309588e-02,-9.090584e-02,9.090584e-02,-8.124564e-02,8.124564e-02,-9.367843e-02,9.367843e-02,-1.011747e-01,1.011747e-01,-9.885045e-02,9.885045e-02,-8.944266e-02,8.944266e-02,-8.453859e-02,8.453859e-02,-8.308847e-02,8.308847e-02,-1.367280e-01,1.367280e-01,-1.295144e-01,1.295144e-01,-1.063965e-01,1.063965e-01,-7.752328e-02,7.752328e-02,-9.681524e-02,9.681524e-02,-7.862345e-02,7.862345e-02,-8.767746e-02,8.767746e-02,-9.198041e-02,9.198041e-02,-9.686489e-02,9.686489e-02]},{"count":564,"threshold":-4.517456e+00,"feature":[{"size":5,"px":[15,9,8,12,11],"py":[3,6,3,0,8],"pz":[0,1,0,0,0],"nx":[6,14,9,22,23],"ny":[8,7,8,17,3],"nz":[1,0,0,0,0]},{"size":5,"px":[12,13,11,14,12],"py":[9,4,4,4,5],"pz":[0,0,0,0,0],"nx":[4,6,10,4,15],"ny":[3,8,7,10,9],"nz":[1,1,0,1,0]},{"size":5,"px":[7,5,6,8,8],"py":[2,13,2,1,1],"pz":[0,0,0,0,-1],"nx":[3,0,4,1,0],"ny":[4,3,10,3,13],"nz":[1,1,1,0,0]},{"size":5,"px":[11,2,2,11,16],"py":[9,4,2,7,11],"pz":[0,2,2,0,0],"nx":[8,4,1,14,0],"ny":[4,4,16,5,13],"nz":[1,1,-1,-1,-1]},{"size":2,"px":[14,14],"py":[18,18],"pz":[0,-1],"nx":[8,13],"ny":[10,16],"nz":[1,0]},{"size":5,"px":[15,17,16,8,18],"py":[1,2,1,0,2],"pz":[0,0,0,1,0],"nx":[21,22,22,22,22],"ny":[1,5,3,4,2],"nz":[0,0,0,0,-1]},{"size":2,"px":[15,4],"py":[23,3],"pz":[0,2],"nx":[7,3],"ny":[10,6],"nz":[1,-1]},{"size":5,"px":[3,6,4,3,11],"py":[10,11,8,3,8],"pz":[1,0,1,1,0],"nx":[3,5,6,3,0],"ny":[4,9,9,9,0],"nz":[1,-1,-1,-1,-1]},{"size":3,"px":[11,11,2],"py":[11,13,16],"pz":[0,0,-1],"nx":[10,10,9],"ny":[10,11,14],"nz":[0,0,0]},{"size":2,"px":[8,4],"py":[12,6],"pz":[0,1],"nx":[4,5],"ny":[11,11],"nz":[1,-1]},{"size":5,"px":[10,11,13,3,12],"py":[3,4,3,0,1],"pz":[0,0,0,2,0],"nx":[14,18,20,19,15],"ny":[13,1,15,2,18],"nz":[0,0,0,0,0]},{"size":5,"px":[20,14,10,12,12],"py":[12,12,4,10,11],"pz":[0,0,1,0,0],"nx":[9,2,9,9,9],"ny":[4,12,5,9,14],"nz":[1,-1,-1,-1,-1]},{"size":5,"px":[3,3,3,4,2],"py":[15,16,14,21,12],"pz":[0,0,0,0,0],"nx":[0,0,0,0,0],"ny":[20,10,5,21,21],"nz":[0,1,2,0,-1]},{"size":2,"px":[18,8],"py":[16,7],"pz":[0,1],"nx":[14,0],"ny":[8,10],"nz":[0,-1]},{"size":4,"px":[12,4,16,1],"py":[14,3,8,3],"pz":[0,-1,-1,-1],"nx":[14,10,20,13],"ny":[13,5,16,9],"nz":[0,1,0,0]},{"size":5,"px":[3,8,2,3,3],"py":[7,2,1,2,4],"pz":[1,-1,-1,-1,-1],"nx":[1,9,2,1,1],"ny":[3,14,9,7,2],"nz":[1,0,1,1,1]},{"size":5,"px":[4,1,3,2,3],"py":[2,1,2,4,3],"pz":[0,1,0,0,0],"nx":[0,0,0,0,0],"ny":[3,1,2,0,0],"nz":[0,1,0,2,-1]},{"size":4,"px":[4,8,7,9],"py":[6,11,11,10],"pz":[1,0,0,0],"nx":[3,10,2,20],"ny":[4,4,4,8],"nz":[1,-1,-1,-1]},{"size":2,"px":[1,8],"py":[3,11],"pz":[2,-1],"nx":[8,2],"ny":[15,5],"nz":[0,2]},{"size":2,"px":[17,0],"py":[13,10],"pz":[0,-1],"nx":[14,14],"ny":[11,10],"nz":[0,0]},{"size":5,"px":[22,22,22,5,22],"py":[16,18,17,2,15],"pz":[0,0,0,2,0],"nx":[8,4,15,6,6],"ny":[4,2,7,11,11],"nz":[1,2,0,1,-1]},{"size":5,"px":[16,9,8,17,15],"py":[12,6,6,22,12],"pz":[0,1,1,0,0],"nx":[11,23,23,23,22],"ny":[11,23,22,21,23],"nz":[1,0,0,0,-1]},{"size":5,"px":[5,2,4,4,9],"py":[22,3,15,20,18],"pz":[0,2,0,0,0],"nx":[9,4,23,7,22],"ny":[8,4,22,19,23],"nz":[0,-1,-1,-1,-1]},{"size":5,"px":[8,6,9,7,3],"py":[3,3,3,3,1],"pz":[0,0,0,0,1],"nx":[5,5,4,4,4],"ny":[0,1,1,2,0],"nz":[0,0,0,0,-1]},{"size":2,"px":[2,3],"py":[3,3],"pz":[2,2],"nx":[3,6],"ny":[4,6],"nz":[1,-1]},{"size":5,"px":[1,1,0,1,0],"py":[17,15,6,16,10],"pz":[0,0,1,0,0],"nx":[4,4,7,4,8],"ny":[2,5,9,4,4],"nz":[2,2,1,2,-1]},{"size":5,"px":[12,12,12,13,13],"py":[10,9,11,13,13],"pz":[0,0,0,0,-1],"nx":[4,3,3,5,3],"ny":[21,18,17,23,16],"nz":[0,0,0,0,0]},{"size":4,"px":[5,6,5,9],"py":[13,7,9,23],"pz":[0,0,1,0],"nx":[6,15,7,5],"ny":[9,20,7,23],"nz":[0,-1,-1,-1]},{"size":2,"px":[6,3],"py":[4,2],"pz":[1,2],"nx":[8,23],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[9,7],"py":[18,0],"pz":[0,0],"nx":[5,7],"ny":[8,10],"nz":[1,1]},{"size":2,"px":[4,6],"py":[11,16],"pz":[1,0],"nx":[10,9],"ny":[16,7],"nz":[0,-1]},{"size":4,"px":[11,11,11,11],"py":[11,10,12,13],"pz":[0,0,0,0],"nx":[13,13,13,9],"ny":[11,9,10,4],"nz":[0,0,0,1]},{"size":4,"px":[12,6,7,6],"py":[7,11,8,4],"pz":[0,1,1,1],"nx":[10,0,19,7],"ny":[21,3,12,11],"nz":[0,-1,-1,-1]},{"size":2,"px":[4,4],"py":[3,4],"pz":[2,2],"nx":[9,1],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[19,19],"py":[21,20],"pz":[0,0],"nx":[7,7],"ny":[3,13],"nz":[1,-1]},{"size":5,"px":[12,9,13,11,5],"py":[0,2,2,0,0],"pz":[0,0,0,0,1],"nx":[6,4,5,5,5],"ny":[1,3,5,2,6],"nz":[0,0,1,0,1]},{"size":5,"px":[4,3,2,5,7],"py":[11,3,3,7,17],"pz":[1,2,2,0,0],"nx":[23,5,11,5,5],"ny":[0,4,10,2,6],"nz":[0,-1,-1,-1,-1]},{"size":2,"px":[20,17],"py":[12,3],"pz":[0,-1],"nx":[20,19],"ny":[21,23],"nz":[0,0]},{"size":2,"px":[2,1],"py":[12,8],"pz":[0,0],"nx":[2,8],"ny":[2,16],"nz":[2,-1]},{"size":2,"px":[16,5],"py":[4,5],"pz":[0,-1],"nx":[7,8],"ny":[9,1],"nz":[1,1]},{"size":2,"px":[2,2],"py":[0,1],"pz":[1,1],"nx":[1,8],"ny":[5,1],"nz":[0,-1]},{"size":2,"px":[1,1],"py":[12,10],"pz":[0,1],"nx":[2,20],"ny":[23,9],"nz":[0,-1]},{"size":4,"px":[11,0,0,2],"py":[14,3,9,22],"pz":[0,-1,-1,-1],"nx":[13,14,7,3],"ny":[6,7,11,1],"nz":[0,0,0,2]},{"size":2,"px":[14,0],"py":[2,3],"pz":[0,-1],"nx":[4,4],"ny":[4,3],"nz":[2,2]},{"size":2,"px":[23,11],"py":[18,11],"pz":[0,1],"nx":[3,2],"ny":[1,21],"nz":[1,-1]},{"size":2,"px":[9,9],"py":[17,14],"pz":[0,-1],"nx":[4,5],"ny":[10,8],"nz":[1,1]},{"size":2,"px":[9,18],"py":[7,14],"pz":[1,0],"nx":[18,9],"ny":[17,8],"nz":[0,-1]},{"size":2,"px":[2,8],"py":[4,22],"pz":[2,0],"nx":[4,3],"ny":[10,1],"nz":[1,-1]},{"size":2,"px":[5,22],"py":[4,9],"pz":[2,-1],"nx":[11,23],"ny":[8,14],"nz":[1,0]},{"size":3,"px":[23,5,5],"py":[8,2,1],"pz":[0,2,2],"nx":[10,10,2],"ny":[4,4,2],"nz":[1,-1,-1]},{"size":2,"px":[11,11],"py":[14,23],"pz":[0,-1],"nx":[3,11],"ny":[4,13],"nz":[1,0]},{"size":2,"px":[3,2],"py":[7,0],"pz":[1,-1],"nx":[4,3],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[12,1],"py":[19,13],"pz":[0,-1],"nx":[9,12],"ny":[10,18],"nz":[1,0]},{"size":2,"px":[10,10],"py":[11,10],"pz":[1,1],"nx":[4,1],"ny":[5,11],"nz":[2,-1]},{"size":5,"px":[9,12,4,8,8],"py":[3,5,2,9,8],"pz":[1,0,2,1,1],"nx":[23,23,23,23,23],"ny":[3,4,6,5,5],"nz":[0,0,0,0,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,9],"ny":[4,6],"nz":[1,-1]},{"size":5,"px":[13,13,13,7,7],"py":[11,10,9,6,6],"pz":[0,0,0,1,-1],"nx":[5,5,15,5,2],"ny":[5,15,9,9,1],"nz":[0,0,0,1,2]},{"size":2,"px":[19,7],"py":[21,7],"pz":[0,1],"nx":[14,10],"ny":[15,4],"nz":[0,-1]},{"size":2,"px":[5,5],"py":[3,4],"pz":[2,2],"nx":[21,0],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[2,0],"py":[0,0],"pz":[1,-1],"nx":[3,2],"ny":[1,2],"nz":[0,0]},{"size":2,"px":[9,0],"py":[4,0],"pz":[0,-1],"nx":[5,12],"ny":[0,1],"nz":[1,0]},{"size":5,"px":[14,16,12,15,13],"py":[0,1,0,0,0],"pz":[0,0,0,0,0],"nx":[4,8,8,4,9],"ny":[2,3,4,1,3],"nz":[2,1,1,2,-1]},{"size":3,"px":[4,17,2],"py":[11,14,1],"pz":[1,-1,-1],"nx":[9,8,17],"ny":[1,4,0],"nz":[1,1,0]},{"size":2,"px":[18,9],"py":[17,7],"pz":[0,1],"nx":[8,4],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[3,0],"pz":[1,2],"nx":[10,11],"ny":[6,5],"nz":[1,-1]},{"size":5,"px":[21,21,21,21,20],"py":[17,16,19,18,21],"pz":[0,0,0,0,0],"nx":[0,0,0,0,0],"ny":[4,9,11,6,6],"nz":[1,0,0,1,-1]},{"size":2,"px":[12,0],"py":[7,1],"pz":[0,-1],"nx":[8,11],"ny":[4,17],"nz":[1,0]},{"size":4,"px":[13,0,0,0],"py":[15,0,0,0],"pz":[0,-1,-1,-1],"nx":[3,7,4,6],"ny":[2,7,5,9],"nz":[2,1,2,1]},{"size":2,"px":[2,9],"py":[3,12],"pz":[2,0],"nx":[2,0],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[10,3],"py":[6,1],"pz":[1,-1],"nx":[20,21],"ny":[19,14],"nz":[0,0]},{"size":5,"px":[5,22,22,11,22],"py":[1,4,3,3,2],"pz":[2,0,0,1,-1],"nx":[7,13,14,8,15],"ny":[3,6,6,3,7],"nz":[1,0,0,1,0]},{"size":2,"px":[12,19],"py":[5,15],"pz":[0,-1],"nx":[16,4],"ny":[8,2],"nz":[0,2]},{"size":2,"px":[1,0],"py":[11,9],"pz":[1,1],"nx":[5,0],"ny":[3,3],"nz":[1,-1]},{"size":4,"px":[8,3,4,2],"py":[6,7,5,3],"pz":[1,-1,-1,-1],"nx":[13,14,11,11],"ny":[11,13,3,5],"nz":[0,0,1,1]},{"size":2,"px":[11,11],"py":[5,6],"pz":[0,0],"nx":[8,4],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[5,9],"py":[6,17],"pz":[1,0],"nx":[9,4],"ny":[15,11],"nz":[0,-1]},{"size":3,"px":[6,3,6],"py":[6,3,5],"pz":[1,2,1],"nx":[11,10,4],"ny":[8,11,5],"nz":[0,0,-1]},{"size":2,"px":[8,16],"py":[0,1],"pz":[1,-1],"nx":[19,17],"ny":[1,0],"nz":[0,0]},{"size":2,"px":[21,20],"py":[4,1],"pz":[0,0],"nx":[11,5],"ny":[0,0],"nz":[1,2]},{"size":2,"px":[8,4],"py":[6,3],"pz":[1,2],"nx":[8,9],"ny":[4,10],"nz":[1,-1]},{"size":2,"px":[10,1],"py":[0,0],"pz":[1,-1],"nx":[13,12],"ny":[6,5],"nz":[0,0]},{"size":2,"px":[5,4],"py":[3,11],"pz":[1,-1],"nx":[3,17],"ny":[1,3],"nz":[2,0]},{"size":2,"px":[12,13],"py":[4,4],"pz":[0,0],"nx":[3,3],"ny":[1,1],"nz":[2,-1]},{"size":2,"px":[3,18],"py":[2,7],"pz":[2,0],"nx":[8,1],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[16,6],"py":[8,2],"pz":[0,1],"nx":[8,9],"ny":[4,19],"nz":[1,-1]},{"size":3,"px":[12,3,14],"py":[13,3,15],"pz":[0,-1,-1],"nx":[0,1,0],"ny":[16,18,15],"nz":[0,0,0]},{"size":2,"px":[3,1],"py":[3,4],"pz":[2,-1],"nx":[7,14],"ny":[10,14],"nz":[1,0]},{"size":2,"px":[9,16],"py":[6,10],"pz":[1,0],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[7,11],"py":[4,4],"pz":[0,0],"nx":[7,23],"ny":[3,11],"nz":[0,-1]},{"size":5,"px":[2,4,3,4,4],"py":[1,2,0,1,1],"pz":[1,0,1,0,-1],"nx":[11,9,4,9,5],"ny":[6,5,3,6,3],"nz":[0,0,1,0,1]},{"size":2,"px":[6,0],"py":[14,1],"pz":[0,-1],"nx":[2,5],"ny":[2,9],"nz":[2,1]},{"size":2,"px":[6,7],"py":[7,12],"pz":[0,0],"nx":[3,22],"ny":[3,16],"nz":[1,-1]},{"size":2,"px":[10,4],"py":[1,1],"pz":[0,1],"nx":[2,6],"ny":[2,21],"nz":[2,-1]},{"size":2,"px":[13,1],"py":[11,6],"pz":[0,-1],"nx":[12,6],"ny":[5,2],"nz":[0,1]},{"size":5,"px":[10,5,11,10,10],"py":[4,3,4,6,5],"pz":[0,1,0,0,0],"nx":[4,7,13,8,4],"ny":[2,8,9,4,4],"nz":[2,1,0,1,-1]},{"size":4,"px":[7,8,7,8],"py":[11,3,4,7],"pz":[1,1,1,1],"nx":[0,7,3,8],"ny":[0,12,2,4],"nz":[0,-1,-1,-1]},{"size":2,"px":[0,0],"py":[4,7],"pz":[2,1],"nx":[10,1],"ny":[7,0],"nz":[0,-1]},{"size":2,"px":[11,5],"py":[19,5],"pz":[0,-1],"nx":[11,5],"ny":[17,10],"nz":[0,1]},{"size":2,"px":[11,12],"py":[4,4],"pz":[0,0],"nx":[7,5],"ny":[8,3],"nz":[0,-1]},{"size":3,"px":[4,8,4],"py":[2,9,4],"pz":[2,1,2],"nx":[3,19,3],"ny":[1,16,5],"nz":[1,-1,-1]},{"size":2,"px":[3,7],"py":[0,1],"pz":[1,0],"nx":[2,3],"ny":[15,2],"nz":[0,-1]},{"size":2,"px":[0,4],"py":[2,0],"pz":[2,-1],"nx":[9,16],"ny":[5,11],"nz":[1,0]},{"size":2,"px":[14,15],"py":[23,16],"pz":[0,0],"nx":[13,3],"ny":[15,1],"nz":[0,-1]},{"size":2,"px":[4,3],"py":[0,1],"pz":[1,-1],"nx":[3,7],"ny":[0,0],"nz":[1,0]},{"size":2,"px":[7,6],"py":[12,12],"pz":[0,0],"nx":[4,8],"ny":[5,4],"nz":[1,-1]},{"size":5,"px":[4,1,2,4,5],"py":[1,0,0,0,6],"pz":[0,2,1,0,1],"nx":[4,8,7,8,6],"ny":[4,10,11,4,4],"nz":[1,0,0,1,1]},{"size":2,"px":[12,12],"py":[15,8],"pz":[0,-1],"nx":[7,15],"ny":[16,14],"nz":[0,0]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[4,6],"ny":[2,8],"nz":[2,-1]},{"size":2,"px":[14,4],"py":[19,23],"pz":[0,-1],"nx":[7,14],"ny":[11,18],"nz":[1,0]},{"size":2,"px":[4,2],"py":[7,4],"pz":[1,2],"nx":[2,22],"ny":[5,19],"nz":[2,-1]},{"size":2,"px":[8,15],"py":[7,17],"pz":[1,0],"nx":[14,4],"ny":[15,5],"nz":[0,2]},{"size":2,"px":[10,11],"py":[9,8],"pz":[1,-1],"nx":[23,5],"ny":[19,4],"nz":[0,2]},{"size":2,"px":[11,1],"py":[7,9],"pz":[0,-1],"nx":[4,4],"ny":[4,5],"nz":[1,1]},{"size":2,"px":[14,7],"py":[6,9],"pz":[0,0],"nx":[4,11],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[5,4],"py":[0,5],"pz":[0,-1],"nx":[2,2],"ny":[0,4],"nz":[1,0]},{"size":2,"px":[10,22],"py":[5,20],"pz":[0,-1],"nx":[3,4],"ny":[1,2],"nz":[2,2]},{"size":3,"px":[23,11,11],"py":[17,9,8],"pz":[0,1,1],"nx":[13,8,8],"ny":[5,3,3],"nz":[0,1,-1]},{"size":2,"px":[18,9],"py":[0,21],"pz":[0,-1],"nx":[10,10],"ny":[2,1],"nz":[1,1]},{"size":5,"px":[11,10,11,11,11],"py":[11,13,10,12,12],"pz":[0,0,0,0,-1],"nx":[11,13,12,3,8],"ny":[5,5,5,1,10],"nz":[0,0,0,2,0]},{"size":2,"px":[7,8],"py":[11,11],"pz":[0,0],"nx":[9,16],"ny":[9,19],"nz":[0,-1]},{"size":2,"px":[9,18],"py":[23,7],"pz":[0,-1],"nx":[21,21],"ny":[7,13],"nz":[0,0]},{"size":2,"px":[8,8],"py":[7,8],"pz":[1,1],"nx":[5,21],"ny":[9,13],"nz":[1,-1]},{"size":2,"px":[17,8],"py":[22,8],"pz":[0,-1],"nx":[4,8],"ny":[5,10],"nz":[2,1]},{"size":5,"px":[2,5,8,8,4],"py":[3,9,13,23,7],"pz":[2,1,0,0,1],"nx":[9,17,18,19,20],"ny":[0,0,0,2,3],"nz":[1,0,0,0,0]},{"size":3,"px":[16,15,2],"py":[3,3,13],"pz":[0,0,-1],"nx":[4,8,4],"ny":[3,6,2],"nz":[2,1,2]},{"size":2,"px":[4,7],"py":[3,7],"pz":[2,1],"nx":[15,1],"ny":[15,0],"nz":[0,-1]},{"size":2,"px":[3,6],"py":[2,3],"pz":[2,1],"nx":[3,18],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[2,4],"py":[2,4],"pz":[2,1],"nx":[3,0],"ny":[5,0],"nz":[1,-1]},{"size":2,"px":[10,0],"py":[10,0],"pz":[0,-1],"nx":[9,4],"ny":[2,0],"nz":[1,2]},{"size":2,"px":[2,0],"py":[8,3],"pz":[1,-1],"nx":[4,8],"ny":[4,14],"nz":[1,0]},{"size":2,"px":[13,18],"py":[14,14],"pz":[0,-1],"nx":[1,1],"ny":[15,13],"nz":[0,0]},{"size":3,"px":[3,2,2],"py":[17,10,15],"pz":[0,1,0],"nx":[13,2,7],"ny":[19,11,0],"nz":[0,-1,-1]},{"size":2,"px":[4,17],"py":[0,2],"pz":[2,0],"nx":[8,5],"ny":[11,3],"nz":[1,-1]},{"size":2,"px":[15,21],"py":[5,4],"pz":[0,-1],"nx":[15,10],"ny":[3,0],"nz":[0,1]},{"size":2,"px":[7,3],"py":[13,8],"pz":[0,-1],"nx":[8,4],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[7,22],"py":[3,4],"pz":[1,-1],"nx":[4,2],"ny":[2,3],"nz":[1,1]},{"size":4,"px":[6,2,6,5],"py":[21,10,22,20],"pz":[0,1,0,0],"nx":[2,3,4,4],"ny":[11,21,23,23],"nz":[1,0,0,-1]},{"size":2,"px":[7,2],"py":[6,8],"pz":[1,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":4,"px":[11,11,5,11],"py":[6,5,2,4],"pz":[1,1,2,1],"nx":[13,7,8,3],"ny":[7,3,5,2],"nz":[0,1,-1,-1]},{"size":2,"px":[3,3],"py":[7,8],"pz":[1,0],"nx":[3,11],"ny":[4,2],"nz":[1,-1]},{"size":3,"px":[16,1,5],"py":[3,3,11],"pz":[0,-1,-1],"nx":[16,4,8],"ny":[2,0,1],"nz":[0,2,1]},{"size":2,"px":[10,0],"py":[8,1],"pz":[0,-1],"nx":[19,18],"ny":[20,23],"nz":[0,0]},{"size":2,"px":[17,4],"py":[10,4],"pz":[0,-1],"nx":[4,14],"ny":[2,9],"nz":[2,0]},{"size":5,"px":[11,12,9,10,11],"py":[2,3,2,2,3],"pz":[0,0,0,0,0],"nx":[6,4,2,2,2],"ny":[18,9,3,2,2],"nz":[0,1,2,2,-1]},{"size":2,"px":[0,1],"py":[6,16],"pz":[1,0],"nx":[8,16],"ny":[5,16],"nz":[0,-1]},{"size":2,"px":[3,3],"py":[2,3],"pz":[2,2],"nx":[8,17],"ny":[4,9],"nz":[1,-1]},{"size":3,"px":[2,5,2],"py":[5,6,4],"pz":[1,-1,-1],"nx":[0,0,0],"ny":[3,5,6],"nz":[2,1,1]},{"size":5,"px":[0,0,0,0,0],"py":[6,15,16,13,14],"pz":[1,0,0,0,0],"nx":[4,5,8,6,8],"ny":[4,16,8,15,4],"nz":[1,0,0,0,-1]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,5],"ny":[4,16],"nz":[1,-1]},{"size":5,"px":[21,19,21,21,21],"py":[17,23,18,19,20],"pz":[0,0,0,0,0],"nx":[5,2,3,6,6],"ny":[12,5,5,12,12],"nz":[0,1,1,0,-1]},{"size":2,"px":[5,2],"py":[11,1],"pz":[1,-1],"nx":[5,11],"ny":[3,5],"nz":[2,1]},{"size":2,"px":[10,5],"py":[5,3],"pz":[0,1],"nx":[6,15],"ny":[11,5],"nz":[1,-1]},{"size":2,"px":[6,2],"py":[4,2],"pz":[1,-1],"nx":[4,3],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[10,6],"py":[20,6],"pz":[0,-1],"nx":[5,10],"ny":[11,17],"nz":[1,0]},{"size":4,"px":[8,4,7,11],"py":[7,4,5,8],"pz":[1,2,1,0],"nx":[13,10,5,21],"ny":[9,3,5,4],"nz":[0,-1,-1,-1]},{"size":2,"px":[7,13],"py":[10,7],"pz":[0,0],"nx":[10,8],"ny":[9,18],"nz":[0,-1]},{"size":2,"px":[3,3],"py":[1,0],"pz":[2,2],"nx":[8,5],"ny":[4,2],"nz":[1,-1]},{"size":5,"px":[5,2,5,8,4],"py":[8,4,14,23,7],"pz":[1,2,0,0,1],"nx":[18,4,16,17,17],"ny":[1,0,0,1,1],"nz":[0,2,0,0,-1]},{"size":2,"px":[6,2],"py":[2,4],"pz":[1,-1],"nx":[8,8],"ny":[4,3],"nz":[1,1]},{"size":2,"px":[6,1],"py":[8,15],"pz":[0,-1],"nx":[8,3],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[10,1],"py":[7,2],"pz":[1,-1],"nx":[6,6],"ny":[9,4],"nz":[1,1]},{"size":2,"px":[4,1],"py":[6,2],"pz":[1,-1],"nx":[1,10],"ny":[16,12],"nz":[0,0]},{"size":2,"px":[8,4],"py":[7,2],"pz":[1,-1],"nx":[8,9],"ny":[8,10],"nz":[1,1]},{"size":5,"px":[4,8,7,6,6],"py":[0,0,0,1,1],"pz":[1,0,0,0,-1],"nx":[11,5,8,4,10],"ny":[5,3,4,4,5],"nz":[0,1,1,1,0]},{"size":2,"px":[5,6],"py":[8,5],"pz":[0,0],"nx":[6,6],"ny":[8,3],"nz":[0,-1]},{"size":2,"px":[18,5],"py":[19,5],"pz":[0,-1],"nx":[4,21],"ny":[5,19],"nz":[2,0]},{"size":2,"px":[9,5],"py":[13,6],"pz":[0,1],"nx":[2,2],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[10,4],"py":[17,6],"pz":[0,1],"nx":[10,2],"ny":[15,4],"nz":[0,-1]},{"size":3,"px":[13,13,19],"py":[11,12,8],"pz":[0,0,-1],"nx":[12,3,8],"ny":[4,1,4],"nz":[0,2,1]},{"size":3,"px":[11,7,4],"py":[5,2,1],"pz":[0,-1,-1],"nx":[9,2,4],"ny":[11,3,6],"nz":[0,2,1]},{"size":2,"px":[10,7],"py":[15,2],"pz":[0,-1],"nx":[4,4],"ny":[0,1],"nz":[2,2]},{"size":5,"px":[8,9,16,18,18],"py":[0,1,1,1,1],"pz":[1,1,0,0,-1],"nx":[5,5,6,4,4],"ny":[21,20,23,17,18],"nz":[0,0,0,0,0]},{"size":2,"px":[6,7],"py":[1,1],"pz":[1,1],"nx":[20,19],"ny":[2,1],"nz":[0,0]},{"size":2,"px":[2,2],"py":[10,11],"pz":[1,1],"nx":[3,3],"ny":[10,10],"nz":[1,-1]},{"size":2,"px":[9,5],"py":[23,1],"pz":[0,-1],"nx":[4,3],"ny":[10,4],"nz":[1,1]},{"size":2,"px":[1,10],"py":[4,7],"pz":[2,-1],"nx":[4,3],"ny":[23,21],"nz":[0,0]},{"size":2,"px":[10,21],"py":[11,18],"pz":[1,0],"nx":[10,4],"ny":[18,1],"nz":[0,-1]},{"size":2,"px":[11,23],"py":[11,15],"pz":[0,-1],"nx":[11,11],"ny":[7,9],"nz":[1,1]},{"size":2,"px":[10,1],"py":[7,7],"pz":[1,-1],"nx":[15,4],"ny":[14,4],"nz":[0,2]},{"size":2,"px":[1,2],"py":[9,20],"pz":[1,0],"nx":[21,3],"ny":[12,20],"nz":[0,-1]},{"size":2,"px":[7,4],"py":[0,0],"pz":[1,2],"nx":[4,2],"ny":[0,19],"nz":[0,-1]},{"size":2,"px":[2,4],"py":[3,6],"pz":[2,1],"nx":[3,0],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[5,1],"py":[5,0],"pz":[1,-1],"nx":[12,10],"ny":[11,4],"nz":[0,1]},{"size":2,"px":[11,12],"py":[11,14],"pz":[1,-1],"nx":[18,16],"ny":[21,15],"nz":[0,0]},{"size":2,"px":[3,18],"py":[1,5],"pz":[2,-1],"nx":[4,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[9,10],"py":[18,7],"pz":[0,-1],"nx":[3,6],"ny":[0,0],"nz":[2,1]},{"size":2,"px":[19,2],"py":[1,4],"pz":[0,-1],"nx":[22,22],"ny":[13,15],"nz":[0,0]},{"size":3,"px":[13,15,20],"py":[14,21,10],"pz":[0,-1,-1],"nx":[15,7,7],"ny":[13,6,8],"nz":[0,1,1]},{"size":2,"px":[9,9],"py":[6,7],"pz":[1,1],"nx":[8,7],"ny":[4,8],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[5,3],"pz":[1,2],"nx":[5,10],"ny":[2,9],"nz":[1,-1]},{"size":2,"px":[14,11],"py":[7,16],"pz":[0,-1],"nx":[1,0],"ny":[17,4],"nz":[0,2]},{"size":2,"px":[14,18],"py":[17,18],"pz":[0,-1],"nx":[8,14],"ny":[10,16],"nz":[1,0]},{"size":2,"px":[6,11],"py":[13,11],"pz":[0,-1],"nx":[8,9],"ny":[12,9],"nz":[0,0]},{"size":2,"px":[8,9],"py":[2,2],"pz":[0,0],"nx":[3,3],"ny":[2,2],"nz":[2,-1]},{"size":3,"px":[21,21,21],"py":[14,16,15],"pz":[0,0,0],"nx":[14,12,0],"ny":[5,12,6],"nz":[0,-1,-1]},{"size":2,"px":[4,21],"py":[6,15],"pz":[1,-1],"nx":[5,1],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[6,3],"py":[2,1],"pz":[1,2],"nx":[8,0],"ny":[4,20],"nz":[1,-1]},{"size":2,"px":[13,2],"py":[9,1],"pz":[0,-1],"nx":[3,5],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[16,1],"py":[5,4],"pz":[0,-1],"nx":[17,8],"ny":[3,2],"nz":[0,1]},{"size":2,"px":[9,2],"py":[7,1],"pz":[1,-1],"nx":[20,20],"ny":[17,16],"nz":[0,0]},{"size":2,"px":[5,7],"py":[3,6],"pz":[2,-1],"nx":[9,9],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[11,17],"py":[4,1],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[15,2],"py":[11,0],"pz":[0,-1],"nx":[5,14],"ny":[1,12],"nz":[2,0]},{"size":2,"px":[22,19],"py":[3,0],"pz":[0,-1],"nx":[9,4],"ny":[6,4],"nz":[1,1]},{"size":2,"px":[1,22],"py":[3,21],"pz":[0,-1],"nx":[0,0],"ny":[1,0],"nz":[2,2]},{"size":2,"px":[11,11],"py":[11,12],"pz":[0,0],"nx":[1,2],"ny":[1,4],"nz":[2,-1]},{"size":2,"px":[18,3],"py":[8,1],"pz":[0,2],"nx":[13,1],"ny":[8,5],"nz":[0,-1]},{"size":2,"px":[13,6],"py":[21,3],"pz":[0,-1],"nx":[11,11],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[15,14],"py":[4,4],"pz":[0,0],"nx":[17,1],"ny":[12,5],"nz":[0,-1]},{"size":2,"px":[11,3],"py":[12,1],"pz":[0,-1],"nx":[1,2],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[3,2],"py":[7,3],"pz":[0,1],"nx":[16,2],"ny":[3,5],"nz":[0,-1]},{"size":2,"px":[10,5],"py":[7,20],"pz":[1,-1],"nx":[9,8],"ny":[4,6],"nz":[1,1]},{"size":2,"px":[19,2],"py":[10,2],"pz":[0,-1],"nx":[9,4],"ny":[3,1],"nz":[1,2]},{"size":2,"px":[14,9],"py":[0,23],"pz":[0,-1],"nx":[4,4],"ny":[3,2],"nz":[2,2]},{"size":2,"px":[6,9],"py":[4,10],"pz":[1,0],"nx":[10,9],"ny":[9,0],"nz":[0,-1]},{"size":4,"px":[6,9,10,8],"py":[20,23,18,23],"pz":[0,0,0,0],"nx":[9,22,1,2],"ny":[21,14,2,5],"nz":[0,-1,-1,-1]},{"size":2,"px":[17,18],"py":[13,6],"pz":[0,-1],"nx":[6,7],"ny":[9,11],"nz":[1,1]},{"size":5,"px":[18,19,20,19,20],"py":[15,19,16,20,17],"pz":[0,0,0,0,0],"nx":[11,22,23,23,23],"ny":[10,22,20,19,19],"nz":[1,0,0,0,-1]},{"size":2,"px":[10,10],"py":[1,0],"pz":[1,1],"nx":[21,11],"ny":[0,4],"nz":[0,-1]},{"size":2,"px":[11,0],"py":[9,3],"pz":[0,-1],"nx":[9,4],"ny":[2,1],"nz":[1,2]},{"size":2,"px":[14,23],"py":[2,18],"pz":[0,-1],"nx":[15,18],"ny":[1,2],"nz":[0,0]},{"size":2,"px":[9,3],"py":[0,0],"pz":[1,-1],"nx":[3,12],"ny":[1,5],"nz":[2,0]},{"size":2,"px":[8,8],"py":[7,8],"pz":[1,1],"nx":[8,8],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[1,0],"py":[1,3],"pz":[2,-1],"nx":[7,19],"ny":[9,15],"nz":[1,0]},{"size":3,"px":[16,6,4],"py":[21,5,4],"pz":[0,-1,-1],"nx":[4,19,8],"ny":[5,21,11],"nz":[2,0,1]},{"size":2,"px":[5,5],"py":[6,6],"pz":[1,-1],"nx":[10,10],"ny":[10,12],"nz":[0,0]},{"size":2,"px":[6,11],"py":[2,5],"pz":[1,0],"nx":[3,4],"ny":[4,7],"nz":[1,-1]},{"size":3,"px":[8,6,2],"py":[4,10,2],"pz":[1,1,2],"nx":[2,18,5],"ny":[0,11,5],"nz":[0,-1,-1]},{"size":2,"px":[11,7],"py":[9,7],"pz":[0,-1],"nx":[12,3],"ny":[9,5],"nz":[0,1]},{"size":2,"px":[14,13],"py":[20,20],"pz":[0,0],"nx":[13,3],"ny":[21,5],"nz":[0,-1]},{"size":2,"px":[13,7],"py":[5,3],"pz":[0,-1],"nx":[3,4],"ny":[1,4],"nz":[2,1]},{"size":2,"px":[6,2],"py":[21,5],"pz":[0,-1],"nx":[2,3],"ny":[5,10],"nz":[2,1]},{"size":2,"px":[23,5],"py":[6,0],"pz":[0,2],"nx":[21,4],"ny":[6,1],"nz":[0,-1]},{"size":2,"px":[9,9],"py":[7,6],"pz":[1,1],"nx":[8,2],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[22,11],"py":[20,9],"pz":[0,1],"nx":[8,8],"ny":[10,10],"nz":[1,-1]},{"size":2,"px":[8,16],"py":[21,12],"pz":[0,-1],"nx":[2,7],"ny":[5,23],"nz":[2,0]},{"size":5,"px":[0,1,1,1,1],"py":[3,1,9,4,7],"pz":[2,2,1,1,1],"nx":[11,22,22,23,23],"ny":[10,21,22,19,20],"nz":[1,0,0,0,-1]},{"size":2,"px":[17,5],"py":[12,4],"pz":[0,-1],"nx":[8,8],"ny":[4,5],"nz":[1,1]},{"size":2,"px":[16,4],"py":[7,10],"pz":[0,-1],"nx":[9,15],"ny":[4,6],"nz":[1,0]},{"size":2,"px":[3,6],"py":[3,5],"pz":[2,1],"nx":[11,12],"ny":[11,23],"nz":[0,-1]},{"size":2,"px":[5,2],"py":[14,7],"pz":[0,1],"nx":[4,17],"ny":[18,16],"nz":[0,-1]},{"size":3,"px":[10,1,1],"py":[12,5,4],"pz":[0,-1,-1],"nx":[7,11,5],"ny":[1,2,1],"nz":[1,0,1]},{"size":2,"px":[7,6],"py":[3,9],"pz":[0,-1],"nx":[2,2],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[13,6],"py":[22,9],"pz":[0,-1],"nx":[8,4],"ny":[4,3],"nz":[1,2]},{"size":5,"px":[12,9,10,11,11],"py":[0,0,0,0,0],"pz":[0,0,0,0,-1],"nx":[16,5,10,4,8],"ny":[10,3,6,4,4],"nz":[0,1,0,1,1]},{"size":2,"px":[18,19],"py":[23,20],"pz":[0,0],"nx":[8,5],"ny":[11,3],"nz":[1,-1]},{"size":2,"px":[8,3],"py":[7,2],"pz":[1,2],"nx":[8,4],"ny":[4,3],"nz":[1,-1]},{"size":5,"px":[8,14,8,7,4],"py":[6,12,8,6,3],"pz":[1,0,1,1,2],"nx":[2,6,6,7,7],"ny":[0,1,2,0,0],"nz":[2,0,0,0,-1]},{"size":3,"px":[1,2,3],"py":[15,18,21],"pz":[0,0,0],"nx":[19,5,18],"ny":[23,5,8],"nz":[0,-1,-1]},{"size":2,"px":[6,2],"py":[6,1],"pz":[1,-1],"nx":[0,0],"ny":[12,4],"nz":[0,1]},{"size":2,"px":[3,5],"py":[5,11],"pz":[2,1],"nx":[14,5],"ny":[19,5],"nz":[0,-1]},{"size":2,"px":[10,4],"py":[4,4],"pz":[1,-1],"nx":[11,5],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[18,4],"py":[6,4],"pz":[0,-1],"nx":[4,8],"ny":[5,4],"nz":[1,1]},{"size":2,"px":[6,12],"py":[2,4],"pz":[1,0],"nx":[8,8],"ny":[3,4],"nz":[1,-1]},{"size":2,"px":[1,0],"py":[1,1],"pz":[1,2],"nx":[7,2],"ny":[4,7],"nz":[0,-1]},{"size":2,"px":[8,0],"py":[20,0],"pz":[0,-1],"nx":[4,5],"ny":[10,11],"nz":[1,1]},{"size":2,"px":[6,14],"py":[5,2],"pz":[1,-1],"nx":[0,0],"ny":[0,2],"nz":[1,0]},{"size":2,"px":[5,15],"py":[4,7],"pz":[1,-1],"nx":[4,7],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[7,5],"py":[2,1],"pz":[0,1],"nx":[3,1],"ny":[4,1],"nz":[1,-1]},{"size":2,"px":[8,9],"py":[4,2],"pz":[0,-1],"nx":[11,9],"ny":[1,3],"nz":[0,0]},{"size":2,"px":[6,3],"py":[2,4],"pz":[1,-1],"nx":[4,8],"ny":[4,4],"nz":[1,1]},{"size":2,"px":[3,7],"py":[3,7],"pz":[2,1],"nx":[6,8],"ny":[14,4],"nz":[0,-1]},{"size":2,"px":[3,0],"py":[21,3],"pz":[0,2],"nx":[20,8],"ny":[10,4],"nz":[0,-1]},{"size":2,"px":[6,3],"py":[5,8],"pz":[0,-1],"nx":[4,3],"ny":[4,2],"nz":[0,1]},{"size":2,"px":[3,6],"py":[7,13],"pz":[1,0],"nx":[3,2],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[16,10],"py":[9,7],"pz":[0,1],"nx":[7,9],"ny":[3,10],"nz":[1,-1]},{"size":2,"px":[13,10],"py":[6,7],"pz":[0,-1],"nx":[8,17],"ny":[4,12],"nz":[1,0]},{"size":2,"px":[5,10],"py":[4,10],"pz":[2,1],"nx":[5,4],"ny":[9,2],"nz":[1,-1]},{"size":4,"px":[15,3,5,0],"py":[12,4,2,3],"pz":[0,-1,-1,-1],"nx":[13,7,5,7],"ny":[12,6,0,7],"nz":[0,1,2,1]},{"size":4,"px":[2,3,16,17],"py":[3,4,6,6],"pz":[2,1,0,0],"nx":[16,16,8,16],"ny":[8,3,10,13],"nz":[0,-1,-1,-1]},{"size":2,"px":[16,8],"py":[1,4],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[9,14],"py":[6,2],"pz":[1,-1],"nx":[8,8],"ny":[6,4],"nz":[1,1]},{"size":2,"px":[8,4],"py":[10,4],"pz":[1,2],"nx":[10,0],"ny":[5,7],"nz":[1,-1]},{"size":2,"px":[9,10],"py":[4,4],"pz":[0,0],"nx":[9,7],"ny":[3,5],"nz":[0,-1]},{"size":5,"px":[11,10,13,6,12],"py":[2,2,2,1,2],"pz":[0,0,0,1,0],"nx":[4,18,18,13,13],"ny":[2,18,19,7,7],"nz":[2,0,0,0,-1]},{"size":4,"px":[13,13,13,2],"py":[13,12,11,3],"pz":[0,0,0,-1],"nx":[4,6,8,11],"ny":[2,2,4,4],"nz":[2,1,1,0]},{"size":2,"px":[4,7],"py":[6,13],"pz":[1,0],"nx":[8,10],"ny":[4,22],"nz":[1,-1]},{"size":2,"px":[0,7],"py":[4,17],"pz":[1,-1],"nx":[0,1],"ny":[5,21],"nz":[2,0]},{"size":2,"px":[12,13],"py":[22,22],"pz":[0,0],"nx":[2,2],"ny":[13,13],"nz":[0,-1]},{"size":3,"px":[4,4,3],"py":[22,23,19],"pz":[0,0,0],"nx":[8,12,3],"ny":[22,15,2],"nz":[0,-1,-1]},{"size":2,"px":[10,12],"py":[3,13],"pz":[0,-1],"nx":[15,2],"ny":[10,2],"nz":[0,2]},{"size":2,"px":[1,1],"py":[3,3],"pz":[2,-1],"nx":[8,4],"ny":[0,0],"nz":[1,2]},{"size":2,"px":[6,12],"py":[6,18],"pz":[1,0],"nx":[12,19],"ny":[17,16],"nz":[0,-1]},{"size":2,"px":[10,5],"py":[2,1],"pz":[0,1],"nx":[5,4],"ny":[4,17],"nz":[0,-1]},{"size":3,"px":[3,12,11],"py":[5,23,23],"pz":[2,0,0],"nx":[12,4,4],"ny":[21,17,1],"nz":[0,-1,-1]},{"size":2,"px":[12,0],"py":[21,5],"pz":[0,-1],"nx":[0,0],"ny":[7,9],"nz":[1,1]},{"size":2,"px":[17,17],"py":[12,11],"pz":[0,0],"nx":[8,11],"ny":[4,11],"nz":[1,-1]},{"size":2,"px":[11,0],"py":[22,1],"pz":[0,-1],"nx":[4,6],"ny":[1,0],"nz":[1,1]},{"size":2,"px":[11,11],"py":[9,5],"pz":[1,1],"nx":[23,11],"ny":[23,20],"nz":[0,-1]},{"size":5,"px":[4,12,11,9,8],"py":[0,1,1,0,1],"pz":[1,0,0,0,0],"nx":[4,17,8,7,7],"ny":[2,13,4,4,4],"nz":[2,0,1,1,-1]},{"size":2,"px":[11,13],"py":[12,12],"pz":[0,-1],"nx":[1,1],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[23,4],"py":[23,2],"pz":[0,-1],"nx":[5,2],"ny":[23,6],"nz":[0,1]},{"size":3,"px":[8,16,0],"py":[5,15,6],"pz":[1,-1,-1],"nx":[23,23,11],"ny":[18,17,8],"nz":[0,0,1]},{"size":2,"px":[1,16],"py":[4,15],"pz":[2,-1],"nx":[2,2],"ny":[3,2],"nz":[2,2]},{"size":2,"px":[3,8],"py":[7,9],"pz":[1,-1],"nx":[4,2],"ny":[10,5],"nz":[1,2]},{"size":3,"px":[22,1,9],"py":[23,2,3],"pz":[0,-1,-1],"nx":[2,2,5],"ny":[5,4,19],"nz":[2,2,0]},{"size":2,"px":[2,20],"py":[5,15],"pz":[1,-1],"nx":[2,1],"ny":[1,2],"nz":[2,2]},{"size":2,"px":[4,8],"py":[1,19],"pz":[1,-1],"nx":[2,2],"ny":[5,4],"nz":[2,2]},{"size":2,"px":[9,10],"py":[21,0],"pz":[0,-1],"nx":[6,5],"ny":[1,1],"nz":[1,1]},{"size":2,"px":[4,8],"py":[3,6],"pz":[2,1],"nx":[9,2],"ny":[4,1],"nz":[1,-1]},{"size":3,"px":[17,3,10],"py":[8,0,2],"pz":[0,2,0],"nx":[13,2,6],"ny":[15,5,1],"nz":[0,-1,-1]},{"size":2,"px":[9,6],"py":[20,21],"pz":[0,-1],"nx":[4,2],"ny":[10,5],"nz":[1,2]},{"size":2,"px":[3,7],"py":[0,1],"pz":[2,1],"nx":[7,20],"ny":[1,19],"nz":[0,-1]},{"size":2,"px":[4,5],"py":[0,1],"pz":[1,0],"nx":[3,2],"ny":[4,2],"nz":[0,-1]},{"size":2,"px":[2,7],"py":[4,19],"pz":[2,0],"nx":[5,2],"ny":[10,2],"nz":[1,-1]},{"size":5,"px":[3,3,4,7,7],"py":[1,0,0,0,1],"pz":[1,1,1,0,0],"nx":[5,4,10,8,8],"ny":[3,3,5,4,4],"nz":[1,1,0,1,-1]},{"size":2,"px":[1,5],"py":[0,3],"pz":[1,-1],"nx":[1,0],"ny":[0,1],"nz":[0,1]},{"size":2,"px":[10,0],"py":[5,5],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[0,9],"py":[0,4],"pz":[2,-1],"nx":[13,10],"ny":[0,0],"nz":[0,0]},{"size":2,"px":[13,4],"py":[14,5],"pz":[0,-1],"nx":[4,2],"ny":[0,0],"nz":[0,1]},{"size":2,"px":[17,4],"py":[13,3],"pz":[0,-1],"nx":[4,2],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[1,0],"py":[6,2],"pz":[1,-1],"nx":[1,6],"ny":[2,12],"nz":[2,0]},{"size":2,"px":[12,4],"py":[6,0],"pz":[0,-1],"nx":[3,3],"ny":[8,9],"nz":[1,1]},{"size":2,"px":[1,5],"py":[1,5],"pz":[1,-1],"nx":[17,17],"ny":[13,7],"nz":[0,0]},{"size":2,"px":[7,3],"py":[12,6],"pz":[0,1],"nx":[3,4],"ny":[4,11],"nz":[1,-1]},{"size":2,"px":[6,17],"py":[2,8],"pz":[1,0],"nx":[3,3],"ny":[1,2],"nz":[1,-1]},{"size":3,"px":[13,6,6],"py":[22,11,10],"pz":[0,1,1],"nx":[13,12,11],"ny":[20,20,20],"nz":[0,0,0]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,12],"ny":[4,20],"nz":[1,-1]},{"size":2,"px":[5,2],"py":[1,1],"pz":[1,-1],"nx":[13,6],"ny":[0,0],"nz":[0,1]},{"size":2,"px":[2,8],"py":[3,9],"pz":[2,0],"nx":[8,16],"ny":[5,17],"nz":[0,-1]},{"size":2,"px":[16,15],"py":[1,1],"pz":[0,0],"nx":[7,11],"ny":[8,0],"nz":[1,-1]},{"size":2,"px":[11,18],"py":[21,23],"pz":[0,-1],"nx":[1,1],"ny":[4,3],"nz":[1,2]},{"size":2,"px":[1,5],"py":[0,2],"pz":[1,-1],"nx":[15,11],"ny":[8,7],"nz":[0,0]},{"size":2,"px":[5,4],"py":[7,8],"pz":[1,-1],"nx":[9,10],"ny":[13,11],"nz":[0,0]},{"size":2,"px":[7,4],"py":[10,4],"pz":[1,2],"nx":[22,4],"ny":[0,2],"nz":[0,-1]},{"size":2,"px":[11,3],"py":[3,1],"pz":[0,2],"nx":[8,0],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[5,21],"py":[11,22],"pz":[0,-1],"nx":[10,11],"ny":[11,9],"nz":[0,0]},{"size":2,"px":[5,5],"py":[0,1],"pz":[2,2],"nx":[2,21],"ny":[6,14],"nz":[0,-1]},{"size":3,"px":[10,10,1],"py":[11,0,5],"pz":[0,-1,-1],"nx":[6,12,5],"ny":[2,5,2],"nz":[1,0,1]},{"size":2,"px":[9,10],"py":[5,6],"pz":[0,0],"nx":[12,19],"ny":[23,5],"nz":[0,-1]},{"size":2,"px":[11,5],"py":[9,6],"pz":[0,1],"nx":[21,0],"ny":[23,0],"nz":[0,-1]},{"size":2,"px":[13,12],"py":[19,15],"pz":[0,0],"nx":[13,0],"ny":[17,0],"nz":[0,-1]},{"size":2,"px":[14,0],"py":[17,3],"pz":[0,-1],"nx":[7,16],"ny":[8,19],"nz":[1,0]},{"size":2,"px":[3,6],"py":[2,4],"pz":[2,1],"nx":[8,1],"ny":[4,4],"nz":[1,-1]},{"size":2,"px":[13,10],"py":[23,20],"pz":[0,-1],"nx":[4,7],"ny":[5,10],"nz":[2,1]},{"size":2,"px":[16,9],"py":[22,5],"pz":[0,-1],"nx":[4,2],"ny":[10,3],"nz":[1,2]},{"size":4,"px":[3,1,1,5],"py":[4,2,1,2],"pz":[0,2,2,1],"nx":[13,5,8,0],"ny":[22,2,9,2],"nz":[0,-1,-1,-1]},{"size":2,"px":[9,9],"py":[0,0],"pz":[1,-1],"nx":[19,20],"ny":[1,2],"nz":[0,0]},{"size":2,"px":[7,22],"py":[6,8],"pz":[1,0],"nx":[4,4],"ny":[2,4],"nz":[2,-1]},{"size":2,"px":[3,6],"py":[4,4],"pz":[2,1],"nx":[10,20],"ny":[10,6],"nz":[0,-1]},{"size":2,"px":[6,12],"py":[6,15],"pz":[1,-1],"nx":[0,0],"ny":[2,5],"nz":[2,1]},{"size":2,"px":[2,7],"py":[4,10],"pz":[2,-1],"nx":[3,6],"ny":[4,8],"nz":[2,1]},{"size":3,"px":[11,11,4],"py":[0,5,7],"pz":[1,-1,-1],"nx":[6,12,12],"ny":[1,1,2],"nz":[1,0,0]},{"size":2,"px":[11,17],"py":[4,18],"pz":[0,-1],"nx":[8,2],"ny":[10,2],"nz":[0,2]},{"size":2,"px":[17,17],"py":[10,18],"pz":[0,-1],"nx":[8,8],"ny":[2,3],"nz":[1,1]},{"size":2,"px":[9,9],"py":[7,7],"pz":[1,-1],"nx":[7,4],"ny":[6,3],"nz":[1,2]},{"size":2,"px":[18,21],"py":[0,0],"pz":[0,-1],"nx":[11,6],"ny":[5,3],"nz":[0,1]},{"size":2,"px":[5,2],"py":[8,4],"pz":[0,2],"nx":[5,8],"ny":[9,16],"nz":[0,-1]},{"size":2,"px":[12,2],"py":[5,4],"pz":[0,-1],"nx":[4,15],"ny":[4,8],"nz":[1,0]},{"size":2,"px":[1,1],"py":[4,6],"pz":[1,1],"nx":[11,3],"ny":[7,9],"nz":[0,-1]},{"size":2,"px":[2,1],"py":[3,3],"pz":[2,2],"nx":[2,2],"ny":[15,16],"nz":[0,0]},{"size":2,"px":[17,18],"py":[5,5],"pz":[0,0],"nx":[9,21],"ny":[2,10],"nz":[1,-1]},{"size":2,"px":[6,3],"py":[14,7],"pz":[0,1],"nx":[3,4],"ny":[4,5],"nz":[1,-1]},{"size":2,"px":[0,3],"py":[3,1],"pz":[1,-1],"nx":[19,10],"ny":[12,4],"nz":[0,1]},{"size":2,"px":[6,16],"py":[3,8],"pz":[1,0],"nx":[8,10],"ny":[20,4],"nz":[0,-1]},{"size":3,"px":[5,5,2],"py":[21,8,4],"pz":[0,1,2],"nx":[10,6,3],"ny":[15,2,1],"nz":[0,-1,-1]},{"size":2,"px":[11,10],"py":[10,12],"pz":[0,0],"nx":[11,11],"ny":[2,1],"nz":[1,-1]},{"size":2,"px":[10,10],"py":[3,2],"pz":[1,1],"nx":[8,11],"ny":[3,5],"nz":[1,-1]},{"size":2,"px":[13,3],"py":[5,8],"pz":[0,-1],"nx":[12,3],"ny":[3,1],"nz":[0,2]},{"size":2,"px":[13,7],"py":[2,1],"pz":[0,1],"nx":[5,5],"ny":[1,1],"nz":[0,-1]},{"size":2,"px":[11,10],"py":[10,8],"pz":[0,-1],"nx":[14,16],"ny":[10,15],"nz":[0,0]},{"size":2,"px":[2,10],"py":[7,8],"pz":[1,-1],"nx":[2,6],"ny":[5,6],"nz":[2,1]},{"size":2,"px":[10,10],"py":[1,8],"pz":[0,-1],"nx":[2,2],"ny":[3,2],"nz":[2,2]},{"size":2,"px":[4,0],"py":[5,2],"pz":[1,-1],"nx":[1,2],"ny":[2,3],"nz":[2,1]},{"size":2,"px":[1,12],"py":[1,9],"pz":[2,-1],"nx":[16,17],"ny":[3,3],"nz":[0,0]},{"size":2,"px":[12,6],"py":[5,8],"pz":[0,-1],"nx":[3,4],"ny":[7,4],"nz":[1,1]},{"size":2,"px":[14,3],"py":[11,5],"pz":[0,-1],"nx":[11,4],"ny":[0,0],"nz":[0,1]},{"size":2,"px":[6,10],"py":[6,6],"pz":[1,-1],"nx":[0,0],"ny":[1,0],"nz":[2,2]},{"size":2,"px":[3,7],"py":[0,7],"pz":[1,-1],"nx":[15,13],"ny":[8,4],"nz":[0,0]},{"size":2,"px":[18,1],"py":[15,0],"pz":[0,-1],"nx":[18,18],"ny":[18,17],"nz":[0,0]},{"size":2,"px":[5,2],"py":[4,4],"pz":[0,-1],"nx":[4,18],"ny":[4,15],"nz":[1,0]},{"size":3,"px":[3,14,13],"py":[2,7,8],"pz":[2,0,0],"nx":[10,0,2],"ny":[8,3,2],"nz":[0,-1,-1]},{"size":2,"px":[16,0],"py":[14,3],"pz":[0,-1],"nx":[18,3],"ny":[12,5],"nz":[0,2]},{"size":2,"px":[5,3],"py":[8,3],"pz":[1,2],"nx":[13,4],"ny":[10,4],"nz":[0,-1]},{"size":2,"px":[3,6],"py":[1,2],"pz":[2,1],"nx":[8,1],"ny":[4,20],"nz":[1,-1]},{"size":2,"px":[10,10],"py":[8,3],"pz":[1,-1],"nx":[12,7],"ny":[2,1],"nz":[0,1]},{"size":2,"px":[17,3],"py":[9,2],"pz":[0,2],"nx":[7,6],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[12,1],"py":[2,1],"pz":[0,-1],"nx":[4,4],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[22,5],"py":[15,3],"pz":[0,2],"nx":[16,17],"ny":[14,2],"nz":[0,-1]},{"size":2,"px":[8,11],"py":[19,13],"pz":[0,-1],"nx":[0,0],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[8,11],"py":[8,1],"pz":[1,-1],"nx":[3,3],"ny":[2,5],"nz":[1,2]},{"size":3,"px":[3,8,0],"py":[7,7,5],"pz":[1,-1,-1],"nx":[11,5,1],"ny":[11,7,5],"nz":[0,1,1]},{"size":2,"px":[12,6],"py":[12,6],"pz":[0,1],"nx":[9,0],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[16,12],"py":[7,1],"pz":[0,-1],"nx":[16,7],"ny":[6,4],"nz":[0,1]},{"size":2,"px":[13,5],"py":[14,0],"pz":[0,-1],"nx":[13,10],"ny":[0,0],"nz":[0,0]},{"size":5,"px":[11,12,13,12,7],"py":[0,1,0,0,0],"pz":[0,0,0,0,1],"nx":[13,16,14,4,4],"ny":[18,23,18,5,5],"nz":[0,0,0,2,-1]},{"size":2,"px":[14,5],"py":[12,4],"pz":[0,-1],"nx":[7,7],"ny":[8,2],"nz":[1,1]},{"size":2,"px":[19,3],"py":[2,5],"pz":[0,-1],"nx":[11,23],"ny":[7,13],"nz":[1,0]},{"size":2,"px":[0,0],"py":[19,20],"pz":[0,0],"nx":[9,4],"ny":[5,2],"nz":[0,-1]},{"size":2,"px":[15,4],"py":[12,3],"pz":[0,2],"nx":[9,5],"ny":[4,5],"nz":[1,-1]},{"size":4,"px":[8,0,1,21],"py":[6,0,7,16],"pz":[1,-1,-1,-1],"nx":[11,6,11,5],"ny":[8,6,4,3],"nz":[1,1,1,2]},{"size":2,"px":[11,11],"py":[7,5],"pz":[0,-1],"nx":[9,10],"ny":[6,7],"nz":[0,0]},{"size":2,"px":[2,4],"py":[1,2],"pz":[2,1],"nx":[16,6],"ny":[0,1],"nz":[0,-1]},{"size":2,"px":[0,0],"py":[5,3],"pz":[1,2],"nx":[1,21],"ny":[23,8],"nz":[0,-1]},{"size":2,"px":[10,0],"py":[7,0],"pz":[0,-1],"nx":[4,13],"ny":[4,10],"nz":[1,0]},{"size":2,"px":[11,4],"py":[0,4],"pz":[1,-1],"nx":[4,2],"ny":[16,8],"nz":[0,1]},{"size":2,"px":[5,3],"py":[12,6],"pz":[0,1],"nx":[3,3],"ny":[4,2],"nz":[1,-1]},{"size":2,"px":[10,0],"py":[19,11],"pz":[0,-1],"nx":[9,5],"ny":[21,9],"nz":[0,1]},{"size":2,"px":[0,0],"py":[17,9],"pz":[0,1],"nx":[0,5],"ny":[0,9],"nz":[2,-1]},{"size":2,"px":[4,5],"py":[2,4],"pz":[0,-1],"nx":[4,4],"ny":[5,6],"nz":[1,1]},{"size":2,"px":[8,4],"py":[1,0],"pz":[1,2],"nx":[4,3],"ny":[3,6],"nz":[0,-1]},{"size":2,"px":[11,0],"py":[7,2],"pz":[1,-1],"nx":[5,5],"ny":[1,0],"nz":[2,2]},{"size":2,"px":[13,0],"py":[17,2],"pz":[0,-1],"nx":[3,6],"ny":[5,8],"nz":[2,1]},{"size":2,"px":[2,1],"py":[0,5],"pz":[2,-1],"nx":[4,9],"ny":[2,7],"nz":[2,1]},{"size":2,"px":[12,5],"py":[13,8],"pz":[0,-1],"nx":[23,11],"ny":[13,7],"nz":[0,1]},{"size":2,"px":[0,0],"py":[0,2],"pz":[1,0],"nx":[3,6],"ny":[11,18],"nz":[0,-1]},{"size":2,"px":[4,3],"py":[6,5],"pz":[0,-1],"nx":[1,1],"ny":[1,3],"nz":[2,1]},{"size":4,"px":[3,6,3,6],"py":[3,6,2,5],"pz":[2,1,2,1],"nx":[0,4,1,1],"ny":[0,22,17,0],"nz":[0,-1,-1,-1]},{"size":2,"px":[8,4],"py":[6,3],"pz":[1,2],"nx":[9,15],"ny":[4,8],"nz":[1,-1]},{"size":2,"px":[8,18],"py":[7,8],"pz":[1,0],"nx":[8,5],"ny":[4,0],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[4,5],"pz":[1,-1],"nx":[5,6],"ny":[0,0],"nz":[1,1]},{"size":2,"px":[13,18],"py":[23,19],"pz":[0,0],"nx":[7,13],"ny":[10,20],"nz":[1,-1]},{"size":2,"px":[10,6],"py":[2,0],"pz":[0,1],"nx":[4,1],"ny":[5,1],"nz":[1,-1]},{"size":2,"px":[1,1],"py":[5,4],"pz":[2,2],"nx":[0,20],"ny":[4,4],"nz":[2,-1]},{"size":2,"px":[5,5],"py":[1,0],"pz":[2,2],"nx":[12,6],"ny":[18,11],"nz":[0,-1]},{"size":5,"px":[2,1,3,1,5],"py":[3,3,7,4,9],"pz":[2,2,1,2,1],"nx":[9,3,8,16,10],"ny":[5,3,10,6,7],"nz":[1,-1,-1,-1,-1]},{"size":2,"px":[4,1],"py":[12,3],"pz":[0,-1],"nx":[10,1],"ny":[11,2],"nz":[0,2]},{"size":2,"px":[19,0],"py":[10,7],"pz":[0,-1],"nx":[14,7],"ny":[6,3],"nz":[0,1]},{"size":2,"px":[7,4],"py":[2,1],"pz":[1,2],"nx":[6,0],"ny":[2,18],"nz":[0,-1]},{"size":2,"px":[14,8],"py":[3,0],"pz":[0,1],"nx":[17,1],"ny":[1,4],"nz":[0,-1]},{"size":2,"px":[18,19],"py":[1,17],"pz":[0,-1],"nx":[5,11],"ny":[2,5],"nz":[2,1]},{"size":5,"px":[12,12,12,6,12],"py":[10,11,12,6,9],"pz":[0,0,0,1,0],"nx":[13,3,12,6,6],"ny":[4,1,4,2,2],"nz":[0,2,0,1,-1]},{"size":2,"px":[11,10],"py":[3,3],"pz":[0,0],"nx":[4,9],"ny":[4,17],"nz":[1,-1]},{"size":2,"px":[11,0],"py":[13,5],"pz":[0,2],"nx":[8,18],"ny":[15,15],"nz":[0,-1]},{"size":2,"px":[3,4],"py":[6,5],"pz":[1,1],"nx":[0,0],"ny":[9,4],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[1,0],"pz":[2,2],"nx":[2,15],"ny":[2,1],"nz":[2,-1]},{"size":3,"px":[2,4,2],"py":[4,9,5],"pz":[2,1,2],"nx":[2,5,14],"ny":[0,1,4],"nz":[0,-1,-1]},{"size":2,"px":[11,12],"py":[20,20],"pz":[0,0],"nx":[6,10],"ny":[9,19],"nz":[1,-1]},{"size":2,"px":[7,0],"py":[16,8],"pz":[0,-1],"nx":[2,3],"ny":[2,4],"nz":[2,1]},{"size":5,"px":[16,17,15,16,15],"py":[1,1,1,0,0],"pz":[0,0,0,0,0],"nx":[8,8,4,12,12],"ny":[8,7,2,23,23],"nz":[1,1,2,0,-1]},{"size":2,"px":[2,4],"py":[6,12],"pz":[1,-1],"nx":[8,13],"ny":[1,1],"nz":[1,0]},{"size":2,"px":[9,2],"py":[3,2],"pz":[0,-1],"nx":[3,4],"ny":[6,5],"nz":[1,1]},{"size":2,"px":[10,8],"py":[6,1],"pz":[1,-1],"nx":[11,8],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[9,3],"py":[7,0],"pz":[1,-1],"nx":[19,19],"ny":[18,16],"nz":[0,0]},{"size":2,"px":[3,2],"py":[1,1],"pz":[2,2],"nx":[22,11],"ny":[4,0],"nz":[0,-1]},{"size":2,"px":[10,10],"py":[9,8],"pz":[1,1],"nx":[4,4],"ny":[10,2],"nz":[1,-1]},{"size":2,"px":[0,1],"py":[0,5],"pz":[0,-1],"nx":[10,8],"ny":[2,2],"nz":[0,0]},{"size":2,"px":[3,3],"py":[8,7],"pz":[1,1],"nx":[8,2],"ny":[8,3],"nz":[0,-1]},{"size":2,"px":[13,5],"py":[21,3],"pz":[0,-1],"nx":[13,3],"ny":[20,5],"nz":[0,2]},{"size":2,"px":[12,5],"py":[11,2],"pz":[0,-1],"nx":[1,0],"ny":[19,9],"nz":[0,1]},{"size":2,"px":[7,10],"py":[9,10],"pz":[1,1],"nx":[8,4],"ny":[10,2],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[5,9],"pz":[2,1],"nx":[2,11],"ny":[9,19],"nz":[1,-1]},{"size":2,"px":[3,5],"py":[1,2],"pz":[2,1],"nx":[8,23],"ny":[4,9],"nz":[1,-1]},{"size":2,"px":[3,4],"py":[2,4],"pz":[2,1],"nx":[5,9],"ny":[2,5],"nz":[2,-1]},{"size":2,"px":[11,11],"py":[2,3],"pz":[1,1],"nx":[19,9],"ny":[6,5],"nz":[0,-1]},{"size":2,"px":[9,4],"py":[5,10],"pz":[1,-1],"nx":[10,22],"ny":[0,16],"nz":[1,0]},{"size":3,"px":[19,9,19],"py":[3,1,2],"pz":[0,1,0],"nx":[6,3,6],"ny":[10,3,0],"nz":[1,-1,-1]},{"size":2,"px":[8,3],"py":[10,3],"pz":[1,2],"nx":[23,14],"ny":[3,18],"nz":[0,-1]},{"size":2,"px":[11,11],"py":[19,0],"pz":[0,-1],"nx":[4,16],"ny":[4,11],"nz":[1,0]},{"size":2,"px":[22,23],"py":[3,22],"pz":[0,-1],"nx":[9,3],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[7,2],"py":[12,4],"pz":[0,-1],"nx":[8,4],"ny":[10,5],"nz":[0,1]},{"size":2,"px":[12,13],"py":[5,13],"pz":[0,-1],"nx":[11,3],"ny":[2,0],"nz":[0,2]},{"size":2,"px":[3,17],"py":[0,16],"pz":[1,-1],"nx":[12,12],"ny":[5,6],"nz":[0,0]},{"size":2,"px":[4,3],"py":[1,0],"pz":[2,2],"nx":[4,3],"ny":[0,3],"nz":[0,-1]},{"size":2,"px":[10,3],"py":[12,0],"pz":[0,-1],"nx":[12,12],"ny":[13,12],"nz":[0,0]},{"size":2,"px":[13,4],"py":[11,14],"pz":[0,-1],"nx":[0,0],"ny":[4,6],"nz":[1,0]},{"size":2,"px":[8,7],"py":[7,8],"pz":[1,1],"nx":[3,0],"ny":[5,21],"nz":[2,-1]},{"size":2,"px":[1,3],"py":[4,14],"pz":[2,0],"nx":[8,8],"ny":[7,7],"nz":[1,-1]},{"size":2,"px":[13,11],"py":[20,7],"pz":[0,-1],"nx":[21,21],"ny":[20,18],"nz":[0,0]},{"size":2,"px":[2,1],"py":[11,0],"pz":[0,-1],"nx":[2,2],"ny":[15,14],"nz":[0,0]},{"size":2,"px":[10,1],"py":[8,0],"pz":[1,-1],"nx":[8,4],"ny":[7,4],"nz":[1,2]},{"size":2,"px":[17,6],"py":[13,1],"pz":[0,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[7,15],"py":[1,3],"pz":[1,0],"nx":[15,5],"ny":[1,8],"nz":[0,-1]},{"size":2,"px":[16,1],"py":[20,10],"pz":[0,-1],"nx":[6,8],"ny":[11,10],"nz":[1,1]},{"size":2,"px":[7,14],"py":[0,0],"pz":[1,0],"nx":[7,8],"ny":[7,3],"nz":[1,-1]},{"size":2,"px":[12,5],"py":[17,4],"pz":[0,-1],"nx":[12,5],"ny":[16,10],"nz":[0,1]},{"size":2,"px":[13,3],"py":[15,0],"pz":[0,-1],"nx":[12,7],"ny":[17,8],"nz":[0,1]},{"size":2,"px":[7,1],"py":[14,1],"pz":[0,-1],"nx":[4,6],"ny":[6,12],"nz":[1,0]},{"size":2,"px":[8,7],"py":[0,0],"pz":[0,0],"nx":[6,20],"ny":[5,5],"nz":[0,-1]},{"size":2,"px":[10,2],"py":[22,5],"pz":[0,-1],"nx":[4,8],"ny":[4,9],"nz":[2,1]},{"size":4,"px":[8,2,2,9],"py":[6,5,3,11],"pz":[1,-1,-1,-1],"nx":[2,7,4,3],"ny":[2,1,0,2],"nz":[2,0,1,2]},{"size":2,"px":[12,6],"py":[12,6],"pz":[0,1],"nx":[8,2],"ny":[4,1],"nz":[1,-1]},{"size":2,"px":[13,11],"py":[19,8],"pz":[0,-1],"nx":[13,13],"ny":[20,17],"nz":[0,0]},{"size":2,"px":[11,19],"py":[5,14],"pz":[0,-1],"nx":[3,4],"ny":[8,4],"nz":[1,1]},{"size":2,"px":[10,0],"py":[8,6],"pz":[1,-1],"nx":[21,21],"ny":[16,15],"nz":[0,0]},{"size":2,"px":[1,12],"py":[7,6],"pz":[1,-1],"nx":[2,7],"ny":[5,14],"nz":[2,0]},{"size":2,"px":[2,9],"py":[7,5],"pz":[1,-1],"nx":[2,5],"ny":[5,9],"nz":[2,1]},{"size":2,"px":[12,5],"py":[15,6],"pz":[0,-1],"nx":[3,12],"ny":[0,2],"nz":[2,0]},{"size":2,"px":[23,22],"py":[23,1],"pz":[0,-1],"nx":[0,0],"ny":[2,3],"nz":[2,2]},{"size":2,"px":[3,6],"py":[1,2],"pz":[2,1],"nx":[8,0],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[5,1],"py":[9,1],"pz":[0,-1],"nx":[4,2],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[0,1],"py":[0,0],"pz":[2,0],"nx":[2,3],"ny":[9,10],"nz":[0,-1]},{"size":2,"px":[6,0],"py":[16,14],"pz":[0,-1],"nx":[6,3],"ny":[23,14],"nz":[0,0]},{"size":2,"px":[3,3],"py":[2,3],"pz":[2,1],"nx":[13,3],"ny":[19,14],"nz":[0,-1]},{"size":2,"px":[11,5],"py":[8,18],"pz":[0,-1],"nx":[4,7],"ny":[1,2],"nz":[2,1]},{"size":2,"px":[4,4],"py":[5,6],"pz":[1,1],"nx":[2,2],"ny":[5,3],"nz":[2,-1]},{"size":2,"px":[7,3],"py":[13,7],"pz":[0,1],"nx":[4,3],"ny":[4,1],"nz":[1,-1]},{"size":2,"px":[0,0],"py":[5,6],"pz":[1,0],"nx":[2,1],"ny":[5,1],"nz":[1,-1]},{"size":2,"px":[7,14],"py":[3,5],"pz":[1,0],"nx":[5,0],"ny":[16,7],"nz":[0,-1]},{"size":2,"px":[11,2],"py":[18,5],"pz":[0,2],"nx":[11,4],"ny":[16,4],"nz":[0,-1]},{"size":2,"px":[6,16],"py":[19,20],"pz":[0,-1],"nx":[3,2],"ny":[10,5],"nz":[1,2]},{"size":2,"px":[5,3],"py":[3,1],"pz":[0,1],"nx":[1,3],"ny":[4,8],"nz":[0,-1]},{"size":2,"px":[12,6],"py":[13,6],"pz":[0,1],"nx":[10,1],"ny":[12,2],"nz":[0,-1]},{"size":2,"px":[8,3],"py":[6,2],"pz":[1,-1],"nx":[4,8],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[9,3],"py":[21,2],"pz":[0,-1],"nx":[8,4],"ny":[1,0],"nz":[1,2]},{"size":2,"px":[8,4],"py":[1,0],"pz":[1,-1],"nx":[8,6],"ny":[4,2],"nz":[1,1]},{"size":2,"px":[2,7],"py":[1,6],"pz":[2,-1],"nx":[7,9],"ny":[6,4],"nz":[1,1]},{"size":2,"px":[6,3],"py":[8,3],"pz":[1,2],"nx":[10,5],"ny":[19,11],"nz":[0,-1]},{"size":2,"px":[2,2],"py":[3,4],"pz":[2,2],"nx":[3,6],"ny":[4,6],"nz":[1,-1]},{"size":2,"px":[3,11],"py":[5,20],"pz":[2,0],"nx":[11,5],"ny":[21,8],"nz":[0,-1]},{"size":3,"px":[5,9,5],"py":[4,7,5],"pz":[2,0,2],"nx":[23,10,4],"ny":[23,3,22],"nz":[0,-1,-1]},{"size":4,"px":[11,9,7,1],"py":[13,8,11,10],"pz":[0,-1,-1,-1],"nx":[8,2,11,12],"ny":[4,2,4,4],"nz":[1,2,0,0]},{"size":2,"px":[0,0],"py":[7,6],"pz":[1,1],"nx":[0,4],"ny":[1,0],"nz":[2,-1]},{"size":2,"px":[19,20],"py":[0,1],"pz":[0,0],"nx":[21,1],"ny":[0,2],"nz":[0,-1]},{"size":2,"px":[8,5],"py":[11,0],"pz":[0,-1],"nx":[11,0],"ny":[12,1],"nz":[0,2]},{"size":2,"px":[11,11],"py":[1,1],"pz":[0,-1],"nx":[4,7],"ny":[5,4],"nz":[1,1]},{"size":2,"px":[5,12],"py":[4,23],"pz":[2,-1],"nx":[13,15],"ny":[5,4],"nz":[0,0]},{"size":2,"px":[12,20],"py":[4,16],"pz":[0,-1],"nx":[9,4],"ny":[2,1],"nz":[0,1]},{"size":2,"px":[12,13],"py":[2,2],"pz":[0,0],"nx":[4,16],"ny":[2,11],"nz":[2,0]},{"size":2,"px":[19,14],"py":[10,17],"pz":[0,-1],"nx":[3,8],"ny":[0,2],"nz":[2,0]},{"size":2,"px":[8,12],"py":[1,2],"pz":[1,0],"nx":[19,10],"ny":[3,1],"nz":[0,-1]},{"size":4,"px":[17,2,3,10],"py":[8,6,2,12],"pz":[0,1,2,0],"nx":[17,9,12,2],"ny":[9,22,13,5],"nz":[0,-1,-1,-1]},{"size":2,"px":[20,10],"py":[15,7],"pz":[0,1],"nx":[13,9],"ny":[7,3],"nz":[0,-1]},{"size":2,"px":[0,0],"py":[1,0],"pz":[2,2],"nx":[10,3],"ny":[9,2],"nz":[1,-1]},{"size":2,"px":[4,3],"py":[1,0],"pz":[2,2],"nx":[0,22],"ny":[14,6],"nz":[0,-1]},{"size":2,"px":[16,3],"py":[4,0],"pz":[0,2],"nx":[16,3],"ny":[2,0],"nz":[0,-1]},{"size":2,"px":[8,16],"py":[6,12],"pz":[1,0],"nx":[8,12],"ny":[4,7],"nz":[1,-1]},{"size":2,"px":[5,11],"py":[0,5],"pz":[2,1],"nx":[10,1],"ny":[5,5],"nz":[1,-1]},{"size":2,"px":[7,4],"py":[5,5],"pz":[0,-1],"nx":[3,6],"ny":[2,3],"nz":[1,0]},{"size":2,"px":[11,11],"py":[11,12],"pz":[0,0],"nx":[23,7],"ny":[20,2],"nz":[0,-1]},{"size":2,"px":[16,8],"py":[12,5],"pz":[0,1],"nx":[8,2],"ny":[2,1],"nz":[1,-1]},{"size":3,"px":[6,11,11],"py":[11,23,20],"pz":[1,0,0],"nx":[11,3,22],"ny":[21,3,16],"nz":[0,-1,-1]},{"size":2,"px":[17,15],"py":[3,2],"pz":[0,-1],"nx":[4,4],"ny":[3,2],"nz":[2,2]},{"size":2,"px":[21,21],"py":[11,10],"pz":[0,0],"nx":[11,3],"ny":[6,2],"nz":[1,-1]},{"size":2,"px":[23,21],"py":[22,10],"pz":[0,-1],"nx":[20,10],"ny":[18,10],"nz":[0,1]},{"size":2,"px":[4,2],"py":[6,3],"pz":[1,2],"nx":[3,2],"ny":[4,3],"nz":[1,-1]},{"size":2,"px":[16,0],"py":[18,11],"pz":[0,-1],"nx":[8,7],"ny":[4,4],"nz":[0,0]},{"size":2,"px":[6,21],"py":[3,16],"pz":[0,-1],"nx":[1,8],"ny":[2,14],"nz":[2,0]},{"size":2,"px":[8,1],"py":[3,0],"pz":[0,-1],"nx":[11,11],"ny":[2,1],"nz":[0,0]},{"size":3,"px":[11,11,11],"py":[9,10,8],"pz":[1,1,1],"nx":[23,1,0],"ny":[23,9,11],"nz":[0,-1,-1]},{"size":2,"px":[6,3],"py":[2,1],"pz":[1,2],"nx":[7,1],"ny":[8,2],"nz":[0,-1]},{"size":2,"px":[10,17],"py":[17,19],"pz":[0,-1],"nx":[10,4],"ny":[16,9],"nz":[0,1]},{"size":2,"px":[3,6],"py":[7,1],"pz":[1,-1],"nx":[11,0],"ny":[11,8],"nz":[0,1]},{"size":2,"px":[10,5],"py":[11,4],"pz":[1,2],"nx":[5,5],"ny":[0,0],"nz":[2,-1]},{"size":2,"px":[3,6],"py":[3,6],"pz":[2,1],"nx":[8,0],"ny":[4,16],"nz":[1,-1]},{"size":2,"px":[14,1],"py":[20,2],"pz":[0,-1],"nx":[7,7],"ny":[11,9],"nz":[1,1]},{"size":3,"px":[11,13,4],"py":[16,21,3],"pz":[0,0,2],"nx":[14,16,5],"ny":[20,14,9],"nz":[0,-1,-1]},{"size":2,"px":[7,0],"py":[1,1],"pz":[1,-1],"nx":[4,7],"ny":[2,4],"nz":[2,1]},{"size":2,"px":[23,11],"py":[9,4],"pz":[0,1],"nx":[11,3],"ny":[1,3],"nz":[0,-1]},{"size":2,"px":[11,13],"py":[23,23],"pz":[0,0],"nx":[13,13],"ny":[20,20],"nz":[0,-1]},{"size":2,"px":[10,8],"py":[5,11],"pz":[0,-1],"nx":[20,19],"ny":[18,20],"nz":[0,0]},{"size":2,"px":[19,5],"py":[22,4],"pz":[0,-1],"nx":[2,9],"ny":[3,17],"nz":[1,0]},{"size":2,"px":[15,2],"py":[13,7],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":2,"px":[14,13],"py":[17,2],"pz":[0,-1],"nx":[15,13],"ny":[19,15],"nz":[0,0]},{"size":2,"px":[12,23],"py":[8,22],"pz":[0,-1],"nx":[7,10],"ny":[5,9],"nz":[1,0]},{"size":2,"px":[2,6],"py":[21,10],"pz":[0,-1],"nx":[3,4],"ny":[3,3],"nz":[1,1]},{"size":2,"px":[15,11],"py":[5,0],"pz":[0,-1],"nx":[3,4],"ny":[17,16],"nz":[0,0]},{"size":2,"px":[3,1],"py":[18,8],"pz":[0,1],"nx":[14,4],"ny":[17,7],"nz":[0,-1]},{"size":2,"px":[15,3],"py":[18,3],"pz":[0,2],"nx":[1,22],"ny":[0,1],"nz":[0,-1]},{"size":2,"px":[13,3],"py":[9,3],"pz":[0,-1],"nx":[0,1],"ny":[9,20],"nz":[1,0]},{"size":2,"px":[1,1],"py":[1,0],"pz":[2,2],"nx":[9,23],"ny":[10,12],"nz":[1,-1]},{"size":4,"px":[9,0,9,1],"py":[8,0,0,10],"pz":[1,-1,-1,-1],"nx":[23,7,5,23],"ny":[20,7,5,19],"nz":[0,1,2,0]},{"size":2,"px":[18,18],"py":[12,12],"pz":[0,-1],"nx":[8,4],"ny":[4,2],"nz":[1,2]},{"size":3,"px":[0,4,1],"py":[3,5,3],"pz":[1,-1,-1],"nx":[16,11,8],"ny":[8,5,6],"nz":[0,0,0]},{"size":5,"px":[9,10,14,11,11],"py":[0,0,0,0,0],"pz":[0,0,0,0,-1],"nx":[8,3,4,6,2],"ny":[22,9,5,4,0],"nz":[0,1,0,0,2]},{"size":2,"px":[6,5],"py":[2,2],"pz":[1,1],"nx":[7,3],"ny":[8,7],"nz":[0,-1]},{"size":2,"px":[11,5],"py":[15,2],"pz":[0,-1],"nx":[3,10],"ny":[0,1],"nz":[2,0]},{"size":2,"px":[0,11],"py":[11,12],"pz":[1,-1],"nx":[22,22],"ny":[14,13],"nz":[0,0]},{"size":2,"px":[2,2],"py":[15,14],"pz":[0,0],"nx":[1,2],"ny":[11,8],"nz":[1,-1]},{"size":2,"px":[11,6],"py":[0,7],"pz":[1,-1],"nx":[19,5],"ny":[3,0],"nz":[0,2]},{"size":2,"px":[2,3],"py":[3,7],"pz":[2,1],"nx":[1,5],"ny":[5,0],"nz":[1,-1]},{"size":2,"px":[10,14],"py":[4,5],"pz":[0,-1],"nx":[4,18],"ny":[2,12],"nz":[2,0]},{"size":2,"px":[19,10],"py":[12,2],"pz":[0,-1],"nx":[13,4],"ny":[10,2],"nz":[0,2]},{"size":2,"px":[6,1],"py":[21,6],"pz":[0,-1],"nx":[6,5],"ny":[0,0],"nz":[1,1]}],"alpha":[-1.044179e+00,1.044179e+00,-6.003138e-01,6.003138e-01,-4.091282e-01,4.091282e-01,-4.590148e-01,4.590148e-01,-4.294004e-01,4.294004e-01,-3.360846e-01,3.360846e-01,-3.054186e-01,3.054186e-01,-2.901743e-01,2.901743e-01,-3.522417e-01,3.522417e-01,-3.195838e-01,3.195838e-01,-2.957309e-01,2.957309e-01,-2.876727e-01,2.876727e-01,-2.637460e-01,2.637460e-01,-2.607900e-01,2.607900e-01,-2.455714e-01,2.455714e-01,-2.749847e-01,2.749847e-01,-2.314217e-01,2.314217e-01,-2.540871e-01,2.540871e-01,-2.143416e-01,2.143416e-01,-2.565697e-01,2.565697e-01,-1.901272e-01,1.901272e-01,-2.259981e-01,2.259981e-01,-2.012333e-01,2.012333e-01,-2.448460e-01,2.448460e-01,-2.192845e-01,2.192845e-01,-2.005951e-01,2.005951e-01,-2.259000e-01,2.259000e-01,-1.955758e-01,1.955758e-01,-2.235332e-01,2.235332e-01,-1.704490e-01,1.704490e-01,-1.584628e-01,1.584628e-01,-2.167710e-01,2.167710e-01,-1.592909e-01,1.592909e-01,-1.967292e-01,1.967292e-01,-1.432268e-01,1.432268e-01,-2.039949e-01,2.039949e-01,-1.404068e-01,1.404068e-01,-1.788201e-01,1.788201e-01,-1.498714e-01,1.498714e-01,-1.282541e-01,1.282541e-01,-1.630182e-01,1.630182e-01,-1.398111e-01,1.398111e-01,-1.464143e-01,1.464143e-01,-1.281712e-01,1.281712e-01,-1.417014e-01,1.417014e-01,-1.779164e-01,1.779164e-01,-2.067174e-01,2.067174e-01,-1.344947e-01,1.344947e-01,-1.357351e-01,1.357351e-01,-1.683191e-01,1.683191e-01,-1.821768e-01,1.821768e-01,-2.158307e-01,2.158307e-01,-1.812857e-01,1.812857e-01,-1.635445e-01,1.635445e-01,-1.474934e-01,1.474934e-01,-1.771993e-01,1.771993e-01,-1.517620e-01,1.517620e-01,-1.283184e-01,1.283184e-01,-1.862675e-01,1.862675e-01,-1.420491e-01,1.420491e-01,-1.232165e-01,1.232165e-01,-1.472696e-01,1.472696e-01,-1.192156e-01,1.192156e-01,-1.602034e-01,1.602034e-01,-1.321473e-01,1.321473e-01,-1.358101e-01,1.358101e-01,-1.295821e-01,1.295821e-01,-1.289102e-01,1.289102e-01,-1.232520e-01,1.232520e-01,-1.332227e-01,1.332227e-01,-1.358887e-01,1.358887e-01,-1.179559e-01,1.179559e-01,-1.263694e-01,1.263694e-01,-1.444876e-01,1.444876e-01,-1.933141e-01,1.933141e-01,-1.917886e-01,1.917886e-01,-1.199760e-01,1.199760e-01,-1.359937e-01,1.359937e-01,-1.690073e-01,1.690073e-01,-1.894222e-01,1.894222e-01,-1.699422e-01,1.699422e-01,-1.340361e-01,1.340361e-01,-1.840622e-01,1.840622e-01,-1.277397e-01,1.277397e-01,-1.381610e-01,1.381610e-01,-1.282241e-01,1.282241e-01,-1.211334e-01,1.211334e-01,-1.264628e-01,1.264628e-01,-1.373010e-01,1.373010e-01,-1.363356e-01,1.363356e-01,-1.562568e-01,1.562568e-01,-1.268735e-01,1.268735e-01,-1.037859e-01,1.037859e-01,-1.394322e-01,1.394322e-01,-1.449225e-01,1.449225e-01,-1.109657e-01,1.109657e-01,-1.086931e-01,1.086931e-01,-1.379135e-01,1.379135e-01,-1.881974e-01,1.881974e-01,-1.304956e-01,1.304956e-01,-9.921777e-02,9.921777e-02,-1.398624e-01,1.398624e-01,-1.216469e-01,1.216469e-01,-1.272741e-01,1.272741e-01,-1.878236e-01,1.878236e-01,-1.336894e-01,1.336894e-01,-1.256289e-01,1.256289e-01,-1.247231e-01,1.247231e-01,-1.853400e-01,1.853400e-01,-1.087805e-01,1.087805e-01,-1.205676e-01,1.205676e-01,-1.023182e-01,1.023182e-01,-1.268422e-01,1.268422e-01,-1.422900e-01,1.422900e-01,-1.098174e-01,1.098174e-01,-1.317018e-01,1.317018e-01,-1.378142e-01,1.378142e-01,-1.274550e-01,1.274550e-01,-1.142944e-01,1.142944e-01,-1.713488e-01,1.713488e-01,-1.103035e-01,1.103035e-01,-1.045221e-01,1.045221e-01,-1.293015e-01,1.293015e-01,-9.763183e-02,9.763183e-02,-1.387213e-01,1.387213e-01,-9.031167e-02,9.031167e-02,-1.283052e-01,1.283052e-01,-1.133462e-01,1.133462e-01,-9.370681e-02,9.370681e-02,-1.079269e-01,1.079269e-01,-1.331913e-01,1.331913e-01,-8.969902e-02,8.969902e-02,-1.044560e-01,1.044560e-01,-9.387466e-02,9.387466e-02,-1.208988e-01,1.208988e-01,-1.252011e-01,1.252011e-01,-1.401277e-01,1.401277e-01,-1.461381e-01,1.461381e-01,-1.323763e-01,1.323763e-01,-9.923889e-02,9.923889e-02,-1.142899e-01,1.142899e-01,-9.110853e-02,9.110853e-02,-1.106607e-01,1.106607e-01,-1.253140e-01,1.253140e-01,-9.657895e-02,9.657895e-02,-1.030010e-01,1.030010e-01,-1.348857e-01,1.348857e-01,-1.237793e-01,1.237793e-01,-1.296943e-01,1.296943e-01,-1.323385e-01,1.323385e-01,-8.331554e-02,8.331554e-02,-8.417589e-02,8.417589e-02,-1.104431e-01,1.104431e-01,-1.170710e-01,1.170710e-01,-1.391725e-01,1.391725e-01,-1.485189e-01,1.485189e-01,-1.840393e-01,1.840393e-01,-1.238250e-01,1.238250e-01,-1.095287e-01,1.095287e-01,-1.177869e-01,1.177869e-01,-1.036409e-01,1.036409e-01,-9.802581e-02,9.802581e-02,-9.364054e-02,9.364054e-02,-9.936022e-02,9.936022e-02,-1.117201e-01,1.117201e-01,-1.081300e-01,1.081300e-01,-1.331861e-01,1.331861e-01,-1.192122e-01,1.192122e-01,-9.889761e-02,9.889761e-02,-1.173456e-01,1.173456e-01,-1.032917e-01,1.032917e-01,-9.268551e-02,9.268551e-02,-1.178563e-01,1.178563e-01,-1.215065e-01,1.215065e-01,-1.060437e-01,1.060437e-01,-1.010044e-01,1.010044e-01,-1.021683e-01,1.021683e-01,-9.974968e-02,9.974968e-02,-1.161528e-01,1.161528e-01,-8.686721e-02,8.686721e-02,-8.145259e-02,8.145259e-02,-9.937060e-02,9.937060e-02,-1.170885e-01,1.170885e-01,-7.693779e-02,7.693779e-02,-9.047233e-02,9.047233e-02,-9.168442e-02,9.168442e-02,-1.054105e-01,1.054105e-01,-9.036177e-02,9.036177e-02,-1.251949e-01,1.251949e-01,-9.523847e-02,9.523847e-02,-1.038930e-01,1.038930e-01,-1.433660e-01,1.433660e-01,-1.489830e-01,1.489830e-01,-8.393174e-02,8.393174e-02,-8.888026e-02,8.888026e-02,-9.347861e-02,9.347861e-02,-1.044838e-01,1.044838e-01,-1.102144e-01,1.102144e-01,-1.383415e-01,1.383415e-01,-1.466476e-01,1.466476e-01,-1.129741e-01,1.129741e-01,-1.310915e-01,1.310915e-01,-1.070648e-01,1.070648e-01,-7.559007e-02,7.559007e-02,-8.812082e-02,8.812082e-02,-1.234272e-01,1.234272e-01,-1.088022e-01,1.088022e-01,-8.388703e-02,8.388703e-02,-7.179593e-02,7.179593e-02,-1.008961e-01,1.008961e-01,-9.030070e-02,9.030070e-02,-8.581345e-02,8.581345e-02,-9.023431e-02,9.023431e-02,-9.807321e-02,9.807321e-02,-9.621402e-02,9.621402e-02,-1.730195e-01,1.730195e-01,-8.984631e-02,8.984631e-02,-9.556661e-02,9.556661e-02,-1.047576e-01,1.047576e-01,-7.854313e-02,7.854313e-02,-8.682118e-02,8.682118e-02,-1.159761e-01,1.159761e-01,-1.339540e-01,1.339540e-01,-1.003048e-01,1.003048e-01,-9.747544e-02,9.747544e-02,-9.501058e-02,9.501058e-02,-1.321566e-01,1.321566e-01,-9.194706e-02,9.194706e-02,-9.359276e-02,9.359276e-02,-1.015916e-01,1.015916e-01,-1.174192e-01,1.174192e-01,-1.039931e-01,1.039931e-01,-9.746733e-02,9.746733e-02,-1.286120e-01,1.286120e-01,-1.044899e-01,1.044899e-01,-1.066385e-01,1.066385e-01,-8.368626e-02,8.368626e-02,-1.271919e-01,1.271919e-01,-1.055946e-01,1.055946e-01,-8.272876e-02,8.272876e-02,-1.370564e-01,1.370564e-01,-8.539379e-02,8.539379e-02,-1.100343e-01,1.100343e-01,-8.102170e-02,8.102170e-02,-1.028728e-01,1.028728e-01,-1.305065e-01,1.305065e-01,-1.059506e-01,1.059506e-01,-1.264646e-01,1.264646e-01,-8.383843e-02,8.383843e-02,-9.357698e-02,9.357698e-02,-7.474400e-02,7.474400e-02,-7.814045e-02,7.814045e-02,-8.600970e-02,8.600970e-02,-1.206090e-01,1.206090e-01,-9.986512e-02,9.986512e-02,-8.516476e-02,8.516476e-02,-7.198783e-02,7.198783e-02,-7.838409e-02,7.838409e-02,-1.005142e-01,1.005142e-01,-9.951857e-02,9.951857e-02,-7.253998e-02,7.253998e-02,-9.913739e-02,9.913739e-02,-7.500360e-02,7.500360e-02,-9.258090e-02,9.258090e-02,-1.400287e-01,1.400287e-01,-1.044404e-01,1.044404e-01,-7.404339e-02,7.404339e-02,-7.256833e-02,7.256833e-02,-1.006995e-01,1.006995e-01,-1.426043e-01,1.426043e-01,-1.036529e-01,1.036529e-01,-1.208443e-01,1.208443e-01,-1.074245e-01,1.074245e-01,-1.141448e-01,1.141448e-01,-1.015809e-01,1.015809e-01,-1.028822e-01,1.028822e-01,-1.055682e-01,1.055682e-01,-9.468699e-02,9.468699e-02,-1.010098e-01,1.010098e-01,-1.205054e-01,1.205054e-01,-8.392956e-02,8.392956e-02,-8.052297e-02,8.052297e-02,-9.576507e-02,9.576507e-02,-9.515692e-02,9.515692e-02,-1.564745e-01,1.564745e-01,-7.357238e-02,7.357238e-02,-1.129262e-01,1.129262e-01,-1.013265e-01,1.013265e-01,-8.760761e-02,8.760761e-02,-8.714771e-02,8.714771e-02,-9.605039e-02,9.605039e-02,-9.064677e-02,9.064677e-02,-8.243857e-02,8.243857e-02,-8.495858e-02,8.495858e-02,-8.350249e-02,8.350249e-02,-7.423234e-02,7.423234e-02,-7.930799e-02,7.930799e-02,-6.620023e-02,6.620023e-02,-7.311919e-02,7.311919e-02,-1.237938e-01,1.237938e-01,-1.086814e-01,1.086814e-01,-6.379798e-02,6.379798e-02,-7.526021e-02,7.526021e-02,-8.297097e-02,8.297097e-02,-8.186337e-02,8.186337e-02,-7.627362e-02,7.627362e-02,-1.061638e-01,1.061638e-01,-8.328494e-02,8.328494e-02,-1.040895e-01,1.040895e-01,-7.649056e-02,7.649056e-02,-7.299058e-02,7.299058e-02,-9.195198e-02,9.195198e-02,-7.990880e-02,7.990880e-02,-7.429346e-02,7.429346e-02,-9.991702e-02,9.991702e-02,-9.755385e-02,9.755385e-02,-1.344138e-01,1.344138e-01,-1.707917e-01,1.707917e-01,-8.325450e-02,8.325450e-02,-8.137793e-02,8.137793e-02,-8.308659e-02,8.308659e-02,-7.440414e-02,7.440414e-02,-7.012744e-02,7.012744e-02,-8.122943e-02,8.122943e-02,-8.845462e-02,8.845462e-02,-8.803450e-02,8.803450e-02,-9.653392e-02,9.653392e-02,-8.795691e-02,8.795691e-02,-1.119045e-01,1.119045e-01,-1.068308e-01,1.068308e-01,-8.406359e-02,8.406359e-02,-1.220414e-01,1.220414e-01,-1.024235e-01,1.024235e-01,-1.252897e-01,1.252897e-01,-1.121234e-01,1.121234e-01,-9.054150e-02,9.054150e-02,-8.974435e-02,8.974435e-02,-1.351578e-01,1.351578e-01,-1.106442e-01,1.106442e-01,-8.093913e-02,8.093913e-02,-9.800762e-02,9.800762e-02,-7.012823e-02,7.012823e-02,-7.434949e-02,7.434949e-02,-8.684816e-02,8.684816e-02,-8.916388e-02,8.916388e-02,-8.773159e-02,8.773159e-02,-7.709608e-02,7.709608e-02,-7.230518e-02,7.230518e-02,-9.662156e-02,9.662156e-02,-7.957632e-02,7.957632e-02,-7.628441e-02,7.628441e-02,-8.050202e-02,8.050202e-02,-1.290593e-01,1.290593e-01,-9.246182e-02,9.246182e-02,-9.703662e-02,9.703662e-02,-7.866445e-02,7.866445e-02,-1.064783e-01,1.064783e-01,-1.012339e-01,1.012339e-01,-6.828389e-02,6.828389e-02,-1.005039e-01,1.005039e-01,-7.559687e-02,7.559687e-02,-6.359878e-02,6.359878e-02,-8.387002e-02,8.387002e-02,-7.851323e-02,7.851323e-02,-8.878569e-02,8.878569e-02,-7.767654e-02,7.767654e-02,-8.033338e-02,8.033338e-02,-9.142797e-02,9.142797e-02,-8.590585e-02,8.590585e-02,-1.052318e-01,1.052318e-01,-8.760062e-02,8.760062e-02,-9.222192e-02,9.222192e-02,-7.548828e-02,7.548828e-02,-8.003344e-02,8.003344e-02,-1.177076e-01,1.177076e-01,-1.064964e-01,1.064964e-01,-8.655553e-02,8.655553e-02,-9.418112e-02,9.418112e-02,-7.248163e-02,7.248163e-02,-7.120974e-02,7.120974e-02,-6.393114e-02,6.393114e-02,-7.997487e-02,7.997487e-02,-1.220941e-01,1.220941e-01,-9.892518e-02,9.892518e-02,-8.270271e-02,8.270271e-02,-1.069400e-01,1.069400e-01,-5.860771e-02,5.860771e-02,-9.126600e-02,9.126600e-02,-6.212559e-02,6.212559e-02,-9.397538e-02,9.397538e-02,-8.070447e-02,8.070447e-02,-8.415587e-02,8.415587e-02,-8.564455e-02,8.564455e-02,-7.791811e-02,7.791811e-02,-6.642259e-02,6.642259e-02,-8.266167e-02,8.266167e-02,-1.134986e-01,1.134986e-01,-1.045267e-01,1.045267e-01,-7.122085e-02,7.122085e-02,-7.979415e-02,7.979415e-02,-7.922347e-02,7.922347e-02,-9.003421e-02,9.003421e-02,-8.796449e-02,8.796449e-02,-7.933279e-02,7.933279e-02,-8.307947e-02,8.307947e-02,-8.946349e-02,8.946349e-02,-7.643384e-02,7.643384e-02,-7.818534e-02,7.818534e-02,-7.990991e-02,7.990991e-02,-9.885664e-02,9.885664e-02,-8.071329e-02,8.071329e-02,-6.952112e-02,6.952112e-02,-6.429706e-02,6.429706e-02,-6.307229e-02,6.307229e-02,-8.100137e-02,8.100137e-02,-7.693623e-02,7.693623e-02,-6.906625e-02,6.906625e-02,-7.390462e-02,7.390462e-02,-6.487217e-02,6.487217e-02,-1.233681e-01,1.233681e-01,-6.979273e-02,6.979273e-02,-8.358669e-02,8.358669e-02,-1.095420e-01,1.095420e-01,-8.519717e-02,8.519717e-02,-7.599857e-02,7.599857e-02,-6.042816e-02,6.042816e-02,-6.546304e-02,6.546304e-02,-1.016245e-01,1.016245e-01,-8.308787e-02,8.308787e-02,-7.385708e-02,7.385708e-02,-6.751630e-02,6.751630e-02,-9.036695e-02,9.036695e-02,-9.371335e-02,9.371335e-02,-1.116088e-01,1.116088e-01,-5.693741e-02,5.693741e-02,-6.383983e-02,6.383983e-02,-5.389843e-02,5.389843e-02,-8.383191e-02,8.383191e-02,-7.820822e-02,7.820822e-02,-7.067557e-02,7.067557e-02,-7.971948e-02,7.971948e-02,-7.360668e-02,7.360668e-02,-7.008027e-02,7.008027e-02,-8.013378e-02,8.013378e-02,-8.331605e-02,8.331605e-02,-7.145702e-02,7.145702e-02,-7.863940e-02,7.863940e-02,-6.992679e-02,6.992679e-02,-5.716495e-02,5.716495e-02,-5.306006e-02,5.306006e-02,-8.855639e-02,8.855639e-02,-7.656397e-02,7.656397e-02,-6.939272e-02,6.939272e-02,-7.523742e-02,7.523742e-02,-8.472299e-02,8.472299e-02,-8.114341e-02,8.114341e-02,-6.795517e-02,6.795517e-02,-7.890130e-02,7.890130e-02,-7.488741e-02,7.488741e-02,-9.281972e-02,9.281972e-02,-9.325498e-02,9.325498e-02,-1.401587e-01,1.401587e-01,-1.176284e-01,1.176284e-01,-8.867597e-02,8.867597e-02,-8.124232e-02,8.124232e-02,-9.441235e-02,9.441235e-02,-8.029452e-02,8.029452e-02,-8.581848e-02,8.581848e-02,-1.029819e-01,1.029819e-01,-9.569118e-02,9.569118e-02,-7.690893e-02,7.690893e-02,-9.018228e-02,9.018228e-02,-1.049209e-01,1.049209e-01,-8.969413e-02,8.969413e-02,-8.651891e-02,8.651891e-02,-8.613331e-02,8.613331e-02,-7.120468e-02,7.120468e-02,-8.743959e-02,8.743959e-02,-7.607158e-02,7.607158e-02,-1.015547e-01,1.015547e-01,-8.090879e-02,8.090879e-02,-7.114079e-02,7.114079e-02,-8.744835e-02,8.744835e-02,-6.074904e-02,6.074904e-02,-6.919871e-02,6.919871e-02,-7.607774e-02,7.607774e-02,-9.444600e-02,9.444600e-02,-7.833429e-02,7.833429e-02,-6.817555e-02,6.817555e-02,-8.997390e-02,8.997390e-02,-9.845223e-02,9.845223e-02,-7.894180e-02,7.894180e-02,-7.921373e-02,7.921373e-02,-7.448032e-02,7.448032e-02,-1.178165e-01,1.178165e-01,-8.216686e-02,8.216686e-02,-8.103286e-02,8.103286e-02,-6.981470e-02,6.981470e-02,-8.709008e-02,8.709008e-02,-8.336259e-02,8.336259e-02,-6.213589e-02,6.213589e-02,-7.068045e-02,7.068045e-02,-6.915676e-02,6.915676e-02,-7.103416e-02,7.103416e-02,-6.523849e-02,6.523849e-02,-7.634760e-02,7.634760e-02,-7.263038e-02,7.263038e-02,-7.164396e-02,7.164396e-02,-8.745559e-02,8.745559e-02,-6.960181e-02,6.960181e-02,-8.500098e-02,8.500098e-02,-6.523260e-02,6.523260e-02,-7.319714e-02,7.319714e-02,-6.268125e-02,6.268125e-02,-7.083135e-02,7.083135e-02,-7.984517e-02,7.984517e-02,-1.256265e-01,1.256265e-01,-1.065412e-01,1.065412e-01,-8.524323e-02,8.524323e-02,-9.291364e-02,9.291364e-02,-7.936567e-02,7.936567e-02,-8.607723e-02,8.607723e-02,-7.583416e-02,7.583416e-02,-7.931928e-02,7.931928e-02,-7.408357e-02,7.408357e-02,-1.034404e-01,1.034404e-01,-1.012127e-01,1.012127e-01,-7.916689e-02,7.916689e-02,-8.753651e-02,8.753651e-02,-6.090366e-02,6.090366e-02,-7.500103e-02,7.500103e-02,-1.228709e-01,1.228709e-01,-6.318201e-02,6.318201e-02,-7.585420e-02,7.585420e-02,-7.089090e-02,7.089090e-02,-1.053542e-01,1.053542e-01,-8.549521e-02,8.549521e-02,-7.906308e-02,7.906308e-02,-6.338780e-02,6.338780e-02,-8.417910e-02,8.417910e-02,-7.115511e-02,7.115511e-02,-7.693949e-02,7.693949e-02,-7.446749e-02,7.446749e-02,-1.037929e-01,1.037929e-01,-7.991005e-02,7.991005e-02,-7.119439e-02,7.119439e-02,-7.071340e-02,7.071340e-02,-8.587362e-02,8.587362e-02,-7.001236e-02,7.001236e-02,-7.567115e-02,7.567115e-02,-7.118930e-02,7.118930e-02,-6.844895e-02,6.844895e-02,-1.035118e-01,1.035118e-01,-8.156618e-02,8.156618e-02,-7.449593e-02,7.449593e-02,-8.154360e-02,8.154360e-02,-9.110878e-02,9.110878e-02,-6.222534e-02,6.222534e-02,-1.033841e-01,1.033841e-01,-6.811687e-02,6.811687e-02,-6.828443e-02,6.828443e-02,-5.769408e-02,5.769408e-02,-5.917684e-02,5.917684e-02,-8.358868e-02,8.358868e-02]}]};

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
