//requires: ccv, numeric.js

var clm = {
	tracker : function(params) {
    
    if (!params) params = {};
		if (params.constantVelocity === undefined) params.constantVelocity = true;
		if (params.searchWindow === undefined) params.searchWindow = 10;
	
		var numPatches, patchSize, numParameters, patchType;
		var gaussianPD;
		var eigenVectors, eigenValues;
		var sketchCC, sketchW, sketchH;
		var candidate;
		
		var currentParameters = [];
		var currentPositions = [];
		var previousParameters = [];
		
		var weightMatrices = [];
		var weightMatricesOld = [];
		var weights = [];
		
		var patches = [];
		var responses = [];
		var meanShape = [];
		
		/*
		It's possible to experiment with the sequence of variances used for the finding the maximum in the KDE.
		This sequence is pretty arbitrary, but was found to be okay using some manual testing.
		In Saragih's paper he used the sequence [20,10,5,1], this was however found to be slower and equally precise.
		*/
		//var varianceSeq = [10,7,5];
		var varianceSeq = [3,1.5,0.75];
		
		/*
		The PDM variance determines how the PDM models parameters can vary when fitting.
		A low variance puts very little constraint on the parameters, a high variance puts a lot of constraint.
		According to the formula in Saragih's paper, this parameter should be around 0.00005, however a higher variance was found to work better for this model.
		*/
		var PDMVariance = 0.05;
		
		var relaxation = 0.1;
		
		var first = true;
		
		var convergenceLimit = 0.01;
		
		var learningRate = [];
		var stepParameter = 1.25;
		var prevCostFunc = []
		
		var searchWindow;
		var modelWidth;
		
		var halfSearchWindow, vecProbs, responsePixels;
		var updatePosition = new Float64Array(2);
		var vecpos = new Float64Array(2);
		
		if (pModel.hints && mosseFilter && left_eye_filter && right_eye_filter && nose_filter) {
		  var mossef_lefteye = new mosseFilter();
      mossef_lefteye.load(left_eye_filter);
      var mossef_righteye = new mosseFilter();
      mossef_righteye.load(right_eye_filter);
      var mossef_nose = new mosseFilter();
      mossef_nose.load(nose_filter);
      var mossef_face = new mosseFilter();
      mossef_face.load(face_filter);
      
      var right_eye_position = [0.0,0.0];
      var left_eye_position = [0.0,0.0];
      var nose_position = [0.0,0.0];
      var face_position = [0.0,0.0];
      var lep, rep, mep;
		} else {
      console.log("MOSSE filters not found, using rough approximation for initialization.");
    }
    var facecheck_count = 0;
    var face_peak = [];
    var face_diff = [];
		
		var webglFi, mosseCalc;
		
		var gidTime1, gidTime2;
		
		this.init = function(canvas) {
			// do all prep stuff
			
			sketchCC = canvas.getContext('2d');
			sketchW = canvas.width;
			sketchH = canvas.height;
			
			//console.log('starting initialization');
			
			
			// load from model
			patchType = pModel.patchModel.patchType;
			numPatches = pModel.patchModel.numPatches;
			patchSize = pModel.patchModel.patchSize[0];
			if (patchType == "MOSSE") {
			  searchWindow = patchSize;
			} else {
			  searchWindow = params.searchWindow;
			}
			numParameters = pModel.shapeModel.numEvalues;
			modelWidth = pModel.patchModel.canvasSize[0];
			
			// load eigenvectors
			eigenVectors = numeric.rep([numPatches*2,numParameters],0.0);
			for (var i = 0;i < numPatches*2;i++) {
				for (var j = 0;j < numParameters;j++) {
          eigenVectors[i][j] = pModel.shapeModel.eigenVectors[i][j];
				}
			}
			
			// load mean shape
			for (var i = 0; i < numPatches;i++) {
			  meanShape[i] = [pModel.shapeModel.meanShape[i][0], pModel.shapeModel.meanShape[i][1]];
			}
			
			// load eigenvalues
			eigenValues = pModel.shapeModel.eigenValues;
			
			// load patchweight matrices for comparing
			var weight;
			/*for (var w = 0;w < numPatches;w++) {
				weightMatricesOld[w] = numeric.rep([patchSize, patchSize],0);
				// insert
				for (var i = 0;i < patchSize;i++) {
					for (var j = 0;j < patchSize;j++) {
						//weight = pModel.patchModel.weights[w+(((i*patchSize)+j)*numPatches)];
						weight = pModel.patchModel.weights[w][(i*patchSize)+j];
						// cut off all values above 1 and below -1
						weight = weight > 1 ? 1 : weight;
						weight = weight < -1 ? -1 : weight;
						weightMatricesOld[w][i][j] = weight;
					}
				}
			}*/
			
			if (patchType == "SVM") {
        for (var w = 0;w < numPatches;w++) {
          for (var i = 0;i < (patchSize*patchSize);i++) {
            weight = pModel.patchModel.weights[w][i];
            weight = weight > 1 ? 1 : weight;
            weight = weight < -1 ? -1 : weight;
            weights[i +(w*patchSize*patchSize)] = weight;
          }
        }
      } else {
        weights = pModel.patchModel.weights;
      }
			
			// precalculate gaussianPriorDiagonal
			gaussianPD = numeric.rep([numParameters+4, numParameters+4],0);
			// set values and append manual inverse
			for (var i = 0;i < numParameters;i++) {
        gaussianPD[i+4][i+4] = 1/eigenValues[i];
			}
				
			for (var i = 0;i < numParameters+4;i++) {
				currentParameters[i] = 0;
			}
			
			// set up webgl filter calculation
			if (patchType == "SVM") {
			  webglFi = new webglFilter();
        webglFi.init(weights, numPatches, searchWindow+patchSize, searchWindow+patchSize, patchSize, patchSize, true);
			} else if (patchType == "MOSSE") {
			  mosseCalc = new mosseFilterResponses();
			  mosseCalc.init(weights, numPatches, patchSize, patchSize);
			}
			
			halfSearchWindow = searchWindow/2;
			responsePixels = searchWindow*searchWindow;
			vecProbs = new Float64Array(responsePixels);
			
			//console.log('ended initialization');
			for (var i = 0;i < numPatches;i++) {
			  learningRate[i] = 1.0;
		    prevCostFunc[i] = 0.0;
			}
		}
		
		var gaussianProb = function(mean, coordinate, variance) {
		  // calculate pdf gaussian probability
		  var dx = coordinate[0] - mean[0];
			var dy = coordinate[1] - mean[1];
			
			return Math.exp(-0.5*((dx*dx)+(dy*dy))/variance);
		}
		
		var createJacobian = function(parameters, eigenVectors) {
			// generates the jacobian matrix
      
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
		
		var calculatePositions = function(parameters, useTransforms) {
			var x, y, a, b;
			var numParameters = parameters.length;
			var positions = [];
			for (var i = 0;i < numPatches;i++) {
				x = pModel.shapeModel.meanShape[i][0];
				y = pModel.shapeModel.meanShape[i][1];
				for (var j = 0;j < numParameters-4;j++) {
					x += pModel.shapeModel.eigenVectors[(i*2)][j]*parameters[j+4];
					y += pModel.shapeModel.eigenVectors[(i*2)+1][j]*parameters[j+4];
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
		
		this.calculatePositions = function(parameters) {
		  return calculatePositions(parameters, true);
		}
		
		this.detectPosition = function(el) {
		  var canvas = document.createElement('canvas');
		  canvas.width = el.width;
		  canvas.height = el.height;
		  var cc = canvas.getContext('2d');
		  cc.drawImage(el, 0, 0, el.width, el.height);
		  
		  // do viola-jones on canvas to get initial guess, if we don't have any points
      /*var comp = ccv.detect_objects(
        ccv.grayscale(canvas), ccv.cascade, 5, 1
      );*/
      var jf = new jsfeat_face(canvas);
      var comp = jf.findFace(1);
      
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
		
		this.getCurrentPosition = function() {
		  return currentPositions;
		}
		
		this.getCurrentParameters = function() {
		  return currentParameters;
		}
		
		var gpopt = function(responseWidth, currentPositionsj, updatePosition, vecProbs, responses, opj0, opj1, i, j, variance) {
      var pos_idx = 0;
      var vpsum = 0;
      var dx, dy;
      for (var k = 0;k < responseWidth;k++) {
        updatePosition[1] = opj1+k;
        for (var l = 0;l < responseWidth;l++) {
          updatePosition[0] = opj0+l;
          
          dx = currentPositionsj[0] - updatePosition[0];
          dy = currentPositionsj[1] - updatePosition[1];
          
          //vecProbs[pos_idx] = responses[j][pos_idx] * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
          vecProbs[pos_idx] = responses[j][pos_idx] * Math.exp(-0.5*((dx*dx)+(dy*dy))/variance);
          vpsum += vecProbs[pos_idx];
          pos_idx++;
        }
      }
      
      return vpsum;
    }
    
    var gpopt2 = function(responseWidth, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1) {
      var pos_idx = 0;
      var vecsum = 0;
      vecpos[0] = 0;
      vecpos[1] = 0;
      for (var k = 0;k < responseWidth;k++) {
        updatePosition[1] = opj1+k;
        for (var l = 0;l < responseWidth;l++) {
          updatePosition[0] = opj0+l;
          vecsum = vecProbs[pos_idx]/vpsum;
          vecpos[0] += vecsum*updatePosition[0];
          vecpos[1] += vecsum*updatePosition[1];
          pos_idx++;
        }
      }
    }
		
		/*
     *  element : canvas or video element
     *  TODO: should be able to take img element as well
     */
		this.track = function(element) {
			
			// for timing only
			var startTime = (new Date).getTime();
			
			var scaling, translateX, translateY, rotation;
			var croppedPatches = [];
			var ptch, px, py, pw;
			
			if (!first) {
        // TODO: check here if the previous tracked position is good enough, every 100th frame?
        facecheck_count += 1;
        if (facecheck_count % 10 == 0) {
          //get centerpoint and approximate width
          var centerpoint = [];
          centerpoint[0] = currentPositions[62][0];
          centerpoint[1] = currentPositions[62][1];
          //get xmin and xmax
          var xmin = 10000000;
          var xmax = 0;
          for (var i = 0;i < currentPositions.length;i++) {
            if (currentPositions[i][0] < xmin) {
              xmin = currentPositions[i][0];
            }
            if (currentPositions[i][0] > xmax) {
              xmax = currentPositions[i][0];
            }
          }
          var modelwidth = xmax-xmin
          var face_result = mossef_face.track(element, Math.round(centerpoint[0]-(modelwidth/2)), Math.round(centerpoint[1]-(modelwidth/2)), modelwidth, modelwidth, false, false, true);
          var peak = mossef_face.peak_prev;
          var psr = mossef_face.psr_prev;
          // TODO : check if peak is within center of model?
          // or check if peak/sidelobe is over some threshold
          // else, set first = true
          
          var peak_avg = 0;
          face_peak.push(peak);
          if (face_peak.length > 5) {
            face_peak.splice(0,1);
            for (var i = 0;i < face_peak.length;i++) {
              peak_avg += face_peak[i]; 
            }
            peak_avg /= face_peak.length;
          }
          var diffx = centerpoint[0]-(face_result[0]+centerpoint[0]-(modelwidth/2));
          var diffy = centerpoint[1]-(face_result[1]+centerpoint[1]-(modelwidth/2));
          var diff = Math.sqrt(diffx*diffx + diffy*diffy);
          var diff_avg = 0;
          face_diff.push(diff);
          if (face_diff.length > 5) {
            face_diff.splice(0,1);
            for (var i = 0;i < face_diff.length;i++) {
              diff_avg += face_diff[i]; 
            }
            diff_avg /= face_diff.length;
          }
          
          //document.getElementById('peak').innerHTML = "peak average :"+peak_avg;
          //document.getElementById('psr').innerHTML = "diff :"+diff_avg;
          
          if ((face_peak.length > 5 && peak_avg < 0.10) || (face_diff.length && diff_avg > 15)) {
            first = true;
            face_diff = [];
            face_peak = [];
            for (var i = 0;i < currentParameters.length;i++) {
              currentParameters[i] = 0;
              previousParameters = [];
            }
          }
        }
			}

			
			if (first) {
				// do viola-jones on canvas to get initial guess, if we don't have any points
				
				var det = this.detectPosition(element);
				if (!det) {
				  // if no face found, stop.
				  return false;
				}
				
				// calculate modelWidth/height from meanshape
				var xmin, ymin, xmax, ymax;
        xmin = ymin = 1000000;
        xmax = ymax = 0;
        for (var i = 0;i < meanShape.length;i++) {
          xmin = meanShape[i][0] < xmin ? meanShape[i][0] : xmin;
          xmax = meanShape[i][0] > xmax ? meanShape[i][0] : xmax;
          ymin = meanShape[i][1] < ymin ? meanShape[i][1] : ymin;
          ymax = meanShape[i][1] > ymax ? meanShape[i][1] : ymax;
        }
        var modelwidth = xmax-xmin;
        var modelheight = ymax-ymin;
        
        if (pModel.hints && mosseFilter && left_eye_filter && right_eye_filter && nose_filter) {
				  var eyeFilterWidth = candidate.width * 4.5/10;
          var noseFilterWidth = candidate.width * 6/10;
          
          // detect position of eyes and nose via mosse filter
          var right_result = mossef_righteye.track(element, Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth, false);
          var left_result = mossef_lefteye.track(element, Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth, false);
          //var nose_result = mossef_nose.track(element, Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(3/4)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth, false);
          var nose_result = mossef_nose.track(element, Math.round(candidate.x+(candidate.width/2)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(5/8)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth, false);
          right_eye_position[0] = Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2))+right_result[0];
          right_eye_position[1] = Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2))+right_result[1];
          left_eye_position[0] = Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2))+left_result[0];
          left_eye_position[1] = Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2))+left_result[1];
          //nose_position[0] = Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2))+nose_result[0];
          //nose_position[1] = Math.round(candidate.y+candidate.height*(3/4)-(noseFilterWidth/2))+nose_result[1];
          nose_position[0] = Math.round(candidate.x+(candidate.width/2)-(eyeFilterWidth/2))+nose_result[0];
          nose_position[1] = Math.round(candidate.y+candidate.height*(5/8)-(eyeFilterWidth/2))+nose_result[1];
          
          /*var canvasContext = document.getElementById('overlay').getContext('2d')
          
          canvasContext.strokeRect(candidate.x, candidate.y, candidate.width, candidate.height);
          canvasContext.strokeRect(Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
          canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
          //canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(3/4)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth);
          canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(5/8)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
          
          //element.pause()
          canvasContext.fillStyle = "rgb(200,200,200)";
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
          canvasContext.fill();*/
          
          // get eye and nose positions of model
          var lep = pModel.hints.leftEye;
          var rep = pModel.hints.rightEye;
          var mep = pModel.hints.nose;
          
          // get scaling, rotation, etc. via procrustes analysis
          var procrustes_params = procrustes([left_eye_position, right_eye_position, nose_position], [lep, rep, mep]);
          translateX = procrustes_params[0];
          translateY = procrustes_params[1];
          scaling = procrustes_params[2];
          rotation = procrustes_params[3];
          
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
          var ccc = document.getElementById('overlay').getContext('2d');
          ccc.strokeRect(candidate.x,candidate.y,candidate.width,candidate.height);
          translateX = candidate.x-(xmin*scaling)+0.1*candidate.width;
          translateY = candidate.y-(ymin*scaling)+0.25*candidate.height;
          currentParameters[0] = scaling-1;
          currentParameters[2] = translateX;
          currentParameters[3] = translateY;
        }
        
        currentPositions = calculatePositions(currentParameters, true);
				
				first = false;
			} else {
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
				scaling = 1+currentParameters[0];
				translateX = currentParameters[2];
				translateY = currentParameters[3];
			}
			
			var secondTime = (new Date).getTime();
			//console.log("detectiontiming: "+(secondTime-startTime)+" ms")
			
			// copy canvas to a new dirty canvas
			sketchCC.save();
			
			// clear canvas
			sketchCC.clearRect(0, 0, sketchW, sketchH);
			
			sketchCC.scale(1/scaling, 1/scaling);
			sketchCC.rotate(-Math.asin(currentParameters[1]/scaling));
			sketchCC.translate(-translateX, -translateY);
			gidTime1 = (new Date).getTime();
			
			sketchCC.drawImage(element, 0, 0, element.width, element.height);
			
			gidTime2 = (new Date).getTime();
			sketchCC.restore();
			//	get cropped images around new points based on model parameters (not scaled and translated)
			var patchPositions = calculatePositions(currentParameters, false);
			
			//console.log("gidtime:"+(gidTime2-gidTime1));
			
			var pw, pl;
			if (patchType == "SVM") {
			  pw = pl = patchSize+searchWindow;
			} else {
			  pw = pl = searchWindow;
			}
			var pdataLength = pw*pl;
			
			var pdata, pmatrix, grayscaleColor;
			for (var i = 0; i < numPatches; i++) {
			  px = patchPositions[i][0]-(pw/2);
			  py = patchPositions[i][1]-(pl/2);
			  ptch = sketchCC.getImageData(px >> 0, py >> 0, pw, pl);
			  // add channels and convert to grayscale
			  // load patchdata to matrix
			  pdata = ptch.data;
			  
			  /*pmatrix = numeric.rep([pw, pl],0);
			  for (var j = 0;j < pdataLength;j++) {
			    grayscaleColor = pdata[j*4]*0.3 + pdata[(j*4)+1]*0.59 + pdata[(j*4)+2]*0.11;
			    pmatrix[j % pw][(j / pw) >> 0] = grayscaleColor;
			  }*/
			  
			  // TODO: preinitialize this
			  pmatrix = new Float64Array(pdataLength);
			  for (var j = 0;j < pdataLength;j++) {
			    grayscaleColor = pdata[j*4]*0.3 + pdata[(j*4)+1]*0.59 + pdata[(j*4)+2]*0.11;
			    pmatrix[j] = grayscaleColor;
			  }
			  patches[i] = pmatrix;
			}
			
			/*testing of printing weights*/
			/*sketchCC.clearRect(0, 0, sketchW, sketchH);
			for (var i = 0; i < numPatches; i++) {
			  px = patchPositions[i][0]-(patchSize/2);
			  py = patchPositions[i][1]-(patchSize/2);
			  var psci = sketchCC.createImageData(patchSize, patchSize);
			  var pscidata = psci.data;
			  for (var j = 0;j < (patchSize)*(patchSize);j++) {
			    //var val = weightMatricesOld[i].getValueAt(j % (patchSize), (j / (patchSize)) >> 0);
			    var val = weights[(i*patchSize*patchSize)+j];
			    val = (val*2000)+127;
			    val = val > 255 ? 255 : val;
			    val = val < 0 ? 0 : val;
			    pscidata[j*4] = val;
			    pscidata[(j*4)+1] = val;
			    pscidata[(j*4)+2] = val;
			    pscidata[(j*4)+3] = 255;
			  }
			  sketchCC.putImageData(psci, px >> 0, py >> 0);
			}*/
			
			/*testing of printing patches*/
			/*sketchCC.clearRect(0, 0, sketchW, sketchH); 
			for (var i = 0; i < numPatches; i++) {
			  px = patchPositions[i][0]-(pw/2);
			  py = patchPositions[i][1]-(pl/2);
			  var psci = sketchCC.createImageData(pw ,pl);
			  var pscidata = psci.data;
			  for (var j = 0;j < (pw)*(pl);j++) {
			    var val = patches[i][(((j / (pw)) >> 0)*(pw))+(j % (pl))];
			    //var val = patches[i].getValueAt(j % (pw), (j / (pl)) >> 0);
			    pscidata[j*4] = val;
			    pscidata[(j*4)+1] = val;
			    pscidata[(j*4)+2] = val;
			    pscidata[(j*4)+3] = 255;
			  }
			  sketchCC.putImageData(psci, px >> 0, py >> 0);
			}*/
			//console.log("starting response calculations");
			//sketchCC.clearRect(0, 0, sketchW, sketchH);
			
			var thirdTime = (new Date).getTime();
			//console.log("getting patches: "+(thirdTime-secondTime)+" ms")
			
			/* start webgl work here */
			/*var patchResponse, submatrix, submsum;
			for (var i = 0;i < numPatches;i++) {
			  // calculate response map
			  console.log("response: "+i);
			  patchResponse = new goog.math.Matrix(searchWindow, searchWindow);
			  for (var k = 0;k < searchWindow;k++) {
			    for (var l = 0;l < searchWindow;l++) {
			      // get submatrix for this patch
			      submatrix = new goog.math.Matrix(patchSize, patchSize);
			      goog.math.Matrix.forEach(submatrix, function(value, a, b) {
			        submatrix.setValueAt(a,b, this.getValueAt(a+k,b+l));
			      }, patches[i]);
            // normalize matrix so that mean pixel value is 0 and variance of 1
            normalizeMatrix(submatrix);
			      // multiply by this patches svm
			      submsum = 0;
			      goog.math.Matrix.forEach(submatrix, function(value, a, b) {
			        submsum += value * this.getValueAt(a,b);
			      }, weightMatrices[i])
			      // calculate exponential stuff
			      //submsum = 1/(1 + Math.exp(submsum));
			      
			      // insert to matrix
			      patchResponse.setValueAt(k,l,submsum);
			    }
			  }
			  
			  goog.math.Matrix.forEach(patchResponse, function(entry, a, b) {
			    patchResponse.setValueAt(a,b, 1-(1/(1+ Math.exp(entry))));
			  });
			  
			  // TODO : is this normalization correct?
			  // try using standard deviation to figure out gain
			  normalizeFilterMatrix(patchResponse);
			  
			  responses[i] = patchResponse;
      }*/
			/* end webgl work here */
			
			if (patchType == "SVM") {
        var responses = webglFi.getResponses(patches, true);
      } else if (patchType == "MOSSE") {
        var responses = mosseCalc.getResponses(patches);
      }
      
      // change responses to matrix form
      /*var respi;
      var responseSize = searchWindow*searchWindow;
      for (var i = 0;i < numPatches;i++) {
        respi = new goog.math.Matrix(searchWindow, searchWindow);
			  for (var j = 0;j < responseSize;j++) {
			    respi.setValueAt((j / searchWindow) >> 0, j % searchWindow, responses[i][j]);
        }
        //responses[i] = respi.getTranspose();
        responses[i] = respi;
      }*/
      
			var secondTime = (new Date).getTime();
			//console.log("calculating responses: "+(secondTime - thirdTime)+" ms")
			
			//testing of printing patches
      /*var canvasppt = document.createElement('canvas')
      canvasppt.setAttribute('width', searchWindow+"px");
      canvasppt.setAttribute('height', (searchWindow*numPatches)+"px");
      canvasppt.setAttribute('id', 'patchprinttest');
      document.body.appendChild(canvasppt);
      var pptcc = canvasppt.getContext('2d');
			for (var i = 0; i < numPatches; i++) {
			  var psci = pptcc.createImageData(searchWindow ,searchWindow );
			  var pscidata = psci.data;
			  for (var j = 0;j < (searchWindow )*(searchWindow );j++) {
			    //var val = responses[i].getValueAt(j % searchWindow, (j / searchWindow) >> 0)*255;
			    var val = responses[i][((j % searchWindow)*searchWindow) + ((j / searchWindow) >> 0)]*255;
			    pscidata[j*4] = val;
			    pscidata[(j*4)+1] = val;
			    pscidata[(j*4)+2] = val;
			    pscidata[(j*4)+3] = 255;
			  }
			  pptcc.putImageData(psci, 0, searchWindow*i);
			}*/
			
			// print responses
			/*sketchCC.clearRect(0, 0, sketchW, sketchH);
			for (var i = 0; i < numPatches; i++) {
			  px = patchPositions[i][0]-searchWindow/2;
			  py = patchPositions[i][1]-searchWindow/2;
			  var psci = sketchCC.createImageData(searchWindow,searchWindow);
			  var pscidata = psci.data;
			  for (var j = 0;j < (searchWindow)*(searchWindow);j++) {
			    //var val = responses[i].getValueAt(j % searchWindow, (j / searchWindow) >> 0);
			    var val = responses[i][((j % searchWindow)*searchWindow) + ((j / searchWindow) >> 0)]
			    val = (val)*255;
			    val = val > 255 ? 255 : val;
			    val = val < 0 ? 0 : val;
			    pscidata[j*4] = val;
			    pscidata[(j*4)+1] = val;
			    pscidata[(j*4)+2] = val;
			    pscidata[(j*4)+3] = 255;
			  }
			  sketchCC.putImageData(psci, px >> 0, py >> 0);
			}*/
			
			// iterate until convergence or max 10, 20 iterations?:
			var originalPositions = currentPositions;
			var jac;
			var meanshiftVectors = [];
			
			var dootholder = 0;
			var partbholder = 0;
			var partaholder = 0;
			
			for (var i = 0; i < varianceSeq.length; i++) {
				
				var partastart = (new Date).getTime();
				
				// calculate jacobian
				jac = createJacobian(currentParameters, eigenVectors);

				/* continue webgl work here? */
				//debugging
				//var debugMVs = [];
				//
				
				var partaend = (new Date).getTime();
				partaholder += (partaend-partastart);
				
				var opj0, opj1;
				
				for (var j = 0;j < numPatches;j++) {
				  // for every point in each response:
					// calculate PI x gaussian
				  
          //var boot = (new Date).getTime();
          
          opj0 = originalPositions[j][0]-halfSearchWindow;
          opj1 = originalPositions[j][1]-halfSearchWindow;
          
				  /*goog.math.Matrix.forEach(responses[j], function(value, a, b) {
				    var pos = (searchWindow*b)+a;
				    var updatePosition = [opj0+b, opj1+a];
				    vecProbs[pos] = value * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
				  });*/
				  
          /*var pos_idx = 0;
          for (var k = 0;k < searchWindow;k++) {
            updatePosition[1] = opj1+k;
            for (var l = 0;l < searchWindow;l++) {
              updatePosition[0] = opj0+l;
              
              vecProbs[pos_idx] = responses[j][pos_idx] * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
              pos_idx++;
            }
          }*/
				  
				  //var doot = (new Date).getTime();
          //dootholder += (doot-boot);
				  
				  // sum PI x gaussians 
				  /*var vpsum = 0;
				  for (var k = 0;k < responsePixels;k++) {
				    vpsum += vecProbs[k];
				  }*/
				  
				  var vpsum = gpopt(searchWindow, currentPositions[j], updatePosition, vecProbs, responses, opj0, opj1, i, j, varianceSeq[i]);
				  
				  //debugging
				  //var debugMatrixMV = numeric.rep([searchWindow,searchWindow],0.0);
				  //
				  
          /*var pos_idx = 0;
				  var vecsum = 0;
				  vecpos[0] = 0;
				  vecpos[1] = 0;
				  for (var k = 0;k < searchWindow;k++) {
				    updatePosition[1] = opj1+k;
				    for (var l = 0;l < searchWindow;l++) {
				      
				      updatePosition[0] = opj0+l;
				      vecsum = vecProbs[pos_idx]/vpsum;
				      
				      //debugging
				      //debugMatrixMV[k][l] = vecsum;
				      //
				      vecpos[0] += vecsum*updatePosition[0];
				      vecpos[1] += vecsum*updatePosition[1];
				      pos_idx++;
				    }
				  }*/
				  
				  gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1);
				  
				  // evaluate here whether to increase/decrease stepSize
				  /*if (vpsum >= prevCostFunc[j]) {
				    learningRate[j] *= stepParameter;
				  } else {
				    learningRate[j] = 1.0;
				  }
				  prevCostFunc[j] = vpsum;*/
				  
				  // compute mean shift vectors
				  // extrapolate meanshiftvectors
				  //var msv = [];
				  //msv[0] = learningRate[j]*(vecpos[0] - currentPositions[j][0]);
				  //msv[1] = learningRate[j]*(vecpos[1] - currentPositions[j][1]);
				  //meanshiftVectors[j] = msv;
				  meanshiftVectors[j] = [vecpos[0] - currentPositions[j][0], vecpos[1] - currentPositions[j][1]];
				  
				  //if (isNaN(msv[0]) || isNaN(msv[1])) debugger;
				  
				  //debugging
				  //debugMVs[j] = debugMatrixMV;
				  //
				}
				
				var partbstart = (new Date).getTime();
				
				// plot the meanshiftvector to see if it's correct
        /*sketchCC.clearRect(0, 0, sketchW, sketchH);
        for (var k = 0; k < numPatches; k++) {
          px = patchPositions[k][0]-searchWindow/2;
          py = patchPositions[k][1]-searchWindow/2;
          var psci = sketchCC.createImageData(searchWindow,searchWindow);
          var pscidata = psci.data;
          for (var j = 0;j < (searchWindow)*(searchWindow);j++) {
            var val = debugMVs[k][(j / searchWindow) >> 0][j % searchWindow];
            val = val*255*500;
            val = val > 255 ? 255 : val;
            val = val < 0 ? 0 : val;
            pscidata[j*4] = val;
            pscidata[(j*4)+1] = val;
            pscidata[(j*4)+2] = val;
            pscidata[(j*4)+3] = 255;
          }
          sketchCC.putImageData(psci, px >> 0, py >> 0);
        }*/
				
				/* end webgl work here */
				var meanShiftVector = numeric.rep([numPatches*2, 1],0.0);
				for (var k = 0;k < numPatches;k++) {
				  meanShiftVector[k*2][0] = meanshiftVectors[k][0];
				  meanShiftVector[(k*2)+1][0] = meanshiftVectors[k][1];
				}
				
				
				//test : draw meanshiftvector diagram
				/*var testMV = [];
				for (var k = 0;k < numPatches;k++) {
				  testMV[k] = []
				  testMV[k][0] = currentPositions[k][0] + meanShiftVector.getValueAt(k*2, 0);
				  testMV[k][1] = currentPositions[k][1] + meanShiftVector.getValueAt((k*2)+1, 0);
				}*/
				
				// compute pdm parameter update
				//var prior = gaussianPD.multiply(PDMVariance);
				var prior = numeric.mul(gaussianPD, varianceSeq[i]);
				var jtj = numeric.dot(numeric.transpose(jac), jac);
				var cpMatrix = numeric.rep([numParameters+4, 1],0.0);
        for (var l = 0;l < (numParameters+4);l++) {
          cpMatrix[l][0] = currentParameters[l];
        }
				var priorP = numeric.dot(prior, cpMatrix);
        var jtv = numeric.dot(numeric.transpose(jac), meanShiftVector);
        var paramUpdateLeft = numeric.add(prior, jtj);
        var paramUpdateRight = numeric.sub(priorP, jtv);
				var paramUpdate = numeric.dot(numeric.inv(paramUpdateLeft), paramUpdateRight);
				
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
				
				var partbend = (new Date).getTime();
				partbholder += (partbend-partbstart);
				
			}
			
			if (params.constantVelocity) {
        // add current parameter to array of previous parameters
        previousParameters.push(currentParameters.slice());
        previousParameters.splice(0, previousParameters.length == 3 ? 1 : 0);
			}
			
			//console.log("doot: "+(dootholder)+" ms")
			//console.log("partaholder: "+(partaholder)+" ms")
			//console.log("partbholder: "+(partbholder)+" ms")
			
      document.getElementById('overlay').getContext('2d').clearRect(0, 0, 720, 576);
			//debugger;
			
			this.draw(document.getElementById('overlay'), currentParameters);
			
			var thirdTime = (new Date).getTime();
			//console.log("optimization of parameters: "+(thirdTime - secondTime)+" ms")
			//console.log("all time: "+(thirdTime - startTime)+" ms")
			
			// return new points
			return currentPositions;
		}
		
		var normalizeMatrix = function(matrix) {
      // normalize matrix 
      // only used in non-webGL code
		  
		  var entry, mean, sd, entries;
		  
		  var msize = matrix.getSize();
		  var msum = 0;
		  var msumsq = 0;
		  for (var i = 0;i < msize.height;i++) {
		    for (var j = 0;j < msize.width;j++) {
		      entry = matrix.getValueAt(i,j);
		      msum += entry;
		      msumsq += entry * entry;
		    }
		  }
		  entries = msize.height * msize.width;
		  mean = msum / entries;
		  sd = Math.sqrt((msumsq/entries) - (mean * mean));
		  
		  //var normMatrix = new goog.math.Matrix(msize.height, msize.width);
		  for (var i = 0;i < msize.height;i++) {
		    for (var j = 0;j < msize.width; j++) {
		      entry = matrix.getValueAt(i,j);
		      matrix.setValueAt(i,j, (entry-mean)/sd);
		      //normMatrix.setValueAt(i,j, (entry-mean)/sd);
		    }
		  }
		  
		  //return normMatrix;
		}
		
		var normalizeFilterMatrix = function(matrix) {
		  // normalize filter matrix 
      // only used in non-webGL code
		  
		  var entry;
		  
		  var msize = matrix.getSize();
		  var max = 0;
		  var min = 1;
		  
		  for (var i = 0;i < msize.height;i++) {
		    for (var j = 0;j < msize.width;j++) {
		      entry = matrix.getValueAt(i,j);
		      max = entry > max ? entry : max;
		      min = entry < min ? entry : min;
		    }
		  }
		  var dist = max-min;
		  
		  //var normMatrix = new goog.math.Matrix(msize.height, msize.width);
		  for (var i = 0;i < msize.height;i++) {
		    for (var j = 0;j < msize.width; j++) {
		      entry = matrix.getValueAt(i,j);
		      matrix.setValueAt(i,j, (entry-min)/dist);
		      //normMatrix.setValueAt(i,j, (entry-mean)/sd);
		    }
		  }
		}
		
		var drawPath = function(canvasContext, path, dp) {
			canvasContext.beginPath();
			var i, x, y, a, b;
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i/2][0];
				y = pModel.shapeModel.meanShape[i/2][1];
				for (var j = 0;j < numParameters;j++) {
					x += pModel.shapeModel.eigenVectors[i][j]*dp[j+4];
					y += pModel.shapeModel.eigenVectors[i+1][j]*dp[j+4];
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
		
		function drawPoint(canvasContext, point, dp) {
		  var i, x, y, a, b;
		  i = point*2;
			x = pModel.shapeModel.meanShape[i/2][0];
      y = pModel.shapeModel.meanShape[i/2][1];
			for (var j = 0;j < numParameters;j++) {
				x += pModel.shapeModel.eigenVectors[i][j]*dp[j+4];
				y += pModel.shapeModel.eigenVectors[i+1][j]*dp[j+4];
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
		
		this.draw = function(canvas, pv) {
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
			cc.save();
			
			var paths = pModel.path.normal;
			for (var i = 0;i < paths.length;i++) {
			  if (typeof(paths[i]) == 'number') {
			    drawPoint(cc, paths[i], params);
			  } else {
			    drawPath(cc, paths[i], params);
			  }
			}
			
			cc.restore()
		}
		
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
		
		return true
	}
}