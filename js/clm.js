  //requires: ccv, closure matrix library

var clm = {
	tracker : function(params) {
    
    if (!params) params = {};
		if (params.constantVelocity === undefined) params.constantVelocity = false;
		if (params.searchWindow === undefined) params.searchWindow = 10;
	
		var numPatches, patchSize, numParameters;
		var gaussianPD;
		var eigenVectors, eigenValues;
		var sketchCC, sketchW, sketchH;
		
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
		var searchWindow;
		var modelWidth;
		
		var webglFi;
		
		var gidTime1, gidTime2;
		
		this.init = function(canvas) {
			// do all prep stuff
			
			sketchCC = canvas.getContext('2d');
			sketchW = canvas.width;
			sketchH = canvas.height;
			
			//console.log('starting initialization');
			
			// load from model
			numPatches = pModel.patchModel.numPatches;
			patchSize = pModel.patchModel.patchSize[0];
			searchWindow = params.searchWindow;
			numParameters = pModel.shapeModel.numEvalues;
			modelWidth = pModel.patchModel.canvasSize[0];
			
			// load eigenvectors
			eigenVectors = new goog.math.Matrix(numPatches*2, numParameters);
			for (var i = 0;i < numPatches*2;i++) {
				for (var j = 0;j < numParameters;j++) {
					eigenVectors.setValueAt(i, j, pModel.shapeModel.eigenVectors[i][j]);
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
			for (var w = 0;w < numPatches;w++) {
				weightMatricesOld[w] = new goog.math.Matrix(patchSize, patchSize);
				// insert
				for (var i = 0;i < patchSize;i++) {
					for (var j = 0;j < patchSize;j++) {
						//weight = pModel.patchModel.weights[w+(((i*patchSize)+j)*numPatches)];
						weight = pModel.patchModel.weights[w][(i*patchSize)+j];
						// cut off all values above 1 and below -1
						weight = weight > 1 ? 1 : weight;
						weight = weight < -1 ? -1 : weight;
						weightMatricesOld[w].setValueAt(i, j, weight);
					}
				}
			}
			
			//var weights = [];
			for (var w = 0;w < numPatches;w++) {
			  for (var i = 0;i < (patchSize*patchSize);i++) {
			    weight = pModel.patchModel.weights[w][i];
			    weight = weight > 1 ? 1 : weight;
          weight = weight < -1 ? -1 : weight;
          weights[i +(w*patchSize*patchSize)] = weight;
			  }
			}
			
			// precalculate gaussianPriorDiagonal
			gaussianPD = new goog.math.Matrix(numParameters+4, numParameters+4);
			// set values and append manual inverse
			for (var i = 0;i < numParameters;i++) {
			  gaussianPD.setValueAt(i+4,i+4,1/eigenValues[i]);
			}
				
			for (var i = 0;i < numParameters+4;i++) {
				currentParameters[i] = 0;
			}
			
			// set up webgl filter calculation
			webglFi = new webglFilter();
      webglFi.init(weights, numPatches, searchWindow+patchSize, searchWindow+patchSize, patchSize, patchSize, true);
			
			//console.log('ended initialization');
		}
		
		var gaussianProb = function(mean, coordinate, variance) {
		  // calculate pdf gaussian probability
		  var dx = coordinate[0] - mean[0];
			var dy = coordinate[1] - mean[1];
			var prob = Math.exp(-0.5*((dx*dx)+(dy*dy))/variance);
			
			return prob;
		}
		
		var createJacobian = function(parameters, eigenVectors) {
			// generates the jacobian matrix

			var jacobian = new goog.math.Matrix(2*numPatches, numParameters+4);
			var j0,j1;
			for (var i = 0;i < numPatches;i ++) {
				// 1
				j0 = meanShape[i][0];
				j1 = meanShape[i][1];
				for (var p = 0;p < numParameters;p++) {
					j0 += parameters[p+4]*eigenVectors.getValueAt((i*2),p);
					j1 += parameters[p+4]*eigenVectors.getValueAt((i*2)+1,p);
				}
				jacobian.setValueAt((i*2), 0, j0);
				jacobian.setValueAt((i*2)+1, 0, j1);
				// 2
				j0 = meanShape[i][1];
				j1 = meanShape[i][0];
				for (var p = 0;p < numParameters;p++) {
					j0 += parameters[p+4]*eigenVectors.getValueAt((i*2)+1,p);
					j1 += parameters[p+4]*eigenVectors.getValueAt((i*2),p);
				}
				jacobian.setValueAt(i*2, 1, -j0);
				jacobian.setValueAt((i*2)+1, 1, j1);
				// 3
				jacobian.setValueAt((i*2), 2, 1)
				jacobian.setValueAt((i*2)+1, 2, 0)
				// 4
				jacobian.setValueAt((i*2), 3, 0)
				jacobian.setValueAt((i*2)+1, 3, 1)
				// the rest
				for (var j = 0;j < numParameters;j++) {
					j0 = parameters[0]*eigenVectors.getValueAt(i*2,j) - parameters[1]*eigenVectors.getValueAt((i*2)+1,j) + eigenVectors.getValueAt(i*2,j);
					j1 = parameters[0]*eigenVectors.getValueAt((i*2)+1,j) + parameters[1]*eigenVectors.getValueAt((i*2),j) + eigenVectors.getValueAt((i*2)+1,j);
					jacobian.setValueAt(i*2,j+4,j0);
					jacobian.setValueAt((i*2)+1,j+4,j1);
				}
			}
			
			return jacobian;
		}
		
		var calculatePositions = function(parameters, useTransforms) {
			var x, y;
			var numParameters = parameters.length;
			positions = [];
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
      var comp = ccv.detect_objects(
        ccv.grayscale(canvas), ccv.cascade, 5, 1
      );
      
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
		
		/*
     *  element : canvas or video element
     *  TODO: should be able to take img element as well
     */
		this.track = function(element) {
			
			// for timing only
			var startTime = (new Date).getTime();
			
			var scaling, translateX, translateY;
			var croppedPatches = [];
			var ptch, px, py, pw;
			
			if (first) {
				// do viola-jones on canvas to get initial guess, if we don't have any points
				
				this.detectPosition(element);
				
				// set translation parameters
				
				// calculate modelWidth/height from meanshape
        var xmin = ymin = 1000000;
        var xmax = ymax = 0;
        for (var i = 0;i < meanShape.length;i++) {
          xmin = meanShape[i][0] < xmin ? meanShape[i][0] : xmin;
          xmax = meanShape[i][0] > xmax ? meanShape[i][0] : xmax;
          ymin = meanShape[i][1] < ymin ? meanShape[i][1] : ymin;
          ymax = meanShape[i][1] > ymax ? meanShape[i][1] : ymax;
        }
        var modelwidth = xmax-xmin;
        var modelheight = ymax-ymin;
        
        scaling = candidate.width/modelheight;
        var ccc = document.getElementById('overlay').getContext('2d');
        ccc.strokeRect(candidate.x,candidate.y,candidate.width,candidate.height);
        translateX = candidate.x-(xmin*scaling)+0.1*candidate.width;
			  translateY = candidate.y-(ymin*scaling)+0.25*candidate.height;
        currentParameters[0] = scaling-1;
        currentParameters[2] = translateX;
        currentParameters[3] = translateY;
        
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
			
			pw = pl = patchSize+searchWindow;
			
			var pdata, pmatrix, grayscaleColor, pdataLength;
			for (var i = 0; i < numPatches; i++) {
			  px = patchPositions[i][0]-(pw/2);
			  py = patchPositions[i][1]-(pl/2);
			  ptch = sketchCC.getImageData(px >> 0, py >> 0, pw, pl);
			  // add channels and convert to grayscale
			  // load patchdata to matrix
			  pdata = ptch.data;
			  pdataLength = pw*pl;
			  /*pmatrix = new goog.math.Matrix(pw, pl);
			  for (var j = 0;j < pdataLength;j++) {
			    grayscaleColor = pdata[j*4]*0.3 + pdata[(j*4)+1]*0.59 + pdata[(j*4)+2]*0.11;
			    pmatrix.setValueAt(j % pw, (j / pw) >> 0, grayscaleColor)
			  }*/
			  pmatrix = [];
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
			  var psci = sketchCC.createImageData(patchSize+searchWindow ,patchSize+searchWindow );
			  var pscidata = psci.data;
			  for (var j = 0;j < (patchSize+searchWindow )*(patchSize+searchWindow );j++) {
			    var val = patches[i][(((j / (patchSize+searchWindow)) >> 0)*(patchSize+searchWindow))+(j % (patchSize+searchWindow))];
			    //var val = patches[i].getValueAt(j % (patchSize+searchWindow ), (j / (patchSize+searchWindow )) >> 0);
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
			
      var responses = webglFi.getResponses(patches, true);
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
			var jac, vecProbs;
			var meanshiftVectors = [];
			
			var dootholder = 0;
			var partbholder = 0;
			var partaholder = 0;
			
			//var maxDiff = 0;
			
			for (var i = 0; i < varianceSeq.length; i++) {
				
				var partastart = (new Date).getTime();
				
				// calculate jacobian
				jac = createJacobian(currentParameters, eigenVectors);

				/* continue webgl work here? */
				
				//debugging
				var debugMVs = [];
				//
				
				var partaend = (new Date).getTime();
				partaholder += (partaend-partastart);
				
				for (var j = 0;j < numPatches;j++) {
				  // for every point in each response:
					// calculate PI x gaussian
				  vecProbs = [];
				  
          var boot = (new Date).getTime();
          
          var halfSearchWindow = searchWindow/2;
          var opj0 = originalPositions[j][0]-halfSearchWindow;
          var opj1 = originalPositions[j][1]-halfSearchWindow;
          
				  /*goog.math.Matrix.forEach(responses[j], function(value, a, b) {
				    var pos = (searchWindow*b)+a;
				    var updatePosition = [opj0+b, opj1+a];
				    vecProbs[pos] = value * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
				  });*/
				  var resplength = searchWindow*searchWindow;
				  var a, b;
				  for (var k = 0;k < resplength;k++) {
            
            a = k % searchWindow;
            b = (k / searchWindow) >> 0;
            //var pos = (searchWindow*b)+a;
				    var updatePosition = [opj0+b, opj1+a];
				    
				    /*if (Math.abs(updatePosition[0]-currentPositions[j][0]) > maxDiff) {
				      maxDiff = Math.abs(updatePosition[0]-currentPositions[j][0]);
				    }
				    if (Math.abs(updatePosition[1]-currentPositions[j][1]) > maxDiff) {
				      maxDiff = Math.abs(updatePosition[1]-currentPositions[j][1]);
				    }*/
				    
				    //vecProbs[pos] = responses[j][k] * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
				    vecProbs[k] = responses[j][k] * gaussianProb(updatePosition, currentPositions[j], varianceSeq[i]);
				  }
				  
				  var doot = (new Date).getTime();
          dootholder += (doot-boot);
				  
				  // sum PI x gaussians 
				  var vpsum = 0;
				  for (var k = 0;k < vecProbs.length;k++) {
				    vpsum += vecProbs[k];
				  }
				  
				  //debugging
				  //var debugMatrixMV = new goog.math.Matrix(searchWindow, searchWindow);
				  //
				  
				  var vecsum = 0;
				  var vecpos = [0,0];
				  for (var k = 0;k < searchWindow;k++) {
				    for (var l = 0;l < searchWindow;l++) {
				      var pos = (searchWindow*l)+k;
				      var updatePosition = [originalPositions[j][0]+l-(searchWindow/2), originalPositions[j][1]+k-(searchWindow/2)];
				      vecsum = (vecProbs[pos]/vpsum)
				      //debugging
				      //debugMatrixMV.setValueAt(k,l,vecsum);
				      //
				      vecpos[0] += vecsum*updatePosition[0];
				      vecpos[1] += vecsum*updatePosition[1];
				    }
				  }
				  
				  // compute mean shift vectors
				  var msv = [];
				  msv[0] = vecpos[0] - currentPositions[j][0];
				  msv[1] = vecpos[1] - currentPositions[j][1];
				  meanshiftVectors[j] = msv;
				  
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
            var val = debugMVs[k].getValueAt((j / searchWindow) >> 0, j % searchWindow);
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
				var meanShiftVector = new goog.math.Matrix(numPatches*2, 1);
				for (var k = 0;k < numPatches;k++) {
				  meanShiftVector.setValueAt(k*2, 0, meanshiftVectors[k][0]);
				  meanShiftVector.setValueAt((k*2)+1, 0, meanshiftVectors[k][1]);
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
				var prior = gaussianPD.multiply(varianceSeq[i]);
        var jtj = jac.getTranspose().multiply(jac);
        var cpMatrix = new goog.math.Matrix((numParameters+4),1);
        for (var l = 0;l < (numParameters+4);l++) {
          cpMatrix.setValueAt(l,0,currentParameters[l]);
        }
				var priorP = prior.multiply(cpMatrix);
				var jtv = jac.getTranspose().multiply(meanShiftVector);
				var paramUpdateLeft = (prior.add(jtj)).getInverse();
				var paramUpdateRight = priorP.subtract(jtv);
				var paramUpdate = paramUpdateLeft.multiply(paramUpdateRight);
				
				var oldPositions = currentPositions;
				
				// update estimated parameters
				for (var k = 0;k < numParameters+4;k++) {
				  currentParameters[k] -= paramUpdate.getValueAt(k,0);
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
			//console.log("maxdiff:"+maxDiff);
			
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
			var i, x, y;
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
		  var i, x, y;
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
			  if (typeof(paths.i) == 'number') {
			    drawPoint(cc, paths[i], params);
			  } else {
			    drawPath(cc, paths[i], params);
			  }
			}
			
			cc.restore()
		}
		
		return true
	}
}