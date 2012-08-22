//requires: ccv, closure matrix library

var clm = {
	tracker : function(params) {
    
    if (!params) params = {};
		if (params.constantVelocity === undefined) params.constantVelocity = false;
	
		var numPatches, patchSize, numParameters;
		var gaussianPD, gaussianTable;
		var eigenVectors, eigenValues;
		var sketchCC, sketchW, sketchH;
		
		var currentParameters = [];
		var currentPositions = [];
		var previousParameters = [];
		
		var weightMatrices = [];
		var weightMatricesOld = [];
		
		var patches = [];
		var responses = [];
		var meanShape = [];
		
		//var varianceSeq = [20,10,5,1,1,1,1,1,1,1];
		//var varianceSeq = [20,20,20,20,20,20,20,20,20,20];
		//var varianceSeq = [10,10,10,10,10,10,10,10,10,10];
		//var varianceSeq = [5,5,1,1,1,1,1,1,1,1];
		//var varianceSeq = [20,10,5,1];
		var varianceSeq = [20,10,5,1];
		//var varianceSeq = [50,25,12.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5];
		var maxSearches = 4;
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
			
			console.log('starting initialization');
			
			// load from model
			numPatches = pModel.patchModel.numPatches;
			patchSize = pModel.patchModel.patchSize[0];
			searchWindow = Math.ceil(patchSize*1.25);
			numParameters = pModel.shapeModel.numEvalues;
			modelWidth = 0.75*pModel.patchModel.canvasSize[0];
			
			// load eigenvectors
			eigenVectors = new goog.math.Matrix(numPatches*2, numParameters);
			for (var i = 0;i < numPatches*2;i++) {
				for (var j = 0;j < numParameters;j++) {
					eigenVectors.setValueAt(i, j, pModel.shapeModel.eVectors[(j*numPatches*2)+i]);
				}
			}
			
			// load mean shape
			for (var i = 0; i < numPatches;i++) {
			  meanShape[i] = [pModel.shapeModel.meanShape[(i*2)], pModel.shapeModel.meanShape[(i*2)+1]];
			}
			
			// load eigenvalues
			eigenValues = pModel.shapeModel.eValues;
			
			// load patchweight matrices for comparing
			var weight;
			for (var w = 0;w < numPatches;w++) {
				weightMatricesOld[w] = new goog.math.Matrix(patchSize, patchSize);
				// insert
				for (var i = 0;i < patchSize;i++) {
					for (var j = 0;j < patchSize;j++) {
						weight = pModel.patchModel.weights[w+(((i*patchSize)+j)*numPatches)];
						// cut off all values above 1 and below -1
						weight = weight > 1 ? 1 : weight;
						weight = weight < -1 ? -1 : weight;
						weightMatricesOld[w].setValueAt(i, j, weight);
					}
				}
			}
			
			var weights = [];
			for (var w = 0;w < numPatches;w++) {
			  for (var i = 0;i < (patchSize*patchSize);i++) {
			    weight = pModel.patchModel.weights[w+(i*numPatches)];
			    weight = weight > 1 ? 1 : weight;
          weight = weight < -1 ? -1 : weight;
          weights[Math.floor(i / patchSize)*patchSize + (i % patchSize) +(w*patchSize*patchSize)] = weight;
			  }
			}
			
			// precalculate gaussianPriorDiagonal
			gaussianPD = new goog.math.Matrix(numParameters+4, numParameters+4);
			// manual inverse
			gaussianPD.setValueAt(4,4,1/eigenValues[0]);
			gaussianPD.setValueAt(5,5,1/eigenValues[1]);
			gaussianPD.setValueAt(6,6,1/eigenValues[2]);
			gaussianPD.setValueAt(7,7,1/eigenValues[3]);
			gaussianPD.setValueAt(8,8,1/eigenValues[4]);
			gaussianPD.setValueAt(9,9,1/eigenValues[5]);
			gaussianPD.setValueAt(10,10,1/eigenValues[6]);
			//gaussianPD = gaussianPD.getInverse();
			
			// load precalculated gaussian bivariate table
			gaussianTable = [];
			var gx, gy;
			var gtdim = Math.sqrt(normProbs[0].length);
			for (var i = 0;i < normProbs.length;i++) {
				/*gaussianTable[i] = new goog.math.Matrix(161,161);
				for (var j = 0; j < 25921;j++) {
					gx = j % 161;
					gy = j / 161 >> 0;
					gaussianTable[i].setValueAt(gx, gy, normProbs[i][j]);
				}*/
				gaussianTable[i] = [];
				for (var j = 0;j < gtdim;j++) {
				  gaussianTable[i][j] = [];
				}
				for (var j = 0;j< normProbs[0].length;j++) {
				  gx = j % gtdim;
					gy = j / gtdim >> 0;
					gaussianTable[i][gx][gy] = normProbs[i][j];
				}
			}
				
			for (var i = 0;i < numParameters+4;i++) {
				currentParameters[i] = 0;
			}
			
			// set up webgl filter calculation
			webglFi = new webglFilter();
      webglFi.init(weights, numPatches, searchWindow+patchSize, searchWindow+patchSize, patchSize, patchSize, true);
			
			console.log('ended initialization');
		}
		
		var gaussianProb = function(mean, coordinate, table) {
			// looks up the relevant gaussian probability in table
			var correctedCoord = [Math.abs(coordinate[0] - mean[0]), Math.abs(coordinate[1] - mean[1])];
			var cx = (correctedCoord[0]+0.1) / 0.1 >> 0;
			var cy = (correctedCoord[1]+0.1) / 0.1 >> 0;
			//var prob = gaussianTable[varianceSeq.indexOf(variance)].getValueAt(cx,cy);
			var prob = table[cx][cy];
			
			return prob;
		}
		
		var createJacobian = function(parameters, eigenVectors) {
			// generates the jacobian matrix

			var jacobian = new goog.math.Matrix(2*numPatches, numParameters+4);
			var j0,j1;
			for (var i = 0;i < numPatches;i ++) {
				// 1
				//j0 = eigenVectors.getValueAt(i,0);
				//j1 = eigenVectors.getValueAt(i+1,0);
				j0 = meanShape[i][0];
				j1 = meanShape[i][1];
				for (var p = 0;p < numParameters;p++) {
					j0 += parameters[p+4]*eigenVectors.getValueAt((i*2),p);
					j1 += parameters[p+4]*eigenVectors.getValueAt((i*2)+1,p);
				}
				jacobian.setValueAt((i*2), 0, j0);
				jacobian.setValueAt((i*2)+1, 0, j1);
				// 2
				//j0 = eigenVectors.getValueAt(i+1,0);
				//j1 = eigenVectors.getValueAt(i,0);
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
				x = pModel.shapeModel.meanShape[(2*i)];
				y = pModel.shapeModel.meanShape[(2*i)+1];
				for (var j = 0;j < numParameters-4;j++) {
					x += pModel.shapeModel.eVectors[(j*68*2)+(2*i)]*parameters[j+4];
					y += pModel.shapeModel.eVectors[(j*68*2)+((2*i)+1)]*parameters[j+4];
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
		  if (el.tagName == 'VIDEO') {
		    var canvas = document.createElement('canvas');
		    canvas.width = el.width;
		    canvas.height = el.height;
		    var cc = canvas.getContext('2d');
		    cc.drawImage(el, 0, 0, el.width, el.height);
		  } else if (el.tagName == 'CANVAS') {
		    var canvas = el;
		    var cc = canvas.getContext('2d');
		  }
		  
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
     *
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
				
				// TODO: can't use this here, need way to pass back scaling etc.
				this.detectPosition(element);
				
				// set translation parameters
        scaling = candidate.width/modelWidth;
        translateX = candidate.x-Math.round(candidate.width*0.08);
			  translateY = candidate.y+Math.round(candidate.height*0.08);
        currentParameters[0] = scaling-1;
        currentParameters[2] = translateX;
        currentParameters[3] = translateY;
        
        currentPositions = calculatePositions(currentParameters, true);
				
				first = false;
			} else {
				// TODO : do cross-correlation or similar to find translation of face (feature detection?)
				
				if (params.constantVelocity) {
          // calculate where to get patches via constant velocity prediction
          if (previousParameters.length >= 3) {
            for (var i = 0;i < currentParameters.length;i++) {
              currentParameters[i] = (2*previousParameters[1][i]) - previousParameters[0][i];
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
			console.log("detectiontiming: "+(secondTime-startTime)+" ms")
			
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
			
			console.log("gidtime:"+(gidTime2-gidTime1));
			
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
			  px = patchPositions[i][0]-(32/2);
			  py = patchPositions[i][1]-(32/2);
			  var psci = sketchCC.createImageData(32 ,32 );
			  var pscidata = psci.data;
			  for (var j = 0;j < (32)*(32);j++) {
			    var val = weightMatricesOld[i].getValueAt(j % (32), (j / (32)) >> 0);
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
			    //var val = patches[i].getValueAt(j % (32+searchWindow ), (j / (32+searchWindow )) >> 0);
			    pscidata[j*4] = val;
			    pscidata[(j*4)+1] = val;
			    pscidata[(j*4)+2] = val;
			    pscidata[(j*4)+3] = 255;
			  }
			  sketchCC.putImageData(psci, px >> 0, py >> 0);
			}*/
			console.log("starting response calculations");
			//sketchCC.clearRect(0, 0, sketchW, sketchH);
			
			var thirdTime = (new Date).getTime();
			console.log("getting patches: "+(thirdTime-secondTime)+" ms")
			
			/* start webgl work here */
			/*var patchResponse, submatrix, submsum;
			for (var i = 0;i < numPatches;i++) {
			  // calculate response map
			  console.log("response: "+i);
			  patchResponse = new goog.math.Matrix(searchWindow, searchWindow);
			  for (var k = 0;k < searchWindow;k++) {
			    for (var l = 0;l < searchWindow;l++) {
			      // get submatrix for this patch
			      submatrix = new goog.math.Matrix(32,32);
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
			console.log("calculating responses: "+(secondTime - thirdTime)+" ms")
			
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
			
			var gaussianLookup;
			//var maxDiff = 0;
			
			for (var i = 0; i < maxSearches; i++) {
				
				var partastart = (new Date).getTime();
				
				// calculate jacobian
				jac = createJacobian(currentParameters, eigenVectors);

				/* continue webgl work here? */
				
				//debugging
				var debugMVs = [];
				//
				
				var partaend = (new Date).getTime();
				partaholder += (partaend-partastart);
				
				gaussianLookUp = gaussianTable[gaussianTableOrder.indexOf(varianceSeq[i])];
				
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
				    vecProbs[pos] = value * gaussianProb(updatePosition, varianceSeq[i], currentPositions[j], gaussianLookUp);
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
				    
				    //vecProbs[pos] = responses[j][k] * gaussianLookUp[Math.abs(currentPositions[j][0] - (opj0+b)) / 0.1 >> 0][Math.abs(currentPositions[j][1] - (opj1+a)) / 0.1 >> 0];
				    vecProbs[k] = responses[j][k] * gaussianProb(updatePosition, currentPositions[j], gaussianLookUp);
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
				var meanShiftVector = new goog.math.Matrix(68*2, 1);
				for (var k = 0;k < 68;k++) {
				  meanShiftVector.setValueAt(k*2, 0, meanshiftVectors[k][0]);
				  meanShiftVector.setValueAt((k*2)+1, 0, meanshiftVectors[k][1]);
				}
				
				
				//test : draw meanshiftvector diagram
				/*var testMV = [];
				for (var k = 0;k < 68;k++) {
				  testMV[k] = []
				  testMV[k][0] = currentPositions[k][0] + meanShiftVector.getValueAt(k*2, 0);
				  testMV[k][1] = currentPositions[k][1] + meanShiftVector.getValueAt((k*2)+1, 0);
				}*/
				//this.drawPos(canvas, testMV, 1);
				//
				
				// compute pdm parameter update
				var prior = gaussianPD.multiply(varianceSeq[i]);
        var jtj = jac.getTranspose().multiply(jac);
        var cpMatrix = new goog.math.Matrix(11,1);
        for (var l = 0;l < 11;l++) {
          cpMatrix.setValueAt(l,0,currentParameters[l]);
        }
				var priorP = prior.multiply(cpMatrix);
				var jtv = jac.getTranspose().multiply(meanShiftVector);
				var paramUpdateLeft = (prior.add(jtj)).getInverse();
				var paramUpdateRight = priorP.subtract(jtv);
				var paramUpdate = paramUpdateLeft.multiply(paramUpdateRight);
				
				// check if converged
        // calculate norm of parameterdifference
				var parameterNorm = 0;
				var pnsq;
				for (var k = 0;k < 4+7;k++) {
				  pnsq = paramUpdate.getValueAt(k,0);
				  parameterNorm += (pnsq * pnsq);
				}
				parameterNorm = Math.sqrt(parameterNorm);
				console.log("parameternorm:"+parameterNorm);
				
				// update estimated parameters
				for (var k = 0;k < 7+4;k++) {
				  currentParameters[k] -= paramUpdate.getValueAt(k,0);
				}
				
				// if norm < limit, then break
				if (parameterNorm < convergenceLimit) {
				  break;
				}
				// update current coordinates
				currentPositions = calculatePositions(currentParameters, true);
				
				var partbend = (new Date).getTime();
				partbholder += (partbend-partbstart);
				
			}
			//console.log("maxdiff:"+maxDiff);
			
			if (params.constantVelocity) {
        // add current parameter to array of previous parameters
        previousParameters.push(currentParameters.slice());
        previousParameters.splice(0, previousParameters.length == 3 ? 1 : 0);
			}
			
			console.log("doot: "+(dootholder)+" ms")
			console.log("partaholder: "+(partaholder)+" ms")
			console.log("partbholder: "+(partbholder)+" ms")
			
      document.getElementById('overlay').getContext('2d').clearRect(0, 0, 720, 576);
			//debugger;
			
			this.draw(document.getElementById('overlay'), currentParameters);
			
			var thirdTime = (new Date).getTime();
			console.log("optimization of parameters: "+(thirdTime - secondTime)+" ms")
			console.log("all time: "+(thirdTime - startTime)+" ms")
			
			// return new points
			return currentPositions;
		}
		
		this.constrain = function(updatePosition) {
		  var jac = createJacobian(currentParameters, eigenVectors);
		  
		  // compute pdm parameter update
      var prior = gaussianPD.multiply(varianceSeq[0]);
      var jtj = jac.getTranspose().multiply(jac);
      var cpMatrix = new goog.math.Matrix(11,1);
      for (var l = 0;l < 11;l++) {
        cpMatrix.setValueAt(l,0,currentParameters[l]);
      }
			var priorP = prior.multiply(cpMatrix);
			
			var up = new goog.math.Matrix(68*2,1);
			for (var l = 0;l < 68;l++) {
			  up.setValueAt(l*2,0,updatePosition[l][0]);
			  up.setValueAt((l*2)+1,0,updatePosition[l][1]);
			}
			
			var jtv = jac.getTranspose().multiply(up);
			var paramUpdateLeft = (prior.add(jtj)).getInverse();
			var paramUpdateRight = priorP.subtract(jtv);
			var paramUpdate = paramUpdateLeft.multiply(paramUpdateRight);
		  
		  return paramUpdate
		}
		
		var normalizeMatrix = function(matrix) {
		  
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
		
		this.draw = function(canvas, pv, scale) {
			// if no previous points, just draw in the middle of canvas
			
			var params;
			if (pv === undefined) {
			  params = currentParameters.slice(0);
			} else {
        params = pv.slice(0);
			}
			if (scale === undefined) {
			  scale = 1;
			}
			
			/*
			if (params[1] == 0) {
			  rotation = 0;
			  var paramscale = (params[0]+1)/Math.cos(0);
			} else {
			  var rotation = 0.5*Math.asin(2*(params[0]+1)*params[1]);
        var paramscale = params[1]/Math.sin(rotation);
			}
			
			params[0] = (paramscale*scale*Math.cos(rotation))-1;
			params[1] = paramscale*scale*Math.sin(rotation);
			*/
			
			var cc = canvas.getContext('2d');
			cc.fillStyle = "rgb(200,200,200)";
			var x, y, i, path;
			
			cc.save();
			
			path = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			//right eyebrow
			path = [15,16,17,18,19,20,15]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// left eyebrow
			path = [21,22,23,24,25,26,21]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// left eye
			path = [27,28,29,30,27,30,31]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// right eye
			path = [32,33,34,35,32,35,36]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// nose
			path = [37,38,39,40,46,41,47,42,43,44,45]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// mouth
			path = [48,49,50,51,52,53,54,55,56,57,58,59,48,60,61,62,54,63,64,65,48]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i];
				y = pModel.shapeModel.meanShape[i+1];
				for (var j = 0;j < 7;j++) {
					x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
					y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
				}
				a = params[0]*x - params[1]*y + params[2];
				b = params[0]*y + params[1]*x + params[3];
				x += a;
				y += b;
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// mid mouth
			i = 66*2;
			x = pModel.shapeModel.meanShape[i];
			y = pModel.shapeModel.meanShape[i+1];
			for (var j = 0;j < 7;j++) {
				x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
				y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
			}
			a = params[0]*x - params[1]*y + params[2];
			b = params[0]*y + params[1]*x + params[3];
			x += a;
			y += b;
			cc.fillRect(x,y,3,3);
			
			// mid nose
			i = 67*2;
			x = pModel.shapeModel.meanShape[i];
			y = pModel.shapeModel.meanShape[i+1];
			for (var j = 0;j < 7;j++) {
				x += pModel.shapeModel.eVectors[((j)*68*2)+i]*params[j+4];
				y += pModel.shapeModel.eVectors[((j)*68*2)+(i+1)]*params[j+4];
			}
			a = params[0]*x - params[1]*y + params[2];
			b = params[0]*y + params[1]*x + params[3];
			x += a;
			y += b;
			cc.fillRect(x,y,3,3);
			
			cc.restore()
		}
	
    this.drawPos = function(canvas, positions) {
			// if no previous points, just draw in the middle of canvas
			
			var cc = canvas.getContext('2d');
			cc.fillStyle = "rgb(200,200,200)";
			var x, y, i, path;
			
			cc.save();
			
			path = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			//right eyebrow
			path = [15,16,17,18,19,20,15]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// left eyebrow
			path = [21,22,23,24,25,26,21]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// left eye
			path = [27,28,29,30,27,30,31]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// right eye
			path = [32,33,34,35,32,35,36]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// nose
			path = [37,38,39,40,46,41,47,42,43,44,45]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// mouth
			path = [48,49,50,51,52,53,54,55,56,57,58,59,48,60,61,62,54,63,64,65,48]
			cc.beginPath();
			for (var p = 0;p < path.length;p++) {
				i = path[p];
				x = positions[i][0];
				y = positions[i][1];
				//cc.fillRect(x,y,3,3);
				if (i == 0) {
					cc.moveTo(x,y);
				} else {
					cc.lineTo(x,y);
				}
			}
			cc.moveTo(0,0);
			cc.closePath();
			cc.stroke();
			
			// mid mouth
			i = 66;
			x = positions[i][0];
			y = positions[i][1];
			cc.fillRect(x,y,3,3);
			
			// mid nose
			i = 67;
			x = positions[i][0];
			y = positions[i][1];
			cc.fillRect(x,y,3,3);
			
			cc.restore()
		}
		
		return true
	}
}