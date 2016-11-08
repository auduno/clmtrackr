import EventEmitter from 'events';

// 3rd party libs
import numeric from 'numeric';

// libs
import MosseFilter from './utils/mosse';
import MosseFilterResponses from './utils/mosseFilterResponses';
import JsfeatFace from './jsfeat/JsfeatFace';
import WebglFilter from './svmfilter/svmfilter_webgl';
import SvmFilter from './svmfilter/svmfilter_fft';

// filters
// import entireFaceFilter from './filters/entire_face_filter.json';
import leftEyeFilter from './filters/left_eye_filter.json';
import rightEyeFilter from './filters/right_eye_filter.json';
import noseFilter from './filters/nose_filter.json';

import modelPca20Svm from '!json!clmtrackr/models/model_pca_20_svm.json';

import {
  requestAnimFrame,
  cancelRequestAnimFrame
} from './utils/anim';
import { gpopt, gpopt2 } from './utils/meanshift';
import procrustes from './utils/procrustes';


const halfPI = Math.PI / 2;

const DEFAULT_MODEL = modelPca20Svm;

const HAS_MOSSE_FILTERS = MosseFilter && leftEyeFilter && rightEyeFilter && noseFilter;
const VALID_RESPONSEMODE_LIST = ['raw', 'sobel', 'lbp'];
const VALID_RESPONSEMODES = ['single', 'blend', 'cycle'];


export default class Tracker extends EventEmitter {
  constructor (params) {
    super();

    if (!params) params = {};
    if (params.constantVelocity === undefined) params.constantVelocity = true;
    if (params.searchWindow === undefined) params.searchWindow = 11;
    if (params.useWebGL === undefined) params.useWebGL = true;
    if (params.scoreThreshold === undefined) params.scoreThreshold = 0.5;
    if (params.stopOnConvergence === undefined) params.stopOnConvergence = false;
    if (params.weightPoints === undefined) params.weightPoints = undefined;
    if (params.sharpenResponse === undefined) params.sharpenResponse = false;

    this.params = params;

    /** @type {Number} Minimum convergence before firing `converged` event. */
    this.convergenceThreshold = 0.5;

    this.numPatches = undefined;
    this.patchSize = undefined;
    this.numParameters = undefined;
    this.patchType = undefined;

    this.gaussianPD = undefined;
    this.eigenVectors = undefined;
    this.eigenValues = undefined;

    this.sketchCC = undefined;
    this.sketchW = undefined;
    this.sketchH = undefined;
    this.sketchCanvas = undefined;

    this.candidate = undefined;
    this.weights = undefined;
    this.model = undefined;
    this.biases = undefined;

    this.sobelInit = false;
    this.lbpInit = false;

    this.currentParameters = [];
    this.currentPositions = [];
    this.previousParameters = [];
    this.previousPositions = [];

    this.patches = [];
    this.responses = [];
    this.meanShape = [];

    this.responseMode = 'single';
    this.responseList = ['raw'];
    this.responseIndex = 0;

    /*
    It's possible to experiment with the sequence of variances used for the finding the maximum in the KDE.
    This sequence is pretty arbitrary, but was found to be okay using some manual testing.
    */
    this.varianceSeq = [10, 5, 1];
    // this.varianceSeq = [3, 1.5, 0.75];
    // this.varianceSeq = [6, 3, 0.75];
    this.PDMVariance = 0.7;

    this.relaxation = 0.1;

    this.first = true;
    this.gettingPosition = false;

    this.convergenceLimit = 0.01;

    this.learningRate = [];
    this.stepParameter = 1.25;
    this.prevCostFunc = []

    this.searchWindow = undefined;

    this.modelWidth = undefined;
    this.modelHeight = undefined;

    this.halfSearchWindow = undefined;
    this.vecProbs = undefined;
    this.responsePixels = undefined;

    if (typeof Float64Array !== 'undefined') {
      this.updatePosition = new Float64Array(2);
      this.vecpos = new Float64Array(2);
    } else {
      this.updatePosition = new Array(2);
      this.vecpos = new Array(2);
    }
    this.pw = undefined;
    this.pl = undefined;
    this.pdataLength = undefined;

    this.facecheck_count = 0;

    this.webglFi = undefined;
    this.svmFi = undefined;
    this.mosseCalc = undefined;

    this.scoringCanvas = document.createElement('canvas');
    // document.body.appendChild(this.scoringCanvas);
    this.scoringContext = this.scoringCanvas.getContext('2d');
    this.msxmin = undefined;
    this.msymin = undefined;
    this.msxmax = undefined;
    this.msymax = undefined;
    this.msmodelwidth = undefined;
    this.msmodelheight = undefined;
    this.scoringWeights = undefined;
    this.scoringBias = undefined;

    this.scoringHistory = [];
    this.meanscore = 0;

    this.mossef_lefteye = undefined;
    this.mossef_righteye = undefined;
    this.mossef_nose = undefined;

    this.right_eye_position = [0.0, 0.0];
    this.left_eye_position = [0.0, 0.0];
    this.nose_position = [0.0, 0.0];
    this.lep = undefined;
    this.rep = undefined;
    this.mep = undefined;
    this.runnerTimeout = undefined;
    this.runnerElement = undefined;
    this.runnerBox = undefined;

    this.pointWeights = undefined;

    this.jsfeatFace = new JsfeatFace();

    this._running = false;

    this._meanshiftVectors = [];
  }

  /*
   *  load model data, initialize filters, etc.
   *
   *  @param  <Object>  pdm model object
   */
  init (pdmmodel = DEFAULT_MODEL) {
    this.model = pdmmodel;

    // load from model
    this.patchType = this.model.patchModel.patchType;
    this.numPatches = this.model.patchModel.numPatches;
    this.patchSize = this.model.patchModel.patchSize[0];
    if (this.patchType === 'MOSSE') {
      this.searchWindow = this.patchSize;
    } else {
      this.searchWindow = this.params.searchWindow;
    }
    this.numParameters = this.model.shapeModel.numEvalues;
    this.modelWidth = this.model.patchModel.canvasSize[0];
    this.modelHeight = this.model.patchModel.canvasSize[1];

    // set up canvas to work on
    this.sketchCanvas = document.createElement('canvas');
    this.sketchCC = this.sketchCanvas.getContext('2d');

    this.sketchW = this.sketchCanvas.width = this.modelWidth + (this.searchWindow - 1) + this.patchSize - 1;
    this.sketchH = this.sketchCanvas.height = this.modelHeight + (this.searchWindow - 1) + this.patchSize - 1;

    if (this.model.hints && HAS_MOSSE_FILTERS) {
      // var mossef_lefteye = new MosseFilter({drawResponse : document.getElementById('overlay2')});
      this.mossef_lefteye = new MosseFilter();
      this.mossef_lefteye.load(leftEyeFilter);
      // var mossef_righteye = new MosseFilter({drawResponse : document.getElementById('overlay2')});
      this.mossef_righteye = new MosseFilter();
      this.mossef_righteye.load(rightEyeFilter);
      // var mossef_nose = new MosseFilter({drawResponse : document.getElementById('overlay2')});
      this.mossef_nose = new MosseFilter();
      this.mossef_nose.load(noseFilter);
    } else {
      console.log('MOSSE filters not found, using rough approximation for initialization.');
    }

    // load eigenvectors
    this.eigenVectors = numeric.rep([this.numPatches * 2, this.numParameters], 0.0);
    for (let i = 0; i < this.numPatches * 2; i++) {
      for (let j = 0; j < this.numParameters; j++) {
        this.eigenVectors[i][j] = this.model.shapeModel.eigenVectors[i][j];
      }
    }

    // load mean shape
    for (let i = 0; i < this.numPatches; i++) {
      this.meanShape[i] = [this.model.shapeModel.meanShape[i][0], this.model.shapeModel.meanShape[i][1]];
    }

    // get max and mins, width and height of meanshape
    this.msxmax = this.msymax = 0;
    this.msxmin = this.msymin = 1000000;
    for (let i = 0; i < this.numPatches; i++) {
      if (this.meanShape[i][0] < this.msxmin) this.msxmin = this.meanShape[i][0];
      if (this.meanShape[i][1] < this.msymin) this.msymin = this.meanShape[i][1];
      if (this.meanShape[i][0] > this.msxmax) this.msxmax = this.meanShape[i][0];
      if (this.meanShape[i][1] > this.msymax) this.msymax = this.meanShape[i][1];
    }
    this.msmodelwidth = this.msxmax - this.msxmin;
    this.msmodelheight = this.msymax - this.msymin;

    // get scoringweights if they exist
    if (this.model.scoring) {
      this.scoringWeights = new Float64Array(this.model.scoring.coef);
      this.scoringBias = this.model.scoring.bias;
      this.scoringCanvas.width = this.model.scoring.size[0];
      this.scoringCanvas.height = this.model.scoring.size[1];
    }

    // load eigenvalues
    this.eigenValues = this.model.shapeModel.eigenValues;

    this.weights = this.model.patchModel.weights;
    this.biases = this.model.patchModel.bias;

    // precalculate gaussianPriorDiagonal
    this.gaussianPD = numeric.rep(
      [this.numParameters + 4, this.numParameters + 4],
      0
    );
    // set values and append manual inverse
    for (let i = 0; i < this.numParameters; i++) {
      if (this.model.shapeModel.nonRegularizedVectors.indexOf(i) >= 0) {
        this.gaussianPD[i + 4][i + 4] = 1 / 10000000;
      } else {
        this.gaussianPD[i + 4][i + 4] = 1 / this.eigenValues[i];
      }
    }

    for (let i = 0; i < this.numParameters + 4; i++) {
      this.currentParameters[i] = 0;
    }

    if (this.patchType === 'SVM') {
      let webGLContext;
      const webGLTestCanvas = document.createElement('canvas');
      if (window.WebGLRenderingContext) {
        webGLContext = (
          webGLTestCanvas.getContext('webgl') ||
          webGLTestCanvas.getContext('experimental-webgl')
        );
        if (!webGLContext || !webGLContext.getExtension('OES_texture_float')) {
          webGLContext = null;
        }
      }

      if (
        webGLContext &&
        this.params.useWebGL &&
        typeof WebglFilter !== 'undefined'
      ) {
        this.webglFi = new WebglFilter();
        try {
          this.webglFi.init(
            this.weights,
            this.biases,
            this.numPatches,
            this.searchWindow + this.patchSize - 1,
            this.searchWindow + this.patchSize - 1,
            this.patchSize,
            this.patchSize
          );
          if ('lbp' in this.weights) this.lbpInit = true;
          if ('sobel' in this.weights) this.sobelInit = true;
        } catch (err) {
          console.error(err);
          alert('There was a problem setting up webGL programs, falling back to slightly slower javascript version. :(');
          this.webglFi = undefined;
          this.svmFi = new SvmFilter();
          this.svmFi.init(
            this.weights['raw'],
            this.biases['raw'],
            this.numPatches,
            this.patchSize,
            this.searchWindow
          );
        }
      } else if (typeof SvmFilter !== 'undefined') {
        // use fft convolution if no webGL is available
        this.svmFi = new SvmFilter();
        this.svmFi.init(
          this.weights['raw'],
          this.biases['raw'],
          this.numPatches,
          this.patchSize,
          this.searchWindow
        );
      } else {
        throw new Error('Could not initiate filters, please make sure that svmfilter.js or svmfilter_conv_js.js is loaded.');
      }
    } else if (this.patchType === 'MOSSE') {
      this.mosseCalc = new MosseFilterResponses();
      this.mosseCalc.init(this.weights, this.numPatches, this.patchSize, this.patchSize);
    }

    if (this.patchType === 'SVM') {
      this.pw = this.pl = this.patchSize + this.searchWindow - 1;
    } else {
      this.pw = this.pl = this.searchWindow;
    }
    this.pdataLength = this.pw * this.pl;
    this.halfSearchWindow = (this.searchWindow - 1) / 2;
    this.responsePixels = this.searchWindow * this.searchWindow;
    if (typeof Float64Array !== 'undefined') {
      this.vecProbs = new Float64Array(this.responsePixels);
      for (let i = 0; i < this.numPatches; i++) {
        this.patches[i] = new Float64Array(this.pdataLength);
      }
    } else {
      this.vecProbs = new Array(this.responsePixels);
      for (let i = 0; i < this.numPatches; i++) {
        this.patches[i] = new Array(this.pdataLength);
      }
    }

    for (let i = 0; i < this.numPatches; i++) {
      this.learningRate[i] = 1.0;
      this.prevCostFunc[i] = 0.0;
    }

    if (this.params.weightPoints) {
      // weighting of points
      this.pointWeights = [];
      for (let i = 0; i < this.numPatches; i++) {
        if (i in this.params.weightPoints) {
          this.pointWeights[(i * 2)] = this.params.weightPoints[i];
          this.pointWeights[(i * 2) + 1] = this.params.weightPoints[i];
        } else {
          this.pointWeights[(i * 2)] = 1;
          this.pointWeights[(i * 2) + 1] = 1;
        }
      }
      this.pointWeights = numeric.diag(this.pointWeights);
    }

    // Populate meanshift vectors
    for (let i = 0; i < this.numPatches; i++) {
      this._meanshiftVectors[i] = [0, 0];
    }
  }

  /*
   *  starts the tracker to run on a regular interval
   */
  start (element, box) {
    // check if model is initalized, else return false
    if (typeof this.model === 'undefined') {
      console.log('tracker needs to be initalized before starting to track.');
      return false;
    }
    // check if a runnerelement already exists, if not, use passed parameters
    if (typeof this.runnerElement === 'undefined') {
      this.runnerElement = element;
      this.runnerBox = box;
    }
    // start named timeout function
    this.runnerTimeout = requestAnimFrame(this._runnerFunction.bind(this));
    this._running = true;
    this.emit('started');
  }

  getIsRunning () {
    return this._running;
  }

  _runnerFunction () {
    this.runnerTimeout = requestAnimFrame(this._runnerFunction.bind(this));
    // // schedule as many iterations as we can during each request
    // var startTime = (new Date()).getTime();
    // while (((new Date()).getTime() - startTime) < 16) {
    //  var tracking = this.track(runnerElement, runnerBox);
    //  if (!tracking) continue;
    // }
    this.track(this.runnerElement, this.runnerBox);
  }

  /*
   *  stop the running tracker
   */
  stop () {
    // stop the running tracker if any exists
    cancelRequestAnimFrame(this.runnerTimeout);
    this._running = false;
    this.emit('stopped');
  }

  recheck () {
    console.log('RECHECKING');
    this.first = true;
  }

  /*
   *  element : canvas or video element
   *  TODO: should be able to take img element as well
   */
  track (element, box, gi) {
    this.emit('beforeTrack');

    var scaling, translateX, translateY, rotation;

    if (gi) {
      scaling = gi[0];
      rotation = gi[1];
      translateX = gi[2];
      translateY = gi[3];

      this.first = false;
      this.gettingPosition = false;
    } else if (this.first) {
      // do viola-jones on canvas to get initial guess, if we don't have any points
      if (!this.gettingPosition) {
        this.gettingPosition = true;
        this._getInitialPosition(element, box, (gi) => {
          this.gettingPosition = false;
          if (!gi) {
            // send an event on no face found
            this.emit('notFound');
            this.first = true;

            return false;
          } else {
            this.track(element, box, gi);
          }
        });
      }
      return;
    } else {
      this.facecheck_count += 1;

      if (this.params.constantVelocity) {
        // calculate where to get patches via constant velocity prediction
        if (this.previousParameters.length >= 2) {
          for (var i = 0; i < this.currentParameters.length; i++) {
            this.currentParameters[i] = (
              this.relaxation *
              this.previousParameters[1][i] +
              (1 - this.relaxation) *
              (
                (2 * this.previousParameters[1][i]) -
                this.previousParameters[0][i]
              )
            );
            // this.currentParameters[i] = (3*this.previousParameters[2][i]) - (3*this.previousParameters[1][i]) + this.previousParameters[0][i];
          }
        }
      }

      // change translation, rotation and scale parameters
      rotation = halfPI - Math.atan(
        (this.currentParameters[0] + 1) / this.currentParameters[1]
      );
      if (rotation > halfPI) {
        rotation -= Math.PI;
      }
      scaling = this.currentParameters[1] / Math.sin(rotation);
      translateX = this.currentParameters[2];
      translateY = this.currentParameters[3];
    }

    // copy canvas to a new dirty canvas
    this.sketchCC.save();

    // clear canvas
    this.sketchCC.clearRect(0, 0, this.sketchW, this.sketchH);

    this.sketchCC.scale(1 / scaling, 1 / scaling);
    this.sketchCC.rotate(-rotation);
    this.sketchCC.translate(-translateX, -translateY);

    this.sketchCC.drawImage(element, 0, 0, element.width, element.height);

    this.sketchCC.restore();
    //  get cropped images around new points based on model parameters (not scaled and translated)
    var patchPositions = this._calculatePositions(this.currentParameters, false);

    // check whether tracking is ok
    if (this.scoringWeights && (this.facecheck_count % 10 === 0)) {
      if (!this._checkTracking()) {
        // reset all parameters
        this._resetParameters();
        // send event to signal that tracking was lost
        this.emit('lost');
        return false;
      }
    }

    for (let i = 0; i < this.numPatches; i++) {
      const px = patchPositions[i][0] - (this.pw / 2);
      const py = patchPositions[i][1] - (this.pl / 2);
      const ptch = this.sketchCC.getImageData(
        Math.round(px),
        Math.round(py),
        this.pw,
        this.pl
      );
      const pdata = ptch.data;

      // convert to grayscale
      const pmatrix = this.patches[i];
      for (let j = 0; j < this.pdataLength; j++) {
        pmatrix[j] = (
          pdata[j * 4] * 0.3 +
          pdata[(j * 4) + 1] * 0.59 +
          pdata[(j * 4) + 2] * 0.11
        );
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
    if (this.patchType === 'SVM') {
      if (typeof this.webglFi !== 'undefined') {
        this.responses = this._getWebGLResponses(this.patches);
      } else if (typeof this.svmFi !== 'undefined') {
        this.responses = this.svmFi.getResponses(this.patches);
      } else {
        throw new Error('SVM-filters do not seem to be initiated properly.');
      }
    } else if (this.patchType === 'MOSSE') {
      this.responses = this.mosseCalc.getResponses(this.patches);
    }

    // option to increase sharpness of responses
    if (this.params.sharpenResponse) {
      for (let i = 0; i < this.numPatches; i++) {
        for (let j = 0; j < this.responses[i].length; j++) {
          this.responses[i][j] = Math.pow(
            this.responses[i][j],
            this.params.sharpenResponse
          );
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

      //if ([27,32,44,50].indexOf(i) > -1) {
      //  drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[i][0]-((searchWindow-1)/2), patchPositions[i][1]-((searchWindow-1)/2));
      //}
      drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[i][0]-((searchWindow-1)/2), patchPositions[i][1]-((searchWindow-1)/2));
    }*/

    // iterate until convergence or max 10, 20 iterations?:
    var jac;

    for (let i = 0; i < this.varianceSeq.length; i++) {
      // calculate jacobian
      jac = this._createJacobian(this.currentParameters, this.eigenVectors);

      // for debugging
      // var debugMVs = [];

      const meanshiftVectors = this._getMeanshiftVectors(scaling, i);

      // draw meanshiftVector
      /*sketchCC.clearRect(0, 0, sketchW, sketchH);
      var nuWeights;
      for (var npidx = 0;npidx < numPatches;npidx++) {
        nuWeights = debugMVs[npidx].map(function(x) {return x*255*500;});
        drawData(sketchCC, nuWeights, searchWindow, searchWindow, false, patchPositions[npidx][0]-((searchWindow-1)/2), patchPositions[npidx][1]-((searchWindow-1)/2));
      }*/

      var meanShiftVector = numeric.rep([this.numPatches * 2, 1], 0.0);
      for (let k = 0; k < this.numPatches; k++) {
        meanShiftVector[k * 2][0] = meanshiftVectors[k][0];
        meanShiftVector[(k * 2) + 1][0] = meanshiftVectors[k][1];
      }

      // compute pdm parameter update
      // var prior = numeric.mul(this.gaussianPD, this.PDMVariance);
      var prior = numeric.mul(this.gaussianPD, this.varianceSeq[i]);
      let jtj;
      if (this.params.weightPoints) {
        jtj = numeric.dot(numeric.transpose(jac), numeric.dot(this.pointWeights, jac));
      } else {
        jtj = numeric.dot(numeric.transpose(jac), jac);
      }
      var cpMatrix = numeric.rep([this.numParameters + 4, 1], 0.0);
      for (let l = 0; l < (this.numParameters + 4); l++) {
        cpMatrix[l][0] = this.currentParameters[l];
      }
      var priorP = numeric.dot(prior, cpMatrix);
      let jtv;
      if (this.params.weightPoints) {
        jtv = numeric.dot(numeric.transpose(jac), numeric.dot(this.pointWeights, meanShiftVector));
      } else {
        jtv = numeric.dot(numeric.transpose(jac), meanShiftVector);
      }
      var paramUpdateLeft = numeric.add(prior, jtj);
      var paramUpdateRight = numeric.sub(priorP, jtv);

      let paramUpdateLeftDet = numeric.det(paramUpdateLeft);
      if (paramUpdateLeftDet === 0) {
        console.warn(paramUpdateLeft);
        throw new Error('paramUpdateLeft is singular (determinate == 0), cannot invert');
      }
      if (isNaN(paramUpdateLeftDet)) {
        console.warn(paramUpdateLeft);
        throw new Error('paramUpdateLeft has invalid determinate (NaN)');
      }

      var paramUpdate = numeric.dot(numeric.inv(paramUpdateLeft), paramUpdateRight);
      // var paramUpdate = numeric.solve(paramUpdateLeft, paramUpdateRight, true);

      var oldPositions = this.currentPositions;

      // update estimated parameters
      for (let k = 0; k < this.numParameters + 4; k++) {
        this.currentParameters[k] -= paramUpdate[k];
      }

      // clipping of parameters if they're too high
      for (let k = 0; k < this.numParameters; k++) {
        const clip = Math.abs(3 * Math.sqrt(this.eigenValues[k]));
        if (Math.abs(this.currentParameters[k + 4]) > clip) {
          if (this.currentParameters[k + 4] > 0) {
            this.currentParameters[k + 4] = clip;
          } else {
            this.currentParameters[k + 4] = -clip;
          }
        }
      }

      // update current coordinates
      this.currentPositions = this._calculatePositions(this.currentParameters, true);

      // check if converged
      // calculate norm of parameterdifference
      let positionNorm = 0;
      let pnsqX;
      let pnsqY;
      for (let k = 0; k < this.currentPositions.length; k++) {
        pnsqX = (this.currentPositions[k][0] - oldPositions[k][0]);
        pnsqY = (this.currentPositions[k][1] - oldPositions[k][1]);
        positionNorm += ((pnsqX * pnsqX) + (pnsqY * pnsqY));
      }
      // console.log("positionnorm:"+positionNorm);

      // if norm < limit, then break
      if (positionNorm < this.convergenceLimit) {
        break;
      }
    }

    if (this.params.constantVelocity) {
      // add current parameter to array of previous parameters
      this.previousParameters.push(this.currentParameters.slice());
      if (this.previousParameters.length === 3) {
        this.previousParameters.shift();
      }
    }

    // store positions, for checking convergence
    if (this.previousPositions.length === 10) {
      this.previousPositions.shift();
    }
    this.previousPositions.push(this.currentPositions.slice(0));

    // send an event on each iteration
    this.emit('iteration');

    // we must get a score before we can say we've converged
    if (
      this.scoringHistory.length >= 5 &&
      this.getConvergence() < this.convergenceThreshold
    ) {
      if (this.params.stopOnConvergence) {
        this.stop();
      }
      this.emit('converged');
    }

    // return new points
    return this.currentPositions;
  }

  _getMeanshiftVectors (scaling, varianceIndex) {
    const originalPositions = this.currentPositions;
    const meanshiftVectors = this._meanshiftVectors;

    for (let j = 0; j < this.numPatches; j++) {
      const opj0 = originalPositions[j][0] - ((this.searchWindow - 1) * scaling / 2);
      const opj1 = originalPositions[j][1] - ((this.searchWindow - 1) * scaling / 2);

      // calculate PI x gaussians
      const vpsum = gpopt(
        this.searchWindow,
        originalPositions[j],
        this.updatePosition,
        this.vecProbs,
        this.responses,
        opj0,
        opj1,
        j,
        this.varianceSeq[varianceIndex],
        scaling
      );

      // calculate meanshift-vector
      gpopt2(
        this.searchWindow,
        this.vecpos,
        this.updatePosition,
        this.vecProbs,
        vpsum,
        opj0,
        opj1,
        scaling
      );

      // for debugging
      // var debugMatrixMV = gpopt2(searchWindow, vecpos, updatePosition, vecProbs, vpsum, opj0, opj1);

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

      meanshiftVectors[j][0] = this.vecpos[0] - originalPositions[j][0];
      meanshiftVectors[j][1] = this.vecpos[1] - originalPositions[j][1];

      // for debugging
      // debugMVs[j] = debugMatrixMV;
    }

    return meanshiftVectors;
  }

  _resetParameters () {
    this.first = true;
    this.scoringHistory = [];
    this.previousParameters = [];
    for (let i = 0, l = this.currentParameters.length; i < l; i++) {
      this.currentParameters[i] = 0;
    }
  }

  /*
   *  reset tracking, so that track() will start a new detection
   */
  reset () {
    this._resetParameters();
    this.runnerElement = undefined;
    this.runnerBox = undefined;
  }

  /*
   *  draw model on given canvas
   */
  draw (canvas, pv, path) {
    // if no previous points, just draw in the middle of canvas
    let params;
    if (pv === undefined) {
      params = this.currentParameters.slice(0);
    } else {
      params = pv.slice(0);
    }

    const cc = canvas.getContext('2d');
    cc.fillStyle = 'rgb(200,200,200)';
    cc.strokeStyle = 'rgb(130,255,50)';
    // cc.lineWidth = 1;

    let paths;
    if (path === undefined) {
      paths = this.model.path.normal;
    } else {
      paths = this.model.path[path];
    }

    for (let i = 0; i < paths.length; i++) {
      if (typeof paths[i] === 'number') {
        this._drawPoint(cc, paths[i], params);
      } else {
        this._drawPath(cc, paths[i], params);
      }
    }
  }

  /*
   *  get the score of the current model fit
   *  (based on svm of face according to current model)
   */
  getScore () {
    return this.meanscore;
  }

  /*
   *  calculate positions based on parameters
   */
  calculatePositions (parameters) {
    return this._calculatePositions(parameters, true);
  }

  /*
   *  get coordinates of current model fit
   */
  getCurrentPosition () {
    if (this.first) {
      return false;
    } else {
      return this.currentPositions;
    }
  }

  /*
   *  get parameters of current model fit
   */
  getCurrentParameters () {
    return this.currentParameters;
  }

  /*
   *  Get the average of recent model movements
   *  Used for checking whether model fit has converged
   */
  getConvergence () {
    if (this.previousPositions.length < 10) return 999999;

    let prevX = 0.0;
    let prevY = 0.0;
    let currX = 0.0;
    let currY = 0.0;

    // average 5 previous positions
    for (let i = 0; i < 5; i++) {
      for (let j = 0; j < this.numPatches; j++) {
        prevX += this.previousPositions[i][j][0];
        prevY += this.previousPositions[i][j][1];
      }
    }
    prevX /= 5;
    prevY /= 5;

    // average 5 positions before that
    for (let i = 5; i < 10; i++) {
      for (let j = 0; j < this.numPatches; j++) {
        currX += this.previousPositions[i][j][0];
        currY += this.previousPositions[i][j][1];
      }
    }
    currX /= 5;
    currY /= 5;

    // calculate difference
    const diffX = currX - prevX;
    const diffY = currY - prevY;
    let msavg = ((diffX * diffX) + (diffY * diffY));
    msavg /= this.previousPositions.length
    return msavg;
  }

  /*
   * Set response mode (only useful if webGL is available)
   * mode : either "single", "blend" or "cycle"
   * list : array of values "raw", "sobel", "lbp"
   */
  setResponseMode (mode, list) {
    // clmtrackr must be initialized with model first
    if (typeof this.model === 'undefined') {
      console.log('Clmtrackr has not been initialized with a model yet. No changes made.');
      return;
    }
    // must check whether webGL or not
    if (typeof this.webglFi === 'undefined') {
      console.log('Responsemodes are only allowed when using webGL. In pure JS, only "raw" mode is available.');
      return;
    }
    if (VALID_RESPONSEMODES.indexOf(mode) < 0) {
      console.log("Tried to set an unknown responsemode : '" + mode + "'. No changes made.");
      return;
    }
    if (!(list instanceof Array)) {
      console.log('List in setResponseMode must be an array of strings! No changes made.');
      return;
    } else {
      for (let i = 0; i < list.length; i++) {
        if (VALID_RESPONSEMODE_LIST.indexOf(list[i]) < 0) {
          console.log("Unknown element in responsemode list : '" + list[i] + "'. No changes made.");
        }
        // check whether filters are initialized
        if (list[i] === 'sobel' && this.sobelInit === false) {
          console.log('The sobel filters have not been initialized! No changes made.');
        }
        if (list[i] === 'lbp' && this.lbpInit === false) {
          console.log('The LBP filters have not been initialized! No changes made.');
        }
      }
    }
    // reset index
    this.responseIndex = 0;
    this.responseMode = mode;
    this.responseList = list;
  }

  _getWebGLResponsesType (type, patches) {
    if (type === 'lbp') {
      return this.webglFi.getLBPResponses(patches);
    } else if (type === 'raw') {
      return this.webglFi.getRawResponses(patches);
    } else if (type === 'sobel') {
      return this.webglFi.getSobelResponses(patches);
    }
  }

  _getWebGLResponses (patches) {
    if (this.responseMode === 'single') {
      return this._getWebGLResponsesType(this.responseList[0], patches);
    } else if (this.responseMode === 'cycle') {
      const response = this._getWebGLResponsesType(this.responseList[this.responseIndex], patches);
      this.responseIndex++;
      if (this.responseIndex >= this.responseList.length) this.responseIndex = 0;
      return response;
    } else {
      // blend
      const responses = [];
      for (let i = 0; i < this.responseList.length; i++) {
        responses[i] = this._getWebGLResponsesType(this.responseList[i], patches);
      }
      const blendedResponses = [];
      const searchWindowSize = this.searchWindow * this.searchWindow;
      for (let i = 0; i < this.numPatches; i++) {
        const response = Array(searchWindowSize);
        for (let k = searchWindowSize - 1; k >= 0; k--) {
          response[k] = 0;
        }

        for (let j = 0; j < this.responseList.length; j++) {
          for (let k = 0; k < searchWindowSize; k++) {
            response[k] += (responses[j][i][k] / this.responseList.length);
          }
        }
        blendedResponses[i] = response;
      }
      return blendedResponses;
    }
  }

  // generates the jacobian matrix used for optimization calculations
  _createJacobian (parameters, eigenVectors) {
    const jacobian = numeric.rep(
      [2 * this.numPatches, this.numParameters + 4],
      0.0
    );
    let j0;
    let j1;
    for (let i = 0; i < this.numPatches; i++) {
      // 1
      j0 = this.meanShape[i][0];
      j1 = this.meanShape[i][1];
      for (let p = 0; p < this.numParameters; p++) {
        j0 += parameters[p + 4] * eigenVectors[i * 2][p];
        j1 += parameters[p + 4] * eigenVectors[(i * 2) + 1][p];
      }
      jacobian[i * 2][0] = j0;
      jacobian[(i * 2) + 1][0] = j1;
      // 2
      j0 = this.meanShape[i][1];
      j1 = this.meanShape[i][0];
      for (let p = 0; p < this.numParameters; p++) {
        j0 += parameters[p + 4] * eigenVectors[(i * 2) + 1][p];
        j1 += parameters[p + 4] * eigenVectors[i * 2][p];
      }
      jacobian[i * 2][1] = -j0;
      jacobian[(i * 2) + 1][1] = j1;
      // 3
      jacobian[i * 2][2] = 1;
      jacobian[(i * 2) + 1][2] = 0;
      // 4
      jacobian[i * 2][3] = 0;
      jacobian[(i * 2) + 1][3] = 1;
      // the rest
      for (let j = 0; j < this.numParameters; j++) {
        j0 = (
          parameters[0] *
          eigenVectors[i * 2][j] - parameters[1] *
          eigenVectors[(i * 2) + 1][j] + eigenVectors[i * 2][j]
        );
        j1 = (
          parameters[0] *
          eigenVectors[(i * 2) + 1][j] + parameters[1] *
          eigenVectors[i * 2][j] + eigenVectors[(i * 2) + 1][j]
        );
        jacobian[i * 2][j + 4] = j0;
        jacobian[(i * 2) + 1][j + 4] = j1;
      }
    }

    return jacobian;
  }

  // calculate positions from parameters
  _calculatePositions (parameters, useTransforms) {
    var x, y, a, b;
    var numParameters = parameters.length;
    var positions = [];
    for (let i = 0; i < this.numPatches; i++) {
      x = this.meanShape[i][0];
      y = this.meanShape[i][1];
      for (let j = 0; j < numParameters - 4; j++) {
        x += this.model.shapeModel.eigenVectors[(i * 2)][j] * parameters[j + 4];
        y += this.model.shapeModel.eigenVectors[(i * 2) + 1][j] * parameters[j + 4];
      }
      if (useTransforms) {
        a = parameters[0] * x - parameters[1] * y + parameters[2];
        b = parameters[0] * y + parameters[1] * x + parameters[3];
        x += a;
        y += b;
      }
      positions[i] = [x, y];
    }

    return positions;
  }

  _faceDetected (comp, callback) {
    if (comp && comp.length > 0) {
      this.candidate = comp[0];
    } else {
      callback(false);
      return false;
    }

    for (var i = 1; i < comp.length; i++) {
      if (comp[i].confidence > this.candidate.confidence) {
        this.candidate = comp[i];
      }
    }

    // return candidate;
    callback(this.candidate);
  }

  // detect position of face on canvas/video element
  _detectPosition (el, callback) {
    const detectedCb = (comp) => {
      this.jsfeatFace.removeListener('faceDetected', detectedCb);
      this._faceDetected(comp, callback);
    };

    this.jsfeatFace.on('faceDetected', detectedCb);

    // TODO Allow option that limit simultaneous trigger of WebWorkers
    this.jsfeatFace.findFace(el);
  }

  // calculate score of current fit
  _checkTracking () {
    const trackingImgW = 20;
    const trackingImgH = 22;

    this.scoringContext.drawImage(
      this.sketchCanvas,
      Math.round(this.msxmin + (this.msmodelwidth / 4.5)),
      Math.round(this.msymin - (this.msmodelheight / 12)),
      Math.round(this.msmodelwidth - (this.msmodelwidth * 2 / 4.5)),
      Math.round(this.msmodelheight - (this.msmodelheight / 12)),
      0,
      0,
      trackingImgW,
      trackingImgH
    );
    // getImageData of canvas
    var imgData = this.scoringContext.getImageData(0, 0, trackingImgW, trackingImgH);
    // convert data to grayscale
    const trackingImgSize = trackingImgW * trackingImgH;
    var scoringData = new Array(trackingImgSize);
    var scdata = imgData.data;
    var scmax = 0;
    for (let i = 0; i < trackingImgSize; i++) {
      scoringData[i] = (
        scdata[i * 4] * 0.3 +
        scdata[(i * 4) + 1] * 0.59 +
        scdata[(i * 4) + 2] * 0.11
      );
      scoringData[i] = Math.log(scoringData[i] + 1);
      if (scoringData[i] > scmax) {
        scmax = scoringData[i];
      }
    }

    if (scmax > 0) {
      // normalize & multiply by svmFilter
      const mean = this._calcMean(scoringData, trackingImgSize);
      const sd = this._calcStandardDeviation(scoringData, trackingImgSize, mean);

      var score = 0;
      for (let i = 0; i < trackingImgSize; i++) {
        scoringData[i] = (scoringData[i] - mean) / sd;
        score += scoringData[i] * this.scoringWeights[i];
      }
      score += this.scoringBias;
      score = 1 / (1 + Math.exp(-score));

      // Keep length === 5
      if (this.scoringHistory.length === 5) {
        this.scoringHistory.shift();
      }
      this.scoringHistory.push(score);

      if (this.scoringHistory.length > 4) {
        // get average
        this.meanscore = this._calcMean(this.scoringHistory, 5);
        // if below threshold, then reset (return false)
        if (this.meanscore < this.params.scoreThreshold) {
          return false;
        }
      }
    }
    return true;
  }

  _calcMean (data, length) {
    let mean = 0;
    for (let i = 0; i < length; i++) {
      mean += data[i];
    }
    return mean / length;
  }

  _calcStandardDeviation (data, length, mean) {
    let avgSquareDiff = 0;
    for (let i = 0; i < length; i++) {
      const diff = data[i] - mean;
      avgSquareDiff += diff * diff;
    }
    avgSquareDiff /= length;
    return Math.sqrt(avgSquareDiff);
  }

  // get initial starting point for model
  _getInitialPosition (element, box, callback, det) {
    var translateX, translateY, scaling, rotation;
    if (box) {
      this.candidate = {
        x: box[0],
        y: box[1],
        width: box[2],
        height: box[3]
      };
    } else {
      if (!det) {
        this._detectPosition(element, (det) => {
          if (!det) {
            // if no face found, stop.
            callback(false);
          } else {
            this._getInitialPosition(element, box, callback, det);
          }
        });
        return;
      }
    }

    const candidate = this.candidate;
    if (this.model.hints && HAS_MOSSE_FILTERS) {
      var noseFilterWidth = candidate.width * 4.5 / 10;
      var eyeFilterWidth = candidate.width * 6 / 10;

      // detect position of eyes and nose via mosse filter
      //
      /*element.pause();

      var canvasContext = document.getElementById('overlay2').getContext('2d')
      canvasContext.clearRect(0,0,500,375);
      canvasContext.strokeRect(candidate.x, candidate.y, candidate.width, candidate.height);*/
      //

      var noseResult = this.mossef_nose.track(
        element,
        Math.round(candidate.x + (candidate.width / 2) - (noseFilterWidth / 2)),
        Math.round(candidate.y + candidate.height * (5 / 8) - (noseFilterWidth / 2)),
        noseFilterWidth,
        noseFilterWidth,
        false
      );
      var rightResult = this.mossef_righteye.track(
        element,
        Math.round(candidate.x + (candidate.width * 3 / 4) - (eyeFilterWidth / 2)),
        Math.round(candidate.y + candidate.height * (2 / 5) - (eyeFilterWidth / 2)),
        eyeFilterWidth,
        eyeFilterWidth,
        false
      );
      var leftResult = this.mossef_lefteye.track(
        element,
        Math.round(candidate.x + (candidate.width / 4) - (eyeFilterWidth / 2)),
        Math.round(candidate.y + candidate.height * (2 / 5) - (eyeFilterWidth / 2)),
        eyeFilterWidth,
        eyeFilterWidth,
        false
      );
      this.right_eye_position[0] = Math.round(candidate.x + (candidate.width * 3 / 4) - (eyeFilterWidth / 2)) + rightResult[0];
      this.right_eye_position[1] = Math.round(candidate.y + candidate.height * (2 / 5) - (eyeFilterWidth / 2)) + rightResult[1];
      this.left_eye_position[0] = Math.round(candidate.x + (candidate.width / 4) - (eyeFilterWidth / 2)) + leftResult[0];
      this.left_eye_position[1] = Math.round(candidate.y + candidate.height * (2 / 5) - (eyeFilterWidth / 2)) + leftResult[1];
      this.nose_position[0] = Math.round(candidate.x + (candidate.width / 2) - (noseFilterWidth / 2)) + noseResult[0];
      this.nose_position[1] = Math.round(candidate.y + candidate.height * (5 / 8) - (noseFilterWidth / 2)) + noseResult[1];

      //
      /*canvasContext.strokeRect(Math.round(candidate.x+(candidate.width*3/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
      canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/4)-(eyeFilterWidth/2)), Math.round(candidate.y+candidate.height*(2/5)-(eyeFilterWidth/2)), eyeFilterWidth, eyeFilterWidth);
      //canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(3/4)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth);
      canvasContext.strokeRect(Math.round(candidate.x+(candidate.width/2)-(noseFilterWidth/2)), Math.round(candidate.y+candidate.height*(5/8)-(noseFilterWidth/2)), noseFilterWidth, noseFilterWidth);

      canvasContext.fillStyle = "rgb(0,0,250)";
      canvasContext.beginPath();
      canvasContext.arc(this.left_eye_position[0], this.left_eye_position[1], 3, 0, Math.PI*2, true);
      canvasContext.closePath();
      canvasContext.fill();

      canvasContext.beginPath();
      canvasContext.arc(this.right_eye_position[0], this.right_eye_position[1], 3, 0, Math.PI*2, true);
      canvasContext.closePath();
      canvasContext.fill();

      canvasContext.beginPath();
      canvasContext.arc(this.nose_position[0], this.nose_position[1], 3, 0, Math.PI*2, true);
      canvasContext.closePath();
      canvasContext.fill();

      debugger;
      element.play()
      canvasContext.clearRect(0,0,element.width,element.height);*/
      //

      // get eye and nose positions of model
      var lep = this.model.hints.leftEye;
      var rep = this.model.hints.rightEye;
      var mep = this.model.hints.nose;

      // get scaling, rotation, etc. via procrustes analysis
      var procrustes_params = procrustes([this.left_eye_position, this.right_eye_position, this.nose_position], [lep, rep, mep]);
      translateX = procrustes_params[0];
      translateY = procrustes_params[1];
      scaling = procrustes_params[2];
      rotation = procrustes_params[3];

      // element.play();

      // var maxscale = 1.10;
      // if ((scaling*modelHeight)/candidate.height < maxscale*0.7) scaling = (maxscale*0.7*candidate.height)/modelHeight;
      // if ((scaling*modelHeight)/candidate.height > maxscale*1.2) scaling = (maxscale*1.2*candidate.height)/modelHeight;

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

      this.currentParameters[0] = (scaling * Math.cos(rotation)) - 1;
      this.currentParameters[1] = (scaling * Math.sin(rotation));
      this.currentParameters[2] = translateX;
      this.currentParameters[3] = translateY;

      // this.draw(document.getElementById('overlay'), this.currentParameters);
    } else {
      scaling = candidate.width / this.modelheight;
      // var ccc = document.getElementById('overlay').getContext('2d');
      // ccc.strokeRect(candidate.x,candidate.y,candidate.width,candidate.height);
      // FIXME: xmin and ymin are undefined
      translateX = candidate.x - (xmin * scaling) + 0.1 * candidate.width;
      translateY = candidate.y - (ymin * scaling) + 0.25 * candidate.height;
      this.currentParameters[0] = scaling - 1;
      this.currentParameters[2] = translateX;
      this.currentParameters[3] = translateY;
    }

    this.currentPositions = this._calculatePositions(this.currentParameters, true);

    callback([scaling, rotation, translateX, translateY]);
  }

  // draw a parametrized line on a canvas
  _drawPath (canvasContext, path, dp) {
    canvasContext.beginPath();
    var i, x, y, a, b;
    for (let p = 0; p < path.length; p++) {
      i = path[p] * 2;
      x = this.meanShape[i / 2][0];
      y = this.meanShape[i / 2][1];
      for (let j = 0; j < this.numParameters; j++) {
        x += this.model.shapeModel.eigenVectors[i][j] * dp[j + 4];
        y += this.model.shapeModel.eigenVectors[i + 1][j] * dp[j + 4];
      }
      a = dp[0] * x - dp[1] * y + dp[2];
      b = dp[0] * y + dp[1] * x + dp[3];
      x += a;
      y += b;

      if (i === 0) {
        canvasContext.moveTo(x, y);
      } else {
        canvasContext.lineTo(x, y);
      }
    }
    canvasContext.moveTo(0, 0);
    canvasContext.closePath();
    canvasContext.stroke();
  }

  // draw a point on a canvas
  _drawPoint (canvasContext, point, dp) {
    var i, x, y, a, b;
    i = point * 2;
    x = this.meanShape[i / 2][0];
    y = this.meanShape[i / 2][1];
    for (let j = 0; j < this.numParameters; j++) {
      x += this.model.shapeModel.eigenVectors[i][j] * dp[j + 4];
      y += this.model.shapeModel.eigenVectors[i + 1][j] * dp[j + 4];
    }
    a = dp[0] * x - dp[1] * y + dp[2];
    b = dp[0] * y + dp[1] * x + dp[3];
    x += a;
    y += b;
    canvasContext.beginPath();
    canvasContext.arc(x, y, 1, 0, Math.PI * 2, true);
    canvasContext.closePath();
    canvasContext.fill();
  }
}
