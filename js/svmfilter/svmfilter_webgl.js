import {
  setupWebGL,
  loadShader,
  loadProgram
} from 'clmtrackr/js/utils/webgl';

import createLbpResponseVS from './shaders/lbpResponse.vert';
import createLbpResponseFS from './shaders/lbpResponse.frag';

import createGradientResponseVS from './shaders/gradientResponse.vert';
import createGradientResponseFS from './shaders/gradientResponse.frag';

import createDrawResponsesVS from './shaders/drawResponses.vert';
import createDrawResponsesFS from './shaders/drawResponses.frag';

import createPatchResponseVS from './shaders/patchResponse.vert';
import createPatchResponseFS from './shaders/patchResponse.frag';


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
export default class WebglFilter {
  constructor () {
    this.gl;
    this.canvas;

    this.filterWidth;
    this.filterHeight;
    this.patchWidth;
    this.patchHeight;
    this.numPatches;
    this.canvasWidth;
    this.canvasHeight;
    //
    this.patchResponseProgram;
    this.patchDrawProgram;

    this.fbo;
    this.numBlocks;
    this.patchTex;

    this.drawRectBuffer;
    this.drawLayerBuffer;
    this.drawImageBuffer;
    this.rttTexture;

    this.texCoordBuffer;
    this.texCoordLocation;
    this.apositionBuffer;

    this.newCanvasWidth;
    this.newCanvasBlockHeight;
    this.newCanvasHeight;

    this.drawOutRectangles;
    this.drawOutImages;
    this.drawOutLayer;

    this.patchCells;
    this.textureWidth;
    this.textureHeight;
    this.patchSize;
    this.patchArray;

    this.biases;
    //
    this.lbpResponseProgram;

    this.lbo;
    this.lbpTexCoordLocation;
    this.lbpTexCoordBuffer;
    this.lbpPositionLocation;
    this.lbpAPositionBuffer;
    //
    this.gradientResponseProgram;

    this.gbo;
    this.gradTexCoordLocation;
    this.gradTexCoordBuffer;
    this.gradPositionLocation;
    this.gradAPositionBuffer;

    this.lbpInit = false;
    this.sobelInit = false;
    this.rawInit = false;

    this.lbpResponseVS = createLbpResponseVS();
    this.lbpResponseFS;

    this.gradientResponseVS = createGradientResponseVS();
    this.gradientResponseFS;

    this.patchResponseVS;
    this.patchResponseFS;

    this.drawResponsesVS = createDrawResponsesVS();
    this.drawResponsesFS = createDrawResponsesFS()
  }

  init (filters, bias, nP, pW, pH, fW, fH) {
    // we assume filterVector goes from left to right, rowwise, i.e. row-major order
    if (fW !== fH) {
      alert('filter width and height must be same size!');
      return;
    }

    // if filter width is not odd, alert
    if (fW % 2 === 0 || fH % 2 === 0) {
      alert('filters used in svm must be of odd dimensions!');
      return;
    }

    // setup variables
    this.biases = bias;
    this.filterWidth = fW;
    this.filterHeight = fH;
    this.patchWidth = pW;
    this.patchHeight = pH;
    this.numPatches = nP;
    this.numBlocks = Math.floor(this.numPatches / 4) + Math.ceil((this.numPatches % 4) / 4);
    this.canvasWidth = this.patchWidth;
    this.canvasHeight = this.patchHeight * this.numBlocks;
    this.newCanvasWidth = this.patchWidth - this.filterWidth + 1;
    this.newCanvasBlockHeight = this.patchHeight - this.filterWidth + 1;
    this.newCanvasHeight = this.newCanvasBlockHeight * this.numPatches;
    this.patchCells = (Math.floor(this.numPatches / 4) + Math.ceil((this.numPatches % 4) / 4));
    this.textureWidth = this.patchWidth;
    this.textureHeight = this.patchHeight * this.patchCells;
    this.patchSize = this.patchWidth * this.patchHeight;
    this.patchArray = new Float32Array(this.patchSize * this.patchCells * 4);
    let opp = [1 / this.patchWidth, 1 / (this.patchHeight * this.numBlocks)];

    // write out shaders
    this.patchResponseFS = createPatchResponseFS({
      onePixelPatchesX: (1 / this.patchWidth).toFixed(10),
      onePixelPatchesY: (1 / (this.patchHeight * this.numBlocks)).toFixed(10),
      onePixelFiltersX: (1 / this.filterWidth).toFixed(10),
      onePixelFiltersY: (1 / (this.filterHeight * this.numBlocks)).toFixed(10),
      halfFilterWidth: ((this.filterWidth - 1.0) / 2).toFixed(1),
      halfFilterHeight: ((this.filterHeight - 1.0) / 2).toFixed(1),
      filterWidth: this.filterWidth,
      filterHeight: this.filterHeight
    });

    this.patchResponseVS = createPatchResponseVS({
      resolutionX: this.canvasWidth.toFixed(1),
      resolutionY: this.canvasHeight.toFixed(1),
      patchHeight: (1 / this.numBlocks).toFixed(10),
      filterHeight: (1 / this.numBlocks).toFixed(10),
      midpointY: (1 / (this.numBlocks * 2)).toFixed(10)
    });

    if ('lbp' in filters) {
      // lbpResponseFragment
      this.lbpResponseFS = createLbpResponseFS({
        opp0: opp[0].toFixed(5),
        opp1: opp[1].toFixed(5)
      });
    }

    if ('sobel' in filters) {
      // gradResponseFragment
      this.gradientResponseFS = createGradientResponseFS({
        opp0: opp[0].toFixed(5),
        opp1: opp[1].toFixed(5)
      });
    }

    // create webglcanvas
    this.canvas = document.createElement('canvas')
    this.canvas.setAttribute('width', (this.patchWidth - this.filterWidth + 1) + 'px');
    this.canvas.setAttribute('height', ((this.patchHeight - this.filterHeight + 1) * this.numPatches) + 'px');
    this.canvas.setAttribute('id', 'renderCanvas');
    this.canvas.setAttribute('style', 'display:none;');
    // document.body.appendChild(this.canvas);
    const gl = setupWebGL(this.canvas, {
      premultipliedAlpha: false,
      preserveDrawingBuffer: true,
      antialias: false
    });
    this.gl = gl;


    // check for float textures support and fail if not
    if (!gl.getExtension('OES_texture_float')) {
      alert('Your graphics card does not support floating point textures! :(');
      return;
    }

    /** insert filters into textures **/
    if ('raw' in filters) {
      this._insertFilter(filters['raw'], gl.TEXTURE0)
      this.rawInit = true;
    }
    if ('sobel' in filters) {
      this._insertFilter(filters['sobel'], gl.TEXTURE4)
      this.sobelInit = true;
    }
    if ('lbp' in filters) {
      this._insertFilter(filters['lbp'], gl.TEXTURE5)
      this.lbpInit = true;
    }

    /** calculate vertices for calculating responses **/

    // vertex rectangles to draw out
    var rectangles = [];
    var halfFilter = (this.filterWidth - 1) / 2;
    let yOffset;
    for (let i = 0; i < this.numBlocks; i++) {
      yOffset = i * this.patchHeight;
      // first triangle
      rectangles = rectangles.concat([
        halfFilter, yOffset + halfFilter,
        this.patchWidth - halfFilter, yOffset + halfFilter,
        halfFilter, yOffset + this.patchHeight - halfFilter
      ]);
      // second triangle
      rectangles = rectangles.concat([
        halfFilter, yOffset + this.patchHeight - halfFilter,
        this.patchWidth - halfFilter, yOffset + halfFilter,
        this.patchWidth - halfFilter, yOffset + this.patchHeight - halfFilter
      ]);
    }
    rectangles = new Float32Array(rectangles);

    // image rectangles to draw out
    var irectangles = [];
    for (let i = 0; i < rectangles.length; i++) {
      if (i % 2 === 0) {
        irectangles[i] = rectangles[i] / this.canvasWidth;
      } else {
        irectangles[i] = rectangles[i] / this.canvasHeight;
      }
    }
    irectangles = new Float32Array(irectangles);

    if ('lbp' in filters || 'sobel' in filters) {
      var topCoord = 1.0 - 2 / (this.patchHeight * this.numBlocks);
      var bottomCoord = 1.0 - 2 / this.numBlocks + 2 / (this.patchHeight * this.numBlocks);
      // calculate position of vertex rectangles for gradient/lbp program
      var gradRectangles = [];
      for (let i = 0; i < this.numBlocks; i++) {
        yOffset = i * (2 / this.numBlocks);
        // first triangle
        gradRectangles = gradRectangles.concat(
          [-1.0, topCoord - yOffset,
          1.0, topCoord - yOffset,
          -1.0, bottomCoord - yOffset]
        );
        // second triangle
        gradRectangles = gradRectangles.concat(
          [-1.0, bottomCoord - yOffset,
          1.0, topCoord - yOffset,
          1.0, bottomCoord - yOffset]
        );
      }
      gradRectangles = new Float32Array(gradRectangles);

      topCoord = 1.0 - 1 / (this.patchHeight * this.numBlocks);
      bottomCoord = 1.0 - 1 / this.numBlocks + 1 / (this.patchHeight * this.numBlocks);
      // calculate position of image rectangles to draw out
      var gradIRectangles = [];
      for (let i = 0; i < this.numBlocks; i++) {
        yOffset = i * (1 / this.numBlocks);
        // first triangle
        gradIRectangles = gradIRectangles.concat(
          [0.0, topCoord - yOffset,
          1.0, topCoord - yOffset,
          0.0, bottomCoord - yOffset]
        );
        // second triangle
        gradIRectangles = gradIRectangles.concat(
          [0.0, bottomCoord - yOffset,
          1.0, topCoord - yOffset,
          1.0, bottomCoord - yOffset]
        );
      }
      gradIRectangles = new Float32Array(gradIRectangles);
    }

    // vertices for drawing out responses

    // drawOutRectangles
    this.drawOutRectangles = new Float32Array(12 * this.numPatches);
    let indexOffset;
    for (var i = 0; i < this.numPatches; i++) {
      yOffset = i * this.newCanvasBlockHeight;
      indexOffset = i * 12;

      // first triangle
      this.drawOutRectangles[indexOffset] = 0.0;
      this.drawOutRectangles[indexOffset + 1] = yOffset;
      this.drawOutRectangles[indexOffset + 2] = this.newCanvasWidth;
      this.drawOutRectangles[indexOffset + 3] = yOffset;
      this.drawOutRectangles[indexOffset + 4] = 0.0;
      this.drawOutRectangles[indexOffset + 5] = yOffset + this.newCanvasBlockHeight;

      // second triangle
      this.drawOutRectangles[indexOffset + 6] = 0.0;
      this.drawOutRectangles[indexOffset + 7] = yOffset + this.newCanvasBlockHeight;
      this.drawOutRectangles[indexOffset + 8] = this.newCanvasWidth;
      this.drawOutRectangles[indexOffset + 9] = yOffset;
      this.drawOutRectangles[indexOffset + 10] = this.newCanvasWidth;
      this.drawOutRectangles[indexOffset + 11] = yOffset + this.newCanvasBlockHeight;
    }

    // images
    this.drawOutImages = new Float32Array(this.numPatches * 12);
    var halfFilterWidth = ((this.filterWidth - 1) / 2) / this.patchWidth;
    var halfFilterHeight = ((this.filterWidth - 1) / 2) / (this.patchHeight * this.patchCells);
    var patchHeightT = this.patchHeight / (this.patchHeight * this.patchCells);
    for (let i = 0; i < this.numPatches; i++) {
      yOffset = Math.floor(i / 4) * patchHeightT;
      indexOffset = i * 12;

      // first triangle
      this.drawOutImages[indexOffset] = halfFilterWidth;
      this.drawOutImages[indexOffset + 1] = yOffset + halfFilterHeight;
      this.drawOutImages[indexOffset + 2] = 1.0 - halfFilterWidth;
      this.drawOutImages[indexOffset + 3] = yOffset + halfFilterHeight;
      this.drawOutImages[indexOffset + 4] = halfFilterWidth;
      this.drawOutImages[indexOffset + 5] = yOffset + patchHeightT - halfFilterHeight;

      // second triangle
      this.drawOutImages[indexOffset + 6] = halfFilterWidth;
      this.drawOutImages[indexOffset + 7] = yOffset + patchHeightT - halfFilterHeight;
      this.drawOutImages[indexOffset + 8] = 1.0 - halfFilterWidth;
      this.drawOutImages[indexOffset + 9] = yOffset + halfFilterHeight;
      this.drawOutImages[indexOffset + 10] = 1.0 - halfFilterWidth;
      this.drawOutImages[indexOffset + 11] = yOffset + patchHeightT - halfFilterHeight;
    }

    // layer
    this.drawOutLayer = new Float32Array(this.numPatches * 6);
    var layernum;
    for (let i = 0; i < this.numPatches; i++) {
      layernum = i % 4;
      indexOffset = i * 6;
      this.drawOutLayer[indexOffset] = layernum;
      this.drawOutLayer[indexOffset + 1] = layernum;
      this.drawOutLayer[indexOffset + 2] = layernum;
      this.drawOutLayer[indexOffset + 3] = layernum;
      this.drawOutLayer[indexOffset + 4] = layernum;
      this.drawOutLayer[indexOffset + 5] = layernum;
    }

    /** set up programs and load attributes etc **/

    if ('sobel' in filters) {
      var grVertexShader = loadShader(gl, this.gradientResponseVS, gl.VERTEX_SHADER);
      var grFragmentShader = loadShader(gl, this.gradientResponseFS, gl.FRAGMENT_SHADER);
      this.gradientResponseProgram = loadProgram(gl, [grVertexShader, grFragmentShader]);
      gl.useProgram(this.gradientResponseProgram);

      // set up vertices with rectangles
      this.gradPositionLocation = gl.getAttribLocation(this.gradientResponseProgram, 'a_position');
      this.gradAPositionBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.gradAPositionBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, gradRectangles, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(this.gradPositionLocation);
      gl.vertexAttribPointer(this.gradPositionLocation, 2, gl.FLOAT, false, 0, 0);

      // set up texture positions
      this.gradTexCoordLocation = gl.getAttribLocation(this.gradientResponseProgram, 'a_texCoord');
      this.gradTexCoordBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.gradTexCoordBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, gradIRectangles, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(this.gradTexCoordLocation);
      gl.vertexAttribPointer(this.gradTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

      // set up patches texture in gradientResponseProgram
      gl.uniform1i(gl.getUniformLocation(this.gradientResponseProgram, 'u_patches'), 1);
    }
    if ('lbp' in filters) {
      var lbpVertexShader = loadShader(gl, this.lbpResponseVS, gl.VERTEX_SHADER);
      var lbpFragmentShader = loadShader(gl, this.lbpResponseFS, gl.FRAGMENT_SHADER);
      this.lbpResponseProgram = loadProgram(gl, [lbpVertexShader, lbpFragmentShader]);
      gl.useProgram(this.lbpResponseProgram);

      // set up vertices with rectangles
      this.lbpPositionLocation = gl.getAttribLocation(this.lbpResponseProgram, 'a_position');
      this.lbpAPositionBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.lbpAPositionBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, gradRectangles, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(this.lbpPositionLocation);
      gl.vertexAttribPointer(this.lbpPositionLocation, 2, gl.FLOAT, false, 0, 0);

      // set up texture positions
      this.gradTexCoordLocation = gl.getAttribLocation(this.lbpResponseProgram, 'a_texCoord');
      this.lbpTexCoordBuffer = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.lbpTexCoordBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, gradIRectangles, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(this.lbpTexCoordLocation);
      gl.vertexAttribPointer(this.lbpTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

      // set up patches texture in lbpResponseProgram
      gl.uniform1i(gl.getUniformLocation(this.lbpResponseProgram, 'u_patches'), 1);
    }

    // setup patchdraw program
    var drVertexShader = loadShader(gl, this.drawResponsesVS, gl.VERTEX_SHADER);
    var drFragmentShader = loadShader(gl, this.drawResponsesFS, gl.FRAGMENT_SHADER);
    this.patchDrawProgram = loadProgram(gl, [drVertexShader, drFragmentShader]);
    gl.useProgram(this.patchDrawProgram);

    // set the resolution/dimension of the canvas
    var resolutionLocation = gl.getUniformLocation(this.patchDrawProgram, 'u_resolutiondraw');
    gl.uniform2f(resolutionLocation, this.newCanvasWidth, this.newCanvasHeight);

    // set u_responses
    var responsesLocation = gl.getUniformLocation(this.patchDrawProgram, 'u_responses');
    gl.uniform1i(responsesLocation, 2);

    // setup patchresponse program
    var prVertexShader = loadShader(gl, this.patchResponseVS, gl.VERTEX_SHADER);
    var prFragmentShader = loadShader(gl, this.patchResponseFS, gl.FRAGMENT_SHADER);
    this.patchResponseProgram = loadProgram(gl, [prVertexShader, prFragmentShader]);
    gl.useProgram(this.patchResponseProgram);

    // set up vertices with rectangles
    var positionLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_position');
    this.apositionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.apositionBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, rectangles, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    this.texCoordLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_texCoord');
    this.texCoordBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.texCoordBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, irectangles, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(this.texCoordLocation);
    gl.vertexAttribPointer(this.texCoordLocation, 2, gl.FLOAT, false, 0, 0);

    if ('lbp' in filters || 'sobel' in filters) {
      // set up gradient/lbp buffer (also used for lbp)
      gl.activeTexture(gl.TEXTURE3);
      var gradients = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, gradients);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, this.patchWidth, this.patchHeight * this.numBlocks, 0, gl.RGBA, gl.FLOAT, null);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

      // set up gradient/lbp framebuffer
      this.gbo = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.gbo);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, gradients, 0);
    }

    // set up buffer to draw to
    gl.activeTexture(gl.TEXTURE2);
    this.rttTexture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.rttTexture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, this.patchWidth, this.patchHeight * this.numBlocks, 0, gl.RGBA, gl.FLOAT, null);

    // set up response framebuffer
    this.fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.rttTexture, 0);

    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    /* initialize some textures and buffers used later on */

    this.patchTex = gl.createTexture();
    this.drawRectBuffer = gl.createBuffer();
    this.drawImageBuffer = gl.createBuffer();
    this.drawLayerBuffer = gl.createBuffer();
  }

  getRawResponses (patches) {
    // TODO: check patches correct length/dimension

    this._insertPatches(patches);
    const gl = this.gl;

    // switch to correct program
    gl.useProgram(this.patchResponseProgram);

    // set u_patches to point to texture 1
    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_patches'), 1);

    // set u_filters to point to correct filter
    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_filters'), 0);

    // set up vertices with rectangles
    var positionLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_position');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.apositionBuffer);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    var texCoordLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_texCoord');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.texCoordBuffer);
    gl.enableVertexAttribArray(texCoordLocation);
    gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

    // set framebuffer to the original one if not already using it
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);

    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, this.patchCells * 6);

    // gl.finish();

    var responses = this._drawOut('raw');

    return responses;
  }

  getSobelResponses (patches) {
    // check that it is initialized
    if (!this.sobelInit) return;

    this._insertPatches(patches);

    /* do sobel filter on patches */

    const gl = this.gl;
    // switch to correct program
    gl.useProgram(this.gradientResponseProgram);

    // set up vertices with rectangles
    var gradPositionLocation = gl.getAttribLocation(this.gradientResponseProgram, 'a_position');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.gradAPositionBuffer);
    gl.enableVertexAttribArray(gradPositionLocation);
    gl.vertexAttribPointer(gradPositionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    var gradTexCoordLocation = gl.getAttribLocation(this.gradientResponseProgram, 'a_texCoord');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.gradTexCoordBuffer);
    gl.enableVertexAttribArray(gradTexCoordLocation);
    gl.vertexAttribPointer(gradTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

    // set framebuffer to the original one if not already using it
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.gbo);

    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, this.patchCells * 6);

    /* calculate responses */

    gl.useProgram(this.patchResponseProgram);

    // set patches and filters to point to correct textures
    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_filters'), 4);
    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_patches'), 3);

    var positionLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_position');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.apositionBuffer);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    var texCoordLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_texCoord');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.texCoordBuffer);
    gl.enableVertexAttribArray(texCoordLocation);
    gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, this.patchCells * 6);

    /* get the responses */

    var responses = this._drawOut('sobel');

    return responses;
  }

  getLBPResponses (patches) {
    // check that it is initialized
    if (!this.lbpInit) return;

    this._insertPatches(patches);

    /* do sobel filter on patches */

    const gl = this.gl;
    // switch to correct program
    gl.useProgram(this.lbpResponseProgram);

    // set up vertices with rectangles
    var lbpPositionLocation = gl.getAttribLocation(this.lbpResponseProgram, 'a_position');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.lbpAPositionBuffer);
    gl.enableVertexAttribArray(lbpPositionLocation);
    gl.vertexAttribPointer(lbpPositionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    var lbpTexCoordLocation = gl.getAttribLocation(this.lbpResponseProgram, 'a_texCoord');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.lbpTexCoordBuffer);
    gl.enableVertexAttribArray(lbpTexCoordLocation);
    gl.vertexAttribPointer(lbpTexCoordLocation, 2, gl.FLOAT, false, 0, 0);

    // set framebuffer to the original one if not already using it
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.gbo);

    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, this.patchCells * 6);

    /* calculate responses */

    gl.useProgram(this.patchResponseProgram);

    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_filters'), 5);
    gl.uniform1i(gl.getUniformLocation(this.patchResponseProgram, 'u_patches'), 3);

    var positionLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_position');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.apositionBuffer);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    // set up texture positions
    var texCoordLocation = gl.getAttribLocation(this.patchResponseProgram, 'a_texCoord');
    gl.bindBuffer(gl.ARRAY_BUFFER, this.texCoordBuffer);
    gl.enableVertexAttribArray(texCoordLocation);
    gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.viewport(0, 0, this.patchWidth, this.patchHeight * this.numBlocks);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    // draw to framebuffer
    gl.drawArrays(gl.TRIANGLES, 0, this.patchCells * 6);

    /* get the responses */

    var responses = this._drawOut('lbp');
    return responses;
  }

  _insertPatches (patches) {
    // pass patches into texture, each patch in either r, g, b or a
    var patchArrayIndex = 0;
    var patchesIndex1 = 0;
    var patchesIndex2 = 0;
    for (var i = 0; i < this.patchCells; i++) {
      for (var j = 0; j < this.patchHeight; j++) {
        for (var k = 0; k < this.patchWidth; k++) {
          patchesIndex1 = i * 4;
          patchesIndex2 = (j * this.patchWidth) + k;
          patchArrayIndex = ((this.patchSize * i) + patchesIndex2) * 4;

          // set r with first patch
          if (patchesIndex1 < this.numPatches) {
            this.patchArray[patchArrayIndex] = patches[patchesIndex1][patchesIndex2];
          } else {
            this.patchArray[patchArrayIndex] = 0;
          }
          // set g with 2nd patch
          if (patchesIndex1 + 1 < this.numPatches) {
            this.patchArray[patchArrayIndex + 1] = patches[patchesIndex1 + 1][patchesIndex2];
          } else {
            this.patchArray[patchArrayIndex + 1] = 0;
          }
          // set b with 3rd patch
          if (patchesIndex1 + 2 < this.numPatches) {
            this.patchArray[patchArrayIndex + 2] = patches[patchesIndex1 + 2][patchesIndex2];
          } else {
            this.patchArray[patchArrayIndex + 2] = 0;
          }
          // set a with 4th patch
          if (patchesIndex1 + 3 < this.numPatches) {
            this.patchArray[patchArrayIndex + 3] = patches[patchesIndex1 + 3][patchesIndex2];
          } else {
            this.patchArray[patchArrayIndex + 3] = 0;
          }
        }
      }
    }

    // pass texture into an uniform
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, this.patchTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, this.textureWidth, this.textureHeight, 0, gl.RGBA, gl.FLOAT, this.patchArray);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  }

  _insertFilter (filter, textureNum) {
    const { filterWidth, filterHeight, numBlocks, gl } = this;

    var filterSize = filterWidth * filterHeight;
    var filterArray = new Float32Array(filterSize * (numBlocks) * 4);
    for (var i = 0; i < numBlocks; i++) {
      for (var j = 0; j < filterHeight; j++) {
        for (var k = 0; k < filterWidth; k++) {
          // set r with first filter
          if (i * 4 < filter.length) {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4] = filter[i * 4][(j * filterWidth) + k];
          } else {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4] = 0;
          }
          // set g with 2nd filter
          if ((i * 4 + 1) < filter.length) {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 1] = filter[(i * 4) + 1][(j * filterWidth) + k];
          } else {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 1] = 0;
          }
          // set b with 3rd filter
          if ((i * 4 + 2) < filter.length) {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 2] = filter[(i * 4) + 2][(j * filterWidth) + k];
          } else {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 2] = 0;
          }
          // set a with 4th filter
          if ((i * 4 + 3) < filter.length) {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 3] = filter[(i * 4) + 3][(j * filterWidth) + k];
          } else {
            filterArray[((filterSize * i) + (j * filterWidth) + k) * 4 + 3] = 0;
          }
        }
      }
    }

    gl.activeTexture(textureNum);
    var filterTexture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, filterTexture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, filterWidth, filterHeight * numBlocks, 0, gl.RGBA, gl.FLOAT, filterArray);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  }

  _drawOut (type) {
    // switch programs
    const gl = this.gl;
    gl.useProgram(this.patchDrawProgram);

    // bind canvas buffer
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, this.newCanvasWidth, this.newCanvasHeight);

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER)

    gl.bindBuffer(gl.ARRAY_BUFFER, this.drawRectBuffer);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      this.drawOutRectangles,
      gl.STATIC_DRAW);
    var positionLocation = gl.getAttribLocation(this.patchDrawProgram, 'a_position_draw');
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, this.drawImageBuffer);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      this.drawOutImages,
      gl.STATIC_DRAW);
    var textureLocation = gl.getAttribLocation(this.patchDrawProgram, 'a_texCoord_draw');
    gl.enableVertexAttribArray(textureLocation);
    gl.vertexAttribPointer(textureLocation, 2, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, this.drawLayerBuffer);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      this.drawOutLayer,
      gl.STATIC_DRAW);
    var layerLocation = gl.getAttribLocation(this.patchDrawProgram, 'a_patchChoice_draw');
    gl.enableVertexAttribArray(layerLocation);
    gl.vertexAttribPointer(layerLocation, 1, gl.FLOAT, false, 0, 0);

    // draw out
    gl.drawArrays(gl.TRIANGLES, 0, this.numPatches * 6);

    let responses = this._getOutput();
    responses = this._unpackToFloat(responses);
    responses = this._splitArray(responses, this.numPatches);
    responses = this._addBias(responses, this.biases[type]);

    // normalize responses to lie within [0,1]
    for (let i = 0, rl = responses.length; i < rl; i++) {
      responses[i] = this._normalizeFilterMatrix(responses[i]);
    }

    return responses;
  }

  _addBias (responses, bias) {
    // do a little trick to add bias in the logit function
    for (let i = 0; i < responses.length; i++) {
      const biasMult = Math.exp(bias[i]);
      for (let j = 0; j < responses[i].length; j++) {
        let response = responses[i][j];
        responses[i][j] = 1 / (1 + ((1 - response) / (response * biasMult)));
      }
    }
    return responses;
  }

  _splitArray (array, parts) {
    var sp = [];
    var al = array.length;
    var splitlength = al / parts;
    var ta = [];
    for (let i = 0; i < al; i++) {
      if (i % splitlength === 0) {
        if (i !== 0) {
          sp.push(ta);
        }
        ta = [];
      }
      ta.push(array[i]);
    }
    sp.push(ta);
    return sp;
  }

  _getOutput () {
    const { gl, canvas } = this;
    // get data
    var pixelValues = new Uint8Array(4 * canvas.width * canvas.height);
    gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, pixelValues);
    return pixelValues;
  }

  _unpackToFloat (array) {
    // convert packed floats to proper floats : see http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work
    const newArray = [];
    const al = array.length;
    for (let i = 0; i < al; i += 4) {
      newArray[(i / 4) >> 0] = (
        array[i] / 4294967296 + // 256 * 256 * 256 * 256
        array[i + 1] / 16777216 + // 256 * 256 * 256
        array[i + 2] / 65536 + // 256 * 256
        array[i + 3] / 256
      );
    }
    return newArray;
  }

  _normalizeFilterMatrix (response) {
    // normalize responses to lie within [0,1]
    const msize = response.length;
    let max = 0;
    let min = 1;

    for (let i = 0; i < msize; i++) {
      let val = response[i];
      if (val > max) {
        max = val;
      }
      if (val < min) {
        min = val;
      }
    }
    const dist = max - min;

    if (dist === 0) {
      console.log('a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.')
      response = response.map(() => 1);
    } else {
      for (let i = 0; i < msize; i++) {
        response[i] = (response[i] - min) / dist;
      }
    }

    return response;
  }
}
