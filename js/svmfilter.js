// requires webgl-utils.js from html5rocks tutorial:
//  https://github.com/html5rocks/www.html5rocks.com/blob/master/content/tutorials/webgl/webgl_fundamentals/static/webgl/resources/webgl-utils.js

var webglFilter = function() {
  
  var gl, canvas;
  var filterWidth, filterHeight, patchWidth, patchHeight, numPatches, canvasWidth, canvasHeight;
  var corrFilterWidth, corrFilterHeight;
  var patchResponseProgram, patchDrawProgram;
  var initiated;
  var fbo, renderBuffer, numBlocks, patchTex;
  var first = true;
  var drawRectBuffer, drawLayerBuffer, drawImageBuffer, rttTexture, filters;
  var texCoordBuffer, texCoordLocation, apositionBuffer;
  
  var newCanvasWidth, newCanvasBlockHeight, newCanvasHeight;
  var drawOutRectangles, drawOutImages, drawOutLayer;
  
  this.init = function(filterVector, nP, pW, pH, fW, fH, drawOut) {
    // we assume filterVector goes from left to right, rowwise
    if (fW != fH) {
      alert("filter width and height must be same size!");
      return;
    }
    
    // if filter width is not odd, we pad it with 0s on bottom and right side
    if (fW % 2 == 0 || fH % 2 == 0) {
      var newFilter = [];
      filterWidth = fW + 1;
      filterHeight = fH + 1;
      var nfsize = filterWidth * filterHeight;
      
      for (var i = nP-1;i > -1;i--) {
        for (var j = 0;j < fW+1;j++) {
          filterVector.splice(((i+1)*fW*fH), 0, 0)
        }
        for (var j = fH;j > 0;j--) {
          filterVector.splice((i*fW*fH)+(j*fW), 0, 0)
        }
      }
      corrFilterWidth = filterWidth-1;
      corrFilterHeight = filterHeight-1;
    } else {
      filterWidth = fW;
      filterHeight = fH;
      corrFilterWidth = filterWidth;
      corrFilterHeight = filterHeight;
    }

    patchWidth = pW;
    patchHeight = pH;
    numPatches = nP;
    
    //testing of printing patches
    
    /*var canvasppt = document.createElement('canvas')
    canvasppt.setAttribute('width', filterWidth+"px");
    canvasppt.setAttribute('height', (filterHeight*numPatches)+"px");
    canvasppt.setAttribute('id', 'patchprinttest2');
    document.body.appendChild(canvasppt);
    var pptcc = canvasppt.getContext('2d');
    for (var i = 0; i < numPatches; i++) {
      var psci = pptcc.createImageData(filterWidth, filterHeight);
      var pscidata = psci.data;
      for (var j = 0;j < (filterWidth)*(filterHeight);j++) {
        var val = filterVector[(i*filterWidth*filterHeight)+j]*2000+127;
        pscidata[j*4] = val;
        pscidata[(j*4)+1] = val;
        pscidata[(j*4)+2] = val;
        pscidata[(j*4)+3] = 255;
      }
      pptcc.putImageData(psci, 0, filterHeight*i);
    }*/
    
    patchResponseFS = [
      "precision mediump float;",
      "",
      "uniform vec2 u_onePixelFilters;",
      "uniform vec2 u_onePixelPatches;",
      "const float u_filterwidth = "+corrFilterWidth.toFixed(1)+";",
      "const float u_filterheight = "+corrFilterHeight.toFixed(1)+";",
      "const float u_halffilterwidth = "+(corrFilterWidth/2).toFixed(1)+";",
      "const float u_halffilterheight = "+(corrFilterHeight/2).toFixed(1)+";",
      "const float u_filtersize = "+(corrFilterWidth*corrFilterHeight)+".0;",
      "const float u_divfiltersize = "+(1/(corrFilterWidth*corrFilterHeight)).toFixed(10)+";",
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
      "  // normalize the selection we take first",
      "  vec4 msum = vec4(0.0, 0.0, 0.0, 0.0);",
      "  vec4 msumsq = vec4(0.0, 0.0, 0.0, 0.0);",
      "  vec4 mean = vec4(0.0, 0.0, 0.0, 0.0);",
      "  for (float w = 0.0;w < u_filterwidth;w++) {",
      "    for (float h = 0.0;h < u_filterheight;h++) {",
      "      mean = texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(w-u_halffilterwidth, h-u_halffilterheight));",
      "      msum += mean;",
      "      msumsq += (mean * mean);",
      "    } ",
      "  }",
      "  mean = msum * u_divfiltersize;",
      "  vec4 sd = sqrt((msumsq * u_divfiltersize) - (mean * mean));",
      "  ",
      "  for (float w = 0.0;w < u_filterwidth;w++) {",
      "    for (float h = 0.0;h < u_filterheight;h++) {",
      "      colorSum += ",
      "        ((texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(w-u_halffilterwidth, h-u_halffilterheight))-mean)/(sd)) * ",
      "        texture2D(u_filters, v_texCoordFilters + u_onePixelFilters * vec2(w-u_halffilterwidth, h-u_halffilterheight)); ",
      "    } ",
      "  }",
      "  // logit",
      "  colorSum = 1.0-(1.0/(1.0 + exp(colorSum)));",
      "  gl_FragColor = colorSum;",
      "}"
    ].join('\n');
    
    numBlocks = Math.floor(numPatches / 4) + Math.ceil(numPatches % 4);
    canvasWidth = patchWidth;
    canvasHeight = patchHeight*numBlocks;
    
    //create canvas
    canvas = document.createElement('canvas')
    canvas.setAttribute('width', (patchWidth-corrFilterWidth)+"px");
    //canvas.setAttribute('width', (patchWidth-(filterWidth))+"px");
    canvas.setAttribute('height', ((patchHeight-corrFilterHeight)*numPatches)+"px");
    //canvas.setAttribute('height', ((patchHeight-(filterHeight))*numPatches)+"px");
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
    var halfFilter = corrFilterWidth/2;
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
    
    // setup patchresponse program
    var prVertexShader = loadShader(gl, patchResponseVS, gl.VERTEX_SHADER);
    var prFragmentShader = loadShader(gl, patchResponseFS, gl.FRAGMENT_SHADER);
    patchResponseProgram = createProgram(gl, [prVertexShader, prFragmentShader]);
    gl.useProgram(patchResponseProgram);
    
    // setup patchdraw program
    var drVertexShader = loadShader(gl, drawResponsesVS, gl.VERTEX_SHADER);
    var drFragmentShader = loadShader(gl, drawResponsesFS, gl.FRAGMENT_SHADER);
    patchDrawProgram = createProgram(gl, [drVertexShader, drFragmentShader]);
    
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
    //gl.uniform2f(onePixelPatchLocation, 1/filterWidth, 1/(filterHeight*numBlocks));
    gl.uniform2f(onePixelPatchLocation, 1/patchWidth, 1/(patchHeight*numBlocks));
    
    // set the onepixel size for filters
    var onePixelFilterLocation = gl.getUniformLocation(patchResponseProgram, "u_onePixelFilters");
    //gl.uniform2f(onePixelFilterLocation, 1/canvasWidth, 1/canvasHeight);
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
    var filterArray = new Float32Array(filterSize*(numBlocks+1)*4);
    for (var i = 0;i < numBlocks;i++) {
      for (var j = 0;j < filterHeight;j++) {
        for (var k = 0;k < filterWidth;k++) {
          //set r with first filter
          if (filterSize*i*4 < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = filterVector[(filterSize*i*4) + (k*filterWidth) + j];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4] = 0;
          }
          //set g with 2nd filter
          if ((filterSize*i*4 + 1) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = filterVector[(filterSize*(i*4 + 1)) + (k*filterWidth) + j];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 1] = 0;
          }
          //set b with 3rd filter
          if ((filterSize*i*4 + 2) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = filterVector[(filterSize*(i*4 + 2)) + (k*filterWidth) + j];
          } else {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 2] = 0;
          }
          //set a with 4th filter
          if ((filterSize*i*4 + 3) < filterVector.length) {
            filterArray[((filterSize*i) + (j*filterWidth) + k)*4 + 3] = filterVector[(filterSize*(i*4 + 3)) + (k*filterWidth) + j];
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
    //renderBuffer = gl.createRenderbuffer();
    //gl.bindRenderbuffer(gl.RENDERBUFFER, renderBuffer);
    
    //gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, patchWidth, patchHeight*numBlocks);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, rttTexture, 0);
    //gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, renderBuffer);
    
    gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);
    
    // preinitialization of drawOut variables
    
    newCanvasWidth = patchWidth-corrFilterWidth;
    newCanvasBlockHeight = patchHeight-corrFilterWidth;
    newCanvasHeight = newCanvasBlockHeight*numPatches;
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
    var patchCells = (Math.floor(numPatches / 4) + Math.ceil(numPatches % 4));
    var halfFilterWidth = (corrFilterWidth/2)/patchWidth;
    var halfFilterHeight = (corrFilterWidth/2)/(patchHeight*patchCells);
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
  }

  this.getResponses = function(patches, drawOut) {
    // TODO: check patches correct length/dimension
    
    // switch to response generation program if we're not already using it
    if (!first) {
      gl.useProgram(patchResponseProgram);
      
      /*
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
      */
      
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

      /*gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, filters);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.uniform1i(gl.getUniformLocation(patchResponseProgram, "u_filters"), 0);

      gl.activeTexture(gl.TEXTURE2);
      gl.bindTexture(gl.TEXTURE_2D, rttTexture);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);*/
      
      // set framebuffer to the original one if not already using it
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      //gl.bindRenderbuffer(gl.RENDERBUFFER, renderBuffer);
      
      gl.viewport(0, 0, patchWidth, patchHeight*numBlocks);
    }
    
    // pass patches into texture, each patch in either r, g, b or a
    var numPatches = patches.length;
    var patchCells = (Math.floor(numPatches / 4) + Math.ceil(numPatches % 4));
    var textureWidth = patchWidth;
    var textureHeight = patchHeight*patchCells;
    var patchSize = patchWidth*patchHeight;
    var patchArray = new Float32Array(patchSize*(patchCells+1)*4);
    
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
    
    gl.finish();
    
    if (drawOut) {
      
      // switch programs
      gl.useProgram(patchDrawProgram);
      
      // bind canvas buffer
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, newCanvasWidth, newCanvasHeight);
      
      gl.clearColor(0.0, 0.0, 0.0, 1.0);
      gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER)
      
      // set the resolution/dimension of the canvas
      var resolutionLocation = gl.getUniformLocation(patchDrawProgram, "u_resolutiondraw");
      gl.uniform2f(resolutionLocation, newCanvasWidth, newCanvasHeight);
      
      // set u_responses
      var responsesLocation = gl.getUniformLocation(patchDrawProgram, "u_responses");
      gl.uniform1i(responsesLocation, 2);
      
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
    var starttime = (new Date).getTime();
    var data = gl.readPixels(0, 0, canvas.width, canvas.height, gl.RGBA, gl.UNSIGNED_BYTE, pixelValues);
    //var data = gl.readPixels(0, 0, 2, 2, gl.RGBA, gl.UNSIGNED_BYTE, pixelValues);
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
    "  if (v_select == 0.0) {",
    "    value = colorSum[0];",
    "  } else if (v_select == 1.0) {",
    "    value = colorSum[1];",
    "  } else if (v_select == 2.0) {",
    "    value = colorSum[2];",
    "  } else if (v_select == 3.0) {",
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