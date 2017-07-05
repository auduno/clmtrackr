// This code is based on webgl-utils.js authored by Gregg Tavares, license below:
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
(function(){

var LOGGING_ENABLED = true;

/**
 * Wrapped logging function.
 * @param {string} msg The message to log.
 */
const log = function (msg) {
  if (!LOGGING_ENABLED) { return; }
  if (window.console && window.console.log) {
    window.console.log(msg);
  }
};

/**
 * Wrapped logging function.
 * @param {string} msg The message to log.
 */
const error = function (msg) {
  if (!LOGGING_ENABLED) { return; }
  if (window.console) {
    if (window.console.error) {
      window.console.error(msg);
    } else if (window.console.log) {
      window.console.log(msg);
    }
  }
  throw msg;
};

/**
 * Turn off all logging.
 */
const loggingOff = function () {
  LOGGING_ENABLED = false;
};

/**
 * Check if the page is embedded.
 * @return {boolean} True of we are in an iframe
 */
const isInIFrame = function () {
  return window !== window.top;
};

/**
 * Converts a WebGL enum to a string
 * @param {!WebGLContext} gl The WebGLContext to use.
 * @param {number} value The enum value.
 * @return {string} The enum as a string.
 */
const glEnumToString = function (gl, value) {
  for (var p in gl) {
    if (gl[p] === value) {
      return p;
    }
  }
  return '0x' + value.toString(16);
};


/**
 * Creates the HTLM for a failure message
 * @param {string} canvasContainerId id of container of th
 *        canvas.
 * @return {string} The html.
 */
const makeFailHTML = function (msg) {
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
// const GET_A_WEBGL_BROWSER = '' +
//   'This page requires a browser that supports WebGL.<br/>' +
//   '<a href="http://get.webgl.org">Click here to upgrade your browser.</a>';


/**
 * Mesasge for need better hardware
 * @type {string}
 */
// const OTHER_PROBLEM = '' +
//   "It doesn't appear your computer can support WebGL.<br/>" +
//   '<a href="http://get.webgl.org/troubleshooting/">Click here for more information.</a>';


/**
 * Creates a webgl context. If creation fails it will
 * change the contents of the container of the <canvas>
 * tag to an error message with the correct links for WebGL.
 * @param {Element} canvas. The canvas element to create a
 *     context from.
 * @param {WebGLContextCreationAttirbutes} optAttribs Any
 *     creation attributes you want to pass in.
 * @return {WebGLRenderingContext} The created context.
 */
const setupWebGL = function (canvas, optAttribs) {
  // const showLink = function (str) {
  //   var container = canvas.parentNode;
  //   if (container) {
  //     container.innerHTML = makeFailHTML(str);
  //   }
  // };

  if (!window.WebGLRenderingContext) {
    // showLink(GET_A_WEBGL_BROWSER);
    return null;
  }

  var context = create3DContext(canvas, optAttribs);
  if (!context) {
    // showLink(OTHER_PROBLEM);
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
const create3DContext = function (canvas, optAttribs) {
  var names = ['webgl', 'experimental-webgl'];
  var context = null;
  for (var ii = 0; ii < names.length; ++ii) {
    try {
      context = canvas.getContext(names[ii], optAttribs);
    } catch (e) {}
    if (context) {
      break;
    }
  }
  return context;
};

const updateCSSIfInIFrame = function () {
  if (isInIFrame()) {
    document.body.className = 'iframe';
  }
};

/**
 * Gets a WebGL context.
 * makes its backing store the size it is displayed.
 */
const getWebGLContext = function (canvas) {
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
 * @param {function(string): void) optErrorCallback callback for errors.
 * @return {!WebGLShader} The created shader.
 */
const loadShader = function (gl, shaderSource, shaderType, optErrorCallback) {
  var errFn = optErrorCallback || error;
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
    var lastError = gl.getShaderInfoLog(shader);
    errFn("*** Error compiling shader '" + shader + "':" + lastError);
    gl.deleteShader(shader);
    return null;
  }

  return shader;
};


/**
 * Creates a program, attaches shaders, binds attrib locations, links the
 * program and calls useProgram.
 * @param {!Array.<!WebGLShader>} shaders The shaders to attach
 * @param {!Array.<string>} optAttribs The attribs names.
 * @param {!Array.<number>} optLocations The locations for the attribs.
 */
const loadProgram = function (gl, shaders, optAttribs, optLocations) {
  var program = gl.createProgram();
  for (var i = 0; i < shaders.length; ++i) {
    gl.attachShader(program, shaders[i]);
  }
  if (optAttribs) {
    for (var i = 0; i < optAttribs.length; ++i) {
      gl.bindAttribLocation(
          program,
          optLocations ? optLocations[i] : i,
          optAttribs[i]);
    }
  }
  gl.linkProgram(program);

  // Check the link status
  const linked = gl.getProgramParameter(program, gl.LINK_STATUS);
  if (!linked) {
    // something went wrong with the link
    const lastError = gl.getProgramInfoLog(program);
    error('Error in program linking:' + lastError);

    gl.deleteProgram(program);
    return null;
  }
  return program;
};


/**
 * Loads a shader from a script tag.
 * @param {!WebGLContext} gl The WebGLContext to use.
 * @param {string} scriptId The id of the script tag.
 * @param {number} optShaderType The type of shader. If not passed in it will
 *     be derived from the type of the script tag.
 * @param {function(string): void) optErrorCallback callback for errors.
 * @return {!WebGLShader} The created shader.
 */
const createShaderFromScript = function (
  gl, scriptId, optShaderType, optErrorCallback
) {
  var shaderSource = '';
  var shaderType;
  var shaderScript = document.getElementById(scriptId);
  if (!shaderScript) {
    throw new Error('*** Error: unknown script element' + scriptId);
  }
  shaderSource = shaderScript.text;

  if (!optShaderType) {
    if (shaderScript.type === 'x-shader/x-vertex') {
      shaderType = gl.VERTEX_SHADER;
    } else if (shaderScript.type === 'x-shader/x-fragment') {
      shaderType = gl.FRAGMENT_SHADER;
    } else if (
      shaderType !== gl.VERTEX_SHADER &&
      shaderType !== gl.FRAGMENT_SHADER
    ) {
      throw new Error('*** Error: unknown shader type');
    }
  }

  return loadShader(
    gl,
    shaderSource,
    optShaderType || shaderType,
    optErrorCallback
  );
};

if (typeof exports === 'object' && typeof module !== 'undefined') {
  module.exports = {
    setupWebGL : setupWebGL,
    createProgram : loadProgram,
    createShaderFromScript : createShaderFromScript,
    getWebGLContext : getWebGLContext,
    loadShader : loadShader
  }
} else {
  // export to global
  window.setupWebGL = setupWebGL;
  window.createProgram = loadProgram;
  window.createShaderFromScriptElement = createShaderFromScript;
  window.getWebGLContext = getWebGLContext;
  window.loadShader = loadShader;
}

}());