var faceDeformer = function() {
  
  var gl, verticeMap, originalPositions;
  var numTriangles;
  var maxx, minx, maxy, miny;
  var first = true;
  
  // set vertice coordinates to draw on a_position
  // set original image vertices on a_texCoord
  
  this.init = function(canvas) {
    // ready a webgl element
    gl = getWebGLContext(canvas); 
  }

  this.load = function(element, sketchCanvasContext, points, triangles, initw, inith) {
    verticeMap = triangles;
    numTriangles = triangles.length;
    
    // get cropping
    maxx = 0;
    minx = element.width;
    maxy = 0;
    miny = element.height;
    for (var i = 0;i < points.length;i++) {
      if (points[i][0] > maxx) maxx = points[i][0];
      if (points[i][0] < minx) minx = points[i][0];
      if (points[i][1] > maxy) maxy = points[i][1];
      if (points[i][1] < miny) miny = points[i][1];
    }
    // get image
    if (element.tagName == 'VIDEO' || element.tagName == 'IMG') {
      var ca = document.createElement('canvas');
      ca.width = element.width;
      ca.height = element.height;
      var cc = ca.getContext('2d');
      cc.drawImage(element, 0, 0, element.width, element.height);
      //sketchCanvasContext.drawImage(element, 0, 0, element.width, element.height);
    } else if (element.tagName == 'CANVAS') {
      var cc = element.getContext('2d');
    }
    var image = cc.getImageData(minx, miny, maxx-minx+1, maxy-miny+1);
    //var image = sketchCanvasContext.getImageData(Math.round(minx), Math.round(miny), Math.round(maxx-minx+1), Math.round(maxy-miny+1));
    
    // correct points
    var nupoints = [];
    for (var i = 0;i < points.length;i++) {
      nupoints[i] = [];
      nupoints[i][0] = points[i][0] - minx;
      nupoints[i][1] = points[i][1] - miny;
    }
    
    originalPositions = points;
    
    // create vertices based on points
    var vertices = [];
    for (var i = 0;i < verticeMap.length;i++) {
      vertices.push(nupoints[verticeMap[i][0]][0]/(maxx-minx));
      vertices.push(nupoints[verticeMap[i][0]][1]/(maxy-miny));
      vertices.push(nupoints[verticeMap[i][1]][0]/(maxx-minx));
      vertices.push(nupoints[verticeMap[i][1]][1]/(maxy-miny));
      vertices.push(nupoints[verticeMap[i][2]][0]/(maxx-minx));
      vertices.push(nupoints[verticeMap[i][2]][1]/(maxy-miny)); 
    }
    
    if (first) {
      var vertexShaderProg = [
        "attribute vec2 a_texCoord;",
        "attribute vec2 a_position;",
        "",
        "varying vec2 v_texCoord;",
        "",
        "uniform vec2 u_resolution;",
        "",
        "void main() {",
        "  vec2 zeroToOne = a_position / u_resolution;",
        "  vec2 zeroToTwo = zeroToOne * 2.0;",
        "  vec2 clipSpace = zeroToTwo - 1.0;",
        "  gl_Position = vec4(clipSpace * vec2(1, -1), 0, 1);",
        "  ",
        "  v_texCoord = a_texCoord;",
        "}"
      ].join('\n');
      
      var fragmentShaderProg = [
        "precision mediump float;",
        "",
        "uniform sampler2D u_image;",
        "",
        "varying vec2 v_texCoord;",
        "",
        "void main() {",
        "  gl_FragColor = texture2D(u_image, v_texCoord);",
        "}"
      ].join('\n');
    
      var vertexShader = loadShader(gl, vertexShaderProg, gl.VERTEX_SHADER);
      var fragmentShader = loadShader(gl, fragmentShaderProg, gl.FRAGMENT_SHADER);
      program = createProgram(gl, [vertexShader, fragmentShader]);
      gl.useProgram(program);
      first = false;
    }
    
    // look up where the vertex data needs to go.
    var texCoordLocation = gl.getAttribLocation(program, "a_texCoord");
    
    // set coordinates to get face texture from
    var texCoordBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(texCoordLocation);
    gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);
    
    // Create a texture.
    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    
    // Set the parameters so we can render any size image.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    
    // Upload the face image into the texture.
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
    
    // set the resolution
    var resolutionLocation = gl.getUniformLocation(program, "u_resolution");
    gl.uniform2f(resolutionLocation, initw, inith);
  }

  this.draw = function(points) {
    // correct points
    var nupoints = [];
    for (var i = 0;i < points.length;i++) {
      nupoints[i] = [];
      nupoints[i][0] = points[i][0]*1.0;
      nupoints[i][1] = points[i][1]*1.0;
      //nupoints[i][0] = (points[i][0] - minx)*0.9;
      //nupoints[i][1] = (points[i][1] - miny)*0.9;
    }
    
    // create drawvertices based on points
    var vertices = [];
    for (var i = 0;i < verticeMap.length;i++) {
      vertices.push(nupoints[verticeMap[i][0]][0]);
      vertices.push(nupoints[verticeMap[i][0]][1]);
      vertices.push(nupoints[verticeMap[i][1]][0]);
      vertices.push(nupoints[verticeMap[i][1]][1]);
      vertices.push(nupoints[verticeMap[i][2]][0]);
      vertices.push(nupoints[verticeMap[i][2]][1]); 
    }
    
    var positionLocation = gl.getAttribLocation(program, "a_position");
    
    // Create a buffer for the position of the rectangle corners.
    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
    
    // Draw the rectangle.
    gl.drawArrays(gl.TRIANGLES, 0, numTriangles*3);
  }
  
  this.calculatePositions = function(parameters, useTransforms) {
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
}
