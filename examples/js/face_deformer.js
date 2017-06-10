"use strict";
var faceDeformer = function() {

	var gl, verticeMap;
	var numTriangles;
	var maxx, minx, maxy, miny;
	var width, height;
	var first = true;
	var texCoordBuffer, gridCoordbuffer;
	var texCoordLocation;
	var pdmModel;

	var usegrid = false;
	var drawProgram, gridProgram;

	this.init = function(canvas) {
		// ready a webgl element
		gl = getWebGLContext(canvas);
		gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, true);

	}

	this.load = function(element, points, pModel, vertices) {
		pdmModel = pModel;
		if (vertices) {
			verticeMap = vertices;
		} else {
			verticeMap = pdmModel.path.vertices;
		}
		numTriangles = verticeMap.length;

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
		minx = Math.floor(minx);
		maxx = Math.ceil(maxx);
		miny = Math.floor(miny);
		maxy = Math.ceil(maxy);
		width = maxx-minx;
		height = maxy-miny;

		if (element.tagName == 'VIDEO' || element.tagName == 'IMG') {
			var ca = document.createElement('canvas');
			ca.width = element.width;
			ca.height = element.height;
			var cc = ca.getContext('2d');
			cc.drawImage(element, 0, 0, element.width, element.height);
		} else if (element.tagName == 'CANVAS') {
			var cc = element.getContext('2d');
		}
		var image = cc.getImageData(minx, miny, width, height);

		// correct points
		var nupoints = [];
		for (var i = 0;i < points.length;i++) {
			nupoints[i] = [];
			nupoints[i][0] = points[i][0] - minx;
			nupoints[i][1] = points[i][1] - miny;
		}

		// create vertices based on points
		var textureVertices = [];
		for (var i = 0;i < verticeMap.length;i++) {
			textureVertices.push(nupoints[verticeMap[i][0]][0]/width);
			textureVertices.push(nupoints[verticeMap[i][0]][1]/height);
			textureVertices.push(nupoints[verticeMap[i][1]][0]/width);
			textureVertices.push(nupoints[verticeMap[i][1]][1]/height);
			textureVertices.push(nupoints[verticeMap[i][2]][0]/width);
			textureVertices.push(nupoints[verticeMap[i][2]][1]/height);
		}

		if (first) {
			// create program for drawing grid
			var gridVertexShaderProg = [
				"attribute vec2 a_position;",
				"",
				"uniform vec2 u_resolution;",
				"",
				"void main() {",
				"  vec2 zeroToOne = a_position / u_resolution;",
				"  vec2 zeroToTwo = zeroToOne * 2.0;",
				"  vec2 clipSpace = zeroToTwo - 1.0;",
				"  gl_Position = vec4(clipSpace * vec2(1, -1), 0, 1);",
				"}"
			].join('\n');

			var gridFragmentShaderProg = [
				"void main() {",
				"  gl_FragColor = vec4(0.2, 0.2, 0.2, 1.0);",
				"}"
			].join('\n');

			var gridVertexShader = loadShader(gl, gridVertexShaderProg, gl.VERTEX_SHADER);
			var gridFragmentShader = loadShader(gl, gridFragmentShaderProg, gl.FRAGMENT_SHADER);
			try {
				gridProgram = createProgram(gl, [gridVertexShader, gridFragmentShader]);
			} catch(err) {
				alert("There was a problem setting up the webGL programs. Maybe you should try it in another browser. :(");
			}

			gridCoordbuffer = gl.createBuffer();

			// create program for drawing deformed face
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
			drawProgram = createProgram(gl, [vertexShader, fragmentShader]);

			texCoordBuffer = gl.createBuffer();

			first = false;
		}

		// load program for drawing grid
		gl.useProgram(gridProgram);

		// set the resolution for grid program
		var resolutionLocation = gl.getUniformLocation(gridProgram, "u_resolution");
		gl.uniform2f(resolutionLocation, gl.drawingBufferWidth, gl.drawingBufferHeight);

		// load program for drawing deformed face
		gl.useProgram(drawProgram);

		// look up where the vertex data needs to go.
		texCoordLocation = gl.getAttribLocation(drawProgram, "a_texCoord");

		// provide texture coordinates for face vertices (i.e. where we're going to copy face vertices from).
		gl.enableVertexAttribArray(texCoordLocation);

		gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(textureVertices), gl.STATIC_DRAW);

		gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

		// Create the texture.
		var texture = gl.createTexture();
		gl.bindTexture(gl.TEXTURE_2D, texture);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

		// Upload the image into the texture.
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);

		// set the resolution for draw program
		resolutionLocation = gl.getUniformLocation(drawProgram, "u_resolution");
		gl.uniform2f(resolutionLocation, gl.drawingBufferWidth, gl.drawingBufferHeight);
	}

	this.draw = function(points) {

		if (usegrid) {
			// switch program if needed
			gl.useProgram(drawProgram);

			//texCoordLocation = gl.getAttribLocation(drawProgram, "a_texCoord");

			gl.enableVertexAttribArray(texCoordLocation);
			gl.bindBuffer(gl.ARRAY_BUFFER, texCoordBuffer);
			gl.vertexAttribPointer(texCoordLocation, 2, gl.FLOAT, false, 0, 0);

			usegrid = false;
		}

		// create drawvertices based on points
		var vertices = [];
		for (var i = 0;i < verticeMap.length;i++) {
			vertices.push(points[verticeMap[i][0]][0]);
			vertices.push(points[verticeMap[i][0]][1]);
			vertices.push(points[verticeMap[i][1]][0]);
			vertices.push(points[verticeMap[i][1]][1]);
			vertices.push(points[verticeMap[i][2]][0]);
			vertices.push(points[verticeMap[i][2]][1]);
		}

		var positionLocation = gl.getAttribLocation(drawProgram, "a_position");

		// Create a buffer for the position of the vertices.
		var drawPosBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ARRAY_BUFFER, drawPosBuffer);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);

		// Draw the face vertices
		gl.drawArrays(gl.TRIANGLES, 0, numTriangles*3);
	}

	this.drawGrid = function(points) {

		if (!usegrid) {
			gl.useProgram(gridProgram);

			usegrid = true;
		}

		// create drawvertices based on points
		var vertices = [];
		// create new texturegrid
		for (var i = 0;i < verticeMap.length;i++) {
			vertices.push(points[verticeMap[i][0]][0]);
			vertices.push(points[verticeMap[i][0]][1]);
			vertices.push(points[verticeMap[i][1]][0]);
			vertices.push(points[verticeMap[i][1]][1]);

			vertices.push(points[verticeMap[i][1]][0]);
			vertices.push(points[verticeMap[i][1]][1]);
			vertices.push(points[verticeMap[i][2]][0]);
			vertices.push(points[verticeMap[i][2]][1]);

			vertices.push(points[verticeMap[i][2]][0]);
			vertices.push(points[verticeMap[i][2]][1]);
			vertices.push(points[verticeMap[i][0]][0]);
			vertices.push(points[verticeMap[i][0]][1]);
		}

		var positionLocation = gl.getAttribLocation(gridProgram, "a_position");

		// Create a buffer for position of the vertices (lines)
		gl.bindBuffer(gl.ARRAY_BUFFER, gridCoordbuffer);
		gl.enableVertexAttribArray(positionLocation);
		gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);

		// Draw the lines
		gl.drawArrays(gl.LINES, 0, numTriangles*6);
	}

	this.clear = function() {
		gl.clear(gl.COLOR_BUFFER_BIT);
	}

	this.calculatePositions = function(parameters, useTransforms) {
		var x, y, a, b;
		var numParameters = parameters.length;
		var positions = [];
		for (var i = 0;i < pdmModel.patchModel.numPatches;i++) {
			x = pdmModel.shapeModel.meanShape[i][0];
			y = pdmModel.shapeModel.meanShape[i][1];
			for (var j = 0;j < numParameters-4;j++) {
				x += pdmModel.shapeModel.eigenVectors[(i*2)][j]*parameters[j+4];
				y += pdmModel.shapeModel.eigenVectors[(i*2)+1][j]*parameters[j+4];
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
