
var drawDot = function(ctx, coordinates) {
	ctx.beginPath();
	ctx.arc(coordinates[0], coordinates[1], 3, 0, Math.PI*2, true);
	ctx.closePath();
	ctx.fill();
}

export const drawBoundingBox = function(ctx, box) {
	ctx.strokeRect( box.x, box.y, box.width, box.height );
}

export const drawDetection = function(ctx, bbox, facePoints) {
	var noseFilterWidth = bbox.width * 4.5/10;
	var eyeFilterWidth = bbox.width * 6/10;

	ctx.clearRect(0, 0, ctx.canvas.clientWidth, ctx.canvas.clientHeight);
	drawBoundingBox(ctx, bbox);

	drawBoundingBox(ctx, [
		Math.round( bbox.x + (bbox.width*3/4)-(eyeFilterWidth/2)),
		Math.round(bbox.y+bbox.height*(2/5)-(eyeFilterWidth/2)),
		eyeFilterWidth,
		eyeFilterWidth
	]);
	drawBoundingBox(ctx, [
		Math.round( bbox.x + (bbox.width/4)-(eyeFilterWidth/2)),
		Math.round(bbox.y+bbox.height*(2/5)-(eyeFilterWidth/2)),
		eyeFilterWidth,
		eyeFilterWidth
	]);
	drawBoundingBox(ctx, [
		Math.round( bbox.x + (bbox.width/2)-(noseFilterWidth/2)),
		Math.round(bbox.y+bbox.height*(5/8)-(noseFilterWidth/2)),
		noseFilterWidth,
		noseFilterWidth
	]);

	drawFacialPoints(ctx, facePoints);
}

export const drawFacialPoints = function(ctx, facePoints, transformParams) {

	// transform facepoints
	if (transformParams) {
		var translateX = transformParams[0];
		var translateY = transformParams[1];
		var scaling = transformParams[2];
		var rotation = transformParams[3];

		facePoints = facePoints.map(function(coord) {
			return [
				(coord[0]*scaling*Math.cos(-rotation) + coord[1]*scaling*Math.sin(-rotation)) + translateX,
				(coord[0]*scaling*(-Math.sin(-rotation)) + coord[1]*scaling*Math.cos(-rotation)) + translateY
			];
		});

		ctx.fillStyle = 'rgb(200,10,100)';
	}

	var leftEye = facePoints[0];
	var rightEye = facePoints[1];
	var nose = facePoints[2];

	drawDot(ctx, leftEye);
	drawDot(ctx, rightEye);
	drawDot(ctx, nose);
}

// function to draw pixeldata on some canvas, only used for debugging
export const drawData = function(canvasContext, data, width, height, transposed, drawX, drawY) {
	const psci = canvasContext.createImageData(width, height);
	const pscidata = psci.data;
	for (var j = 0; j < width * height; j++) {
		var val;
		if (!transposed) {
			val = data[(j % width) + ((j / width) >> 0) * width];
		} else {
			val = data[(j % height) * height + ((j / height) >> 0)];
		}
		val = val > 255 ? 255 : val;
		val = val < 0 ? 0 : val;
		pscidata[j * 4] = val;
		pscidata[(j * 4) + 1] = val;
		pscidata[(j * 4) + 2] = val;
		pscidata[(j * 4) + 3] = 255;
	}
	canvasContext.putImageData(psci, drawX, drawY);
}

export const drawPatches = function(ctx, patches, patchSize, patchPositions, mapFunction, subset) {
	var halfPatch = Math.floor(patchSize/2);
	ctx.clearRect(0, 0, ctx.canvas.clientWidth, ctx.canvas.clientHeight);
	for (var i = 0;i < patches.length;i++) {
		if (!subset || (subset && subset.indexOf(i) > -1)) {
			var patch = mapFunction ? patches[i].map(mapFunction) : patches[i];
			drawData(ctx, patch, patchSize, patchSize, false, patchPositions[i][0]-halfPatch, patchPositions[i][1]-halfPatch);
		}
	}
}
