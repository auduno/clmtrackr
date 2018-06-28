var paramsRange = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

      function rangeInput(val, num) {
        // get rangeinput
        paramsRange[num] = val;
        drawNew2(paramsRange);
      }

      var drawNew2 = function(params) {
        cc.clearRect(0,0,600,600);
        draw(document.getElementById('compare'), similarityTransforms.concat(params));
      }

			var canvasInput = document.getElementById('compare');
			var cc = canvasInput.getContext('2d');
			cc.fillStyle = "rgb(200,0,0)";

      var similarityTransforms = [2,0,0,-50];
      var paramslength = paramsRange.length;
      var num_points = pModel.shapeModel.numPtsPerSample;

			var x,y;

			var i, path;

			var drawPath = function(canvasContext, path, dp) {
			canvasContext.beginPath();
			var i, x, y;
			for (var p = 0;p < path.length;p++) {
				i = path[p]*2;
				x = pModel.shapeModel.meanShape[i/2][0];
				y = pModel.shapeModel.meanShape[i/2][1];
				for (var j = 0;j < pModel.shapeModel.numEvalues;j++) {
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
			for (var j = 0;j < pModel.shapeModel.numEvalues;j++) {
				x += pModel.shapeModel.eigenVectors[i][j]*dp[j+4];
				y += pModel.shapeModel.eigenVectors[i+1][j]*dp[j+4];
			}
			a = dp[0]*x - dp[1]*y + dp[2];
			b = dp[0]*y + dp[1]*x + dp[3];
			x += a;
			y += b;
			canvasContext.beginPath();
		  canvasContext.arc(x, y, 2, 0, Math.PI*2, true);
			canvasContext.closePath();
			canvasContext.fill();
		}


		var draw = function(canvas, pv) {
      // if no previous points, just draw in the middle of canvas

      var params;
      if (pv === undefined) {
        params = currentParameters.slice(0);
      } else {
        params = pv.slice(0);
      }

      var cc = canvas.getContext('2d');
      cc.fillStyle = "rgb(50,50,50)";
      cc.strokeStyle = "rgb(50,50,50)";
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
	  draw(document.getElementById('compare'), similarityTransforms.concat(paramsRange))
			