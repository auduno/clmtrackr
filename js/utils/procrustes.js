const procrustes = (template, shape) => {
  // assume template and shape is a vector of x,y-coordinates
  //i.e. template = [[x1,y1], [x2,y2], [x3,y3]];
  var templateClone = [];
  var shapeClone = [];
  for (var i = 0;i < template.length;i++) {
    templateClone[i] = [template[i][0], template[i][1]];
  }
  for (var i = 0;i < shape.length;i++) {
    shapeClone[i] = [shape[i][0], shape[i][1]];
  }
  shape = shapeClone;
  template = templateClone;

  // calculate translation
  var templateMean = [0.0, 0.0];
  for (var i = 0;i < template.length;i++) {
    templateMean[0] += template[i][0];
    templateMean[1] += template[i][1];
  }
  templateMean[0] /= template.length;
  templateMean[1] /= template.length;

  var shapeMean = [0.0, 0.0];
  for (var i = 0;i < shape.length;i++) {
    shapeMean[0] += shape[i][0];
    shapeMean[1] += shape[i][1];
  }
  shapeMean[0] /= shape.length;
  shapeMean[1] /= shape.length;

  var translationX = templateMean[0] - shapeMean[0];
  var translationY = templateMean[1] - shapeMean[1];

  // centralize
  for (var i = 0;i < shape.length;i++) {
    shape[i][0] -= shapeMean[0];
    shape[i][1] -= shapeMean[1];
  }
  for (var i = 0;i < template.length;i++) {
    template[i][0] -= templateMean[0];
    template[i][1] -= templateMean[1];
  }

  // scaling

  var scaleS = 0.0;
  for (var i = 0;i < shape.length;i++) {
    scaleS += ((shape[i][0])*(shape[i][0]));
    scaleS += ((shape[i][1])*(shape[i][1]));
  }
  scaleS = Math.sqrt(scaleS/shape.length);

  var scaleT = 0.0;
  for (var i = 0;i < template.length;i++) {
    scaleT += ((template[i][0])*(template[i][0]));
    scaleT += ((template[i][1])*(template[i][1]));
  }
  scaleT = Math.sqrt(scaleT/template.length);

  var scaling = scaleT/scaleS;

  for (var i = 0;i < shape.length;i++) {
    shape[i][0] *= scaling;
    shape[i][1] *= scaling;
  }

  // rotation

  var top = 0.0;
  var bottom = 0.0;
  for (var i = 0;i < shape.length;i++) {
    top += (shape[i][0]*template[i][1] - shape[i][1]*template[i][0]);
    bottom += (shape[i][0]*template[i][0] + shape[i][1]*template[i][1]);
  }
  var rotation = Math.atan(top/bottom);

  translationX += (shapeMean[0]-(scaling*Math.cos(-rotation)*shapeMean[0])-(scaling*shapeMean[1]*Math.sin(-rotation)));
  translationY += (shapeMean[1]+(scaling*Math.sin(-rotation)*shapeMean[0])-(scaling*shapeMean[1]*Math.cos(-rotation)));

  //returns rotation, scaling, transformx and transformx
  return [translationX, translationY, scaling, rotation];
}

export default procrustes;
