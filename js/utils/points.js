export const getBoundingBox = (points) => {
  let maxX = 0;
  let minX = 0;
  let maxY = 0;
  let minY = 0;

  for (let i = 0; i < points.length; i++) {
    const x = points[i][0];
    if (x > maxX) maxX = x;
    else if (x < minX) minX = x;

    const y = points[i][1];
    if (y > maxY) maxY = y;
    else if (y < minY) minY = y;
  }

  minX = Math.floor(minX);
  maxX = Math.ceil(maxX);
  minY = Math.floor(minY);
  maxY = Math.ceil(maxY);

  return {
    maxX,
    minX,
    maxY,
    minY,
    width: maxX - minX,
    height: maxY - minY
  };
};


export const generateTextureVertices = (points, vertMap, scaleX = 1, scaleY = 1) => {
  const tvs = new Float32Array(vertMap.length * 6);
  for (let i = 0; i < vertMap.length; i++) {
    const j = i * 6;
    tvs[j] = points[vertMap[i][0]][0] * scaleX;
    tvs[j + 1] = points[vertMap[i][0]][1] * scaleY;
    tvs[j + 2] = points[vertMap[i][1]][0] * scaleX;
    tvs[j + 3] = points[vertMap[i][1]][1] * scaleY;
    tvs[j + 4] = points[vertMap[i][2]][0] * scaleX;
    tvs[j + 5] = points[vertMap[i][2]][1] * scaleY;
  }
  return tvs;
};
