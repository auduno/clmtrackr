// part one of meanshift calculation
export const gpopt = (
  responseWidth,
  currentPositionsj,
  updatePosition,
  vecProbs,
  responses,
  opj0,
  opj1,
  j,
  variance,
  scaling
) => {
  let posIdx = 0;
  let vpsum = 0;
  let dx;
  let dy;
  for (let k = 0; k < responseWidth; k++) {
    updatePosition[1] = opj1 + (k * scaling);
    for (let l = 0; l < responseWidth; l++) {
      updatePosition[0] = opj0 + (l * scaling);

      dx = currentPositionsj[0] - updatePosition[0];
      dy = currentPositionsj[1] - updatePosition[1];
      vecProbs[posIdx] = responses[j][posIdx] * Math.exp(
        -0.5 * ((dx * dx) + (dy * dy)) / (variance * scaling)
      );

      vpsum += vecProbs[posIdx];
      posIdx++;
    }
  }
  return vpsum;
}

// part two of meanshift calculation
export const gpopt2 = (
  responseWidth,
  vecpos,
  updatePosition,
  vecProbs,
  vpsum,
  opj0,
  opj1,
  scaling
) => {
  // for debugging
  // const vecmatrix = [];
  let posIdx = 0;
  let vecsum = 0;
  vecpos[0] = 0;
  vecpos[1] = 0;
  for (let k = 0; k < responseWidth; k++) {
    updatePosition[1] = opj1 + (k * scaling);
    for (let l = 0; l < responseWidth; l++) {
      updatePosition[0] = opj0 + (l * scaling);
      vecsum = vecProbs[posIdx] / vpsum;

      // for debugging
      // vecmatrix[k*responseWidth + l] = vecsum;

      vecpos[0] += vecsum * updatePosition[0];
      vecpos[1] += vecsum * updatePosition[1];
      posIdx++;
    }
  }
  // for debugging
  // return vecmatrix;
}
