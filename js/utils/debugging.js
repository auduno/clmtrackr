// function to draw pixeldata on some canvas, only used for debugging
export const drawData = (canvasContext, data, width, height, transposed, drawX, drawY) => {
  const psci = canvasContext.createImageData(width, height);
  const pscidata = psci.data;
  for (let j = 0; j < width * height; j++) {
    let val;
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
