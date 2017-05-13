import jsfeat from './jsfeatWithFrontalface';

onmessage = function (e) {
  const buffer = e.data.imageDataBuffer;
  const imageData = new Uint8Array(buffer, 0, e.data.imageDataLength);
  const w = e.data.w;
  const h = e.data.h;

  const imgU8 = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
  const iiSum = new Int32Array((w + 1) * (h + 1));
  const iiSqsum = new Int32Array((w + 1) * (h + 1));
  const iiTilted = new Int32Array((w + 1) * (h + 1));

  const classifier = jsfeat.haar.frontalface;

  // Old findFace
  jsfeat.imgproc.grayscale(imageData, imgU8.data);
  jsfeat.imgproc.equalize_histogram(imgU8, imgU8);
  jsfeat.imgproc.compute_integral_image(imgU8, iiSum, iiSqsum, null);
  let rects = jsfeat.haar.detect_multi_scale(
    iiSum,
    iiSqsum,
    iiTilted,
    null,
    imgU8.cols,
    imgU8.rows,
    classifier,
    1.15,
    2
  );
  rects = jsfeat.haar.group_rectangles(rects, 1);

  const rl = rects.length;

  // console.timeEnd('findFace');
  if (rl <= 0) {
    self.postMessage({ comp: false });
    return;
  }

  let best = rects[0];
  for (let i = 1; i < rl; i++) {
    if (rects[i].neighbors > best.neighbors) {
      best = rects[i]
    } else if (rects[i].neighbors === best.neighbors) {
      if (rects[i].confidence > best.confidence) {
        best = rects[i];
      }
    }
  }

  self.postMessage({ comp: [best] });
  // END Old findFace
};
