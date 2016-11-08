import EventEmitter from 'events';
// simple wrapper for jsfeat face detector
import findFaceWorker from './findFace.worker';

// import frontalface from '../filters/frontalface.json';


export default class JsfeatFace extends EventEmitter {
  /**
   * @param  {Canvas|Image|Video} image
   */
  constructor () {
    super();

    this._waitingForResponse = false;
    this.worker = findFaceWorker();

    // Listen for messages coming out of the web worker
    this.worker.addEventListener('message', (e) => {
      if (e.data.type === 'console') {
        console[e.data.func].apply(window, e.data.args);
        return;
      }

      this._waitingForResponse = false;
      this.emit('faceDetected', e.data.comp);
    }, false);
  }

  findFace (image) {
    if (!image) {
      throw new Error('Image is falsey');
    }
    if (this._waitingForResponse) {
      throw new Error('Already finding face');
    }
    this._waitingForResponse = true;

    let workCtx;
    const w = image.width;
    const h = image.height;

    if (image.tagName === 'VIDEO' || image.tagName === 'IMG') {
      const workCanvas = document.createElement('canvas');
      workCanvas.height = image.width;
      workCanvas.width = image.height;
      workCtx = workCanvas.getContext('2d');
      // Draw a single frame
      workCtx.drawImage(image, 0, 0);
    } else if (image.tagName === 'CANVAS') {
      workCtx = image.getContext('2d');
    } else {
      throw new Error('unknown image tagName: ' + image.tagName);
    }

    // img_u8 = new jsfeat.matrix_t(w, h, jsfeat.U8_t | jsfeat.C1_t);
    // ii_sum = new Int32Array((w+1)*(h+1));
    // ii_sqsum = new Int32Array((w+1)*(h+1));
    // ii_tilted = new Int32Array((w+1)*(h+1));

    // var classifier = frontalface;

    let imageData;
    try {
      imageData = workCtx.getImageData(0, 0, w, h);
    } catch (e) {
      console.warn(
        'Could not getImageData, is your element too large?',
        `w= ${w} h= ${h}`
      );
      console.error(e);
    }

    // console.time('findFace');
    // Send the underlying ArrayBuffer to worker
    const imageDataBuffer = imageData.data.buffer;
    const message = {
      w: w,
      h: h,
      imageDataBuffer: imageDataBuffer,
      imageDataLength: imageData.length
    };
    this.worker.postMessage(message, [message.imageDataBuffer]);
  }

}
