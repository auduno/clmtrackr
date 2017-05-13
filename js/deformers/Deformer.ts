import { EventEmitter } from 'events';

import { getWebGLContext } from 'clmtrackr/js/utils/webgl';
import {
  getBoundingBox,
  generateTextureVertices
} from 'clmtrackr/js/utils/points';
import { getImageData } from 'clmtrackr/js/utils/image';


abstract class Deformer extends EventEmitter {
  protected _gl: WebGLRenderingContext;

  protected _tracker;
  protected _verticeMap;

  protected _dynamicMaskTexture: boolean;
  protected _maskTextureSrcElement: HTMLElement;
  protected _maskTextureCanvas: HTMLCanvasElement;

  protected _pointBB;
  protected _maskTextureCoord;


  abstract setBackground (element: HTMLElement): void;
  abstract draw (points: number[][]): void;
  abstract drawGrid (): void;


  constructor () {
    super();

    this._dynamicMaskTexture = false;
    this._maskTextureSrcElement = null;
    this._maskTextureCanvas = document.createElement('canvas');
    this._pointBB = null;

    this._maskTextureCoord = null;
  }


  public init (canvas: HTMLCanvasElement): void {
    if (!canvas) {
      throw new Error('canvas parameter is falsey');
    }
    this._gl = getWebGLContext(canvas);
    if (!this._gl) {
      throw new Error('Could not get a webgl context; have you already tried getting a 2d context on this canvas?');
    }
  }

  public setTracker (tracker): void {
    this._tracker = tracker;
    // Set verts for this mask
    this._verticeMap = tracker.model.path.vertices;
  }

  public setMaskTexture (element: HTMLElement): void {
    this._maskTextureSrcElement = element;

    const tagName = this._maskTextureSrcElement.tagName;
    if (tagName === 'CANVAS') {
      // Use the element as texture (its dynamic!)
      this._dynamicMaskTexture = true;
    } else {
      // We need a still frame from it
      this._dynamicMaskTexture = false;
      this.updateMaskTexture();
    }
  }

  public setPoints (points: number[][]): void {
    if (!points) {
      throw new Error('points is falsey');
    }

    // Find texture cropping from mask points
    this._pointBB = getBoundingBox(points);

    // offset points by bounding box
    const nupoints = points.map(p => [
      p[0] - this._pointBB.minX,
      p[1] - this._pointBB.minY
    ]);

    // create UVs based on map points
    this._maskTextureCoord = generateTextureVertices(
      nupoints,
      this._verticeMap,
      1 / this._pointBB.width,
      1 / this._pointBB.height
    );

    this.updateMaskTexture();
  }

  protected updateMaskTexture (): HTMLElement {
    if (
      !this._maskTextureSrcElement ||
      !this._pointBB
    ) {
      return null;
    }

    this.emit('maskReady');

    if (!this._dynamicMaskTexture) {
      // Draw the srcElement to the mask texture canvas
      const {
        minX, minY, width, height
      } = this._pointBB;

      const maskImage = getImageData(
        this._maskTextureSrcElement,
        minX,
        minY,
        width,
        height
      );

      const canvas = this._maskTextureCanvas;
      const ctx = canvas.getContext('2d');
      canvas.width = width;
      canvas.height = height;
      ctx.putImageData(maskImage, 0, 0);

      return canvas;
    } else {
      return this._maskTextureSrcElement;
    }
  }

  public getGLContext (): WebGLRenderingContext {
    return this._gl;
  }

  public clear (): void {
    const gl = this.getGLContext();
    gl.clear(gl.COLOR_BUFFER_BIT);
  }
}

export default Deformer;
