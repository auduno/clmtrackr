import twgl from 'twgl.js/dist/twgl';

import { generateTextureVertices } from 'clmtrackr/js/utils/points';

import Background from './Background';

import createDeformVert from './shaders/deform.vert';
import createDeformFrag from './shaders/deform.frag';

import Deformer from '../Deformer';


export default class TwglDeformer extends Deformer {
  private background: Background;
  private debug: Background;

  private _maskProgramInfo;
  private _maskTextures;
  private _maskUniforms;
  private _maskBufferInfo;

  constructor (params = {}) {
    super();

    twgl.setDefaults({ attribPrefix: 'a_' });

    this.background = new Background(this);
    this.debug = new Background(this);

    this._maskProgramInfo = null;
    this._maskTextures = null;
    this._maskUniforms = null;
    this._maskBufferInfo = null;
  }

  public init (canvas: HTMLCanvasElement): void {
    super.init(canvas);

    const gl = this.getGLContext();
    this._maskProgramInfo = twgl.createProgramInfo(gl, [
      createDeformVert(),
      createDeformFrag()
    ]);
  }

  protected updateMaskTexture (): HTMLElement {
    const srcElement = super.updateMaskTexture();
    if (!srcElement) { return; }

    // Update the webgl stuff
    const gl = this.getGLContext();
    this._maskTextures = twgl.createTextures(gl, {
      mask: {
        mag: gl.LINEAR,
        min: gl.LINEAR,
        src: srcElement
      }
    });

    this._maskUniforms = {
      u_resolution: [gl.drawingBufferWidth, gl.drawingBufferHeight],
      u_sampler: this._maskTextures.mask
    };
  }

  public setBackground (bgElement: HTMLElement): void {
    this.background.setElement(bgElement);
  }

  public draw (points: number[][]): void {
    const gl = this.getGLContext();
    twgl.resizeCanvasToDisplaySize(gl.canvas);
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

    this.background.draw();
    this.drawMask(points);
    this.debug.draw();
  }

  private drawMask (points: number[][]): void {
    const gl = this.getGLContext();

    // create drawvertices based on points
    let vertices = generateTextureVertices(points, this._verticeMap);

    this._maskBufferInfo = twgl.createBufferInfoFromArrays(gl, {
      position: {
        numComponents: 2,
        data: vertices
      },
      texcoord: {
        numComponents: 2,
        data: this._maskTextureCoord
      }
    });

    gl.useProgram(this._maskProgramInfo.program);

    twgl.setBuffersAndAttributes(gl, this._maskProgramInfo, this._maskBufferInfo);
    twgl.setUniforms(this._maskProgramInfo, this._maskUniforms);
    twgl.drawBufferInfo(gl, gl.TRIANGLES, this._maskBufferInfo);
  }

  public drawGrid () {
    // TODO: implement (what is drawGrid?)
  }
}
