import {
  WebGLRenderer,
  Scene,
  PerspectiveCamera,
  Mesh,
  MeshBasicMaterial,
  PlaneGeometry,
  Texture,
  LinearFilter,
  NearestFilter,
  DoubleSide,
  ShaderMaterial,
  BufferGeometry,
  BufferAttribute,
  WebGLRenderTarget
} from 'three';

import { generateTextureVertices } from 'clmtrackr/js/utils/points';
import Deformer from '../Deformer';

import createMaskVS from './shaders/mask.vert';
import createMaskFS from './shaders/mask.frag';


const RAD_TO_DEG = 180 / Math.PI;
const DEG_TO_RAD = Math.PI / 180;


export default class ThreeDeformer extends Deformer {

  private scene: Scene;

  private camera: PerspectiveCamera;

  private renderer: WebGLRenderer;

  private maskMesh: Mesh;
  private bgMesh: Mesh;

  private bgScaleX: number;
  private bgScaleY: number;


  constructor () {
    super();
  }

  public init (canvas: HTMLCanvasElement): void {
    super.init(canvas);

    this.renderer = new WebGLRenderer({
      canvas,
      preserveDrawingBuffer: true
    });
    this.renderer.autoClear = false;
    this.renderer.setSize(canvas.width, canvas.height);

    this.scene = new Scene();

    this.camera = new PerspectiveCamera(
      75,
      canvas.width / canvas.height,
      1,
      10000
    );
    this.camera.position.z = 1000;

    // Make the background
    const tan = Math.tan(this.camera.fov / 2 * DEG_TO_RAD) * (2 * this.camera.position.z);
    const bgGeom = new PlaneGeometry(
      tan * this.camera.aspect,
      tan
    );
    const bgMat = new MeshBasicMaterial({
      color: 0x00ff00,
      wireframe: true
    });
    this.bgMesh = new Mesh(bgGeom, bgMat);
    this.scene.add(this.bgMesh);

    this.bgScaleX = bgGeom.parameters.width / canvas.width;
    this.bgScaleY = bgGeom.parameters.height / canvas.height;

    // Mask the mask geometry
    const maskGeom = new BufferGeometry();
    const maskMat = new ShaderMaterial({
      uniforms: {
        texture: { value: null },
        bgTexture: { value: null },
        bgWidth: { value: canvas.width },
        bgHeight: { value: canvas.height }
      },
      vertexShader: createMaskVS(),
      fragmentShader: createMaskFS()
    });

    this.maskMesh = new Mesh(maskGeom, maskMat);

    // Dont add mask to scene until it is ready
    this.once('maskReady', () => {
      this.scene.add(this.maskMesh);
    })
  }

  public setBackground (element: HTMLElement): void {
    const texture = new Texture(element);
    texture.minFilter = LinearFilter;
    const bgMaterial = this.bgMesh.material;
    bgMaterial.map = texture;

    const maskBgTexture = this.maskMesh.material.uniforms.bgTexture;
    maskBgTexture.value = texture;
    maskBgTexture.needsUpdate = true;

    // Un-set the defaults
    bgMaterial.wireframe = false;
    bgMaterial.color.set(0xffffff);
  }

  protected updateMaskTexture (): HTMLElement {
    const srcElement = super.updateMaskTexture();
    if (!srcElement) { return; }
    // Update mask texture
    const texture = new Texture(srcElement);
    texture.minFilter = LinearFilter;
    texture.needsUpdate = true;

    const maskMaterial = this.maskMesh.material;
    maskMaterial.map = texture;
    maskMaterial.side = DoubleSide;
    // Un-set the defaults
    maskMaterial.wireframe = false;

    // Update the shader uniform
    const uTexture = maskMaterial.uniforms.texture;
    uTexture.value = texture;
    uTexture.needsUpdate = true;

    return;
  }

  public setPoints (points: number[][]): void {
    super.setPoints(points);

    const geom = this.maskMesh.geometry;

    const faceCount = Math.floor(this._maskTextureCoord.length / 6 * 3)
    // Initialize the verts
    geom.addAttribute(
      'position',
      new BufferAttribute(new Float32Array(faceCount * 3), 3)
    );

    // Initialize the UVs
    const faceVertexUvs = new Float32Array(faceCount * 3);
    for (let i = 0; i < this._maskTextureCoord.length; i += 2) {
      faceVertexUvs[i] = this._maskTextureCoord[i];
      faceVertexUvs[i + 1] = 1 - this._maskTextureCoord[i + 1];
    }

    geom.addAttribute(
      'uv',
      new BufferAttribute(faceVertexUvs, 2)
    );
    geom.attributes.uv.needsUpdate = true;
  }

  private updateMaskGeom (points: number[][]): void {
    const maskVertices = generateTextureVertices(points, this._verticeMap);

    const geom = this.maskMesh.geometry;
    const position = geom.attributes.position;

    const bgW = this.bgMesh.geometry.parameters.width;
    const bgH = this.bgMesh.geometry.parameters.height;
    const offsetX = bgW * -0.5;
    const offsetY = bgH * -0.5;

    const verts = position.array;
    let vertIndex = 0;
    for (let i = 0; i < maskVertices.length; i += 2) {
      verts[vertIndex++] = (maskVertices[i] * this.bgScaleX) + offsetX;
      verts[vertIndex++] = (bgH - (maskVertices[i + 1] * this.bgScaleY)) + offsetY;
      verts[vertIndex++] = 1;
    }

    position.needsUpdate = true;
  }

  public draw (points: number[][]): void {
    // Update the scene
    // TODO: this should move to a separate tick function to avoid rendering
    // hiccups
    this.updateMaskGeom(points);

    // Update bg texture
    const bgTex = this.bgMesh.material.map;
    if (bgTex) {
      bgTex.needsUpdate = true;
      this.maskMesh.material.uniforms.bgTexture.needsUpdate = true;
    }

    this.renderer.render(this.scene, this.camera);
  }

  public drawGrid (): void { }
}
