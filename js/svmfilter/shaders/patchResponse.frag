precision mediump float;

const vec2 u_onePixelPatches = vec2($onePixelPatchesX, $onePixelPatchesY);
const vec2 u_onePixelFilters = vec2($onePixelFiltersX, $onePixelFiltersY);
const float u_halffilterwidth = $halfFilterWidth;
const float u_halffilterheight = $halfFilterHeight;

// our patches
uniform sampler2D u_patches;
// our filters
uniform sampler2D u_filters;

// the texCoords passed in from the vertex shader.
varying vec2 v_texCoord;
varying vec2 v_texCoordFilters; // this should give us correct filter

void main() {
  vec4 colorSum = vec4(0.0, 0.0, 0.0, 0.0);
  vec4 maxn = vec4(0.0, 0.0, 0.0, 0.0);
  vec4 minn = vec4(256.0, 256.0, 256.0, 256.0);
  vec4 scale = vec4(0.0, 0.0, 0.0, 0.0);
  vec4 patchValue = vec4(0.0, 0.0, 0.0, 0.0);
  vec4 filterValue = vec4(0.0, 0.0, 0.0, 0.0);
  vec4 filterTemp = vec4(0.0, 0.0, 0.0, 0.0);
  for (int w = 0; w < $filterWidth; w++) {
    for (int h = 0; h < $filterHeight; h++) {
      patchValue = texture2D(u_patches, v_texCoord + u_onePixelPatches * vec2(float(w)-u_halffilterwidth, float(h)-u_halffilterheight));
      filterValue = texture2D(u_filters, v_texCoordFilters + u_onePixelFilters * vec2(float(w)-u_halffilterwidth, float(h)-u_halffilterheight));
      maxn = max(patchValue, maxn);
      minn = min(patchValue, minn);
      colorSum += patchValue*filterValue;
      filterTemp += filterValue;
    }
  }
  scale = maxn-minn;
  colorSum = (colorSum-(minn*filterTemp))/scale;
  // logistic transformation
  colorSum = 1.0/(1.0 + exp(- (colorSum) ));
  gl_FragColor = colorSum;
}
