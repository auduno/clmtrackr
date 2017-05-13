precision mediump float;

uniform vec2 u_onePixelPatches;

// our patches
uniform sampler2D u_patches;

// the texCoords passed in from the vertex shader.
varying vec2 v_texCoord;

void main() {
  vec4 topLeft = texture2D(u_patches, v_texCoord + vec2(-$opp0, -$opp1));
  vec4 topMid = texture2D(u_patches, v_texCoord + vec2(0.0, -$opp1));
  vec4 topRight = texture2D(u_patches, v_texCoord + vec2($opp0, -$opp1));
  vec4 midLeft = texture2D(u_patches, v_texCoord + vec2(-$opp0, 0.0));
  vec4 midMid = texture2D(u_patches, v_texCoord);
  vec4 midRight = texture2D(u_patches, v_texCoord + vec2($opp0, 0.0));
  vec4 bottomLeft = texture2D(u_patches, v_texCoord + vec2(-$opp0, $opp1));
  vec4 bottomMid = texture2D(u_patches, v_texCoord + vec2(0.0, $opp1));
  vec4 bottomRight = texture2D(u_patches, v_texCoord + vec2($opp0, $opp1));
  vec4 lbp = step(midMid, midRight) * 1.0 + step(midMid, topRight) * 2.0 + step(midMid, topMid) * 4.0;
  lbp = lbp + step(midMid, topLeft) * 8.0 + step(midMid, midLeft) * 16.0 + step(midMid, bottomLeft) * 32.0;
  lbp = lbp + step(midMid, bottomMid) * 64.0 + step(midMid, bottomRight) * 128.0;
  gl_FragColor = lbp;
}
