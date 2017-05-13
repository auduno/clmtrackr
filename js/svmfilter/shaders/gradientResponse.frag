precision mediump float;

uniform vec2 u_onePixelPatches;

// our patches
uniform sampler2D u_patches;

// the texCoords passed in from the vertex shader.
varying vec2 v_texCoord;

void main() {
  vec4 bottomLeft = texture2D(u_patches, v_texCoord + vec2(-$opp0, $opp1));
  vec4 bottomRight = texture2D(u_patches, v_texCoord + vec2($opp0, $opp1));
  vec4 topLeft = texture2D(u_patches, v_texCoord + vec2(-$opp0, -$opp1));
  vec4 topRight = texture2D(u_patches, v_texCoord + vec2($opp0, -$opp1));
  vec4 dx = (
    bottomLeft +
    (texture2D(u_patches, v_texCoord + vec2(-$opp0, 0.0))*vec4(2.0,2.0,2.0,2.0)) +
    topLeft -
    bottomRight -
    (texture2D(u_patches, v_texCoord + vec2($opp0, 0.0))*vec4(2.0,2.0,2.0,2.0)) -
    topRight)/4.0;
  vec4 dy = (
    bottomLeft +
    (texture2D(u_patches, v_texCoord + vec2(0.0, $opp1))*vec4(2.0,2.0,2.0,2.0)) +
    bottomRight -
    topLeft -
    (texture2D(u_patches, v_texCoord + vec2(0.0, -$opp1))*vec4(2.0,2.0,2.0,2.0)) -
    topRight)/4.0;
  vec4 gradient = sqrt((dx*dx) + (dy*dy));
  gl_FragColor = gradient;
}
