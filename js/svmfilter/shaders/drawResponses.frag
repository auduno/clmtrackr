precision mediump float;

// our responses
uniform sampler2D u_responses;

// the texCoords passed in from the vertex shader.
varying vec2 v_texCoord;
varying float v_select;

const vec4 bit_shift = vec4(256.0 * 256.0 * 256.0, 256.0 * 256.0, 256.0, 1.0);
const vec4 bit_mask  = vec4(0.0, 1.0 / 256.0, 1.0 / 256.0, 1.0 / 256.0);

// packing code from here http://stackoverflow.com/questions/9882716/packing-float-into-vec4-how-does-this-code-work
void main() {
  vec4 colorSum = texture2D(u_responses, v_texCoord);
  float value = 0.0;
  if (v_select < 0.1) {
    value = colorSum[0];
  } else if (v_select > 0.9 && v_select < 1.1) {
    value = colorSum[1];
  } else if (v_select > 1.9 && v_select < 2.1) {
    value = colorSum[2];
  } else if (v_select > 2.9 && v_select < 3.1) {
    value = colorSum[3];
  } else {
    value = 1.0;
  }

  vec4 res = fract(value * bit_shift);
  res -= res.xxyz * bit_mask;

  //gl_FragColor = vec4(value, value, value, value);
  //gl_FragColor = vec4(1.0, value, 1.0, 1.0);
  gl_FragColor = res;
}
