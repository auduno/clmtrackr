attribute vec4 a_position;
attribute vec2 a_texcoord;

varying vec2 v_texCoord;

void main() {
  v_texCoord = a_texcoord;
  gl_Position = a_position;
}
