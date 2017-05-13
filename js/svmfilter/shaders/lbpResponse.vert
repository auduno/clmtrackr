attribute vec2 a_texCoord;
attribute vec2 a_position;

varying vec2 v_texCoord;

void main() {
   // transform coordinates to regular coordinates
   gl_Position = vec4(a_position, 0.0, 1.0);

   // pass the texCoord to the fragment shader
   v_texCoord = a_texCoord;
}
