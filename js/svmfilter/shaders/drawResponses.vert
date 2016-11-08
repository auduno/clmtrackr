attribute vec2 a_texCoord_draw;
attribute vec2 a_position_draw;
attribute float a_patchChoice_draw;

uniform vec2 u_resolutiondraw;

varying vec2 v_texCoord;
varying float v_select;

void main() {
   // convert the rectangle from pixels to 0.0 to 1.0
   vec2 zeroToOne = a_position_draw / u_resolutiondraw;

   // convert from 0->1 to 0->2
   vec2 zeroToTwo = zeroToOne * 2.0;

   // convert from 0->2 to -1->+1 (clipspace)
   vec2 clipSpace = zeroToTwo - 1.0;

   // transform coordinates to regular coordinates
   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);

   // pass the texCoord to the fragment shader
   v_texCoord = a_texCoord_draw;

   v_select = a_patchChoice_draw;
}
