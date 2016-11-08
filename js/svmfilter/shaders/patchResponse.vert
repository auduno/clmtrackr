attribute vec2 a_texCoord;
attribute vec2 a_position;

const vec2 u_resolution = vec2($resolutionX, $resolutionY);
const float u_patchHeight = $patchHeight;
const float u_filterHeight = $filterHeight;
const vec2 u_midpoint = vec2(0.5, $midpointY);

varying vec2 v_texCoord;
varying vec2 v_texCoordFilters;

void main() {
   // convert the rectangle from pixels to 0.0 to 1.0
   vec2 zeroToOne = a_position / u_resolution;

   // convert from 0->1 to 0->2
   vec2 zeroToTwo = zeroToOne * 2.0;

   // convert from 0->2 to -1->+1 (clipspace)
   vec2 clipSpace = zeroToTwo - 1.0;

   // transform coordinates to regular coordinates
   gl_Position = vec4(clipSpace * vec2(1.0, 1.0), 0, 1);

   // pass the texCoord to the fragment shader
   v_texCoord = a_texCoord;

   // set the filtertexture coordinate based on number filter to use
   v_texCoordFilters = u_midpoint + vec2(0.0, u_filterHeight * floor(a_texCoord[1]/u_patchHeight));
}
