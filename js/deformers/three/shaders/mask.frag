uniform sampler2D texture;
uniform sampler2D bgTexture;

uniform float bgWidth;
uniform float bgHeight;

varying vec2 vUv;


#include ./ColorSpaces;


void main() {
  vec4 bgColor = texture2D(
    bgTexture,
    vec2(
      gl_FragCoord.x / bgWidth,
      gl_FragCoord.y / bgHeight
    )
  );
  vec4 texColor = texture2D(texture, vUv);

  vec3 texHSV = rgb_to_hsv(vec3(texColor.x, texColor.y, texColor.z));
  vec3 bgHSV = rgb_to_hsv(vec3(bgColor.x, bgColor.y, bgColor.z));

  float vDiff = bgHSV.z - texHSV.z;

  gl_FragColor = vec4(hsv_to_rgb(vec3(
    texHSV.x,
    texHSV.y,
    texHSV.z + vDiff * 0.4
  )), 1.0);
}
