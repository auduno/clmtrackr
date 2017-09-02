// webgl setup tests

/*
 * Test whether we can render to floating point texture
 */
var canRenderToFloatTexture = function(webGLContext) {
    var renderingSupported = false;
    var gl = webGLContext;

    // setup the texture
    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 2, 2, 0, gl.RGBA, gl.FLOAT, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);

    // setup the framebuffer
    var framebuffer = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
    gl.framebufferTexture2D(
        gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);

    // check the framebuffer
    var check = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (check == gl.FRAMEBUFFER_COMPLETE) {
        renderingSupported = true;
    }

    // cleanup
    gl.deleteTexture(texture);
    gl.deleteFramebuffer(framebuffer);
    gl.bindTexture(gl.TEXTURE_2D, null);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);

    return renderingSupported
}

export default canRenderToFloatTexture;