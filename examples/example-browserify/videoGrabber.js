var getUserMedia = require('getusermedia');

function videoGrabber(width, height) {
	var v = document.createElement('video');
	var constraints = {
		audio: false
		, video: {
			width: { min: width, ideal: width }
			, height: { min: height, ideal: height}
		}
	}

	v.width = width;
	v.height = height;


	getUserMedia(constraints, function(err, stream) {
		if(err) {
			throw new Error(err);
		}

		// http://www.kirupa.com/html5/accessing_your_webcam_in_html5.htm
		v.src = window.URL.createObjectURL(stream);
	});

	return v;
}

module.exports = videoGrabber;