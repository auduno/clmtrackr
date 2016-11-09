var clmtrackr = require('../../js/index.js');
var model = require('../../models/model_pca_20_svm.json');
var videoGrabber = require('./videoGrabber.js');
var raf = require('raf');

var width = 600;
var height = 400;

var tracker = new clmtrackr.tracker();
var video = videoGrabber(width, height);
var canvas = document.createElement('canvas');
canvas.width = width;
canvas.height = height;
var cc = canvas.getContext('2d');

function setup() {
	video.play();
		
	document.body.appendChild(canvas);

	tracker.init(model);
	tracker.start(video);
}

function update() {
	raf(update);
	cc.clearRect(0, 0, width, height);
	cc.drawImage(video, 0, 0, width, height);
	if(tracker.getCurrentPosition()) {
		tracker.draw(canvas);
	}
}



setup();
update();