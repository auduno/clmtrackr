module.exports = function(grunt) {
	grunt.initConfig({
		concat: {
			dist: {
				src: [
						'../examples/ext_js/jsfeat_detect.js',
						'../examples/ext_js/mosse.js',
						'../examples/ext_js/left_eye_filter.js',
						'../examples/ext_js/right_eye_filter.js',
						'../examples/ext_js/nose_filter.js',
						'./clm.js',
						'./svmfilter_webgl.js',
						'./svmfilter_fft.js',
						'./svmfilter_mosse.js',
						],
				dest: './clmtrackr.js'
			}
		}, 
		min: {
			dist: {
				src: ['./clmtrackr.js'],
				dest: './clmtrackr.min.js'
			}
		}
	});

	// Default task.
	grunt.registerTask('default', 'concat min');
};