module.exports = function(grunt) {
	grunt.initConfig({
		concat: {
			dist: {
				src: [
						'./license.js',
						'./clm.js',
						'./svmfilter_webgl.js',
						'./svmfilter_fft.js',
						'./mossefilter.js',
						'../examples/ext_js/left_eye_filter.js',
						'../examples/ext_js/right_eye_filter.js',
						'../examples/ext_js/nose_filter.js',
						'../examples/ext_js/numeric-1.2.6.js',
						'../examples/ext_js/ccv.js',
						'../examples/ext_js/cascade.js',
						'../examples/ext_js/mosse.js',
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