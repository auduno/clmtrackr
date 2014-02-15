module.exports = function(grunt) {
	grunt.loadNpmTasks('grunt-contrib-concat');
	grunt.loadNpmTasks('grunt-contrib-uglify');

	grunt.initConfig({
		concat: {
			dist: {
				src: [
						'js/license.js',
						'js/clm.js',
						'js/svmfilter_webgl.js',
						'js/svmfilter_fft.js',
						'js/mossefilter.js',
						'examples/ext_js/left_eye_filter.js',
						'examples/ext_js/right_eye_filter.js',
						'examples/ext_js/nose_filter.js',
						'examples/ext_js/numeric-1.2.6.js',
						'examples/ext_js/jsfeat-min.js',
						'examples/ext_js/frontalface.js',
						'examples/ext_js/jsfeat_detect.js',
						'examples/ext_js/mosse.js',
						],
				dest: './clmtrackr.js'
			}
		},
		uglify: {
			options: {
				report: 'gzip',
				preserveComments: 'false',
				mangle: {
					except: ['clmtrackr']
				}
			},
			dist: {
				src: ['./clmtrackr.js'],
				dest: './clmtrackr.min.js'
			}
		}
	});

	// Default task.
	grunt.registerTask('default', [
		'concat',
		'uglify'
	]);
};
