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
						'examples/js/libs/left_eye_filter.js',
						'examples/js/libs/right_eye_filter.js',
						'examples/js/libs/nose_filter.js',
						'examples/js/libs/numeric-1.2.6.js',
						'examples/js/libs/jsfeat-min.js',
						'examples/js/libs/frontalface.js',
						'examples/js/libs/jsfeat_detect.js',
						'examples/js/libs/mosse.js',
						],
				dest: './clmtrackr.js'
			}
		},
		uglify: {
			options: {
				report: 'min',
				preserveComments: 'false',
				mangle: {
					except: ['clmtrackr']
				},
				compress: {
					drop_console: true
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
