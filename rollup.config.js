import resolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';

export default {
	entry: 'js/clm.js',
	targets: [
		{
			format: 'umd',
			moduleName: 'clm',
			dest: 'dist/clmtrackr.js'
		},
		{
			format: 'es',
			dest: 'dist/clmtrackr.module.js'
		}
	],
	plugins: [
		resolve({
			module: true,
			main: true
		}),
		commonjs()
	]
};