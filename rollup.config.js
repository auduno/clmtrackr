import resolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';
import json from 'rollup-plugin-json';

export default {
	entry: 'src/clm.js',
	targets: [
		{
			format: 'umd',
			moduleName: 'clm',
			dest: 'build/clmtrackr.js'
		},
		{
			format: 'es',
			dest: 'build/clmtrackr.module.js'
		}
	],
	plugins: [
		resolve({
			module: true,
			main: true
		}),
		commonjs(),
		json()
	]
};