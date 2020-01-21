import {terser} from 'rollup-plugin-terser';
import buble from 'rollup-plugin-buble';
import resolve from '@rollup/plugin-node-resolve';

const config = (file, plugins) => ({
    input: 'constrainment.js',
    output: {
        name: 'ConstrainoDelaunato',
        format: 'umd',
        file
    },
	plugins: [
		resolve()
	]
});

const bubleConfig = {transforms: {dangerousForOf: true}};

export default [
    config('constrainodelaunato.js', [buble(bubleConfig)]),
    config('constrainodelaunato.min.js', [terser(), buble(bubleConfig)])
];
