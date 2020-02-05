import { terser } from 'rollup-plugin-terser'
import resolve from '@rollup/plugin-node-resolve'

const config = (file, plugins) => ({
  input: 'constrainment.js',
  plugins,
  output: {
    name: 'ConstrainoDelaunato',
    format: 'umd',
    file
  }
})

export default [
  config('constrainodelaunato.js', [resolve()]),
  config('constrainodelaunato.min.js', [resolve(), terser()])
]
