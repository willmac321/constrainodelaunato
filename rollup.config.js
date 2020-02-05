import { terser } from 'rollup-plugin-terser'
import buble from 'rollup-plugin-buble'
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

const bubleConfig = { transforms: { dangerousForOf: true } }

export default [
  config('constrainodelaunato.js', [resolve(), buble(bubleConfig)]),
  config('constrainodelaunato.min.js', [resolve(), terser(), buble(bubleConfig)])
]
