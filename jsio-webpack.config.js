'use strict';
const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');

const NODE_ENV = process.env.NODE_ENV || 'development';


const configure = (configurator, options) => {
  configurator.merge({
    entry: {
      clmtrackr: path.resolve(__dirname, 'js'),
      examples: path.resolve(__dirname, 'examples')
    },
    output: {
      filename: '[name].js',
      path: path.resolve(__dirname, 'dist'),
      publicPath: NODE_ENV === 'production' ? '/clmtrackr/' : '/',
      libraryTarget: 'umd',
      library: 'clm'
    },
    resolve: {
      alias: {
        'clmtrackr': path.resolve(__dirname),
        'stats.js': path.resolve(__dirname, 'lib', 'stats.js')
      }
    }
  });

  configurator.plugin('htmlWebpack', HtmlWebpackPlugin, [{
    title: 'Clmtrackr Examples',
    template: path.resolve(__dirname, 'examples', 'index.html'),
    chunks: ['examples']
  }]);

  // Use normal style-loader for stylus
  options.noStylusExtractText = true;

  return configurator;
};


const postConfigure = (configurator, options) => {
  configurator.removePreLoader('eslint');
};


module.exports = {
  configure: configure,
  postConfigure: postConfigure
};
