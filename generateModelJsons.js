var fs = require('fs')
var path = require('path')

var modelsPath = path.resolve(__dirname, 'models')

var files = fs.readdirSync(modelsPath).filter(function (file) {
  return path.extname(file) === '.js';
}).map(function (file) {
  return path.join(modelsPath, file)
})

console.log('Found', files)

files.forEach(function (fpath) {
  var contents = fs.readFileSync(fpath, 'utf8');
  eval(contents);
  fs.writeFileSync(
    fpath + 'on',
    JSON.stringify(pModel)
  )
})

console.log('Done')
