// requires mosse.js

var mosseFilterResponses = function() {
  
  var filters = [];
  var responses = [];
  var num_Patches = 0;
  
  this.init = function(filter_input, numPatches, filterWidth, filterHeight) {
    // load filters, make fft ready
    
    for (var i = 0;i < numPatches;i++) {
      var temp = {};
      temp.width = filterWidth;
      temp.height = filterHeight;
      var filterLength = filterWidth*filterHeight
      var flar_fi0 = new Float64Array(filterLength);
      var flar_fi1 = new Float64Array(filterLength);
      for (var j = 0;j < filterLength;j++) {
        flar_fi0[j] = filter_input[i][0][j];
        flar_fi1[j] = filter_input[i][1][j];
      }
      temp.real = flar_fi0;
      temp.imag = flar_fi1;
      filters[i] = new mosseFilter();
      filters[i].load(temp);
    }
    
    num_Patches = numPatches;
  }
  
  this.getResponses = function(patches) {
    for (var i = 0;i < num_Patches;i++) {
      responses[i] = filters[i].getResponse(patches[i]);
      //responses[i] = logisticResponse(responses[i]);
      responses[i] = normalizeFilterMatrix(responses[i]);
    }
    
    return responses;
  }
  
  var logisticResponse = function(response) {
    // create probability by doing logistic transformation
    var filter_size = response.length;
    for (var j = 0;j < filter_size;j++) {
      response[j] = 1.0/(1.0 + Math.exp(- (response[j]-1.0) ));
    }
    return response;
  }
  
  var normalizeFilterMatrix = function(response) {
    // normalize responses to lie within [0,1]
    var msize = response.length;
    var max = 0;
    var min = 1;
    
    for (var i = 0;i < msize;i++) {
      max = response[i] > max ? response[i] : max;
      min = response[i] < min ? response[i] : min;
    }
    var dist = max-min;
    
    if (dist == 0) {
      console.log("a patchresponse was monotone, causing normalization to fail. Leaving it unchanged.")
      response = response.map(function() {return 1});
    } else {
      for (var i = 0;i < msize;i++) {
        response[i] = (response[i]-min)/dist;
      }
    }
    
    return response
  }
}
