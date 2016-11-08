// Action types
export const SET_EXAMPLE = 'examples/SET_EXAMPLE';

export const HEADER_TITLES = {
  example: 'Examples',
  tool: 'Tools'
};

export const EXAMPLES = [
  {
    id: 'simple',
    name: 'Simple'
  },
  {
    id: 'clmImage',
    name: 'Clm Image'
  },
  {
    id: 'clmVideo',
    name: 'Clm Video'
  },
  {
    id: 'clmVideoResponses',
    name: 'Clm Video Responses',
    original: 'https://auduno.github.io/clmtrackr/clm_video_responses.html'
  },
  {
    id: 'clmEmotionDetection',
    name: 'Emotion Detection',
    original: 'https://auduno.github.io/clmtrackr/examples/clm_emotiondetection.html'
  },
  {
    id: 'clmGenderDetection',
    name: 'Gender Detection',
    original: 'https://auduno.github.io/clmtrackr/examples/clm_genderdetection.html'
  },
  {
    id: 'faceDeformStill',
    name: 'Face Deformation Still'
  },
  {
    id: 'faceMask',
    name: 'Face Mask'
  },
  {
    id: 'caricature',
    name: 'Caricature',
    original: 'https://auduno.github.io/clmtrackr/examples/caricature.html'
  },
  {
    id: 'faceDeformVideo',
    name: 'Face Deformation Video',
    original: 'https://auduno.github.io/clmtrackr/examples/facedeform.html'
  },
  {
    id: 'classViewer',
    name: 'Class Viewer',
    type: 'tool',
    original: 'https://auduno.github.io/clmtrackr/examples/classviewer.html'
  },
  {
    id: 'modelViewerPca',
    name: 'Model Viewer: pca',
    type: 'tool',
    original: 'https://auduno.github.io/clmtrackr/examples/modelviewer_pca.html'
  },
  {
    id: 'modelViewerSpca',
    name: 'Model Viewer: spca',
    type: 'tool',
    original: 'https://auduno.github.io/clmtrackr/examples/modelviewer_spca.html'
  },
  {
    id: 'paramModel',
    name: 'Model Preview: clm pca',
    type: 'tool',
    original: 'https://auduno.github.io/clmtrackr/docs/param_model/clm_pca.html'
  }
];


export const getExample = (id) => {
  return _.find(EXAMPLES, { id });
};


// Action generators
export const setExample = (example) => {
  if (typeof example === 'object') {
    example = example.id;
  }
  return { type: SET_EXAMPLE, value: example };
};


// The reducers
const DEFAULT_STATE = {
  activeExample: EXAMPLES[7].id
};

export default function (state = DEFAULT_STATE, action) {
  switch (action.type) {
    case SET_EXAMPLE:
      return Object.assign({}, state, {
        activeExample: action.value
      });
    default:
      return state
  }
};
