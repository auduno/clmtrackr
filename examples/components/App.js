import React, { PropTypes } from 'react';
import { connect } from 'react-redux';

import { getExample } from 'clmtrackr/examples/reducers/examples';

import MenuDrawer from './MenuDrawer';
import SimpleExample from './SimpleExample';
import ClmImageExample from './ClmImageExample';
import ClmVideoExample from './ClmVideoExample';
import FaceDeformationVideoExample from './FaceDeformationVideoExample';
import FaceMaskExample from './FaceMaskExample';

import './App.styl'


const EXAMPLE_TO_CMPT = {
  simple: SimpleExample,
  clmImage: ClmImageExample,
  clmVideo: ClmVideoExample,
  faceDeformStill: FaceDeformationVideoExample,
  faceMask: FaceMaskExample
}


class App extends React.Component {
  render () {
    const exampleCtor = EXAMPLE_TO_CMPT[this.props.activeExample];
    let example;
    if (exampleCtor) {
      example = React.createElement(exampleCtor);
    } else {
      const curExample = getExample(this.props.activeExample);
      example = (
        <div className='coming-soon-wrapper'>
          <h1>Coming Soon!</h1>
          <p>The <a href={curExample.original}>original demo</a> is presented below:</p>
          <iframe src={curExample.original} />
        </div>
      );
    }

    return (
      <div className='app-cmpt'>
        <MenuDrawer />
        <div className='example-wrapper'>
          {example}
        </div>

        <a href='https://github.com/auduno/clmtrackr'>
          <img
            style={{ position: 'absolute', top: 0, right: 0, border: 0 }}
            src='https://s3.amazonaws.com/github/ribbons/forkme_right_green_007200.png'
            alt='Fork me on GitHub'
          />
        </a>
      </div>
    );
  }
}

App.propTypes = {
  activeExample: PropTypes.string
};

const mapStateToProps = state => {
  return {
    activeExample: state.examples.activeExample
  };
};

const mapDispatchToProps = dispatch => {
  return { };
};

export default connect(mapStateToProps, mapDispatchToProps)(App);
