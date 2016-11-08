import React from 'react';
import Dat, { DatNumber, DatBoolean } from 'react-dat-gui';
import '!style!css!react-dat-gui/lib/Dat.css';
import _ from 'lodash';

import VideoExample from './VideoExample';

import Tracker from 'clmtrackr/js/Tracker';
import { ParamHolder } from 'clmtrackr/js/deformers/presets';
import Deformer from 'clmtrackr/js/deformers/twgl';

import './FaceDeformationVideoExample.styl';


import WEBM_VIDEO from '../media/franck.webm';
import OGG_VIDEO from '../media/franck.ogv';


export default class FaceDeformationVideoExample extends VideoExample {
  constructor () {
    super();
    this.state = _.merge(this.state, {
      deformer: null,
      data: new ParamHolder(),
      positions: null,
      title: 'Face deformation'
    });

    this.mediaSize = { width: 370, height: 288 };
    this.webmVideoSrc = WEBM_VIDEO;
    this.oggVideoSrc = OGG_VIDEO;
    this.mp4VideoSrc = null;
  }

  newTracker () {
    return new Tracker({ stopOnConvergence: true });
  }

  componentDidMount () {
    super.componentDidMount();

    const deformer = new Deformer();
    const trackerContainer = this.refs.trackerContainer;
    deformer.init(trackerContainer.refs.deformerCanvas);

    this.setState({ deformer });

    setTimeout(() => {
      const tracker = this.state.tracker;
      tracker.on('converged', () => {
        this._setupFaceDeformation();
      });
      tracker.convergenceThreshold = 1.8;

      this._setupVideoStream();
    });
  }

  _setupFaceDeformation () {
    // hide message element
    this.setState({ showScore: false });

    const trackerContainer = this.refs.trackerContainer;
    const video = trackerContainer.refs.media;
    video.pause();

    // draw face deformation model
    const tracker = this.state.tracker;
    const deformer = this.state.deformer;
    deformer.load(
      video,
      tracker.getCurrentPosition(),
      tracker,
      null
    );

    // hide video
    this.setState({ hideMedia: true });

    setTimeout(() => { this._drawDeformedFace(); });
  }

  _drawDeformedFace () {
    const rotation = 0;
    const scale = 3;
    const xOffset = -10;
    const yOffset = 0;

    const tracker = this.state.tracker;
    const parameters = tracker.getCurrentParameters();
    parameters[0] = scale * Math.cos(rotation) - 1;
    parameters[1] = scale * Math.sin(rotation);
    parameters[2] = xOffset;
    parameters[3] = yOffset;

    const data = this.state.data;
    for (let i = 0; i < 20; i++) {
      parameters[i + 4] = data['param' + (i + 1)];
    }

    const positions = tracker.calculatePositions(parameters);
    const deformer = this.state.deformer;
    deformer.clear();
    if (data.drawFace) {
      deformer.draw(positions);
    }
    if (data.drawGrid) {
      deformer.drawGrid(positions);
    }
  }

  startVideo () {
    super.startVideo();
    // start loop to draw face
    this._onFrame();
  }

  _handleDatUpdate (dataUpdated) {
    this.setState({ data: dataUpdated });
    setTimeout(() => this._drawDeformedFace());
  }

  renderInstructions () {
    return (
      <ol>
        <li>allow the page to your use webcamera</li>
        <li>make sure that your face is clearly visible in the video, and click start</li>
        <li>keep your head still until the model fits your face properly</li>
        <li>choose a preset or fiddle with the parameters on the right hand side, to deform the captured face</li>
      </ol>
    );
  }

  renderAltInstructions () {
    return (
      <ol>
        <li>click start</li>
        <li>wait till model fits the face and a mask is applied</li>
        <li>try out different masks from the dropdown</li>
      </ol>
    );
  }

  render () {
    const datControls = [];
    for (let i = 1; i <= 20; i++) {
      const min = -20;
      const max = 20;
      datControls.push(
        <DatNumber
          path={'param' + i}
          label={'Param ' + i}
          min={min}
          max={max}
          step={0.1}
          key={i}
        />
      );
    }
    datControls.push(<DatBoolean path='drawGrid' label='Draw Grid' key={50} />);
    datControls.push(<DatBoolean path='drawFace' label='Draw Face' key={51} />);

    const children = (
      <Dat
        data={this.state.data}
        onUpdate={this._handleDatUpdate.bind(this)}
      >
        {datControls}
      </Dat>
    );

    return super.render({
      children,
      className: 'face-deformation-video-example-cmpt'
    });
  }
}
