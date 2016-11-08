import React from 'react';
import _ from 'lodash';
import { disposableEvent } from 'jsio-event-kit';

import DropDownMenu from 'material-ui/DropDownMenu';
import MenuItem from 'material-ui/MenuItem';

import VideoExample from './VideoExample';

import Tracker from 'clmtrackr/js/Tracker';
import {
  requestAnimFrame,
  cancelRequestAnimFrame
} from 'clmtrackr/js/utils/anim';
import { MASKS, getMask } from 'clmtrackr/examples/masks';
import Deformer from 'clmtrackr/js/deformers/twgl';

import './FaceMaskExample.styl';


import WEBM_VIDEO from '../media/cap13_edit2.webm';
import OGG_VIDEO from '../media/cap13_edit2.ogv';
import MP4_VIDEO from '../media/cap13_edit2.mp4';


export default class FaceMaskExample extends VideoExample {
  constructor () {
    super();
    this.state = _.merge(this.state, {
      selectedMask: null,
      mask: null
    });

    this.mediaSize = { width: 370, height: 288 };
    this._animateRequestId = null;

    this.webmVideoSrc = WEBM_VIDEO;
    this.oggVideoSrc = OGG_VIDEO;
    this.mp4VideoSrc = MP4_VIDEO;
  }

  newTracker () {
    return new Tracker();
  }

  componentDidMount () {
    super.componentDidMount();

    const deformer = new Deformer();
    const trackerContainer = this.refs.trackerContainer;
    deformer.init(trackerContainer.refs.deformerCanvas);

    this.setState({ deformer });

    setTimeout(() => {
      this.state.tracker.convergenceThreshold = 1.8;
      this._handleMaskChange(MASKS[0].id, true);
    });
  }

  componentWillUnmount () {
    super.componentWillUnmount();
    if (this._animateRequestId) {
      cancelRequestAnimFrame(this._animateRequestId);
      this._animateRequestId = null;
    }
  }

  _setupFaceDeformation () {
    const trackerContainer = this.refs.trackerContainer;
    const video = trackerContainer.refs.media;

    // draw face deformation model
    const tracker = this.state.tracker;
    const deformer = this.state.deformer;
    deformer.setTracker(tracker);
    deformer.setMaskTexture(video);
    deformer.setBackground(video);

    const disposable = disposableEvent(tracker, 'converged', () => {
      // hide video and score
      this.setState({ hideMedia: true, showScore: false });
      disposable.dispose();
    });
  }

  _drawGridLoop () {
    // get position of face
    const trackerContainer = this.refs.trackerContainer;
    const media = trackerContainer.refs.media;
    const tracker = this.state.tracker;
    const positions = tracker.getCurrentPosition(media);

    this.clearOverlayCanvas();

    if (positions) {
      // draw current grid
      tracker.draw(trackerContainer.refs.canvas);
    }
    // check whether mask has converged
    var pn = tracker.getConvergence();
    if (pn < 0.4) {
      this._switchMasks();
      this._animateRequestId = requestAnimFrame(this._drawMaskLoop.bind(this));
    } else {
      this._animateRequestId = requestAnimFrame(this._drawGridLoop.bind(this));
    }
  }

  _drawMaskLoop () {
    this.clearOverlayCanvas();
    // get position of face
    const tracker = this.state.tracker;
    const positions = tracker.getCurrentPosition();
    if (positions) {
      // draw mask on top of face
      const deformer = this.state.deformer;
      deformer.draw(positions);
    }
    this._animateRequestId = requestAnimFrame(this._drawMaskLoop.bind(this));
  }

  startVideo () {
    super.startVideo();
    // start loop to draw face
    // this._onFrame();
    this._setupFaceDeformation();
    this._drawGridLoop();
  }

  _switchMasks () {
    const mask = this.state.mask;
    const deformer = this.state.deformer;
    deformer.setPoints(mask.map);

    const image = new Image();
    image.onload = () => {
      deformer.setMaskTexture(image)
    };
    image.src = mask.texture;
  }

  _handleMaskChange (value, skipUpdate) {
    this.setState({
      selectedMask: value,
      mask: getMask(value)
    });
    if (!skipUpdate) {
      setTimeout(() => this._switchMasks());
    }
  }

  renderInstructions () {
    return (
      <ol>
        <li>allow the page to use your webcamera</li>
        <li>make sure that your face is clearly visible in the video, and click start</li>
        <li>keep your face still, and wait till model fits your face and mask is applied</li>
        <li>try out different masks from the dropdown</li>
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
    const dropdownItems = MASKS.map((mask, i) => <MenuItem
      value={mask.id}
      primaryText={mask.name}
      key={i}
    />);

    const childControls = (
      <div className='mask-select'>
        <div>Select Mask</div>
        <DropDownMenu
          value={this.state.selectedMask}
          onChange={(event, index, value) => this._handleMaskChange(value)}
        >
          {dropdownItems}
        </DropDownMenu>
      </div>
    );

    return super.render({
      childControls,
      className: 'face-mask-example-cmpt'
    });
  }
}
