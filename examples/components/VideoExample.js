import React from 'react';
import classNames from 'classnames';

import RaisedButton from 'material-ui/RaisedButton';
import Toggle from 'material-ui/Toggle';

import {
  SUPPORTED_CODECS,
  supportsVideo,
  supportsUserMedia,
  loadVideo
} from 'clmtrackr/js/utils/video';
import {
  requestAnimFrame,
  cancelRequestAnimFrame
} from 'clmtrackr/js/utils/anim';

import TrackerContainer from 'clmtrackr/ui/container/TrackerContainer';

import './VideoExample.styl';


export default class VideoExample extends React.Component {
  constructor () {
    super();
    this.state = {
      tracker: null,

      showAltInstructions: false,
      startButtonDisabled: true,
      showScore: false,

      mediaSrc: null,
      hideMedia: false,
      useStockVideo: true,

      title: 'Cool Example'
    };

    this.mediaSize = { width: 400, height: 300 };
    this.webmVideoSrc = null;
    this.oggVideoSrc = null;
    this.mp4VideoSrc = null;

    this._onFrameAnimId = null;
  }

  newTracker () {
    throw new Error('STUB');
  }

  componentDidMount () {
    const tracker = this.newTracker();
    tracker.init();

    this.setState({ tracker });

    const trackerContainer = this.refs.trackerContainer;
    const video = trackerContainer.refs.media;
    video.addEventListener('canplay', () => {
      this.setState({ startButtonDisabled: false });
    });

    setTimeout(() => {
      this._setupVideoStream();
    });
  }

  componentWillUnmount () {
    // Stop tracker
    const tracker = this.state.tracker;
    tracker.stop();

    if (this._onFrameAnimId) {
      cancelRequestAnimFrame(this._onFrameAnimId);
      this._onFrameAnimId = null;
    }
  }

  _setupVideoStream () {
    if (this.state.useStockVideo) {
      this._insertAltVideo();
      return;
    }

    // check for camerasupport
    if (supportsUserMedia()) {
      // set up stream
      loadVideo((err, stream) => {
        if (err) {
          this._insertAltVideo();
          this.setState({ showAltInstructions: true });
          alert('There was some problem trying to fetch video from your webcam, using a fallback video instead.');
          return;
        }

        this.setState({ mediaSrc: stream });
      });
    } else {
      this._insertAltVideo();
      this.setState({ showAltInstructions: true });
      alert('Your browser does not seem to support getUserMedia, using a fallback video instead.');
    }
  }

  _onFrame () {
    // Update overlay
    const trackerContainer = this.refs.trackerContainer;
    if (trackerContainer) {
      // TODO: This should really draw to the same canvas as the deformer,
      // Need to add webgl support to `tracker#draw`
      this.clearOverlayCanvas();
      const tracker = this.state.tracker;
      tracker.draw(trackerContainer.refs.canvas);
    }
    this._onFrameAnimId = requestAnimFrame(this._onFrame.bind(this));
  }

  _resetTracker () {
    const tracker = this.state.tracker;
    this.clearOverlayCanvas();
    tracker.stop();
    tracker.reset();
  }

  clearOverlayCanvas () {
    const trackerContainer = this.refs.trackerContainer;
    const overlayCC = trackerContainer.refs.canvas.getContext('2d');
    overlayCC.clearRect(0, 0, this.mediaSize.width, this.mediaSize.height);
  }

  _insertAltVideo () {
    if (!supportsVideo()) {
      alert('Your browser does not support video');
      return;
    }
    if (this.webmVideoSrc && SUPPORTED_CODECS.webm) {
      this.setState({ mediaSrc: this.webmVideoSrc });
    } else if (this.oggVideoSrc && SUPPORTED_CODECS.ogg) {
      this.setState({ mediaSrc: this.oggVideoSrc });
    } else if (this.mp4VideoSrc && SUPPORTED_CODECS.h264) {
      this.setState({ mediaSrc: this.mp4VideoSrc });
    } else {
      alert('No stock video available for your browser');
    }
  }

  startVideo () {
    // start video
    const trackerContainer = this.refs.trackerContainer;
    const video = trackerContainer.refs.media;
    // start tracking
    this._resetTracker();
    this.state.tracker.start(video);
    this.setState({ showScore: true });
    // start loop to draw face
    // this._onFrame();
  }

  _handleStockVideoToggle (flag) {
    this.setState({ useStockVideo: flag });
    setTimeout(() => this._setupVideoStream());
  }

  renderInstructions () {
    throw new Error('STUB');
  }

  renderAltInstructions () {
    throw new Error('STUB');
  }

  render (opts = {}) {
    const { childControls, children, className } = opts;

    let instructions;
    if (!this.state.showAltInstructions) {
      instructions = (
        <div>
          <p>To try it out:</p>
          {this.renderInstructions()}
        </div>
      );
    } else {
      instructions = (
        <div>
          <p>There was some problem trying to capture your webcamera, please check that your browser supports WebRTC. Using a fallback video instead. try it out:</p>
          {this.renderAltInstructions()}
        </div>
      );
    }

    const startButtonDisabled = this.state.startButtonDisabled;
    let scoreText;
    if (this.state.showScore) {
      scoreText = <p>Please keep head still while model fits</p>;
    }

    return (
      <div
        className={classNames('video-example-cmpt', className, {
          'hide-media': this.state.hideMedia
        })}
      >
        <h1>{this.state.title}</h1>

        <TrackerContainer
          ref='trackerContainer'
          mediaType={'video'}
          mediaSrc={this.state.mediaSrc}
          mediaSize={this.mediaSize}
          showStats={true}
          tracker={this.state.tracker}
          hasDeformer={true}
        />

        <div className='control-row'>
          <RaisedButton
            label={startButtonDisabled ? 'wait, loading video' : 'start'}
            disabled={startButtonDisabled}
            onClick={this.startVideo.bind(this)}
          />

          <Toggle
            label='Use stock video'
            toggled={this.state.useStockVideo}
            style={{ width: 'auto' }}
            onToggle={(event, flag) => this._handleStockVideoToggle(flag)}
          />

          {childControls}
        </div>
        {scoreText}
        <div>
          <p>This is an example of face deformation using the javascript library <a href='https://github.com/auduno/clmtrackr'><em>clmtrackr</em></a>.</p>
          <p>Note that this example needs support for WebGL, and works best in Google Chrome.</p>

          {instructions}
        </div>

        {children}
      </div>
    );
  }
}
