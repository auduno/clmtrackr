import React, { PropTypes } from 'react';

import Stats from 'stats.js/src/Stats';

import { setVideoSrc } from 'clmtrackr/js/utils/video';

import './TrackerContainer.styl';


export default class TrackerContainer extends React.Component {
  constructor () {
    super();
    this._stats = null;
    this._boundUpdateStats = this._updateStats.bind(this);

    this._mediaSrcStream = null;
  }

  _updateStats () {
    const stats = this._stats;
    if (!stats) { return; }
    stats.update();
  }

  componentDidMount () {
    // Initialize the tracker
    this.componentDidUpdate({ tracker: null });
  }

  componentWillUnmount () {
    this._shutdownMediaStream();
  }

  _addTrackerListeners (tracker) {
    tracker.on('iteration', this._boundUpdateStats);
  }

  _removeTrackerListeners (tracker) {
    if (!tracker) { return; }
    tracker.removeListener('iteration', this._boundUpdateStats);
  }

  shouldComponentUpdate (nextProps) {
    return (
      nextProps.tracker !== this.props.tracker ||
      nextProps.showStats !== this.props.showStats ||
      nextProps.mediaSrc !== this.props.mediaSrc
    );
  }

  componentDidUpdate (prevProps) {
    // Update tracker listeners for stats
    const oldTracker = prevProps.tracker;
    const tracker = this.props.tracker;
    if (!this.props.showStats) {
      this._removeTrackerListeners(oldTracker);
    } else {
      this._removeTrackerListeners(oldTracker);
      if (tracker) {
        this._addTrackerListeners(tracker);
        // new stats!
        const statsContainer = this.refs.statsContainer;
        if (statsContainer) {
          this._stats = new Stats({ container: statsContainer, cssText: '' });
        } else {
          console.error('statsContainer is not available yet');
        }
      }
    }
    // Check the mediaSrc
    this._updateVideoSrc();
  }

  _shutdownMediaStream () {
    if (!this._mediaSrcStream) { return; }
    // Shut down video stream (if it has been set)
    const tracks = this._mediaSrcStream.getVideoTracks();
    tracks.forEach(track => track.stop());
    this._mediaSrcStream = null;
  }

  _updateVideoSrc () {
    const mediaType = this.props.mediaType;
    if (mediaType !== 'video') { return; }
    const videoEl = this.refs.media;
    const mediaSrc = this.props.mediaSrc;
    // Lazy check to see if it is a MediaStream
    if (videoEl && mediaSrc && mediaSrc.getVideoTracks) {
      this._mediaSrcStream = mediaSrc;
      setVideoSrc(videoEl, mediaSrc);
    } else {
      this._shutdownMediaStream();
      if (videoEl.src) {
        videoEl.removeAttribute('src');
      }
      if (mediaSrc) {
        videoEl.src = mediaSrc;
      }
    }
  }

  render () {
    const mediaType = this.props.mediaType;
    let media;
    if (mediaType === 'video') {
      media = (
        <video
          className='media'
          ref='media'
          autoPlay
          loop
          {...this.props.mediaSize}
        />
      );
    } else if (mediaType === 'image') {
      media = (
        <canvas
          className='media'
          ref='media'
          {...this.props.mediaSize}
        />
      );
    } else {
      throw new Error('unknown mediaType: ' + mediaType);
    }

    let deformerCanvas;
    if (this.props.hasDeformer) {
      deformerCanvas = (
        <canvas
          ref='deformerCanvas'
          className='overlay'
          {...this.props.mediaSize}
        />
      );
    }

    return (
      <div className='tracker-container-cmpt'>
        {media}

        {deformerCanvas}

        <canvas
          ref='canvas'
          className='overlay'
          {...this.props.mediaSize}
        />

        <div className='stats-container' ref='statsContainer' />
      </div>
    );
  }
}

TrackerContainer.propTypes = {
  mediaType: PropTypes.string.isRequired,
  mediaSrc: PropTypes.any,
  mediaSize: PropTypes.object.isRequired,
  showStats: PropTypes.bool,
  tracker: PropTypes.any,
  hasDeformer: PropTypes.bool
};
