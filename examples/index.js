// See: http://stackoverflow.com/questions/33541949/react-material-ui-dropdown-menu-not-working
import injectTapEventPlugin from 'react-tap-event-plugin';
injectTapEventPlugin();
// ----

import React from 'react';
import ReactDOM from 'react-dom';

import App from './components/App';

import { createStore, applyMiddleware } from 'redux';
import { Provider } from 'react-redux';
import createLogger from 'redux-logger';
import reducers from './reducers';

import MuiThemeProvider from 'material-ui/styles/MuiThemeProvider';


let store;
if (process.env.NODE_ENV === 'production') {
  store = createStore(reducers);
} else {
  store = createStore(reducers, applyMiddleware(createLogger()));
}


// getUserMedia only works over https in Chrome 47+, so we redirect to https.
// Also notify user if running from file.
if (window.location.protocol === 'file:') {
  alert('You seem to be running this example directly from a file. Note that these examples only work when served from a server or localhost due to canvas cross-domain restrictions.');
} else if (
  window.location.hostname !== 'localhost' &&
  window.location.protocol !== 'https:'
) {
  window.location.protocol = 'https';
}


ReactDOM.render(
  <Provider store={store}>
    <MuiThemeProvider>
      <App />
    </MuiThemeProvider>
  </Provider>,
  document.getElementById('main')
);
