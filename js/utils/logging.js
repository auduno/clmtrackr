let LOGGING_ENABLED = true;


/**
 * Wrapped logging function.
 * @param {string} msg The message to log.
 */
export const log = function (msg) {
  if (!LOGGING_ENABLED) { return; }
  if (window.console && window.console.log) {
    window.console.log(msg);
  }
};


/**
 * Wrapped logging function.
 * @param {string} msg The message to log.
 */
export const error = function (msg) {
  if (!LOGGING_ENABLED) { return; }
  if (window.console) {
    if (window.console.error) {
      window.console.error(msg);
    } else if (window.console.log) {
      window.console.log(msg);
    }
  }
  throw msg;
};


/**
 * Turn off all logging.
 */
export const loggingOff = function () {
  LOGGING_ENABLED = false;
};
