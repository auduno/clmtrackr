/**
 * Check if the page is embedded.
 * @return {boolean} True of we are in an iframe
 */
export const isInIFrame = function () {
  return window !== window.top;
};


export const updateCSSIfInIFrame = function () {
  if (isInIFrame()) {
    document.body.className = 'iframe';
  }
};
