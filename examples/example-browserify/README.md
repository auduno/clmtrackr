Browserify Example
===================

This example will include the source files in "js" using browserify, rather than including the fully built file.

I set this up to use a prototyping tool called [budo](https://github.com/mattdesl/budo), but you could also use similar tools like [beefy](https://github.com/chrisdickinson/beefy) and [wzrd](https://github.com/maxogden/wzrd).

To run this example using budo, first make sure you have it install globally:

    npm install budo -g

Then, with the "example-browserify" folder as your working directory, run:

    budo index.js --live

Finally head to [http://127.0.0.1:9966/](http://127.0.0.1:9966/) in your browser.