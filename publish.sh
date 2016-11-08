#!/bin/bash
set -e
# Remove existing dist
npm run clean
# Run build
npm run build

# Publish using gh-pages tool
gh-pages -d dist
