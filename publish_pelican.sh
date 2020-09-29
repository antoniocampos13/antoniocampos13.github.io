#!/usr/bin/env bash

pelican content -o output -s publishconf.py
ghp-import -m "${1}" --no-jekyll -b master output
git push origin master

git add content
git commit -m "${1}"
git push origin content
