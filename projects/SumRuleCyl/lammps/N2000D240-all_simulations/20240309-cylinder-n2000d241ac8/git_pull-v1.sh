#!/bin/bash

git pull git@github.com:amirhs1/PROBING_PACKAGE.git master # pull from repository
git clone git@github.com:amirhs1/PROBING_PACKAGE.git master # clone a repo for the first time
rsync -axvH --no-g --no-p --exclude='.*' PROBING_PACKAGE /destination # copy directory exclude hidden files