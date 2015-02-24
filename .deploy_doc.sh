#!/bin/bash

git clone git@github.com:cosmicpy/cosmicpy.github.io.git

(cd cosmicpy.github.io
	git config user.name "CosmicPy Travis CI Worker"
	git config user.email "travis@cosmicpy.github.io"
	cp -r $HOME/build/cosmicpy/cosmicpy/doc/build/html/* .
	git add *
	git commit -m "Deployed documentation from Travis"
	git push --quiet git@github.com:cosmicpy/cosmicpy.github.io.git > /dev/null 2>&1
)
rm -rf cosmicpy.github.io
