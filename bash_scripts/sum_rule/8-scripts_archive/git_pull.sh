git pull git@github.com:amirhs1/sumrule_pipeline.git main # pull from repository
git clone git@github.com:amirhs1/sumrule_pipeline.git main # clone a repo for the first time
rsync -axvH --no-g --no-p --exclude='.*' PACKAGE_NAME /destination # copy directory exclude hidden files