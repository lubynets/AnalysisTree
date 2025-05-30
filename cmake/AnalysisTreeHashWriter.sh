#!/bin/bash

FILE_HASH=AnalysisTreeHash.sh
FILE_DIFF=AnalysisTreeDiff.patch
for F in $FILE_HASH $FILE_DIFF; do
  if [ -f $F ]; then
    rm $F
  fi
done

SRC_DIR=${1}

cd $SRC_DIR
if [ -d ".git" ]; then
  GITTAG=$(git describe --tags)
  GITCOMMIT=$(git rev-parse HEAD)
  GITSTATUS=$(git status --porcelain)
  echo "export ANALYSIS_TREE_TAG=\"${GITTAG}\"" >> $FILE_HASH
  echo "export ANALYSIS_TREE_COMMIT_HASH=${GITCOMMIT}" >> $FILE_HASH
  if [ -z "${GITSTATUS}" ]; then
    echo "export ANALYSIS_TREE_COMMIT_ORIGINAL=TRUE" >> $FILE_HASH
  else
    echo "export ANALYSIS_TREE_COMMIT_ORIGINAL=FALSE" >> $FILE_HASH
    git diff >> $FILE_DIFF
  fi
else
  echo "export ANALYSIS_TREE_TAG=NOT_A_GIT_REPO" >> $FILE_HASH
  echo "export ANALYSIS_TREE_COMMIT_HASH=NOT_A_GIT_REPO" >> $FILE_HASH
  echo "export ANALYSIS_TREE_COMMIT_ORIGINAL=NOT_A_GIT_REPO" >> $FILE_HASH
fi
cd -
mv $SRC_DIR/$FILE_HASH .
if [ -f $SRC_DIR/$FILE_DIFF ]; then
mv $SRC_DIR/$FILE_DIFF .
fi
