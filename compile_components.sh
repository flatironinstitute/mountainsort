#!/bin/bash

# verify if we're running Qt4 or Qt5 qmake
QT_VERSION=`qmake -query QT_VERSION`
if [[ "$QT_VERSION" =~ ^4.*$ ]]
then
  echo "You're trying to build MountainSort with Qt4 ($QT_VERSION) but MountainSort only supports Qt5"
  echo "Please make sure qmake from Qt5 installation is first in your PATH."
  
  # try to find a Qt5 install in $HOME or /opt
  if [ -d "/opt/Qt" ]
  then
    ls -1 /opt/Qt/|grep -sqe '^5\.[[:digit:]]$'
    if [[ "$?" == "0" ]]
    then
      echo "It looks like Qt5 installation(s) might be present in /opt/Qt/"
      echo
      exit 1
    fi
  fi
  if [ -d "$HOME/Qt" ]
  then
    ls -1 "$HOME/Qt"|grep -sqe '^5\.[[:digit:]]$'
    if [[ "$?" == "0" ]]
    then
      echo "It looks like Qt5 installation(s) might be present in $HOME/Qt"
      echo
      exit 1
    fi
  fi
  echo "If you don't have Qt5 installed, please go to http://qt.io and download Qt5"
  echo
  exit 1
fi

args0=""
components0=""
for var in "$@"; do
	if [[ "$var" == "nogui" ]]; then
		args0="$args0 \"GUI = off\""
	elif [[ "$var" == "default" ]]; then
		echo ""
	else
		components0="$components0 $var"
	fi

    [ "$var" == 'nogui' ] && ARGS+=("$var")
done

qmake -recursive

if [ -z "$components0" ]; then
eval qmake $args0
else
eval qmake \"COMPONENTS = $components0\" $args0
fi

NPROCCMD=$(which nproc)

if [ -z "$NPROCCMD" ]
then
  $NPROCCMD="lscpu|awk '/^CPU\(s\)/ { print $2; }'"
fi

NPROC=$($NPROCCMD)

if [ -z "$NPROC" ]
then
  NPROC="2"
fi
echo "Building with $NPROC parallel jobs"
make -j $NPROC
EXIT_CODE=$?
if [[ $EXIT_CODE -ne 0 ]]; then
	echo "Problem in compilation."
else
	make install #not strictly needed at this point
	echo ""
	echo "Compilation successful."
fi

