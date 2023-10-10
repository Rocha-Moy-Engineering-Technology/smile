#!/bin/sh

# In a script, history expansion is turned off by default, enable it with
set -o history -o histexpand

check_error() {
    local retval=$?
    if [ $retval -ne 0 ]; then
        echo "'$1' returns code $retval"
        exit $retval
    fi
}

if [[ "$OSTYPE" == "darwin"* ]]; then
    JAVA_HOME=$(/usr/libexec/java_home -v 1.8)
    export JAVA_HOME
fi

sbt clean
rm -rf doc/api/*
sbt unidoc
check_error "!!"
mv target/javaunidoc doc/api/java

sbt json/doc
check_error "!!"

sbt scala/doc
check_error "!!"

#cd kotlin
#gradle dokkaHtml
#check_error "!!"

#cd ../clojure
#lein codox
#check_error "!!"

#cd ../web
#npx @11ty/eleventy
#check_error "!!"

# cd ..
sbt universal:packageBin
check_error "!!"
