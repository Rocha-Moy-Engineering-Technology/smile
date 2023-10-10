#!/bin/sh

while true; do
    read -r -p "Do you want to publish it? " ans
    case $ans in
    [Yy]*)
        sbt publishSigned
        check_error "!!"

        git checkout scala3
        check_error "!!"
        git pull
        check_error "!!"
        git merge master
        check_error "!!"
        sbt ++3.3.0 scala/publishSigned
        check_error "!!"
        sbt ++3.3.0 json/publishSigned
        check_error "!!"
        # sbt ++3.3.0 spark/publishSigned
        # check_error "!!"

        # git checkout master
        break
        ;;
    [Nn]*) break ;;
    *) echo "Please answer yes or no." ;;
    esac
done

# while true; do
#     read -p "Do you want to publish smile-kotlin? " ans
#     case $ans in
#     [Yy]*)
#         cd kotlin
#         gradle publishMavenJavaPublicationToMavenRepository
#         check_error "gradle publish"
#         cd ..
#         break
#         ;;
#     [Nn]*) break ;;
#     *) echo "Please answer yes or no." ;;
#     esac
# done

# while true; do
#     read -p "Do you want to publish smile-clojure? " ans
#     case $ans in
#     [Yy]*)
#         cd ../clojure
#         lein deploy clojars
#         check_error "lein deploy clojars"
#         cd ..
#         break
#         ;;
#     [Nn]*) break ;;
#     *) echo "Please answer yes or no." ;;
#     esac
# done
