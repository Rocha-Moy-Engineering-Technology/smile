name := "smile-core"

packageOptions += Package.ManifestAttributes("Automatic-Module-Name" -> "smile.core")

libraryDependencies ++= Seq(
  "org.bytedeco" % "javacpp"   % "1.5.9"        % "provided" classifier "macosx-x86_64" classifier "windows-x86_64" classifier "linux-x86_64",
  "org.bytedeco" % "openblas"  % "0.3.23-1.5.9" % "provided" classifier "macosx-x86_64" classifier "windows-x86_64" classifier "linux-x86_64",
  "org.bytedeco" % "arpack-ng" % "3.9.0-1.5.9"  % "provided" classifier "macosx-x86_64" classifier "windows-x86_64" classifier "linux-x86_64" classifier ""
)
