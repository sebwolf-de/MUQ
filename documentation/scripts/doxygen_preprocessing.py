import sys
import subprocess
import os

"""
This script is run during the "make doc" build just before calling doxygen.  It
can be used for generating doxygen input files.
"""

# Read the root and build directoriess
assert (len(sys.argv)>2)

baseDir = sys.argv[1]
if(baseDir[-1]=='/'):
    baseDir = baseDir[0:-1]

buildDir = sys.argv[2]
if(buildDir[-1]=='/'):
    buildDir = buildDir[0:-1]

# Create the output directory if it doesn't exists
if not os.path.exists(buildDir+'/doxygen_prep'):
    os.makedirs(buildDir+'/doxygen_prep')


def MakeInstallation():
    """
    Combines markdown files in "doxFiles/Installation" with the bash script in
    "SupportScripts" to create a single set of installation instructions.
    """
    md1 = baseDir + '/documentation/doxFiles/excluded/source_install1.md'
    md2 = baseDir + '/documentation/doxFiles/excluded/source_install2.md'
    bash = baseDir + '/SupportScripts/install-instructions.sh'

    newText = subprocess.run(["cat", md1, bash, md2], stdout=subprocess.PIPE)
    print(newText.stdout)
    with open(buildDir+'/doxygen_prep/source_install.md', 'wb') as f:
        f.write(newText.stdout)


MakeInstallation()
