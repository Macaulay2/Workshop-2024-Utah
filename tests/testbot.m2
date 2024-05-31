-- This is the default testbot.m2 provided for this workshop.

-- It contains Macaulay2 code and is automatically run every time you open
-- a pull request or push a new commit to a project branch on the GitHub
-- respository, unless you add the string "[skip ci]" to the commit message.

-- If you have the software Docker installed and would like to run the same
-- tests locally or run M2 inside a container, use one of these commands from
-- the top directory of the workshop repository:
--  docker run -v "`pwd`":"/home/macaulay" ghcr.io/macaulay2/testbot --script tests/testbot.m2
--  docker run -v "`pwd`":"/home/macaulay" -it --entrypoint bash ghcr.io/macaulay2/testbot:latest

-- Uncomment and edit the following line to add your project directories
-- containing Macaulay2 source code files to the load path. Terminate each
-- directory name with a "/".
path = join( { currentDirectory() | "CYToolsM2" }, path )

-- Uncomment and edit the following lines to preload and check your package or
-- to run a series of examples with every push on GitHub.
--needsPackage "LocalRings"
--check LocalRings
--load "tests/example.m2"
--capture get "tests/example.m2"

installPackage("StringTorics", FileName => currentDirectory() | "CYToolsM2/StringTorics.m2")
check StringTorics

-- The following lines automatically run every file in the "tests" directory.
-- If you wish, you can change testDir to any other directory.
testDir = currentDirectory() | "tests/"

testFiles = select(readDirectory testDir,
    file -> match("\\.m2$", file) and file != "testbot.m2")
printerr("Found ", toString(#testFiles), " test file(s) matching '", testDir, "*.m2'.")
TEST(testFiles / (filename -> testDir | "/" | filename), FileName => true)

-- workaround for https://github.com/Macaulay2/M2/issues/2835
importFrom_Core {"PackageIsLoaded"}
User.PackageIsLoaded = true

check(User, Verbose => true)
