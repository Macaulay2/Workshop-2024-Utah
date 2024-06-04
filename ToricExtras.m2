-- -*- coding: utf-8 -*-
------------------------------------------------------------------------------
-- Copyright 2009--2020 Gregory G. Smith
--
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- more details.
--
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------------
newPackage(
    "ToricExtras",
    AuxiliaryFiles => true,
    Version => "0.1",
    Date => "1 June 2024",
    Authors => {{
            Name => "Gregory G. Smith",
            Email => "ggsmith@mast.queensu.ca",
            HomePage => "http://www.mast.queensu.ca/~ggsmith"},
        {
            Name => "Zhengning Hu",
            Email => "zhengninghu@arizona.edu"},
        {
            Name => "Prajwal Samal",
            Email => "psamal@impan.pl"},
        {
            Name => "Griffin Edwards",
            Email => "griffinedwards@gatech.edu"},
        {
            Name => "Noah Solomon",
            Email => "nsolomon33@gatech.edu"},
        {
            Name => "Jay Yang",
            Email => "jayy@wustl.edu"},
        {
            Name => "Anna Natalie Chlopecki",
            Email => "achlopec@purdue.edu",
	    HomePage => "https://anchlopecki.github.io"},
        {
            Name => "Eduardo Torres Dávila",
            Email => "torre680@umn.edu"},
        {
            Name => "Nate Gallup",
            Email => "npgallup@gmail.com",
	    HomePage => "https://sites.google.com/view/nathanielgallup/home"},
        {
            Name => "Sean Grate",
            Email => "sean.grate@auburn.edu",
            HomePage => "https://seangrate.com/"},
        {
            Name => "Thiago Holleben",
            Email => "hollebenthiago@dal.ca",
            HomePage => "https://hollebenthiago.github.io/"},
	{
            Name => "Will Gilroy",
            Email => "wlg38@cornell.edu"
	},
        {
	    Name => "Rohan Joshi",
	    Email => "rohansjoshi@math.ucla.edu"}
        },

    Headline => "new routines for working with normal toric varieties",
    Keywords => {"Toric Geometry"},
    PackageExports => {"NormalToricVarieties"},
    PackageImports => {"PrimaryDecomposition"},
    DebuggingMode => true
    )

export {
    "ToricLinearSeries",
    "toricLinearSeries",
    "baseLocusIdeal",
    "isBasepointFree",
    "batyrevConstructor",
    "toricMap", 
    "projectivizationOfBundle"
    }


------------------------------------------------------------------------------
-- CODE
------------------------------------------------------------------------------

load "./ToricExtras/ToricLinearSeries.m2"

load "./ToricExtras/BatyrevConstructions.m2"

load "./ToricExtras/ProjectiveBundlesDivisors.m2"

------------------------------------------------------------------------------
-- DOCUMENTATION
------------------------------------------------------------------------------
beginDocumentation ()
doc ///
    Key
        ToricExtras
    Headline
        new features for normal toric varieties
    Description
    	Text
	    This temporary package implements several new features that will
	    be incorporated into the existing NormalToricVarieties package.
///

doc ///
    Key
        ToricLinearSeries
    Headline
        Linear series on a toric variety
    Description
    	Text
	        Work in progress implementation of linear series on a toric variety
///

load "./ToricExtras/BatyrevDocs.m2"

------------------------------------------------------------------------------
-- TESTS
------------------------------------------------------------------------------

-- test 0
TEST ///
    X = toricProjectiveSpace 1;
    assert isWellDefined X
///

load "./ToricExtras/ToricLinearSeriesTests.m2"
load "./ToricExtras/BatyrevTests.m2"
load "./ToricExtras/ProjectiveBundlesDivisorsTests.m2"

end---------------------------------------------------------------------------

------------------------------------------------------------------------------
-- SCRATCH SPACE
------------------------------------------------------------------------------

-- XXX
uninstallPackage "ToricExtras";
restart
installPackage "ToricExtras"
check ToricExtras
check (ToricExtras, Verbose => true)

needsPackage "ToricExtras";

-----------------------------------------
--------- ANNA'S SCRATCH SPACE ----------
-----------------------------------------

restart
needsPackage "NormalToricVarieties";
load "Desktop/Projects/Workshop-2024-Utah/ToricExtras/BatyrevConstructions.m2";

allGood = (V) -> (
    i := isWellDefined V;
    j := isSmooth V;
    k := isProjective V;
    l := rank picardGroup V == 3;
    i and j and k and l
    )
