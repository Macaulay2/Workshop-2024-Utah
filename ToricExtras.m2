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
    AuxiliaryFiles => false,
    Version => "0.1",
    Date => "1 June 2024",
    Authors => {{
        Name => "Gregory G. Smith", 
        Email => "ggsmith@mast.queensu.ca", 
        HomePage => "http://www.mast.queensu.ca/~ggsmith"}},
    Headline => "new routines for working with normal toric varieties",
    Keywords => {"Toric Geometry"},
    PackageExports => {"NormalToricVarieties"},
    PackageImports => {},
    DebuggingMode => true
    )

export {   
    }


------------------------------------------------------------------------------
-- CODE
------------------------------------------------------------------------------


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

------------------------------------------------------------------------------
-- TESTS
------------------------------------------------------------------------------


end---------------------------------------------------------------------------     

------------------------------------------------------------------------------
-- SCRATCH SPACE
------------------------------------------------------------------------------

-- XXX
uninstallPackage "ToricExtras";
restart
installPackage "ToricExtras"
check NormalToricVarieties

needsPackage "ToricExtras";
