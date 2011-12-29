To the user of this Library,

Once one have obtained this library by

  svn checkout http://num4lcp.googlecode.com/svn/trunk/ num4lcp

or similar. Then this should create a num4lcp folder on ones local machine. Inside this one should find a matlab folder containg the whole library.

When running the code one matlab function/script may call another function inside the svn repository. Thus one should make sure that Matlab currents folder is changed to the newly checked out Matlab folder or that the newly checked out Matlab folder is added to the search path of Matlab. If this is unfamiliar the one may try to write the following in the Matlab prompt

   addpath /my/full/path/to/check-out/destination/num4lcp/matlab

That should take care of the path issues for the current matlab session.

The repository is organized in a flat structure where functionality is divided into individual files. All files named:

* make_XXX.m are factory functions used to generate matrices and LCP problems

* test_XXX.m are all test scripts and the XXX part of the name indicates which method or properties that are being tested.

* All remaining matlab files are functions implementing the numerical methods of the library. There should be plenty of inline code comments to help get an overview of the individual implementations.

The test scripts are a good place to start looking for boiler-plate code for how one can make calls for the numerical methods in the library in ones own Matlab code.

Enjoy

Kenny Erleben (2012-12-29)
