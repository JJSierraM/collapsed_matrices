Collapsed matrices is a methodology that I developed for the niche case of having to calculate a great number of instances of a symmetrical function.
Say you have the N-body problem: you have to calculate the gravity casted between N celestial objects. You can loop over each of them and, for each one, loop thorugh every other object to calculate said gravity. That would mean N^2 calculations, but the gravity between A and B and between B and A is the same, so you are indeed doing twice 
as many operations than you could. Also, if you store said gravities in a 2-dimension array (NxN items) you are allocating over twice the needed memory (only (NxN-N)/2 is needed). 

Collapsed matrices makes the conversion between that square matrix and "collapses" it to a 1-dimension array through a custom data-structure that lets you:
-Allocate memory on a contiguous array, favoring data access by making it cache-friendly.
-Use exactly the amount of memory needed. 
-Easily distributable among threads, making it very easy to parallelize.

It works for any amount of items and with relations of any number of dimensions. Gravity is a funciton between 2 objects, but it can be used to calculate functions with any amount of dimensions, as long as the function is symmetrical between all of them (e.g. mean point of a set of three points).  

Under the Presentation/Presentation.pptx there is a PowerPoint presentation explaining the gist of it, although is thought to be _talked_, so written info on it is limited. Also, it is recommended to see it on presentation mode, as it has some animations that should make it a bit easier to catch the idea. 

The propper version is the one in the C folder. C++, CUDA and Rust are only poor implementations that I am developing as an exercise to practice in those languages. It is not advised to used those versions yet. The first one that I will trully develop to an end is the CUDA one, the rest will be left as academic exercises for myself. 

Right now, for 1024 objects (which means 523 776 gravitational relations), the N-body problem is solved in under 70 ms on an Intel© Core™ i7-6700HQ CPU @ 2.60GHz × 4 (2015 processor))
