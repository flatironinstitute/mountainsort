# MountainSort

Spike sorting software

## Which version to use? 

Depending on your situation, you will want to use one of the following two versions of MountainSort/MountainLab:

**Option 1. The current, stable (and now frozen) version.**

For documentation and installation instructions for this version of MountainSort, MountainLab, and related tools, please visit [MountainSort on *Read the Docs*](https://mountainsort.readthedocs.org).

**Option 2. The new beta version, under development**

See [mountainlab-js](https://github.com/flatironinstitute/mountainlab-js)

Here is a list of pros and cons for the two methods, which can help you decide which to use.

Pros of option 1:
* Tested and stable
* Relatively complete documentation and tutorial
* Stable user interface (MountainView) for exploring results
* Install via apt (ppa) if on Ubuntu 16.04

Cons of option 1:
* Sorting algorithm not optimized for large electrode arrays
* Does not run on Mac
* Will no longer have bug fixes

Pros of option 2
* Better implementation of MountainSort -- can handle large electrode arrays
* More portable -- install on any Linux flavor or Mac
* Better scripting of pipelines (although not fully documented)
* Active bug fixes
* Cleaner code base, ultimately better documentation, and overall better software.

Cons of option 2
* In flux, under development, not fully documented
* Curation step not yet implemented for spike sorting algorithm
* User interface (MountainView) not yet available (although one can simply use the MountainView from option 1 -- it is compatible)

**So, in summary, which one should you use?**

Answer: ultimately you'll want to go with option 2, but depending on your situation, you may want to wait until it is more mature. If you have relatively small electrode arrays (say 32 channels or less per shank/probe), then the first version is your best bet at this point.

**Note:** Many of the components between the two options can be mixed and matched. But I don't want to make the docs too confusing.

