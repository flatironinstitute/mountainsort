# MountainSort

Spike sorting software

Please use the [new version of MountainSort](https://github.com/flatironinstitute/mountainsort_examples/blob/master/README.md).

## Which version to use? 

Depending on your situation, you will want to use one of the following two versions of MountainSort/MountainLab:

**Option 1. The old (and now frozen) version.**

For documentation and installation instructions for the old version of MountainSort, MountainLab, and related tools, please visit [MountainSort on *Read the Docs*](https://mountainsort.readthedocs.org).

**Option 2. The new version, under active development**

See [the new version of MountainSort](https://github.com/flatironinstitute/mountainsort_examples/blob/master/README.md)

Here is a list of pros and cons for the two methods, which can help you decide which to use.

Pros of option 1:
* Relatively complete documentation and tutorial
* Install via apt (ppa) if on Ubuntu 16.04

Cons of option 1:
* Will no longer have bug fixes
* Sorting algorithm not optimized for large electrode arrays
* Does not run on Mac

Pros of option 2
* Better implementation of MountainSort -- can handle large electrode arrays
* More portable -- install on any Linux flavor or Mac
* Integration with python and JupyterLab
* Active bug fixes
* Cleaner code base, ultimately better documentation, and overall better software.

Cons of option 2
* It is under development

**So, in summary, which one should you use?**

Answer: ultimately you'll want to go with option 2.

**Note:** Many of the components between the two options can be mixed and matched. But I don't want to make the docs too confusing.

