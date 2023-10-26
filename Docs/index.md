<!-- Phys 581 - The Standard Model of Particle Physics
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
-->


<!-- Include ../README.md
     If you would like to use the contents of your top-level README.md file here, then
     you can literally include it here with the following:

```{include} ../README.md
``` 

    Note that this may will break `sphinx-autobuild` (`make doc-server`) which will not rebuild
    this index file when ../README.md changes.  See the note at the bottom of the file
    if you want to do this while using sphinx-autobuild.
--> 


Welcome to Phys 581 - The Standard Model of Particle Physics!  This is the main
documentation page for the course.  For more class information, please see the
{ref}`sec:sylabus`.

The following topics are relatively complete and should be reviewed.  More material in a
somewhat incomplete state is available under Class Notes in the sidebar.

```{toctree}
---
maxdepth: 2
caption: "Renormalization Group:"
titlesonly:
---
ClassNotes/RenormalizationGroup
ClassNotes/RandomWalks
```

```{toctree}
---
maxdepth: 2
caption: "Resources:"
titlesonly:
hidden:
glob:
---
Syllabus
Prerequisites/*
References
```

```{toctree}
---
maxdepth: 2
caption: "Class Notes:"
titlesonly:
hidden:
glob:
---
ClassNotes/*
```

```{toctree}
---
maxdepth: 2
caption: "Miscellaneous:"
hidden:
---
Demonstration
CoCalc
ClassLog
../InstructorNotes

README.md <../README>
```

<!-- If you opt to literally include files like ../README.md and would like to be able
     to take advantage of `sphinx-autobuild` (`make doc-server`), then you must make
     sure that you pass the name of any of these files to `sphinx-autobuild` in the
     `Makefile` so that those files will be regenerated.  We do this already for
     `index.md` but leave this note in case you want to do this elsewhere.
     
     Alternatively, you can include them separately and view these directly when editing.
     We do not include this extra toc when we build on RTD or on CoCalc.  We do this
     using the `sphinx.ext.ifconfig extension`:
     
     https://www.sphinx-doc.org/en/master/usage/extensions/ifconfig.html

```{eval-rst}
.. ifconfig:: not on_rtd and not on_cocalc

   .. toctree::
      :maxdepth: 0
      :caption: Top-level Files:
      :titlesonly:
      :hidden:

      README.md <../README>
      InstructorNotes.md <../InstructorNotes>
```
-->
