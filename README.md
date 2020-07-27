# dtsa2Scripts

This repository contains exemplar j/python scripts demonstrating how
to accomplish common tasks with
[DTSA-II](http://www.cstl.nist.gov/div837/837.02/epq/dtsa2/index.html),
developed by Nicholas Ritchie at NIST. He describes DTSA-II:

> DTSA-II is a multiplatform software package for quantitative x-ray
> microanalysis. ... DTSA-II has been designed with the goal of making
> standards-based microanalysis more accessible for the novice
> microanalyst. We want to encourage standards-based analysis by
> making it as easy as possible to get reliable results.

P/Jython scripts are useful for automating common tasks, making them
reproducible. Because these files are plain text, they work well with
version control systems like `git`.

The `mc3Scripts` folders contain exemplar scripts for modelling spectrum
generation using the classes/functions in
`gov.nist.microanalysis.NISTMonte.Gen3`. Many convenient wrapper
functions are exported in `mcSimulate3.py` in the `$DTSA_ROOT/Lib/dtsa2`
folder. These are imported in a script with:

```
import dtsa2.mcSimulate3 as mc3
```

