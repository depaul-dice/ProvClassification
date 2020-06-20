Using prov2dot and DFGen to generate a first-pass scan of dependencies:

prov2dot.py:
```
usage: prov2dot.py [-h] [--nosub] [--nofilter] [--withfork] [-f FIN_NAME]
                   [-d DIR_NAME] [--withgraph]

Process provenance log file.

optional arguments:
  -h, --help   show this help message and exit
  --nosub
  --nofilter
  --withfork
  -f FIN_NAME
  -d DIR_NAME
  --withgraph

```

DFGen:
```
USAGE: dfgen /path/to/prov.log
```
