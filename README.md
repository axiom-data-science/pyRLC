pyRLC
=====
[![Build Status](https://travis-ci.org/axiom-data-science/pyRLC.svg?branch=master)](https://travis-ci.org/axiom-data-science/pyRLC)

Read and convert HF Radar RLC files to PNG images.

## Install

#### From source

Installs library and command line tool `rec_to_png`

```bash
git clone https://github.com/axiom-data-science/pyRLC
cd pyRLC
python setup.py install
```

## Test
```bash
python setup.py test
```

With coverage (missing bits are due to untested exception paths):

```bash
python -m pytest tests --cov=pyRLC
```

## Usage

Convert a single .rec file to .png and add overlay

```python
from pyRLC import pyRLC
with pyRLC.RLCfile('test.rec') as f:
    f.rec_to_png()

radar = pyRLC.RadarOverlay('radar.png', 'overlay.png', 'new.png')
radar.add_overlay()
```

Using command line script to convert single .rec to .png:

```rec_to_png radar.png output/```

Using command line script to convert single .rec to .png with overlay:

```rec_to_png radar.png output/ -o overlay.png```

