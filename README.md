pyRCL
=====

Read and convert .rcl HF Radar files to .png images.

## Install

#### From source

```bash
git clone https://github.com/axiomdatascience/pyRLC
cd pyRLC
python setup.py install
```

## Test
```bash
python setup.py test
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
