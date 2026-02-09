
import numpy as np

xval =
yval =
rotation =
zstart = 1 # in um
zend = 100 # in um
zstep = 1 # um

nreps = 100 # number of repetitions

filename = "/path/to/file"


# open file as writeable
with open(filename, 'w'):
    for zval in np.arange(zstart, zend, zstep):
        string2write = strrep("\{
    "type": "operation",\n "data": {\n\t  "device": "stacks",
      "command": "add",
      "name": "stack_2",
      "elements": [
        {
          "name": "x",
          "start": xval,
          "end": xval
        },
        {
          "name": "y",
          "start": yval,
          "end": yval
        },
        {
          "name": "z",
          "start": zval,
          "end": zval
        },
        {
          "name": "r",
          "start": rot,
          "end": rot
        }
      ],
      "n": 1,
      "reps": %d
    }
    })", nreps)
        #write the string by appending it into the file

