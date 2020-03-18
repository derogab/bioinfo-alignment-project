# bioinfo-alignment-project
### Program to print reads from a SAM/BAM file filtered for start positions.

```shell
$ python app.py -h
usage: app.py [-h] -f FILE [-p POS [POS ...]] [--only-pos]

Print reads from a SAM/BAM file filtered for start positions.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  path to input file
  -p POS [POS ...], --pos POS [POS ...]
                        list of positions
  --only-pos            get list of all align start positions
```

##### W/ Docker 
```shell
docker build -t bioinfo-alignment-project .
```
```shell
docker run \
  -v /host/path/to/examples/:/usr/src/app/examples/ \
  bioinfo-alignment-project
```